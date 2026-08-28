# =============================================================================
# Meta-analysis: effect of LLM assistance on clinicians' diagnostic accuracy
# Primary analysis script (metafor + clubSandwich - "harmonized" effect size pooling)
# run the same script with  AGGREGATION_MODE = "granular" to get Fig. S5)
# -----------------------------------------------------------------------------
# Estimand: Hedges' g_av (paired rows via metafor SMCRPH; independent rows via
#   SMD/SMDH). Positive g = benefit of LLM assistance.
# Working model: 3-level random effects, effects within experiments
#   ( ~ 1 | Experiment_ID/Effect_ID ), REML.
# Inference: cluster-robust (CR2 / RVE) at the PAPER / shared-case-corpus level
#   (ID_2); the overlapping TUM brain-MRI studies are merged into one ID_2 cluster.
#
# Pipeline (top -> bottom):
#   1. Packages, paths, data import, ID handling + TUM cluster merge
#   2. Harmonized pre-ES pooling of within-paper reader/strata rows
#      (raw-arm reconstruction) + candidate audit + homogeneity guard
#   3. Effect-size + sampling-variance construction via metafor (escalc/conv.*),
#      within-pair r handling, reader-clustering design effect (ICC)
#   4. Main multilevel model + CR2 robustification; forest plots
#   5. r- and ICC-sensitivity; heterogeneity (multilevel I^2)
#   6. Robustness: low-RoB subset, publication status, study design,
#      4-level working model
#   7. Moderators: effect-coded joint model (supplement) + no-intercept
#      cell-means with CR2 + HTZ omnibus (main reporting)
#   8. Small-study effects (funnel, Egger/PET), leave-one-out (effect /
#      experiment / study), LaTeX export
# NOTE: raw coding file is treated as read-only source of truth; every pooling,
#   transformation and exclusion is audited to a CSV in output_dir.
# =============================================================================

# Prepare workspace
rm(list = ls())

# clubSandwich underpins all CR2/Satterthwaite inference (vcovCR, coef_test,
# Wald_test) and is also what metafor::robust(clubSandwich = TRUE) calls
# internally; readxl reads the PROBAST workbook. knitr provides kable(),
# kableExtra the LaTeX styling.
packages <- c(
  "metafor", "clubSandwich", "tidyverse", "knitr", "kableExtra",
  "rstudioapi", "readr", "rlang", "dplyr", "readxl"
)

install_and_load <- function(pkg) {
  # paired d_av-style effect sizes are computed by metafor rather than manually.
  if (pkg == "metafor") {
    need_install <- !requireNamespace(pkg, quietly = TRUE) ||
      utils::packageVersion("metafor") < "4.6-0"
    if (need_install) install.packages("metafor")
  } else if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(packages, install_and_load))

# Project directories
if (rstudioapi::isAvailable()) {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  project_root <- dirname(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

data_dir   <- file.path(project_root, "Data")
output_dir <- file.path(project_root, "Plots")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Load data
CSV_ENCODING <- "UTF-8"
df_performance <- read_delim(
  file.path(data_dir, "Coding_Scheme_Final.csv"),
  delim = ";",
  quote = "\"",
  escape_double = TRUE,
  locale = locale(decimal_mark = ".", grouping_mark = "", encoding = CSV_ENCODING),
  na = c("NR", "n.d.", "", "NA"),
  show_col_types = FALSE
)

# Remove completely empty CSV records, such as trailing delimiter-only rows
df_performance <- df_performance %>%
  dplyr::filter(
    dplyr::if_any(dplyr::everything(), ~ !is.na(.x))
  )

# Clean names / standard IDs --------------------------------------------------
names(df_performance) <- trimws(names(df_performance))

df_performance <- df_performance %>%
  dplyr::rename(Experiment_ID = Exp_ID)

# Inner grouping factor of the ~ 1 | Experiment_ID/Effect_ID model. The unique
# effect-size id (ES_ID) is used as the inner factor, so each effect size is its
# own innermost random-effect unit; ES_ID is also retained for the forest-plot
# slabs. (The coding sheet's native 'Effect_ID' is the treatment-arm number and
# is deliberately not used here to avoid collapsing distinct effects within an
# experiment.)
df_performance <- df_performance %>%
  dplyr::mutate(Effect_ID = ES_ID)

# Merge the overlapping TUM brain-MRI corpus into ONE analysis cluster (ID_2 = 11):
#   ID_2 = 11  Kim, "Human-AI Collaboration ... Usability Study"
#   ID_2 = 12  Kim, "Boosting LLM-assisted diagnosis ..." (shares 20/30 cases)
#   ID_2 = 21  Schramm (shares the same TUM brain-MRI case corpus)
# All three draw on the same (partly shared) case corpus and are therefore treated
# as a single cluster for CR2 robustification. Only ID_2 is merged; Experiment_ID 
# and ES_ID stay separate so the within-cluster structure is preserved.
#
# Study_ID retains the distinct study identity for leave-one-STUDY-out analyses.
# ID_2 represents the paper/shared-case-corpus cluster used for CR2 inference
# and the PROBAST join. It enters the random-effects structure only in the
# separate four-level sensitivity model.
df_performance$Study_ID <- df_performance$ID_2
df_performance$ID_2[df_performance$ID_2 %in% c(11, 12, 21)] <- 11

# Helpers

norm <- function(x) {
  str_to_lower(str_squish(str_replace_all(as.character(x), "\u00A0", " ")))
}

# Cast numeric columns (only if present)
numeric_cols <- intersect(
  c(
    "M_treatment","SD_treatment","SE_treatment",
    "M_control","SD_control","SE_control",
    "n_treatment","n_control","n_total",
    # Additional raw-count fields used by the harmonized pooling step.
    "N_total","N_physicians_total","N_cases_total","N_assessments_total",
    "n_control_case","n_control_physician","n_treatment_case","n_treatment_physician",
    "h_n_pairs","h_r_used","h_b","h_c",
    "Discordant_b","Discordant_c",
    "Performance_Difference", "Performance_Difference_SD", "Performance_Difference_SE",
    "95CI_lower_performance_diff", "95CI_upper_performance_diff", "Performance_difference_p_value",
    "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment",
    "F_value","Chi_square","unstandardised_beta",
    "SE_beta","Z_value","t_value","Odds_Ratio",
    # OR fields are retained for audit purposes only and are not used by the primary effect-size router; 
    # h_icc is the optional per-study ICC column (absent -> anchor ICC used).
    "95CI_lower_odds_ratio","95CI_upper_odds_ratio","h_icc"
  ),
  names(df_performance)
)

df_performance <- df_performance %>% mutate(across(all_of(numeric_cols), as.numeric))

# =============================================================================
# SEMI-AUTOMATED HARMONIZED AGGREGATION
# =============================================================================
# Purpose:
#   - Keep the granular coding file unchanged as the source of truth.
#   - Automatically FLAG additional within-paper rows that share all moderators.
#   - Apply pooling ONLY to the study-specific rules agreed a priori:
#       Wang/RetiZero: 19 reader rows -> 2 career-stage rows
#       Ji:            14 reader/time rows -> 4 career-stage x time rows
#       McDuff:         8 specialty/outcome rows -> 4 medical-field x outcome rows
#       Chen:           3 case-subset rows -> 1 overall row
#       Siepmann:       6 experience/outcome rows -> 2 outcome rows
#     Further shared-case reader-stratification rules (same mechanism: readers
#     evaluate the SAME cases in the SAME condition, rows differ ONLY by reader
#     identity). Grouping keys use content columns only (Title + Career_stage /
#     Type_of_Results / N_cases_total / Time_lag_sessions), so they are robust to
#     any ID_2 renumbering:
#       Wang_ovarian:    5 reader rows        -> 1 row              (Doctors A-E, 77 cases)
#       Xu_pGGN:         6 reader rows        -> 2 career-stage rows (R7-R12, 110 cases)
#       Jin:             6 reader rows        -> 2 career-stage rows (junior/senior 1-3, 103 cases)
#       Sheng:          12 reader rows        -> 4 career x outcome rows (R1-R6 x top-1/top-3, 228 cases)
#       Camur:           4 reader rows        -> 2 condition rows   (R1+R2 within browser / optimal, 86 cases)
#       Zhang_ovarianCT: 4 reader rows        -> 3 rows: R4+R5 pooled (attending radiologists, 230 cases);
#                        G4 (attending, 224) and G5 (resident, 224) stay separate -- pooling them would
#                        collapse the Career_stage moderator and disjoint case cohorts (230 vs 224).
#   - Pool RAW arm-level information before Hedges' g / sampling variances are
#     calculated. No already calculated effect sizes are averaged.
#
# Toggle between the granular (one row per source observation - only used for robustness check in Figure S5) and the harmonized
# (pooled) analysis (main model, shown in Figure 1). The candidate checks below run in both modes.
AGGREGATION_MODE <- "harmonized"  # allowed: "granular", "harmonized"
stopifnot(AGGREGATION_MODE %in% c("granular", "harmonized"))

# Analysis moderators. Defined once here and reused by both the pooling
# homogeneity check and the meta-regression, so the two can never drift apart.
ANALYSIS_MODERATORS <- c(
  "Task_Format", "Explanation_Format_of_LLM", "LLM_Model",
  "Input_Form_Patient_Information_in_LLM", "Career_stage",
  "Type_of_Results", "Medical_Field", "Follow_Up_Prompts_Allowed",
  "Case_setting"
)
mods_pool_check <- intersect(ANALYSIS_MODERATORS, names(df_performance))

# Variables that are NOT moderators but can still indicate substantively
# different observations that should not be pooled blindly (e.g. validation vs
# test set, or Ji 2-week vs 3-month measurements).
pool_guard_vars <- intersect(
  c(
    "Test_Data", "Time_lag_sessions", 
    "h_design", "Same_Cases_Both_Conditions",
    "Control_group_resources", "Treatment_group_resources",
    # Distinguishes, for example, the two Cesur experiments (multimodal vs text-only clinician input).
    "Input_Form_Patient_Information_clinician"
  ),
  names(df_performance)
)

.pool_as_num <- function(x) suppressWarnings(as.numeric(x))

.pool_norm_title <- function(x) {
  stringr::str_to_lower(
    stringr::str_squish(
      stringr::str_replace_all(as.character(x), "\u00A0", " ")
    )
  )
}

# Study-specific rules. Matching is done on normalized titles so that minor
# whitespace / capitalization differences in the CSV do not matter.
.assign_pool_rule <- function(df) {
  title_norm <- .pool_norm_title(df$Title)

  df %>%
    dplyr::mutate(
      .pool_rule = dplyr::case_when(
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("Enhancing diagnostic accuracy in rare and common fundus diseases with a knowledge-rich vision-language model"))
        ) ~ "Wang_RetiZero",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("Early diagnosis of axial spondyloarthritis in primary care using multi-agent systems"))
        ) ~ "Ji",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("Towards accurate differential diagnosis with large language models"))
        ) ~ "McDuff",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("Grounding large language models in clinical diagnostics"))
        ) ~ "Chen",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("The virtual reference radiologist: comprehensive AI assistance for clinical image reading and interpretation"))
        ) ~ "Siepmann",
        # Shared-case reader-stratification rules. Matched on apostrophe-free,
        # unique title substrings to avoid the Wang/Xu/Zhang author collisions and
        # any apostrophe-encoding pitfalls.
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("ultrasound-based deep learning model as an assistant improves the diagnosis of ovarian tumors"))
        ) ~ "Wang_ovarian",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("predicting invasiveness of lung adenocarcinoma from chest ct with few-shot vision-language ternary"))
        ) ~ "Xu_pGGN",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("cystic renal masses to senior-level accuracy"))
        ) ~ "Jin",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("diagnosing focal liver lesions from ct/mri reports"))
        ) ~ "Sheng",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("hyperparameter and input optimization in gpt-5"))
        ) ~ "Camur",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("gpt-4o-driven automated feature recognition and validation in clinical settings"))
        ) ~ "Zhang_ovarianCT",
        # Rules that reuse the distinct-strata (McDuff) and reader-strata
        # (Wang/RetiZero) modes:
        #   Liu      -> distinct_strata_pool  (same harmonised field "general medicine";
        #               disjoint pulmonology vs endocrinology case sets, different
        #               physicians -> independent-subgroup combination).
        #   Schramm  -> shared_cases_reader_pool (neurology vs radiology resident groups
        #               read the SAME 40 brain-MRI cases; pooled within Type_of_Results).
        #   Ver Berne-> shared_cases_reader_pool (general dentists vs DMF radiologists read
        #               the SAME 25 jawbone cases).
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("a generalist medical language model for disease diagnosis assistance"))
        ) ~ "Liu",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("reader experience shapes accuracy gains"))
        ) ~ "Schramm",
        stringr::str_detect(
          title_norm,
          stringr::fixed(.pool_norm_title("radiological diagnosis of jawbone lesions"))
        ) ~ "Ver_Berne",
        TRUE ~ "none"
      )
    )
}

# Preserve an untouched granular data object for direct sensitivity comparison.
df_performance_granular <- df_performance %>%
  dplyr::mutate(
    .raw_order = dplyr::row_number(),
    # IDs are identifiers, not quantitative variables. Character type also permits
    # auditable pooled IDs without bind_rows() type conflicts.
    Experiment_ID = as.character(Experiment_ID),
    Effect_ID = as.character(Effect_ID),
    ES_ID = as.character(ES_ID)
  ) %>%
  .assign_pool_rule() %>%
  dplyr::mutate(
    .study_pool_id = dplyr::coalesce(
      as.character(.data$ID),
      as.character(.data$Title)
    )
  )

# -----------------------------------------------------------------------------
# Candidate checks
# -----------------------------------------------------------------------------
# Use the same whitespace normalization later applied to moderators in the
# meta-regression, but only inside the audit copy (the granular source stays intact).
candidate_check_vars <- unique(c(mods_pool_check, pool_guard_vars))
df_aggregation_candidate_check <- df_performance_granular %>%
  # Exclude empty/placeholder CSV rows from the audit only; the analysis source
  # object itself is left unchanged.
  dplyr::filter(!is.na(ES_ID), !is.na(Title), !is.na(.study_pool_id)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(candidate_check_vars),
      ~ if (is.character(.x) || is.factor(.x)) norm(.x) else .x
    )
  )

# Check 1 follows the proposed rule literally: within a paper, which sets of
# rows have identical values across ALL moderators used in the analyses?
strict_candidate_cols <- c(".study_pool_id", mods_pool_check)

aggregation_candidates_same_moderators <- df_aggregation_candidate_check %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(strict_candidate_cols))) %>%
  dplyr::summarise(
    Author = dplyr::first(Author),
    Title = dplyr::first(Title),
    ID = dplyr::first(ID),
    n_rows = dplyr::n(),
    source_ES_IDs = paste(ES_ID, collapse = " | "),
    covered_by_explicit_rule = any(.pool_rule != "none"),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_rows > 1) %>%
  dplyr::arrange(dplyr::desc(n_rows), Author)

readr::write_csv(
  aggregation_candidates_same_moderators,
  file.path(output_dir, "aggregation_candidates_same_moderators.csv")
)

# Check 2 is deliberately more conservative. It additionally requires equality
# on non-moderator guard variables that can distinguish separate measurements.
safe_candidate_cols <- unique(c(".study_pool_id", mods_pool_check, pool_guard_vars))

aggregation_candidates_safe <- df_aggregation_candidate_check %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(safe_candidate_cols))) %>%
  dplyr::summarise(
    Author = dplyr::first(Author),
    Title = dplyr::first(Title),
    ID = dplyr::first(ID),
    n_rows = dplyr::n(),
    source_ES_IDs = paste(ES_ID, collapse = " | "),
    covered_by_explicit_rule = any(.pool_rule != "none"),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_rows > 1) %>%
  dplyr::arrange(dplyr::desc(n_rows), Author)

readr::write_csv(
  aggregation_candidates_safe,
  file.path(output_dir, "aggregation_candidates_safe_same_moderators.csv")
)

# This is the high-value review file: rows that satisfy the conservative check
# but are NOT part of the study-specific pooling rules above. They are only
# flagged; they are never pooled automatically.
aggregation_candidates_additional <- aggregation_candidates_safe %>%
  dplyr::filter(!covered_by_explicit_rule)

readr::write_csv(
  aggregation_candidates_additional,
  file.path(output_dir, "aggregation_candidates_additional_review.csv")
)

message(
  "[aggregation] Candidate check: ",
  nrow(aggregation_candidates_same_moderators),
  " within-paper groups share all moderators; ",
  nrow(aggregation_candidates_additional),
  " conservative candidate groups are NOT covered by the five explicit rules."
)
if (nrow(aggregation_candidates_additional) > 0) {
  message("[aggregation] Additional pooling candidates require manual review; they are NOT changed automatically:")
  print(as.data.frame(
    aggregation_candidates_additional %>%
      dplyr::select(dplyr::any_of(c("Author", "Title", "ID", "n_rows", "source_ES_IDs")))
  ))
}

# -----------------------------------------------------------------------------
# Explicit pooling keys
# -----------------------------------------------------------------------------
# These rules encode only the study-specific pooling decisions used in  the
# harmonized (primary) analysis. All other rows receive a unique KEEP key.
df_performance_granular <- df_performance_granular %>%
  dplyr::mutate(
    .pool_mode = dplyr::case_when(
      .pool_rule %in% c("Wang_RetiZero", "Ji", "Siepmann",
                        # same mechanism throughout: distinct readers evaluate
                        # the SAME cases in the SAME condition.
                        "Wang_ovarian", "Xu_pGGN", "Jin", "Sheng",
                        "Camur", "Zhang_ovarianCT",
                        # Schramm/Ver Berne: distinct reader GROUPS on shared cases.
                        "Schramm", "Ver_Berne") ~ "shared_cases_reader_pool",
      # Liu: independent disjoint-case specialty strata, same harmonised field.
      .pool_rule %in% c("McDuff", "Liu") ~ "distinct_strata_pool",
      .pool_rule == "Chen" ~ "disjoint_case_subsets_same_physicians",
      TRUE ~ "keep"
    ),
    .pool_group = dplyr::case_when(
      .pool_rule == "Wang_RetiZero" ~ paste0("Wang_RetiZero|", norm(Career_stage)),
      .pool_rule == "Ji" ~ paste0("Ji|", norm(Career_stage), "|", norm(Time_lag_sessions)),
      .pool_rule == "McDuff" ~ paste0("McDuff|", norm(Medical_Field), "|", norm(Type_of_Results)),
      .pool_rule == "Chen" ~ "Chen|overall",
      .pool_rule == "Siepmann" ~ paste0("Siepmann|", norm(Type_of_Results)),
      # Effect-level pool keys (one pooled effect per key). Content columns
      # only -> robust to the ID_2/Exp_ID/ES_ID renumbering.
      .pool_rule == "Wang_ovarian" ~ paste0("Wang_ovarian|", norm(Career_stage)),
      .pool_rule == "Xu_pGGN" ~ paste0("Xu_pGGN|", norm(Career_stage)),
      .pool_rule == "Jin" ~ paste0("Jin|", norm(Career_stage)),
      .pool_rule == "Sheng" ~ paste0("Sheng|", norm(Career_stage), "|", norm(Type_of_Results)),
      # Camur: browser (13-month washout) vs optimal-settings (16-month washout) are
      # distinguished by Time_lag_sessions; pooling R1+R2 WITHIN each condition.
      .pool_rule == "Camur" ~ paste0("Camur|", norm(Time_lag_sessions)),
      # Zhang: key on career x case-cohort so ONLY the two attending radiologists
      # (230 cases) pool; G4 (attending, 224) and G5 (resident, 224) fall into
      # singleton groups and pass through unchanged.
      .pool_rule == "Zhang_ovarianCT" ~ paste0("Zhang_ovarianCT|", norm(Career_stage), "|", norm(N_cases_total)),
      # Liu: pool the two disjoint specialty strata WITHIN each career level.
      .pool_rule == "Liu" ~ paste0("Liu|", norm(Career_stage)),
      # Schramm: 6 rows over TWO career stages. Pool the two resident reader groups
      # (neuro + radiology) within each metric; the attending neuroradiologist group
      # is a singleton. Key on career x outcome (Sheng pattern) so Career_stage is
      # never collapsed. -> 6 rows become 4 (2 pooled resident + 2 attending singletons).
      .pool_rule == "Schramm" ~ paste0("Schramm|", norm(Career_stage), "|", norm(Type_of_Results)),
      # Ver Berne: single top-1 outcome -> one pooled group across both reader groups.
      .pool_rule == "Ver_Berne" ~ paste0("Ver_Berne|", norm(Type_of_Results)),
      TRUE ~ paste0("KEEP|", ES_ID)
    ),
    # Experiment groups are slightly broader than effect groups when multiple
    # accuracy outcomes arise from the same pooled experimental condition.
    .pool_experiment_group = dplyr::case_when(
      .pool_rule == "Wang_RetiZero" ~ paste0("Wang_RetiZero|", norm(Career_stage)),
      .pool_rule == "Ji" ~ paste0("Ji|", norm(Career_stage), "|", norm(Time_lag_sessions)),
      .pool_rule == "McDuff" ~ paste0("McDuff|", norm(Medical_Field)),
      .pool_rule == "Chen" ~ "Chen|overall",
      .pool_rule == "Siepmann" ~ "Siepmann|overall",
      # Experiment-level grouping. Sheng and Camur mirror Siepmann: several
      # accuracy outcomes / conditions from the SAME pooled experimental condition
      # sit as separate effects inside ONE experiment.
      .pool_rule == "Wang_ovarian" ~ paste0("Wang_ovarian|", norm(Career_stage)),
      .pool_rule == "Xu_pGGN" ~ paste0("Xu_pGGN|", norm(Career_stage)),
      .pool_rule == "Jin" ~ paste0("Jin|", norm(Career_stage)),
      .pool_rule == "Sheng" ~ paste0("Sheng|", norm(Career_stage)),
      .pool_rule == "Camur" ~ "Camur|overall",
      .pool_rule == "Zhang_ovarianCT" ~ paste0("Zhang_ovarianCT|", norm(Career_stage), "|", norm(N_cases_total)),
      # Liu: one experiment per career level. Schramm/Ver Berne mirror Siepmann
      # (the accuracy metrics are separate effects inside ONE experiment).
      .pool_rule == "Liu" ~ paste0("Liu|", norm(Career_stage)),
      .pool_rule == "Schramm" ~ paste0("Schramm|", norm(Career_stage)),
      .pool_rule == "Ver_Berne" ~ "Ver_Berne|overall",
      TRUE ~ paste0("KEEP|", Experiment_ID)
    )
  )

# Give every pooled experimental condition a stable Experiment_ID by reusing
# the smallest source Experiment_ID within that condition. This changes no
# study-level ID and creates no new cross-paper cluster.
pool_experiment_id_map <- df_performance_granular %>%
  dplyr::filter(.pool_rule != "none") %>%
  dplyr::group_by(.pool_experiment_group) %>%
  dplyr::summarise(
    # Reuse the first source Experiment_ID as a stable CHARACTER identifier.
    # No numerical operation on IDs is needed for the meta-analytic model.
    .pooled_Experiment_ID = {
      ord <- order(.raw_order)
      z <- as.character(Experiment_ID[ord])
      z <- z[!is.na(z) & z != ""]
      if (length(z) == 0) NA_character_ else z[[1]]
    },
    .groups = "drop"
  )

df_performance_granular <- df_performance_granular %>%
  dplyr::left_join(pool_experiment_id_map, by = ".pool_experiment_group")

# -----------------------------------------------------------------------------
# Pooling helpers
# -----------------------------------------------------------------------------
.pool_sum_available <- function(x) {
  z <- .pool_as_num(x)
  if (all(is.na(z))) return(NA_real_)
  sum(z, na.rm = TRUE)
}

.pool_sum_complete <- function(x) {
  z <- .pool_as_num(x)
  if (length(z) == 0 || any(is.na(z))) return(NA_real_)
  sum(z)
}

.pool_common_num <- function(x) {
  z <- unique(.pool_as_num(x))
  z <- z[!is.na(z)]
  if (length(z) == 1) z else NA_real_
}

.pool_common_chr <- function(x) {
  z <- unique(as.character(x))
  z <- z[!is.na(z) & z != ""]
  if (length(z) == 1) z else NA_character_
}

.pool_weighted_mean <- function(x, w) {
  xx <- .pool_as_num(x)
  ww <- .pool_as_num(w)
  ok <- is.finite(xx) & is.finite(ww) & ww > 0
  if (!any(ok)) {
    ok_x <- is.finite(xx)
    if (!any(ok_x)) return(NA_real_)
    return(mean(xx[ok_x]))
  }
  stats::weighted.mean(xx[ok], ww[ok])
}

.pool_clear_numeric <- function(df, cols) {
  for (nm in intersect(cols, names(df))) df[[nm]] <- NA_real_
  df
}

.pool_bernoulli_sd_same_scale <- function(m) {
  m <- .pool_as_num(m)
  if (!is.finite(m)) return(NA_real_)
  scale <- if (m > 1) 100 else 1
  p <- m / scale
  if (p < 0 || p > 1) return(NA_real_)
  sqrt(p * (1 - p)) * scale
}

.pool_one_group <- function(g) {
  # No transformation is made when a rule leaves a singleton group.
  if (nrow(g) == 1 || identical(g$.pool_mode[[1]], "keep")) {
    return(g)
  }

  mode <- g$.pool_mode[[1]]
  rule <- g$.pool_rule[[1]]
  source_ids <- paste(g$ES_ID, collapse = " | ")

  # Start from the first row so all non-analytic study descriptors remain intact.
  out <- g[1, , drop = FALSE]
  out$.raw_order <- min(g$.raw_order, na.rm = TRUE)

  # Arm means are pooled on their own observation denominators. This is exact
  # for the binary accuracy data used in the study-specific pooling rules.
  out$M_control   <- .pool_weighted_mean(g$M_control,   g$n_control)
  out$M_treatment <- .pool_weighted_mean(g$M_treatment, g$n_treatment)

  # Raw analysis denominators are additive across source rows in all study-specific rules.
  if ("n_control" %in% names(out))   out$n_control   <- .pool_sum_available(g$n_control)
  if ("n_treatment" %in% names(out)) out$n_treatment <- .pool_sum_available(g$n_treatment)

  # Case counts and physician counts depend on WHY rows are being pooled.
  if (mode == "shared_cases_reader_pool") {
    # Same cases are repeatedly evaluated by distinct readers/reader strata.
    if ("n_control_case" %in% names(out))   out$n_control_case   <- .pool_common_num(g$n_control_case)
    if ("n_treatment_case" %in% names(out)) out$n_treatment_case <- .pool_common_num(g$n_treatment_case)
    if ("N_cases_total" %in% names(out))    out$N_cases_total    <- .pool_common_num(g$N_cases_total)

    if ("n_control_physician" %in% names(out))
      out$n_control_physician <- .pool_sum_available(g$n_control_physician)
    if ("n_treatment_physician" %in% names(out))
      out$n_treatment_physician <- .pool_sum_available(g$n_treatment_physician)

  } else if (mode == "distinct_strata_pool") {
    # McDuff: pooled strata represent different specialty case/reader strata
    # that map to the same harmonized Medical_Field.
    if ("n_control_case" %in% names(out))   out$n_control_case   <- .pool_sum_available(g$n_control_case)
    if ("n_treatment_case" %in% names(out)) out$n_treatment_case <- .pool_sum_available(g$n_treatment_case)
    if ("N_cases_total" %in% names(out))    out$N_cases_total    <- .pool_sum_available(g$n_control_case)

    if ("n_control_physician" %in% names(out))
      out$n_control_physician <- .pool_sum_available(g$n_control_physician)
    if ("n_treatment_physician" %in% names(out))
      out$n_treatment_physician <- .pool_sum_available(g$n_treatment_physician)

  } else if (mode == "disjoint_case_subsets_same_physicians") {
    # Chen: the 20 + 20 + 20 case subsets are disjoint, while the same physician
    # groups are used across subsets. Cases add; physician counts do not.
    if ("n_control_case" %in% names(out))   out$n_control_case   <- .pool_sum_available(g$n_control_case)
    if ("n_treatment_case" %in% names(out)) out$n_treatment_case <- .pool_sum_available(g$n_treatment_case)
    if ("N_cases_total" %in% names(out))    out$N_cases_total    <- .pool_sum_available(g$n_control_case)

    if ("n_control_physician" %in% names(out))
      out$n_control_physician <- .pool_common_num(g$n_control_physician)
    if ("n_treatment_physician" %in% names(out))
      out$n_treatment_physician <- .pool_common_num(g$n_treatment_physician)
  }

  # h_n_pairs is additive only when it was actually defined in every source row.
  # For within-subject rows with missing h_n_pairs, fall back to the pooled arm N
  # only when both arm Ns are complete and equal.
  if ("h_n_pairs" %in% names(out)) {
    hp <- .pool_sum_complete(g$h_n_pairs)
    is_within <- any(norm(g$h_design) == "within_subject", na.rm = TRUE)
    nc <- if ("n_control" %in% names(out)) .pool_as_num(out$n_control) else NA_real_
    nt <- if ("n_treatment" %in% names(out)) .pool_as_num(out$n_treatment) else NA_real_
    if (!is.finite(hp) && is_within && is.finite(nc) && is.finite(nt) && nc == nt) hp <- nc
    out$h_n_pairs <- hp
  }

  # Discordant cells are counts and can be added across the explicitly pooled
  # reader strata, but ONLY if every contributing row supplies the cell.
  if ("h_b" %in% names(out)) out$h_b <- .pool_sum_complete(g$h_b)
  if ("h_c" %in% names(out)) out$h_c <- .pool_sum_complete(g$h_c)
  if ("Discordant_b" %in% names(out)) out$Discordant_b <- .pool_sum_complete(g$Discordant_b)
  if ("Discordant_c" %in% names(out)) out$Discordant_c <- .pool_sum_complete(g$Discordant_c)

  # A reader-specific r is not itself additive. After pooling, r is therefore
  # re-derived from pooled b/c when available; otherwise the existing donor-r
  # procedure is used later in the script.
  if ("h_r_used" %in% names(out)) out$h_r_used <- NA_real_
  if ("h_r_source" %in% names(out)) {
    if (is.finite(.pool_as_num(out$h_b)) && is.finite(.pool_as_num(out$h_c))) {
      out$h_r_source <- "pooled discordant cells (aggregation step)"
    } else {
      out$h_r_source <- NA_character_
    }
  }

  # Reconstruct study-level N fields from the pooled arm-specific counts.
  nc_phys <- if ("n_control_physician" %in% names(out)) .pool_as_num(out$n_control_physician) else NA_real_
  nt_phys <- if ("n_treatment_physician" %in% names(out)) .pool_as_num(out$n_treatment_physician) else NA_real_
  is_parallel <- any(norm(g$h_design) == "parallel", na.rm = TRUE) &&
    !any(norm(g$h_design) == "within_subject", na.rm = TRUE)

  if ("N_physicians_total" %in% names(out)) {
    out$N_physicians_total <- if (is_parallel && is.finite(nc_phys) && is.finite(nt_phys)) {
      nc_phys + nt_phys
    } else {
      max(c(nc_phys, nt_phys), na.rm = TRUE)
    }
    if (!is.finite(out$N_physicians_total)) out$N_physicians_total <- NA_real_
  }

  nc <- if ("n_control" %in% names(out)) .pool_as_num(out$n_control) else NA_real_
  nt <- if ("n_treatment" %in% names(out)) .pool_as_num(out$n_treatment) else NA_real_
  if ("N_assessments_total" %in% names(out)) {
    out$N_assessments_total <- if (is.finite(nc) && is.finite(nt)) nc + nt else NA_real_
  }

  # N_total / n_total are retained as the per-comparison analysis-unit count,
  # not as the sum across the two arms. This mirrors the source coding in the
  # paired/shared-case studies and gives Chen 60 rather than 120.
  n_comp <- if (is.finite(nc) && is.finite(nt) && nc == nt) nc else {
    if (is.finite(nc) && is.finite(nt)) nc + nt else dplyr::coalesce(nc, nt)
  }
  if ("N_total" %in% names(out)) out$N_total <- n_comp
  if ("n_total" %in% names(out)) out$n_total <- n_comp

  # The pooling rules concern diagnostic accuracy. Once the arm means have been
  # recomputed, source-row SD/SE/CI values are no longer valid for the pooled row.
  # For binary/proportion-like accuracy, reconstruct the Bernoulli raw-score SD on
  # the SAME scale as the pooled mean, so the downstream metafor code standardizes
  # the pooled effect exactly as it does for non-pooled rows.
  measure_class <- unique(norm(as.character(g$h_measure_class)))
  measure_class <- measure_class[!is.na(measure_class) & nzchar(measure_class)]
  if (length(measure_class) != 1)
    stop("[aggregation] Pooling rule '", rule, "' mixes h_measure_class values (",
         paste(measure_class, collapse = " / "), ").")
  is_proportion_like <- identical(measure_class, "proportion")

  if (is_proportion_like) {
    if ("SD_control" %in% names(out))   out$SD_control   <- .pool_bernoulli_sd_same_scale(out$M_control)
    if ("SD_treatment" %in% names(out)) out$SD_treatment <- .pool_bernoulli_sd_same_scale(out$M_treatment)
  } else {
    if ("SD_control" %in% names(out))   out$SD_control   <- NA_real_
    if ("SD_treatment" %in% names(out)) out$SD_treatment <- NA_real_
  }

  out <- .pool_clear_numeric(
    out,
    c(
      "SE_control", "95CI_lower_control", "95CI_upper_control",
      "SE_treatment", "95CI_lower_treatment", "95CI_upper_treatment",
      "Performance_Difference_SD", "Performance_Difference_SE",
      "95CI_lower_performance_diff", "95CI_upper_performance_diff",
      "Performance_difference_p_value", "Odds_Ratio",
      "95CI_lower_odds_ratio", "95CI_upper_odds_ratio",
      "F_value", "Chi_square", "unstandardised_beta", "SE_beta",
      "Z_value", "t_value"
    )
  )

  if ("Performance_Difference" %in% names(out) &&
      is.finite(.pool_as_num(out$M_control)) && is.finite(.pool_as_num(out$M_treatment))) {
    out$Performance_Difference <- .pool_as_num(out$M_treatment) - .pool_as_num(out$M_control)
  }

  # Preserve study-level variance assumptions only when they are common across
  # all source rows. Otherwise force later code to make the variance route
  # explicit rather than inheriting an arbitrary first-row value.
  if ("h_icc" %in% names(out)) out$h_icc <- .pool_common_num(g$h_icc)
  if ("variance_source" %in% names(out)) out$variance_source <- .pool_common_chr(g$variance_source)
  if ("Reader_cluster_adjusted_SD_SE_CI" %in% names(out))
    out$Reader_cluster_adjusted_SD_SE_CI <- .pool_common_chr(g$Reader_cluster_adjusted_SD_SE_CI)

  # The pooled row's unit of analysis depends on
  # WHAT the pooling did, not on the mode name:
  #   * shared_cases_reader_pool STACKS distinct readers onto the SAME cases. The
  #     pooled unit becomes physician-case even if each single-reader source row sat
  #     at "case" level, so we FORCE "physician-case" -> the reader-clustering design
  #     effect fires (.compute_de keys on low(h_var_level) == "physician-case").
  #     n_pairs, b/c and the physician count were aggregated to exactly this level
  #     above, so m_bar = n_pairs / n_physician is mutually consistent.
  #   * distinct_strata_pool (McDuff, Liu) combines INDEPENDENT disjoint-case strata,
  #     and disjoint_case_subsets_same_physicians (Chen) concatenates disjoint case
  #     batches. Neither stacks readers onto shared cases, so the faithful unit is the
  #     SOURCE level, which is inherited (case -> DE = 1 for the reader-aggregated
  #     McDuff/Chen rows; physician-case -> DE fires with the correct m_bar for the
  #     per-reader-resolved Liu rows). Inheriting the source level avoids imposing a
  #     spurious reader design effect on case-level source data and keeps the McDuff
  #     physician-count sum irrelevant to the variance.
  if (identical(mode, "shared_cases_reader_pool")) {
    if ("h_var_level" %in% names(out))  out$h_var_level  <- "physician-case"
    if ("h_unit_level" %in% names(out)) out$h_unit_level <- "physician-case"
  } else {
    if ("h_var_level" %in% names(out))  out$h_var_level  <- .pool_common_chr(g$h_var_level)
    if ("h_unit_level" %in% names(out)) out$h_unit_level <- .pool_common_chr(g$h_unit_level)
  }

  # New auditable IDs. The study cluster (ID_2) is intentionally unchanged.
  pooled_exp_id <- as.character(g$.pooled_Experiment_ID[[1]])
  if (!is.na(pooled_exp_id) && pooled_exp_id != "") {
    out$Experiment_ID <- pooled_exp_id
  }
  pooled_es_id <- paste0("POOL_", rule, "_", sprintf("%03d", min(g$.raw_order, na.rm = TRUE)))
  out$ES_ID <- pooled_es_id
  out$Effect_ID <- pooled_es_id

  if ("Note" %in% names(out)) {
    old_note <- as.character(out$Note)
    if (is.na(old_note)) old_note <- ""
    agg_note <- paste0(
      "[aggregation] Harmonized pooling applied in R; rule=", rule,
      "; source rows=", nrow(g), "; source ES_IDs=", source_ids,
      ". Raw coding remains unchanged in the source file."
    )
    out$Note <- stringr::str_squish(paste(old_note, agg_note))
  }

  out
}

.apply_harmonized_pooling <- function(df) {
  keep_rows <- df %>% dplyr::filter(.pool_rule == "none")
  pool_rows <- df %>% dplyr::filter(.pool_rule != "none")

  if (nrow(pool_rows) == 0) return(df)

  pooled_rows <- pool_rows %>%
    dplyr::group_split(.pool_group, .keep = TRUE) %>%
    purrr::map_dfr(.pool_one_group)

  dplyr::bind_rows(keep_rows, pooled_rows) %>%
    dplyr::arrange(.raw_order)
}

# -----------------------------------------------------------------------------
# PRE-POOLING HOMOGENEITY GUARD
# -----------------------------------------------------------------------------
# Every set of rows collapsed into ONE pooled row MUST be constant on all analysis
# moderators. Otherwise the pooled row silently inherits the first source row's
# moderator value (out <- g[1, ]) and any later meta-regression on that moderator
# would run on a mislabeled aggregate. This guard runs BEFORE pooling and hard-stops
# on any within-group moderator conflict.
#
# Severity per (pool group x variable):
#   "conflict"   : >1 distinct non-NA value in the group  -> pooling would mislabel.
#   "partial_NA" : exactly one non-NA value, but some rows NA -> reported, not fatal.
# Only a moderator "conflict" aborts the run. Guard variables are additionally
# reported (for review) but do not, by themselves, stop the pipeline.
.check_pool_group_homogeneity <- function(df, mod_vars, guard_vars = character(0),
                                          stop_on_violation = TRUE) {
  pool_only <- df %>% dplyr::filter(.pool_rule != "none")
  if (nrow(pool_only) == 0) {
    message("[homogeneity guard] Homogeneity guard: no rows flagged for pooling; nothing to check.")
    return(invisible(NULL))
  }

  chk_vars <- intersect(unique(c(mod_vars, guard_vars)), names(pool_only))

  viol <- purrr::map_dfr(chk_vars, function(v) {
    pool_only %>%
      dplyr::group_by(.pool_rule, .pool_group) %>%
      dplyr::summarise(
        variable      = v,
        n_rows        = dplyr::n(),
        distinct_vals = dplyr::n_distinct(norm(.data[[v]][!is.na(.data[[v]])])),
        any_na        = any(is.na(.data[[v]])),
        values        = paste(sort(unique(norm(.data[[v]][!is.na(.data[[v]])]))),
                              collapse = " | "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        is_moderator = v %in% mod_vars,
        severity = dplyr::case_when(
          distinct_vals > 1           ~ "conflict",
          distinct_vals == 1 & any_na ~ "partial_NA",
          TRUE                        ~ "ok"
        )
      ) %>%
      dplyr::filter(severity != "ok")
  })

  if (nrow(viol) == 0) {
    message("[homogeneity guard] Homogeneity guard passed: all pooled groups are constant on ",
            length(chk_vars), " checked variables (", length(mod_vars),
            " moderators + ", length(guard_vars), " guard variables).")
    return(invisible(NULL))
  }

  viol <- viol %>%
    dplyr::arrange(dplyr::desc(is_moderator), severity, .pool_group, variable)
  readr::write_csv(viol, file.path(output_dir, "aggregation_homogeneity_violations.csv"))

  message("[homogeneity guard] Homogeneity guard found issues (written to aggregation_homogeneity_violations.csv):")
  print(as.data.frame(viol))

  conflicts <- viol %>% dplyr::filter(severity == "conflict", is_moderator)
  if (stop_on_violation && nrow(conflicts) > 0) {
    stop(
      "[homogeneity guard] Pooling aborted: ", nrow(conflicts),
      " moderator/pool-group combination(s) are NOT constant within the group. ",
      "Collapsing those rows would mislabel the aggregate on a moderator you test. ",
      "Fix the coding or refine the .pool_group key, then rerun. ",
      "(partial_NA rows are reported but do not stop the run.)"
    )
  }
  invisible(viol)
}

# Guard is only meaningful in harmonized mode (granular mode pools nothing).
if (AGGREGATION_MODE == "harmonized") {
  .check_pool_group_homogeneity(
    df_performance_granular,
    mod_vars   = mods_pool_check,   # the moderators we test repeatedly -> hard stop
    guard_vars = pool_guard_vars,   # reported for review only
    stop_on_violation = TRUE
  )
}

# Apply only when requested. In granular mode the original observations are
# passed to the remainder of the script unchanged.
if (AGGREGATION_MODE == "harmonized") {
  df_performance <- .apply_harmonized_pooling(df_performance_granular)
} else {
  df_performance <- df_performance_granular
}

# -----------------------------------------------------------------------------
# Audit the five explicit decisions and verify expected row counts
# -----------------------------------------------------------------------------
expected_pool_counts <- tibble::tibble(
  .pool_rule = c("Wang_RetiZero", "Ji", "McDuff", "Chen", "Siepmann"),
  expected_rows_after = c(2L, 4L, 4L, 1L, 2L)
)

aggregation_applied_audit <- df_performance_granular %>%
  dplyr::filter(.pool_rule != "none") %>%
  dplyr::group_by(.pool_rule) %>%
  dplyr::summarise(
    rows_before = dplyr::n(),
    groups_defined = dplyr::n_distinct(.pool_group),
    source_ES_IDs = paste(ES_ID, collapse = " | "),
    .groups = "drop"
  ) %>%
  dplyr::left_join(expected_pool_counts, by = ".pool_rule") %>%
  dplyr::mutate(
    expected_structure_ok = groups_defined == expected_rows_after,
    aggregation_mode = AGGREGATION_MODE
  )

readr::write_csv(
  aggregation_applied_audit,
  file.path(output_dir, "aggregation_applied_audit.csv")
)

message("[aggregation] Explicit aggregation audit:")
print(as.data.frame(aggregation_applied_audit))

bad_expected <- aggregation_applied_audit %>%
  dplyr::filter(!expected_structure_ok)
if (nrow(bad_expected) > 0) {
  warning(
    "[aggregation] At least one explicit study no longer matches the expected pooling structure. ",
    "Review aggregation_applied_audit.csv before interpreting results."
  )
}

message(
  "[aggregation] Analysis mode = ", AGGREGATION_MODE,
  "; rows entering the effect-size pipeline: ", nrow(df_performance),
  " (granular source rows: ", nrow(df_performance_granular), ")."
)

# Remove internal aggregation-helper columns before the downstream analysis code.
df_performance <- df_performance %>%
  dplyr::select(-dplyr::any_of(c(
    ".raw_order", ".study_pool_id", ".pool_rule", ".pool_mode",
    ".pool_group", ".pool_experiment_group", ".pooled_Experiment_ID"
  )))


# Quick structure checks
dplyr::glimpse(df_performance[, c(
  "M_treatment","SD_treatment","SE_treatment",
  "M_control","SD_control","SE_control",
  "n_treatment","n_control","Odds_Ratio",
  "Performance_Difference", "95CI_lower_performance_diff", "95CI_upper_performance_diff",
  "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment"
)])

print(names(df_performance))

sapply(
  df_performance[, c(
    "M_treatment","SD_treatment","SE_treatment",
    "M_control","SD_control","SE_control",
    "n_treatment","n_control","Odds_Ratio",
    "Performance_Difference", "95CI_lower_performance_diff", "95CI_upper_performance_diff",
    "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment"
  )],
  function(x) sum(!is.na(x))
)

# =============================================================================
# Effect-size and sampling-variance computation (metafor-based)
# -----------------------------------------------------------------------------
# Effect sizes and their sampling variances are computed with metafor's own
# effect-size machinery wherever the information reported by the primary study
# permits it.
#
# PRINCIPLE:
#   1) Point estimate / basic sampling variance are computed by metafor.
#   2) Paired designs use metafor::escalc(measure = "SMCRPH"), so r enters the
#      SAMPLING VARIANCE while the standardized point estimate is based on the
#      average raw-score variance. This preserves the intended d_av logic.
#   3) Binary (i.e., h_measure_class == "proportion") propoaccuracy is represented as a 0/1 outcome. Hence p is the mean and
#      sqrt(p*(1-p)) is its raw-score SD. This SD standardizes the point estimate;
#      r / ICC / design effects concern the precision of that estimate.
#   4) The reader-clustering design effect is applied AFTER metafor has produced
#      the row-level sampling variance. It never changes the point estimate.
#
# IMPORTANT: metafor::escalc() already applies the relevant small-sample bias
# correction for SMD/SMDH/SMCRPH by default. Therefore there is NO separate
# manual Hedges-J step in this module.
# =============================================================================

# ---- small helpers ----------------------------------------------------------
ensure_col <- function(df, nm, default = NA) {
  if (!nm %in% names(df)) df[[nm]] <- default
  df
}
low     <- function(x) stringr::str_squish(stringr::str_to_lower(as.character(x)))
to_prop <- function(x) dplyr::if_else(!is.na(x) & x > 1, x / 100, x)
as_num  <- function(x) suppressWarnings(as.numeric(x))

# Make sure every referenced column exists (defensive; only adds columns if absent)
for (nm in c(
  "h_design","h_unit_level","h_var_level","h_n_pairs","h_r_used",
  "h_r_source","h_b","h_c","h_measure_class",
  "Parallel_or_Crossover","Same_Cases_Both_Conditions",
  "Unit_of_Analysis_for_ES",
  "Discordant_b","Discordant_c",
  "n_control_case","n_control_physician","n_treatment_case",
  "N_cases_total","N_physicians_total","n_treatment_physician",
  "95CI_lower_odds_ratio","95CI_upper_odds_ratio",
  "Reader_cluster_adjusted_SD_SE_CI","variance_source","h_icc",
  "cohens_d"
)) df_performance <- ensure_col(df_performance, nm)

# Preserve any Cohen's d that came from the coding sheet. Later sensitivity
# refits call safe_calc_effects() on an already processed data frame, so we must
# not accidentally treat our own calculated `cohens_d` as a reported input.
if (!".cohens_d_reported" %in% names(df_performance)) {
  df_performance$.cohens_d_reported <- as_num(df_performance$cohens_d)
}

# ---- proportion normalization audit ----------------------------------------
prop_cols <- intersect(c("M_control", "M_treatment"), names(df_performance))
norm_audit <- df_performance %>%
  dplyr::mutate(.row = dplyr::row_number()) %>%
  dplyr::filter(dplyr::if_any(dplyr::all_of(prop_cols),
                              ~ !is.na(as_num(.x)) & as_num(.x) > 1)) %>%
  dplyr::select(.row, dplyr::any_of(c("Author", "ES_ID", prop_cols)))
if (nrow(norm_audit) > 0) {
  message("[effect sizes] Values > 1 in proportion columns; interpreted as percentages and divided by 100 for proportion rows:")
  print(as.data.frame(norm_audit))
}

# ---- harmonized design / measure / pairing fields --------------------------
stopifnot(all(df_performance$h_design %in% c("within_subject", "parallel")))
stopifnot(all(df_performance$h_measure_class %in% c("proportion", "continuous")))

.prepare_es_inputs <- function(df) {
  df %>%
    dplyr::mutate(
      .design = as.character(h_design),
      .paired = dplyr::case_when(
        .design == "within_subject" &
          (is.na(Same_Cases_Both_Conditions) |
             as_num(Same_Cases_Both_Conditions) == 1) ~ TRUE,
        TRUE ~ FALSE
      ),
      .measure = as.character(h_measure_class),
      .n_pairs = dplyr::case_when(
        !is.na(as_num(h_n_pairs))        ~ as_num(h_n_pairs),
        low(h_unit_level) == "physician" ~ as_num(n_control_physician),
        TRUE                              ~ as_num(n_control_case)
      ),
      .p_c = dplyr::if_else(.measure == "proportion",
                            to_prop(as_num(M_control)), NA_real_),
      .p_t = dplyr::if_else(.measure == "proportion",
                            to_prop(as_num(M_treatment)), NA_real_),
      .b = dplyr::coalesce(as_num(h_b), as_num(Discordant_b)),
      .c = dplyr::coalesce(as_num(h_c), as_num(Discordant_c)),
      .r_precomputed     = as_num(h_r_used),
      .r_source_text     = low(dplyr::coalesce(as.character(h_r_source), "")),

      # Binary 0/1 outcome SDs. These standardize the POINT ESTIMATE.
      .sd_c_bern = dplyr::if_else(.measure == "proportion" &
                                    !is.na(.p_c) & .p_c >= 0 & .p_c <= 1,
                                  sqrt(.p_c * (1 - .p_c)), NA_real_),
      .sd_t_bern = dplyr::if_else(.measure == "proportion" &
                                    !is.na(.p_t) & .p_t >= 0 & .p_t <= 1,
                                  sqrt(.p_t * (1 - .p_t)), NA_real_),

      # If a continuous arm SD is absent, reconstruct it from a reported arm SE
      # or a symmetric Wald-type arm CI. This is only an input conversion; the
      # actual SMD and sampling variance are still computed by metafor.
      .n_c_for_sd = dplyr::coalesce(as_num(n_control), .n_pairs),
      .n_t_for_sd = dplyr::coalesce(as_num(n_treatment), .n_pairs),
      .se_c_from_ci = dplyr::if_else(
        .measure != "proportion" &
          !is.na(as_num(`95CI_lower_control`)) & !is.na(as_num(`95CI_upper_control`)),
        (as_num(`95CI_upper_control`) - as_num(`95CI_lower_control`)) /
          (2 * stats::qnorm(0.975)), NA_real_),
      .se_t_from_ci = dplyr::if_else(
        .measure != "proportion" &
          !is.na(as_num(`95CI_lower_treatment`)) & !is.na(as_num(`95CI_upper_treatment`)),
        (as_num(`95CI_upper_treatment`) - as_num(`95CI_lower_treatment`)) /
          (2 * stats::qnorm(0.975)), NA_real_),

      .sd_c_use = dplyr::case_when(
        .measure == "proportion" ~ .sd_c_bern,
        !is.na(as_num(SD_control)) ~ as_num(SD_control),
        !is.na(as_num(SE_control)) & !is.na(.n_c_for_sd) & .n_c_for_sd > 0 ~
          as_num(SE_control) * sqrt(.n_c_for_sd),
        !is.na(.se_c_from_ci) & !is.na(.n_c_for_sd) & .n_c_for_sd > 0 ~
          .se_c_from_ci * sqrt(.n_c_for_sd),
        TRUE ~ NA_real_
      ),
      .sd_t_use = dplyr::case_when(
        .measure == "proportion" ~ .sd_t_bern,
        !is.na(as_num(SD_treatment)) ~ as_num(SD_treatment),
        !is.na(as_num(SE_treatment)) & !is.na(.n_t_for_sd) & .n_t_for_sd > 0 ~
          as_num(SE_treatment) * sqrt(.n_t_for_sd),
        !is.na(.se_t_from_ci) & !is.na(.n_t_for_sd) & .n_t_for_sd > 0 ~
          .se_t_from_ci * sqrt(.n_t_for_sd),
        TRUE ~ NA_real_
      ),
      .sd_c_source = dplyr::case_when(
        .measure == "proportion" & !is.na(.sd_c_bern) ~ "Bernoulli p(1-p)",
        !is.na(as_num(SD_control)) ~ "reported SD",
        !is.na(as_num(SE_control)) ~ "SE -> SD",
        !is.na(.se_c_from_ci) ~ "arm CI -> SE -> SD",
        TRUE ~ NA_character_
      ),
      .sd_t_source = dplyr::case_when(
        .measure == "proportion" & !is.na(.sd_t_bern) ~ "Bernoulli p(1-p)",
        !is.na(as_num(SD_treatment)) ~ "reported SD",
        !is.na(as_num(SE_treatment)) ~ "SE -> SD",
        !is.na(.se_t_from_ci) ~ "arm CI -> SE -> SD",
        TRUE ~ NA_character_
      )
    )
}

df_performance <- .prepare_es_inputs(df_performance)

# Audit boundary proportions: the standardized mean difference is undefined if
# the average raw-score SD is zero (e.g., both arms exactly 0% or 100%).
boundary_prop <- df_performance %>%
  dplyr::filter(.measure == "proportion",
                (!is.na(.p_c) & .p_c %in% c(0, 1)) |
                  (!is.na(.p_t) & .p_t %in% c(0, 1))) %>%
  dplyr::select(dplyr::any_of(c("Author", "ES_ID", ".p_c", ".p_t")))
if (nrow(boundary_prop) > 0) {
  message("[effect sizes] Proportion rows at 0 or 1 detected. metafor may return NA/very unstable SMDs if the standardizing SD is zero; review these rows:")
  print(as.data.frame(boundary_prop))
}

# ---- r from paired binary discordant counts --------------------------------
# For a paired 0/1 outcome, the discordant cells imply the variance of the
# individual-level difference. Together with the two marginal Bernoulli
# variances, this identifies the within-pair correlation used by SMCRPH.
.r_from_bc_one <- function(p_t, p_c, b, c, n) {
  vals <- c(p_t, p_c, b, c, n)
  if (any(!is.finite(vals)) || n <= 1) return(NA_real_)
  sd_t <- sqrt(p_t * (1 - p_t))
  sd_c <- sqrt(p_c * (1 - p_c))
  if (!is.finite(sd_t) || !is.finite(sd_c) || sd_t <= 0 || sd_c <= 0) return(NA_real_)

  # Var(Y_t - Y_c) for an individual paired Bernoulli observation.
  var_diff_unit <- (b + c) / n - ((c - b) / n)^2
  cov_tc <- (sd_t^2 + sd_c^2 - var_diff_unit) / 2
  r <- cov_tc / (sd_t * sd_c)
  if (!is.finite(r)) return(NA_real_)
  # Tiny excursions can arise from rounded percentages/counts.
  max(min(r, 0.999), -0.999)
}

# ---- observed r + donor pool / imputation -----------------------------------

# Derive the best available empirical within-pair correlation for each row.
#
# Priority:
#   1) r recalculated from discordant cells b/c
#   2) empirical/precomputed h_r_used
#
# Values explicitly labelled as imputed/assumed/fallback/sensitivity values in
# h_r_source are NOT allowed to become donors, avoiding circular imputation.
.derive_observed_r <- function(df) {

  df %>%
    dplyr::mutate(

      # Recalculate r from discordant cells wherever possible.
      .r_from_bc = mapply(
        .r_from_bc_one,
        .p_t, .p_c, .b, .c, .n_pairs
      ),

      # h_r_used is eligible as an empirical value unless its source explicitly
      # indicates that it was itself imputed or assumed.
      .r_precomputed_empirical = dplyr::if_else(
        is.finite(.r_precomputed) &
          !grepl(
            "imput|assum|fallback|sensitivity",
            .r_source_text
          ),
        .r_precomputed,
        NA_real_
      ),

      # Best empirical r available for this row.
      .r_observed = dplyr::case_when(

        !.paired ~ NA_real_,

        is.finite(.r_from_bc) &
          .r_from_bc > -1 & .r_from_bc < 1 ~
          .r_from_bc,

        is.finite(.r_precomputed_empirical) &
          .r_precomputed_empirical > -1 &
          .r_precomputed_empirical < 1 ~
          .r_precomputed_empirical,

        TRUE ~ NA_real_
      ),

      # Audit source of the observed r.
      .r_observed_source = dplyr::case_when(

        !.paired ~ NA_character_,

        is.finite(.r_from_bc) &
          .r_from_bc > -1 & .r_from_bc < 1 ~
          "from_bc",

        is.finite(.r_precomputed_empirical) &
          grepl("from_bc|discord", .r_source_text) ~
          "from_bc_precomputed",

        is.finite(.r_precomputed_empirical) &
          grepl("report", .r_source_text) ~
          "reported_h_r_used",

        is.finite(.r_precomputed_empirical) ~
          "h_r_used_other_empirical",

        TRUE ~ NA_character_
      )
    )
}


# Apply once before construction of the donor pool.
df_performance <- .derive_observed_r(df_performance)


# ---------------------------------------------------------------------------
# DONOR UNIT
# ---------------------------------------------------------------------------
# ID_2, because this corresponds to the independent paper /
# shared-case-corpus cluster used for CR2 inference. The three overlapping TUM
# papers therefore contribute only one donor unit.
#
# Optionally: If each of the 42 papers should count once instead, change
# this ONE line to:
# R_DONOR_UNIT <- "Study_ID"

R_DONOR_UNIT <- "ID_2"


# ---------------------------------------------------------------------------
# Build empirical donor pool
# ---------------------------------------------------------------------------
# First collapse multiple r estimates within the same donor unit and variance
# level. Therefore, a study/cluster with many reader-level rows cannot dominate
# the donor median simply because it reported more granular data.
#
# For the overall fallback, each donor unit contributes exactly ONE r value,
# irrespective of how many variance levels or effect-size rows it contains.
build_r_pool <- function(df, donor_unit = R_DONOR_UNIT) {

  stopifnot(donor_unit %in% names(df))

  d <- df %>%
    dplyr::filter(
      .paired,
      is.finite(.r_observed),
      .r_observed > -1,
      .r_observed < 1,
      !is.na(.data[[donor_unit]])
    ) %>%
    dplyr::mutate(
      .donor_unit = as.character(.data[[donor_unit]]),
      .var_level  = low(h_var_level)
    )


  # ONE donor value per study/cluster AND variance level.
  donors_by_level <- d %>%
    dplyr::filter(
      .var_level %in% c("case", "physician", "physician-case")
    ) %>%
    dplyr::group_by(.donor_unit, .var_level) %>%
    dplyr::summarise(
      r_donor = stats::median(.r_observed, na.rm = TRUE),
      n_source_rows = dplyr::n(),
      Study_ID = paste(sort(unique(Study_ID)), collapse = " | "),
      Author = paste(sort(unique(Author)), collapse = " | "),
      sources = paste(
        sort(unique(.r_observed_source)),
        collapse = " | "
      ),
      .groups = "drop"
    )


  # ONE donor value per study/cluster for the overall fallback.
  donors_overall <- d %>%
    dplyr::group_by(.donor_unit) %>%
    dplyr::summarise(
      r_donor = stats::median(.r_observed, na.rm = TRUE),
      n_source_rows = dplyr::n(),
      Study_ID = paste(sort(unique(Study_ID)), collapse = " | "),
      Author = paste(sort(unique(Author)), collapse = " | "),
      sources = paste(
        sort(unique(.r_observed_source)),
        collapse = " | "
      ),
      .groups = "drop"
    )


  list(

    case = donors_by_level %>%
      dplyr::filter(.var_level == "case") %>%
      dplyr::pull(r_donor),

    physician = donors_by_level %>%
      dplyr::filter(.var_level == "physician") %>%
      dplyr::pull(r_donor),

    physician_case = donors_by_level %>%
      dplyr::filter(.var_level == "physician-case") %>%
      dplyr::pull(r_donor),

    all = donors_overall$r_donor,

    # Keep complete audit tables as part of the object.
    donors_by_level = donors_by_level,
    donors_overall  = donors_overall
  )
}


# ---------------------------------------------------------------------------
# Imputation from study/cluster-aggregated donor pool
# ---------------------------------------------------------------------------
impute_r_value <- function(level, pool, override = NA_real_) {

  if (!is.na(override))
    return(override)

  level <- low(level)

  cand <- switch(
    level,
    "physician"      = pool$physician,
    "case"           = pool$case,
    "physician-case" = pool$physician_case,
    pool$all
  )

  # If no level-specific donors exist, use the overall study/cluster-level pool.
  if (length(cand) == 0)
    cand <- pool$all

  # Absolute last resort; should normally never be reached.
  if (length(cand) == 0)
    return(0.5)

  stats::median(cand, na.rm = TRUE)
}


# ---------------------------------------------------------------------------
# Assign final r to each paired effect
# ---------------------------------------------------------------------------
.resolve_r <- function(df, r_pool, r_override = NA_real_) {

  df %>%
    dplyr::mutate(

      .r_used_final = dplyr::case_when(

        !.paired ~ NA_real_,

        # Any empirically available r is retained.
        is.finite(.r_observed) ~
          .r_observed,

        # Only genuinely missing r values are imputed.
        TRUE ~ mapply(
          impute_r_value,
          h_var_level,
          MoreArgs = list(
            pool = r_pool,
            override = r_override
          )
        )
      ),

      .r_flag = dplyr::case_when(

        !.paired ~ "independent",

        .r_observed_source %in%
          c("from_bc", "from_bc_precomputed") ~
          "from_bc",

        .r_observed_source %in%
          c("reported", "reported_h_r_used") ~
          "reported",

        .r_observed_source ==
          "h_r_used_other_empirical" ~
          "h_r_used_empirical",

        TRUE ~ "imputed"
      )
    )
}

# ---- ICC / design-effect correction ----------------------------------------
# This is an approximation for the additional reader clustering axis. It is
# intentionally kept separate from r, which concerns the paired conditions.
ICC_ANCHOR <- 0.027
ICC_GRID   <- c(0, 0.027, 0.05, 0.10, 0.20)

# Every effect size uses the marginal d_av estimand: point estimate and sampling
# variance are built from the raw arm data, and reader clustering is handled by
# the design effect. Reported reader-cluster-adjusted odds ratios are not used as
# effect sizes -- they sit on the logistic latent scale and are frequently
# conditional rather than marginal estimates. The Reader_cluster_adjusted_SD_SE_CI
# notes are retained as documentation but no longer drive any calculation.
.derive_var_source <- function(df) {
  rep("raw_n", nrow(df))
}

# icc_override_scope controls which rows a sensitivity ICC is forced onto:
#   "imputed_only" (default): rows WITH a study-specific h_icc keep it, so the
#      grid value ICC = ICC_ANCHOR reproduces the primary analysis exactly.
#   "all": the override replaces every ICC, including study-specific ones. Use
#      this only for the strict lower bound (ICC = 0 everywhere).
.compute_de <- function(df, icc_override = NA_real_,
                        icc_override_scope = c("imputed_only", "all")) {
  icc_override_scope <- match.arg(icc_override_scope)
  df %>%
    dplyr::mutate(
      .var_source = .derive_var_source(.),
      .icc_used = dplyr::case_when(
        !is.na(icc_override) & icc_override_scope == "all" ~ icc_override,
        !is.na(as_num(h_icc))                              ~ as_num(h_icc),
        !is.na(icc_override)                               ~ icc_override,
        TRUE                                               ~ ICC_ANCHOR
      ),
      .n_phys_pair = as_num(n_control_physician),
      .m_bar = dplyr::case_when(
        .paired & !is.na(.n_pairs) & !is.na(.n_phys_pair) & .n_phys_pair > 0 ~
          .n_pairs / .n_phys_pair,
        !.paired & !is.na(as_num(n_treatment)) & !is.na(as_num(n_control)) &
          !is.na(as_num(n_treatment_physician)) & !is.na(as_num(n_control_physician)) &
          (as_num(n_treatment_physician) + as_num(n_control_physician)) > 0 ~
          (as_num(n_treatment) + as_num(n_control)) /
            (as_num(n_treatment_physician) + as_num(n_control_physician)),
        TRUE ~ NA_real_
      ),
      .DE = dplyr::case_when(
        low(h_var_level) == "physician-case" & .var_source == "raw_n" &
          !is.na(.m_bar) & .m_bar > 1 ~ 1 + (.m_bar - 1) * .icc_used,
        TRUE ~ 1
      ),
      .de_applied = .DE > 1
    )
}

# ---- metafor wrappers -------------------------------------------------------
.mf_escalc <- function(measure, args, correct = TRUE) {
  out <- tryCatch(
    do.call(metafor::escalc,
            c(list(measure = measure, correct = correct, append = FALSE), args)),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(list(yi = NA_real_, vi = NA_real_, error = conditionMessage(out)))
  }
  list(yi = as.numeric(out$yi[1]), vi = as.numeric(out$vi[1]), error = NA_character_)
}

# Assemble the standard effect-size row from a single escalc measure: bias-
# corrected g, uncorrected d, its sampling variance, and the audit labels.
.emit_es <- function(measure, args, source_label, measure_label = measure) {
  g <- .mf_escalc(measure, args, correct = TRUE)
  d <- .mf_escalc(measure, args, correct = FALSE)
  tibble::tibble(
    cohens_d        = d$yi,
    hedges_g        = g$yi,
    vi_metafor      = g$vi,
    metafor_measure = measure_label,
    es_source       = source_label,
    es_error        = dplyr::coalesce(g$error, d$error)
  )
}

# ---- one-row effect-size router --------------------------------------------
# Priority intentionally favors direct means/SDs over converted statistics.
.calc_metafor_es_one <- function(row) {
  paired  <- isTRUE(row$.paired[[1]])
  measure <- as.character(row$.measure[[1]])

  mt <- if (measure == "proportion") row$.p_t[[1]] else as_num(row$M_treatment[[1]])
  mc <- if (measure == "proportion") row$.p_c[[1]] else as_num(row$M_control[[1]])
  st <- as_num(row$.sd_t_use[[1]])
  sc <- as_num(row$.sd_c_use[[1]])
  nt <- as_num(row$n_treatment[[1]])
  nc <- as_num(row$n_control[[1]])
  np <- as_num(row$.n_pairs[[1]])
  rr <- as_num(row$.r_used_final[[1]])

  # 1) Paired / within-subject: raw-score standardized mean change.
  # SMCRPH uses sqrt((sd_t^2 + sd_c^2)/2) as denominator and incorporates r in vi.
  if (paired && all(is.finite(c(mt, mc, st, sc, np, rr))) &&
      st > 0 && sc > 0 && np > 1 && rr > -1 && rr < 1) {
    args <- list(m1i = mt, m2i = mc, sd1i = st, sd2i = sc,
                 ri = rr, ni = np)
    return(.emit_es("SMCRPH", args, paste0("paired_", measure, "_means_SD")))
  }

  # 2) Independent groups from means + SDs.
  if (!paired && all(is.finite(c(mt, mc, st, sc, nt, nc))) &&
      st > 0 && sc > 0 && nt > 1 && nc > 1) {
    # For binary 0/1 accuracy, SMDH exactly preserves the intended average-
    # variance denominator sqrt((p_t(1-p_t)+p_c(1-p_c))/2).
    mf_measure <- if (measure == "proportion") "SMDH" else "SMD"
    args <- list(m1i = mt, m2i = mc, sd1i = st, sd2i = sc,
                 n1i = nt, n2i = nc)
    return(.emit_es(mf_measure, args, paste0("independent_", measure, "_means_SD")))
  }

  # 3) Directly reported Cohen's d (independent-group interpretation only).
  d_rep <- as_num(row$.cohens_d_reported[[1]])
  if (!paired && is.finite(d_rep) && is.finite(nt) && is.finite(nc) && nt > 1 && nc > 1) {
    args <- list(di = d_rep, n1i = nt, n2i = nc)
    return(.emit_es("SMD", args, "independent_reported_cohens_d"))
  }

  # 4) Independent-samples t statistic -> SMD (metafor performs conversion + J).
  tv <- as_num(row$t_value[[1]])
  if (!paired && is.finite(tv) && is.finite(nt) && is.finite(nc) && nt > 1 && nc > 1) {
    args <- list(ti = tv, n1i = nt, n2i = nc)
    return(.emit_es("SMD", args, "independent_t_statistic"))
  }

  # 5) Performance difference + 95% CI in an independent design.
  # conv.wald reconstructs the sampling variance; the implied t statistic is then
  # passed to escalc(SMD) for the standardized effect and its variance.
  pd <- as_num(row$Performance_Difference[[1]])
  pl <- as_num(row$`95CI_lower_performance_diff`[[1]])
  pu <- as_num(row$`95CI_upper_performance_diff`[[1]])
  if (!paired && all(is.finite(c(pd, pl, pu, nt, nc))) && nt > 1 && nc > 1) {
    wd <- tryCatch(
      metafor::conv.wald(out = pd, ci.lb = pl, ci.ub = pu, append = FALSE),
      error = function(e) e
    )
    if (!inherits(wd, "error") && is.finite(wd$vi[1]) && wd$vi[1] > 0) {
      t_imp <- as.numeric(wd$yi[1]) / sqrt(as.numeric(wd$vi[1]))
      args <- list(ti = t_imp, n1i = nt, n2i = nc)
      return(.emit_es("SMD", args, "independent_difference_CI_via_t"))
    }
    }

  # No scientifically compatible metafor route available for this row.
  # Paired t-only / F-only rows are intentionally NOT converted to SMCC here:
  # SMCC uses the SD of change scores and is a different standardization from
  # the d_av / SMCRPH estimand used in the primary analysis.
  tibble::tibble(
    cohens_d = NA_real_,
    hedges_g = NA_real_,
    vi_metafor = NA_real_,
    metafor_measure = NA_character_,
    es_source = "NO_VALID_METAFOR_ROUTE",
    es_error = NA_character_
  )
}

# ---- main wrapper -----------------------------------------------------------
safe_calc_effects <- function(df, r_pool, r_override = NA_real_, icc_override = NA_real_,
                              icc_override_scope = "imputed_only") {
  # Remove stale calculated outputs before a sensitivity refit. Input fields and
  # .cohens_d_reported are retained.
  stale_cols <- intersect(
    c("cohens_d", "hedges_g", "vi_metafor", "metafor_measure", "es_source",
      "es_error", "vi", "vi_g"),
    names(df)
  )
  if (length(stale_cols) > 0) df <- dplyr::select(df, -dplyr::all_of(stale_cols))

  df <- .prepare_es_inputs(df)
  df <- .derive_observed_r(df)
  df <- .resolve_r(df, r_pool, r_override = r_override)
  df <- .compute_de(df, icc_override = icc_override,
                  icc_override_scope = icc_override_scope)

  # Point estimate + ordinary row-level sampling variance from metafor::escalc().
  es_tab <- purrr::map_dfr(seq_len(nrow(df)), function(i) {
    .calc_metafor_es_one(df[i, , drop = FALSE])
  })
  out <- dplyr::bind_cols(df, es_tab)

  out %>%
    dplyr::mutate(
      vi = vi_metafor,
      # Final sampling variance used by the meta-analysis. The design effect
      # corrects precision only and never changes the point estimate.
      vi_g = vi * .DE
    )
}

r_pool <- build_r_pool(
  df_performance,
  donor_unit = R_DONOR_UNIT
)

df_performance <- safe_calc_effects(
  df_performance,
  r_pool
)


# -----------------------------------------------------------------------------
# Diagnostic: which within-pair r is actually imputed?
# -----------------------------------------------------------------------------
local({

  # Final imputation medians.
  pool_tbl <- tibble::tibble(
    var_level = c(
      "case",
      "physician",
      "physician-case",
      "all (fallback)"
    ),

    n_donors = c(
      length(r_pool$case),
      length(r_pool$physician),
      length(r_pool$physician_case),
      length(r_pool$all)
    ),

    imputed_r = c(
      impute_r_value("case", r_pool),
      impute_r_value("physician", r_pool),
      impute_r_value("physician-case", r_pool),
      impute_r_value("__fallback__", r_pool)
    )
  )

  message(
    "\nImputed r per variance level ",
    "(median of study/cluster-aggregated empirical r donors):"
  )

  print(as.data.frame(pool_tbl))


  # Which study/cluster contributes which r to each variance-level pool?
  message("\nEmpirical r donors after within-study/cluster aggregation:")

  print(
    as.data.frame(
      r_pool$donors_by_level %>%
        dplyr::arrange(.var_level, .donor_unit)
    )
  )


  # Overall fallback donors.
  message("\nOverall donor values: one value per study/cluster:")

  print(
    as.data.frame(
      r_pool$donors_overall %>%
        dplyr::arrange(.donor_unit)
    )
  )


  # Which effects actually required imputation?
  imp_rows <- df_performance %>%
    dplyr::filter(.r_flag == "imputed") %>%
    dplyr::select(
      dplyr::any_of(
        c(
          "Author",
          "Study_ID",
          "ID_2",
          "ES_ID",
          "h_var_level",
          ".n_pairs",
          ".r_used_final"
        )
      )
    )

  message("\nRows with imputed r:")

  print(as.data.frame(imp_rows))


  # Save complete audit trail.
  readr::write_csv(
    pool_tbl,
    file.path(output_dir, "r_imputation_medians.csv")
  )

  readr::write_csv(
    r_pool$donors_by_level,
    file.path(output_dir, "r_donor_pool_by_level.csv")
  )

  readr::write_csv(
    r_pool$donors_overall,
    file.path(output_dir, "r_donor_pool_overall.csv")
  )

  readr::write_csv(
    imp_rows,
    file.path(output_dir, "r_imputed_effects.csv")
  )
})


# ---- validation summary -----------------------------------------------------
val <- df_performance %>% dplyr::summarise(
  n_rows      = dplyr::n(),
  paired      = sum(.paired, na.rm = TRUE),
  independent = sum(!.paired, na.rm = TRUE),
  r_from_bc         = sum(.r_flag == "from_bc", na.rm = TRUE),
  r_reported        = sum(.r_flag == "reported", na.rm = TRUE),
  r_h_used_empirical = sum(.r_flag == "h_r_used_empirical", na.rm = TRUE),
  r_imputed         = sum(.r_flag == "imputed", na.rm = TRUE),
  g_missing   = sum(is.na(hedges_g)),
  vi_missing  = sum(is.na(vi_g))
)
message("[effect sizes] Metafor effect-size / variance validation summary:")
print(as.data.frame(val))

# Duplicate inner-cluster identifiers (flag only; not silently altered)
dup_ids <- df_performance %>%
  dplyr::filter(!is.na(Experiment_ID), !is.na(Effect_ID)) %>%
  dplyr::count(Experiment_ID, Effect_ID, name = "n") %>%
  dplyr::filter(n > 1)
if (nrow(dup_ids) > 0) {
  message("[effect sizes] Duplicate (Experiment_ID, Effect_ID) combinations (review):")
  print(as.data.frame(dup_ids))
}

message("Imputed-r rows (covered by r-sensitivity analysis):")
print(df_performance %>%
        dplyr::filter(.r_flag == "imputed") %>%
        dplyr::select(dplyr::any_of(c("Author", "ES_ID", "h_var_level",
                                      ".n_pairs", ".r_used_final"))))

# Rows for which no scientifically compatible metafor route was available.
no_route <- df_performance %>%
  dplyr::filter(!is.finite(hedges_g) | !is.finite(vi_g)) %>%
  dplyr::select(dplyr::any_of(c(
    "Author", "ES_ID", ".design", ".paired", ".measure", "es_source",
    "M_treatment", "M_control", "SD_treatment", "SD_control",
    "SE_treatment", "SE_control", "n_treatment", "n_control", "h_n_pairs",
    "t_value", "F_value", "Odds_Ratio", "Performance_Difference", "es_error"
  )))
if (nrow(no_route) > 0) {
  message("[effect sizes] Rows with no valid/complete metafor SMD route (not silently used):")
  print(as.data.frame(no_route))
}

# Effect-size route audit ------------------------------------------------------
route_summary <- df_performance %>%
  dplyr::count(es_source, metafor_measure, name = "n_rows") %>%
  dplyr::arrange(dplyr::desc(n_rows))
message("[effect sizes] Effect-size calculation routes:")
print(as.data.frame(route_summary))
readr::write_csv(route_summary, file.path(output_dir, "effect_size_route_summary.csv"))

readr::write_csv(
  df_performance %>%
    dplyr::select(dplyr::any_of(c(
      "Author", "ES_ID", ".design", ".paired", ".measure",
      "es_source", "metafor_measure", "cohens_d", "hedges_g", "vi_metafor",
      ".sd_c_source", ".sd_t_source", ".r_from_bc", ".r_reported_direct", ".r_precomputed",
      ".r_observed", ".r_observed_source", ".r_flag", ".r_used_final", "h_var_level",
      ".var_source", ".m_bar", ".icc_used", ".DE", "vi_g", "es_error"
    ))),
  file.path(output_dir, "effect_size_metafor_audit.csv")
)

# ---- design-effect / variance-source validation ----------------------------
de_summary <- df_performance %>% dplyr::summarise(
  physician_case   = sum(low(h_var_level) == "physician-case", na.rm = TRUE),
  de_applied       = sum(.de_applied, na.rm = TRUE),
  icc_from_h_icc   = sum(!is.na(as_num(h_icc)), na.rm = TRUE),
  icc_anchor_used  = sum(.de_applied & is.na(as_num(h_icc)), na.rm = TRUE),
  DE_median        = round(stats::median(.DE[.de_applied], na.rm = TRUE), 3),
  DE_max           = round(max(c(.DE[.de_applied], 1), na.rm = TRUE), 3)
)
message("[design effect] Design-effect / variance-source summary (ICC anchor = ", ICC_ANCHOR, "):")
print(as.data.frame(de_summary))

mbar_missing <- df_performance %>%
  dplyr::filter(low(h_var_level) == "physician-case", .var_source == "raw_n",
                (is.na(.m_bar) | .m_bar <= 1)) %>%
  dplyr::select(dplyr::any_of(c("Author", "ES_ID", ".n_pairs",
                                "n_control_physician", "n_treatment", "n_control", ".m_bar")))
if (nrow(mbar_missing) > 0) {
  message("[design effect] physician-case raw_n rows with unusable m_bar (DE = 1 fallback; verify coding):")
  print(as.data.frame(mbar_missing))
}

level_mismatch <- df_performance %>%
  dplyr::filter(low(h_var_level) != "physician-case",
                !is.na(as_num(N_physicians_total)) & as_num(N_physicians_total) > 1,
                !is.na(as_num(N_cases_total)) & as_num(N_cases_total) > 1) %>%
  dplyr::select(dplyr::any_of(c("Author", "ES_ID", "h_var_level",
                                "N_physicians_total", "N_cases_total")))
if (nrow(level_mismatch) > 0) {
  message("[design effect] Rows with >1 physician AND >1 case but h_var_level != physician-case (review):")
  print(as.data.frame(level_mismatch))
}

readr::write_csv(
  df_performance %>%
    dplyr::select(dplyr::any_of(c(
      "Author", "ES_ID", "es_source", "metafor_measure", "h_var_level",
      ".var_source", ".m_bar", ".icc_used", ".DE", ".de_applied",
      "vi_metafor", "vi_g"
    ))),
  file.path(output_dir, "design_effect_audit.csv")
)


# Simple checks
"hedges_g" %in% names(df_performance)
if ("hedges_g" %in% names(df_performance)) sum(!is.na(df_performance$hedges_g))

intersect(c("cohens_d", "hedges_g", "vi", "vi_g"), names(df_performance))
colSums(!is.na(df_performance[, c("cohens_d","hedges_g","vi","vi_g")]))

# Inspect ES table (optional CSV)
inspect_effects <- dplyr::bind_rows(
  df_performance %>%
    dplyr::select(dplyr::any_of(c("ID_2","cohens_d","hedges_g"))) %>%
    dplyr::mutate(dataset = "performance")
)

write_csv(inspect_effects, file.path(output_dir, "all_effect_sizes.csv"))

# Main multilevel model (REML) + CR2-robust SEs
required_cols <- c("hedges_g","vi_g","Experiment_ID","Effect_ID")
stopifnot(all(required_cols %in% names(df_performance)))

df_es <- df_performance %>%
  dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(Experiment_ID), !is.na(Effect_ID)) %>%
  dplyr::mutate(ES_ID = if ("ES_ID" %in% names(.)) ES_ID else paste0(Experiment_ID, "_", Effect_ID))

# Collapse rare levels per moderator (for stable inference)
collapse_by_min_exp <- function(df, var, id_cluster = "ID_2",
                                min_exp = 2, other_label = "Other") {
  v  <- rlang::sym(var)
  vn <- rlang::as_name(v)

  tab <- df |> dplyr::distinct(!!v, .data[[id_cluster]]) |> dplyr::count(!!v, name = "n_exp")
  keep_levels <- tab |> dplyr::filter(n_exp >= min_exp) |> dplyr::pull(!!v) |> as.character()

  df |>
    dplyr::mutate(
      "{vn}_desc"  := factor(!!v),
      "{vn}_infer" := forcats::fct_other(
        factor(!!v), keep = keep_levels, other_level = other_label
      ) |> droplevels()
    )
}

# =============================================================================
# Shared model-fitting and extraction helpers
# -----------------------------------------------------------------------------
# The whole analysis uses ONE model specification: a REML multilevel model with
# effects nested in experiments, made cluster-robust (CR2 / RVE) at the paper /
# shared-case-corpus level (ID_2). These helpers hold that specification in one
# place so every model in the script (main, sensitivity, moderator, LOO, Egger)
# is fit identically (exception: 4-level model used as sensitivity analysis). 
# All reported inference is taken from the CR2 objects;
# the base rma.mv object is only queried for point estimates, variance
# components, the design matrix and weights.
cr2 <- function(mv, cluster) {
  metafor::robust(mv, cluster = cluster, clubSandwich = TRUE, adjust = TRUE)
}

# First matching column name from a set of candidates (handles the differing
# clubSandwich output labels across versions, e.g. df_Satt vs d.f. vs df).
pick_col <- function(x, cands) intersect(cands, colnames(x))[1]

fit_mv <- function(data, mods = NULL,
                   random = ~ 1 | Experiment_ID/Effect_ID,
                   slab   = if ("ES_ID" %in% names(data)) data$ES_ID else NULL) {
  # metafor uses non-standard evaluation for `mods`/`slab` and errors if they are
  # passed as a variable that resolves to NULL. Intercept-only models must OMIT
  # `mods` entirely (as the original code did), so build the call arguments and
  # add mods/slab only when supplied.
  args <- list(
    yi = data$hedges_g, V = data$vi_g, data = data,
    random = random, method = "REML", tdist = TRUE, level = 95
  )
  if (!is.null(mods)) args$mods <- mods
  if (!is.null(slab)) args$slab <- slab
  do.call(metafor::rma.mv, args)
}

# Base model plus its CR2 robustification as a pair (mv = base, rb = robust).
fit_cr2 <- function(data, mods = NULL,
                    random = ~ 1 | Experiment_ID/Effect_ID,
                    cluster = data$ID_2) {
  mv <- fit_mv(data, mods = mods, random = random)
  list(mv = mv, rb = cr2(mv, cluster))
}

# Pull a scalar estimate + CR2/Satterthwaite (or normal-fallback) CI and p from a
# robust object. Returns a list; NULL input yields all-NA (used by the LOO loops).
.extract_rb <- function(rb_obj) {
  if (is.null(rb_obj)) return(list(beta=NA_real_, se=NA_real_, stat=NA_real_, df=NA_real_, p=NA_real_, lo=NA_real_, hi=NA_real_))
  beta <- as.numeric(coef(rb_obj))[1]
  Vrb  <- as.matrix(vcov(rb_obj))
  se   <- sqrt(pmax(diag(Vrb)[1], 0))
  ddf  <- rb_obj$ddf
  has_t <- !is.null(ddf) && is.finite(ddf[1]) && ddf[1] > 0
  stat <- beta / se
  if (has_t) { p <- 2 * stats::pt(abs(stat), df = ddf[1], lower.tail = FALSE); crt <- stats::qt(0.975, df = ddf[1]) }
  else       { p <- 2 * stats::pnorm(abs(stat), lower.tail = FALSE);            crt <- stats::qnorm(0.975); ddf <- NA_real_ }
  lo <- beta - crt * se
  hi <- beta + crt * se
  list(beta=beta, se=se, stat=stat, df=if(is.list(ddf)) ddf[[1]] else ddf, p=p, lo=lo, hi=hi)
}

# Typical within-study variance term for Viechtbauer's multilevel I^2, i.e.
# (k - p) / tr(P). Shared by the heterogeneity summary and the LOO I^2.
.i2_typical_var <- function(fit, vi_vec) {
  W    <- diag(1 / vi_vec)
  X    <- model.matrix(fit)
  XtWX <- t(X) %*% W %*% X
  P    <- W - W %*% X %*% solve(XtWX) %*% t(X) %*% W
  (fit$k - fit$p) / sum(diag(P))
}

.i2_total <- function(fit, vi_vec) {
  tau2_total <- sum(fit$sigma2)
  out <- try(100 * tau2_total / (tau2_total + .i2_typical_var(fit, vi_vec)), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_) else return(as.numeric(out))
}

# =============================================================================
# Base meta-analytic model and Forest plot
# =============================================================================

main_model        <- fit_mv(df_es)
main_model_robust <- cr2(main_model, df_es$ID_2)

# Summary snippets for plotting
k        <- main_model_robust$k
n_pos_es <- sum(df_es$hedges_g > 0, na.rm = TRUE)
n_neg_es <- sum(df_es$hedges_g <= 0, na.rm = TRUE)
es_color <- ifelse(df_es$hedges_g > 0, "darkgreen", "darkred")

neg_label_part <- "≤ 0 effects: "
pos_label_part <- "> 0 effects: "
neg_label <- paste0(neg_label_part, n_neg_es, ", ", round(n_neg_es/k * 100, 1), "%)")
pos_label <- paste0(pos_label_part, n_pos_es, ", ", round(n_pos_es/k * 100, 1), "%)")
neg_label_pos <- min((n_pos_es + k)/2, k - 22)
pos_label_pos <- n_pos_es/2

estimate  <- as.numeric(coef(main_model_robust))
est_color <- ifelse(estimate > 0, "darkgreen", "darkred")
caption   <- paste0(
  "Model Estimate: g = ", round(estimate, 2),
  " [", round(main_model_robust$ci.lb, 2), ", ", round(main_model_robust$ci.ub, 2), "]"
)

print(summary(main_model_robust))

# =============================================================================
# r-SENSITIVITY ANALYSIS over the imputed paired rows
#   Refits the main multilevel + CR2 model forcing a fixed within-pair r for the
#   IMPUTED rows across a grid. Reported-r and b/c-derived rows are untouched.
#   Shows whether the pooled estimate is robust to the assumed correlation.
# =============================================================================
.fit_main_cr2 <- function(data) fit_cr2(data)$rb
r_grid        <- c(0.0, 0.3, 0.5, 0.7, 0.9)
sensitivity_r <- purrr::map_dfr(r_grid, function(rr) {
  d <- safe_calc_effects(df_performance, r_pool, r_override = rr) %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g),
                  !is.na(Experiment_ID), !is.na(Effect_ID))
  fit <- .fit_main_cr2(d)
  tibble::tibble(r_imputed = rr,
                 k         = fit$k,
                 estimate  = as.numeric(coef(fit)),
                 ci_lb     = as.numeric(fit$ci.lb),
                 ci_ub     = as.numeric(fit$ci.ub),
                 tau2      = sum(fit$sigma2))
})
message("[r-sensitivity] Sensitivity of pooled g to the imputed within-pair correlation r:")
print(as.data.frame(sensitivity_r))
readr::write_csv(sensitivity_r, file.path(output_dir, "r_imputation_sensitivity.csv"))

# =============================================================================
# ICC / DESIGN-EFFECT SENSITIVITY ANALYSIS
#   Refits the main multilevel + CR2 model across an ICC grid
#   he ICC grid replaces the anchor ICC only for rows without a study-specific h_icc.
#   Study-specific ICCs remain fixed. The final ICC = 0 scenario removes the design-effect correction from all rows.
#   r-handling is held at its normal per-row value. Shows robustness of pooled g to the assumed ICC.
#   NB: it is (m_bar - 1) * ICC that matters -- a "small" ICC can still halve the
#   effective N when m_bar is large (e.g. m_bar = 500 -> DE ~ 14.5 at ICC = 0.027).
# =============================================================================
.run_icc_scenario <- function(ii, scope, label) {
  d <- safe_calc_effects(df_performance, r_pool,
                         icc_override = ii, icc_override_scope = scope) %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g),
                  !is.na(Experiment_ID), !is.na(Effect_ID))
  fit <- .fit_main_cr2(d)
  tibble::tibble(icc        = label,
                 scope      = scope,
                 k          = fit$k,
                 n_deflated = sum(d$.de_applied, na.rm = TRUE),
                 estimate   = as.numeric(coef(fit)),
                 ci_lb      = as.numeric(fit$ci.lb),
                 ci_ub      = as.numeric(fit$ci.ub),
                 tau2       = sum(fit$sigma2))
}

sensitivity_icc <- dplyr::bind_rows(
  # Grid: study-specific ICCs are held fixed, so the ICC_ANCHOR row reproduces
  # the primary analysis exactly.
  purrr::map_dfr(ICC_GRID, ~ .run_icc_scenario(.x, "imputed_only", as.character(.x))),
  # Strict lower bound: no reader-clustering correction anywhere.
  .run_icc_scenario(0, "all", "0 (all rows, strict lower bound)")
)
message("[ICC sensitivity] Sensitivity of pooled g to the assumed ICC (base case = ", ICC_ANCHOR, "):")
print(as.data.frame(sensitivity_icc))
readr::write_csv(sensitivity_icc, file.path(output_dir, "icc_designeffect_sensitivity.csv"))

# Guard: the anchor row of the grid MUST equal the primary estimate.
local({
  anchor <- sensitivity_icc %>% dplyr::filter(icc == as.character(ICC_ANCHOR),
                                              scope == "imputed_only")
  if (nrow(anchor) == 1 &&
      abs(anchor$estimate - as.numeric(coef(main_model_robust))) > 1e-6)
    warning("[ICC sensitivity] The ICC = ", ICC_ANCHOR,
            " row does not reproduce the primary estimate -- investigate before reporting.")
})

# Forest plot data
author_col <- "Author"

model_est <- as.numeric(coef(main_model_robust))
model_lo  <- as.numeric(main_model_robust$ci.lb)
model_hi  <- as.numeric(main_model_robust$ci.ub)

# For an intercept-only rma.mv model, row-sum weights are the actual weights
# assigned to the individual estimates when estimating the model intercept.
# metafor returns them directly in percent (sum = 100).

model_rowsum_weights <- as.numeric(
  stats::weights(main_model, type = "rowsum")
)
stopifnot(length(model_rowsum_weights) == nrow(df_es))

df_es <- df_es %>%
  dplyr::mutate(w_raw = NA_real_, w_pct = model_rowsum_weights)

# =============================================================================
# Sensitivity analyses: publication status + study design
# =============================================================================

fit_sens <- function(dat, dimension, subgroup) {

  rb <- fit_cr2(dat)$rb

  tibble::tibble(
    dimension = dimension,
    subgroup = subgroup,
    k_effects = nrow(dat),
    k_experiments = dplyr::n_distinct(dat$Experiment_ID),
    k_clusters = dplyr::n_distinct(dat$ID_2),
    g = as.numeric(coef(rb))[1],
    ci_low = as.numeric(rb$ci.lb)[1],
    ci_high = as.numeric(rb$ci.ub)[1]
  )
}

# Publication status
df_es <- df_es %>%
  dplyr::mutate(
    .pub = dplyr::case_when(
      stringr::str_detect(low(Type_of_study), "pre.?print|medrxiv|arxiv|biorxiv") ~ "Preprint",
      stringr::str_detect(low(Type_of_study), "peer.?review") ~ "Peer-reviewed",
      TRUE ~ NA_character_
    )
  )

if (anyNA(df_es$.pub))
  stop("Unclassified Type_of_study value(s): ",
       paste(unique(df_es$Type_of_study[is.na(df_es$.pub)]), collapse = ", "))

# ONE combined output
sensitivity_subgroups <- dplyr::bind_rows(
  fit_sens(dplyr::filter(df_es, .pub == "Peer-reviewed"),
           "Publication status", "Peer-reviewed"),
  fit_sens(dplyr::filter(df_es, .pub == "Preprint"),
           "Publication status", "Preprint"),
  fit_sens(dplyr::filter(df_es, .design == "parallel"),
           "Study design", "Parallel"),
  fit_sens(dplyr::filter(df_es, .design == "within_subject"),
           "Study design", "Within-subject")
)

print(sensitivity_subgroups)

readr::write_csv(
  sensitivity_subgroups,
  file.path(output_dir, "sensitivity_publication_and_design.csv")
)

# =============================================================================
# CONSISTENT FOREST-PLOT LABELS
# =============================================================================

first_author_fun <- function(x) {
  x <- ifelse(is.na(x) | x == "", "Unlabelled", x)
  x <- stringr::str_squish(as.character(x))
  x <- stringr::str_remove(x, "\\bet al\\.?$")
  x <- stringr::str_split_fixed(x, pattern = ",|;| and | & ", n = 2)[, 1]
  stringr::str_split_fixed(x, pattern = "\\s+", n = 2)[, 1]
}

letter_suffix <- function(k) {
  labs <- letters
  if (k <= length(labs)) return(labs[seq_len(k)])

  vapply(seq_len(k), function(i) {
    q <- (i - 1) %/% length(labs)
    r <- (i - 1) %% length(labs) + 1
    paste0(if (q == 0) "" else labs[q], labs[r])
  }, character(1))
}

pick_year_column <- function(nms) {
  cand <- c(
    "Publication_Year", "Year", "Publication.Year",
    "publication_year", "Pub_Year", "Study_Year"
  )
  hit <- cand[cand %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

year_col_labels <- pick_year_column(names(df_es))

if (is.na(year_col_labels)) {
  year_vec_labels <- rep("", nrow(df_es))
} else {
  year_vec_labels <- trimws(as.character(df_es[[year_col_labels]]))
  year_vec_labels[
    is.na(year_vec_labels) |
      year_vec_labels %in% c("NA", "NR", "NaN", "na", "n/a", "")
  ] <- ""
}

forest_label_map <- df_es %>%
  dplyr::mutate(
    First_Author = first_author_fun(Author),
    Study_Year   = year_vec_labels
  ) %>%
  dplyr::arrange(
    First_Author, Study_Year,
    Experiment_ID, ES_ID
  ) %>%
  dplyr::group_by(First_Author, Study_Year) %>%
  dplyr::mutate(
    short_suffix = letter_suffix(dplyr::n())[dplyr::row_number()]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    short_base = ifelse(
      Study_Year == "",
      First_Author,
      paste(First_Author, Study_Year)
    ),
    short_label = paste0(short_base, short_suffix)
  ) %>%
  dplyr::select(ES_ID, short_label)

# -----------------------------------------------------------------------------
# Shared forest-plot building blocks (visualisation only; no model objects are
# touched). Used by the main forest (V1 compact + V2 labelled) and by the
# low-RoB forest, so their styling stays in one place.
# -----------------------------------------------------------------------------
FOREST_COL <- list(
  pos = "#2E7D32",  # CI entirely > 0 (benefit)
  neg = "#C2662E",  # CI entirely < 0 (harm)
  ns  = "#9AA0A6",  # CI crosses 0
  acc = "#C8860D",  # pooled band / dashed line
  gold_txt = "#B8860B", gold_fill = "#FFF7E0", zebra = "grey50"
)
FOREST_ARROW <- grid::arrow(angle = 22, length = grid::unit(0.9, "mm"),
                            ends = "last", type = "closed")

forest_theme <- function(tick_length = grid::unit(1.1, "pt")) {
  th <- theme_classic(base_size = 8, base_family = "") +
    theme(
      axis.line.y  = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y  = element_blank(), axis.title.y = element_blank(),
      axis.line.x  = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.x = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x  = element_text(size = 6.5, colour = "black"),
      axis.title.x = element_text(size = 7.5, margin = margin(t = 3)),
      panel.grid   = element_blank(),
      plot.margin  = margin(2, 4, 2, 4, "mm")
    )
  if (!is.null(tick_length)) th <- th + theme(axis.ticks.length = tick_length)
  th
}

# Data-driven x-limits with 2nd/98th-percentile trimming and pretty breaks.
forest_xlims <- function(plot_df, lo, hi) {
  q_lo <- as.numeric(quantile(plot_df$ci.lb, 0.02, na.rm = TRUE))
  q_hi <- as.numeric(quantile(plot_df$ci.ub, 0.98, na.rm = TRUE))
  pad  <- 0.05 * (q_hi - q_lo)
  xlim_lo <- min(q_lo, lo) - pad
  xlim_hi <- max(q_hi, hi) + pad
  xb <- scales::breaks_pretty(7)(c(xlim_lo, xlim_hi))
  xb <- xb[xb >= xlim_lo & xb <= xlim_hi]
  list(xlim_lo = xlim_lo, xlim_hi = xlim_hi, xspan = xlim_hi - xlim_lo, x_breaks = xb)
}

# Clamp CIs to the plot window and flag capped (arrowed) rows.
forest_cap <- function(plot_df, xlim_lo, xlim_hi) {
  plot_df %>% dplyr::mutate(
    seg_lo = pmax(ci.lb, xlim_lo), seg_hi = pmin(ci.ub, xlim_hi),
    cap_lo = ci.lb < xlim_lo,      cap_hi = ci.ub > xlim_hi,
    pt_x   = pmin(pmax(effect, xlim_lo), xlim_hi)
  )
}

# Cairo-aware exporters. The base pdf() device cannot draw the Turkish glyphs in
# Atakir/Camur (U+0131, c-cedilla); cairo_pdf embeds a font that carries them,
# and ragg/Cairo do the same for PNG. Falls back to the base device otherwise.
FOREST_USE_CAIRO <- isTRUE(capabilities("cairo")[["cairo"]])
forest_save_pdf <- function(file, plot, ...) {
  if (FOREST_USE_CAIRO) ggsave(file, plot, device = grDevices::cairo_pdf, ...)
  else                  ggsave(file, plot, ...)
}
forest_save_png <- function(file, plot, ...) {
  if (requireNamespace("ragg", quietly = TRUE)) ggsave(file, plot, ...)
  else if (FOREST_USE_CAIRO)                     ggsave(file, plot, type = "cairo", ...)
  else                                           ggsave(file, plot, ...)
}

# Labelled forest panel (study labels + weight column + pooled gold label).
# Shared by the main appendix figure (Figure 1 (with labels: Figure S3)) and the low-RoB figure (Figure S2); the two differ
# only in the pooled-label text/size and the axis tick length, which are passed
# in. plot_df must already carry y, effect, ci.lb/ub, seg_lo/hi, cap_lo/hi, pt_x,
# col, short_label and weight_lab.
draw_labeled_forest <- function(plot_df, est, lo, hi, lims, summary_text,
                                gold_size, theme = forest_theme(),
                                zebra = TRUE) {
  left_bound  <- lims$xlim_lo - 0.36 * lims$xspan
  right_bound <- lims$xlim_hi + 0.15 * lims$xspan
  y_first <- min(plot_df$y); y_last <- max(plot_df$y)
  head_y  <- y_first - 1.5;  rule_y <- y_first - 0.55
  top_y   <- y_first - 0.5;  bot_y  <- y_last + 0.5

  p <- ggplot(plot_df, aes(y = y))
  if (zebra)
    p <- p + geom_rect(data = dplyr::filter(plot_df, y %% 2 == 0),
                       aes(xmin = left_bound, xmax = right_bound,
                           ymin = y - 0.5, ymax = y + 0.5),
                       inherit.aes = FALSE, fill = FOREST_COL$zebra, alpha = 0.06)
  p +
    annotate("rect", xmin = lo, xmax = hi, ymin = top_y, ymax = bot_y,
             fill = FOREST_COL$acc, alpha = 0.08) +
    annotate("segment", x = 0, xend = 0, y = top_y, yend = bot_y,
             colour = "grey25", linewidth = 0.35) +
    annotate("segment", x = est, xend = est, y = top_y, yend = bot_y,
             colour = FOREST_COL$acc, linetype = "22", linewidth = 0.48) +
    geom_segment(aes(x = seg_lo, xend = seg_hi, yend = y, colour = col),
                 linewidth = 0.48, alpha = 0.68, lineend = "round") +
    geom_segment(data = dplyr::filter(plot_df, cap_hi),
                 aes(x = lims$xlim_hi - 1e-3, xend = lims$xlim_hi, y = y, yend = y, colour = col),
                 linewidth = 0.48, alpha = 0.68, arrow = FOREST_ARROW) +
    geom_segment(data = dplyr::filter(plot_df, cap_lo),
                 aes(x = lims$xlim_lo + 1e-3, xend = lims$xlim_lo, y = y, yend = y, colour = col),
                 linewidth = 0.48, alpha = 0.68, arrow = FOREST_ARROW) +
    geom_point(aes(x = pt_x, colour = col), size = 1.05, alpha = 0.95) +
    geom_text(aes(x = left_bound,  label = short_label), hjust = 0, size = 2.1, colour = "black") +
    geom_text(aes(x = right_bound, label = weight_lab),  hjust = 1, size = 2.1, colour = "black") +
    annotate("text", x = left_bound,  y = head_y, label = "Study",  hjust = 0, size = 2.9, fontface = "bold") +
    annotate("text", x = right_bound, y = head_y, label = "Weight", hjust = 1, size = 2.9, fontface = "bold") +
    annotate("segment", x = left_bound, xend = right_bound, y = rule_y, yend = rule_y,
             linewidth = 0.4, colour = "black") +
    annotate("label", x = est, y = y_last + 1.8, label = summary_text,
             colour = FOREST_COL$gold_txt, fill = FOREST_COL$gold_fill,
             label.size = 0.35, size = gold_size) +
    scale_colour_identity(guide = "none") +
    scale_y_reverse(expand = expansion(add = c(2, 2.5))) +
    scale_x_continuous(breaks = lims$x_breaks) +
    coord_cartesian(xlim = c(left_bound, right_bound), clip = "off") +
    labs(x = expression("Effect size (Hedges' " * italic(g) * ")")) +
    theme
}

# =============================================================================
# Forest plots (visualisation only)
# -----------------------------------------------------------------------------
# Two figures are drawn from the same fitted objects (main_model /
# main_model_robust); the model, data and pooled estimate are not touched:
#   V1  forest_multilevel_overall_gg_simple.pdf/png  -> compact main figure,
#       no per-row labels. Row -> study lookup is exported to forest_es_map.csv.
#   V2  forest_appendix_labels_weights.pdf/png       -> appendix figure with
#       short study labels and per-row weight %.
# =============================================================================

plot_filename2 <- "forest_multilevel_overall_gg_simple"

# Wrapped in local() so none of the plotting temporaries leak into the global
# environment or collide with objects used later in the script. The figures are
# side effects written to disk, so nothing downstream needs these objects.
local({

  ## ---- style constants (V1 compact figure; V2 uses draw_labeled_forest) -----
  COLOR_BY_SIG     <- TRUE   # colour points/CIs by CI significance (FALSE = grey)
  ZEBRA            <- TRUE   # faint alternating row shading in the labelled V2
  SUMMARY_DECIMALS <- 3      # decimals in the pooled gold label. Kept at 3 to
                             # match the current figure/manuscript (g = 0.213...).
                             # Set to 2 to follow the 2-dp reporting policy.

  PT_SIZE <- 0.70; CI_LW <- 0.7; CI_ALPHA <- 0.60
  ZERO_LW <- 0.35; POOLED_LW <- 0.48

  col_pos <- FOREST_COL$pos; col_neg <- FOREST_COL$neg; col_ns <- FOREST_COL$ns
  col_acc <- FOREST_COL$acc; gold_txt <- FOREST_COL$gold_txt; gold_fill <- FOREST_COL$gold_fill

  ## ---- pooled estimate (reuse the objects computed above) -------------------
  est <- model_est[1]; lo <- model_lo[1]; hi <- model_hi[1]
  summary_text <- sprintf(
    paste0("Overall Hedges' g = %.", SUMMARY_DECIMALS, "f [%.",
           SUMMARY_DECIMALS, "f, %.", SUMMARY_DECIMALS, "f]"),
    est, lo, hi)

  ## ---- per-row weights (self-diagnosing; prefers df_es$w_pct) ---------------
  ## Source order (first that yields finite, positive values wins):
  ##   (1) df_es$w_pct               = main model's own row-sum weights (set above)
  ##   (2) weights(main_model,"rowsum")
  ##   (3) weights(main_model,"diagonal")
  ##   (4) APPROX 1/(vi + tau^2)     = clearly flagged as NOT the model weights
  ## NB: the local vector is called `w_source` (NOT `w_raw`). df_es
  ## already carries an all-NA placeholder column `w_raw` (set upstream). Inside
  ## dplyr::mutate() a data column shadows an identically named env variable, so
  ## a local `w_raw` would be silently overwritten by that NA column -> blank
  ## weights. `w_source` has no column counterpart, and .env$ is used below too.
  n_rows <- nrow(df_es); w_source <- NULL; weight_source <- NA_character_
  if ("w_pct" %in% names(df_es)) {
    cand <- suppressWarnings(as.numeric(df_es$w_pct))
    if (any(is.finite(cand)) && sum(cand, na.rm = TRUE) > 0) {
      w_source <- cand; weight_source <- "df_es$w_pct (main-model row-sum weights)"
    }
  }
  if (is.null(w_source)) {
    try_w <- function(type) tryCatch({
      ww <- as.numeric(stats::weights(main_model, type = type))
      if (length(ww) == n_rows && all(is.finite(ww))) ww else NULL
    }, error = function(e) NULL, warning = function(w) NULL)
    ww <- try_w("rowsum")
    if (!is.null(ww)) { w_source <- ww; weight_source <- "weights(main_model,'rowsum')" }
    if (is.null(w_source)) {
      ww <- try_w("diagonal")
      if (!is.null(ww)) { w_source <- ww; weight_source <- "weights(main_model,'diagonal')" }
    }
  }
  if (is.null(w_source)) {
    tau2_tot <- sum(c(main_model$sigma2, main_model$tau2), na.rm = TRUE)
    w_source <- 1 / (df_es$vi_g + tau2_tot)
    weight_source <- sprintf("APPROX 1/(vi+tau^2=%.4f) -- NOT model weights", tau2_tot)
  }
  message("[forest] weights from: ", weight_source)

  ## ---- year vector (real NA and the strings NA/NR/NaN -> "") ----------------
  year_col <- pick_year_column(names(df_es))
  if (is.na(year_col)) {
    year_vec <- rep("", nrow(df_es))
  } else {
    year_vec <- trimws(as.character(df_es[[year_col]]))
    year_vec[is.na(year_vec) | year_vec %in% c("NA", "NR", "NaN", "na", "n/a", "")] <- ""
  }

  ## ---- assemble the plotting frame ------------------------------------------
  stopifnot(all(c("hedges_g", "vi_g", "Experiment_ID", "ES_ID", "Author") %in% names(df_es)))
  plot_df <- df_es %>%
    dplyr::mutate(
      effect       = hedges_g,
      se           = sqrt(vi_g),
      ci.lb        = effect - 1.96 * se,
      ci.ub        = effect + 1.96 * se,
      weight       = 100 * .env$w_source / sum(.env$w_source, na.rm = TRUE),
      weight_lab   = ifelse(is.finite(weight), sprintf("%.1f%%", weight), ""),
      Author_plot  = ifelse(is.na(Author) | Author == "", "Unlabelled",
                            stringr::str_squish(as.character(Author))),
      First_Author = first_author_fun(Author),
      Study_Year   = year_vec
    ) %>%
    dplyr::arrange(dplyr::desc(effect), Experiment_ID, ES_ID) %>%
    dplyr::mutate(y = dplyr::row_number())

  ## short labels: "AuthorYear" + a/b/c within the same author-year block
  plot_df <- plot_df %>%
    dplyr::arrange(First_Author, Study_Year, Experiment_ID, ES_ID) %>%
    dplyr::group_by(First_Author, Study_Year) %>%
    dplyr::mutate(short_suffix = letter_suffix(dplyr::n())[dplyr::row_number()]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      short_base  = ifelse(Study_Year == "", First_Author, paste(First_Author, Study_Year)),
      short_label = paste0(short_base, short_suffix),
      short_label = dplyr::coalesce(short_label, Author_plot, paste0("Study ", ES_ID))
    ) %>%
    dplyr::arrange(y) %>%
    dplyr::mutate(
      sign = dplyr::case_when(ci.lb > 0 ~ "pos", ci.ub < 0 ~ "neg", TRUE ~ "ns"),
      col  = if (COLOR_BY_SIG)
               dplyr::case_when(sign == "pos" ~ col_pos, sign == "neg" ~ col_neg, TRUE ~ col_ns)
             else col_ns
    )

  if (all(!nzchar(plot_df$weight_lab)))
    warning("[forest] all weight labels are empty; weight source was: ", weight_source)

  ## row -> study lookup companion for the label-free V1 (audit trail)
  readr::write_csv(
    plot_df %>% dplyr::transmute(row = y, ES_ID, Author, short_label,
                                 g = round(effect, 3), ci_lb = round(ci.lb, 3),
                                 ci_ub = round(ci.ub, 3), weight_pct = round(weight, 2)),
    file.path(output_dir, "forest_es_map.csv")
  )

  ## ---- x-limits + arrow-capping (shared helpers) ----------------------------
  lims    <- forest_xlims(plot_df, lo, hi)
  plot_df <- forest_cap(plot_df, lims$xlim_lo, lims$xlim_hi)
  n_es    <- nrow(plot_df)

  ## =========================================================================
  ## V1 - compact main figure (no per-row labels)
  ## =========================================================================
  p1 <- ggplot(plot_df, aes(y = y)) +
    annotate("rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
             fill = col_acc, alpha = 0.08) +
    geom_vline(xintercept = 0,   colour = "grey25", linewidth = ZERO_LW) +
    geom_vline(xintercept = est, colour = col_acc, linetype = "22", linewidth = POOLED_LW) +
    geom_segment(aes(x = seg_lo, xend = seg_hi, yend = y, colour = col),
                 linewidth = CI_LW, alpha = CI_ALPHA, lineend = "round") +
    geom_segment(data = dplyr::filter(plot_df, cap_hi),
                 aes(x = lims$xlim_hi - 1e-3, xend = lims$xlim_hi, y = y, yend = y, colour = col),
                 linewidth = CI_LW, alpha = CI_ALPHA, arrow = FOREST_ARROW) +
    geom_segment(data = dplyr::filter(plot_df, cap_lo),
                 aes(x = lims$xlim_lo + 1e-3, xend = lims$xlim_lo, y = y, yend = y, colour = col),
                 linewidth = CI_LW, alpha = CI_ALPHA, arrow = FOREST_ARROW) +
    geom_point(aes(x = pt_x, colour = col), size = PT_SIZE, alpha = 0.95) +
    annotate("label", x = est, y = max(plot_df$y) + 1.8, label = summary_text,
             colour = gold_txt, fill = gold_fill, label.size = 0.35, size = 2.9) +
    scale_colour_identity(guide = "none") +
    scale_y_reverse(expand = expansion(add = c(1, 4))) +
    scale_x_continuous(breaks = lims$x_breaks) +
    coord_cartesian(xlim = c(lims$xlim_lo, lims$xlim_hi), clip = "off") +
    labs(x = expression("Effect size (Hedges' " * italic(g) * ")")) +
    forest_theme()

  H1 <- max(6.0, min(11.2, n_es * 0.045 + 1.5))   # data-scaled, capped to A4 height
  forest_save_pdf(file.path(output_dir, paste0(plot_filename2, ".pdf")), p1,
                  width = 7.5, height = H1, units = "in", dpi = 600)
  forest_save_png(file.path(output_dir, paste0(plot_filename2, ".png")), p1,
                  width = 7.5, height = H1, units = "in", dpi = 300)

  ## =========================================================================
  ## V2 - appendix figure WITH short labels + weights
  ## =========================================================================
  p2 <- draw_labeled_forest(plot_df, est, lo, hi, lims, summary_text,
                            gold_size = 2.9, theme = forest_theme(), zebra = ZEBRA)
  H2 <- max(7.0, n_es * 0.115 + 1.6)
  forest_save_pdf(file.path(output_dir, "forest_appendix_labels_weights.pdf"), p2,
                  width = 8.2, height = H2, units = "in")
  forest_save_png(file.path(output_dir, "forest_appendix_labels_weights.png"), p2,
                  width = 8.2, height = H2, units = "in", dpi = 400)

  message("[forest] saved to: ", output_dir)
  message("  V1 (main):     ", plot_filename2, ".pdf/png  (", n_es, " ES)")
  message("  V2 (appendix): forest_appendix_labels_weights.pdf/png")
  message("  row map:       forest_es_map.csv")
})


# Heterogeneity (tau^2, I^2 total and per level) ------------------------------
# Based on Viechtbauer's multilevel I^2

tau2_vec <- main_model$sigma2
if (length(tau2_vec) == 2) {
  names(tau2_vec) <- c("between_experiments", "within_experiments")
} else {
  names(tau2_vec) <- paste0("component_", seq_along(tau2_vec))
}

tau2_total <- sum(tau2_vec)

denom    <- tau2_total + .i2_typical_var(main_model, df_es$vi_g)
I2_total <- 100 * tau2_total / denom
I2_lvls  <- 100 * tau2_vec   / denom

cat("\n--- Heterogeneity summary ---\n")
cat("tau^2 components:\n"); print(round(tau2_vec, 6))
cat("tau^2 total: ", round(tau2_total, 6), "\n\n")
cat("I^2 total:   ", round(I2_total, 2), "%\n", sep = "")
cat("I^2 by level:\n"); print(round(I2_lvls, 2))

hetero_out <- tibble::tibble(
  component = c(names(tau2_vec), "total"),
  tau2      = c(as.numeric(tau2_vec), tau2_total),
  I2_pct    = c(as.numeric(I2_lvls), I2_total)
)
readr::write_csv(hetero_out, file.path(output_dir, "heterogeneity_summary.csv"))

# =============================================================================
# Low-risk-of-bias sensitivity analysis  (SECONDARY / ROBUSTNESS)
# -----------------------------------------------------------------------------
# Refits the IDENTICAL main model on the subset of effect sizes that come from
# LOW-RoB experiments only, and draws a matching forest plot.
#
#   * Model spec unchanged: rma.mv(~ 1 | Experiment_ID/Effect_ID, REML) with CR2
#     robustification clustered at the PAPER level (ID_2).
#   * Exclusion is at the EXPERIMENT level: high-RoB experiments are dropped and
#     low-RoB experiments are kept -- even inside a paper whose experiments differ
#     in RoB (the eyecare Wu paper keeps Exp 3 and Exp 2b, drops Exp 2a).
#   * Exclusion (which rows enter the data) and CR2 clustering (how the retained
#     rows are grouped for robust SEs) are ORTHOGONAL, so the surviving Wu
#     experiments still sit inside the single Wu ID_2 cluster.
#
# SMALL-SAMPLE CAVEAT: the low-RoB subset has few paper clusters, so the CR2 /
# Satterthwaite CIs are wide and approximate. Report as robustness, not as a
# confirmatory estimate.
#
# RoB source: PROBAST workbook, sheet "Extraction form".
# The relevant columns are identified by their header text rather than by fixed
# column positions. The TUM merge (ID_2 11/12/21 -> 11) is replicated here so
# that the PROBAST study IDs align with df_es.
# =============================================================================

rob_path <- file.path(data_dir, "20260828_PROBAST_vFinal.xlsm")
if (!file.exists(rob_path))
  stop("[RoB sensitivity] PROBAST file not found next to the coding sheet: ", rob_path)

# Read the WHOLE sheet raw and locate columns by HEADER TEXT (robust to any
# column shift from re-saving). Header layout: Excel row 2 carries the group
# labels ("Overall: RoB judgement" vs "...Applicability judgement"), Excel row 3
# carries the field headers ("Study ID", "Author", "Year"); data starts row 4.
rob_all <- tryCatch(
  readxl::read_xlsx(rob_path, sheet = "Extraction form", col_names = FALSE),
  error = function(e) stop("[RoB sensitivity] Could not read the .xlsm (", conditionMessage(e),
                           "). If readxl rejects it, re-save the sheet as .xlsx."))
if (nrow(rob_all) < 4)
  stop("[RoB sensitivity] 'Extraction form' has too few rows (", nrow(rob_all), ").")

hdr2 <- tolower(trimws(as.character(unlist(rob_all[2, ], use.names = FALSE))))
hdr3 <- tolower(trimws(as.character(unlist(rob_all[3, ], use.names = FALSE))))
hdr2[is.na(hdr2)] <- ""; hdr3[is.na(hdr3)] <- ""

find_col <- function(hit, what) {
  w <- which(hit)
  if (length(w) == 0) stop("[RoB sensitivity] Could not find the ", what, " column in the PROBAST header.")
  w[1]
}
# Overall RoB = row-2 label with 'overall' + 'rob' but NOT 'applicab' (that is AR)
col_rob  <- find_col(grepl("overall", hdr2) & grepl("rob", hdr2) & !grepl("applicab", hdr2),
                     "'Overall: RoB judgement'")
col_sid  <- find_col(hdr3 == "study id", "'Study ID'")
col_auth <- find_col(hdr3 == "author",   "'Author'")
col_year <- find_col(hdr3 == "year",     "'Year'")

rob_dat <- rob_all[-(1:3), , drop = FALSE]     # data rows only (Excel row 4+)
rob_tbl <- tibble::tibble(
  ID_2_raw    = suppressWarnings(as.numeric(rob_dat[[col_sid]])),
  Author_rob  = as.character(rob_dat[[col_auth]]),
  Year_rob    = suppressWarnings(as.numeric(rob_dat[[col_year]])),
  rob_overall = stringr::str_squish(as.character(rob_dat[[col_rob]]))
) %>%
  dplyr::filter(!is.na(rob_overall), !rob_overall %in% c("", "NA")) %>%
  dplyr::mutate(rob_clean = dplyr::case_when(
    stringr::str_to_lower(rob_overall) == "low"  ~ "Low",
    stringr::str_to_lower(rob_overall) == "high" ~ "High",
    TRUE ~ NA_character_))
bad_rob <- unique(rob_tbl$rob_overall[is.na(rob_tbl$rob_clean)])
if (length(bad_rob) > 0)
  stop("[RoB sensitivity] Unrecognised overall-RoB values (expected Low/High): ",
       paste(shQuote(bad_rob), collapse = ", "),
       "  [detected RoB column index = ", col_rob, "].")
rob_tbl <- rob_tbl %>% dplyr::mutate(rob_overall = rob_clean) %>% dplyr::select(-rob_clean)

# replicate the TUM merge on the RoB study ids so they match df_es$ID_2.
# df_es$ID_2 and the PROBAST ids may be read as different R types; coerce BOTH
# sides of the joins below to character so a numeric/character mismatch cannot
# make the join fail or silently drop rows.
rob_tbl <- rob_tbl %>%
  dplyr::mutate(ID_2 = as.character(ifelse(ID_2_raw %in% c(11, 12, 21), 11, ID_2_raw)))

rob_no_id <- rob_tbl %>% dplyr::filter(is.na(ID_2))
if (nrow(rob_no_id) > 0)
  warning("[RoB sensitivity] ", nrow(rob_no_id), " PROBAST row(s) still have no Study ID and are ignored: ",
          paste(unique(rob_no_id$Author_rob), collapse = "; "))
rob_tbl <- rob_tbl %>% dplyr::filter(!is.na(ID_2))

# per-paper (ID_2) verdict; papers whose experiments disagree are 'mixed' (NA)
paper_rob <- rob_tbl %>%
  dplyr::group_by(ID_2) %>%
  dplyr::summarise(n_low  = sum(rob_overall == "Low"),
                   n_high = sum(rob_overall == "High"),
                   .groups = "drop") %>%
  dplyr::mutate(paper_verdict = dplyr::case_when(
    n_high == 0 & n_low  > 0 ~ "Low",
    n_low  == 0 & n_high > 0 ~ "High",
    TRUE ~ NA_character_))          # NA = mixed -> resolve per experiment below

# -----------------------------------------------------------------------------
#   Manual experiment-level RoB overrides used only for merged ID_2 clusters with 
#   mixed experiment-level PROBAST judgements.
# -----------------------------------------------------------------------------
rob_experiment_override <- tibble::tribble(
  ~Experiment_ID, ~rob_overall,
  "341", "Low",
  "342", "High",
  "343", "High",
  "344", "Low",
  "345", "Low",
  "346", "Low",
  "111", "Low",
  "121", "Low",
  "211", "High",
  "212", "High",
  "213", "High"
) %>% dplyr::mutate(Experiment_ID = as.character(Experiment_ID),
                    rob_overall   = as.character(rob_overall))

# resolve a Low/High verdict for every (ID_2, Experiment_ID) present in df_es
es_keys <- df_es %>%
  dplyr::distinct(ID_2, Experiment_ID) %>%
  dplyr::mutate(ID_2 = as.character(ID_2),
                Experiment_ID_chr = as.character(Experiment_ID)) %>%
  dplyr::left_join(paper_rob %>% dplyr::select(ID_2, paper_verdict), by = "ID_2") %>%
  dplyr::left_join(rob_experiment_override %>%
                     dplyr::rename(Experiment_ID_chr = Experiment_ID, exp_verdict = rob_overall),
                   by = "Experiment_ID_chr") %>%
  dplyr::mutate(rob_final = dplyr::coalesce(exp_verdict, paper_verdict))

# VALIDATION 1 (join completeness): every df_es paper must have a RoB verdict row
missing_paper <- setdiff(as.character(unique(df_es$ID_2)), paper_rob$ID_2)
if (length(missing_paper) > 0) {
  print(as.data.frame(df_es %>% dplyr::filter(ID_2 %in% missing_paper) %>%
                        dplyr::distinct(ID_2, Author) %>% dplyr::arrange(ID_2)))
  stop("[RoB sensitivity] No RoB verdict for the ID_2 value(s) above. Complete the Study ID ",
       "column in the PROBAST file (remember the 11/12/21 -> 11 TUM merge).")
}

# VALIDATION 2 (mixed papers): every experiment of a mixed paper needs an override
unresolved <- es_keys %>% dplyr::filter(is.na(rob_final))
if (nrow(unresolved) > 0) {
  message("[RoB sensitivity] These experiments belong to mixed-RoB papers and need an ",
          "explicit verdict in rob_experiment_override:")
  print(as.data.frame(
    df_es %>% dplyr::semi_join(unresolved, by = c("ID_2", "Experiment_ID")) %>%
      dplyr::group_by(ID_2, Experiment_ID, Author) %>%
      dplyr::summarise(n_ES = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(ID_2, Experiment_ID)))
  stop("[RoB sensitivity] Fill rob_experiment_override for the Experiment_ID(s) above, then rerun.")
}

# attach verdict to every ES row and split (ID_2 as character on both sides)
df_es_rob <- df_es %>%
  dplyr::mutate(ID_2 = as.character(ID_2)) %>%
  dplyr::left_join(es_keys %>% dplyr::select(ID_2, Experiment_ID, rob_final),
                   by = c("ID_2", "Experiment_ID"))
stopifnot(!anyNA(df_es_rob$rob_final))
df_es_low <- df_es_rob %>% dplyr::filter(rob_final == "Low")

# exclusion audit + classification export
rob_audit <- tibble::tibble(
  level        = c("effect sizes", "experiments", "paper clusters (ID_2)"),
  total        = c(nrow(df_es_rob), dplyr::n_distinct(df_es_rob$Experiment_ID),
                   dplyr::n_distinct(df_es_rob$ID_2)),
  low_rob_kept = c(nrow(df_es_low), dplyr::n_distinct(df_es_low$Experiment_ID),
                   dplyr::n_distinct(df_es_low$ID_2)))
message("[RoB sensitivity] Low-RoB subset:"); print(as.data.frame(rob_audit))
readr::write_csv(
  df_es_rob %>% dplyr::distinct(ID_2, Experiment_ID, Author, rob_final) %>%
    dplyr::arrange(rob_final, ID_2, Experiment_ID),
  file.path(output_dir, "rob_experiment_classification.csv"))

n_clusters_low <- dplyr::n_distinct(df_es_low$ID_2)
if (n_clusters_low < 5)
  warning("[RoB sensitivity] Only ", n_clusters_low, " paper clusters remain -- ",
          "CR2/Satterthwaite inference is unreliable at this cluster count.")

# -----------------------------------------------------------------------------
# Refit the IDENTICAL main model on the low-RoB subset (+ CR2 at ID_2)
# -----------------------------------------------------------------------------
lowrob_model  <- fit_mv(df_es_low)
lowrob_robust <- cr2(lowrob_model, df_es_low$ID_2)
message("[RoB sensitivity] Low-RoB pooled model (CR2 clustered at ID_2):")
print(summary(lowrob_robust))

readr::write_csv(
  tibble::tibble(
    subset        = "low_RoB_only",
    k_ES          = lowrob_robust$k,
    n_experiments = dplyr::n_distinct(df_es_low$Experiment_ID),
    n_clusters    = n_clusters_low,
    estimate      = as.numeric(coef(lowrob_robust)),
    ci_lb         = as.numeric(lowrob_robust$ci.lb),
    ci_ub         = as.numeric(lowrob_robust$ci.ub),
    tau2_total    = sum(lowrob_model$sigma2)),
  file.path(output_dir, "lowrob_main_effect.csv"))

# -----------------------------------------------------------------------------
# Forest plot for the low-RoB subset (same styling as the main forest)
# weights come from the LOW-RoB model itself, not the full-model w_pct.
# -----------------------------------------------------------------------------
local({
  w_src <- tryCatch(as.numeric(stats::weights(lowrob_model, type = "rowsum")),
                    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(w_src) || length(w_src) != nrow(df_es_low) || !all(is.finite(w_src))) {
    tau2_tot <- sum(c(lowrob_model$sigma2, lowrob_model$tau2), na.rm = TRUE)
    w_src <- 1 / (df_es_low$vi_g + tau2_tot)
    message("[RoB sensitivity] using APPROX 1/(vi+tau^2) weights")
  }
  est <- as.numeric(coef(lowrob_robust))[1]
  lo  <- as.numeric(lowrob_robust$ci.lb)[1]
  hi  <- as.numeric(lowrob_robust$ci.ub)[1]

  plot_df <- df_es_low %>%
    dplyr::left_join(forest_label_map, by = "ES_ID") %>%
    dplyr::mutate(
      effect     = hedges_g,
      se         = sqrt(vi_g),
      ci.lb      = effect - 1.96 * se,
      ci.ub      = effect + 1.96 * se,
      weight     = 100 * .env$w_src / sum(.env$w_src, na.rm = TRUE),
      weight_lab = ifelse(is.finite(weight), sprintf("%.1f%%", weight), "")
    ) %>%
    dplyr::arrange(dplyr::desc(effect), Experiment_ID, ES_ID) %>%
    dplyr::mutate(
      y    = dplyr::row_number(),
      sign = dplyr::case_when(ci.lb > 0 ~ "pos", ci.ub < 0 ~ "neg", TRUE ~ "ns"),
      col  = dplyr::case_when(sign == "pos" ~ FOREST_COL$pos,
                              sign == "neg" ~ FOREST_COL$neg, TRUE ~ FOREST_COL$ns)
    )

  lims    <- forest_xlims(plot_df, lo, hi)
  plot_df <- forest_cap(plot_df, lims$xlim_lo, lims$xlim_hi)

  summary_text <- sprintf("Low-RoB Hedges' g = %.3f [%.3f, %.3f]  (k = %d, clusters = %d)",
                          est, lo, hi, nrow(plot_df), dplyr::n_distinct(df_es_low$ID_2))

  # gold_size 2.7 and default tick length reproduce the original low-RoB figure.
  p <- draw_labeled_forest(plot_df, est, lo, hi, lims, summary_text,
                           gold_size = 2.7, theme = forest_theme(tick_length = NULL))

  H <- max(5.0, nrow(plot_df) * 0.16 + 1.6)
  forest_save_pdf(file.path(output_dir, "forest_lowRoB_labels_weights.pdf"), p,
                  width = 8.2, height = H, units = "in")
  forest_save_png(file.path(output_dir, "forest_lowRoB_labels_weights.png"), p,
                  width = 8.2, height = H, units = "in", dpi = 400)
  message("[RoB sensitivity] forest saved: forest_lowRoB_labels_weights.pdf/png  (",
          nrow(plot_df), " ES, ", dplyr::n_distinct(df_es_low$ID_2), " clusters)")
})


# =============================================================================
# Multivariable meta-regression with effect coding (global intercept + all subgroups) - Results shown in Table S3
# Reported level coefficients are adjusted deviations from the effect-coded intercept, not subgroup-specific pooled means.
# =============================================================================

# Moderator analysis: set up factors (desc vs infer)
mods <- intersect(ANALYSIS_MODERATORS, names(df_es))

# Clean moderator values: remove leading/trailing and duplicate whitespace
df_es <- df_es %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(mods),
      ~ if (is.character(.x) || is.factor(.x)) {
        stringr::str_squish(
          stringr::str_replace_all(as.character(.x), "\u00A0", " ")
        )
      } else {
        .x
      }
    )
  )

for (m in mods) df_es <- collapse_by_min_exp(df_es, m, min_exp = 2, other_label = "Other")


# 1) Prepare moderator list and modeling data
mmods <- paste0(mods, "_infer")
mmods <- mmods[mmods %in% names(df_es)]

df_mm <- df_es %>%
  dplyr::filter(is.finite(hedges_g), is.finite(vi_g),
                !is.na(Experiment_ID), !is.na(Effect_ID)) %>%
  dplyr::select(hedges_g, vi_g, ID_2, Experiment_ID, Effect_ID, dplyr::all_of(mmods)) %>%
  dplyr::mutate(across(all_of(mmods), ~ droplevels(factor(.x))))

# Drop moderators that ended up with <2 levels after the joint filter
mmods <- mmods[sapply(mmods, function(v) nlevels(df_mm[[v]]) >= 2)]
if (length(mmods) == 0) stop("No usable moderators with >=2 levels for effect-coded model.")

# 2) Apply effect coding (deviation/contr.sum) per moderator
for (v in mmods) {
  levs <- levels(df_mm[[v]])
  if (length(levs) >= 2) contrasts(df_mm[[v]]) <- contr.sum(n = length(levs))
}

# 3) Build formula with a global intercept and all moderators (DO NOT wrap in factor())
form_eff <- as.formula(paste("~", paste(mmods, collapse = " + ")))

# 4) Fit multivariable meta-regression (multilevel) with CR2 robust inference
multi_fit_eff <- fit_mv(df_mm, mods = form_eff)

## ==== Inspect which moderator blocks actually entered the model ====
beta_names <- names(coef(multi_fit_eff))
Xnames     <- colnames(multi_fit_eff$X)  # design colnames

cols_for <- function(v) grep(paste0("^`?", v, "`?\\d+$"), beta_names, value = TRUE)

block_map <- setNames(lapply(mmods, cols_for), mmods)
print(block_map)  # which columns per moderator (empty => dropped)

# keep only moderators that actually contributed at least one column
mmods_keep <- names(block_map)[lengths(block_map) > 0]
if (length(mmods_keep) == 0) stop("All moderator blocks were aliased/dropped in the joint model.")

# If we dropped any blocks, refit a lean model that includes only the estimable ones
if (!setequal(mmods_keep, mmods)) {
  message("Refitting model with estimable moderator blocks only: ", paste(mmods_keep, collapse = ", "))
  form_eff <- as.formula(paste("~", paste(mmods_keep, collapse = " + ")))
  multi_fit_eff <- fit_mv(df_mm, mods = form_eff)
  beta_names <- names(coef(multi_fit_eff))
  Xnames     <- colnames(multi_fit_eff$X)
}

# CR2 vcov
Vcr2_eff <- clubSandwich::vcovCR(multi_fit_eff, type = "CR2", cluster = df_mm$ID_2)

# Intercept (grand mean)
ct <- clubSandwich::coef_test(multi_fit_eff, vcov = Vcr2_eff, test = "Satterthwaite")
df_col <- pick_col(ct, c("df_Satt","d.f.","df"))
row_int <- rownames(ct) %in% c("intrcpt","(Intercept)")
grand_mean_tbl <- tibble::tibble(
  Moderator="(Global)", Level="(Intercept / Grand mean)",
  estimate = ct$beta[row_int],
  SE       = ct$SE[row_int],
  df       = ct[[df_col]][row_int],
  t        = ct$tstat[row_int],
  p        = ct[[pick_col(ct, c("p_Satt","p_val","p_value","p"))]][row_int],
  conf.low = ct$beta[row_int] - qt(0.975, df = ct[[df_col]][row_int]) * ct$SE[row_int],
  conf.high= ct$beta[row_int] + qt(0.975, df = ct[[df_col]][row_int]) * ct$SE[row_int]
)

## ==== Reconstruct ALL subgroup deviations per moderator (incl. the omitted level) ====
b   <- coef(multi_fit_eff)
Vb  <- as.matrix(Vcr2_eff)

make_all_levels_table <- function(v) {
  levs <- levels(df_mm[[v]])
  k    <- length(levs)
  # columns for this block in the fitted model
  cols <- cols_for(v)
  if (length(cols) == 0) return(tibble::tibble())  # block not estimable

  # Build L: first k-1 pickers, last row = negative sum (effect coding)
  p <- length(b)
  L <- matrix(0, nrow = k, ncol = p, dimnames = list(levs, names(b)))
  if (length(cols) >= 1) {
    for (j in seq_len(min(k-1, length(cols)))) L[j, cols[j]] <- 1
    if (k >= 2) L[k, cols] <- -1
  }

  est <- as.numeric(L %*% b)
  Var <- L %*% Vb %*% t(L)
  se  <- sqrt(pmax(diag(Var), 0))

  # Satterthwaite denominator df PER LEVEL: one single-row constraint per level.
  # A one-row constraint is always full rank, so the test cannot fail because the
  # block matrix L (which contains the reconstructed omitted level) is rank
  # deficient. A block-level df rolled out across levels would not be the
  # Satterthwaite df of the individual contrast.
  dfS <- vapply(seq_len(k), function(j) {
    wt1 <- try(
      clubSandwich::Wald_test(multi_fit_eff,
                              constraints = L[j, , drop = FALSE],
                              vcov = Vcr2_eff, test = "HTZ"),
      silent = TRUE
    )
    if (inherits(wt1, "try-error")) NA_real_ else as.numeric(wt1$df_denom)[1]
  }, numeric(1))

  # t / p / CI with dfS if available, else normal
  tval <- est / se
  pval <- ifelse(is.finite(dfS),
                 2*pt(abs(tval), df = dfS, lower.tail = FALSE),
                 2*pnorm(abs(tval), lower.tail = FALSE))
  crit <- ifelse(is.finite(dfS), qt(0.975, df = dfS), qnorm(0.975))
  lo   <- est - crit * se
  hi   <- est + crit * se

  tibble::tibble(
    Moderator = v, Level = levs,
    estimate = est, SE = se, df = dfS, t = tval, p = pval,
    conf.low = lo, conf.high = hi
  )
}

all_levels <- purrr::map_dfr(mmods_keep, make_all_levels_table)

## ==== Per-moderator omnibus (only blocks that exist). If CR2 fails, show NA rather than error ====
omni_eff <- purrr::map_dfr(mmods_keep, function(v) {
  cols <- cols_for(v)
  if (length(cols) == 0) {
    return(tibble::tibble(Moderator=v, df1=NA_real_, df2=NA_real_, F=NA_real_, p=NA_real_,
                          k=multi_fit_eff$k, clusters=dplyr::n_distinct(df_mm$ID_2)))
  }
  R <- matrix(0, nrow = length(cols), ncol = length(b), dimnames = list(cols, names(b)))
  for (i in seq_along(cols)) R[i, cols[i]] <- 1
  wt <- try(clubSandwich::Wald_test(multi_fit_eff, constraints = R, vcov = Vcr2_eff, test = "HTZ"), silent = TRUE)
  if (inherits(wt, "try-error")) {
    tibble::tibble(Moderator=v, df1=NA_real_, df2=NA_real_, F=NA_real_, p=NA_real_,
                   k=multi_fit_eff$k, clusters=dplyr::n_distinct(df_mm$ID_2))
  } else {
    tibble::tibble(
      Moderator=v, df1=unname(wt$df_num), df2=unname(wt$df_denom),
      F=unname(wt$Fstat), p=unname(wt$p_val),
      k=multi_fit_eff$k, clusters=dplyr::n_distinct(df_mm$ID_2)
    )
  }
})

## ==== Final table and write ====
final_tbl <- dplyr::bind_rows(grand_mean_tbl, all_levels) |>
  dplyr::group_by(Moderator) |>
  dplyr::arrange(Moderator, estimate, .by_group = TRUE) |>
  dplyr::ungroup()

readr::write_csv(final_tbl, file.path(output_dir, "meta_regression_effect_coded_ALL_LEVELS_CR2.csv"))
readr::write_csv(omni_eff,  file.path(output_dir, "meta_regression_effect_coded_OMNIBUS_HTZ.csv"))

message("Blocks in model:"); print(block_map)
message("Kept blocks: ", paste(mmods_keep, collapse = ", "))


# =============================================================================
# Moderator analyses - Results shown in Figure 2 (main text) and Table S2 (Appendix)
# =============================================================================


# NO-INTERCEPT: per-level subgroup means (CR2 where possible) to determine the subgroup means -----------------
means_out <- dplyr::tibble()

for (mod in mods) {
  message("\n[NO-INT/DESC] Moderator: ", mod)
  fac_desc <- paste0(mod, "_desc")

  df_tmp <- df_es %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(.data[[fac_desc]])) %>%
    dplyr::mutate(.mod_fac = droplevels(factor(.data[[fac_desc]])))

  if (nrow(df_tmp) == 0) { message("  -> no rows after filtering."); next }

  lev_counts <- df_tmp %>% dplyr::count(.mod_fac, name = "k_effects")
  message("  levels: ", nlevels(df_tmp$.mod_fac), " | ",
          paste(paste0(as.character(lev_counts$.mod_fac), "=", lev_counts$k_effects), collapse = ", "))

  nlev <- nlevels(df_tmp$.mod_fac)

  if (nlev == 1) {
    # Single level: intercept-only fits that subgroup mean
    lev_name <- levels(df_tmp$.mod_fac)[1]

    fit1 <- fit_mv(df_tmp)

    per_level <- try({
      Vcr2 <- clubSandwich::vcovCR(fit1, type = "CR2", cluster = df_tmp$ID_2)
      out  <- clubSandwich::coef_test(fit1, vcov = Vcr2, test = "Satterthwaite") %>% as.data.frame()
      df_col <- pick_col(out, c("df_Satt","d.f.","df"))
      out$conf.low  <- out$beta - qt(0.975, df = out[[df_col]]) * out$SE
      out$conf.high <- out$beta + qt(0.975, df = out[[df_col]]) * out$SE
      out$Level <- lev_name
      out$inference <- "CR2 (Satterthwaite)"
      out
    }, silent = TRUE)

    if (inherits(per_level, "try-error")) {
      beta  <- as.numeric(coef(fit1))
      se    <- sqrt(diag(vcov(fit1)))[1]
      zcrit <- qnorm(0.975)
      per_level <- data.frame(
        beta      = beta,
        SE        = se,
        df_Satt   = NA_real_,
        tstat     = beta / se,
        p_Satt    = NA_real_,
        conf.low  = beta - zcrit * se,
        conf.high = beta + zcrit * se,
        Level     = lev_name,
        inference = "REML (descriptive)",
        check.names = FALSE
      )
    }

    k_eff <- nrow(df_tmp)
    k_exp <- dplyr::n_distinct(df_tmp$Experiment_ID)

    means_tab <- tibble::tibble(
      Moderator   = mod,
      Level       = per_level$Level,
      estimate    = per_level$beta,
      std.error   = per_level$SE,
      df          = per_level[[pick_col(per_level, c("df_Satt","d.f.","df"))]],
      t.value     = per_level$tstat,
      p.value     = per_level[[pick_col(per_level, c("p_Satt","p_val","p_value","p"))]],
      conf.low    = per_level$conf.low,
      conf.high   = per_level$conf.high,
      inference   = per_level$inference,
      k_effects   = k_eff,
      k_experiments = k_exp
    )
    means_out <- dplyr::bind_rows(means_out, means_tab)
    next
  }

  # Multiple levels: no-intercept model gives each level mean
  fit <- fit_mv(df_tmp, mods = ~ 0 + .mod_fac)

  per_level <- try({
    Vcr2 <- clubSandwich::vcovCR(fit, type = "CR2", cluster = df_tmp$ID_2)
    out  <- clubSandwich::coef_test(fit, vcov = Vcr2, test = "Satterthwaite") %>% as.data.frame()
    out$term <- rownames(out)
    df_col <- pick_col(out, c("df_Satt","d.f.","df"))
    out$conf.low  <- out$beta - qt(0.975, df = out[[df_col]]) * out$SE
    out$conf.high <- out$beta + qt(0.975, df = out[[df_col]]) * out$SE
    out$inference <- "CR2 (Satterthwaite)"
    out
  }, silent = TRUE)

  if (inherits(per_level, "try-error")) {
    beta  <- coef(fit); V <- vcov(fit); se <- sqrt(diag(V)); zcrit <- qnorm(0.975)
    per_level <- data.frame(
      beta      = as.numeric(beta),
      SE        = as.numeric(se),
      term      = names(beta),
      df_Satt   = NA_real_,
      tstat     = as.numeric(beta) / as.numeric(se),
      p_Satt    = NA_real_,
      conf.low  = as.numeric(beta) - zcrit * as.numeric(se),
      conf.high = as.numeric(beta) + zcrit * as.numeric(se),
      inference = "REML (descriptive)",
      check.names = FALSE
    )
  }

  # Clean term → level names
  lvl <- per_level$term
  lvl <- gsub("`", "", lvl)
  lvl <- sub("^.*\\.mod_fac", "", lvl)
  lvl <- sub("^mod_fac", "", lvl)
  per_level$Level <- lvl

  # Counts per level
  k_eff <- df_tmp %>% dplyr::count(.mod_fac, name = "k_effects") %>%
    dplyr::mutate(Level = as.character(.mod_fac)) %>% dplyr::select(Level, k_effects)
  k_exp <- df_tmp %>% dplyr::group_by(.mod_fac) %>%
    dplyr::summarise(k_experiments = dplyr::n_distinct(Experiment_ID), .groups = "drop") %>%
    dplyr::mutate(Level = as.character(.mod_fac)) %>% dplyr::select(Level, k_experiments)

  means_tab <- dplyr::tibble(
    Moderator  = mod,
    Level      = per_level$Level,
    estimate   = per_level$beta,
    std.error  = per_level$SE,
    df         = per_level[[pick_col(per_level, c("df_Satt","d.f.","df"))]],
    t.value    = per_level[[pick_col(per_level, c("tstat","t","Tstat"))]],
    p.value    = per_level[[pick_col(per_level, c("p_Satt","p_val","p_value","p"))]],
    conf.low   = per_level$conf.low,
    conf.high  = per_level$conf.high,
    inference  = per_level$inference
  ) %>%
    dplyr::left_join(k_eff, by = "Level") %>%
    dplyr::left_join(k_exp, by = "Level")

  means_out <- dplyr::bind_rows(means_out, means_tab)
}

print(means_out %>% dplyr::group_by(Moderator) %>% dplyr::summarise(n_levels = dplyr::n(), .groups = "drop"))
readr::write_csv(means_out, file.path(output_dir, "moderator_NOINT_MEANS_CR2_OR_REML.csv"))

# -----------------------------------------------------------------------------
# Helpers shared by the collapsed-level (_infer) omnibus tests below.
#   .moderator_frame_infer(): modelling frame for a collapsed moderator, keeping
#     only levels observed in >= min_clusters paper clusters (ID_2).
#   .htz_cols(): the (df1, df2, F, p) columns from a clubSandwich Wald_test.
# The no-intercept HTZ test below assesses equality of the collapsed level means.
.moderator_frame_infer <- function(mod, min_clusters) {
  fac_infer <- paste0(mod, "_infer")
  df_tmp <- df_es %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(.data[[fac_infer]])) %>%
    dplyr::mutate(.mod_fac = droplevels(factor(.data[[fac_infer]])))
  keep_levels <- df_tmp %>%
    dplyr::group_by(.mod_fac) %>%
    dplyr::summarise(n_exp = dplyr::n_distinct(ID_2), .groups = "drop") %>%
    dplyr::filter(n_exp >= min_clusters) %>%
    dplyr::pull(.mod_fac) %>% as.character()
  df_tmp %>%
    dplyr::filter(as.character(.mod_fac) %in% keep_levels) %>%
    dplyr::mutate(.mod_fac = droplevels(.mod_fac))
}

.htz_cols <- function(wt) {
  dplyr::tibble(df1 = unname(wt$df_num), df2 = unname(wt$df_denom),
                F = unname(wt$Fstat), p = unname(wt$p_val))
}

# NO-INTERCEPT omnibus: equality of means across levels (collapsed), here we test if all sub-group means are equal (Appendix Table S2) -----------
omni_noint <- dplyr::tibble()
min_clusters_per_level_infer <- 2

for (mod in mods) {
  message("\n[NO-INT/INFER OMNIBUS] Moderator: ", mod)

  df_tmp <- .moderator_frame_infer(mod, min_clusters_per_level_infer)

  if (nlevels(df_tmp$.mod_fac) < 2) { message("  -> Skip omnibus (no-intercept): too few analyzable levels after sparsity filter."); next }

  df_tmp$.mod_fac <- factor(df_tmp$.mod_fac, levels = sort(levels(df_tmp$.mod_fac)))
  levs <- levels(df_tmp$.mod_fac)

  fit_noint <- fit_mv(df_tmp, mods = ~ 0 + .mod_fac)

  Vcr2_noint <- clubSandwich::vcovCR(fit_noint, type = "CR2", cluster = df_tmp$ID_2)

  beta_names <- names(coef(fit_noint))
  vc_names   <- colnames(Vcr2_noint); if (is.null(vc_names)) vc_names <- beta_names

  level_coefs <- setNames(rep(NA_character_, length(levs)), levs)
  for (lv in levs) {
    plain <- paste0(".mod_fac", lv)
    tick  <- paste0("`.mod_fac", lv, "`")
    found <- intersect(vc_names, c(plain, tick))
    if (length(found) == 1) level_coefs[lv] <- found
  }
  if (any(is.na(level_coefs))) { message("  -> Skip omnibus (no-intercept): could not map all level coefficients."); next }

  make_R_row <- function(A, B, vc_names, level_coefs) {
    r <- rep(0, length(vc_names)); names(r) <- vc_names
    r[level_coefs[A]] <-  1
    r[level_coefs[B]] <- -1
    matrix(r, nrow = 1)
  }

  R_list <- lapply(levs[-1], function(lv) make_R_row(levs[1], lv, vc_names, level_coefs))
  R <- do.call(rbind, R_list)
  rownames(R) <- paste0(levs[1], " - ", levs[-1])
  colnames(R) <- vc_names

  omni_ni <- tryCatch(
    clubSandwich::Wald_test(fit_noint, constraints = R, vcov = Vcr2_noint, test = "HTZ"),
    error = function(e) { message("  -> equality-of-means omnibus not estimable (CR2 vcov not PD): ", mod); NULL }
  )
  if (is.null(omni_ni)) next

  omni_noint <- dplyr::bind_rows(
    omni_noint,
    dplyr::bind_cols(
      dplyr::tibble(
        Moderator = mod,
        test      = "CR2 HTZ (equality of level means) on COLLAPSED levels (no-intercept)"
      ),
      .htz_cols(omni_ni),
      dplyr::tibble(levels = paste(levs, collapse = "|"))
    )
  )
}

readr::write_csv(omni_noint, file.path(output_dir, "moderator_NOINTERCEPT_OMNIBUS_CR2_EQUAL_MEANS.csv"))


# =============================================================================
#Robustness and sensitivity checks
# =============================================================================

# Funnel plot
fig_name <- file.path(output_dir, "funnel_main_model.png")

estimate <- as.numeric(coef(main_model_robust))  # robust center line

se_vals <- sqrt(df_es$vi_g)
ylim_se <- c(max(se_vals, na.rm = TRUE) * 1.05, 0)

title    <- "Funnel Plot (Hedges' g)"
#subtitle <- paste0("k = ", main_model$k)

if (!exists(".dev_w", inherits = FALSE))  .dev_w <- 1000L
if (!exists(".dev_h", inherits = FALSE))  .dev_h <- 1000L
if (!exists(".dev_res", inherits = FALSE)) .dev_res <- 150L
png(fig_name, width = .dev_w, height = .dev_h, res = .dev_res)

metafor::funnel(
  x       = main_model,
  yaxis   = "sei",
   xlab    = "Effect size (Hedges' g)",
  ylab    = "Standard Error",
  ylim    = ylim_se,
  steps   = 6,
  digits  = c(1, 2),
  level   = 95,
  shade   = c("white", "gray85"),
  refline = estimate,
  main    = title,
  pch     = 19,
  col     = "#0072B2"
)
dev.off()
message("Saved funnel plot: ", fig_name)

# Egger / PET regression (multilevel, CR2-robust)
# PRIMARY: the precision proxy is built from sample size only. For a standardized
# mean difference the sampling variance is itself a function of the effect size,
# so regressing g on sqrt(vi) induces an artefactual association and inflates the
# false-positive rate of the asymmetry test (Pustejovsky & Rodgers 2019).
# vi_g additionally carries the reader design effect, which is a design feature
# rather than a measure of precision.
df_es <- df_es %>%
  dplyr::mutate(
    n_unit_egger  = dplyr::if_else(.paired,
                                   as_num(.n_pairs),
                                   as_num(n_treatment) + as_num(n_control)),
    n_eff_egger   = n_unit_egger / .DE,
    se_size_egger = sqrt(1 / n_eff_egger)
  )
if (any(!is.finite(df_es$se_size_egger)))
  stop("[Egger] se_size_egger is not finite for ",
       sum(!is.finite(df_es$se_size_egger)),
       " row(s); check .n_pairs / n_treatment / n_control.")

eggers_model        <- fit_mv(df_es, mods = ~ se_size_egger)
eggers_model_robust <- cr2(eggers_model, df_es$ID_2)
print(summary(eggers_model_robust))

eg_coefs <- coef(eggers_model_robust)
eg_se    <- sqrt(diag(vcov(eggers_model_robust)))
eg_stat  <- eg_coefs / eg_se
eg_df    <- eggers_model_robust$ddf
if (!is.null(eg_df) && all(is.finite(eg_df) & eg_df > 0)) {
  eg_p <- 2 * pt(abs(eg_stat), df = eg_df, lower.tail = FALSE)
} else {
  eg_p <- 2 * pnorm(abs(eg_stat), lower.tail = FALSE)
}

# The INTERCEPT is the small-study-adjusted (PET) effect estimate -- a substantive
# robustness result. The SLOPE is the asymmetry parameter that the test refers to.
egger_tab <- data.frame(
  term     = c("intercept (PET-adjusted g)", "slope (small-study asymmetry)"),
  estimate = as.numeric(eg_coefs),
  se       = as.numeric(eg_se),
  stat     = as.numeric(eg_stat),
  df       = as.numeric(rep(eg_df, length.out = length(eg_coefs))),
  p        = as.numeric(eg_p),
  row.names = NULL
)
cat("\nEgger / PET (size-based predictor, CR2-robust):\n"); print(egger_tab)
readr::write_csv(egger_tab, file.path(output_dir, "egger_pet_size_predictor.csv"))

# SECONDARY: the SE-based predictor, retained for comparability with the previous
# version of the analysis. Reported as a supplementary check only.
eggers_model_se    <- fit_mv(df_es, mods = ~ sqrt(vi_g))
eggers_model_se_rb <- cr2(eggers_model_se, df_es$ID_2)
cat("\nEgger (SE-based predictor, secondary check):\n")
print(summary(eggers_model_se_rb))

# =============================================================================
# Sensitivity analysis: four-level working model
#   The primary model is a conventional three-level model:
#     level 1 = sampling error,
#     level 2 = effect sizes,
#     level 3 = experiments.
#
#   This sensitivity model adds a paper/shared-case-corpus random intercept
#   as level 4: ~ 1 | ID_2/Experiment_ID/Effect_ID.
#   CR2 inference remains clustered at ID_2;
#   It only asks whether adding an explicit between-paper variance
#   component to the working model moves the pooled estimate.
# =============================================================================
stopifnot("ID_2" %in% names(df_es))
stopifnot(all(c("hedges_g","vi_g","Experiment_ID","Effect_ID") %in% names(df_es)))

df_paper <- df_es %>%
  dplyr::filter(is.finite(hedges_g), is.finite(vi_g),
                !is.na(ID_2), !is.na(Experiment_ID), !is.na(Effect_ID))

if (nrow(df_paper) < 2) stop("Not enough rows after filtering.")
if (dplyr::n_distinct(df_paper$ID_2) < 2) stop("Need at least 2 paper/case-corpus clusters (ID_2) for CR2-robust inference.")

main_model_4lvl    <- fit_mv(df_paper, random = ~ 1 | ID_2/Experiment_ID/Effect_ID)
main_model_4lvl_rb <- cr2(main_model_4lvl, df_paper$ID_2)

print(summary(main_model_4lvl_rb))

# Variance components (outer -> inner): paper/case-corpus, experiment, effect size
s2           <- main_model_4lvl$sigma2
sigma2_paper <- if (length(s2) >= 1) s2[1] else NA_real_
sigma2_exp   <- if (length(s2) >= 2) s2[2] else NA_real_
sigma2_es    <- if (length(s2) >= 3) s2[3] else NA_real_

rb4 <- .extract_rb(main_model_4lvl_rb)

sens_summary <- tibble::tibble(
  clusters     = dplyr::n_distinct(df_paper$ID_2),
  k            = nrow(df_paper),
  estimate     = rb4$beta,
  se           = rb4$se,
  stat         = rb4$stat,
  df           = rb4$df,
  p            = rb4$p,
  ci.lb        = rb4$lo,
  ci.ub        = rb4$hi,
  sigma2_paper = sigma2_paper,
  sigma2_exp   = sigma2_exp,
  sigma2_es    = sigma2_es
)

print(sens_summary)
readr::write_csv(sens_summary, file.path(output_dir, "sensitivity_4level_working_model_summary.csv"))

# =============================================================================
# Leave-one-out Analysis
# =============================================================================

# Leave-one-EFFECT-out. The working model and CR2 clustering are IDENTICAL to the
# main model: random = ~ 1 | Experiment_ID/Effect_ID with CR2 clustered at the
# paper/case-corpus level (ID_2). Only the unit that is dropped changes here (one
# ES_ID at a time); the inferential structure is unchanged.
stopifnot(all(c("hedges_g","vi_g","Experiment_ID","Effect_ID","ES_ID") %in% names(df_es)))

df_use <- df_es %>% dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(Experiment_ID), !is.na(Effect_ID), !is.na(ES_ID))
if (nrow(df_use) < 3 || dplyr::n_distinct(df_use$Experiment_ID) < 2) stop("Need ≥3 effects and ≥2 experiments after filtering for LOO.")

loo_ids <- unique(df_use$ES_ID)

df_loo_es <- tibble::tibble(
  left_out_ES_ID = loo_ids,
  k              = NA_integer_,
  clusters       = NA_integer_,
  estimate       = NA_real_,
  se             = NA_real_,
  stat           = NA_real_,
  df             = NA_real_,
  p_val          = NA_real_,
  ci.lb          = NA_real_,
  ci.ub          = NA_real_,
  I2_total       = NA_real_
)

for (i in seq_along(loo_ids)) {
  es_id <- loo_ids[i]
  dat   <- dplyr::filter(df_use, ES_ID != es_id)
  if (nrow(dat) < 2 || dplyr::n_distinct(dat$Experiment_ID) < 2) next

  loo_model <- fit_mv(dat)

  loo_rb <- tryCatch(cr2(loo_model, dat$ID_2), error = function(e) NULL)

  ex       <- .extract_rb(loo_rb)
  I2_total <- .i2_total(loo_model, dat$vi_g)

  df_loo_es$k[i]        <- loo_model$k
  df_loo_es$clusters[i]  <- dplyr::n_distinct(dat$ID_2)
  df_loo_es$estimate[i] <- ex$beta
  df_loo_es$se[i]       <- ex$se
  df_loo_es$stat[i]     <- ex$stat
  df_loo_es$df[i]       <- ex$df
  df_loo_es$p_val[i]    <- ex$p
  df_loo_es$ci.lb[i]    <- ex$lo
  df_loo_es$ci.ub[i]    <- ex$hi
  df_loo_es$I2_total[i] <- I2_total
}

print(range(df_loo_es$estimate,   na.rm = TRUE))
print(range(df_loo_es$p_val,      na.rm = TRUE))
print(range(df_loo_es$I2_total,   na.rm = TRUE))

readr::write_csv(df_loo_es, file.path(output_dir, "loo_effect_level_experiment_cluster.csv"))

# LOO helpers -----------------------------------------------------------------
stopifnot(all(c("hedges_g","vi_g","Experiment_ID","Effect_ID","ES_ID","ID_2") %in% names(df_es)))

df_use <- df_es %>% dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(Experiment_ID), !is.na(Effect_ID), !is.na(ID_2))
if (nrow(df_use) < 3) stop("Need ≥3 effects after filtering for LOO.")

# LOO at Experiment level -----------------------------------------------------
exp_ids <- sort(unique(df_use$Experiment_ID))

df_loo_exp <- tibble::tibble(
  left_out_Experiment_ID = exp_ids,
  k              = NA_integer_,
  clusters       = NA_integer_,
  estimate       = NA_real_,
  se             = NA_real_,
  stat           = NA_real_,
  df             = NA_real_,
  p_val          = NA_real_,
  ci.lb          = NA_real_,
  ci.ub          = NA_real_,
  I2_total       = NA_real_
)

for (i in seq_along(exp_ids)) {
  eid  <- exp_ids[i]
  dat  <- dplyr::filter(df_use, Experiment_ID != eid)
  if (nrow(dat) < 2 || dplyr::n_distinct(dat$Experiment_ID) < 2) next

  fit <- fit_mv(dat)

  rb  <- tryCatch(cr2(fit, dat$ID_2), error = function(e) NULL)

  inf <- .extract_rb(rb)
  I2t <- .i2_total(fit, dat$vi_g)

  df_loo_exp$k[i]        <- fit$k
  df_loo_exp$clusters[i] <- dplyr::n_distinct(dat$ID_2)
  df_loo_exp$estimate[i] <- inf$beta
  df_loo_exp$se[i]       <- inf$se
  df_loo_exp$stat[i]     <- inf$stat
  df_loo_exp$df[i]       <- inf$df
  df_loo_exp$p_val[i]    <- inf$p
  df_loo_exp$ci.lb[i]    <- inf$lo
  df_loo_exp$ci.ub[i]    <- inf$hi
  df_loo_exp$I2_total[i] <- I2t
}

print(range(df_loo_exp$estimate, na.rm = TRUE))
print(range(df_loo_exp$p_val,    na.rm = TRUE))
print(range(df_loo_exp$I2_total, na.rm = TRUE))

readr::write_csv(df_loo_exp, file.path(output_dir, "loo_experiment_level_cluster_Experiment_ID.csv"))

# LOO at Study level ----------------------------------------------------------
# Iterate over the DISTINCT Study_ID (42) so a single study can be dropped,
# including one TUM brain-MRI paper at a time. The CR2 cluster stays ID_2 (merged),
# so the remaining TUM studies are still treated as one shared-corpus cluster.
paper_ids <- sort(unique(df_use$Study_ID))

df_loo_paper <- tibble::tibble(
  left_out_Study_ID = paper_ids,
  k              = NA_integer_,
  clusters       = NA_integer_,
  estimate       = NA_real_,
  se             = NA_real_,
  stat           = NA_real_,
  df             = NA_real_,
  p_val          = NA_real_,
  ci.lb          = NA_real_,
  ci.ub          = NA_real_,
  I2_total       = NA_real_
)

for (i in seq_along(paper_ids)) {
  pid  <- paper_ids[i]
  dat  <- dplyr::filter(df_use, Study_ID != pid)
  if (nrow(dat) < 2 || dplyr::n_distinct(dat$ID_2) < 2) next

  fit <- fit_mv(dat)

  rb  <- tryCatch(cr2(fit, dat$ID_2), error = function(e) NULL)

  inf <- .extract_rb(rb)
  I2t <- .i2_total(fit, dat$vi_g)

  df_loo_paper$k[i]        <- fit$k
  df_loo_paper$clusters[i] <- dplyr::n_distinct(dat$ID_2)
  df_loo_paper$estimate[i] <- inf$beta
  df_loo_paper$se[i]       <- inf$se
  df_loo_paper$stat[i]     <- inf$stat
  df_loo_paper$df[i]       <- inf$df
  df_loo_paper$p_val[i]    <- inf$p
  df_loo_paper$ci.lb[i]    <- inf$lo
  df_loo_paper$ci.ub[i]    <- inf$hi
  df_loo_paper$I2_total[i] <- I2t
}

print(range(df_loo_paper$estimate, na.rm = TRUE))
print(range(df_loo_paper$p_val,    na.rm = TRUE))
print(range(df_loo_paper$I2_total, na.rm = TRUE))

readr::write_csv(df_loo_paper, file.path(output_dir, "loo_study_level_leave_one_study_out.csv"))

# =============================================================================
# === Format meta-regression results (ALL_LEVELS_CR2_TABLE) for LaTeX ===
# =============================================================================

# Load data
meta_tbl <- read_csv(
  file.path(output_dir, "meta_regression_effect_coded_ALL_LEVELS_CR2.csv"),
  show_col_types = FALSE
)


# ---- Step 1: Keep p for significance stars, drop other columns we don't report ----
meta_tbl <- meta_tbl %>%
  dplyr::select(-dplyr::any_of(c("t", "df", "test")))

# ---- Step 2: Rename columns ----
meta_tbl <- meta_tbl %>%
  dplyr::rename(
    Subgroup = Level,
    Beta_raw = estimate,
    SE_raw = SE,
    CI_low = conf.low,
    CI_high = conf.high,
    p_value = p
  )

# ---- Step 3: Add significance stars from p_value ----
add_stars <- function(p) {
  dplyr::case_when(
    is.na(p)        ~ "",
    p < 0.001       ~ "***",
    p < 0.01        ~ "**",
    p < 0.05        ~ "*",
    TRUE            ~ ""
  )
}

meta_tbl <- meta_tbl %>%
  dplyr::mutate(stars = add_stars(p_value))

# ---- Step 4: Format numeric cols and append stars to Beta ----
meta_tbl <- meta_tbl %>%
  dplyr::mutate(
    Beta = sprintf("%.3f%s", Beta_raw, stars),
    SE   = sprintf("%.3f", SE_raw),
    CI   = sprintf("[%.2f, %.2f]", CI_low, CI_high)
  )

# ---- Step 5: Capitalize first letter in Subgroup ----
meta_tbl <- meta_tbl %>%
  dplyr::mutate(Subgroup = stringr::str_to_sentence(Subgroup))

# ---- Step 6: Rename moderator labels to nicer display names ----
rename_map <- c(
  "LLM_Model_infer" = "LLM Model",
  "Medical_Field_infer" = "Medical field",
  "Career_stage_infer" = "Career stage",
  "Explanation_Format_of_LLM_infer" = "Response format",
  "Task_Format_infer" = "Task format",
  "Input_Form_Patient_Information_in_LLM_infer" = "Input form",
  "Follow_Up_Prompts_Allowed_infer" = "Customized prompts",
  "Type_of_Results_infer" = "Accuracy measure",
  "Case_setting_infer" = "Setting of patient cases",
  "(Global)" = "NA"  # in case your intercept row had "(Global)"
)

meta_tbl <- meta_tbl %>%
  dplyr::mutate(
    Moderator = dplyr::recode(Moderator, !!!rename_map, .default = Moderator)
  )

# ---- Step 7: Sort moderators in required order ----
desired_order <- c(
  "LLM Model",
  "Medical field",
  "Career stage",
  "Response format",
  "Task format",
  "Input form",
  "Customized prompts",
  "Accuracy measure",
  "Setting of patient cases",
  "NA"  # ensure intercept is last
)

meta_tbl$Moderator <- factor(meta_tbl$Moderator, levels = desired_order)
meta_tbl <- meta_tbl %>%
  arrange(Moderator, Subgroup)

# ---- Step 8: Keep only final display columns ----
meta_tbl <- meta_tbl %>%
  dplyr::select(
    Moderator,
    Subgroup,
    Beta,
    SE,
    CI
  ) %>%
  dplyr::rename(`95\\% CI` = CI)

# ---- Step 9: Build LaTeX with left alignment and single header ----
latex_tbl <- meta_tbl %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Meta-regression results across all moderators (CR2 robust SEs)",
    label = "tab:meta_regression_all_levels",
    # force all columns left-aligned
    align = c("l","l","l","l","l"),
    col.names = c("Moderator", "Subgroup", "Beta", "SE", "95\\% CI")
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 9
  ) %>%
  row_spec(0, bold = TRUE)

# ---- Step 10: reduce row spacing in table ----
latex_tbl <- paste0(
  "% smaller row spacing\n",
  "\\renewcommand{\\arraystretch}{0.8}\n",
  latex_tbl
)

# ---- Step 11: write to .tex ----
latex_path <- file.path(output_dir, "meta_regression_effect_coded_ALL_LEVELS_CR2_TABLE.tex")
cat(latex_tbl, file = latex_path)

message("Clean LaTeX table saved to: ", latex_path)