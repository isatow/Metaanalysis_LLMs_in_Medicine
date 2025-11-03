# Prepare workspace 
rm(list = ls())

# Packages
packages <- c(
  "metafor", "tidyverse", "ggtext", "gridExtra", "broom", "knitr", "xtable",
  "rstudioapi", "readr", "rlang", "dplyr", "kableExtra"
)

install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
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

mods_dir <- file.path(project_root, "moderators")
if (!dir.exists(mods_dir)) dir.create(mods_dir, recursive = TRUE)

# Load data 
df_performance <- read_delim(
  file.path(data_dir, "Coding_Observations.csv"),
  delim = ";",
  quote = "\"",
  escape_double = TRUE,
  locale = locale(decimal_mark = ".", grouping_mark = ""),
  na = c("NR", "n.d.", "", "NA"),
  show_col_types = FALSE
)

# Clean names / standard IDs 
names(df_performance) <- trimws(names(df_performance))

df_performance <- df_performance %>%
  dplyr::rename(
    Experiment_ID = Exp_ID,
    Effect_ID     = ES_ID
  )

# Helpers
library(stringr)

norm <- function(x) {
  str_to_lower(str_squish(str_replace_all(as.character(x), "\u00A0", " ")))
}

# Cast numeric columns (only if present)
numeric_cols <- intersect(
  c(
    "M_treatment","SD_treatment","SE_treatment",
    "M_control","SD_control","SE_control",
    "n_treatment","n_control","n_total",
    "Performance_Difference", "95CI_lower_performance_diff", "95CI_upper_performance_diff",
    "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment",
    "F_value","Chi_square","unstandardised_beta",
    "SE_beta","Z_value","t_value","Odds_Ratio"
  ),
  names(df_performance)
)

df_performance <- df_performance %>% mutate(across(all_of(numeric_cols), as.numeric))

# Quick structure checks 
dplyr::glimpse(df_performance[, c(
  "M_treatment","SD_treatment","SE_treatment",
  "M_control","SD_control","SE_control",
  "n_treatment","n_control","Odds_Ratio",
  "Performance_Difference", "95CI_lower_performance_diff", "95CI_upper_performance_diff",
  "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment"
)])

print(names(df_performance))

dplyr::glimpse(df_performance[, c(
  "M_treatment","SD_treatment","SE_treatment",
  "M_control","SD_control","SE_control",
  "n_treatment","n_control","Odds_Ratio",
  "Performance_Difference", "95CI_lower_performance_diff", "95CI_upper_performance_diff",
  "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment"
)])

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

# Effect-size computation 
# Compute Cohen's d from whatever stats are available; keep existing d if present
compute_cohens_d <- function(df) {
  all.possible <- c(
    "M_treatment","SD_treatment","SE_treatment",
    "M_control","SD_control","SE_control",
    "n_treatment","n_control",
    "F_value","Chi_square","unstandardised_beta",
    "SE_beta","Z_value","n_total",
    "cohens_d","t_value", "Odds_Ratio",
    "Performance_Difference", "95CI_lower_performance_diff", "95CI_upper_performance_diff",
    "95CI_lower_control", "95CI_upper_control", "95CI_lower_treatment", "95CI_upper_treatment"
  )

  # Ensure referenced columns exist
  missing.cols <- setdiff(all.possible, names(df))
  for (col in missing.cols) df[[col]] <- NA_real_

  df <- df %>%
    mutate(
      cohens_d = if_else(
        !is.na(cohens_d),
        cohens_d,
        case_when(
          # Means & SDs (Likert/Accuracy only)
          !is.na(M_treatment) & !is.na(M_control) &
            !is.na(SD_treatment) & !is.na(SD_control) &
            Measure_type %in% c("Likert", "Accuracy") ~
              (M_treatment - M_control) /
                sqrt(((n_treatment - 1)*SD_treatment^2 +
                      (n_control  - 1)*SD_control^2) /
                     (n_treatment + n_control - 2)),

          # Means & SEs (Likert/Accuracy only)
          !is.na(M_treatment) & !is.na(M_control) &
            !is.na(SE_treatment) & !is.na(SE_control) &
            Measure_type %in% c("Likert", "Accuracy") ~
              (M_treatment - M_control) /
                sqrt(((n_treatment - 1)*(SE_treatment*sqrt(n_treatment))^2 +
                      (n_control  - 1)*(SE_control*sqrt(n_control))^2) /
                     (n_treatment + n_control - 2)),

          # Conversions not tied to scale type
          !is.na(F_value) ~ sqrt(F_value) * sqrt(1/n_treatment + 1/n_control),

          !is.na(unstandardised_beta) & !is.na(SE_beta) ~
            (unstandardised_beta / SE_beta) * sqrt(1/n_control + 1/n_treatment),

          !is.na(t_value) ~ t_value * sqrt(1/n_control + 1/n_treatment),

          !is.na(Odds_Ratio) & Odds_Ratio > 0 ~ log(Odds_Ratio) / 1.81,

          # Difference in performance with CIs
          !is.na(Performance_Difference) & 
            !is.na(`95CI_lower_performance_diff`) & 
            !is.na(`95CI_upper_performance_diff`) ~ {
              se_diff <- (`95CI_upper_performance_diff` - `95CI_lower_performance_diff`) / (2 * 1.96)
              n_pooled <- n_treatment + n_control
              ifelse(n_pooled > 0,
                     Performance_Difference / (se_diff * sqrt(2 * (n_treatment * n_control) / n_pooled)),
                     NA_real_)
            },

          # Means + CIs per group (Likert/Accuracy only)
          !is.na(M_treatment) & !is.na(M_control) &
            !is.na(`95CI_lower_treatment`) & !is.na(`95CI_upper_treatment`) &
            !is.na(`95CI_lower_control`) & !is.na(`95CI_upper_control`) &
            Measure_type %in% c("Likert", "Accuracy") ~ {
              SE_treat <- (`95CI_upper_treatment` - `95CI_lower_treatment`) / (2 * 1.96)
              SE_ctrl  <- (`95CI_upper_control` - `95CI_lower_control`) / (2 * 1.96)
              SD_treat <- SE_treat * sqrt(n_treatment)
              SD_ctrl  <- SE_ctrl  * sqrt(n_control)
              SD_pooled <- sqrt(((n_treatment - 1) * SD_treat^2 + (n_control - 1) * SD_ctrl^2) /
                                  (n_treatment + n_control - 2))
              (M_treatment - M_control) / SD_pooled
            },

          TRUE ~ NA_real_
        )
      )
    )
}

# Sampling variance of d 
compute_vi <- function(cohens_d, n_treatment, n_control) {
  ((n_treatment + n_control) / (n_treatment * n_control)) +
    (cohens_d^2 / (2 * (n_treatment + n_control)))
}



# Wrap-up: compute Hedges' g and its variance 
safe_calc_effects <- function(df) {
  if (!"cohens_d" %in% names(df) || any(is.na(df$cohens_d))) df <- compute_cohens_d(df)
  if (!"vi" %in% names(df) || any(is.na(df$vi))) df$vi <- compute_vi(df$cohens_d, df$n_treatment, df$n_control)

  df %>%
    mutate(
      df_total  = n_treatment + n_control - 2,
      J         = 1 - (3 / (4 * df_total - 1)),
      hedges_g  = cohens_d * J,
      vi_g      = vi * J^2
    )
}

# Apply effect-size pipeline 
df_performance <- safe_calc_effects(df_performance)

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
if (!requireNamespace("clubSandwich", quietly = TRUE)) install.packages("clubSandwich")
library(metafor)

required_cols <- c("hedges_g","vi_g","Experiment_ID","Effect_ID")
stopifnot(all(required_cols %in% names(df_performance)))

df_es <- df_performance %>%
  dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(Experiment_ID), !is.na(Effect_ID)) %>%
  dplyr::mutate(ES_ID = if ("ES_ID" %in% names(.)) ES_ID else paste0(Experiment_ID, "_", Effect_ID))

# Collapse rare levels per moderator (for stable inference) 
collapse_by_min_exp <- function(df, var, id_cluster = "Experiment_ID",
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

# =====================================================================
# Base meta-analytic model and Forest plot
# =====================================================================

main_model <- metafor::rma.mv(
  yi     = hedges_g,
  V      = vi_g,
  data   = df_es,
  random = ~ 1 | Experiment_ID/Effect_ID,
  method = "REML",
  tdist  = TRUE,
  level  = 95,
  digits = 7,
  slab   = ES_ID
)

# CR2-robustification (clustered at experiment) 
main_model_robust <- robust(main_model, cluster = df_es$Experiment_ID,
                            clubSandwich = TRUE, adjust = TRUE)

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

# Forest plot data 
author_col <- "Author"

model_est <- as.numeric(coef(main_model_robust))
model_lo  <- as.numeric(main_model_robust$ci.lb)
model_hi  <- as.numeric(main_model_robust$ci.ub)

df_es <- df_es %>%
  dplyr::mutate(w_raw = 1 / vi_g, w_pct = 100 * w_raw / sum(w_raw, na.rm = TRUE))

plot_data <- df_es %>%
  dplyr::mutate(
    effect = hedges_g,
    se     = sqrt(vi_g),
    ci.lb  = effect - 1.96 * se,
    ci.ub  = effect + 1.96 * se,
    weight = w_pct,
    study  = paste0(.data[[author_col]], " \u2014 ", ES_ID)
  ) %>%
  dplyr::arrange(dplyr::desc(effect)) %>%
  dplyr::group_by(Experiment_ID) %>%
  dplyr::mutate(label_grouped = paste0("  ", study)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    order = dplyr::row_number(),
    sign  = dplyr::case_when(ci.lb > 0 ~ "positive", ci.ub < 0 ~ "negative", TRUE ~ "ns"),
    col   = dplyr::case_when(sign == "positive" ~ "#006400", sign == "negative" ~ "#8B0000", TRUE ~ "grey40")
  )

x_rng          <- range(c(plot_data$ci.lb, plot_data$ci.ub, model_lo, model_hi), na.rm = TRUE)
x_pad_val      <- 0.15 * diff(x_rng)
x_pad          <- x_rng[2] + x_pad_val
plot_height_in <- nrow(plot_data) * 0.25 + 3.5

label_text <- sprintf("Overall Hedges' g = %.3f  [%0.3f, %0.3f]", model_est, model_lo, model_hi)

# Forest plot – using Author and numbers as labels 
plot_data_simple <- plot_data %>%
  arrange(desc(effect)) %>%
  dplyr::group_by(!!rlang::sym(author_col)) %>%
  dplyr::mutate(es_num = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(study_simple = paste0(.data[[author_col]], " — ", es_num),
                order_simple = dplyr::row_number())

readr::write_csv(
  plot_data_simple %>% dplyr::select(ES_ID, Author = !!rlang::sym(author_col), es_num),
  file.path(output_dir, "forest_labels_map_author_num.csv")
)

label2          <- "Overall Multilevel Model (Hedges' g)"
plot_filename2  <- "forest_multilevel_overall_gg_simple"
plot_height_in2 <- nrow(plot_data_simple) * 0.15 + 3.5

p2 <- ggplot(plot_data_simple, aes(x = effect, y = order_simple)) +
  annotate("rect", xmin = model_lo, xmax = model_hi, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.25) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6, alpha = 0.6) +
  geom_segment(aes(x = ci.lb, xend = ci.ub, y = order_simple, yend = order_simple, color = col), linewidth = 3, alpha = 0.35) +
  geom_point(aes(color = col), size = 1.8) +
  scale_color_identity(guide = "none") +
  geom_vline(xintercept = model_est, linetype = "dashed", color = "orange", linewidth = 0.7) +
  annotate("label", x = model_est, y = max(plot_data_simple$order_simple) + 1.2, label = label_text, fill = "white", colour = "orange", label.size = 0.2, size = 3.8) +
  geom_text(aes(x = x_pad, y = order_simple, label = sprintf("%.1f%%", weight)), hjust = 0, size = 3.5) +
  scale_y_reverse(breaks = plot_data_simple$order_simple, labels = plot_data_simple$study_simple, expand = expansion(add = c(0.3, 0.3))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y.left = element_text(size = 11, hjust = 0), axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), plot.margin = margin(5, 14, 5, 5, "mm")) +
  coord_cartesian(clip = "off") +
  labs(title = NULL, x = expression("Effect size (Hedges'" ~ italic(g) ~ ")"), y = NULL)

ggsave(file.path(output_dir, paste0(plot_filename2, ".pdf")), p2, width = 6.5, height = plot_height_in2, units = "in", dpi = 600)

ggsave(file.path(output_dir, paste0(plot_filename2, ".png")), p2, width = 6.5, height = plot_height_in2, units = "in", dpi = 300)


# Heterogeneity (tau^2, I^2 total and per level) ------------------------------
# Based on Viechtbauer's multilevel I^2

tau2_vec <- main_model$sigma2
if (length(tau2_vec) == 2) {
  names(tau2_vec) <- c("between_experiments", "within_experiments")
} else {
  names(tau2_vec) <- paste0("component_", seq_along(tau2_vec))
}

tau2_total <- sum(tau2_vec)

W    <- diag(1 / df_es$vi_g)
X    <- model.matrix(main_model)
XtWX <- t(X) %*% W %*% X
P    <- W - W %*% X %*% solve(XtWX) %*% t(X) %*% W

denom    <- tau2_total + (main_model$k - main_model$p) / sum(diag(P))
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


# =========================
# Multivariable meta-regression with effect coding (global intercept + all subgroups)
# =========================

# Moderator analysis: set up factors (desc vs infer) 
if (!requireNamespace("clubSandwich", quietly = TRUE)) install.packages("clubSandwich")
library(clubSandwich)

mods <- c("Task_Format","Explanation_Format_of_LLM","LLM_Model",
          "Input_Form_Patient_Information_in_LLM","Career_stage",
          "Type_of_Results","Medical_Field",
          "Follow_Up_Prompts_Allowed","Case_setting")
mods <- intersect(mods, names(df_es))

for (m in mods) df_es <- collapse_by_min_exp(df_es, m, min_exp = 2, other_label = "Other")


# 1) Prepare moderator list and modeling data
mmods <- paste0(mods, "_infer")
mmods <- mmods[mmods %in% names(df_es)]

df_mm <- df_es %>%
  dplyr::filter(is.finite(hedges_g), is.finite(vi_g),
                !is.na(Experiment_ID), !is.na(Effect_ID)) %>%
  dplyr::select(hedges_g, vi_g, Experiment_ID, Effect_ID, dplyr::all_of(mmods)) %>%
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
multi_fit_eff <- metafor::rma.mv(
  yi     = hedges_g,
  V      = vi_g,
  data   = df_mm,
  random = ~ 1 | Experiment_ID/Effect_ID,
  mods   = form_eff,
  method = "REML",
  tdist  = TRUE,
  level  = 95
)

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
  multi_fit_eff <- metafor::rma.mv(
    yi = hedges_g, V = vi_g, data = df_mm,
    random = ~ 1 | Experiment_ID/Effect_ID,
    mods   = form_eff, method = "REML", tdist = TRUE, level = 95
  )
  beta_names <- names(coef(multi_fit_eff))
  Xnames     <- colnames(multi_fit_eff$X)
}

# CR2 vcov
Vcr2_eff <- clubSandwich::vcovCR(multi_fit_eff, type = "CR2", cluster = df_mm$Experiment_ID)

# Intercept (grand mean)
ct <- clubSandwich::coef_test(multi_fit_eff, vcov = Vcr2_eff, test = "Satterthwaite")
df_col <- intersect(c("df_Satt","d.f.","df"), colnames(ct))[1]
row_int <- rownames(ct) %in% c("intrcpt","(Intercept)")
grand_mean_tbl <- tibble::tibble(
  Moderator="(Global)", Level="(Intercept / Grand mean)",
  estimate = ct$beta[row_int],
  SE       = ct$SE[row_int],
  df       = ct[[df_col]][row_int],
  t        = ct$tstat[row_int],
  p        = ct[[intersect(c("p_Satt","p_val","p_value","p"), names(ct))[1]]][row_int],
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

  # Try Satterthwaite dfs from CR2 for each single-coefficient contrast; fallback to normal if needed
  wt  <- try(clubSandwich::Wald_test(multi_fit_eff, constraints = L, vcov = Vcr2_eff, test = "EDT"), silent = TRUE)
  if (!inherits(wt, "try-error") && !anyNA(wt$df_denom)) {
    dfS <- as.numeric(wt$df_denom)
  } else {
    # conservative fallback: NA dfs => use normal-based CI/p (you’ll see df = NA)
    dfS <- rep(NA_real_, k)
  }

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
                          k=multi_fit_eff$k, clusters=dplyr::n_distinct(df_mm$Experiment_ID)))
  }
  R <- matrix(0, nrow = length(cols), ncol = length(b), dimnames = list(cols, names(b)))
  for (i in seq_along(cols)) R[i, cols[i]] <- 1
  wt <- try(clubSandwich::Wald_test(multi_fit_eff, constraints = R, vcov = Vcr2_eff, test = "HTZ"), silent = TRUE)
  if (inherits(wt, "try-error")) {
    tibble::tibble(Moderator=v, df1=NA_real_, df2=NA_real_, F=NA_real_, p=NA_real_,
                   k=multi_fit_eff$k, clusters=dplyr::n_distinct(df_mm$Experiment_ID))
  } else {
    tibble::tibble(
      Moderator=v, df1=unname(wt$df_num), df2=unname(wt$df_denom),
      F=unname(wt$Fstat), p=unname(wt$p_val),
      k=multi_fit_eff$k, clusters=dplyr::n_distinct(df_mm$Experiment_ID)
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


# =====================================================================
# Moderator analyses
# =====================================================================


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

    fit1 <- metafor::rma.mv(yi = hedges_g, V = vi_g,
                             random = ~ 1 | Experiment_ID/Effect_ID,
                             data = df_tmp, method = "REML")

    per_level <- try({
      Vcr2 <- clubSandwich::vcovCR(fit1, type = "CR2", cluster = df_tmp$Experiment_ID)
      out  <- clubSandwich::coef_test(fit1, vcov = Vcr2, test = "Satterthwaite") %>% as.data.frame()
      df_col <- intersect(c("df_Satt","d.f.","df"), names(out))[1]
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
      df          = per_level[[intersect(c("df_Satt","d.f.","df"), names(per_level))[1]]],
      t.value     = per_level$tstat,
      p.value     = per_level[[intersect(c("p_Satt","p_val","p_value","p"), names(per_level))[1]]],
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
  fit <- metafor::rma.mv(
    yi = hedges_g, V = vi_g,
    mods   = ~ 0 + .mod_fac,
    random = ~ 1 | Experiment_ID/Effect_ID,
    data   = df_tmp, method = "REML"
  )

  per_level <- try({
    Vcr2 <- clubSandwich::vcovCR(fit, type = "CR2", cluster = df_tmp$Experiment_ID)
    out  <- clubSandwich::coef_test(fit, vcov = Vcr2, test = "Satterthwaite") %>% as.data.frame()
    out$term <- rownames(out)
    df_col <- intersect(c("df_Satt","d.f.","df"), names(out))[1]
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
    df         = per_level[[intersect(c("df_Satt","d.f.","df"), names(per_level))[1]]],
    t.value    = per_level[[intersect(c("tstat","t","Tstat"), names(per_level))[1]]],
    p.value    = per_level[[intersect(c("p_Satt","p_val","p_value","p"), names(per_level))[1]]],
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

# NO-INTERCEPT omnibus: equality of means across levels (collapsed), here we test if all sub-group means are equal -----------
omni_noint <- dplyr::tibble()
min_clusters_per_level_infer <- 2

for (mod in mods) {
  message("\n[NO-INT/INFER OMNIBUS] Moderator: ", mod)
  fac_infer <- paste0(mod, "_infer")

  df_tmp <- df_es %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(.data[[fac_infer]])) %>%
    dplyr::mutate(.mod_fac = droplevels(factor(.data[[fac_infer]])))

  lev_exp <- df_tmp %>%
    dplyr::group_by(.mod_fac) %>%
    dplyr::summarise(n_exp = dplyr::n_distinct(Experiment_ID), .groups = "drop")

  keep_levels <- lev_exp %>% dplyr::filter(n_exp >= min_clusters_per_level_infer) %>%
    dplyr::pull(.mod_fac) %>% as.character()

  df_tmp <- df_tmp %>%
    dplyr::filter(as.character(.mod_fac) %in% keep_levels) %>%
    dplyr::mutate(.mod_fac = droplevels(.mod_fac))

  if (nlevels(df_tmp$.mod_fac) < 2) { message("  -> Skip omnibus (no-intercept): too few analyzable levels after sparsity filter."); next }

  df_tmp$.mod_fac <- factor(df_tmp$.mod_fac, levels = sort(levels(df_tmp$.mod_fac)))
  levs <- levels(df_tmp$.mod_fac)

  fit_noint <- metafor::rma.mv(
    yi = hedges_g, V = vi_g,
    mods   = ~ 0 + .mod_fac,
    random = ~ 1 | Experiment_ID/Effect_ID,
    data   = df_tmp, method = "REML"
  )

  Vcr2_noint <- clubSandwich::vcovCR(fit_noint, type = "CR2", cluster = df_tmp$Experiment_ID)

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

  omni_ni <- clubSandwich::Wald_test(fit_noint, constraints = R, vcov = Vcr2_noint, test = "HTZ")

  omni_noint <- dplyr::bind_rows(
    omni_noint,
    dplyr::tibble(
      Moderator = mod,
      test      = "CR2 HTZ (equality of level means) on COLLAPSED levels (no-intercept)",
      df1       = unname(omni_ni$df_num),
      df2       = unname(omni_ni$df_denom),
      F         = unname(omni_ni$Fstat),
      p         = unname(omni_ni$p_val),
      levels    = paste(levs, collapse = "|")
    )
  )
}

readr::write_csv(omni_noint, file.path(output_dir, "moderator_NOINTERCEPT_OMNIBUS_CR2_EQUAL_MEANS.csv"))

# NO-INTERCEPT omnibus: to test whether all subgroup means are equal to zero ("is there any effect at all?")  ----------------------
omni_noint_allzero <- dplyr::tibble()
min_clusters_per_level_infer <- 2

for (mod in mods) {
  message("\n[NO-INT/INFER ALL-ZERO] Moderator: ", mod)
  fac_infer <- paste0(mod, "_infer")

  df_tmp <- df_es %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(.data[[fac_infer]])) %>%
    dplyr::mutate(.mod_fac = droplevels(factor(.data[[fac_infer]])))

  lev_exp <- df_tmp %>% dplyr::group_by(.mod_fac) %>% dplyr::summarise(n_exp = dplyr::n_distinct(Experiment_ID), .groups = "drop")
  keep_levels <- lev_exp %>% dplyr::filter(n_exp >= min_clusters_per_level_infer) %>% dplyr::pull(.mod_fac) %>% as.character()

  df_tmp <- df_tmp %>% dplyr::filter(as.character(.mod_fac) %in% keep_levels) %>% dplyr::mutate(.mod_fac = droplevels(.mod_fac))

  if (nlevels(df_tmp$.mod_fac) < 1) { message("  -> Skip all-zero omnibus: too few analyzable levels after sparsity filter."); next }

  df_tmp$.mod_fac <- factor(df_tmp$.mod_fac, levels = sort(levels(df_tmp$.mod_fac)))

  fit0 <- metafor::rma.mv(
    yi = hedges_g, V = vi_g,
    mods   = ~ 0 + .mod_fac,
    random = ~ 1 | Experiment_ID/Effect_ID,
    data   = df_tmp, method = "REML"
  )

  Vcr2  <- clubSandwich::vcovCR(fit0, type = "CR2", cluster = df_tmp$Experiment_ID)

  vc_names <- colnames(Vcr2); if (is.null(vc_names)) vc_names <- names(coef(fit0))
  R <- diag(length(vc_names)); rownames(R) <- colnames(R) <- vc_names

  wt <- clubSandwich::Wald_test(fit0, constraints = R, vcov = Vcr2, test = "HTZ")

  omni_noint_allzero <- dplyr::bind_rows(
    omni_noint_allzero,
    dplyr::tibble(
      Moderator = mod,
      test      = "CR2 HTZ (all subgroup means = 0) on COLLAPSED levels (no-intercept)",
      df1       = unname(wt$df_num),
      df2       = unname(wt$df_denom),
      F         = unname(wt$Fstat),
      p         = unname(wt$p_val),
      levels    = paste(levels(df_tmp$.mod_fac), collapse = "|")
    )
  )
}

readr::write_csv(omni_noint_allzero, file.path(output_dir, "moderator_NOINTERCEPT_OMNIBUS_CR2_ALLZERO.csv"))

# WITH INTERCEPT: tests whether subgroups means are different from one other in relation to one subgroup mean (which is the reference category) ---------------
omni_with_int <- dplyr::tibble()
min_clusters_per_level_infer <- 2

for (mod in mods) {
  message("\n[WITH INT/INFER] Moderator: ", mod)
  fac_infer <- paste0(mod, "_infer")

  df_tmp <- df_es %>%
    dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(.data[[fac_infer]])) %>%
    dplyr::mutate(.mod_fac = factor(.data[[fac_infer]]))

  lev_exp <- df_tmp %>% dplyr::group_by(.mod_fac) %>% dplyr::summarise(n_exp = dplyr::n_distinct(.data$Experiment_ID), .groups = "drop")

  if (nlevels(df_tmp$.mod_fac) < 2 || any(lev_exp$n_exp < min_clusters_per_level_infer)) {
    message("  -> Skip omnibus: collapsed factor still too sparse."); next
  }

  df_tmp$.mod_fac <- factor(df_tmp$.mod_fac, levels = sort(levels(df_tmp$.mod_fac)))
  ref_level <- levels(df_tmp$.mod_fac)[1]

  fit_wi <- metafor::rma.mv(
    yi = hedges_g, V = vi_g,
    mods   = ~ .mod_fac,
    random = ~ 1 | Experiment_ID/Effect_ID,
    data   = df_tmp, method = "REML"
  )

  Vcr2_wi   <- clubSandwich::vcovCR(fit_wi, type = "CR2", cluster = df_tmp$Experiment_ID) 
  beta_names <- names(coef(fit_wi))
  int_name   <- if ("intrcpt" %in% beta_names) "intrcpt" else "(Intercept)"
  non_int    <- setdiff(beta_names, int_name)
  R          <- diag(length(beta_names))[match(non_int, beta_names), , drop = FALSE]
  rownames(R) <- non_int; colnames(R) <- beta_names

  omni_wi <- clubSandwich::Wald_test(fit_wi, constraints = R, vcov = Vcr2_wi, test = "HTZ")

  omni_with_int <- dplyr::bind_rows(
    omni_with_int,
    dplyr::tibble(
      Moderator = mod,
      test      = "CR2 HTZ (contrasts vs intercept) on COLLAPSED levels",
      df1       = unname(omni_wi$df_num),
      df2       = unname(omni_wi$df_denom),
      F         = unname(omni_wi$Fstat),
      p         = unname(omni_wi$p_val),
      ref_level = ref_level,
      levels    = paste(levels(df_tmp$.mod_fac), collapse = "|")
    )
  )
}

readr::write_csv(omni_with_int, file.path(output_dir, "moderator_WITHINT_OMNIBUS_CR2_COLLAPSED.csv"))


# =====================================================================
#Robustness and sensitivity checks
# =====================================================================

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

# Egger's regression (multilevel, CR2-robust) 

eggers_model <- metafor::rma.mv(
  yi     = hedges_g,
  V      = vi_g,
  data   = df_es,
  mods   = ~ sqrt(vi_g),
  random = ~ 1 | Experiment_ID/Effect_ID,
  method = "REML",
  tdist  = TRUE,
  level  = 95,
  digits = 7,
  slab   = ES_ID
)

eggers_model_robust <- metafor::robust(
  eggers_model,
  cluster      = df_es$Experiment_ID,
  clubSandwich = TRUE,
  adjust       = TRUE
)

print(summary(eggers_model_robust))

eg_coefs <- coef(eggers_model_robust)
eg_vcov  <- vcov(eggers_model_robust)
eg_se    <- sqrt(diag(eg_vcov))
eg_stat  <- eg_coefs / eg_se
eg_df    <- eggers_model_robust$ddf
if (!is.null(eg_df) && all(is.finite(eg_df) & eg_df > 0)) {
  eg_p <- 2 * pt(abs(eg_stat), df = eg_df, lower.tail = FALSE)
} else {
  eg_p <- 2 * pnorm(abs(eg_stat), lower.tail = FALSE)
}
cat("\nEgger (robust) coefficients:\n"); print(data.frame(estimate = eg_coefs, se = eg_se, stat = eg_stat, p = eg_p))

# Begg–Mazumdar rank correlation test 
if (!exists("main_model")) {
  main_model <- metafor::rma.mv(
    yi     = hedges_g,
    V      = vi_g,
    data   = df_es,
    random = ~ 1 | Experiment_ID/Effect_ID,
    method = "REML",
    tdist  = TRUE,
    level  = 95,
    digits = 7,
    slab   = ES_ID
  )
}

ranktest_res <- tryCatch(
  metafor::ranktest(main_model),
  error = function(e) {
    message("ranktest(main_model) not available for this object; using yi/vi fallback.")
    metafor::ranktest(yi = df_es$hedges_g, vi = df_es$vi_g)
  }
)
print(ranktest_res)

# Sensitivity Analysis: Clustering at PAPER (used ID_2 as variable) 
stopifnot("ID_2" %in% names(df_es))
stopifnot(all(c("hedges_g","vi_g","Effect_ID") %in% names(df_es)))

df_paper <- df_es %>%
  dplyr::filter(is.finite(hedges_g), is.finite(vi_g), !is.na(ID_2), !is.na(Effect_ID))

if (nrow(df_paper) < 2) stop("Not enough rows after filtering.")
if (dplyr::n_distinct(df_paper$ID_2) < 2) stop("Need at least 2 papers (ID_2) to compute CR2-robust inference.")

main_model_paper <- metafor::rma.mv(
  yi     = hedges_g,
  V      = vi_g,
  data   = df_paper,
  random = ~ 1 | ID_2/Effect_ID,
  method = "REML",
  tdist  = TRUE,
  level  = 95,
  digits = 7,
  slab   = if ("ES_ID" %in% names(df_paper)) df_paper$ES_ID else NULL
)

main_model_paper_rb <- metafor::robust(
  main_model_paper,
  cluster      = df_paper$ID_2,
  clubSandwich = TRUE,
  adjust       = TRUE
)

print(summary(main_model_paper_rb))

beta <- as.numeric(coef(main_model_paper_rb))[1]
Vrb  <- as.matrix(vcov(main_model_paper_rb))
se   <- sqrt(pmax(diag(Vrb)[1], 0))

ddf  <- main_model_paper_rb$ddf
has_t <- !is.null(ddf) && is.finite(ddf[1]) && ddf[1] > 0

stat <- beta / se
if (has_t) {
  pval <- 2 * stats::pt(abs(stat), df = ddf[1], lower.tail = FALSE)
  crit <- stats::qt(0.975, df = ddf[1])
  df_out <- ddf[1]
} else {
  pval <- 2 * stats::pnorm(abs(stat), lower.tail = FALSE)
  crit <- stats::qnorm(0.975)
  df_out <- NA_real_
}
ci_lo <- beta - crit * se
ci_hi <- beta + crit * se

sens_summary <- tibble::tibble(
  clusters = dplyr::n_distinct(df_paper$ID_2),
  k        = nrow(df_paper),
  estimate = beta,
  se       = se,
  stat     = stat,
  df       = df_out,
  p        = pval,
  ci.lb    = ci_lo,
  ci.ub    = ci_hi
)

print(sens_summary)
readr::write_csv(sens_summary, file.path(output_dir, "sensitivity_paper_cluster_summary.csv"))

# =====================================================================
# Leave-one-out Analysis
# =====================================================================

# Effect Size level (now again clustering at experiment) 
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

  loo_model <- metafor::rma.mv(
    yi     = hedges_g,
    V      = vi_g,
    data   = dat,
    random = ~ 1 | Experiment_ID/Effect_ID,
    method = "REML",
    tdist  = TRUE,
    level  = 95,
    slab   = dat$ES_ID
  )

  loo_rb <- tryCatch(
    metafor::robust(loo_model, cluster = dat$Experiment_ID, clubSandwich = TRUE, adjust = TRUE),
    error = function(e) NULL
  )

  if (!is.null(loo_rb)) {
    beta <- as.numeric(coef(loo_rb))[1]
    Vrb  <- as.matrix(vcov(loo_rb))
    se   <- sqrt(pmax(diag(Vrb)[1], 0))
    ddf  <- loo_rb$ddf
    has_t <- !is.null(ddf) && is.finite(ddf[1]) && ddf[1] > 0

    stat <- beta / se
    if (has_t) {
      pval <- 2 * stats::pt(abs(stat), df = ddf[1], lower.tail = FALSE)
      crit <- stats::qt(0.975, df = ddf[1])
      df_out <- ddf[1]
    } else {
      pval <- 2 * stats::pnorm(abs(stat), lower.tail = FALSE)
      crit <- stats::qnorm(0.975)
      df_out <- NA_real_
    }
    ci_lo <- beta - crit * se
    ci_hi <- beta + crit * se
  } else {
    beta <- se <- stat <- pval <- ci_lo <- ci_hi <- df_out <- NA_real_
  }

  I2_total <- NA_real_
  tau2_total <- sum(loo_model$sigma2)
  I2_try <- try({
    W   <- diag(1 / dat$vi_g)
    X   <- model.matrix(loo_model)
    XtWX <- t(X) %*% W %*% X
    P   <- W - W %*% X %*% solve(XtWX) %*% t(X) %*% W
    denom <- tau2_total + (loo_model$k - loo_model$p) / sum(diag(P))
    I2_total <- 100 * tau2_total / denom
  }, silent = TRUE)
  if (inherits(I2_try, "try-error")) I2_total <- NA_real_

  df_loo_es$k[i]        <- loo_model$k
  df_loo_es$clusters[i] <- dplyr::n_distinct(dat$Experiment_ID)
  df_loo_es$estimate[i] <- beta
  df_loo_es$se[i]       <- se
  df_loo_es$stat[i]     <- stat
  df_loo_es$df[i]       <- df_out
  df_loo_es$p_val[i]    <- pval
  df_loo_es$ci.lb[i]    <- ci_lo
  df_loo_es$ci.ub[i]    <- ci_hi
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

.i2_total <- function(fit, vi_vec) {
  tau2_total <- sum(fit$sigma2)
  out <- try({
    W    <- diag(1 / vi_vec)
    X    <- model.matrix(fit)
    XtWX <- t(X) %*% W %*% X
    P    <- W - W %*% X %*% solve(XtWX) %*% t(X) %*% W
    denom <- tau2_total + (fit$k - fit$p) / sum(diag(P))
    100 * tau2_total / denom
  }, silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_) else return(as.numeric(out))
}

# LOO at Experiment level 
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

  fit <- metafor::rma.mv(
    yi     = hedges_g, V = vi_g,
    data   = dat,
    random = ~ 1 | Experiment_ID/Effect_ID,
    method = "REML",
    tdist  = TRUE,
    level  = 95,
    slab   = dat$ES_ID
  )

  rb  <- tryCatch(metafor::robust(fit, cluster = dat$Experiment_ID, clubSandwich = TRUE, adjust = TRUE), error = function(e) NULL)

  inf <- .extract_rb(rb)
  I2t <- .i2_total(fit, dat$vi_g)

  df_loo_exp$k[i]        <- fit$k
  df_loo_exp$clusters[i] <- dplyr::n_distinct(dat$Experiment_ID)
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

# LOO at Paper level 
paper_ids <- sort(unique(df_use$ID_2))

df_loo_paper <- tibble::tibble(
  left_out_ID_2 = paper_ids,
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
  dat  <- dplyr::filter(df_use, ID_2 != pid)
  if (nrow(dat) < 2 || dplyr::n_distinct(dat$ID_2) < 2) next

  fit <- metafor::rma.mv(
    yi     = hedges_g, V = vi_g,
    data   = dat,
    random = ~ 1 | Experiment_ID/Effect_ID,
    method = "REML",
    tdist  = TRUE,
    level  = 95,
    slab   = dat$ES_ID
  )

  rb  <- tryCatch(metafor::robust(fit, cluster = dat$ID_2, clubSandwich = TRUE, adjust = TRUE), error = function(e) NULL)

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

readr::write_csv(df_loo_paper, file.path(output_dir, "loo_paper_level_cluster_ID_2.csv"))

loo_dir <- output_dir

# Formatting helpers for LaTeX tables 
fmt_p <- function(p) {
  p <- as.numeric(p)
  ifelse(is.na(p), "", ifelse(p < 0.001, "$<0.001$", sprintf("$= %.3f$", p)))
}

fmt_num <- function(x, digits = 2) {
  x <- as.numeric(x)
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_g_ci <- function(g, lo, hi, digits = 2) {
  g  <- as.numeric(g); lo <- as.numeric(lo); hi <- as.numeric(hi)
  has_all <- !(is.na(g) | is.na(lo) | is.na(hi))
  out <- rep("", length(g))
  out[has_all] <- paste0(
    "$", fmt_num(g[has_all],  digits), "$ (",
    "$", fmt_num(lo[has_all], digits), "$–$",
    "$", fmt_num(hi[has_all], digits), "$)"
  )
  out
}

# Read LOO outputs 
loo_es  <- read_csv(file.path(loo_dir, "loo_effect_level_experiment_cluster.csv"), show_col_types = FALSE)
loo_exp <- read_csv(file.path(loo_dir, "loo_experiment_level_cluster_Experiment_ID.csv"), show_col_types = FALSE)
loo_pap <- read_csv(file.path(loo_dir, "loo_paper_level_cluster_ID_2.csv"), show_col_types = FALSE)


# =====================================================================
# === Format meta-regression results (ALL_LEVELS_CR2_TABLE) for LaTeX ===
# =====================================================================

library(kableExtra)

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

message("✅ Clean LaTeX table saved to: ", latex_path)