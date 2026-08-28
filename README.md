# The Effect of LLM Assistance on Diagnostic Accuracy: A Systematic Review and Meta-Analysis


## Project Description

This systematic review and meta-analysis examines whether the assistance of large language models (LLMs) affects the diagnostic accuracy of physicians.

The review includes comparative studies in which physicians performed diagnostic tasks with and without LLM assistance. The literature search covered MEDLINE, Web of Science, EMBASE, Scopus, and medRxiv and included studies published between January 2020 and June 2026.

The repository provides the coded data and R code used for the statistical analyses reported in the manuscript.

## Repository Contents

- `Data/`
  - `Coding_Scheme_Final.csv`: The data extracted from the studies included in the meta-analysis.
  - `20260828_PROBAST_vFinal`: Results of the risk-of-bias assessment.

- `Code/`
  - `Meta_Analysis_Code_vFINAL.r`: Main R script used for data processing, effect-size calculation, meta-analytic models, moderator analyses, sensitivity analyses, and generation of figures and tables.

## Statistical Analysis

Effect sizes are expressed as Hedges' *g*, with positive values indicating higher diagnostic accuracy with LLM assistance.

The main analysis uses a three-level random-effects meta-analytic model estimated using restricted maximum likelihood (REML), accounting for dependencies among multiple effect sizes originating from the same experiment.

Statistical inference is based on cluster-robust variance estimation with small-sample corrections (CR2), clustered at the paper/shared-case-corpus level.

The analysis script additionally includes:

- moderator and subgroup analyses;
- heterogeneity analyses;
- sensitivity analyses for assumptions regarding within-subject correlations and physician-level clustering;
- robustness analyses for study design, publication status, and risk of bias;
- small-study-effect analyses;
- leave-one-out analyses; and
- generation of the figures and supplementary analysis outputs reported in the manuscript.

## Reproducing the Analysis

The analyses were conducted in R.

To reproduce the analyses:

1. Clone or download this repository.
2. Keep the `Code/` and `Data/` directories in their current structure.
3. Open `Meta_Analysis_Code_vFINAL.r` from the `Code/` directory.
4. Run the script from top to bottom.

The script imports the coded data from the `Data/` directory and creates the required output directories automatically.

The primary analyses were conducted using R 4.5.0 and the `metafor` and `clubSandwich` packages. Additional R packages used for data processing, visualization, and table generation are loaded within the analysis script.


## Citation

If you use the data or code from this repository, please cite the associated manuscript.

Citation details will be added upon publication.
