# Can Targeting Nudges Improve Enrollment in California's ACA Health Insurance Marketplace?

Replication code for the heterogeneity and policy-learning analysis: multi-arm
causal forests, CATE (conditional average treatment effect) estimation, and
policy-tree targeting rules, applied to a 5-arm RCT (n = 87,394) testing
outreach letters that encourage ACA marketplace enrollment.

## Project overview

The pipeline takes the trial's public-use analysis file and:

1. Confirms randomization balance across the five arms (Table 1).
2. Fits a multi-arm causal forest to estimate household-level treatment
   effects (CATEs) for each letter vs. Control, and reports these overall,
   by subgroup, and by their full distribution (Figures 1–2, Appendix A2).
3. Estimates arm-specific variable importance from one binary causal forest
   per arm (Appendix A1).
4. Fits doubly-robust (AIPW) policy-learning inputs and cross-fitted policy
   trees (depth 2–4) plus a plug-in targeting rule, and compares their
   out-of-sample value against "treat everyone" / "treat no one" baselines
   (Table 2).
5. Translates the fitted policy trees' binned splits back into real-world
   units (age in years, income in % FPL, subsidy in $/month) and renders
   them as diagrams (Figure 3, + depth-3/4 trees).

All figures are exported as **PDF** (vector, publication-quality); tables
are exported as **CSV** and **LaTeX/HTML**.

## Software and dependencies

- **R** (developed/tested under a recent R 4.x release)
- CRAN packages:
  `policytree`, `grf`, `dplyr`, `tidyr`, `haven`, `sandwich`, `lmtest`,
  `kableExtra`, `ggplot2`, `scales`, `cowplot`, `patchwork`, `forcats`,
  `future.apply`, `progressr`, `DiagrammeR`, `DiagrammeRsvg`, `rsvg`,
  `digest`
- Install everything with:
  ```r
  install.packages(c(
    "policytree", "grf", "dplyr", "tidyr", "haven", "sandwich", "lmtest",
    "kableExtra", "ggplot2", "scales", "cowplot", "patchwork", "forcats",
    "future.apply", "progressr", "DiagrammeR", "DiagrammeRsvg", "rsvg",
    "digest"
  ))
  ```
- **Data:** `AER_2019_0823_analysisfile_public.dta` (Stata `.dta` format,
  read via `haven`) — not included in this repo; you must supply your own
  copy (see Replication Instructions).

## Replication instructions

1. **Get the folder structure right.** The pipeline expects:
   ```
   PROJECT_DIR/
     master.R
     code/                 <- the 11 numbered .R scripts
     data/
       AER_2019_0823_analysisfile_public.dta
     output/                <- created automatically; figures/tables land here
   ```
2. **Edit one line.** Open `master.R` and set `PROJECT_DIR` to the folder
   above (it must contain `code/` and `data/`).
3. **Place the data file** inside `PROJECT_DIR/data/`.
4. **Run it.** Either source `master.R` in RStudio (Session > Source), or
   from a terminal:
   ```bash
   Rscript master.R
   ```
   This runs all 11 scripts in order and writes every figure/table to
   `PROJECT_DIR/output/`.
5. **Runtime note:** step `09_policy_evaluation.R` fits 15 cross-fitted
   policy/hybrid-policy trees (depths 2–4 × 5 folds) plus 3 final full-sample
   trees. This is the slow step; it uses `future.apply` for parallelism
   (defaults to `min(8, cores - 1)` workers) and caches its results to
   `output/tree_fits_cache_*.rds` and `output/final_trees_cache_*.rds`, so
   re-running `master.R` after the first successful run will load from
   cache rather than refitting, unless `policy_covs`, `SPLIT`,
   `SEARCH_DEPTH`, or the fold assignment change.
6. **Reproducibility:** the seed is fixed once, in `01_setup.R`
   (`set.seed(21)`), and every forest/tree call also passes `seed = 21`
   directly, so results should replicate exactly.

## Output map: figures and tables by script

| Output | File(s) | Produced by | Description |
|---|---|---|---|
| Table 1 | `balance_table.csv`, `balance_table.tex` | `03_balance_table.R` | By-arm covariate means, overall mean, and a robust (HC1) joint F-test that arm means are equal — confirms randomization balance. |
| Figure 1 | `fig2_cate_plot_only.pdf`, `fig2_cate_legend_only.pdf` | `06_cate_distribution.R` | Kernel-density plot of household-level predicted treatment effects for each letter, with each arm's ATE marked. Plot and legend saved as separate PDFs for layout flexibility. |
| Figure 2 | `cate_fpl.pdf`, `cate_aptc.pdf`, `cate_age.pdf`, `cate_funnel_sws.pdf` | `05_cate_subgroup_plots.R` | CATE-by-subgroup point-range panels: by income (% FPL), by subsidy amount, by household age, and by enrollment funnel (open enrollment vs. county referral). |
| Appendix Figure A1 | `variable_importance_pooled_top10.pdf`, `variable_importance_by_arm_labeled.pdf`, `variable_importance_by_arm.csv` | `07_variable_importance_by_arm.R` | Panel A: pooled top-10 covariate importance from the multi-arm forest. Panel B: top-10 importance per arm, from separate binary causal forests (arm vs. Control). |
| Appendix Figure A2 (+ related panels) | `cate_race.pdf`, `cate_language.pdf`, `cate_hh_members.pdf`, `cate_email.pdf`, `cate_regions_1.pdf`, `cate_regions_2.pdf`, `cate_marriage.pdf`, `cate_gender.pdf`, `cate_subgroup_legend_horizontal.pdf` | `05_cate_subgroup_plots.R` | Additional CATE-by-subgroup panels (race, language, household size, e-mail availability, California rating region, marital status, gender) and one shared horizontal legend used across all subgroup figures. |
| Table 2 | `policy_comparison_table.html`, `policy_comparison_table.tex`, `policy_comparison_long_holm.csv` | `09_policy_evaluation.R` | Pairwise out-of-sample (AIPW) value comparison across policies: doing nothing, always-treat-with-arm-X, cross-fitted policy trees (depth 2–4), and the fairness-restricted plug-in rule. HTML/TeX show unadjusted stars; the CSV also carries Holm-adjusted p-values for every pair. |
| Figure 3 (+ depth-3/4 trees) | `policy_tree_depth2_5arms.pdf`, `policy_tree_depth3_5arms.pdf`, `policy_tree_depth4_5arms.pdf` | `11_plot_policy_trees.R` | Graphviz diagrams of the final, full-sample-fitted policy trees (depths 2–4), with arm-matched leaf colors and splits labeled in real-world units (via `10_translate_tree_splits.R`'s bin lookup). |
| — (console only) | *(no file — printed output)* | `10_translate_tree_splits.R` | Prints each final policy tree with binned splits translated back into real units (e.g. "Subsidy ≤ $129/mth"), used to write up the targeting-rule narrative. |

### Scripts with no direct output (intermediate/setup steps)

| Script | Role |
|---|---|
| `01_setup.R` | Loads packages, sets the seed (21), and defines shared tuning constants (`N_TREES`, `SPLIT`, `SEARCH_DEPTH`, `PP_SCALE`). |
| `02_prepare_data.R` | Loads the `.dta` analysis file, recodes outcome/treatment, builds the covariate matrix `X`, and defines display names/colors reused by every later figure. |
| `04_heterogeneity_forest.R` | Fits the multi-arm causal forest (`forest_est`) used for CATE plots and subgroup tests; feeds Figures 1–2, Appendix A2, and Appendix A1 (pooled importance). |
| `08_policy_scoring_targeting_forests.R` | Fits the cross-fitted scoring forest (produces the AIPW reward matrix `Gamma`) and the coarse, fairness-restricted targeting forest used by the plug-in rule; both feed `09_policy_evaluation.R`. |

## Notes on data availability

The underlying trial microdata (`AER_2019_0823_analysisfile_public.dta`) is sourced from:
Domurat, Richard, Menashe, Isaac , and Yin, Wesley. Data and code for “The Role of Behavioral Frictions in Health Insurance Marketplace Enrollment and Risk: Evidence from a Field Experiment.” Nashville, TN: American Economic Association [publisher], 2021. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2021-04-09. https://doi.org/10.3886/E125801V1
