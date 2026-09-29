# =====================================================================
# master.R
# Can Targeting Nudges Improve Enrollment in California's ACA Health
# Insurance Marketplace? -- full replication pipeline
# =====================================================================
#
# HOW TO RUN
# ----------
# 1. Edit ONE line below: PROJECT_DIR.
# 2. Put the data file, AER_2019_0823_analysisfile_public.dta, inside
#    PROJECT_DIR/data/  (create the "data" folder if it doesn't exist).
# 3. Source this file (in RStudio: Session > Source, or from a terminal:
#    Rscript master.R). Everything else runs automatically -- figures and
#    tables are written to PROJECT_DIR/output/, which is created for you.
#
# Expected folder layout:
#   PROJECT_DIR/
#     master.R
#     code/                 <- the 11 scripts
#     data/
#       AER_2019_0823_analysisfile_public.dta   <-add this file
#     output/            <- created automatically; figures/tables land here


# OUTPUT
# ---------------------------------------------------------------------
#   01_setup.R                       no output (packages, seed, constants)
#   02_prepare_data.R                no output (builds analysis sample)
#   03_balance_table.R               Table 1
#   04_heterogeneity_forest.R        no standalone figure (feeds 05-07)
#   05_cate_subgroup_plots.R         Figure 2 + Appendix Figure A2 panels
#   06_cate_distribution.R           Figure 1
#   07_variable_importance_by_arm.R  Appendix Figure A1
#   08_policy_scoring_targeting_forests.R   no output (feeds 09)
#   09_policy_evaluation.R           Table 2 (+ final_trees for 10-11)
#   10_translate_tree_splits.R       console output only (real-unit splits)
#   11_plot_policy_trees.R           Figure 3 (+ depth-3/4 trees)
# =====================================================================

# ---- 1. EDIT THIS LINE --------------------------
PROJECT_DIR <- "/PATH/policy_learning_code"   # edit with the path to the folder containing this file

# ---- 2. Nothing below this line needs to change ---------------------------
DATA_PATH <- file.path(PROJECT_DIR, "data", "AER_2019_0823_analysisfile_public.dta")
OUT_DIR   <- file.path(PROJECT_DIR, "output/")
CODE_DIR     <- file.path(PROJECT_DIR, "code")

stopifnot(
  "Can't find the R/ scripts folder inside PROJECT_DIR -- check the path." = dir.exists(CODE_DIR),
  "Can't find the data file. Put AER_2019_0823_analysisfile_public.dta inside PROJECT_DIR/data/." = file.exists(DATA_PATH)
)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

run <- function(script, produces = "no direct output") {
  cat("\n============================================================\n")
  cat("Running:", script, "\n")
  cat("Produces:", produces, "\n")
  cat("============================================================\n")
  source(file.path(CODE_DIR, script))
}

run("01_setup.R",
    "no output -- loads packages, sets seed 21, defines N_TREES/SPLIT/SEARCH_DEPTH/PP_SCALE")

run("02_prepare_data.R",
    "no output -- loads the .dta file (n = 87,394), builds df, X, Y, W, W.hat, arm colours/labels, covariate list")

run("03_balance_table.R",
    "Table 1 (Summary Statistics, by arm) -> balance_table.csv, balance_table.tex")

run("04_heterogeneity_forest.R",
    "no standalone figure -- fits forest_est/tau_hat (used by 05, 06) and prints pooled variable importance (used by 07)")

run("05_cate_subgroup_plots.R",
    paste("Figure 2 (CATE by income/subsidy/age/funnel round) + Appendix Figure A2 panels",
          "(race, language, region, household size, e-mail, marital status, gender)",
          "-> cate_*.pdf, cate_subgroup_legend_horizontal.pdf"))

run("06_cate_distribution.R",
    "Figure 1 (Distribution of household-level predicted treatment effects) -> fig2_cate_plot_only.pdf, fig2_cate_legend_only.pdf")

run("07_variable_importance_by_arm.R",
    "Appendix Figure A1 (pooled + by-arm variable importance) -> variable_importance_pooled_top10.pdf, variable_importance_by_arm_labeled.pdf, variable_importance_by_arm.csv")

run("08_policy_scoring_targeting_forests.R",
    "no output -- fits the cross-fitted scoring forest (Gamma) and the coarse fairness-restricted targeting forest used by the plug-in rule")

run("09_policy_evaluation.R",
    paste("Table 2 (Estimated reward of the policy rules vs. alternative approaches)",
          "-> policy_comparison_table.html, policy_comparison_table.tex, policy_comparison_long_holm.csv.",
          "Also fits final_trees (depth 2-4), used by scripts 10-11."))

run("10_translate_tree_splits.R",
    "no figures -- console output only: each final policy tree with bin indices translated to real-world units ($, % FPL, years)")

run("11_plot_policy_trees.R",
    "Figure 3 (Fitted depth-2 policy tree) + depth-3/depth-4 trees -> policy_tree_depth2_5arms.pdf, policy_tree_depth3_5arms.pdf, policy_tree_depth4_5arms.pdf")

cat("\nAll scripts finished. Outputs written to:", OUT_DIR, "\n")
