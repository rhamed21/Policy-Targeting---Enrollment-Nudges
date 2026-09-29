# =====================================================================
# 04_heterogeneity_forest.R
# Fits the descriptive multi-arm causal forest used for CATE plots and
# subgroup difference tests, and reports pooled (cross-arm) variable
# importance.
#
# Produces: no standalone figure/table here.
#   - forest_est / tau_hat feed Figure 1 (06) and Figure 2 + Appendix A2 (05)
#   - var_imp feeds Appendix Figure A1, Panel A (07)
# Depends on: 02_prepare_data.R (X, Y, W, W.hat, active_arms, lbl)
# =====================================================================

# =============================================================================
# ESTIMATION FOREST (for CATE plots + difference tests only)
# =============================================================================
# I fit a default-tuned forest for descriptive heterogeneity, which gives
# more precise DR scores / subgroup ATEs than the coarse targeting forest
# used later for policy assignment.
# NOTE: I fit on X (fpl_imp, no NAs), not a version that swaps in raw
# subsidy_fpl_percent (NA for the 8,885 fpl_missing==1 households), purely
# so make_cate_plot()/the subgroup tables below can bin on the real FPL
# scale. The forest itself must see the same clean covariates as the
# targeting forest further down, or its splits/tau_hat silently diverge.
forest_est <- multi_arm_causal_forest(
  X, Y, W,
  W.hat     = W.hat,
  num.trees = N_TREES,
  seed      = 21          # default min.node.size = 5
)

tau_hat <- average_treatment_effect(forest_est)
stopifnot(identical(sub(" - Control$", "", tau_hat$contrast), active_arms))
tau_hat$contrast <- factor(active_arms, levels = active_arms)

# =============================================================================
# Pooled variable importance
# =============================================================================
var_imp <- variable_importance(forest_est)
rownames(var_imp) <- lbl(colnames(X))
top10_vars <- rownames(var_imp)[order(var_imp, decreasing = TRUE)][1:10]
cat("\nTop 10 variables by importance:\n")
print(round(sort(var_imp[top10_vars, ], decreasing = TRUE), 4))
