# =====================================================================
# 08_policy_scoring_targeting_forests.R
# Fits the two nuisance forests policy learning depends on: the
# cross-fitted scoring forest that supplies the doubly-robust reward
# matrix Gamma, and a coarser targeting forest (fairness-restricted,
# non-sensitive covariates only) used solely by the plug-in rule.
#
# Produces: nothing directly (no figures/tables). Gamma and
# forest_targeting_restricted feed 09_policy_evaluation.R (Table 2).
# Depends on: 02_prepare_data.R (X, Y, W, W.hat, all_covs, n)
# =====================================================================

# -----------------------------------------------------------------------------
# STEP 1: SCORING FOREST
# -----------------------------------------------------------------------------
# `forest` supplies Gamma, the AIPW reward matrix every policy in the
# comparison table is scored against. This is a NUISANCE object, so it runs
# at grf's default min.node.size = 5 -- AKS likewise use the default node
# size for nuisance components (Appendix C.1, C.2) and reserve coarsening
# for targeting.
#
# clusters = fold_id implements cross-fitting: grf samples whole clusters, so
# a tree contains either all of fold k or none of it, and the OOB score for a
# fold-k unit therefore comes only from trees that never saw fold k. This is
# AKS's own implementation of cross-fitting in grf (Appendix C.1: "We
# implement this procedure by providing the folds as clusters to the grf
# algorithm").
fold_id <- sample(rep(1:5, length.out = n))

forest <- multi_arm_causal_forest(
  X, Y, W,
  W.hat     = W.hat,
  num.trees = N_TREES,
  clusters  = fold_id,
  seed      = 21
)
Gamma <- double_robust_scores(forest)

# -----------------------------------------------------------------------------
# STEP 2: TARGETING FOREST FOR THE PLUG-IN RULE
# -----------------------------------------------------------------------------
# I use an AKS-style coarse forest (min.node.size = 300, "to avoid excess
# noise in targeting"), clustered on fold_id so each unit's OOB prediction
# comes only from trees that never saw its fold -- the same held-out
# standard as the tree rows, and matching Ahrens, who predicts onto his test
# split before taking the argmax. Its predictions drive the plug-in argmax
# rule and NOTHING else -- it never touches Gamma.
#
# I fit this forest on the NON-SENSITIVE covariates only, so the plug-in row
# in the comparison table is the fairness-constrained ceiling: plug-in vs.
# tree is then the cost of the tree's shallow, interpretable functional form
# (Ahrens' "85% of the plug-in rule" comparison).
#
# It uses RAW CONTINUOUS covariates. The 20-bin discretization used by the
# policy tree below exists only to make policy_tree's exact search
# tractable; a forest has no such constraint, so binning here would
# handicap the plug-in rule for no reason.
sensitive_vars <- c("curr_gender_num1", "latino", "black", "white", "asian", "otherrace")
plugin_covs_restricted <- setdiff(all_covs, sensitive_vars)
X_plugin_restricted    <- X[, plugin_covs_restricted, drop = FALSE]

forest_targeting_restricted <- multi_arm_causal_forest(
  X_plugin_restricted, Y, W,
  W.hat         = W.hat,
  num.trees     = N_TREES,
  clusters      = fold_id,
  min.node.size = 300,
  seed          = 21
)
