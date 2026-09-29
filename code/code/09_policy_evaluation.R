# =====================================================================
# 09_policy_evaluation.R
# Builds the policy-tree covariate set, fits cross-validated policy trees
# (depth 2-4) and the plug-in rule, and compares every candidate policy's
# out-of-sample value in a pairwise table. Also fits the final,
# full-sample deployable trees used later for the tree plots.
#
# Produces: Table 2 (Estimated reward of the policy rules compared to
# alternative approaches)
#   -> policy_comparison_table.html, policy_comparison_table.tex,
#      policy_comparison_long_holm.csv
#   Also produces final_trees (d2/d3/d4), consumed by
#   10_translate_tree_splits.R and 11_plot_policy_trees.R.
# Depends on: 02_prepare_data.R (X, n, arms, active_arms, arm_lbl)
#             08_policy_scoring_targeting_forests.R (Gamma, fold_id,
#               forest_targeting_restricted)
# =====================================================================

# =============================================================================
# POLICY EVALUATION: 5-FOLD CROSS-FITTED CAIPWL (Zhou et al. 2023),
# COMPARISON TABLE (Ahrens Table 2 style)
# =============================================================================
# Gamma (doubly-robust/AIPW scores) comes from the fold-clustered forest
# fit in 08_policy_scoring_targeting_forests.R, so scores for fold k use no
# fold-k outcomes. The 5-fold split trains each policy tree on 4/5 of the
# Gamma rows and evaluates it on the held-out 1/5.

# --- Covariates for the decision rules ---------------------------------------
# I use a reduced, low-correlation covariate set (dropping avg_risk_factor:
# r = 0.949 with age_mean, and hh_risk_factor: r = 0.689 with hh_members) to
# cut exact-search runtime across all tree depths, and I exclude sensitive
# attributes (curr_gender_num1, race/ethnicity) from the decision rules --
# nuisance models behind Gamma may use them, but the policy may not.
#
# fpl_imp (~51k distinct values, well above the package's 20k exact-search
# warning threshold), indv_max_aptc_amt (~1.6k values), and age_mean (~513
# values) dominate policy_tree runtime relative to the other covariates. I
# discretize all three to ordinal bin indices (still single numeric columns
# -- policy_tree only ever makes threshold splits, so this is not a
# one-hot/dummy expansion) JUST for the policy tree's decision-rule inputs.
# The causal forest / Gamma above still use the raw continuous variables.
#
# I use quantile bins with duplicate breakpoints dropped, rather than
# ntile()'s equal-COUNT rank bins: ntile(x, 100) wastes candidates on mass
# points (e.g. indv_max_aptc_amt has ~44% of households at exactly $0),
# producing many bins that are literally identical thresholds once
# split.step downsamples them. unique() on the quantile breaks collapses
# those duplicates so every remaining bin is a genuinely distinct
# real-world cutoff.
make_policy_bin <- function(x, n_bins = 20) {
  breaks <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  as.integer(cut(x, breaks = breaks, include.lowest = TRUE))
}

X_policy_df <- as.data.frame(X)

# policy_tree() cannot handle NaN covariates (unlike grf's causal forests,
# which use MIA). I impute missing FPL to 0 ONLY for this bin -- fpl_missing
# stays available as its own covariate, so the tree can still isolate that
# subpopulation explicitly rather than silently merging it into the low end
# of the income distribution.
X_policy_df$fpl_for_bin <- ifelse(is.na(X_policy_df$fpl_imp), 0, X_policy_df$fpl_imp)
X_policy_df$fpl_bin  <- make_policy_bin(X_policy_df$fpl_for_bin)
X_policy_df$aptc_bin <- make_policy_bin(X_policy_df$indv_max_aptc_amt)
X_policy_df$age_bin  <- make_policy_bin(X_policy_df$age_mean)

policy_covs <- c("age_bin", "fpl_bin", "aptc_bin", "hoh_email_ok",
                 "sws", "hh_members", "English", "fpl_missing", "fpl_1000p",
                 "region3", "region18", "region15")
stopifnot(all(policy_covs %in% colnames(X_policy_df)))
X_policy <- as.matrix(X_policy_df[, policy_covs])

# =============================================================================
# Cross-fitted trees: (depth, fold) grid, fit in parallel
# =============================================================================
# I cache this to disk -- it's the expensive step (15 policy_tree/
# hybrid_policy_tree fits, depth 2-4 x 5 folds) and doesn't need to be redone
# just because the session restarted or a later chunk broke.
#
# IMPORTANT: the cache filename encodes the covariate set. If policy_covs,
# SPLIT, or the fold assignment changes, the name below changes with it (via
# the digest of policy_covs), so a stale cache can never be silently reused
# against a different X_policy.
#
# Depth 4 uses hybrid_policy_tree() instead of the exact policy_tree(): an
# exact search over 4 levels of splits is combinatorially infeasible at this
# n and covariate count. hybrid_policy_tree does an exact search down to
# search.depth levels, then a greedy search for the remaining levels -- an
# approximation, so the d=4 tree is not guaranteed globally optimal the way
# the d=2/3 exact trees are.
handlers(handler_cli(format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} ({cli::pb_percent}) | eta: {cli::pb_eta}"))
plan(multisession, workers = min(8, parallel::detectCores() - 1))

cache_tag <- substr(digest::digest(list(policy_covs, SPLIT, SEARCH_DEPTH, fold_id)), 1, 8)
TREE_FITS_CACHE   <- paste0(OUT_DIR, "tree_fits_cache_d2to4_",   cache_tag, ".rds")
FINAL_TREES_CACHE <- paste0(OUT_DIR, "final_trees_cache_d2to4_", cache_tag, ".rds")

if (file.exists(TREE_FITS_CACHE)) {
  cat("Loading cached tree_fits from", TREE_FITS_CACHE, "\n")
  tree_fits <- readRDS(TREE_FITS_CACHE)
} else {
  grid <- expand.grid(d = 2:4, k = 1:5)

  with_progress({
    p <- progressor(steps = nrow(grid))
    tree_fits <- future_lapply(seq_len(nrow(grid)), function(g) {
      d <- grid$d[g]; k <- grid$k[g]
      tr <- which(fold_id != k); te <- which(fold_id == k)

      tree_k <- if (d <= 3) {
        policy_tree(X_policy[tr, , drop = FALSE], Gamma[tr, , drop = FALSE],
                    depth = d, split.step = SPLIT)
      } else {
        hybrid_policy_tree(X_policy[tr, , drop = FALSE], Gamma[tr, , drop = FALSE],
                           depth = d, search.depth = SEARCH_DEPTH, split.step = SPLIT)
      }

      p(sprintf("done depth=%d fold=%d", d, k))
      list(d = d, k = k, te = te, a = predict(tree_k, X_policy[te, , drop = FALSE]), tree = tree_k)
    }, future.seed = TRUE)
  })

  saveRDS(tree_fits, TREE_FITS_CACHE)
}

# =============================================================================
# Per-unit AIPW value under each policy
# =============================================================================
policy_values <- list()

## One-size-fits-all (incl. "Nothing" = always Control)
policy_values[["Nothing"]] <- Gamma[, "Control"]
for (a in active_arms)
  policy_values[[paste0("Always: ", arm_lbl(a))]] <- Gamma[, a]

## Policy trees d = 2, 3, 4 (cross-fitted, non-sensitive covariates)
for (d in 2:4) {
  vals <- rep(NA_real_, n)
  for (f in Filter(function(x) x$d == d, tree_fits))
    vals[f$te] <- Gamma[cbind(f$te, f$a)]
  stopifnot(!anyNA(vals))          # every unit scored exactly once, out-of-fold
  policy_values[[paste0("Policy tree (d=", d, ")")]] <- vals
}

## Plug-in rule: arm with the highest OOB-predicted CATE, Control if all < 0,
## from the coarse targeting forest (AKS min.node.size = 300, Appendix C.2).
## Predictions are OOB under the fold clustering, so each unit's assignment
## comes from trees that excluded its fold -- the same held-out standard as
## the tree rows.
make_plugin_values <- function(forest_obj) {
  tp <- forest_obj$predictions
  if (length(dim(tp)) == 3) tp <- tp[, , 1, drop = TRUE]
  stopifnot(nrow(tp) == n, ncol(tp) == length(active_arms))
  # column order of tp follows the factor levels of W, minus the base level
  stopifnot(identical(levels(W)[-1], active_arms))

  best_arm <- max.col(tp)
  best_tau <- tp[cbind(1:n, best_arm)]
  action   <- ifelse(best_tau > 0, active_arms[best_arm], "Control")
  Gamma[cbind(1:n, match(action, arms))]
}
policy_values[["Plug-in rule (non-sensitive)"]] <- make_plugin_values(forest_targeting_restricted)

# =============================================================================
# Pairwise paired comparisons (column beats row if > 0)
# =============================================================================
# NOTE: these p-values are UNADJUSTED across the K(K-1)/2 cells of the
# table, unlike the Holm-adjusted blocks elsewhere in this document. I write
# a Holm-adjusted copy alongside so the write-up can cite either; the
# caption states which is displayed.
pv_names <- names(policy_values)
K <- length(pv_names)
cell <- matrix("", K, K, dimnames = list(pv_names, pv_names))
stars <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))

comp_rows <- list()
for (i in 1:(K - 1)) for (j in (i + 1):K) {
  dif <- policy_values[[j]] - policy_values[[i]]
  est <- mean(dif) * 100                      # percentage points
  se  <- sd(dif) / sqrt(n) * 100
  p   <- 2 * pnorm(-abs(est / se))
  cell[i, j] <- sprintf("%.3f%s (%.3f)", est, stars(p), se)
  comp_rows[[length(comp_rows) + 1]] <- data.frame(
    row_policy = pv_names[i], col_policy = pv_names[j], est_pp = est, se_pp = se, p = p
  )
}
comp_long <- bind_rows(comp_rows) %>% mutate(p_holm = p.adjust(p, method = "holm"))
write.csv(comp_long, paste0(OUT_DIR, "policy_comparison_long_holm.csv"), row.names = FALSE)

comp_df <- as.data.frame(cell[-K, -1, drop = FALSE])
print(comp_df)

comp_caption <- paste(
  "Estimated reward of column policy relative to row policy, in pp of take-up.",
  "Paired AIPW differences; SEs in parentheses. Stars are unadjusted;",
  "Holm-adjusted p-values for all pairs are in policy_comparison_long_holm.csv.",
  "* p<0.1, ** p<0.05, *** p<0.01."
)
comp_df %>%
  kbl(align = "c", caption = comp_caption) %>%
  kable_classic(full_width = FALSE) %>%
  save_kable(paste0(OUT_DIR, "policy_comparison_table.html"))
comp_df %>%
  kbl(format = "latex", booktabs = TRUE, align = "c", caption = comp_caption) %>%
  save_kable(paste0(OUT_DIR, "policy_comparison_table.tex"))

# =============================================================================
# Final deployable trees (fit once on ALL data; report these), and headline
# numbers vs. doing nothing
# =============================================================================
if (file.exists(FINAL_TREES_CACHE)) {
  cat("Loading cached final_trees from", FINAL_TREES_CACHE, "\n")
  final_trees <- readRDS(FINAL_TREES_CACHE)
} else {
  with_progress({
    p <- progressor(steps = 3)
    final_trees <- future_lapply(2:4, function(d) {
      ft <- if (d <= 3) {
        policy_tree(X_policy, Gamma, depth = d, split.step = SPLIT)
      } else {
        hybrid_policy_tree(X_policy, Gamma, depth = d, search.depth = SEARCH_DEPTH, split.step = SPLIT)
      }
      p(sprintf("final tree d=%d done", d))
      ft
    }, future.seed = TRUE)
  })
  names(final_trees) <- paste0("d", 2:4)
  saveRDS(final_trees, FINAL_TREES_CACHE)
}

for (d in 2:4) {
  nm <- paste0("Policy tree (d=", d, ")")
  dif <- policy_values[[nm]] - policy_values[["Nothing"]]
  cat(nm, "vs Nothing:", sprintf("%.3f pp (SE %.3f)\n", mean(dif) * 100, sd(dif) / sqrt(n) * 100))
}
