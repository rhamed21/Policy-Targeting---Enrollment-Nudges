# =====================================================================
# 07_variable_importance_by_arm.R
# Fits one binary causal forest per active arm (vs. Control) to get
# arm-specific variable importance (Mulungu et al. 2024 JAE, Sec 3.3 /
# Eq. 6 style), since the pooled multi-arm forest only yields one score
# per covariate.
#
# Produces: Appendix Figure A1 (pooled top-10 + by-arm top-10 variable
# importance)
#   -> variable_importance_pooled_top10.pdf, variable_importance_by_arm_labeled.pdf,
#      variable_importance_by_arm.csv
# Depends on: 02_prepare_data.R (X, Y, df, all_covs, lbl)
#             04_heterogeneity_forest.R (var_imp, pooled importance)
# =====================================================================

# =============================================================================
# PER-ARM VARIABLE IMPORTANCE (Mulungu et al. 2024 JAE, Sec 3.3 / Eq. 6 style)
# =============================================================================
# The pooled multi_arm_causal_forest (forest_est, from 04) splits jointly
# across all K-1 arm-contrasts, so variable_importance() on it yields ONE
# importance score per covariate, not one per arm. To get arm-specific
# importance (as in Mulungu et al., who fit a separate causal forest per
# treatment vs. control), I instead fit K-1 independent BINARY causal_forest
# models, one per active arm, each restricted to {Control, that arm}.
# Propensity within each subsample is exactly 0.5 by the RCT's
# equal-allocation design.
imp_list       <- list()
forests_by_arm <- list()

for (a in active_arms) {
  keep <- df$W %in% c("Control", a)
  Xa <- X[keep, ]
  Ya <- Y[keep]
  Wa <- as.integer(df$W[keep] == a)   # 1 = arm, 0 = Control

  cf_a <- causal_forest(
    Xa, Ya, Wa,
    W.hat     = rep(0.5, length(Wa)),
    num.trees = N_TREES,
    seed      = 21
  )
  forests_by_arm[[a]] <- cf_a

  vi_a <- variable_importance(cf_a)
  rownames(vi_a) <- colnames(Xa)
  imp_list[[a]] <- vi_a[, 1]
}

imp_df <- as.data.frame(imp_list, check.names = FALSE)
imp_df$covariate <- lbl(all_covs)
imp_df <- imp_df[, c("covariate", active_arms)]
stopifnot(!any(duplicated(imp_df$covariate)))   # display names must be unique

# CSV keeps arm keys as column names; a display-name header row can be added
# separately for the write-up if needed.
write.csv(imp_df, paste0(OUT_DIR, "variable_importance_by_arm.csv"), row.names = FALSE)

cat("\nTop 8 variables by importance, per arm:\n")
for (a in active_arms) {
  cat("\n---", arm_lbl(a), "vs", arm_lbl("Control"), "---\n")
  top <- imp_df[order(-imp_df[[a]]), c("covariate", a)][1:8, ]
  print(top, row.names = FALSE)
}

# =============================================================================
# Variable importance figures: pooled top-10 and by-arm top-10
# =============================================================================
N_TOP_VARS <- 10

# --- Panel A: pooled variable importance (single grey bars) -----------------
var_imp_df <- tibble(covariate = rownames(var_imp), importance = var_imp[, 1]) %>%
  slice_max(order_by = importance, n = N_TOP_VARS) %>%
  arrange(importance) %>%
  mutate(covariate = factor(covariate, levels = covariate))

p_var_imp_pooled <- ggplot(var_imp_df, aes(x = covariate, y = importance)) +
  geom_col(fill = "grey45") +
  coord_flip() +
  theme_classic(base_size = 13) +
  theme(text = element_text(size = 13)) +
  labs(x = NULL, y = "Variable importance")

ggsave(paste0(OUT_DIR, "variable_importance_pooled_top10.pdf"),
       plot = p_var_imp_pooled, width = 18, height = 14, unit = "cm", dpi = 300)

# --- Panel B: by-arm variable importance (top N per arm, arm-coloured) ------
imp_long <- imp_df %>%
  pivot_longer(cols = all_of(active_arms), names_to = "arm", values_to = "importance") %>%
  mutate(arm = factor(arm, levels = active_arms))

top_per_arm <- imp_long %>%
  group_by(arm) %>%
  slice_max(order_by = importance, n = N_TOP_VARS) %>%
  ungroup() %>%
  # I encode arm into the factor key so each facet's bars sort
  # independently ("reorder within facet"), then strip the suffix for
  # display in scale_x_discrete() below.
  mutate(cov_key = paste(covariate, arm, sep = "___")) %>%
  arrange(arm, importance) %>%
  mutate(cov_key = factor(cov_key, levels = cov_key))

p_var_imp_by_arm <- ggplot(top_per_arm, aes(x = cov_key, y = importance, fill = arm)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = arm_cols) +
  scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
  coord_flip() +
  facet_wrap(~ arm, labeller = as_labeller(arm_labels), scales = "free_y", ncol = 2) +
  theme_classic(base_size = 13) +
  theme(text = element_text(size = 13), strip.text = element_text(face = "bold")) +
  labs(x = NULL, y = "Variable importance")

ggsave(paste0(OUT_DIR, "variable_importance_by_arm_labeled.pdf"),
       plot = p_var_imp_by_arm, width = 18, height = 14, unit = "cm", dpi = 300)

p_var_imp_by_arm   # preview to confirm the arm-coloured facets before moving on
