# =====================================================================
# 10_translate_tree_splits.R
# make_policy_bin() splits report the bin index (<= 20), not the raw
# dollar / percent-of-FPL / age value. This builds a bin -> real-value
# lookup for each binned covariate and prints each final tree with splits
# translated into real-world units.
#
# Produces: no figures/tables -- console output only (used to write the
# "Learning the Optimal Policy Rules" narrative, e.g. "Subsidy <= $129/mth").
# Depends on: 09_policy_evaluation.R (X_policy_df, final_trees)
#             02_prepare_data.R (arm_lbl, lbl)
# =====================================================================

# I build a bin -> real-value lookup for each binned covariate
# (lookup[[var]][k] = max of the raw variable among obs in bin <= k, i.e.
# the actual real-world cutoff a printed "var <= k" split corresponds to)
# and a pretty-printer that substitutes it into the tree output.
make_bin_lookup <- function(bin, raw, n_bins = max(bin, na.rm = TRUE)) {
  by_bin <- tapply(raw, bin, max)
  full <- rep(NA_real_, n_bins)
  full[as.integer(names(by_bin))] <- by_bin
  # I fill empty bins forward; seq_len() rather than 2:n_bins keeps this
  # safe when n_bins == 1.
  for (k in seq_len(n_bins)[-1]) if (is.na(full[k])) full[k] <- full[k - 1]
  cummax(full)  # guard against any non-monotonicity
}

bin_lookup <- list(
  fpl_bin  = make_bin_lookup(X_policy_df$fpl_bin,  X_policy_df$fpl_for_bin),
  aptc_bin = make_bin_lookup(X_policy_df$aptc_bin, X_policy_df$indv_max_aptc_amt),
  age_bin  = make_bin_lookup(X_policy_df$age_bin,  X_policy_df$age_mean)
)
bin_units <- c(fpl_bin = "% of FPL", aptc_bin = "$/month subsidy", age_bin = "years")

print_translated_tree <- function(tree, bin_lookup, bin_units, node_id = 1, depth = 0) {
  node <- tree$nodes[[node_id]]
  indent <- strrep("  ", depth)
  if (isTRUE(node$is_leaf)) {
    cat(indent, "* action:", arm_lbl(tree$action.names[node$action]), "\n")
    return(invisible())
  }
  var <- tree$columns[node$split_variable]
  val <- node$split_value
  if (var %in% names(bin_lookup)) {
    real_val <- bin_lookup[[var]][val]
    cat(sprintf("%s[%s <= bin %d]  ->  %s <= %s %s\n",
                indent, var, val, lbl(var), format(round(real_val, 1), big.mark = ","), bin_units[[var]]))
  } else {
    cat(sprintf("%s[%s <= %s]\n", indent, lbl(var), val))
  }
  print_translated_tree(tree, bin_lookup, bin_units, node$left_child, depth + 1)
  print_translated_tree(tree, bin_lookup, bin_units, node$right_child, depth + 1)
}

for (nm in names(final_trees)) {
  cat("\n===", nm, "===\n")
  print_translated_tree(final_trees[[nm]], bin_lookup, bin_units)
}
