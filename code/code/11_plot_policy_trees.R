# =====================================================================
# 11_plot_policy_trees.R
# Renders each final policy tree (depth 2-4) as a graphviz diagram with
# arm-matched leaf colours and human-readable split labels (dollars,
# % FPL, years, True/False for indicators), then exports as vector PDFs.
#
# Produces: Figure 3 (Fitted depth-2 policy tree), plus the depth-3 and
# depth-4 trees discussed in the text / appendix
#   -> policy_tree_depth2_5arms.pdf, policy_tree_depth3_5arms.pdf,
#      policy_tree_depth4_5arms.pdf
# Depends on: 09_policy_evaluation.R (final_trees)
#             10_translate_tree_splits.R (bin_lookup)
#             02_prepare_data.R (arm_cols, arm_labels)
# =====================================================================

# =============================================================================
# Plot policy trees with (a) arm-matched leaf colours and (b) custom node text
# =============================================================================
# --- 1. Leaf colours and labels ----------------------------------------------
# I reuse the same hue_pal() colours (`arm_cols`, defined once in
# 02_prepare_data.R) that the CATE density / CATE pointrange figures use, so
# leaf colour is consistent across the whole document. Control gets a
# neutral grey.
leaf_fill  <- c(Control = "#BDBDBD", arm_cols)
leaf_label <- arm_labels      # single source of truth, same as every figure

# --- 2. Node text ------------------------------------------------------------
# For each policy covariate, I define how to render "var <= value":
#   fmt  : function(value_on_raw_scale) -> label text
#   flip : TRUE for 0/1 indicators phrased as a positive statement, where the
#          <= 0 branch is the FALSE branch, so the True/False edge labels swap.
#   bin  : name in bin_lookup, if the column is a quantile-bin index
node_spec <- list(
  age_bin      = list(bin = "age_bin",
                      fmt = function(v) sprintf("%s <= %s yrs", lbl("age_bin"), format(round(v), big.mark = ","))),
  aptc_bin     = list(bin = "aptc_bin",
                      fmt = function(v) sprintf("%s <= $%s/mth", lbl("aptc_bin"), format(round(v), big.mark = ","))),
  fpl_bin      = list(bin = "fpl_bin", fmt = function(v) sprintf("%s <= %.0f%%", lbl("fpl_bin"), v)),
  hh_members   = list(fmt = function(v) sprintf("%s <= %g", lbl("hh_members"), v)),
  English      = list(flip = TRUE, fmt = function(v) "Language: English"),
  sws          = list(flip = TRUE, fmt = function(v) "Funnel: County Referrals"),
  hoh_email_ok = list(flip = TRUE, fmt = function(v) lbl("hoh_email_ok")),
  fpl_missing  = list(flip = TRUE, fmt = function(v) lbl("fpl_missing")),
  fpl_1000p    = list(flip = TRUE, fmt = function(v) lbl("fpl_1000p"))
)
for (r in 1:19)
  node_spec[[paste0("region", r)]] <-
    list(flip = TRUE, fmt = local({ rr <- r; function(v) lbl(paste0("region", rr)) }))

# Every covariate the tree can split on must have a rendering rule, so
# nothing ever falls through to the raw-variable-name fallback in
# describe_split().
stopifnot(all(policy_covs %in% names(node_spec)))

describe_split <- function(var, val, bin_lookup) {
  sp <- node_spec[[var]]
  if (is.null(sp))
    return(list(label = sprintf("%s <= %g", lbl(var), val), flip = FALSE))
  v <- if (!is.null(sp$bin)) bin_lookup[[sp$bin]][val] else val
  list(label = sp$fmt(v), flip = isTRUE(sp$flip))
}

# =============================================================================
# DOT builder, then render + export each depth's tree
# =============================================================================
esc <- function(s) gsub('"', '\\\\"', s)

policy_tree_dot <- function(tree, bin_lookup, leaf_fill, leaf_label,
                            fontname = "Helvetica", fontsize = 12) {
  nodes <- character(0); edges <- character(0)

  walk <- function(id) {
    nd <- tree$nodes[[id]]
    if (isTRUE(nd$is_leaf)) {
      arm <- tree$action.names[nd$action]
      lab <- if (arm %in% names(leaf_label)) leaf_label[[arm]] else arm
      col <- if (arm %in% names(leaf_fill))  leaf_fill[[arm]]  else "#EEEEEE"
      nodes <<- c(nodes, sprintf(
        '  n%d [label="%s", shape=box, style="filled,rounded", fillcolor="%s", color="%s", penwidth=1.2];',
        id, esc(lab), col, col))
      return(invisible())
    }
    var <- tree$columns[nd$split_variable]
    d   <- describe_split(var, nd$split_value, bin_lookup)
    nodes <<- c(nodes, sprintf(
      '  n%d [label="%s", shape=box, style="rounded", color="#444444", penwidth=1.2];',
      id, esc(d$label)))

    # The tree's own left child is the "x <= value" branch. For flipped
    # nodes (0/1 indicators phrased positively) that branch is the FALSE
    # one, so I emit the children in swapped order -- graphviz lays
    # siblings out in declaration order, which keeps True on the left in
    # every node.
    true_child  <- if (d$flip) nd$right_child else nd$left_child
    false_child <- if (d$flip) nd$left_child  else nd$right_child

    edges <<- c(edges,
      sprintf('  n%d -> n%d [label="  True  "];',  id, true_child),
      sprintf('  n%d -> n%d [label="  False  "];', id, false_child))
    walk(true_child); walk(false_child)
  }
  walk(1)

  paste0(
    'digraph policy_tree {\n',
    sprintf('  graph [nodesep=0.35, ranksep=0.55, fontname="%s"];\n', fontname),
    sprintf('  node  [fontname="%s", fontsize=%d, fontcolor="#000000", margin="0.14,0.07"];\n', fontname, fontsize),
    sprintf('  edge  [fontname="%s", fontsize=%d, color="#666666", arrowsize=0.7];\n', fontname, max(9, fontsize - 2)),
    paste(nodes, collapse = "\n"), "\n",
    paste(edges, collapse = "\n"), "\n}"
  )
}

for (d in 2:4) {
  nm  <- paste0("d", d)
  dot <- policy_tree_dot(final_trees[[nm]], bin_lookup, leaf_fill, leaf_label,
                         fontsize = if (d == 4) 11 else 12)
  g   <- DiagrammeR::grViz(dot)
  svg_code <- DiagrammeRsvg::export_svg(g)

  pdf_path <- paste0(OUT_DIR, "policy_tree_depth", d, "_5arms.pdf")
  # rsvg_pdf() writes vector PDF directly from the SVG, so (unlike rsvg_png())
  # no pixel width/dpi is needed -- the output stays crisp at any zoom.
  rsvg::rsvg_pdf(charToRaw(svg_code), file = pdf_path)
  cat("wrote", pdf_path, "\n")
}

# I preview the depth-2 tree (Figure 3 in the paper).
DiagrammeR::grViz(policy_tree_dot(final_trees$d2, bin_lookup, leaf_fill, leaf_label))
