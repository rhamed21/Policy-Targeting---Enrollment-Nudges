# =====================================================================
# 01_setup.R
# Loads packages, fixes the seed, and defines constants used throughout
# every later script (grf/policytree tuning parameters, plotting scale).
#
# Produces: nothing (no figures/tables) -- run first, always.
# Depends on: nothing
# =====================================================================

# ---- Packages ---------------------------------------------------------
library(policytree)
library(grf)
library(dplyr)
library(tidyr)
library(haven)
library(sandwich)
library(lmtest)
library(kableExtra)
library(ggplot2)
library(scales)
library(cowplot)
library(patchwork)
library(forcats)
library(future.apply)
library(progressr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ---- Seed ---------------------------------------------------------------
set.seed(21)

# ---- grf / policytree tuning ---------------------------------------------
# Used by every forest and every policy_tree() call from
# "04_heterogeneity_forest.R" through "11_plot_policy_trees.R".
N_TREES <- 3000
# policy_tree()/hybrid_policy_tree() split.step -- with quantile-deduplicated
# bins (see make_policy_bin() in "09_policy_evaluation.R"), every bin is
# already a distinct, meaningful threshold, so an exact search (no
# downsampling) is cheap.
SPLIT        <- 1
SEARCH_DEPTH <- 2   # how far hybrid_policy_tree() looks ahead when splitting

# Percentage-point scale factor, used by every CATE figure (subgroup panels,
# density plot) to convert model estimates into the units the write-up
# reports in.
PP_SCALE <- 100

# NOTE: DATA_PATH and OUT_DIR are set once in master.R (not here), so every
# script below can be reused unchanged on any machine -- just point
# master.R at your files and source it.
