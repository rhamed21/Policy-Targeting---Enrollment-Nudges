# =====================================================================
# 06_cate_distribution.R
# Kernel-density plot of household-level predicted treatment effects for
# each letter, with each arm's ATE marked and annotated.
#
# Produces: Figure 1 (Distribution of household-level predicted treatment
# effects of each letter, using multi-arm causal forests)
#   -> fig2_cate_plot_only.pdf, fig2_cate_legend_only.pdf
# Depends on: 04_heterogeneity_forest.R (forest_est, tau_hat)
#             02_prepare_data.R (active_arms, arm_cols, arm_lbl)
# =====================================================================

# =============================================================================
# CATE DISTRIBUTION -- plot and legend saved separately
# =============================================================================
# The legend text already carries the display arm names (via arm_lbl) plus
# each arm's ATE, so this plot uses ONE colour scale only -- scale_colour_manual.
# Adding scale_color_discrete() on top would override both the arm colours
# and the ATE annotations.
cate_pred_raw <- forest_est$predictions
if (length(dim(cate_pred_raw)) == 3) cate_pred_raw <- cate_pred_raw[, , 1, drop = TRUE]
cate_long <- as.data.frame(cate_pred_raw) %>%
  setNames(active_arms) %>%
  pivot_longer(everything(), names_to = "arm", values_to = "cate") %>%
  mutate(arm = factor(arm, levels = active_arms), cate_pp = cate * PP_SCALE)

ate_lines <- tau_hat %>%
  transmute(arm = factor(contrast, levels = active_arms),
            ate = estimate * PP_SCALE,
            se  = std.err * PP_SCALE)

ate_lines$y_end <- vapply(seq_len(nrow(ate_lines)), function(i) {
  a <- as.character(ate_lines$arm[i])
  d <- density(cate_long$cate_pp[cate_long$arm == a], kernel = "epanechnikov")
  approx(d$x, d$y, xout = ate_lines$ate[i])$y
}, numeric(1))
stopifnot(all(c("ate", "se", "y_end") %in% names(ate_lines)), !anyNA(ate_lines$y_end))

arm_labs <- setNames(
  sprintf("%s; ATE = %.2f (%.2f)", arm_lbl(as.character(ate_lines$arm)), ate_lines$ate, ate_lines$se),
  as.character(ate_lines$arm)
)[active_arms]

rng <- range(cate_long$cate_pp)
brk <- seq(floor(rng[1] / 2) * 2, ceiling(rng[2] / 2) * 2, by = 2)

# --- Callout geometry (the "ATEs" bracket over the vertical ATE lines) -------
y_max <- max(density(cate_long$cate_pp, kernel = "epanechnikov")$y)
y_upper_limit <- ceiling(y_max * 10) / 10

ell_x  <- mean(range(ate_lines$ate))
ell_y  <- 0.17 * y_max
ell_rx <- 0.55 * diff(range(ate_lines$ate)) + 0.35
ell_ry <- 0.035 * y_max
tt  <- seq(0, 2 * pi, length.out = 100)
ell <- data.frame(x = ell_x + ell_rx * cos(tt), y = ell_y + ell_ry * sin(tt))
lab_x <- ell_x + 4.5
lab_y <- ell_y + 0.09 * y_max

# --- Base plot (legend on, so it can be harvested) ---------------------------
p_base <- ggplot(cate_long, aes(x = cate_pp, colour = arm)) +
  geom_density(kernel = "epanechnikov", linewidth = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(data = ate_lines,
               aes(x = ate, xend = ate, y = 0, yend = y_end, colour = arm),
               linetype = "dashed", linewidth = 0.9,
               inherit.aes = FALSE, show.legend = FALSE) +
  annotate("path", x = ell$x, y = ell$y, colour = "black", linewidth = 0.5) +
  annotate("segment", x = ell_x + ell_rx, y = ell_y + ell_ry * 0.5,
           xend = lab_x - 0.4, yend = lab_y, colour = "black", linewidth = 0.5) +
  annotate("text", x = lab_x, y = lab_y, label = "ATEs", hjust = 0, fontface = "bold", size = 5.5) +
  scale_colour_manual(values = arm_cols, labels = arm_labs, name = "Arm") +
  scale_x_continuous(breaks = brk, labels = scales::number_format(accuracy = 1),
                     minor_breaks = NULL, expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, y_upper_limit),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 18) +
  theme(
    legend.title     = element_text(size = 20),
    legend.text      = element_text(size = 17),
    legend.key       = element_blank(),
    legend.key.width = unit(1.4, "lines"),
    axis.text        = element_text(size = 16, colour = "black"),
    axis.title       = element_text(size = 20)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.6))) +
  labs(x = "Effect of treatment letter on enrollment rate (percentage points)", y = "Density")

# --- 1. Plot only ------------------------------------------------------------
p_only <- p_base + theme(legend.position = "none")
ggsave(paste0(OUT_DIR, "fig2_cate_plot_only.pdf"),
       plot = p_only, width = 26, height = 16, unit = "cm", dpi = 300)

# --- 2. Legend only ----------------------------------------------------------
p_with_leg <- p_base + theme(legend.position = "right")
g <- ggplot2::ggplotGrob(p_with_leg)
# I grep for "guide-box" because it catches both the modern "guide-box" and
# the legacy "guide-box-right" grob names across ggplot2 versions.
leg_idx <- grep("guide-box", g$layout$name)
p_leg <- if (length(leg_idx) > 0) g$grobs[[leg_idx[1]]] else NULL
stopifnot(!is.null(p_leg))

ggsave(paste0(OUT_DIR, "fig2_cate_legend_only.pdf"),
       plot = cowplot::ggdraw(p_leg),
       width = 14, height = 6, unit = "cm", dpi = 300, bg = "white")
p_only
