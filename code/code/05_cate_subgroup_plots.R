# =====================================================================
# 05_cate_subgroup_plots.R
# Builds CATE-by-subgroup panels for every covariate group of interest.
#
# Produces:
#   - Figure 2 (By income / By subsidy / By household age / By funnel round)
#     -> cate_fpl.pdf, cate_aptc.pdf, cate_age.pdf, cate_funnel_sws.pdf
#   - Appendix Figure A2 and related panels (race, language, household size,
#     e-mail, region, marital status, gender)
#     -> cate_race.pdf, cate_language.pdf, cate_hh_members.pdf, cate_email.pdf,
#        cate_regions_1.pdf, cate_regions_2.pdf, cate_marriage.pdf, cate_gender.pdf
#   - Shared horizontal legend -> cate_subgroup_legend_horizontal.pdf
#   - cate_table (all subgroup estimates, in memory; not written to disk here)
# Depends on: 04_heterogeneity_forest.R (forest_est, tau_hat)
#             02_prepare_data.R (X, df, race_cols/lang_cols, arm_cols, arm_lbl, lbl, cat_labels)
# =====================================================================

# =============================================================================
# CATE PLOTS BY SUBGROUP -- helper functions
# =============================================================================
# I use one generic pattern for every figure: a NAMED list of logical
# vectors (length n), so dummy sets, custom cutoffs, and integer categories
# all go through the same helper. Subgroup ATEs come from forest_est (the
# fine-tuned estimation forest); tau_hat supplies the "All" reference row.
#
# CAUTION: these are descriptive. Nothing below is corrected for
# multiplicity, and the region panels alone are 19 x 4 = 76 subgroup ATEs.

cate_by_groups <- function(groups, min_n = 300) {
  sizes <- vapply(groups, function(g) sum(g, na.rm = TRUE), integer(1))
  ok    <- sizes >= min_n
  if (any(!ok))
    message("Dropped (n < ", min_n, "): ",
            paste0(names(groups)[!ok], " (n=", sizes[!ok], ")", collapse = ", "))
  groups <- groups[ok]
  stopifnot(length(groups) > 0)

  bind_rows(lapply(names(groups), function(lab) {
    g <- groups[[lab]] & !is.na(groups[[lab]])
    a <- average_treatment_effect(forest_est, subset = g)
    a$contrast <- factor(active_arms, levels = active_arms)
    a$group    <- lab
    a$n_group  <- sum(g)
    a
  }))
}

# I compute a shared y-range across several subgroup panels (+ the "All"
# reference row), so every figure in the same scale-group uses identical
# axis limits.
shared_ylim <- function(group_list, min_n = 300, pad = 0.08) {
  d_all <- bind_rows(lapply(group_list, cate_by_groups, min_n = min_n))
  d_all <- bind_rows(d_all, tau_hat)
  lo <- min(d_all$estimate - 1.96 * d_all$std.err) * 100
  hi <- max(d_all$estimate + 1.96 * d_all$std.err) * 100
  rng <- hi - lo
  c(lo - pad * rng, hi + pad * rng)
}

plot_cate_groups <- function(groups, xlab, file, min_n = 300,
                             flip = FALSE, facet_arm = FALSE,
                             show_n = FALSE, width = 18, height = 15,
                             dodge_w = 0.5, ylim = NULL, y_breaks = NULL, show_ylab = TRUE) {
  d <- cate_by_groups(groups, min_n)
  lev <- unique(d$group)

  lab_map <- setNames(
    if (show_n) sprintf("%s\n(n=%s)", lev, format(d$n_group[match(lev, d$group)],
                                                  big.mark = "", trim = TRUE))
    else lev, lev)

  plot_df <- tau_hat %>%
    mutate(group = "All", n_group = nrow(X)) %>%
    bind_rows(d) %>%
    mutate(group = factor(group, levels = c("All", lev)))

  p <- ggplot(plot_df) +
    geom_pointrange(aes(x = group, y = estimate * 100,
                        ymin = (estimate - 1.96 * std.err) * 100,
                        ymax = (estimate + 1.96 * std.err) * 100,
                        col = contrast, linetype = contrast, shape = contrast),
                    linewidth = 1.2, size = 0.85,
                    position = position_dodge(width = dodge_w)) +
    geom_hline(yintercept = 0, alpha = 0.3) +
    scale_x_discrete(labels = function(x) ifelse(x == "All", "All", lab_map[x])) +
    theme_classic(base_size = 15) +
    theme(text = element_text(size = 14), legend.position = "none") +
    scale_color_manual("Arm", values = arm_cols, labels = arm_lbl) +
    scale_shape_discrete("Arm", labels = arm_lbl) +
    scale_linetype_discrete("Arm", labels = arm_lbl) +
    xlab(xlab) + ylab("Effect of treatment letter on enrollment rate (percentage points)")

  if (!is.null(y_breaks)) p <- p + scale_y_continuous(breaks = y_breaks)
  if (!show_ylab) p <- p + theme(axis.title.y = element_blank())
  if (flip) {
    p <- p + coord_flip(ylim = ylim, clip = "off")
  } else if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim, clip = "off")
  }
  if (facet_arm) p <- p + facet_wrap(~ contrast, ncol = 2, labeller = as_labeller(arm_labels))

  ggsave(paste0(OUT_DIR, file), plot = p, width = width, height = height, unit = "cm")
  invisible(list(plot = p, table = d))
}

dummy_groups <- function(cols, labels = cols) {
  stopifnot(all(cols %in% colnames(X)))
  setNames(lapply(cols, function(c) X[, c] == 1), labels)
}

# =============================================================================
# CATE PLOTS BY SUBGROUP -- group definitions
# =============================================================================
# --- Pull variables out once -------------------------------------------------
fpl     <- X[, "fpl_imp"]
aptc    <- X[, "indv_max_aptc_amt"]
age     <- X[, "age_mean"]
county_referral <- X[, "sws"] == 1
email   <- X[, "hoh_email_ok"] == 1
hhm     <- X[, "hh_members"]
miss    <- X[, "fpl_missing"] == 1
obs     <- !miss
gender  <- X[, "curr_gender_num1"]
married <- X[, "married_infer"] == 1
no_adult <- df$no_eligible_adult == 1

# --- Group definitions (built once, reused for both plotting and ylim) ------
race_groups <- dummy_groups(race_cols, race_short)
lang_groups <- dummy_groups(lang_cols, lang_short)

fpl_groups <- list(
  "<200%"    = obs & fpl <  200,
  "200-400%" = obs & fpl >= 200 & fpl < 400
)

age_groups <- list(
  "< 26"           = age < 26 & !no_adult,
  "26-40"          = age > 26 & age <= 40 & !no_adult,
  "40+"            = age > 40,
  "No adult in HH" = no_adult
)

hhm_groups     <- list("1" = hhm == 1, "2" = hhm == 2, "3+" = hhm >= 3)
email_groups   <- list("No" = !email, "Yes" = email)
sws_groups     <- list("Open Enrollment" = !county_referral, "County Referrals" = county_referral)
region1_groups <- dummy_groups(paste0("region", 1:10), paste0(1:10))
region2_groups <- dummy_groups(paste0("region", 11:19), paste0(11:19))
married_groups <- list("Married" = married, "Single" = !married)
gender_groups  <- list("Female" = gender == 2, "Male" = gender == 1)

# I define subsidy ($/month) bins directly on real dollar cutoffs (rather
# than data-driven quantiles) so the bins are interpretable in the write-up.
aptc_groups <- list(
  "<$100"    = aptc > 0   & aptc <= 100,
  "$101-400" = aptc > 100 & aptc <= 400,
  "$401-700" = aptc > 400 & aptc <= 700,
  ">$700"    = aptc > 700
)
# Cell-size sanity check (per arm; all should clear the ~300/arm floor).
sapply(aptc_groups, sum, na.rm = TRUE)
lapply(names(aptc_groups), function(lab) table(df$W[aptc_groups[[lab]]])) %>%
  setNames(names(aptc_groups))

# --- Shared y-axis ranges -----------------------------------------------
# Income (fpl) and subsidy (aptc) share a manually set -2 to 4 scale for the
# income panel; the subsidy panel gets its own data-driven range below.
ylim_income          <- c(-2, 4)
ylim_subsidy         <- shared_ylim(list(aptc_groups))
ylim_age_funnel      <- c(-4, 4)
ylim_race            <- shared_ylim(list(race_groups))
ylim_funnel_email    <- shared_ylim(list(sws_groups, email_groups))
ylim_regions         <- shared_ylim(list(region1_groups, region2_groups), min_n = 500)
ylim_hhsize_married  <- shared_ylim(list(hhm_groups, married_groups, lang_groups))

# =============================================================================
# CATE PLOTS BY SUBGROUP -- figures
# =============================================================================
res_race <- plot_cate_groups(race_groups, xlab = cat_labels[["race"]],
  file = "cate_race.pdf", width = 16, height = 15, ylim = ylim_race, y_breaks = seq(-4, 6, by = 1))

res_lang <- plot_cate_groups(lang_groups, xlab = cat_labels[["language"]],
  file = "cate_language.pdf", width = 16, height = 15, ylim = c(-2, 3), y_breaks = seq(-2, 3, by = 1))

res_fpl <- plot_cate_groups(fpl_groups, xlab = lbl("fpl_imp"), file = "cate_fpl.pdf",
  width = 16, height = 15, ylim = ylim_income, y_breaks = seq(-3, 4, by = 1))

res_aptc <- plot_cate_groups(aptc_groups,
  xlab = paste0(lbl("indv_max_aptc_amt"), " ($/month)"), file = "cate_aptc.pdf",
  width = 16, height = 15, ylim = ylim_subsidy, show_ylab = TRUE, y_breaks = seq(-5, 6, by = 1))

res_age <- plot_cate_groups(age_groups,
  xlab = paste0(lbl("age_mean"), " (years)"), file = "cate_age.pdf",
  width = 16, height = 15, ylim = ylim_age_funnel, show_ylab = TRUE, y_breaks = seq(-4, 4, by = 1))

res_hhm <- plot_cate_groups(hhm_groups, xlab = lbl("hh_members"),
  file = "cate_hh_members.pdf", width = 16, height = 15, ylim = ylim_hhsize_married)

res_email <- plot_cate_groups(email_groups, xlab = lbl("hoh_email_ok"),
  file = "cate_email.pdf", width = 16, height = 15, ylim = ylim_funnel_email,
  y_breaks = seq(-1, 3, by = 1), show_ylab = FALSE)

res_sws <- plot_cate_groups(sws_groups, xlab = cat_labels[["funnel"]],
  file = "cate_funnel_sws.pdf", width = 16, height = 15, ylim = c(-1, 3),
  y_breaks = seq(-1, 3, by = 1), show_ylab = TRUE)

res_region_1 <- plot_cate_groups(region1_groups,
  xlab = paste0(cat_labels[["region"]], " (1-10)"), file = "cate_regions_1.pdf",
  min_n = 500, width = 34, height = 16, dodge_w = 0.7, ylim = ylim_regions)

res_region_2 <- plot_cate_groups(region2_groups,
  xlab = paste0(cat_labels[["region"]], " (11-19)"), file = "cate_regions_2.pdf",
  min_n = 500, width = 34, height = 16, dodge_w = 0.7, ylim = ylim_regions, show_ylab = TRUE)

res_married <- plot_cate_groups(married_groups, xlab = lbl("married_infer"),
  file = "cate_marriage.pdf", width = 14, height = 13, ylim = ylim_hhsize_married)

res_gender <- plot_cate_groups(gender_groups, xlab = lbl("curr_gender_num1"),
  file = "cate_gender.pdf", width = 16, height = 15, ylim = ylim_funnel_email)

# --- Single, shared, horizontal legend for all subgroup CATE figures --------
legend_source <- res_race$plot +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)),
         shape = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1))

subgroup_legend <- tryCatch(
  cowplot::get_plot_component(legend_source, "guide-box-bottom", return_all = TRUE),
  error = function(e) cowplot::get_legend(legend_source)
)
stopifnot(!is.null(subgroup_legend))

ggsave(paste0(OUT_DIR, "cate_subgroup_legend_horizontal.pdf"),
       plot = cowplot::ggdraw(subgroup_legend),
       width = 18, height = 1.6, unit = "cm", dpi = 300, bg = "white")

# --- Export every subgroup estimate ------------------------------------------
# NB: I name the region entries by the PANEL (1-10 / 11-19), not by the
# covariates region1 / region2, which would collide in cate_table$variable.
res_all <- list(
  race = res_race, language = res_lang, fpl = res_fpl, aptc = res_aptc, age = res_age,
  hh_members = res_hhm, hoh_email_ok = res_email, sws = res_sws,
  regions_1_10 = res_region_1, regions_11_19 = res_region_2,
  married = res_married, gender = res_gender
)

cate_table <- bind_rows(lapply(names(res_all), function(v) {
  res_all[[v]]$table %>% mutate(variable = v, .before = 1)
})) %>%
  mutate(arm_label = arm_lbl(contrast),
         ci_lo  = estimate - 1.96 * std.err,
         ci_hi  = estimate + 1.96 * std.err,
         est_pp = estimate * 100,
         se_pp  = std.err  * 100)
