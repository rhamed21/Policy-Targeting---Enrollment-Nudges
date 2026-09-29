# =====================================================================
# 03_balance_table.R
# Checks that randomization succeeded: by-arm means, overall mean, and a
# heteroskedasticity-robust joint F-test that the arm means are equal.
#
# Produces: Table 1 (Summary Statistics, by arm)
#   -> balance_table.csv, balance_table.tex
# Depends on: 02_prepare_data.R (df, arms, arm_lbl, lbl)
# =====================================================================

df$fpl_b400 <- ifelse(df$subsidy_fpl_percent < 400, df$subsidy_fpl_percent, NA)

balance_vars <- c(
  "age_mean", "curr_gender_num1", "married_infer", "hh_members", "hh_kids", "hoh_email_ok",
  "fpl_b400", "indv_max_aptc_amt", "fpl_missing",
  "latino", "black", "white", "asian", "otherrace",
  "English", "Spanish", "Other"
)

balance_labels <- case_when(
  balance_vars == "fpl_b400"          ~ "Income (% of FPL)",
  balance_vars == "indv_max_aptc_amt" ~ "Subsidy ($/month)",
  balance_vars == "curr_gender_num1"  ~ "Gender",
  balance_vars == "hoh_email_ok"      ~ "E-mail Available",
  balance_vars == "hh_members"        ~ "Household Size",
  balance_vars == "hh_kids"           ~ "Number of Kids",
  balance_vars == "fpl_missing"       ~ "Missing Income",
  balance_vars == "otherrace"         ~ "Other Race",
  balance_vars == "Other"             ~ "Other Language",
  TRUE                                ~ lbl(balance_vars)
)

balance_rows <- lapply(seq_along(balance_vars), function(i) {
  v <- balance_vars[i]
  x <- df[[v]]
  if (!is.numeric(x)) x <- as.numeric(x)

  means      <- tapply(x, df$W, mean, na.rm = TRUE)
  total_mean <- mean(x, na.rm = TRUE)

  keep <- !is.na(x)
  m_full <- lm(x[keep] ~ factor(df$W[keep]))
  m_null <- lm(x[keep] ~ 1)
  p_joint <- waldtest(m_null, m_full, vcov = vcovHC(m_full, type = "HC1"))$`Pr(>F)`[2]

  tibble(
    covariate = balance_labels[i],
    !!!setNames(as.list(sprintf("%.2f", means[arms])), arm_lbl(arms)),
    Total   = sprintf("%.2f", total_mean),
    p_joint = sprintf("%.2f", p_joint)
  )
})
balance_table <- bind_rows(balance_rows)
write.csv(balance_table, paste0(OUT_DIR, "balance_table.csv"), row.names = FALSE)

n_row <- tibble(
  covariate = "Observations",
  !!!setNames(as.list(format(table(df$W)[arms], big.mark = ",")), arm_lbl(arms)),
  Total   = format(sum(table(df$W)[arms]), big.mark = ","),
  p_joint = ""
)
display_table <- bind_rows(balance_table, n_row)
n_row_idx <- nrow(display_table)
n_data    <- ncol(display_table) - 1

notes_text <- paste(
  "Table 1 reports income and demographic characteristics reported in the administrative data for the final study",
  "sample, by treatment arm. Column 7 reports the $p$-value on the null that the mean across all arms are equal.",
  "$p$-values are computed from heteroskedasticity-robust (HC1) standard errors."
)

col_headers_tex <- c(
  "",
  linebreak(paste0(arm_lbl(arms), "\n(", seq_along(arms), ")"), align = "c"),
  linebreak(paste0("Total\n(", length(arms) + 1, ")"), align = "c"),
  linebreak(paste0("p-value\nall arms\n(", length(arms) + 2, ")"), align = "c")
)

kbl_out_tex <- display_table %>%
  kbl(align = c("l", rep("c", n_data)), booktabs = TRUE, format = "latex", linesep = "",
      col.names = col_headers_tex, escape = FALSE,
      caption = "Summary Statistics, by arm") %>%
  kable_styling(font_size = NULL, latex_options = "hold_position", position = "center") %>%
  add_header_above(c(" " = 1, "Treatment Arm" = length(arms), " " = 2)) %>%
  row_spec(n_row_idx - 1, extra_css = "border-bottom: 1px solid black;")

writeLines(as.character(kbl_out_tex), paste0(OUT_DIR, "balance_table.tex"))
