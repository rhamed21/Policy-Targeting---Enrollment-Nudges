# =====================================================================
# 02_prepare_data.R
# Loads the public-use analysis file, recodes the outcome/treatment,
# builds the covariate matrix, and defines display names/colours reused
# by every figure in the later scripts.
#
# Produces: nothing (no figures/tables)
# Depends on: 01_setup.R (packages, DATA_PATH from master.R)
# =====================================================================

# ---- Treatment arms -----------------------------------------------------
# I define the treatment arms in the order the RCT used them:
#   arm 1: control
#   arm 2: basic reminder
#   arm 3: subsidy + penalty letter (and arm 2)
#   arm 4: price comparison letter (and arms 2-3)
#   arm 5: price + quality comparison letter (and arms 2-4)
arms        <- c("Control", "Basic Reminder", "SubsidyPenalty", "PriceCompare", "PriceQuality")
n_arms      <- length(arms)
active_arms <- setdiff(arms, "Control")

# I build one shared colour palette for the active arms here, so every plot
# later in the document (CATE panels, density plot, variable importance,
# trees) references the same `arm_cols` instead of re-deriving it.
arm_cols <- setNames(scales::hue_pal()(length(active_arms)), active_arms)

# `arms` stays the internal key; arm_labels are the display names.
arm_labels <- c(
  Control          = "Control",
  "Basic Reminder" = "Basic Letter",
  SubsidyPenalty   = "Subsidy + Penalty",
  PriceCompare     = "Price Compare",
  PriceQuality     = "Price-Quality Compare"
)
stopifnot(setequal(names(arm_labels), arms))
arm_lbl <- function(a) unname(arm_labels[a])

# I load the paper's post-exclusion analysis sample (n = 87,394) and recode
# the outcome, treatment index, and treatment factor.
raw <- read_dta(DATA_PATH)
df <- raw %>%
  filter(samplea == 1) %>%
  mutate(
    Y  = as.numeric(any_takeup),
    Wn = as.integer(rctarm),
    W  = factor(rctarm, levels = 1:5, labels = arms)
  )
n <- nrow(df)

# =====================================================================
# Covariates
# =====================================================================

# --- Income --------------------------------------------------------------
# I keep fpl_missing / fpl_1000p as explicit indicator covariates, while
# fpl_imp itself stays RAW subsidy_fpl_percent, unmodified:
#   - Missing income stays genuine NA -> NaN, so grf's causal forests can
#     split on it via "missing incorporated in attributes" (MIA) rather
#     than my guessing an imputed value.
#   - >1000% FPL values are left at their raw (large) magnitude; fpl_1000p
#     flags them for the forest/tree to use, but fpl_imp is untouched.
#     (min observed FPL is ~100, so 0 sits safely below the support.)
df <- df %>%
  mutate(
    fpl_missing = as.numeric(fpl_missing == 1),
    fpl_1000p   = as.numeric(subsidy_fpl_percent > 1000 & !is.na(subsidy_fpl_percent)),
    fpl_imp     = subsidy_fpl_percent
  )

# --- Written-language dummies --------------------------------------------
# I one-hot encode the language variable myself instead of trusting
# model.matrix() blindly, because model.matrix() silently DROPS rows with
# NA (which would misalign the cbind below) and the labels assume the
# numeric codes are exactly 1/2/3 sorting as English/Spanish/Other. I check
# both assumptions rather than assume them.
stopifnot(!anyNA(df$prefer_written_lang_num1),
          setequal(sort(unique(df$prefer_written_lang_num1)), 1:3))

lang_dummies <- model.matrix(
  ~ factor(prefer_written_lang_num1, labels = c("English", "Spanish", "Other")) - 1,
  data = df
)
stopifnot(nrow(lang_dummies) == nrow(df))
colnames(lang_dummies) <- c("English", "Spanish", "Other")
df <- cbind(df, lang_dummies)

# --- Covariate list --------------------------------------------------------
all_covs <- c(
  "age_mean", "curr_gender_num1",
  "fpl_missing", "fpl_1000p", "fpl_imp",
  "hh_members", "hh_kids", "hoh_email_ok", "indv_max_aptc_amt", "married_infer",
  "latino", "black", "white", "asian", "otherrace",
  "English", "Spanish", "Other",
  paste0("region", 1:19),
  "sws"
)

var_labels <- c(
  age_mean = "Household Age (Mean)", curr_gender_num1 = "Gender of Head of Household",
  sws = "Funnel", fpl_missing = "Missing Income", fpl_1000p = "Income Over 1000% FPL",
  fpl_imp = "Income (% of FPL)", hh_members = "Household Size", hh_kids = "Number of Kids",
  hoh_email_ok = "E-mail Available", indv_max_aptc_amt = "Subsidy", married_infer = "Married",
  x7 = "Subsidy-Ineligible Aid Code",
  white = "White", black = "Black", asian = "Asian", latino = "Latino", otherrace = "Other Race",
  English = "English", Spanish = "Spanish", Other = "Other Language",
  # Binned analogues used only by the policy tree
  age_bin = "Household Age (Mean)", fpl_bin = "Income (% of FPL)", aptc_bin = "Subsidy"
)
var_labels[paste0("region", 1:19)] <- paste("Region", 1:19)

# I keep short labels for the standalone race / language figures, where the
# category is already named on the axis so "Other" alone is unambiguous.
race_cols  <- c("white", "black", "asian", "latino", "otherrace")
race_short <- c("White", "Black", "Asian", "Latino", "Other")
lang_cols  <- c("English", "Spanish", "Other")
lang_short <- c("English", "Spanish", "Other")

cat_labels <- c(race = "Race", language = "Language", funnel = "Funnel",
                region = "California Rating Region")

# I unname() so downstream paste()/sprintf() calls never inherit stray names.
lbl <- function(v) unname(ifelse(v %in% names(var_labels), var_labels[v], v))

# --- Estimation inputs -----------------------------------------------------
X <- df %>% select(all_of(all_covs)) %>% as.matrix()
Y <- df$Y
W <- df$W

# Propensities are known by RCT design (equal allocation across 5 arms).
W.hat <- matrix(0.2, nrow = nrow(X), ncol = n_arms, dimnames = list(NULL, arms))
