# =============================================================================
# Jenkins & Van Kerm (2017) Attrition Diagnostics
#
# This script implements a three-step attrition analysis following the
# Jenkins & Van Kerm (2017) methodology for poverty persistence studies:
#
#   Step 1 — Descriptive attrition patterns:
#     Attrition flow diagram (wave-to-wave and cumulative), baseline
#     characteristic comparison (chi-square / t-test), and attrition
#     rates stratified by income decile, poverty status, employment,
#     formality, and household structure.
#
#   Step 2 — Cross-sectional vs. longitudinal representativeness test:
#     For each wave, compare the poverty rate from the full cross-section
#     with the poverty rate from the balanced longitudinal panel. The null
#     hypothesis (ΔP = 0) is tested with a two-proportion z-test.
#     An extended version repeats the comparison for persistent poverty.
#
#   Step 3 — Subgroup representativeness assessment:
#     Repeat Step 2 stratified by employment status, formality, baseline
#     poverty status, household size, head education, and head age group.
#
# Outputs:  tables/attrition_flow.csv
#           tables/attrition_baseline_comparison.csv
#           tables/attrition_by_subgroup.csv
#           tables/jvk_representativeness_test.csv
#           tables/jvk_persistent_poverty_test.csv
#           tables/jvk_subgroup_representativeness.csv
#           (+ gt HTML/LaTeX versions when gt is available)
#
# Reference:
#   Jenkins, S. P. & Van Kerm, P. (2017). How does attrition affect
#   estimates of persistent poverty rates? The case of EU-SILC.
#   In: Ferriss, A. L. (ed.) Advances in Quality-of-Life Theory and
#   Research, pp. 401–461. Springer.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(haven)
  library(purrr)
  library(readr)
  library(survey)
  library(tibble)
  library(tidyr)
})

# ---- Locate project root and source configuration / functions ---------------

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  sourced_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )
  if (!is.na(sourced_file)) return(sourced_file)
  NA_character_
}

script_path <- get_script_path()
project_root <- if (is.na(script_path)) {
  getwd()
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

source(file.path(project_root, "00_config.R"))
source(file.path(project_root, "01_functions.R"))

project$data_path  <- file.path(project_root, project$data_path)
project$out_dir    <- file.path(project_root, project$out_dir)
project$table_dir  <- file.path(project_root, project$table_dir)
project$figure_dir <- file.path(project_root, project$figure_dir)
project$model_dir  <- file.path(project_root, project$model_dir)

dir.create(project$out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(project$table_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read data and build analysis panel -------------------------------------

cat("Reading raw SILC panel data ...\n")
raw <- read_dta(project$data_path)

cat("Constructing balanced panel ...\n")
panel_balanced <- construct_balanced_panel(
  raw = raw, vars = vars,
  panel_years = project$panel_years,
  reference_year = project$reference_year
)

cat("Computing equivalised income and household context ...\n")
panel_income <- panel_balanced %>%
  add_equivalised_income(vars) %>%
  add_household_context(vars, codes) %>%
  add_employment_stability(vars)

cat("Adding poverty status ...\n")
poverty_lines <- compute_poverty_lines(panel_income, vars, project$thresholds)
panel_poverty <- add_poverty_status(panel_income, poverty_lines, vars, project$thresholds)

cat("Classifying poverty spells ...\n")
classified_main <- classify_poverty_spells(
  panel       = panel_poverty,
  vars        = vars,
  panel_years = project$panel_years,
  reference_year = project$reference_year,
  threshold   = project$main_threshold
)

# ---- Shorthand aliases -------------------------------------------------------

id          <- vars$person_id
hh          <- vars$household_id
year_var    <- vars$year
age_var     <- vars$age
ref_var     <- vars$household_reference_person
income_var  <- vars$household_income
weight_long <- vars$longitudinal_weight
sex_var     <- vars$sex
edu_var     <- vars$education
labour_var  <- vars$labour_status
ss_var      <- vars$social_security

threshold   <- project$main_threshold / 100
panel_years <- project$panel_years
first_year  <- min(panel_years)
ref_year    <- project$reference_year

balanced_ids <- panel_balanced %>%
  distinct(.data[[id]]) %>%
  pull(.data[[id]])

# =============================================================================
# STEP 1: Descriptive Comparison of Attrition Patterns
# =============================================================================
cat("\n=== STEP 1: Descriptive attrition patterns ===\n")

# ---- 1a. Attrition flow diagram (wave-to-wave and cumulative) ---------------

# Identify who is present in each wave
wave_presence <- raw %>%
  filter(.data[[year_var]] %in% panel_years) %>%
  distinct(.data[[id]], .data[[year_var]])

# Build cumulative presence matrix
ids_by_wave <- lapply(panel_years, function(y) {
  wave_presence %>%
    filter(.data[[year_var]] == y) %>%
    pull(.data[[id]])
})
names(ids_by_wave) <- panel_years

# Attrition flow: cumulative survivors across consecutive waves
flow_rows <- list()
cumulative_ids <- ids_by_wave[[1]]
n_wave1 <- length(cumulative_ids)

for (i in seq_along(panel_years)) {
  y <- panel_years[i]
  current_ids <- ids_by_wave[[as.character(y)]]

  if (i == 1) {
    flow_rows[[i]] <- tibble(
      wave              = y,
      n_in_wave         = length(current_ids),
      stayed_from_prev  = NA_integer_,
      dropped_from_prev = NA_integer_,
      wave_to_wave_retention = NA_real_,
      wave_to_wave_attrition = NA_real_,
      cumulative_retention_from_wave1 = 1.0,
      cumulative_attrition_from_wave1 = 0.0
    )
  } else {
    stayed  <- length(intersect(cumulative_ids, current_ids))
    dropped <- length(cumulative_ids) - stayed
    cumulative_ids <- intersect(cumulative_ids, current_ids)

    flow_rows[[i]] <- tibble(
      wave              = y,
      n_in_wave         = length(current_ids),
      stayed_from_prev  = stayed,
      dropped_from_prev = dropped,
      wave_to_wave_retention = stayed / (stayed + dropped),
      wave_to_wave_attrition = dropped / (stayed + dropped),
      cumulative_retention_from_wave1 = length(cumulative_ids) / n_wave1,
      cumulative_attrition_from_wave1 = 1 - length(cumulative_ids) / n_wave1
    )
  }
}

attrition_flow <- bind_rows(flow_rows)

cat("\nAttrition flow diagram:\n")
print(attrition_flow)

# ---- 1b. Baseline comparison: stayers vs. attritors -------------------------

# Compute equivalised income at baseline from the raw data
baseline_raw <- raw %>%
  filter(.data[[year_var]] == first_year) %>%
  mutate(
    oecd_weight_bl = case_when(
      .data[[ref_var]] == 1                           ~ 1.0,
      .data[[ref_var]] != 1 & .data[[age_var]] >= 14  ~ 0.5,
      .data[[ref_var]] != 1 & .data[[age_var]] < 14   ~ 0.3,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(.data[[hh]]) %>%
  mutate(
    hh_eq_size_bl = sum(oecd_weight_bl, na.rm = TRUE),
    eq_income_bl  = .data[[income_var]] / hh_eq_size_bl,
    hh_size_bl    = n_distinct(.data[[id]])
  ) %>%
  ungroup() %>%
  filter(!is.na(eq_income_bl), hh_eq_size_bl > 0)

# Baseline poverty line (60% of weighted median at wave 1, equal weights)
baseline_poverty_line <- weighted_median(baseline_raw$eq_income_bl, rep(1, nrow(baseline_raw))) * threshold

baseline_attrition <- baseline_raw %>%
  mutate(
    stayer  = .data[[id]] %in% balanced_ids,
    attritor = !stayer,
    stayer_label = if_else(stayer, "Stayer", "Attritor"),
    poor_baseline = as.integer(eq_income_bl < baseline_poverty_line),
    female = if_else(
      .data[[sex_var]] == codes$sex[["female"]], 1L, 0L, missing = NA_integer_),
    employed = if_else(
      .data[[labour_var]] %in% codes$employed_labour_status, 1L, 0L, missing = NA_integer_),
    informal = if_else(
      .data[[ss_var]] %in% codes$likely_informal_social_security_values,
      1L, 0L, missing = 0L),
    edu_primary   = if_else(.data[[edu_var]] %in% c(0, 1, 2), 1L, 0L, missing = NA_integer_),
    edu_secondary = if_else(.data[[edu_var]] %in% c(3, 4, 5), 1L, 0L, missing = NA_integer_),
    edu_tertiary  = if_else(.data[[edu_var]] %in% c(6, 7, 8), 1L, 0L, missing = NA_integer_),
    income_decile = ntile(eq_income_bl, 10),
    age_group_bl  = cut(
      .data[[age_var]],
      breaks = c(-Inf, 17, 34, 54, 64, Inf),
      labels = c("0-17", "18-34", "35-54", "55-64", "65+")
    )
  )

# ---- 1b-i. Descriptive table: mean characteristics by stayer/attritor -------

fmt_v <- function(x) sprintf("%.3f", x)
fmt_p <- function(p) {
  p <- as.numeric(p)
  stars <- add_significance_stars(p)
  paste0(ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)), stars)
}

# Survey-weighted t-test p-value
wt_ttest_p <- function(var, data = baseline_attrition) {
  df <- data %>% filter(!is.na(.data[[var]]))
  des <- svydesign(ids = ~1, weights = ~1, data = df)
  tt <- tryCatch(
    svyttest(as.formula(paste0(var, " ~ stayer")), design = des),
    error = function(e) NULL
  )
  if (is.null(tt)) return(NA_real_)
  tt$p.value
}

# Chi-square p-value for categorical variable
chisq_p <- function(var, data = baseline_attrition) {
  df <- data %>% filter(!is.na(.data[[var]]))
  tt <- tryCatch(
    chisq.test(table(df[[var]], df$stayer_label)),
    error = function(e) NULL
  )
  if (is.null(tt)) return(NA_real_)
  tt$p.value
}

# Weighted mean by stayer/attritor group
wmean_by_group <- function(var, data = baseline_attrition) {
  data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(stayer_label) %>%
    summarise(m = mean(.data[[var]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = stayer_label, values_from = m)
}

build_comparison_row <- function(label, var, test_type = "t") {
  vals <- wmean_by_group(var)
  p <- if (test_type == "t") wt_ttest_p(var) else chisq_p(var)
  tibble(
    Characteristic = label,
    Stayer   = fmt_v(vals$Stayer),
    Attritor = fmt_v(vals$Attritor),
    Difference = fmt_v(vals$Stayer - vals$Attritor),
    p_value  = fmt_p(p)
  )
}

n_stayers   <- sum(baseline_attrition$stayer)
n_attritors <- sum(baseline_attrition$attritor)

baseline_comparison <- bind_rows(
  tibble(Characteristic = "N (unweighted)",
         Stayer = as.character(n_stayers),
         Attritor = as.character(n_attritors),
         Difference = "", p_value = ""),
  build_comparison_row("Age, mean",                             age_var,    "t"),
  build_comparison_row("Female, proportion",                    "female",   "t"),
  build_comparison_row("Household income, mean",                income_var, "t"),
  build_comparison_row("Equivalised income, mean",              "eq_income_bl", "t"),
  build_comparison_row("Primary or below, proportion",          "edu_primary",   "t"),
  build_comparison_row("Secondary, proportion",                 "edu_secondary", "t"),
  build_comparison_row("Tertiary, proportion",                  "edu_tertiary",  "t"),
  build_comparison_row("Employed, proportion",                  "employed",      "t"),
  build_comparison_row("Informal (not registered), proportion", "informal",      "t"),
  build_comparison_row("Poor at baseline (60%), proportion",    "poor_baseline", "t")
)

cat("\nBaseline comparison (stayers vs. attritors):\n")
print(baseline_comparison)

# ---- 1c. Attrition rates by subgroups --------------------------------------

attrition_rate_by <- function(data, group_var, group_label) {
  data %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(group = group_label, category = as.character(.data[[group_var]])) %>%
    summarise(
      n = n(),
      n_attritors = sum(attritor),
      attrition_rate = mean(attritor),
      .groups = "drop"
    )
}

attrition_by_subgroup <- bind_rows(
  # By income decile
  attrition_rate_by(baseline_attrition, "income_decile", "Income decile"),
  # By baseline poverty status
  baseline_attrition %>%
    mutate(poverty_label = if_else(poor_baseline == 1L, "Poor", "Non-poor")) %>%
    attrition_rate_by("poverty_label", "Baseline poverty status"),
  # By employment status
  baseline_attrition %>%
    mutate(emp_label = if_else(employed == 1L, "Employed", "Not employed")) %>%
    attrition_rate_by("emp_label", "Employment status"),
  # By formality
  baseline_attrition %>%
    mutate(form_label = if_else(informal == 1L, "Informal", "Formal")) %>%
    attrition_rate_by("form_label", "Formality"),
  # By sex
  baseline_attrition %>%
    mutate(sex_label = if_else(female == 1L, "Female", "Male")) %>%
    attrition_rate_by("sex_label", "Sex"),
  # By age group
  attrition_rate_by(baseline_attrition, "age_group_bl", "Age group"),
  # By household size bracket
  baseline_attrition %>%
    mutate(hh_size_bracket = case_when(
      hh_size_bl == 1 ~ "1",
      hh_size_bl <= 3 ~ "2-3",
      hh_size_bl <= 5 ~ "4-5",
      TRUE            ~ "6+"
    )) %>%
    attrition_rate_by("hh_size_bracket", "Household size"),
  # By education level
  baseline_attrition %>%
    mutate(edu_label = case_when(
      edu_primary   == 1L ~ "Primary or below",
      edu_secondary == 1L ~ "Secondary",
      edu_tertiary  == 1L ~ "Tertiary",
      TRUE ~ NA_character_
    )) %>%
    attrition_rate_by("edu_label", "Education")
)

cat("\nAttrition rates by subgroup:\n")
print(attrition_by_subgroup, n = 40)

# ---- 1d. Chi-square and t-test significance for key stratifications ---------

chi_poverty   <- chisq_p("poor_baseline")
chi_employed  <- chisq_p("employed")
chi_informal  <- chisq_p("informal")
chi_female    <- chisq_p("female")
chi_education <- tryCatch({
  df <- baseline_attrition %>%
    mutate(edu_label = case_when(
      edu_primary   == 1L ~ "Primary",
      edu_secondary == 1L ~ "Secondary",
      edu_tertiary  == 1L ~ "Tertiary"
    )) %>%
    filter(!is.na(edu_label))
  chisq.test(table(df$edu_label, df$stayer_label))$p.value
}, error = function(e) NA_real_)

attrition_significance <- tibble(
  variable = c("Baseline poverty", "Employment", "Informal", "Female", "Education"),
  test = "Chi-square",
  p_value = c(chi_poverty, chi_employed, chi_informal, chi_female, chi_education),
  significance = add_significance_stars(p_value)
)

cat("\nChi-square tests for differential attrition:\n")
print(attrition_significance)

# =============================================================================
# STEP 2: Cross-sectional vs. Longitudinal Representativeness Test
#         (Core Jenkins & Van Kerm innovation)
# =============================================================================
cat("\n=== STEP 2: Cross-sectional vs. longitudinal representativeness test ===\n")

# Helper: compute OECD equivalised income for any subset of the raw data
compute_eq_income_raw <- function(data) {
  data %>%
    mutate(
      oecd_w = case_when(
        .data[[ref_var]] == 1                          ~ 1.0,
        .data[[ref_var]] != 1 & .data[[age_var]] >= 14 ~ 0.5,
        .data[[ref_var]] != 1 & .data[[age_var]] < 14  ~ 0.3,
        TRUE ~ NA_real_
      )
    ) %>%
    group_by(.data[[hh]], .data[[year_var]]) %>%
    mutate(
      hh_eq_size = sum(oecd_w, na.rm = TRUE),
      eq_income_raw = .data[[income_var]] / hh_eq_size
    ) %>%
    ungroup() %>%
    filter(!is.na(eq_income_raw), hh_eq_size > 0)
}

# Two-proportion z-test for difference in poverty rates
prop_z_test <- function(n_poor_1, n_total_1, n_poor_2, n_total_2) {
  p1 <- n_poor_1 / n_total_1
  p2 <- n_poor_2 / n_total_2
  p_pool <- (n_poor_1 + n_poor_2) / (n_total_1 + n_total_2)
  se <- sqrt(p_pool * (1 - p_pool) * (1 / n_total_1 + 1 / n_total_2))
  if (se == 0) return(list(z = NA_real_, p = NA_real_))
  z <- (p1 - p2) / se
  p <- 2 * pnorm(-abs(z))
  list(z = z, p = p)
}

# ---- 2a. Poverty rate comparison by wave ------------------------------------

jvk_results <- purrr::map_dfr(panel_years, function(y) {
  # Cross-sectional sample: all individuals present in wave y
  xs_wave <- raw %>%
    filter(.data[[year_var]] == y) %>%
    compute_eq_income_raw()

  xs_poverty_line <- weighted_median(xs_wave$eq_income_raw, rep(1, nrow(xs_wave))) * threshold
  xs_n_poor  <- sum(xs_wave$eq_income_raw < xs_poverty_line)
  xs_n_total <- nrow(xs_wave)
  xs_rate    <- xs_n_poor / xs_n_total

  # Longitudinal sample: balanced panel members present in wave y
  long_wave <- panel_income %>%
    filter(.data[[year_var]] == y, .data[[id]] %in% balanced_ids)

  long_poverty_line <- weighted_median(
    long_wave$eq_income, long_wave[[weight_long]]
  ) * threshold
  long_n_poor  <- sum(long_wave$eq_income < long_poverty_line)
  long_n_total <- nrow(long_wave)
  long_rate    <- long_n_poor / long_n_total

  # ΔP and z-test
  delta_p <- xs_rate - long_rate
  z_test  <- prop_z_test(xs_n_poor, xs_n_total, long_n_poor, long_n_total)

  # Confidence intervals (normal approximation)
  xs_se   <- sqrt(xs_rate * (1 - xs_rate) / xs_n_total)
  long_se <- sqrt(long_rate * (1 - long_rate) / long_n_total)

  tibble(
    wave = y,
    cross_sectional_n     = xs_n_total,
    cross_sectional_rate  = xs_rate,
    cross_sectional_ci_lo = xs_rate - 1.96 * xs_se,
    cross_sectional_ci_hi = xs_rate + 1.96 * xs_se,
    longitudinal_n        = long_n_total,
    longitudinal_rate     = long_rate,
    longitudinal_ci_lo    = long_rate - 1.96 * long_se,
    longitudinal_ci_hi    = long_rate + 1.96 * long_se,
    delta_p               = delta_p,
    z_statistic           = z_test$z,
    p_value               = z_test$p,
    significance          = add_significance_stars(z_test$p),
    interpretation = case_when(
      abs(delta_p) < 0.005               ~ "No meaningful difference",
      delta_p > 0 & z_test$p < 0.05      ~ "Poor more likely to stay (negative selection into panel)",
      delta_p < 0 & z_test$p < 0.05      ~ "Poor more likely to leave (positive selection out of panel)",
      TRUE                                ~ "Difference not statistically significant"
    )
  )
})

cat("\nJ&VK representativeness test (poverty rates by wave):\n")
print(jvk_results %>% select(wave, cross_sectional_rate, longitudinal_rate,
                              delta_p, z_statistic, p_value, significance, interpretation))

# ---- 2b. Extended test for persistent poverty (Core J&VK) -------------------

# Persistent poverty = poor in current year AND poor in >= 2 of preceding 3 years
# For year 4 (2019): poor in 2019 AND poor in >= 2 of {2016, 2017, 2018}

# Cross-sectional persistent poverty (all respondents present in year 4)
# We can only compute persistent poverty for individuals observed in all 4 waves
# (need prior-year status), so the cross-sectional version uses everyone
# present in year 4 for whom we can compute at least 3 preceding observations.

xs_year4 <- raw %>%
  filter(.data[[year_var]] %in% panel_years) %>%
  compute_eq_income_raw()

# Build per-year poverty lines from the cross-sectional sample
xs_poverty_lines <- xs_year4 %>%
  group_by(.data[[year_var]]) %>%
  summarise(
    poverty_line_xs = weighted_median(eq_income_raw, rep(1, n())) * threshold,
    .groups = "drop"
  )

xs_year4 <- xs_year4 %>%
  left_join(xs_poverty_lines, by = year_var) %>%
  mutate(poor_xs = as.integer(eq_income_raw < poverty_line_xs))

# Identify individuals present in year 4
ids_year4 <- xs_year4 %>%
  filter(.data[[year_var]] == ref_year) %>%
  pull(.data[[id]])

# Build poverty history for everyone present in year 4
xs_poverty_wide <- xs_year4 %>%
  filter(.data[[id]] %in% ids_year4) %>%
  select(all_of(c(id, year_var)), poor_xs) %>%
  distinct(.data[[id]], .data[[year_var]], .keep_all = TRUE) %>%
  pivot_wider(
    names_from  = all_of(year_var),
    values_from = poor_xs,
    names_prefix = "poor_"
  )

previous_years <- setdiff(panel_years, ref_year)
prev_poor_vars <- paste0("poor_", previous_years)
current_var    <- paste0("poor_", ref_year)

xs_persistent <- xs_poverty_wide %>%
  rowwise() %>%
  mutate(
    n_prev_years_observed = sum(!is.na(c_across(any_of(prev_poor_vars)))),
    n_prev_poor = sum(c_across(any_of(prev_poor_vars)), na.rm = TRUE),
    current_poor = .data[[current_var]],
    persistent_poor_xs = !is.na(current_poor) & current_poor == 1 &
                         n_prev_years_observed >= 2 & n_prev_poor >= 2
  ) %>%
  ungroup() %>%
  filter(!is.na(current_poor), n_prev_years_observed >= 2)

xs_persist_n_poor  <- sum(xs_persistent$persistent_poor_xs)
xs_persist_n_total <- nrow(xs_persistent)
xs_persist_rate    <- xs_persist_n_poor / xs_persist_n_total

# Longitudinal persistent poverty (balanced panel only)
long_persistent <- classified_main %>%
  filter(!is.na(persistent_poor))
long_persist_n_poor  <- sum(long_persistent$persistent_poor)
long_persist_n_total <- nrow(long_persistent)
long_persist_rate    <- long_persist_n_poor / long_persist_n_total

# Test
delta_persist <- xs_persist_rate - long_persist_rate
z_persist     <- prop_z_test(xs_persist_n_poor, xs_persist_n_total,
                              long_persist_n_poor, long_persist_n_total)

jvk_persistent <- tibble(
  metric = "Persistent poverty (year 4, anchored)",
  cross_sectional_n     = xs_persist_n_total,
  cross_sectional_rate  = xs_persist_rate,
  longitudinal_n        = long_persist_n_total,
  longitudinal_rate     = long_persist_rate,
  delta_p               = delta_persist,
  z_statistic           = z_persist$z,
  p_value               = z_persist$p,
  significance          = add_significance_stars(z_persist$p),
  interpretation = case_when(
    abs(delta_persist) < 0.005                  ~ "No meaningful difference",
    delta_persist > 0 & z_persist$p < 0.05      ~ "Cross-section shows higher persistent poverty (panel under-represents persistently poor)",
    delta_persist < 0 & z_persist$p < 0.05      ~ "Panel shows higher persistent poverty (panel over-represents persistently poor)",
    TRUE                                         ~ "Difference not statistically significant"
  )
)

cat("\nJ&VK persistent poverty representativeness test:\n")
print(jvk_persistent)

# =============================================================================
# STEP 3: Representativeness Assessment by Subgroup
# =============================================================================
cat("\n=== STEP 3: Subgroup representativeness assessment ===\n")

# Recode baseline characteristics for subgroup stratification
baseline_chars <- baseline_raw %>%
  mutate(
    female_label = if_else(
      .data[[sex_var]] == codes$sex[["female"]], "Female", "Male",
      missing = NA_character_),
    employed_label = if_else(
      .data[[labour_var]] %in% codes$employed_labour_status,
      "Employed", "Not employed", missing = NA_character_),
    formal_label = if_else(
      .data[[ss_var]] %in% codes$likely_informal_social_security_values,
      "Informal", "Formal", missing = "Formal"),
    poor_baseline_label = if_else(
      eq_income_bl < baseline_poverty_line, "Poor at baseline", "Non-poor at baseline"),
    hh_size_group = case_when(
      hh_size_bl <= 2 ~ "1-2 members",
      hh_size_bl <= 4 ~ "3-4 members",
      TRUE            ~ "5+ members"
    ),
    edu_head_label = case_when(
      .data[[edu_var]] %in% c(0, 1, 2) ~ "Primary or below",
      .data[[edu_var]] %in% c(3, 4, 5) ~ "Secondary",
      .data[[edu_var]] %in% c(6, 7, 8) ~ "Tertiary",
      TRUE ~ NA_character_
    ),
    age_head_group = cut(
      .data[[age_var]],
      breaks = c(-Inf, 34, 54, Inf),
      labels = c("Young (<=34)", "Middle (35-54)", "Older (55+)")
    )
  )

# Merge baseline characteristics onto year-4 cross-sectional data
xs_year4_with_chars <- xs_year4 %>%
  filter(.data[[year_var]] == ref_year) %>%
  left_join(
    baseline_chars %>% select(all_of(id), female_label, employed_label,
                               formal_label, poor_baseline_label,
                               hh_size_group, edu_head_label, age_head_group),
    by = id
  )

# Merge baseline characteristics onto longitudinal panel year 4
long_year4 <- panel_income %>%
  filter(.data[[year_var]] == ref_year, .data[[id]] %in% balanced_ids)

long_poverty_line_y4 <- weighted_median(
  long_year4$eq_income, long_year4[[weight_long]]
) * threshold

long_year4_with_chars <- long_year4 %>%
  mutate(poor_long = as.integer(eq_income < long_poverty_line_y4)) %>%
  left_join(
    baseline_chars %>% select(all_of(id), female_label, employed_label,
                               formal_label, poor_baseline_label,
                               hh_size_group, edu_head_label, age_head_group),
    by = id
  )

# Compute xs poverty status for year 4
xs_poverty_line_y4 <- xs_poverty_lines %>%
  filter(.data[[year_var]] == ref_year) %>%
  pull(poverty_line_xs)

xs_year4_with_chars <- xs_year4_with_chars %>%
  mutate(poor_xs_y4 = as.integer(eq_income_raw < xs_poverty_line_y4))

# Function: subgroup representativeness test
subgroup_jvk <- function(subgroup_var, subgroup_label) {
  xs_data   <- xs_year4_with_chars %>% filter(!is.na(.data[[subgroup_var]]))
  long_data <- long_year4_with_chars %>% filter(!is.na(.data[[subgroup_var]]))

  subgroups <- sort(unique(c(xs_data[[subgroup_var]], long_data[[subgroup_var]])))

  purrr::map_dfr(subgroups, function(sg) {
    xs_sub   <- xs_data %>% filter(.data[[subgroup_var]] == sg)
    long_sub <- long_data %>% filter(.data[[subgroup_var]] == sg)

    xs_n_poor  <- sum(xs_sub$poor_xs_y4, na.rm = TRUE)
    xs_n_total <- nrow(xs_sub)
    xs_rate    <- if (xs_n_total > 0) xs_n_poor / xs_n_total else NA_real_

    long_n_poor  <- sum(long_sub$poor_long, na.rm = TRUE)
    long_n_total <- nrow(long_sub)
    long_rate    <- if (long_n_total > 0) long_n_poor / long_n_total else NA_real_

    if (xs_n_total < 5 || long_n_total < 5) {
      return(tibble(
        stratification = subgroup_label,
        subgroup = sg,
        xs_n = xs_n_total, xs_poverty_rate = xs_rate,
        long_n = long_n_total, long_poverty_rate = long_rate,
        delta_p = NA_real_, z_statistic = NA_real_,
        p_value = NA_real_, significance = "",
        interpretation = "Insufficient observations"
      ))
    }

    delta_p <- xs_rate - long_rate
    z_test  <- prop_z_test(xs_n_poor, xs_n_total, long_n_poor, long_n_total)

    tibble(
      stratification    = subgroup_label,
      subgroup          = sg,
      xs_n              = xs_n_total,
      xs_poverty_rate   = xs_rate,
      long_n            = long_n_total,
      long_poverty_rate = long_rate,
      delta_p           = delta_p,
      z_statistic       = z_test$z,
      p_value           = z_test$p,
      significance      = add_significance_stars(z_test$p),
      interpretation = case_when(
        abs(delta_p) < 0.005               ~ "No meaningful difference",
        delta_p > 0 & z_test$p < 0.05      ~ "Panel under-represents poor in this subgroup",
        delta_p < 0 & z_test$p < 0.05      ~ "Panel over-represents poor in this subgroup",
        TRUE                                ~ "Difference not statistically significant"
      )
    )
  })
}

jvk_subgroup <- bind_rows(
  subgroup_jvk("employed_label",       "Employment status"),
  subgroup_jvk("formal_label",         "Formality (formal/informal)"),
  subgroup_jvk("poor_baseline_label",  "Initial poverty status"),
  subgroup_jvk("hh_size_group",        "Household size"),
  subgroup_jvk("edu_head_label",       "Head education"),
  subgroup_jvk("age_head_group",       "Head age group"),
  subgroup_jvk("female_label",         "Sex")
)

cat("\nSubgroup representativeness (ΔP = cross-sectional - longitudinal poverty rate, year 4):\n")
print(jvk_subgroup %>% select(stratification, subgroup, xs_poverty_rate,
                               long_poverty_rate, delta_p, p_value,
                               significance, interpretation), n = 40)

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
cat("\n=== Saving outputs ===\n")

write_csv(attrition_flow,       file.path(project$table_dir, "attrition_flow.csv"))
write_csv(baseline_comparison,  file.path(project$table_dir, "attrition_baseline_comparison.csv"))
write_csv(attrition_by_subgroup, file.path(project$table_dir, "attrition_by_subgroup.csv"))
write_csv(jvk_results,          file.path(project$table_dir, "jvk_representativeness_test.csv"))
write_csv(jvk_persistent,       file.path(project$table_dir, "jvk_persistent_poverty_test.csv"))
write_csv(jvk_subgroup,         file.path(project$table_dir, "jvk_subgroup_representativeness.csv"))

# ---- gt formatted tables (HTML + LaTeX) when gt is available -----------------

if (requireNamespace("gt", quietly = TRUE)) {
  gt_theme <- function(x) {
    x %>%
      gt::tab_options(
        table.font.names = "Times New Roman",
        table.font.size  = gt::px(12),
        table.border.top.color    = "black",
        table.border.top.width    = gt::px(1.5),
        table.border.bottom.color = "black",
        table.border.bottom.width = gt::px(1.5),
        column_labels.border.bottom.color = "black",
        column_labels.border.bottom.width = gt::px(1)
      )
  }

  save_gt <- function(gt_tbl, stem) {
    gt::gtsave(gt_tbl, file.path(project$table_dir, paste0(stem, ".html")))
    gt::gtsave(gt_tbl, file.path(project$table_dir, paste0(stem, ".tex")))
  }

  # --- Attrition flow table ---
  flow_gt <- attrition_flow %>%
    gt::gt() %>%
    gt::tab_header(
      title    = gt::md("**Attrition flow diagram**"),
      subtitle = gt::md("TR-SILC balanced panel, 2016--2019")
    ) %>%
    gt::fmt_number(columns = c(n_in_wave, stayed_from_prev, dropped_from_prev),
                   decimals = 0, use_seps = TRUE) %>%
    gt::fmt_percent(columns = c(wave_to_wave_retention, wave_to_wave_attrition,
                                cumulative_retention_from_wave1,
                                cumulative_attrition_from_wave1),
                    decimals = 1) %>%
    gt::cols_label(
      wave = "Wave",
      n_in_wave = "N in wave",
      stayed_from_prev = "Stayed",
      dropped_from_prev = "Dropped",
      wave_to_wave_retention = "W-to-W retention",
      wave_to_wave_attrition = "W-to-W attrition",
      cumulative_retention_from_wave1 = "Cum. retention",
      cumulative_attrition_from_wave1 = "Cum. attrition"
    ) %>%
    gt::sub_missing(missing_text = "—") %>%
    gt_theme()
  save_gt(flow_gt, "attrition_flow")

  # --- Baseline comparison table ---
  baseline_gt <- baseline_comparison %>%
    gt::gt() %>%
    gt::tab_header(
      title    = gt::md("**Baseline characteristics: stayers vs. attritors**"),
      subtitle = gt::md("Wave 1 sample. t-tests for mean differences.")
    ) %>%
    gt::cols_label(
      Characteristic = "Characteristic",
      Stayer = "Stayer",
      Attritor = "Attritor",
      Difference = "Difference",
      p_value = gt::md("p-value")
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        "Significance: \\*\\*\\* p<0.01, \\*\\* p<0.05, \\* p<0.10."
      )
    ) %>%
    gt_theme()
  save_gt(baseline_gt, "attrition_baseline_comparison")

  # --- J&VK representativeness test ---
  jvk_gt <- jvk_results %>%
    select(wave, cross_sectional_n, cross_sectional_rate,
           longitudinal_n, longitudinal_rate, delta_p,
           z_statistic, p_value, significance, interpretation) %>%
    gt::gt() %>%
    gt::tab_header(
      title    = gt::md("**Jenkins & Van Kerm representativeness test**"),
      subtitle = gt::md("Cross-sectional vs. longitudinal poverty rates by wave. 60% threshold.")
    ) %>%
    gt::fmt_number(columns = c(cross_sectional_n, longitudinal_n),
                   decimals = 0, use_seps = TRUE) %>%
    gt::fmt_percent(columns = c(cross_sectional_rate, longitudinal_rate, delta_p),
                    decimals = 2) %>%
    gt::fmt_number(columns = c(z_statistic, p_value), decimals = 3) %>%
    gt::cols_label(
      wave = "Wave",
      cross_sectional_n = "XS N",
      cross_sectional_rate = "XS poverty rate",
      longitudinal_n = "Long. N",
      longitudinal_rate = "Long. poverty rate",
      delta_p = gt::md("ΔP"),
      z_statistic = "z",
      p_value = gt::md("p-value"),
      significance = "",
      interpretation = "Interpretation"
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        paste0(
          "ΔP = cross-sectional rate − longitudinal rate. ",
          "Two-proportion z-test, H₀: ΔP = 0. ",
          "Significance: \\*\\*\\* p<0.01, \\*\\* p<0.05, \\* p<0.10. ",
          "Jenkins & Van Kerm (2017)."
        )
      )
    ) %>%
    gt_theme()
  save_gt(jvk_gt, "jvk_representativeness_test")

  # --- J&VK persistent poverty test ---
  persist_gt <- jvk_persistent %>%
    gt::gt() %>%
    gt::tab_header(
      title    = gt::md("**J&VK persistent poverty representativeness test**"),
      subtitle = gt::md("Cross-sectional vs. longitudinal persistent poverty at year 4 (2019).")
    ) %>%
    gt::fmt_number(columns = c(cross_sectional_n, longitudinal_n),
                   decimals = 0, use_seps = TRUE) %>%
    gt::fmt_percent(columns = c(cross_sectional_rate, longitudinal_rate, delta_p),
                    decimals = 2) %>%
    gt::fmt_number(columns = c(z_statistic, p_value), decimals = 3) %>%
    gt::cols_label(
      metric = "Metric",
      cross_sectional_n = "XS N",
      cross_sectional_rate = "XS rate",
      longitudinal_n = "Long. N",
      longitudinal_rate = "Long. rate",
      delta_p = gt::md("ΔP"),
      z_statistic = "z",
      p_value = gt::md("p-value"),
      significance = "",
      interpretation = "Interpretation"
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        paste0(
          "Persistent poverty = poor in 2019 AND poor in ≥ 2 of 2016--2018. ",
          "ΔP = cross-sectional rate − longitudinal rate. ",
          "Two-proportion z-test, H₀: ΔP = 0. ",
          "Significance: \\*\\*\\* p<0.01, \\*\\* p<0.05, \\* p<0.10. ",
          "Jenkins & Van Kerm (2017)."
        )
      )
    ) %>%
    gt_theme()
  save_gt(persist_gt, "jvk_persistent_poverty_test")

  # --- Subgroup representativeness ---
  subgroup_gt <- jvk_subgroup %>%
    select(stratification, subgroup, xs_n, xs_poverty_rate,
           long_n, long_poverty_rate, delta_p, p_value, significance,
           interpretation) %>%
    gt::gt(groupname_col = "stratification") %>%
    gt::tab_header(
      title    = gt::md("**Subgroup representativeness assessment (J&VK Step 3)**"),
      subtitle = gt::md("Cross-sectional vs. longitudinal poverty rates at year 4, by subgroup. 60% threshold.")
    ) %>%
    gt::fmt_number(columns = c(xs_n, long_n), decimals = 0, use_seps = TRUE) %>%
    gt::fmt_percent(columns = c(xs_poverty_rate, long_poverty_rate, delta_p),
                    decimals = 2) %>%
    gt::fmt_number(columns = p_value, decimals = 3) %>%
    gt::cols_label(
      subgroup = "Subgroup",
      xs_n = "XS N",
      xs_poverty_rate = "XS poverty rate",
      long_n = "Long. N",
      long_poverty_rate = "Long. poverty rate",
      delta_p = gt::md("ΔP"),
      p_value = gt::md("p-value"),
      significance = "",
      interpretation = "Interpretation"
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        paste0(
          "ΔP = cross-sectional rate − longitudinal rate. ",
          "Two-proportion z-test. ",
          "Significance: \\*\\*\\* p<0.01, \\*\\* p<0.05, \\* p<0.10. ",
          "Jenkins & Van Kerm (2017)."
        )
      )
    ) %>%
    gt_theme()
  save_gt(subgroup_gt, "jvk_subgroup_representativeness")

  cat("gt tables (HTML + LaTeX) saved.\n")
}

# ---- LaTeX input snippet file ------------------------------------------------

attrition_latex <- c(
  "% Auto-generated LaTeX snippets for J&VK attrition diagnostics.",
  "% Requires packages: booktabs, longtable, caption.",
  "",
  "\\input{tables/attrition_flow.tex}",
  "\\input{tables/attrition_baseline_comparison.tex}",
  "\\input{tables/jvk_representativeness_test.tex}",
  "\\input{tables/jvk_persistent_poverty_test.tex}",
  "\\input{tables/jvk_subgroup_representativeness.tex}"
)
writeLines(attrition_latex, file.path(project$table_dir, "attrition_latex_inputs.tex"))
cat("LaTeX input file saved to:", file.path(project$table_dir, "attrition_latex_inputs.tex"), "\n")

cat("\nJenkins & Van Kerm (2017) attrition diagnostics completed.\n")
