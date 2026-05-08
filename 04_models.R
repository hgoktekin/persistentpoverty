# ============================================================================
# PANEL POVERTY REGRESSION ANALYSIS
# Three-Model Approach: Random Effects + Fixed Effects + Mundlak Correction
#
# Purpose: Quantify effects of informality and household composition on poverty
# Data: 4-year balanced panel (TR-SILC 2016-2019)
# Author: Hatice Göktekin
#
# Input:  panel_16_19.dta (raw SILC micro-data)
# Output: 04_models_table.txt, 04_models_table.html, 04_models_table.tex,
#         04_model_re.rds, 04_model_fe.rds, 04_model_mundlak.rds
#
# Main models:
#   Model 1: Random Effects (Primary – comprehensive, efficient)
#   Model 2: Fixed Effects  (Robustness – checks for selection bias)
#   Model 3: Mundlak Correction (Tests RE assumptions directly)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(haven)
  library(readr)
  library(plm)
  library(lmtest)
  library(sandwich)
  library(tibble)
  library(broom)
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

source(file.path(project_root, "R", "00_config.R"))
source(file.path(project_root, "R", "01_functions.R"))

project$data_path  <- file.path(project_root, project$data_path)
project$out_dir    <- file.path(project_root, project$out_dir)
project$table_dir  <- file.path(project_root, project$table_dir)
project$model_dir  <- file.path(project_root, project$model_dir)

dir.create(project$out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(project$table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(project$model_dir, showWarnings = FALSE, recursive = TRUE)

# ---- NUTS-2 → 7 geographic regions mapping ----------------------------------
# Turkey has 26 NUTS-2 (İBBS Düzey-2) regions; we collapse them into the
# traditional 7 geographic regions for the regression specification.

nuts2_to_region <- c(
  "TR10" = "Marmara",    "TR21" = "Marmara",    "TR22" = "Marmara",
  "TR41" = "Marmara",    "TR42" = "Marmara",
  "TR31" = "Ege",        "TR32" = "Ege",        "TR33" = "Ege",
  "TR61" = "Akdeniz",    "TR62" = "Akdeniz",    "TR63" = "Akdeniz",
  "TR51" = "Ic Anadolu", "TR52" = "Ic Anadolu",
  "TR71" = "Ic Anadolu", "TR72" = "Ic Anadolu",
  "TR81" = "Karadeniz",  "TR82" = "Karadeniz",
  "TR83" = "Karadeniz",  "TR90" = "Karadeniz",
  "TRA1" = "Dogu Anadolu",  "TRA2" = "Dogu Anadolu",
  "TRB1" = "Dogu Anadolu",  "TRB2" = "Dogu Anadolu",
  "TRC1" = "Guneydogu Anadolu", "TRC2" = "Guneydogu Anadolu",
  "TRC3" = "Guneydogu Anadolu"
)

# ============================================================================
# TASK 1: DATA PREPARATION & SUMMARY STATISTICS
# ============================================================================

cat("="," TASK 1: DATA PREPARATION ", "=", "\n")
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
  add_household_context(vars, codes)

cat("Adding poverty status ...\n")
poverty_lines <- compute_poverty_lines(panel_income, vars, project$thresholds)
panel_poverty <- add_poverty_status(panel_income, poverty_lines, vars, project$thresholds)

# ---- Recode variables for regression ----------------------------------------

cat("Recoding variables for regression ...\n")

model_panel <- panel_poverty %>%
  mutate(
    id   = .data[[vars$person_id]],
    year = .data[[vars$year]],

    # Dependent variable: poverty status at 60% threshold
    poverty_status = as.integer(poor_60),

    # Informality: 1 = not registered with social security, 0 = registered
    informality = case_when(
      .data[[vars$social_security]] %in% codes$likely_informal_social_security_values ~ 1L,
      !is.na(.data[[vars$social_security]]) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Dependents = household members who are not earners
    number_of_dependents = hh_size - hh_earners_proxy,
    number_of_earners    = hh_earners_proxy,

    # Dependency ratio (+ 0.1 to avoid division by zero)
    dependency_ratio = number_of_dependents / (number_of_earners + 0.1),

    # Employment type
    employment_type = factor(
      case_when(
        .data[[vars$labour_status]] %in% c(1, 2) ~ "Employee",
        .data[[vars$labour_status]] %in% c(3, 4) ~ "Self-employed",
        .data[[vars$labour_status]] == 5          ~ "Unemployed",
        .data[[vars$labour_status]] == 7          ~ "Retired",
        .data[[vars$labour_status]] %in% c(6, 8, 9, 10) ~ "Inactive",
        TRUE ~ NA_character_
      ),
      levels = c("Employee", "Self-employed", "Unemployed", "Retired", "Inactive")
    ),

    # Education level
    education_level = factor(
      case_when(
        .data[[vars$education]] %in% c(0, 1, 2) ~ "Primary",
        .data[[vars$education]] %in% c(3, 4, 5) ~ "Secondary",
        .data[[vars$education]] %in% c(6, 7, 8) ~ "Tertiary",
        TRUE ~ NA_character_
      ),
      levels = c("Primary", "Secondary", "Tertiary")
    ),

    # Age group
    age_group = factor(
      case_when(
        .data[[vars$age]] >= 15 & .data[[vars$age]] <= 24 ~ "15-24",
        .data[[vars$age]] >= 25 & .data[[vars$age]] <= 34 ~ "25-34",
        .data[[vars$age]] >= 35 & .data[[vars$age]] <= 44 ~ "35-44",
        .data[[vars$age]] >= 45 & .data[[vars$age]] <= 54 ~ "45-54",
        .data[[vars$age]] >= 55 & .data[[vars$age]] <= 64 ~ "55-64",
        .data[[vars$age]] >= 65                            ~ "65+",
        TRUE ~ NA_character_
      ),
      levels = c("15-24", "25-34", "35-44", "45-54", "55-64", "65+")
    ),

    # Sex: 1 = female, 0 = male
    sex = factor(
      case_when(
        .data[[vars$sex]] == codes$sex[["female"]] ~ "Female",
        .data[[vars$sex]] == codes$sex[["male"]]   ~ "Male",
        TRUE ~ NA_character_
      ),
      levels = c("Male", "Female")
    ),

    # NUTS-2 → 7 geographic regions
    nuts2_code   = as.character(.data[[vars$nuts2]]),
    nuts_region  = factor(
      nuts2_to_region[nuts2_code],
      levels = c("Marmara", "Ege", "Akdeniz", "Ic Anadolu",
                 "Karadeniz", "Dogu Anadolu", "Guneydogu Anadolu")
    )
  ) %>%
  select(id, year, poverty_status, informality, dependency_ratio,
         number_of_dependents, number_of_earners,
         employment_type, education_level, age_group, sex, nuts_region) %>%
  drop_na()

# ---- Data verification ------------------------------------------------------

cat("\n--- Data Verification ---\n")
cat("Total observations:", nrow(model_panel), "\n")
cat("Unique individuals:", n_distinct(model_panel$id), "\n")
cat("Years:", paste(sort(unique(model_panel$year)), collapse = ", "), "\n\n")

cat("Poverty status distribution:\n")
print(table(model_panel$poverty_status, dnn = "poor"))

cat("\nInformality distribution:\n")
print(table(model_panel$informality, dnn = "informal"))

cat("\nEmployment type:\n")
print(table(model_panel$employment_type))

cat("\nEducation level:\n")
print(table(model_panel$education_level))

cat("\nAge group:\n")
print(table(model_panel$age_group))

cat("\nSex:\n")
print(table(model_panel$sex))

cat("\nNUTS region (7 geographic regions):\n")
print(table(model_panel$nuts_region))

# ---- Summary statistics by poverty status ------------------------------------

cat("\n--- Summary Statistics by Poverty Status ---\n")

summary_stats <- model_panel %>%
  group_by(poverty_status) %>%
  summarise(
    N = n(),
    informality_mean    = mean(informality, na.rm = TRUE),
    informality_sd      = sd(informality, na.rm = TRUE),
    dependency_ratio_mean = mean(dependency_ratio, na.rm = TRUE),
    dependency_ratio_sd   = sd(dependency_ratio, na.rm = TRUE),
    mean_age_group      = NA_real_,
    prop_female         = mean(sex == "Female", na.rm = TRUE),
    prop_primary_edu    = mean(education_level == "Primary", na.rm = TRUE),
    prop_secondary_edu  = mean(education_level == "Secondary", na.rm = TRUE),
    prop_tertiary_edu   = mean(education_level == "Tertiary", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(poverty_status = ifelse(poverty_status == 1, "Poor", "Non-poor"))

cat("\n")
print(as.data.frame(summary_stats), row.names = FALSE)

write_csv(summary_stats, file.path(project$table_dir, "04_summary_stats_by_poverty.csv"))
cat("\nSummary statistics saved.\n")

# ---- Set panel structure -----------------------------------------------------

cat("\nSetting panel structure (plm) ...\n")

pdata <- pdata.frame(model_panel, index = c("id", "year"), drop.index = FALSE)

cat("Panel dimensions: N =", pdim(pdata)$nT$n, " T =", pdim(pdata)$nT$T,
    " total obs =", pdim(pdata)$nT$N, "\n")
cat("Balanced panel:", ifelse(pdim(pdata)$balanced, "Yes", "No"), "\n\n")

# ============================================================================
# TASK 2: MODEL 1 – RANDOM EFFECTS (PRIMARY MODEL)
# ============================================================================

cat("=", " TASK 2: RANDOM EFFECTS MODEL ", "=", "\n")

re_formula <- poverty_status ~ informality + dependency_ratio +
  factor(employment_type) + factor(education_level) +
  factor(age_group) + factor(sex) +
  factor(nuts_region) + factor(year)

model_re <- plm(
  re_formula,
  data   = pdata,
  model  = "random",
  effect = "individual",
  random.method = "swar"
)

# Cluster-robust standard errors
re_robust_vcov <- vcovHC(model_re, method = "arellano", type = "HC1")
re_robust <- coeftest(model_re, vcov. = re_robust_vcov)

cat("\n--- Model 1: Random Effects (Swamy-Arora) ---\n")
cat("--- Cluster-robust standard errors (Arellano) ---\n\n")
print(re_robust)

cat("\nModel summary:\n")
print(summary(model_re))

# Extract rho (proportion of variance from individual effects)
re_summary <- summary(model_re)
re_ercomp  <- re_summary$ercomp
rho_re     <- re_ercomp$sigma2["id"] / (re_ercomp$sigma2["id"] + re_ercomp$sigma2["idios"])
cat("\nRho (proportion of variance from u_i):", round(rho_re, 4), "\n")
if (rho_re > 0.3) {
  cat("Interpretation: Rho > 0.3 — unobserved individual heterogeneity is important;\n")
  cat("  the panel structure is informative and pooled OLS would be inappropriate.\n")
} else {
  cat("Interpretation: Rho <= 0.3 — unobserved heterogeneity is moderate.\n")
}

# ============================================================================
# TASK 3: MODEL 2 – FIXED EFFECTS (ROBUSTNESS CHECK)
# ============================================================================

cat("\n=", " TASK 3: FIXED EFFECTS MODEL ", "=", "\n")

# Time-invariant variables (sex, education, age_group, nuts_region) drop automatically.
fe_formula <- poverty_status ~ informality + dependency_ratio +
  factor(employment_type) + factor(year)

model_fe <- plm(
  fe_formula,
  data   = pdata,
  model  = "within",
  effect = "individual"
)

fe_robust_vcov <- vcovHC(model_fe, method = "arellano", type = "HC1")
fe_robust <- coeftest(model_fe, vcov. = fe_robust_vcov)

cat("\n--- Model 2: Fixed Effects (Within estimator) ---\n")
cat("--- Cluster-robust standard errors (Arellano) ---\n\n")
print(fe_robust)

cat("\nModel summary:\n")
print(summary(model_fe))

fe_r2_within <- summary(model_fe)$r.squared["rsq"]
cat("\nR-squared (within):", round(fe_r2_within, 4), "\n")
cat("Note: Within R-squared uses only time-variation within individuals.\n")

# F-test for individual effects
fe_ftest <- tryCatch(pFtest(model_fe, plm(fe_formula, data = pdata, model = "pooling")),
                     error = function(e) NULL)
if (!is.null(fe_ftest)) {
  cat("\nF-test for individual effects:\n")
  print(fe_ftest)
  if (fe_ftest$p.value < 0.05) {
    cat("Interpretation: Significant (p < 0.05) — individual fixed effects are jointly significant.\n")
  }
}

# ---- Coefficient comparison: RE vs FE ---------------------------------------

cat("\n--- Coefficient Comparison: Informality ---\n")
re_coef_inf <- coef(model_re)["informality"]
fe_coef_inf <- coef(model_fe)["informality"]
diff_inf    <- re_coef_inf - fe_coef_inf

cat("  RE coefficient:", round(re_coef_inf, 4), "\n")
cat("  FE coefficient:", round(fe_coef_inf, 4), "\n")
cat("  Difference (RE - FE):", round(diff_inf, 4), "\n")

if (abs(fe_coef_inf) < abs(re_coef_inf) * 0.5) {
  cat("  The FE estimate is substantially smaller than RE, suggesting positive\n")
  cat("  selection bias: ambitious/able individuals are more likely to formalise\n")
  cat("  and also less likely to be poor for other reasons.\n")
} else if (abs(diff_inf) < abs(re_coef_inf) * 0.2) {
  cat("  The FE and RE estimates are similar, suggesting minimal selection bias.\n")
  cat("  The RE estimator appears credible for informality.\n")
} else {
  cat("  The difference of", round(diff_inf, 4), "suggests moderate selection bias.\n")
  cat("  Part of the RE association comes from unobserved individual characteristics.\n")
}

cat("\n--- Coefficient Comparison: Dependency Ratio ---\n")
re_coef_dep <- coef(model_re)["dependency_ratio"]
fe_coef_dep <- coef(model_fe)["dependency_ratio"]
diff_dep    <- re_coef_dep - fe_coef_dep

cat("  RE coefficient:", round(re_coef_dep, 4), "\n")
cat("  FE coefficient:", round(fe_coef_dep, 4), "\n")
cat("  Difference (RE - FE):", round(diff_dep, 4), "\n")

# ============================================================================
# TASK 4: MODEL 3 – MUNDLAK CORRECTION (ASSUMPTION TESTING)
# ============================================================================

cat("\n=", " TASK 4: MUNDLAK CORRECTION ", "=", "\n")

# Create within-means (time-averages) for each individual
cat("Computing Mundlak correction terms (within-person means) ...\n")

mundlak_means <- model_panel %>%
  mutate(employment_numeric = as.numeric(employment_type)) %>%
  group_by(id) %>%
  summarise(
    avg_informality      = mean(informality, na.rm = TRUE),
    avg_dependency_ratio = mean(dependency_ratio, na.rm = TRUE),
    avg_employment_type  = mean(employment_numeric, na.rm = TRUE),
    .groups = "drop"
  )

model_panel_mundlak <- model_panel %>%
  left_join(mundlak_means, by = "id")

pdata_mundlak <- pdata.frame(model_panel_mundlak, index = c("id", "year"), drop.index = FALSE)

mundlak_formula <- poverty_status ~ informality + dependency_ratio +
  factor(employment_type) + factor(education_level) +
  factor(age_group) + factor(sex) +
  factor(nuts_region) + factor(year) +
  avg_informality + avg_dependency_ratio + avg_employment_type

model_mundlak <- plm(
  mundlak_formula,
  data   = pdata_mundlak,
  model  = "random",
  effect = "individual",
  random.method = "swar"
)

mundlak_robust_vcov <- vcovHC(model_mundlak, method = "arellano", type = "HC1")
mundlak_robust <- coeftest(model_mundlak, vcov. = mundlak_robust_vcov)

cat("\n--- Model 3: Random Effects with Mundlak Correction ---\n")
cat("--- Cluster-robust standard errors (Arellano) ---\n\n")
print(mundlak_robust)

cat("\nModel summary:\n")
print(summary(model_mundlak))

# ---- Interpret Mundlak terms -------------------------------------------------

cat("\n--- Mundlak Correction Terms ---\n")

mundlak_vars <- c("avg_informality", "avg_dependency_ratio", "avg_employment_type")
mundlak_labels <- c("Average informality", "Average dependency ratio", "Average employment type")

for (i in seq_along(mundlak_vars)) {
  v <- mundlak_vars[i]
  lab <- mundlak_labels[i]

  idx <- which(rownames(mundlak_robust) == v)
  if (length(idx) == 0) {
    cat("\n", lab, ": not found in model (may have been dropped).\n")
    next
  }

  coef_val <- mundlak_robust[idx, 1]
  p_val    <- mundlak_robust[idx, 4]

  cat("\n", lab, ":\n")
  cat("  Coefficient:", round(coef_val, 4), "\n")
  cat("  p-value:", round(p_val, 4), "\n")

  if (p_val < 0.05) {
    cat("  => SIGNIFICANT: Individuals with higher average", tolower(lab),
        "are more poverty-prone\n")
    cat("     even controlling for current values. This suggests selection bias:\n")
    cat("     the RE estimate for this variable may overstate the true causal effect.\n")
  } else {
    cat("  => NOT significant: No evidence of correlation between unobserved\n")
    cat("     heterogeneity and this variable. RE assumption appears reasonable.\n")
  }
}

# ============================================================================
# TASK 5: STATISTICAL TESTS FOR MODEL COMPARISON
# ============================================================================

cat("\n=", " TASK 5: HAUSMAN TEST ", "=", "\n")

# Re-estimate both models on the same (reduced) formula for valid comparison.
# The Hausman test requires that both models use the same specification for
# the variables that appear in both.
hausman_formula <- poverty_status ~ informality + dependency_ratio +
  factor(employment_type) + factor(year)

model_re_hausman <- plm(hausman_formula, data = pdata, model = "random",
                        effect = "individual", random.method = "swar")
model_fe_hausman <- plm(hausman_formula, data = pdata, model = "within",
                        effect = "individual")

hausman_test <- tryCatch(
  phtest(model_fe_hausman, model_re_hausman),
  error = function(e) {
    cat("Hausman test error:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(hausman_test)) {
  cat("\n--- Hausman Test: FE vs RE ---\n")
  print(hausman_test)

  cat("\nTest statistic:", round(hausman_test$statistic, 4), "\n")
  cat("Degrees of freedom:", hausman_test$parameter, "\n")
  cat("p-value:", format.pval(hausman_test$p.value, digits = 4), "\n\n")

  if (hausman_test$p.value < 0.05) {
    cat("Interpretation: Reject H0 (p < 0.05). The RE assumption that unobserved\n")
    cat("  heterogeneity is uncorrelated with the regressors is violated.\n")
    cat("  FE estimates are safer for causal interpretation.\n")
    cat("  However, RE remains useful for estimating time-invariant effects.\n")
  } else {
    cat("Interpretation: Cannot reject H0 (p > 0.05). No evidence that the RE\n")
    cat("  assumption is violated. RE is more efficient and preferred.\n")
  }
}

# ============================================================================
# TASK 6: PUBLICATION-QUALITY RESULTS TABLE
# ============================================================================

cat("\n=", " TASK 6: REGRESSION TABLES ", "=", "\n")

# Stargazer requires the stargazer package; fall back to manual table if absent.
if (requireNamespace("stargazer", quietly = TRUE)) {

  cat("Generating publication tables with stargazer ...\n")

  # Plain text
  stargazer::stargazer(
    model_re, model_fe, model_mundlak,
    type = "text",
    title = "Panel Regression Models: Poverty Determinants (TR-SILC 2016-2019)",
    column.labels = c("Random Effects", "Fixed Effects", "Mundlak Correction"),
    dep.var.labels = "Poverty Status (1 = poor)",
    covariate.labels = c(
      "Informality (not registered)",
      "Dependency Ratio",
      "Self-employed", "Unemployed", "Retired", "Inactive",
      "Secondary education", "Tertiary education",
      "25-34", "35-44", "45-54", "55-64", "65+",
      "Female",
      "Ege", "Akdeniz", "Ic Anadolu", "Karadeniz",
      "Dogu Anadolu", "Guneydogu Anadolu",
      "2017", "2018", "2019",
      "Avg Informality (Mundlak)", "Avg Dependency Ratio (Mundlak)",
      "Avg Employment Type (Mundlak)"
    ),
    se = list(
      sqrt(diag(re_robust_vcov)),
      sqrt(diag(fe_robust_vcov)),
      sqrt(diag(mundlak_robust_vcov))
    ),
    omit.stat = c("f", "ser"),
    add.lines = list(
      c("Individual effects", "Random", "Fixed", "Random + Mundlak"),
      c("Time effects", "Yes", "Yes", "Yes"),
      c("Cluster-robust SE", "Yes", "Yes", "Yes")
    ),
    notes = "Robust standard errors clustered at individual level (Arellano).",
    notes.align = "l",
    out = file.path(project$table_dir, "04_models_table.txt")
  )

  # HTML
  stargazer::stargazer(
    model_re, model_fe, model_mundlak,
    type = "html",
    title = "Panel Regression Models: Poverty Determinants (TR-SILC 2016-2019)",
    column.labels = c("Random Effects", "Fixed Effects", "Mundlak Correction"),
    dep.var.labels = "Poverty Status (1 = poor)",
    covariate.labels = c(
      "Informality (not registered)",
      "Dependency Ratio",
      "Self-employed", "Unemployed", "Retired", "Inactive",
      "Secondary education", "Tertiary education",
      "25-34", "35-44", "45-54", "55-64", "65+",
      "Female",
      "Ege", "Akdeniz", "Ic Anadolu", "Karadeniz",
      "Dogu Anadolu", "Guneydogu Anadolu",
      "2017", "2018", "2019",
      "Avg Informality (Mundlak)", "Avg Dependency Ratio (Mundlak)",
      "Avg Employment Type (Mundlak)"
    ),
    se = list(
      sqrt(diag(re_robust_vcov)),
      sqrt(diag(fe_robust_vcov)),
      sqrt(diag(mundlak_robust_vcov))
    ),
    omit.stat = c("f", "ser"),
    add.lines = list(
      c("Individual effects", "Random", "Fixed", "Random + Mundlak"),
      c("Time effects", "Yes", "Yes", "Yes"),
      c("Cluster-robust SE", "Yes", "Yes", "Yes")
    ),
    notes = "Robust standard errors clustered at individual level (Arellano).",
    notes.align = "l",
    out = file.path(project$table_dir, "04_models_table.html")
  )

  # LaTeX
  stargazer::stargazer(
    model_re, model_fe, model_mundlak,
    type = "latex",
    title = "Panel Regression Models: Poverty Determinants (TR-SILC 2016--2019)",
    label = "tab:panel-models",
    column.labels = c("Random Effects", "Fixed Effects", "Mundlak Correction"),
    dep.var.labels = "Poverty Status (1 = poor)",
    covariate.labels = c(
      "Informality (not registered)",
      "Dependency Ratio",
      "Self-employed", "Unemployed", "Retired", "Inactive",
      "Secondary education", "Tertiary education",
      "25-34", "35-44", "45-54", "55-64", "65+",
      "Female",
      "Ege", "Akdeniz", "\\Ic{} Anadolu", "Karadeniz",
      "Do\\u{g}u Anadolu", "G\\\"uneydogu Anadolu",
      "2017", "2018", "2019",
      "Avg Informality (Mundlak)", "Avg Dependency Ratio (Mundlak)",
      "Avg Employment Type (Mundlak)"
    ),
    se = list(
      sqrt(diag(re_robust_vcov)),
      sqrt(diag(fe_robust_vcov)),
      sqrt(diag(mundlak_robust_vcov))
    ),
    omit.stat = c("f", "ser"),
    add.lines = list(
      c("Individual effects", "Random", "Fixed", "Random + Mundlak"),
      c("Time effects", "Yes", "Yes", "Yes"),
      c("Cluster-robust SE", "Yes", "Yes", "Yes")
    ),
    notes = "Robust standard errors clustered at individual level (Arellano).",
    notes.align = "l",
    out = file.path(project$table_dir, "04_models_table.tex")
  )

  cat("Tables saved: 04_models_table.txt, .html, .tex\n")

} else {
  cat("stargazer package not available. Saving manual coefficient table ...\n")

  # Fallback: tidy coefficient table
  tidy_model <- function(model, robust_vcov, model_name) {
    cf  <- coef(model)
    se  <- sqrt(diag(robust_vcov))
    z   <- cf / se
    p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
    tibble(
      model = model_name,
      term  = names(cf),
      estimate  = cf,
      std_error = se,
      z_value   = z,
      p_value   = p,
      ci_low    = cf - 1.96 * se,
      ci_high   = cf + 1.96 * se
    )
  }

  coef_table <- bind_rows(
    tidy_model(model_re,      re_robust_vcov,      "Random Effects"),
    tidy_model(model_fe,      fe_robust_vcov,      "Fixed Effects"),
    tidy_model(model_mundlak, mundlak_robust_vcov, "Mundlak Correction")
  )

  write_csv(coef_table, file.path(project$table_dir, "04_models_coefficients.csv"))
  cat("Coefficient table saved: 04_models_coefficients.csv\n")
}

# ============================================================================
# TASK 7: COEFFICIENT COMPARISON & INTERPRETATION
# ============================================================================

cat("\n=", " TASK 7: COEFFICIENT COMPARISON TABLE ", "=", "\n")

make_comparison_row <- function(var_name, label, model_re, re_vcov,
                                model_fe, fe_vcov,
                                model_mk, mk_vcov) {
  extract <- function(model, vcov) {
    cf <- coef(model)
    if (var_name %in% names(cf)) {
      se <- sqrt(diag(vcov))[var_name]
      return(list(est = cf[var_name], se = se,
                  ci_low = cf[var_name] - 1.96 * se,
                  ci_high = cf[var_name] + 1.96 * se))
    }
    list(est = NA, se = NA, ci_low = NA, ci_high = NA)
  }

  re <- extract(model_re, re_vcov)
  fe <- extract(model_fe, fe_vcov)
  mk <- extract(model_mk, mk_vcov)

  tibble(
    variable    = label,
    RE_coef     = re$est,  RE_se = re$se,
    RE_ci       = ifelse(is.na(re$est), NA_character_,
                         sprintf("[%.4f, %.4f]", re$ci_low, re$ci_high)),
    FE_coef     = fe$est,  FE_se = fe$se,
    FE_ci       = ifelse(is.na(fe$est), NA_character_,
                         sprintf("[%.4f, %.4f]", fe$ci_low, fe$ci_high)),
    Mundlak_coef = mk$est, Mundlak_se = mk$se,
    Mundlak_ci   = ifelse(is.na(mk$est), NA_character_,
                          sprintf("[%.4f, %.4f]", mk$ci_low, mk$ci_high)),
    diff_RE_FE  = ifelse(is.na(re$est) | is.na(fe$est), NA_real_,
                         re$est - fe$est)
  )
}

comparison_table <- bind_rows(
  make_comparison_row("informality",      "Informality (not registered)",
                      model_re, re_robust_vcov, model_fe, fe_robust_vcov,
                      model_mundlak, mundlak_robust_vcov),
  make_comparison_row("dependency_ratio", "Dependency Ratio",
                      model_re, re_robust_vcov, model_fe, fe_robust_vcov,
                      model_mundlak, mundlak_robust_vcov)
)

cat("\n")
print(as.data.frame(comparison_table), row.names = FALSE)

write_csv(comparison_table, file.path(project$table_dir, "04_coefficient_comparison.csv"))

# ---- Interpretation ----------------------------------------------------------

cat("\n--- Interpretation ---\n")

re_inf <- coef(model_re)["informality"]
fe_inf <- coef(model_fe)["informality"]
re_dep <- coef(model_re)["dependency_ratio"]
fe_dep <- coef(model_fe)["dependency_ratio"]

cat("\nInformality:\n")
cat(sprintf(
  "  Workers not registered with social security are %.1f percentage points\n",
  re_inf * 100
))
cat("  more likely to be poor in the RE model. The FE estimate (",
    sprintf("%.1f pp", fe_inf * 100), ") ",
    ifelse(abs(fe_inf) < abs(re_inf),
           "is smaller, suggesting that part of the association comes from\n  selection (unobserved characteristics) rather than causation alone.\n",
           "is similar, suggesting minimal selection bias.\n"),
    sep = "")

bias_pct <- ifelse(re_inf != 0, (1 - fe_inf / re_inf) * 100, NA)
if (!is.na(bias_pct) && abs(bias_pct) > 5) {
  cat(sprintf("  Approximately %.0f%% of the RE effect may reflect selection.\n", abs(bias_pct)))
}

cat("\nDependency Ratio:\n")
cat(sprintf(
  "  Each additional dependent increases poverty risk by %.1f percentage points (RE).\n",
  re_dep * 100
))
cat(sprintf(
  "  The FE estimate is %.1f pp. ", fe_dep * 100
))
if (abs(re_dep - fe_dep) < abs(re_dep) * 0.2) {
  cat("The similarity suggests little selection on household composition.\n")
} else {
  cat("The difference suggests some selection on household composition.\n")
}

# ============================================================================
# TASK 8: EFFECT SIZE DOCUMENTATION
# ============================================================================

cat("\n=", " TASK 8: EFFECT SIZE DOCUMENTATION ", "=", "\n")

cat("\n--- RANDOM EFFECTS (Primary for policy) ---\n")

re_cf <- coef(model_re)

cat(sprintf(
  "  Non-registration (informality) increases poverty probability by %.1f pp,\n  holding other factors constant.\n",
  re_cf["informality"] * 100
))
cat(sprintf(
  "  Each additional household dependent increases poverty risk by %.1f pp.\n",
  re_cf["dependency_ratio"] * 100
))

if ("factor(sex)Female" %in% names(re_cf)) {
  cat(sprintf(
    "  Female individuals are %.1f pp %s likely to be poor.\n",
    abs(re_cf["factor(sex)Female"]) * 100,
    ifelse(re_cf["factor(sex)Female"] > 0, "more", "less")
  ))
}

if ("factor(education_level)Secondary" %in% names(re_cf)) {
  cat(sprintf(
    "  Secondary education reduces poverty probability by %.1f pp\n  compared to primary education.\n",
    abs(re_cf["factor(education_level)Secondary"]) * 100
  ))
}

if ("factor(education_level)Tertiary" %in% names(re_cf)) {
  cat(sprintf(
    "  Tertiary education reduces poverty probability by %.1f pp\n  compared to primary education.\n",
    abs(re_cf["factor(education_level)Tertiary"]) * 100
  ))
}

cat("\n--- FIXED EFFECTS (Causal interpretation) ---\n")

fe_cf <- coef(model_fe)
cat(sprintf(
  "  The FE specification controls for all time-invariant individual\n  characteristics. When a worker transitions from formal to informal\n  status, their poverty probability changes by %.1f pp.\n",
  fe_cf["informality"] * 100
))
cat("  This suggests formalization itself — not just selection of different\n")
cat("  types into formal work — provides poverty protection.\n")

cat("\n--- MUNDLAK TEST (Assumption validity) ---\n")

for (v in mundlak_vars) {
  idx <- which(rownames(mundlak_robust) == v)
  if (length(idx) == 0) next

  p_val <- mundlak_robust[idx, 4]
  if (p_val < 0.05) {
    cat(sprintf("  %s: SIGNIFICANT (p = %.4f)\n", v, p_val))
    cat("    Individuals with persistently high values also have unobserved\n")
    cat("    poverty-increasing characteristics. RE may overstate the causal effect.\n")
  } else {
    cat(sprintf("  %s: NOT significant (p = %.4f)\n", v, p_val))
    cat("    RE assumption appears reasonable for this variable.\n")
  }
}

# ============================================================================
# TASK 9: DIAGNOSTIC CHECKS & ROBUSTNESS
# ============================================================================

cat("\n=", " TASK 9: DIAGNOSTIC CHECKS ", "=", "\n")

cat("\n1. Rho from RE model:", round(rho_re, 4), "\n")
cat("   (Proportion of total variance from individual-specific u_i)\n")
if (rho_re > 0.3) {
  cat("   Panel structure is important; pooled OLS is inappropriate.\n")
} else {
  cat("   Moderate heterogeneity; pooling may be acceptable but panel is preferred.\n")
}

cat("\n2. Within R-squared from FE model:", round(fe_r2_within, 4), "\n")
cat("   FE explains", sprintf("%.1f%%", fe_r2_within * 100),
    "of within-person variation in poverty status.\n")

if (!is.null(hausman_test)) {
  cat("\n3. Hausman test p-value:", format.pval(hausman_test$p.value, digits = 4), "\n")
  if (hausman_test$p.value < 0.05) {
    cat("   FE assumption (E(u_i|X) != 0) is supported; FE results are safer.\n")
  } else {
    cat("   No evidence against RE; RE is more efficient and preferred.\n")
  }
}

cat("\n4. Number of observations used:\n")
cat("   RE model:", nobs(model_re), "\n")
cat("   FE model:", nobs(model_fe), "\n")
cat("   Mundlak model:", nobs(model_mundlak), "\n")

cat("\n5. Panel balance:\n")
pdims <- pdim(pdata)
cat("   Balanced:", ifelse(pdims$balanced, "Yes", "No"), "\n")
cat("   Individuals:", pdims$nT$n, "  Periods:", pdims$nT$T,
    "  Total:", pdims$nT$N, "\n")

if (!pdims$balanced) {
  obs_per_id <- model_panel %>% count(id, name = "n_years")
  cat("   Distribution of observations per individual:\n")
  print(table(obs_per_id$n_years))
}

# ============================================================================
# TASK 10: SAVE OUTPUTS FOR FURTHER ANALYSIS
# ============================================================================

cat("\n=", " TASK 10: SAVING OUTPUTS ", "=", "\n")

saveRDS(model_re,      file.path(project$model_dir, "04_model_re.rds"))
saveRDS(model_fe,      file.path(project$model_dir, "04_model_fe.rds"))
saveRDS(model_mundlak, file.path(project$model_dir, "04_model_mundlak.rds"))
cat("Model objects saved: 04_model_re.rds, 04_model_fe.rds, 04_model_mundlak.rds\n")

# Predicted poverty probabilities from RE model
pdata$pred_poverty_re <- fitted(model_re)
pred_data <- as.data.frame(pdata) %>%
  select(id, year, poverty_status, pred_poverty_re)
write_csv(pred_data, file.path(project$out_dir, "04_predictions_re.csv"))
cat("RE predictions saved: 04_predictions_re.csv\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Panel poverty regression analysis completed.\n")
cat("Models:   ", normalizePath(project$model_dir), "\n")
cat("Tables:   ", normalizePath(project$table_dir), "\n")
cat("Outputs:  ", normalizePath(project$out_dir), "\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# To reload models later:
# model_re      <- readRDS("models/04_model_re.rds")
# model_fe      <- readRDS("models/04_model_fe.rds")
# model_mundlak <- readRDS("models/04_model_mundlak.rds")
