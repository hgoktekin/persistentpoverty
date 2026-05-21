# ============================================================================
# PANEL POVERTY REGRESSION ANALYSIS
# Three-Model Approach: Random Effects + Fixed Effects + Mundlak Correction
#
# Purpose: Quantify effects of informality and household composition on poverty
# Data: 4-year balanced panel (TR-SILC 2016-2019)
# Author: Hatice Göktekin
#
# Input:  panel_16_19.dta (raw SILC micro-data)
# Output: 04_models_table.html, 04_models_table.tex,
#         04_model_re.rds, 04_model_fe.rds, 04_model_mundlak.rds
#
# Dependent variable: poor_60 (binary: 1 = poor at 60% threshold, 0 = non-poor)
#
# Main models:
#   Model 1: Random Effects (Primary – comprehensive, efficient)
#   Model 2: Fixed Effects  (Robustness – controls for all time-invariant heterogeneity)
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

source(file.path(project_root,  "00_config.R"))
source(file.path(project_root,  "01_functions.R"))

project$data_path  <- file.path(project_root, project$data_path)
project$out_dir    <- file.path(project_root, project$out_dir)
project$table_dir  <- file.path(project_root, project$table_dir)
project$model_dir  <- file.path(project_root, project$model_dir)

dir.create(project$out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(project$table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(project$model_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# DATA PREPARATION
# ============================================================================

cat("Reading raw SILC panel data ...\n")
raw <- read_dta('/Users/haticegoktekin/Desktop/phd application/lisans tez/panel_16_19.dta')

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

# ---- Prepare model data frame -----------------------------------------------

cat("Recoding variables for regression ...\n")

model_panel <- panel_poverty %>%
  mutate(
    id   = .data[[vars$person_id]],
    year = .data[[vars$year]],

    # Dependent variable: poverty status at 60% threshold
    poverty_status = as.integer(poor_60),

    # Binary variables (already 0/1 from add_household_context)
    # informal_status, female, earner_loss are already in the data

    # Employment type
    labour_recoded = factor(
      case_when(
        .data[[vars$labour_status]] %in% c(1, 2)          ~ "Employee",
        .data[[vars$labour_status]] %in% c(3, 4)          ~ "Self-employed",
        .data[[vars$labour_status]] == 5                   ~ "Unemployed",
        .data[[vars$labour_status]] == 7                   ~ "Retired",
        .data[[vars$labour_status]] %in% c(6, 8, 9, 10) ~ "Inactive",
        TRUE ~ NA_character_ ),
      levels = c("Employee", "Self-employed","Unemployed", "Retired", "Inactive")
    ),
    # Education level
    education_recoded = factor(
      case_when(
        .data[[vars$education]] %in% c(0, 1, 2) ~ "Primary or below",
        .data[[vars$education]] %in% c(3, 4, 5) ~ "Secondary",
        .data[[vars$education]] == 6 ~ "Tertiary",
        TRUE ~ NA_character_),
      levels = c("Primary or below", "Secondary", "Tertiary")),
    # Age group
    age_group = cut(
      .data[[vars$age]],
      breaks = c(-Inf, 17, 24, 34, 44, 54, 64, Inf),
      labels = c("0-17", "18-24", "25-34", "35-44", "45-54", "55-64", "65+")
    ) ) %>%
  select(id, year, poverty_status,
         informal_status, female, dependency_ratio,
         dependency_ratio_oecd, log_dependency_ratio_oecd,
         labour_recoded, education_recoded, age_group)

# ---- Diagnose missing values before dropping --------------------------------

cat("\n--- Missing values per variable (before drop_na) ---\n")
na_counts <- colSums(is.na(model_panel))
for (v in names(na_counts)) {
  if (na_counts[v] > 0) {
    cat(sprintf("  %-30s %d NAs (%.1f%%)\n", v, na_counts[v],
                100 * na_counts[v] / nrow(model_panel)))
  }
}
if (all(na_counts == 0)) cat("  No missing values.\n")
cat("  Total rows before drop_na:", nrow(model_panel), "\n")

model_panel <- model_panel %>% drop_na() %>% droplevels()
cat("  Total rows after  drop_na:", nrow(model_panel), "\n")

# ---- Check factor levels ----------------------------------------------------

cat("\n--- Factor levels after drop_na ---\n")
factor_vars <- c("employment_type", "education_level", "age_group")
usable_factors <- character(0)

for (v in factor_vars) {
  nlev <- nlevels(model_panel[[v]])
  cat(sprintf("  %-20s %d levels: %s\n", v, nlev,
              paste(levels(model_panel[[v]]), collapse = ", ")))
  if (nlev >= 2) {
    usable_factors <- c(usable_factors, v)
  } else {
    cat("    ** WARNING: only", nlev, "level — excluded from regression **\n")
  }
}

# ---- Data verification ------------------------------------------------------

cat("\n--- Data Verification ---\n")
cat("Total observations:", nrow(model_panel), "\n")
cat("Unique individuals:", n_distinct(model_panel$id), "\n")
cat("Years:", paste(sort(unique(model_panel$year)), collapse = ", "), "\n\n")

cat("Poverty status distribution:\n")
print(table(model_panel$poverty_status, dnn = "poor"))

cat("\nInformal status distribution:\n")
print(table(model_panel$informal_status, dnn = "informal"))

cat("\nFemale distribution:\n")
print(table(model_panel$female, dnn = "female"))

cat("\nEmployment type:\n")
print(table(model_panel$labour_recoded))

cat("\nEducation level:\n")
print(table(model_panel$education_recoded))

cat("\n--- Summary Statistics by Poverty Status ---\n")

summary_stats <- model_panel %>%
  group_by(poverty_status) %>%
  summarise(
    N = n(),
    informal_mean    = mean(informal_status, na.rm = TRUE),
    informal_sd      = sd(informal_status, na.rm = TRUE),
    dep_ratio_oecd_mean = mean(dependency_ratio_oecd, na.rm = TRUE),
    dep_ratio_oecd_sd   = sd(dependency_ratio_oecd, na.rm = TRUE),
    prop_female         = mean(female, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    poverty_status = ifelse(poverty_status == 1, "Poor", "Non-poor"),
    dep_ratio_oecd_mean = round(dep_ratio_oecd_mean, 2),
    dep_ratio_oecd_sd   = round(dep_ratio_oecd_sd, 2) )
cat("\n")
print(as.data.frame(summary_stats), row.names = FALSE)
write_csv(summary_stats, file.path(project$table_dir, "04_summary_stats_by_poverty.csv"))

# ---- Set panel structure -----------------------------------------------------

cat("\nSetting panel structure (plm) ...\n")
pdata <- pdata.frame(model_panel, index = c("id", "year"), drop.index = FALSE)

cat("Panel dimensions: N =", pdim(pdata)$nT$n, " T =", pdim(pdata)$nT$T,
    " total obs =", pdim(pdata)$nT$N, "\n")
cat("Balanced panel:", ifelse(pdim(pdata)$balanced, "Yes", "No"), "\n\n")

# ---- Build formulas dynamically ----------------------------------------------

# Core time-varying covariates (appear in all models)
core_vars <- c("informal_status", "dependency_ratio")

# Factor terms (only include factors with 2+ levels)
re_factor_terms <- paste0("factor(", usable_factors, ")")

# RE formula: core + factors + female + year
re_rhs <- paste(c(core_vars, re_factor_terms, "female", "factor(year)"),
                collapse = " + ")

# FE formula: only time-varying variables; time-invariant (education, age, female) drop
fe_time_varying <- intersect(c("employment_type"), usable_factors)
fe_factor_terms <- if (length(fe_time_varying) > 0) {
  paste0("factor(", fe_time_varying, ")")
} else {
  character(0)
}
fe_rhs <- paste(c(core_vars, fe_factor_terms, "factor(year)"), collapse = " + ")

# ============================================================================
# MODEL 1: RANDOM EFFECTS
# ============================================================================

cat("="," MODEL 1: RANDOM EFFECTS ", "=", "\n")
# no gender in the formula 
re_formula <- as.formula(paste("poverty_status ~", re_rhs))
cat("RE formula:", deparse(re_formula), "\n")

model_re <- plm(
  re_formula,
  data   = pdata,
  model  = "random",
  effect = "individual",
  random.method = "swar"
)

re_robust_vcov <- vcovHC(model_re, method = "arellano", type = "HC1")
re_robust <- coeftest(model_re, vcov. = re_robust_vcov)

cat("\n--- Model 1: Random Effects (Swamy-Arora) ---\n")
cat("--- Cluster-robust standard errors (Arellano) ---\n\n")
print(re_robust)

re_summary <- summary(model_re)
re_ercomp  <- re_summary$ercomp
rho_re     <- re_ercomp$sigma2["id"] / (re_ercomp$sigma2["id"] + re_ercomp$sigma2["idios"])
cat("\nRho (proportion of variance from u_i):", round(rho_re, 4), "\n")

# ============================================================================
# MODEL 2: FIXED EFFECTS
# ============================================================================

cat("\n=", " MODEL 2: FIXED EFFECTS ", "=", "\n")

fe_formula <- as.formula(paste("poverty_status ~", fe_rhs))
cat("FE formula:", deparse(fe_formula), "\n")

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

fe_r2_within <- summary(model_fe)$r.squared["rsq"]
cat("\nR-squared (within):", round(fe_r2_within, 4), "\n")

# F-test for individual effects
fe_ftest <- tryCatch(
  pFtest(model_fe, plm(fe_formula, data = pdata, model = "pooling")),
  error = function(e) NULL)
if (!is.null(fe_ftest)) {
  cat("\nF-test for individual effects:\n")
  print(fe_ftest)
}

# ============================================================================
# MODEL 3: MUNDLAK CORRECTION
# ============================================================================

cat("\n=", " MODEL 3: MUNDLAK CORRECTION ", "=", "\n")
# controlling terms 
# Try estimation; if singular, progressively drop Mundlak terms
mundlak_terms <- c("avg_informal", "avg_dep_ratio", "avg_employment")

# Check collinearity before estimating
cat("\nMundlak term correlations with originals:\n")
cor_inf <- cor(model_panel_mundlak$informal_status, model_panel_mundlak$avg_informal, use = "complete.obs")
cor_dep <- cor(model_panel_mundlak$dependency_ratio, model_panel_mundlak$avg_dep_ratio, use = "complete.obs")
emp_num <- as.numeric(model_panel_mundlak$labour_recoded)
cor_emp <- cor(emp_num, model_panel_mundlak$avg_employment, use = "complete.obs")
cat(sprintf("  informal_status vs avg_informal: %.4f\n", cor_inf))
cat(sprintf("  dependency_ratio vs avg_dep_ratio: %.4f\n", cor_dep))
cat(sprintf("  employment vs avg_employment: %.4f\n", cor_emp))

# corr coeff : 0.93 ; 0.92 ; 0.94

# Mundlak correction terms: within-person means of time-varying variables
mundlak_means <- model_panel %>%
  mutate(employment_numeric = as.numeric(labour_recoded)) %>%
  group_by(id) %>%
  summarise(
    avg_informal       = mean(informal_status, na.rm = TRUE),
    avg_dep_ratio = mean(dependency_ratio, na.rm = TRUE),
    avg_employment     = mean(employment_numeric, na.rm = TRUE),
    .groups = "drop")

model_panel_mundlak <- model_panel %>% left_join(mundlak_means, by = "id")
pdata_mundlak <- pdata.frame(model_panel_mundlak, index = c("id", "year"), drop.index = FALSE)

#mundlak_terms <- c("avg_informal", "avg_dep_ratio", "avg_employment")
mundlak_rhs <- paste(c(re_rhs, mundlak_terms), collapse = " + ")
mundlak_formula <- as.formula(paste("poverty_status ~", mundlak_rhs))
cat("Mundlak formula:", deparse(mundlak_formula), "\n")

model_mundlak <- tryCatch(
  plm(mundlak_formula, data = pdata_mundlak,
      model = "random", effect = "individual", random.method = "swar"),
  error = function(e) {
    cat("  Swar failed, trying walhus...\n")
    tryCatch(
      plm(mundlak_formula, data = pdata_mundlak,
          model = "random", effect = "individual", random.method = "walhus"),
      error = function(e2) {
        cat("  Walhus failed, trying amemiya...\n")
        plm(mundlak_formula, data = pdata_mundlak,
            model = "random", effect = "individual", random.method = "amemiya")
      }
    )
  }
)

mundlak_robust_vcov <- vcovHC(model_mundlak, method = "arellano", type = "HC1")
mundlak_robust <- coeftest(model_mundlak, vcov. = mundlak_robust_vcov)

cat("\n--- Model 3: Random Effects with Mundlak Correction ---\n")
cat("--- Cluster-robust standard errors (Arellano) ---\n\n")
print(mundlak_robust)

# ---- Mundlak terms interpretation -------------------------------------------

cat("\n--- Mundlak Correction Terms ---\n")
for (v in mundlak_terms) {
  idx <- which(rownames(mundlak_robust) == v)
  if (length(idx) == 0) next
  p_val <- mundlak_robust[idx, 4]
  stars <- add_significance_stars(p_val)
  cat(sprintf("  %-30s coef = %7.4f  p = %s%s\n",
              v, mundlak_robust[idx, 1],
              ifelse(p_val < 0.001, "<0.001", sprintf("%.4f", p_val)), stars))
}

# ============================================================================
# HAUSMAN TEST: FE vs RE
# ============================================================================

cat("\n=", " HAUSMAN TEST ", "=", "\n")

hausman_formula <- as.formula(paste("poverty_status ~", fe_rhs))
model_re_h <- plm(hausman_formula, data = pdata, model = "random",
                   effect = "individual", random.method = "walhus")
model_fe_h <- plm(hausman_formula, data = pdata, model = "within",
                   effect = "individual")

hausman_test <- tryCatch(phtest(model_fe_h, model_re_h), error = function(e) NULL)

if (!is.null(hausman_test)) {
  cat("Test statistic:", round(hausman_test$statistic, 4), "\n")
  cat("p-value:", format.pval(hausman_test$p.value, digits = 4), "\n")
  if (hausman_test$p.value < 0.05) {
    cat("Reject H0: FE is safer (RE assumption violated).\n")
  } else {
    cat("Cannot reject H0: RE is appropriate.\n")
  }
}

# ============================================================================
# PUBLICATION TABLE: Three-model comparison (stargazer + gt fallback)
# ============================================================================

cat("\n=", " PUBLICATION TABLE ", "=", "\n")

if (requireNamespace("stargazer", quietly = TRUE)) {

  cat("Generating publication tables with stargazer ...\n")

  sg_labels <- c(
    "Informal (not registered)",
    "Modified dep. ratio (OECD)",
    "Earner loss (t vs t-1)",
    "No. formal earners",
    "No. informal earners"
  )

  # Dynamically build labels for employment_type factors if present
  if ("employment_type" %in% usable_factors) {
    emp_ref <- levels(model_panel$employment_type)[1]
    emp_labs <- setdiff(levels(model_panel$employment_type), emp_ref)
    sg_labels <- c(sg_labels, emp_labs)
  }

  # Education, age, sex labels for RE / Mundlak
  if ("education_level" %in% usable_factors) {
    edu_ref <- levels(model_panel$education_level)[1]
    edu_labs <- setdiff(levels(model_panel$education_level), edu_ref)
    sg_labels <- c(sg_labels, edu_labs)
  }
  if ("age_group" %in% usable_factors) {
    age_ref <- levels(model_panel$age_group)[1]
    age_labs <- setdiff(levels(model_panel$age_group), age_ref)
    sg_labels <- c(sg_labels, age_labs)
  }
  sg_labels <- c(sg_labels, "Female")

  # Year dummies
  year_labs <- paste0("Year: ", sort(unique(model_panel$year))[-1])
  sg_labels <- c(sg_labels, year_labs)

  # Mundlak terms
  mundlak_labels <- c(
    "Avg. informal (Mundlak)",
    "Avg. dep. ratio OECD (Mundlak)",
    "Avg. employment (Mundlak)",
    "Avg. earner loss (Mundlak)",
    "Avg. formal earners (Mundlak)",
    "Avg. informal earners (Mundlak)"
  )
  sg_labels <- c(sg_labels, mundlak_labels)

  # HTML
  stargazer::stargazer(
    model_re, model_fe, model_mundlak,
    type = "html",
    title = "Panel Regression Models: Poverty Determinants (TR-SILC 2016-2019)",
    column.labels = c("Random Effects", "Fixed Effects", "Mundlak Correction"),
    dep.var.labels = "Poverty Status (1 = poor at 60% threshold)",
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
    notes = "Robust standard errors clustered at individual level (Arellano). *** p<0.01, ** p<0.05, * p<0.10.",
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
    dep.var.labels = "Poverty Status (1 = poor at 60\\% threshold)",
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

  # Plain text to console
  stargazer::stargazer(
    model_re, model_fe, model_mundlak,
    type = "text",
    title = "Panel Regression Models: Poverty Determinants",
    column.labels = c("Random Effects", "Fixed Effects", "Mundlak Correction"),
    dep.var.labels = "Poverty Status (1 = poor)",
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
    notes = "Robust SE clustered at individual level.",
    notes.align = "l",
    out = file.path(project$table_dir, "04_models_table.txt")
  )

  cat("Tables saved: 04_models_table.html, .tex, .txt\n")

} else {
  cat("stargazer not available. Saving manual coefficient table ...\n")

  tidy_model <- function(model, robust_vcov, model_name) {
    cf  <- coef(model)
    se  <- sqrt(diag(robust_vcov))
    z   <- cf / se
    p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
    tibble(
      model = model_name, term = names(cf),
      estimate = cf, std_error = se, z_value = z, p_value = p,
      stars = add_significance_stars(p),
      ci_low = cf - 1.96 * se, ci_high = cf + 1.96 * se
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
# KEY COEFFICIENT COMPARISON
# ============================================================================

cat("\n=", " COEFFICIENT COMPARISON ", "=", "\n")

make_comparison_row <- function(var_name, label) {
  extract <- function(model, vcov) {
    cf <- coef(model)
    if (var_name %in% names(cf)) {
      se <- sqrt(diag(vcov))[var_name]
      p  <- 2 * pnorm(abs(cf[var_name] / se), lower.tail = FALSE)
      return(list(est = cf[var_name], se = se, p = p))
    }
    list(est = NA, se = NA, p = NA)
  }

  re <- extract(model_re, re_robust_vcov)
  fe <- extract(model_fe, fe_robust_vcov)
  mk <- extract(model_mundlak, mundlak_robust_vcov)

  tibble(
    Variable   = label,
    RE_coef    = ifelse(is.na(re$est), "", sprintf("%.4f%s", re$est, add_significance_stars(re$p))),
    RE_se      = ifelse(is.na(re$se),  "", sprintf("(%.4f)", re$se)),
    FE_coef    = ifelse(is.na(fe$est), "", sprintf("%.4f%s", fe$est, add_significance_stars(fe$p))),
    FE_se      = ifelse(is.na(fe$se),  "", sprintf("(%.4f)", fe$se)),
    Mundlak_coef = ifelse(is.na(mk$est), "", sprintf("%.4f%s", mk$est, add_significance_stars(mk$p))),
    Mundlak_se   = ifelse(is.na(mk$se),  "", sprintf("(%.4f)", mk$se)),
    Diff_RE_FE = ifelse(is.na(re$est) | is.na(fe$est), "", sprintf("%.4f", re$est - fe$est))
  )
}

comparison_table <- bind_rows(
  make_comparison_row("informal_status",       "Informal (not registered)"),
  make_comparison_row("dependency_ratio_oecd", "Modified dep. ratio (OECD)"),
  make_comparison_row("dependency_ratio", "Modified dep. ratio"),
  make_comparison_row("female",                "Female")
)

cat("\n")
print(as.data.frame(comparison_table), row.names = FALSE)
write_csv(comparison_table, file.path(project$table_dir, "04_coefficient_comparison.csv"))

# ============================================================================
# INTERPRETATION
# ============================================================================

cat("\n--- Interpretation ---\n")

re_cf <- coef(model_re)
fe_cf <- coef(model_fe)

if ("informal_status" %in% names(re_cf)) {
  re_inf <- re_cf["informal_status"]
  fe_inf <- fe_cf["informal_status"]
  cat(sprintf("\nInformality:\n  RE: %.1f pp | FE: %.1f pp\n", re_inf * 100, fe_inf * 100))
  if (abs(fe_inf) < abs(re_inf)) {
    bias_pct <- (1 - fe_inf / re_inf) * 100
    cat(sprintf("  ~%.0f%% of the RE effect may reflect selection.\n", abs(bias_pct)))
  }
}

if ("dependency_ratio_oecd" %in% names(re_cf)) {
  re_dep <- re_cf["dependency_ratio_oecd"]
  fe_dep <- fe_cf["dependency_ratio_oecd"]
  cat(sprintf("\nModified dependency ratio:\n  RE: %.1f pp | FE: %.1f pp\n",
              re_dep * 100, fe_dep * 100))}

# ============================================================================
# DIAGNOSTIC SUMMARY
# ============================================================================
cat("\n--- Diagnostics ---\n")
cat("Rho (RE):", round(rho_re, 4), "\n")
cat("R² within (FE):", round(fe_r2_within, 4), "\n")
if (!is.null(hausman_test)) {
  cat("Hausman p-value:", format.pval(hausman_test$p.value, digits = 4), "\n")
}
cat("Observations — RE:", nobs(model_re), "| FE:", nobs(model_fe),
    "| Mundlak:", nobs(model_mundlak), "\n")

# ============================================================================
# SAVE MODEL OBJECTS
# ============================================================================
saveRDS(model_re,      file.path(project$model_dir, "04_model_re.rds"))
saveRDS(model_fe,      file.path(project$model_dir, "04_model_fe.rds"))
saveRDS(model_mundlak, file.path(project$model_dir, "04_model_mundlak.rds"))
cat("Model objects saved.\n")

pdata$pred_poverty_re <- fitted(model_re)
pred_data <- as.data.frame(pdata) %>% select(id, year, poverty_status, pred_poverty_re)
write_csv(pred_data, file.path(project$out_dir, "04_predictions_re.csv"))
cat("RE predictions saved.\n")

cat("\nPanel poverty regression analysis completed.\n")
cat("Models: ", normalizePath(project$model_dir), "\n")
cat("Tables: ", normalizePath(project$table_dir), "\n")

