# ============================================================
# 06_model.R  –  Dynamic CRE probit (Wooldridge 2005)
#
# Sources 00_config.R (paths, variable map) and 01_functions.R
# (data-pipeline helpers).  Builds the balanced analysis panel
# using the same pipeline as 02, then estimates a dynamic
# correlated random-effects probit.
#
# ★ To test a different specification, edit the SPECIFICATION
#   block below.  Everything downstream adjusts automatically.
# ============================================================

suppressPackageStartupMessages({
  library(glmmTMB)   # RE probit — more robust than lme4::glmer for
                     # near-singular designs and survey weights
})

# --- Project setup -----------------------------------------------------------

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f    <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(normalizePath(sub("^--file=", "", f[1]), mustWork = TRUE))
  s <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
                error = function(e) NA_character_)
  if (!is.na(s)) return(s)
  NA_character_
}
script_path  <- get_script_path()
project_root <- if (is.na(script_path)) getwd() else
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(project_root, "00_config.R"))
source(file.path(project_root, "01_functions.R"))

project$data_path <- file.path(project_root, project$data_path)
project$out_dir   <- file.path(project_root, project$out_dir)
project$table_dir <- file.path(project_root, project$table_dir)
project$model_dir <- file.path(project_root, project$model_dir)
for (d in c(project$out_dir, project$table_dir, project$model_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- Data pipeline (functions from 01, config from 00) -----------------------

cat("Step 1 – reading raw SILC panel data\n")
raw <- read_dta(project$data_path)

cat("Step 2 – balanced panel, equivalised income, household context\n")
panel_balanced <- construct_balanced_panel(raw, vars,
                                           project$panel_years,
                                           project$reference_year)
panel_income <- panel_balanced %>%
  add_equivalised_income(vars) %>%
  add_household_context(vars, codes) %>%
  add_employment_stability(vars)

poverty_lines <- compute_poverty_lines(panel_income, vars, project$thresholds)
panel_poverty <- add_poverty_status(panel_income, poverty_lines, vars,
                                    project$thresholds)

# ============================================================
# ★  SPECIFICATION  –  change these to try different models
# ============================================================

spec <- list(
  threshold    = 60,     # poverty line (50, 60, or 70)
  heads_only   = TRUE,   # restrict sample to household reference persons
  use_mundlak  = TRUE,   # Mundlak (CRE) means of time-varying regressors
  use_initial  = TRUE,   # Wooldridge initial-condition term y_{i0}
  use_year_fe  = TRUE,   # year dummies
  use_weights  = TRUE,   # longitudinal survey weights (normalised)
  nAGQ         = 1       # quadrature points (1 = Laplace; >1 not supported by glmmTMB)
)

# X_it  – time-varying covariates of the household head
x_tv <- c("age_young",         # 0/1  (age < 30; ref = 30-64)
           "age_elderly",        # 0/1  (age > 64; ref = 30-64)
           "informal_status",   # 0/1  (social-security registration)
           "labour_recoded")    # factor

# H_it  – time-varying household-level covariates
h_tv <- c("hh_size",                # integer
           "dependency_ratio_oecd",  # continuous (OECD-weighted)
           "other_earners")          # integer (earners excl. head)

# Z_i   – time-invariant covariates
z_ti <- c("female",              # 0/1
           "education_recoded")  # factor

# ============================================================
# DATA PREPARATION
# ============================================================

cat("Step 3 – preparing estimation data\n")

poor_var   <- paste0("poor_", spec$threshold)
first_year <- min(project$panel_years)

# -- a) restrict to household heads (if requested) and recode -----------------

df <- panel_poverty
if (spec$heads_only)
  df <- df %>% filter(.data[[vars$household_reference_person]] == 1)

df <- df %>%
  mutate(
    pid  = .data[[vars$person_id]],
    year = .data[[vars$year]],
    wt   = .data[[vars$longitudinal_weight]],
    y    = as.integer(.data[[poor_var]]),

    # X_it  (Andriopoulou 2011: ref = aged 30-64)
    age_young   = as.integer(.data[[vars$age]] < 30),
    age_elderly = as.integer(.data[[vars$age]] > 64),
    labour_recoded = factor(
      case_when(
        .data[[vars$labour_status]] %in% 1:2         ~ "Employee",
        .data[[vars$labour_status]] %in% 3:4         ~ "Self-employed",
        .data[[vars$labour_status]] == 5             ~ "Unemployed",
        .data[[vars$labour_status]] == 7             ~ "Retired",
        .data[[vars$labour_status]] %in% c(6, 8:10) ~ "Inactive"),
      levels = c("Employee", "Self-employed", "Unemployed",
                 "Retired", "Inactive")),

    # H_it
    other_earners = pmax(hh_earners_proxy - is_earner_proxy, 0L),

    # Z_i
    education_recoded = factor(
      case_when(
        .data[[vars$education]] %in% 0:2 ~ "Primary or below",
        .data[[vars$education]] %in% 3:5 ~ "Secondary",
        .data[[vars$education]] == 6     ~ "Tertiary"),
      levels = c("Primary or below", "Secondary", "Tertiary"))
  )

# -- b) keep only individuals with complete data in every year ----------------

model_vars <- c("y", x_tv, h_tv, z_ti)
df <- df %>%
  filter(if_all(all_of(model_vars), ~ !is.na(.))) %>%
  group_by(pid) %>%
  filter(n_distinct(year) == length(project$panel_years)) %>%
  ungroup()

# ============================================================
# DIAGNOSTIC: time-varying vs time-invariant classification
# ============================================================
# For each candidate variable, compute the share of individuals
# whose value changes across the panel.  A variable is
# "time-invariant" if it never (or almost never) changes within
# person; "time-varying" if it does.

candidate_vars <- intersect(c("y", x_tv, h_tv, z_ti), names(df))

variation_check <- df %>%
  group_by(pid) %>%
  summarise(
    across(all_of(candidate_vars),
           ~ n_distinct(.) > 1,
           .names = "{.col}"),
    .groups = "drop")

variation_summary <- tibble(
  variable = candidate_vars,
  n_individuals     = nrow(variation_check),
  n_who_change      = sapply(candidate_vars,
                             function(v) sum(variation_check[[v]])),
  pct_who_change    = round(100 * n_who_change / n_individuals, 1),
  classification    = if_else(pct_who_change > 5,
                              "TIME-VARYING", "TIME-INVARIANT")
)

cat("\n--- Within-individual variation diagnostic ---\n")
cat("  (% of individuals whose value changes across 2016-2019)\n\n")
print(as.data.frame(variation_summary), row.names = FALSE)
cat("\n")

# -- c) lagged poverty y_{it-1} and initial condition y_{i0} -----------------

df <- df %>%
  arrange(pid, year) %>%
  group_by(pid) %>%
  mutate(y_lag = lag(y, order_by = year)) %>%
  ungroup()

y0 <- df %>%
  filter(year == first_year) %>%
  transmute(pid, y_init = y)
df <- left_join(df, y0, by = "pid")

# -- d) Mundlak means (computed over the full 2016-2019 window) --------------
#    Factors are expanded to K-1 treatment-contrast dummies before averaging.

mundlak_names <- character(0)
if (spec$use_mundlak) {
  tv_all  <- c(x_tv, h_tv)
  tv_mat  <- model.matrix(reformulate(tv_all), data = df)[, -1, drop = FALSE]
  colnames(tv_mat) <- make.names(colnames(tv_mat), unique = TRUE)

  mundlak_df <- as_tibble(tv_mat) %>%
    mutate(pid = df$pid) %>%
    group_by(pid) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") %>%
    rename_with(\(nm) paste0("m_", nm), .cols = -pid)

  df <- left_join(df, mundlak_df, by = "pid")
  mundlak_names <- setdiff(names(mundlak_df), "pid")
}

# -- e) estimation sample: 2017-2019 (drop initial period) -------------------

est <- df %>%
  filter(year > first_year, !is.na(y_lag)) %>%
  droplevels()

if (spec$use_weights)
  est <- est %>% mutate(wt = wt / mean(wt))

# -- f) scale continuous covariates for numerical stability --------------------
# Keeps binary / factor columns untouched; only centres & scales numeric
# variables with SD > 0.  Scaling improves optimizer conditioning and avoids
# false-convergence warnings.

continuous_vars <- c("hh_size", "dependency_ratio_oecd", "other_earners")
scale_info <- list()  # store centres & scales if you need to back-transform
for (v in continuous_vars) {
  if (v %in% names(est) && is.numeric(est[[v]])) {
    mu <- mean(est[[v]], na.rm = TRUE)
    sd <- sd(est[[v]], na.rm = TRUE)
    if (sd > 0) {
      est[[v]] <- (est[[v]] - mu) / sd
      scale_info[[v]] <- list(center = mu, scale = sd)
      # also scale the corresponding Mundlak mean (same transformation)
      m_v <- paste0("m_", v)
      if (m_v %in% names(est))
        est[[m_v]] <- (est[[m_v]] - mu) / sd
    }
  }
}
if (length(scale_info) > 0)
  cat(sprintf("  Scaled %d continuous covariate(s): %s\n",
              length(scale_info), paste(names(scale_info), collapse = ", ")))

cat(sprintf("  Sample: %d obs, %d individuals, T = %s\n",
            nrow(est), n_distinct(est$pid),
            paste(sort(unique(est$year)), collapse = ", ")))

# ============================================================
# BUILD FORMULA
# ============================================================

rhs <- c(
  "y_lag",                                      # state dependence (rho)
  x_tv, h_tv,                                   # X_it, H_it
  z_ti,                                          # Z_i
  if (spec$use_initial) "y_init",                # initial condition
  if (spec$use_mundlak) mundlak_names,           # Mundlak means
  if (spec$use_year_fe) "factor(year)"           # year dummies
)
frm <- as.formula(paste("y ~", paste(rhs, collapse = " + "), "+ (1 | pid)"))
cat("Formula:\n "); cat(deparse(frm, width.cutoff = 300), "\n\n")

# -- collinearity check -------------------------------------------------------
X_check <- model.matrix(update(frm, . ~ . - (1 | pid)), data = est)
qr_rank <- qr(X_check)$rank
if (qr_rank < ncol(X_check)) {
  cat(sprintf("  WARNING: design matrix rank = %d but has %d columns.\n",
              qr_rank, ncol(X_check)))
  cat("  Near-collinear terms may cause estimation problems.\n")
  cat("  Consider dropping redundant Mundlak means or covariates.\n\n")
}

# ============================================================
# ESTIMATION
# ============================================================
# Random-effects probit via glmmTMB (Laplace approximation).
# Survey weights are normalised to mean 1 and passed as prior
# weights; this gives consistent point estimates but standard
# errors may be conservative.  For strict design-based inference
# consider bootstrap or pseudo-likelihood alternatives.

cat("Estimating dynamic CRE probit ...\n")
fit <- glmmTMB(
  frm,
  data    = est,
  family  = binomial(link = "probit"),
  weights = if (spec$use_weights) est$wt,
  REML    = FALSE,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS"),
    optCtrl   = list(maxit = 1e4)
  )
)

cat("\n"); print(summary(fit))

# ============================================================
# COEFFICIENT TABLE & AVERAGE PARTIAL EFFECTS
# ============================================================

fe <- fixef(fit)$cond
se <- sqrt(diag(vcov(fit)$cond))
z  <- fe / se
p  <- 2 * pnorm(abs(z), lower.tail = FALSE)

coef_tbl <- tibble(
  term      = names(fe),
  estimate  = round(fe, 5),
  std_error = round(se, 5),
  z_value   = round(z, 3),
  p_value   = round(p, 4),
  stars     = add_significance_stars(p)
)

# APE = mean[phi(X*beta)] * beta  (continuous / binary regressors)
Xb      <- predict(fit, type = "link", re.form = NA)
avg_phi <- mean(dnorm(Xb))
coef_tbl <- coef_tbl %>%
  mutate(ape = round(avg_phi * estimate, 5))

cat("\n--- Coefficients and average partial effects ---\n")
print(coef_tbl, n = Inf)

# Variance components
vc      <- as.data.frame(VarCorr(fit))
sigma_a <- vc$sdcor[vc$grp == "pid" & vc$component == "cond"]
rho     <- sigma_a^2 / (sigma_a^2 + 1)
cat(sprintf("\nsigma_a = %.4f,  rho = sigma_a^2/(sigma_a^2+1) = %.4f\n",
            sigma_a, rho))

if ("y_lag" %in% names(fe))
  cat(sprintf("State dependence: rho_hat = %.4f  (p %s)\n",
              fe["y_lag"],
              ifelse(p["y_lag"] < 0.001, "< 0.001",
                     sprintf("= %.4f", p["y_lag"]))))

# ============================================================
# SAVE
# ============================================================

write_csv(coef_tbl, file.path(project$table_dir, "06_cre_probit_coefficients.csv"))
saveRDS(fit,        file.path(project$model_dir, "06_cre_probit.rds"))

# --- LaTeX table -------------------------------------------------------------

tex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Dynamic Correlated Random-Effects Probit (Wooldridge 2005)}",
  "\\label{tab:cre-probit}",
  sprintf("\\begin{tabular}{l%s}", paste(rep("r", 4), collapse = "")),
  "\\hline\\hline",
  "& Estimate & Std.~Error & APE & \\\\",
  "\\hline"
)

for (i in seq_len(nrow(coef_tbl))) {
  row <- coef_tbl[i, ]
  lab <- gsub("_", "\\\\_", row$term)
  tex_lines <- c(tex_lines, sprintf(
    "%s & %.4f%s & (%.4f) & %.4f \\\\",
    lab, row$estimate, row$stars, row$std_error, row$ape))
}

tex_lines <- c(tex_lines,
  "\\hline",
  sprintf("$\\sigma_a$ & \\multicolumn{3}{l}{%.4f} \\\\", sigma_a),
  sprintf("$\\rho$ & \\multicolumn{3}{l}{%.4f} \\\\", rho),
  sprintf("Observations & \\multicolumn{3}{l}{%s} \\\\",
          formatC(nrow(est), format = "d", big.mark = ",")),
  sprintf("Individuals & \\multicolumn{3}{l}{%s} \\\\",
          formatC(n_distinct(est$pid), format = "d", big.mark = ",")),
  "\\hline\\hline",
  "\\multicolumn{4}{l}{\\footnotesize *** p<0.01, ** p<0.05, * p<0.10.} \\\\",
  "\\multicolumn{4}{l}{\\footnotesize Cluster-robust SE at individual level.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(tex_lines, file.path(project$table_dir, "06_cre_probit.tex"))

cat("\nOutputs:\n",
    " ", file.path(project$table_dir, "06_cre_probit_coefficients.csv"), "\n",
    " ", file.path(project$table_dir, "06_cre_probit.tex"), "\n",
    " ", file.path(project$model_dir, "06_cre_probit.rds"), "\n",
    "\nDone.\n")
