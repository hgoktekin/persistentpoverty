# ============================================================
# 06_model.R  –  Dynamic CRE probit (Wooldridge 2005)
#
# Sources 00_config.R (paths, variable map) and 01_functions.R
# (data-pipeline helpers).  Builds the balanced analysis panel
# using the same pipeline as 02, then estimates three
# progressively augmented dynamic random-effects probit models:
#
#   Model 1 – base:  y_lag + y_init + X_it + H_it + Z_i
#   Model 2 – CRE:   + Mundlak means
#   Model 3 – full:  + wave dummies
#
# ★ To change the specification, edit the SPECIFICATION block.
# ============================================================

suppressPackageStartupMessages({
  library(glmmTMB)
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
  use_weights  = TRUE    # longitudinal survey weights (normalised)
)

# X_it  – time-varying covariates of the household head
x_tv <- c("age_young",         # 0/1  (age < 30; ref = 30-64)
           "age_elderly",        # 0/1  (age > 64; ref = 30-64)
           "informal_status",   # 0/1  (social-security registration)
           "labour_recoded")    # factor

# H_it  – time-varying household-level covariates
h_tv <- c("hh_size",                # integer
           "dependency_ratio_oecd",  # continuous (OECD-weighted)
           "other_earners",          # integer (earners excl. head)
           "has_child_u14")          # 0/1  (any child < 14 in HH)

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
    other_earners  = pmax(hh_earners_proxy - is_earner_proxy, 0L),
    has_child_u14  = as.integer(hh_n_children_u14 > 0),

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

# -- e) estimation sample: 2017-2019 (drop initial period) -------------------

est <- df %>%
  filter(year > first_year, !is.na(y_lag)) %>%
  droplevels()

if (spec$use_weights)
  est <- est %>% mutate(wt = wt / mean(wt))

cat(sprintf("  Sample: %d obs, %d individuals, T = %s\n",
            nrow(est), n_distinct(est$pid),
            paste(sort(unique(est$year)), collapse = ", ")))

# ============================================================
# ESTIMATION  –  three progressive specifications
# ============================================================
# Model 1 (base):  y_lag + y_init + X + H + Z
# Model 2 (CRE):   + Mundlak means
# Model 3 (full):  + wave dummies

base_rhs <- c("y_lag", "y_init", x_tv, h_tv, z_ti)

model_specs <- list(
  list(name = "Model 1 (base)",
       rhs  = base_rhs),
  list(name = "Model 2 (+Mundlak)",
       rhs  = c(base_rhs, mundlak_names)),
  list(name = "Model 3 (+wave FE)",
       rhs  = c(base_rhs, mundlak_names, "factor(year)"))
)

fits      <- list()
coef_tbls <- list()

for (m in seq_along(model_specs)) {
  sp  <- model_specs[[m]]
  frm <- as.formula(paste("y ~", paste(sp$rhs, collapse = " + "),
                          "+ (1 | pid)"))

  cat(sprintf("\n========== %s ==========\n", sp$name))
  cat("Formula:\n "); cat(deparse(frm, width.cutoff = 300), "\n\n")

  # collinearity check
  X_check <- model.matrix(update(frm, . ~ . - (1 | pid)), data = est)
  qr_rank <- qr(X_check)$rank
  if (qr_rank < ncol(X_check))
    cat(sprintf("  WARNING: rank = %d < %d columns (near-collinearity)\n\n",
                qr_rank, ncol(X_check)))

  cat("Estimating ...\n")
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

  # coefficients and APE
  fe <- fixef(fit)$cond
  se <- sqrt(diag(vcov(fit)$cond))
  z  <- fe / se
  p  <- 2 * pnorm(abs(z), lower.tail = FALSE)

  Xb      <- predict(fit, type = "link", re.form = NA)
  avg_phi <- mean(dnorm(Xb))

  coef_tbl <- tibble(
    term      = names(fe),
    estimate  = round(fe, 5),
    std_error = round(se, 5),
    z_value   = round(z, 3),
    p_value   = round(p, 4),
    stars     = add_significance_stars(p),
    ape       = round(avg_phi * fe, 5)
  )

  cat("\n--- Coefficients and APEs ---\n")
  print(coef_tbl, n = Inf)

  # variance components
  vc      <- as.data.frame(VarCorr(fit))
  sigma_a <- vc$sdcor[vc$grp == "pid" & vc$component == "cond"]
  rho_re  <- sigma_a^2 / (sigma_a^2 + 1)
  cat(sprintf("\nsigma_a = %.4f,  rho = %.4f\n", sigma_a, rho_re))

  fits[[m]]      <- fit
  coef_tbls[[m]] <- coef_tbl %>%
    mutate(sigma_a = sigma_a, rho_re = rho_re,
           n_obs = nrow(est), n_pid = n_distinct(est$pid))
}

# ============================================================
# SAVE
# ============================================================

# --- CSV (one per model) -----------------------------------------------------
for (m in seq_along(model_specs)) {
  fname <- sprintf("06_cre_probit_m%d_coefficients.csv", m)
  write_csv(coef_tbls[[m]], file.path(project$table_dir, fname))
}
for (m in seq_along(fits)) {
  fname <- sprintf("06_cre_probit_m%d.rds", m)
  saveRDS(fits[[m]], file.path(project$model_dir, fname))
}

# --- LaTeX table (all 3 models side by side) ---------------------------------

# collect all unique terms across models, in the order they first appear
all_terms <- unique(unlist(lapply(coef_tbls, \(ct) ct$term)))

tex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Dynamic Random-Effects Probit: Progressive Specifications}",
  "\\label{tab:cre-probit}",
  "\\begin{tabular}{lccc}",
  "\\hline\\hline",
  "& Model 1 & Model 2 & Model 3 \\\\",
  "& (Base) & (+Mundlak) & (+Wave FE) \\\\",
  "\\hline"
)

for (trm in all_terms) {
  lab <- gsub("_", "\\\\_", trm)
  cells <- character(3)
  se_cells <- character(3)
  for (m in 1:3) {
    ct  <- coef_tbls[[m]]
    row <- ct[ct$term == trm, ]
    if (nrow(row) == 1) {
      cells[m]    <- sprintf("%.4f%s", row$estimate, row$stars)
      se_cells[m] <- sprintf("(%.4f)", row$std_error)
    } else {
      cells[m]    <- ""
      se_cells[m] <- ""
    }
  }
  tex_lines <- c(tex_lines,
    sprintf("%s & %s & %s & %s \\\\", lab, cells[1], cells[2], cells[3]),
    sprintf(" & %s & %s & %s \\\\", se_cells[1], se_cells[2], se_cells[3])
  )
}

# footer: sigma_a, rho, N
sa <- sapply(coef_tbls, \(ct) ct$sigma_a[1])
rr <- sapply(coef_tbls, \(ct) ct$rho_re[1])
nn <- sapply(coef_tbls, \(ct) ct$n_obs[1])
np <- sapply(coef_tbls, \(ct) ct$n_pid[1])

tex_lines <- c(tex_lines,
  "\\hline",
  sprintf("$\\sigma_a$ & %.4f & %.4f & %.4f \\\\", sa[1], sa[2], sa[3]),
  sprintf("$\\rho$ & %.4f & %.4f & %.4f \\\\", rr[1], rr[2], rr[3]),
  sprintf("Observations & %s & %s & %s \\\\",
          formatC(nn[1], big.mark = ","),
          formatC(nn[2], big.mark = ","),
          formatC(nn[3], big.mark = ",")),
  sprintf("Individuals & %s & %s & %s \\\\",
          formatC(np[1], big.mark = ","),
          formatC(np[2], big.mark = ","),
          formatC(np[3], big.mark = ",")),
  "Mundlak means & No & Yes & Yes \\\\",
  "Wave dummies & No & No & Yes \\\\",
  "\\hline\\hline",
  "\\multicolumn{4}{l}{\\footnotesize *** p$<$0.01, ** p$<$0.05, * p$<$0.10.} \\\\",
  "\\multicolumn{4}{l}{\\footnotesize Standard errors in parentheses.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(tex_lines, file.path(project$table_dir, "06_cre_probit.tex"))

cat("\nOutputs:\n")
for (m in 1:3) {
  cat(sprintf("  %s\n", file.path(project$table_dir,
              sprintf("06_cre_probit_m%d_coefficients.csv", m))))
  cat(sprintf("  %s\n", file.path(project$model_dir,
              sprintf("06_cre_probit_m%d.rds", m))))
}
cat(sprintf("  %s\n", file.path(project$table_dir, "06_cre_probit.tex")))
cat("\nDone.\n")
