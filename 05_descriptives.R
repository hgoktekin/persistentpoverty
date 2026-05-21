# =============================================================================
# Descriptive statistics for all variables used in the analysis
#
# This script produces two publication-ready descriptive tables:
#
#   Table A — Continuous variables:
#     N, Mean, SD, Min, Median, Max for age, household income,
#     equivalised income, household composition, and dependency ratios.
#
#   Table B — Categorical variables:
#     Unweighted count (N) and share (%) for sex, education, employment
#     status, social security registration, age group, and poverty status.
#
# Both tables are computed on the balanced panel (2016-2019, all person-years)
# and saved as LaTeX (.tex) files via the gt package.  The underlying data
# frames (table_a, table_b) remain attached in the R environment for
# interactive inspection.
#
# Outputs:  tables/table_descriptives_continuous.tex
#           tables/table_descriptives_categorical.tex
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(haven)
  library(readr)
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

source(file.path(project_root,  "00_config.R"))
source(file.path(project_root,  "01_functions.R"))

project$data_path  <- file.path(project_root, project$data_path)
project$out_dir    <- file.path(project_root, project$out_dir)
project$table_dir  <- file.path(project_root, project$table_dir)

dir.create(project$out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(project$table_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read data and build analysis panel -------------------------------------

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

# ---- Recode categorical variables for display ------------------------------

panel_desc <- panel_poverty %>%
  mutate(
    sex_recoded = factor(
      case_when(
        .data[[vars$sex]] == codes$sex[["female"]] ~ "Female",
        .data[[vars$sex]] == codes$sex[["male"]]   ~ "Male",
        TRUE ~ NA_character_
      ),
      levels = c("Female", "Male")
    ),
    age_group = cut(
      .data[[vars$age]],
      breaks = c(-Inf, 17, 24, 34, 44, 54, 64, Inf),
      labels = c("0-17", "18-24", "25-34", "35-44", "45-54", "55-64", "65+")
    ),
    education_recoded = factor(
      case_when(
        .data[[vars$education]] %in% c(0, 1, 2) ~ "Primary or below",
        .data[[vars$education]] %in% c(3, 4, 5) ~ "Secondary",
        .data[[vars$education]] %in% c(6, 7, 8) ~ "Tertiary",
        TRUE ~ NA_character_
      ),
      levels = c("Primary or below", "Secondary", "Tertiary")
    ),
    labour_recoded = factor(
      case_when(
        .data[[vars$labour_status]] %in% c(1, 2)          ~ "Employee",
        .data[[vars$labour_status]] %in% c(3, 4)          ~ "Self-employed",
        .data[[vars$labour_status]] == 7                   ~ "Retired",
        .data[[vars$labour_status]] %in% c(5, 6, 8, 9, 10) ~ "Inactive",
        TRUE ~ NA_character_
      ),
      levels = c("Employee", "Self-employed", "Retired", "Inactive")
    ),
    social_security_recoded = factor(
      if_else(informal_status == 1L, "Not registered", "Registered"),
      levels = c("Registered", "Not registered")
    ),
    poverty_status_60 = factor(
      if_else(poor_60 == 1L, "Poor", "Non-poor"),
      levels = c("Non-poor", "Poor")
    )
  )

# =============================================================================
# TABLE A — Continuous variables: N, Mean, SD, Min, Median, Max
# =============================================================================

cat("Building Table A (continuous variables) ...\n")

continuous_vars <- list(
  list(label = "Age",                               var = vars$age),
  list(label = "Household income (TL)",              var = vars$household_income),
  list(label = "Equivalised income (TL)",            var = "eq_income"),
  list(label = "Household size",                     var = "hh_size"),
  list(label = "No. children under 14",              var = "hh_children_u14"),
  list(label = "No. children under 18",              var = "hh_children_u18"),
  list(label = "No. elderly (65+)",                  var = "hh_elderly_65plus"),
  list(label = "No. earners (proxy)",                var = "hh_earners_proxy"),
  list(label = "No. formal earners",                 var = "hh_formal_earners"),
  list(label = "No. informal earners",               var = "hh_informal_earners"),
  list(label = "Dependency ratio",                   var = "dependency_ratio"),
  list(label = "Modified dep. ratio (OECD)",         var = "dependency_ratio_oecd")
)

fmt_num <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

table_a <- do.call(rbind, lapply(continuous_vars, function(v) {
  x <- panel_desc[[v$var]]
  x <- x[!is.na(x)]
  data.frame(
    Variable  = v$label,
    N         = formatC(length(x), format = "d", big.mark = ","),
    Mean      = fmt_num(mean(x)),
    SD        = fmt_num(sd(x)),
    Min       = fmt_num(min(x)),
    Median    = fmt_num(median(x)),
    Max       = fmt_num(max(x)),
    stringsAsFactors = FALSE
  )
}))
rownames(table_a) <- NULL

cat("\n--- Table A: Continuous variables ---\n")
print(table_a)

# =============================================================================
# TABLE B — Categorical variables: unweighted N and share (%)
# =============================================================================

cat("\nBuilding Table B (categorical variables) ...\n")

categorical_vars <- list(
  list(section = "Sex",                    var = "sex_recoded"),
  list(section = "Age group",              var = "age_group"),
  list(section = "Education",      var = "education_recoded"),
  list(section = "Employment status",      var = "labour_recoded"),
  list(section = "Social security",        var = "social_security_recoded"),
  list(section = "Poverty status (60%)",   var = "poverty_status_60")
)

table_b_rows <- list()
for (cv in categorical_vars) {
  x <- panel_desc[[cv$var]]
  x <- x[!is.na(x)]
  total <- length(x)
  tab <- table(x)
  # Section header row
  table_b_rows[[length(table_b_rows) + 1]] <- data.frame(
    Variable = cv$section,
    N        = "",
    Share    = "",
    row_type = "section",
    stringsAsFactors = FALSE
  )
  for (lvl in names(tab)) {
    n <- as.integer(tab[lvl])
    pct <- 100 * n / total
    table_b_rows[[length(table_b_rows) + 1]] <- data.frame(
      Variable = paste0("  ", lvl),
      N        = formatC(n, format = "d", big.mark = ","),
      Share    = sprintf("%.1f%%", pct),
      row_type = "category",
      stringsAsFactors = FALSE
    )
  }
}
table_b <- do.call(rbind, table_b_rows)
rownames(table_b) <- NULL

cat("\n--- Table B: Categorical variables ---\n")
print(table_b[, c("Variable", "N", "Share")])

# =============================================================================
# SAVE .tex TABLES (gt)
# =============================================================================

cat("\nGenerating LaTeX tables ...\n")

if (!requireNamespace("gt", quietly = TRUE)) {
  stop("The gt package is required to produce .tex output. Install it with install.packages('gt').")
}

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

# ---- Table A: continuous variables ------------------------------------------

table_a_gt <- table_a %>%
  gt::gt() %>%
  gt::tab_header(
    title    = gt::md("**Descriptive Statistics: Continuous Variables**"),
    subtitle = "Balanced panel, TR-SILC 2016--2019 (all person-years)"
  ) %>%
  gt::cols_label(
    Variable = "",
    N        = "N",
    Mean     = "Mean",
    SD       = "SD",
    Min      = "Min",
    Median   = "Median",
    Max      = "Max"
  ) %>%
  gt::cols_align(align = "left",  columns = Variable) %>%
  gt::cols_align(align = "right", columns = c(N, Mean, SD, Min, Median, Max)) %>%
  gt::tab_source_note(
    source_note = gt::md(
      "Unweighted statistics computed on the balanced panel (all 4 waves). Source: TR-SILC 2016--2019."
    )
  ) %>%
  gt_theme()

gt::gtsave(table_a_gt, file.path(project$table_dir, "table_descriptives_continuous.tex"))
cat("Saved:", file.path(project$table_dir, "table_descriptives_continuous.tex"), "\n")

# ---- Table B: categorical variables ----------------------------------------

table_b_gt <- table_b %>%
  gt::gt() %>%
  gt::cols_hide(columns = row_type) %>%
  gt::tab_header(
    title    = gt::md("**Descriptive Statistics: Categorical Variables**"),
    subtitle = "Balanced panel, TR-SILC 2016--2019 (all person-years)"
  ) %>%
  gt::cols_label(
    Variable = "",
    N        = "N (unweighted)",
    Share    = "Share (%)"
  ) %>%
  gt::cols_align(align = "left",  columns = Variable) %>%
  gt::cols_align(align = "right", columns = c(N, Share)) %>%
  gt::tab_style(
    style     = list(gt::cell_text(weight = "bold")),
    locations = gt::cells_body(rows = row_type == "section")
  ) %>%
  gt::tab_source_note(
    source_note = gt::md(
      "Unweighted counts and shares computed on the balanced panel (all 4 waves). Source: TR-SILC 2016--2019."
    )
  ) %>%
  gt_theme()

gt::gtsave(table_b_gt, file.path(project$table_dir, "table_descriptives_categorical.tex"))
cat("Saved:", file.path(project$table_dir, "table_descriptives_categorical.tex"), "\n")

# =============================================================================
# LATEX INPUT SNIPPET
# =============================================================================

desc_tex <- c(
  "% LaTeX input snippet for descriptive statistics tables.",
  "% Add these lines to your main .tex document.",
  "",
  "\\input{tables/table_descriptives_continuous.tex}",
  "\\input{tables/table_descriptives_categorical.tex}"
)
writeLines(desc_tex, file.path(project$table_dir, "05_descriptives_inputs.tex"))
cat("LaTeX input snippet saved:", file.path(project$table_dir, "05_descriptives_inputs.tex"), "\n")

cat("\n05_descriptives completed.\n")
cat("Objects in R environment: table_a (continuous), table_b (categorical)\n")
