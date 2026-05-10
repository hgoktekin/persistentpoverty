# =============================================================================
# Master script: poverty dynamics in Türkiye using SILC longitudinal micro-data
#
# Outputs:
# - balanced panel with equivalised income and annual poverty status
# - poverty lines and FGT indices
# - poverty spell typologies
# - duration and transition tables
# - socio-demographic profile table
# - probit model with robust standard errors and marginal effects
# - publication-quality poverty trend figure
#
# Attrition note:
# The analysis uses the balanced 2016-2019 longitudinal panel. This reduces
# attrition-related missingness in transition estimates because every included
# individual is observed in all four waves. It does not remove all selection
# concerns: people who leave the panel, migrate, die, or otherwise become
# unobserved may differ systematically from stayers. Recommended robustness
# checks include comparing baseline characteristics of stayers and attritors
# in the original rotational sample, repeating descriptive estimates with
# cross-sectional weights where appropriate, and checking whether results are
# sensitive to alternative poverty thresholds and spell definitions.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(readr)
  library(tidyr)
})

# Locate the project root robustly. This lets the script run either from the
# project directory with `Rscript R/02_run_silc_poverty_analysis.R` or from any
# other working directory with an absolute script path.
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }

  sourced_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
  if (!is.na(sourced_file)) {
    return(sourced_file)
  }

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

project$data_path <- file.path(project_root, project$data_path)
project$out_dir <- file.path(project_root, project$out_dir)
project$table_dir <- file.path(project_root, project$table_dir)
project$figure_dir <- file.path(project_root, project$figure_dir)
project$model_dir <- file.path(project_root, project$model_dir)

dir.create(project$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(project$table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(project$figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(project$model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project$out_dir, "r_cache"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project$out_dir, "r_cache", "sass"), showWarnings = FALSE, recursive = TRUE)
Sys.setenv(
  R_USER_CACHE_DIR = normalizePath(file.path(project$out_dir, "r_cache")),
  XDG_CACHE_HOME = normalizePath(file.path(project$out_dir, "r_cache")),
  XDG_CONFIG_HOME = normalizePath(file.path(project$out_dir, "r_cache")),
  SASS_CACHE_PATH = normalizePath(file.path(project$out_dir, "r_cache", "sass"))
)

cat("Step 1: reading raw SILC panel data\n")
raw <- read_dta(project$data_path)

# 02 dosyasında Step 1 kısmını buna çevir:
panel_audit <- calculate_attrition_rate(raw, vars, project$panel_years, project$reference_year)

# Artık bunlar çalışacaktır:
cat("Attrition rate:", panel_audit$attrition_rate_percent, "%\n")
write_csv(panel_audit$baseline_comparison, file.path(project$out_dir, "attrition_comparison.csv"))

cat("Step 2: constructing balanced panel and propagating longitudinal weights\n")
panel_balanced <- construct_balanced_panel(
  raw = raw,
  vars = vars,
  panel_years = project$panel_years,
  reference_year = project$reference_year
)

panel_audit <- panel_balanced %>%
  summarise(
    individuals = n_distinct(.data[[vars$person_id]]),
    person_years = n(),
    first_year = min(.data[[vars$year]]),
    last_year = max(.data[[vars$year]]),
    duplicated_person_years = sum(duplicated(paste(.data[[vars$person_id]], .data[[vars$year]]))),
    missing_weights = sum(is.na(.data[[vars$longitudinal_weight]]))
  )

cat("Step 3: calculating modified OECD equivalence scale and equivalised income\n")
panel_income <- panel_balanced %>%
  add_equivalised_income(vars) %>%
  add_household_context(vars, codes)

cat("Step 4: computing annual within-sample poverty thresholds\n")
poverty_lines <- compute_poverty_lines(panel_income, vars, project$thresholds)

panel_poverty <- add_poverty_status(
  panel = panel_income,
  poverty_lines = poverty_lines,
  vars = vars,
  thresholds = project$thresholds)

cat("Step 5: classifying poverty spells and Eurostat persistent poverty\n")
classified_main <- classify_poverty_spells(
  panel = panel_poverty,
  vars = vars,
  panel_years = project$panel_years,
  reference_year = project$reference_year,
  threshold = project$main_threshold)

classified_all <- purrr::map_dfr(project$thresholds, function(threshold) {
  classify_poverty_spells(
    panel = panel_poverty,
    vars = vars,
    panel_years = project$panel_years,
    reference_year = project$reference_year,
    threshold = threshold) %>%
    select(all_of(vars$person_id), all_of(vars$longitudinal_weight), threshold,
           n_poor_4yr, n_poor_previous3, current_poor, persistent_poor, poverty_group)
})

poverty_definitions <- tibble::tribble(
  ~category, ~logical_definition,
  "Never poor", "poor in zero of the four observed years",
  "Transient poor", "poor in exactly one of the four observed years",
  "Frequently poor", "poor in at least two years, but not persistent by the Eurostat current-year rule",
  "Current poor (2019)", "equivalised income in 2019 is below the selected poverty threshold",
  "Persistent poor", "current poor in 2019 and poor in at least two of 2016, 2017, and 2018"
)

cat("Step 6: creating descriptive poverty tables and FGT indices\n")
table1_poverty_rates <- make_table1_poverty_rates(panel_poverty, vars, project$thresholds)
table2_group_distribution <- make_poverty_group_distribution(classified_main, vars)

duration_outputs <- make_poverty_duration_table(
  panel = panel_poverty,
  classified = classified_main,
  vars = vars,
  panel_years = project$panel_years,
  threshold = project$main_threshold)

table3_duration <- duration_outputs$duration
table3_episode_summary <- duration_outputs$episodes %>%
  summarise(
    mean_episodes = weighted_mean_safe(n_episodes, panel_weight),
    mean_max_spell_duration = weighted_mean_safe(max_spell_duration, panel_weight),
    mean_spell_duration = weighted_mean_safe(mean_spell_duration, panel_weight)
  )

table4_profile <- make_profile_table(
  panel = panel_poverty,
  classified = classified_main,
  vars = vars,
  codes = codes,
  reference_year = project$reference_year,
  threshold = project$main_threshold)

# 1. Formal çalışanlar için matris
table5_formal <- make_transition_matrices(
  panel = panel_poverty,
  vars = vars,
  panel_years = project$panel_years,
  threshold = project$main_threshold,
  filter_expr = "informal_status == 0" 
) %>% mutate(group = "Formal")

# 2. Informal çalışanlar için matris
table5_informal <- make_transition_matrices(
  panel = panel_poverty,
  vars = vars,
  panel_years = project$panel_years,
  threshold = project$main_threshold,
  filter_expr = "informal_status == 1"
) %>% mutate(group = "Informal")

# 3. İkisini birleştir
table5_comparison <- bind_rows(table5_formal, table5_informal)

# İstatistiksel Test (Chi-Square): 
# Formal ve informal gruplar arasında yoksulluktan çıkış (Exit) veya giriş (Entry) 
# oranları farklı mı?
test_data <- table5_comparison %>%
  filter(transition_type == "Exit") # Örnek: Yoksulluktan çıkış oranları farkı

# Not: Geçiş oranları arasındaki farkın testi genellikle 'prop.test' veya 
# survey paketindeki 'svychisq' ile panel verisi üzerinden yapılır.

table5_transitions <- make_transition_matrices(
  panel = panel_poverty,
  vars = vars,
  panel_years = project$panel_years,
  threshold = project$main_threshold)

mobility_summary <- table5_transitions %>%
  filter(transition_type %in% c("Entry", "Exit", "Poverty persistence")) %>%
  select(transition, transition_type, from_status, to_status, weighted_n, row_probability)

cat("Step 8: exporting tables, model outputs, and figure\n")
write_csv(panel_audit, file.path(project$out_dir, "panel_audit.csv"))
write_csv(poverty_lines, file.path(project$out_dir, "poverty_lines.csv"))
write_csv(panel_poverty, file.path(project$out_dir, "analysis_panel_2016_2019.csv"))
write_csv(classified_main, file.path(project$out_dir, "poverty_typology_60_2019_anchor.csv"))
write_csv(classified_all, file.path(project$out_dir, "poverty_typology_all_thresholds.csv"))
write_csv(poverty_definitions, file.path(project$table_dir, "poverty_group_definitions.csv"))
write_csv(table1_poverty_rates, file.path(project$table_dir, "table1_poverty_rates_fgt.csv"))
write_csv(table2_group_distribution, file.path(project$table_dir, "table2_poverty_group_distribution.csv"))
write_csv(table3_duration, file.path(project$table_dir, "table3_poverty_duration.csv"))
write_csv(table3_episode_summary, file.path(project$table_dir, "table3_episode_summary.csv"))
write_csv(table4_profile, file.path(project$table_dir, "table4_sociodemographic_profile.csv"))
write_csv(table5_transitions, file.path(project$table_dir, "table5_transition_matrices.csv"))
write_csv(mobility_summary, file.path(project$table_dir, "transition_mobility_summary.csv"))

if (requireNamespace("gt", quietly = TRUE)) {
  gt_theme <- function(x) {
    x %>%
      gt::tab_options(
        table.font.names = "Times New Roman",
        table.font.size = gt::px(12),
        table.border.top.color = "black",
        table.border.top.width = gt::px(1.5),
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

  table1_gt <- table1_poverty_rates %>%
      mutate(threshold = paste0(threshold, "%")) %>%
      gt::gt(groupname_col = "threshold") %>%
      gt::tab_header(
        title = gt::md("**Table 1. Poverty rates and FGT indices**"),
        subtitle = "Weighted estimates, Türkiye SILC balanced panel, 2016-2019"
      ) %>%
      gt::fmt_percent(
        columns = c(poverty_rate, poverty_rate_ci_low, poverty_rate_ci_high,
                    FGT0_headcount, FGT1_gap, FGT2_severity),
        decimals = 1
      ) %>%
      gt::fmt_number(columns = c(poverty_line, weighted_n), decimals = 0, use_seps = TRUE) %>%
      gt_theme()
  save_gt(table1_gt, "table1_poverty_rates_fgt")

  table2_gt <- table2_group_distribution %>%
      gt::gt(groupname_col = "table_family") %>%
      gt::tab_header(
        title = gt::md("**Table 2. Poverty group distribution**"),
        subtitle = "Persistent poor is shown both in the mutually exclusive typology and as a subset of current poverty."
      ) %>%
      gt::fmt_number(columns = weighted_n, decimals = 0, use_seps = TRUE) %>%
      gt::fmt_percent(columns = population_share, decimals = 1) %>%
      gt_theme()
  save_gt(table2_gt, "table2_poverty_group_distribution")

  table3_gt <- table3_duration %>%
      gt::gt() %>%
      gt::tab_header(
        title = gt::md("**Table 3. Poverty duration**"),
        subtitle = "Number of years below the 60% poverty threshold."
      ) %>%
      gt::fmt_number(columns = weighted_n, decimals = 0, use_seps = TRUE) %>%
      gt::fmt_percent(columns = population_share, decimals = 1) %>%
      gt_theme()
  save_gt(table3_gt, "table3_poverty_duration")

  table4_headers <- attr(table4_profile, "headers")
  table4_gt <- table4_profile %>%
    gt::gt() %>%
    gt::cols_hide(columns = row_type) %>%
    gt::cols_label(
      variable = "",
      p_value = gt::md("p-value^1^"),
      overall = gt::md(table4_headers$label[table4_headers$col == "overall"]),
      never_poor = gt::md(table4_headers$label[table4_headers$col == "never_poor"]),
      transient_poor = gt::md(table4_headers$label[table4_headers$col == "transient_poor"]),
      currently_poor = gt::md(table4_headers$label[table4_headers$col == "currently_poor"]),
      persistent_poor = gt::md(table4_headers$label[table4_headers$col == "persistent_poor"])
    ) %>%
    gt::tab_header(
      title = gt::md("**Table 4. Socio-demographic profile by poverty type, 2019**"),
      subtitle = gt::md(
        "*Note:* Weighted column percentages shown as N (%). p-values: weighted chi-squared (categorical); weighted Kruskal-Wallis rank test (continuous). Currently poor (2019) = poor in 2019 regardless of spell history. Persistent poor is a strict subset of currently poor (2019). Source: TR-SILC 2016--2019."
      )
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        "^1^ N (%) = unweighted count (weighted column percentage). Mean (SD) = weighted mean (weighted standard deviation). For this profile table, transient poor means poor in one or two years but not persistently poor. P-values are computed over the mutually exclusive four-year profile typology: never poor, transient poor, frequently poor, persistent poor."
      )
    ) %>%
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold")),
      locations = gt::cells_body(rows = row_type == "section")
    ) %>%
    gt::tab_style(
      style = gt::cell_borders(sides = "left", color = "#808080", weight = gt::px(1)),
      locations = list(
        gt::cells_body(columns = p_value),
        gt::cells_column_labels(columns = p_value)
      )
    ) %>%
    gt::tab_style(
      style = gt::cell_borders(sides = "left", color = "#d9d9d9", weight = gt::px(1)),
      locations = list(
        gt::cells_body(columns = currently_poor),
        gt::cells_column_labels(columns = currently_poor)
      )
    ) %>%
    gt::cols_align(align = "left", columns = everything()) %>%
    gt::cols_width(
      variable ~ gt::px(190),
      p_value ~ gt::px(90),
      overall ~ gt::px(145),
      never_poor ~ gt::px(145),
      transient_poor ~ gt::px(145),
      currently_poor ~ gt::px(145),
      persistent_poor ~ gt::px(145)
    ) %>%
    gt_theme()
  save_gt(table4_gt, "table4_sociodemographic_profile")

  table5_gt <- table5_transitions %>%
      gt::gt(groupname_col = "transition") %>%
      gt::tab_header(
        title = gt::md("**Table 5. Poverty transition matrix**"),
        subtitle = "Weighted row probabilities; 60% poverty threshold."
      ) %>%
      gt::fmt_number(columns = weighted_n, decimals = 0, use_seps = TRUE) %>%
      gt::fmt_percent(columns = row_probability, decimals = 1) %>%
      gt_theme()
  save_gt(table5_gt, "table5_transition_matrices") }

poverty_trend <- make_poverty_trend_figure(table1_poverty_rates)
ggsave(
  filename = file.path(project$figure_dir, "figure1_poverty_trends_ci.png"),
  plot = poverty_trend,
  width = 7.2,
  height = 4.6,
  dpi = 320 )
ggsave(
  filename = file.path(project$figure_dir, "figure1_poverty_trends_ci.pdf"),
  plot = poverty_trend,
  width = 7.2,
  height = 4.6 )

figure1_tex <- c(
  "\\begin{figure}[htbp]",
  "  \\centering",
  "  \\includegraphics[width=0.85\\textwidth]{figures/figure1_poverty_trends_ci.png}",
  "  \\caption{At-risk-of-poverty trends by poverty threshold, Türkiye SILC balanced panel, 2016--2019}",
  "  \\label{fig:poverty-trends}",
  "\\end{figure}"
)
writeLines(figure1_tex, file.path(project$figure_dir, "figure1_poverty_trends_ci.tex"))

latex_snippets <- c(
  "% Auto-generated LaTeX snippets for SILC poverty dynamics outputs.",
  "% In the main thesis preamble, useful packages include: booktabs, graphicx, longtable, caption.",
  "",
  "\\input{tables/table1_poverty_rates_fgt.tex}",
  "\\input{tables/table2_poverty_group_distribution.tex}",
  "\\input{tables/table3_poverty_duration.tex}",
  "\\input{tables/table4_sociodemographic_profile.tex}",
  "\\input{tables/table5_transition_matrices.tex}",
  "",
  "\\input{figures/figure1_poverty_trends_ci.tex}" )
writeLines(latex_snippets, file.path(project_root, "latex_inputs.tex"))

cat("\nCompleted SILC poverty dynamics workflow.\n")
cat("Individuals in balanced panel:", panel_audit$individuals, "\n")
cat("Main threshold:", project$main_threshold, "% of annual median equivalised income\n")
cat("Outputs written to:", normalizePath(project$out_dir), normalizePath(project$table_dir),
    normalizePath(project$figure_dir), normalizePath(project$model_dir), "\n")

