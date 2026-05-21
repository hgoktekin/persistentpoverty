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

source(file.path(project_root, "00_config.R"))
source(file.path(project_root, "01_functions.R"))

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
raw <- read_dta('/Users/haticegoktekin/Desktop/phd application/lisans tez/panel_16_19.dta')

#summary(raw[, c("fk070", "fk090", "hg110", "fg140", "fi190")])

#attrition_results <- calculate_attrition_rate(raw, vars, codes, project$panel_years, project$reference_year)
#cat("Attrition rate:", attrition_results$attrition_rate_percent, "%\n")
#write_csv(attrition_results$attrition_table, file.path(project$table_dir, "attrition_comparison.csv"))

cat("Step 2: constructing balanced panel and propagating longitudinal weights\n")
panel_balanced <- construct_balanced_panel(
  raw = raw,
  vars = vars,
  panel_years = project$panel_years,
  reference_year = project$reference_year)

panel_balance_check <- panel_balanced %>%
  summarise(
    individuals = n_distinct(.data[[vars$person_id]]),
    person_years = n(),
    first_year = min(.data[[vars$year]]),
    last_year = max(.data[[vars$year]]),
    duplicated_person_years = sum(duplicated(paste(.data[[vars$person_id]], .data[[vars$year]]))),
    missing_weights = sum(is.na(.data[[vars$longitudinal_weight]])))

cat("Step 3: calculating modified OECD equivalence scale and equivalised income\n")
panel_income <- panel_balanced %>%
  add_equivalised_income(vars) %>%
  add_household_context(vars, codes) %>%
  add_employment_stability(vars)

cat("Step 4: computing annual within-sample poverty thresholds\n")
poverty_lines <- compute_poverty_lines(panel_income, vars, project$thresholds)

panel_poverty <- add_poverty_status(
  panel = panel_income,
  poverty_lines = poverty_lines,
  vars = vars,
  thresholds = project$thresholds)


# =============================================================================
# Attrition diagnostics: Step 1 and Step 2
# Jenkins & Van Kerm-style representativeness check
# =============================================================================

xs_weight_var <- NA_character_   # replace with a cross-sectional weight if you have one

id <- vars$person_id
hh <- vars$household_id
year <- vars$year
age <- vars$age
ref <- vars$household_reference_person
income <- vars$household_income
weight_long <- vars$longitudinal_weight
threshold <- project$main_threshold / 100
first_year <- min(project$panel_years)
ref_year <- project$reference_year

# ---- Step 1: Describe attrition ----

balanced_ids <- panel_balanced %>%
  distinct(.data[[id]]) %>%
  pull(.data[[id]])

sample_flow <- raw %>%
  filter(.data[[year]] %in% project$panel_years) %>%
  group_by(.data[[year]]) %>%
  summarise(
    n_individuals = n_distinct(.data[[id]]),
    .groups = "drop"
  ) %>%
  mutate(
    retention_from_wave1 = n_individuals / first(n_individuals),
    attrition_from_wave1 = 1 - retention_from_wave1
  )

sample_flow <- sample_flow %>%
  mutate(
    wave_to_wave_retention = n_individuals / lag(n_individuals),
    wave_to_wave_attrition = 1 - wave_to_wave_retention
  )

baseline_raw <- panel_income %>%
  filter(.data[[year]] == first_year)



baseline_poverty_line <- weighted_median(
  baseline_raw$eq_income_raw,
  baseline_raw$xs_weight
) * threshold

baseline_attrition <- baseline_raw %>%
  mutate(
    stayed_balanced_panel = .data[[id]] %in% balanced_ids,
    attrited = !stayed_balanced_panel,
    poor_baseline = eq_income_raw < baseline_poverty_line,
    employment_group = if_else(
      .data[[vars$labour_status]] %in% codes$employed_labour_status,
      "Employed", "Unemployed", missing = NA_character_
    ) )

attrition_by_group <- bind_rows(
  baseline_attrition %>%
    group_by(group = "Baseline poverty", category = if_else(poor_baseline, "Poor", "Non-poor")) %>%
    summarise(n = n(), attrition_rate = mean(attrited), .groups = "drop"),
  
  baseline_attrition %>%
    filter(!is.na(employment_group)) %>%
    group_by(group = "Employment", category = employment_group) %>%
    summarise(n = n(), attrition_rate = mean(attrited), .groups = "drop")
)

sample_flow
attrition_by_group

baseline_attrition %>%
  summarise(
  attrition_unweighted = mean(attrited),
  attrition_weighted = weighted_mean_safe(attrited, xs_weight)
)
# ---- Step 2: Representativeness test ----

wave4_raw <- panel_income %>%
  filter(.data[[year]] == ref_year) %>%
  mutate(
    xs_weight = if (!is.na(xs_weight_var) && xs_weight_var %in% names(raw)) {
      .data[[xs_weight_var]]
    } else {
      1
    },
    oecd_weight = case_when(
      .data[[ref]] == 1 ~ 1.0,
      .data[[ref]] != 1 & .data[[age]] >= 14 ~ 0.5,
      .data[[ref]] != 1 & .data[[age]] < 14 ~ 0.3,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(.data[[hh]], .data[[year]]) %>%
  mutate(
    hh_eq_size = sum(oecd_weight, na.rm = TRUE),
    eq_income_raw = .data[[income]] / hh_eq_size
  ) %>%
  ungroup() %>%
  filter(!is.na(eq_income_raw), hh_eq_size > 0)

wave4_poverty_line <- weighted_median(
  wave4_raw$eq_income_raw,
  wave4_raw$xs_weight ) * threshold

cross_sectional_rate <- weighted_mean_safe(
  wave4_raw$eq_income_raw < wave4_poverty_line,
  wave4_raw$xs_weight)

longitudinal_rate <- panel_income %>%
  filter(
    .data[[year]] == ref_year,
    .data[[id]] %in% balanced_ids
  ) %>%
  summarise(
    poverty_rate = weighted_mean_safe(
      eq_income < wave4_poverty_line,
      .data[[weight_long]]
    )
  ) %>%
  pull(poverty_rate)


representativeness_test <- tibble(
  reference_year = ref_year,
  poverty_threshold = project$main_threshold,
  cross_sectional_rate = cross_sectional_rate,
  longitudinal_rate = longitudinal_rate,
  difference = cross_sectional_rate - longitudinal_rate,
  interpretation = case_when(
    abs(difference) < 0.005 ~ "Balanced panel is close to representative",
    difference > 0 ~ "Balanced panel may under-represent poorer people",
    difference < 0 ~ "Balanced panel may over-represent poorer people"
  )
)

representativeness_test

chisq.test(
  table(
    baseline_attrition$poor_baseline,
    baseline_attrition$attrited
  )
)
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
  "Frequently poor", "poor in at least two years, but not in the current-year",
  "Current poor (2019)", "equivalised income in 2019 is below the selected poverty threshold",
  "Persistent poor", "current poor in 2019 and poor in at least two of 2016, 2017, and 2018")

cat("Step 6: creating descriptive poverty tables and FGT indices\n")

# poverty rates with 60% threshold on BALANCED PANEL 
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
    mean_spell_duration = weighted_mean_safe(mean_spell_duration, panel_weight))

table4_profile <- make_profile_table(
  panel = panel_poverty,
  classified = classified_main,
  vars = vars,
  codes = codes,
  reference_year = project$reference_year,
  threshold = project$main_threshold)

table5_transitions <- make_transition_matrices(
  panel = panel_poverty,
  vars = vars,
  panel_years = project$panel_years,
  threshold = project$main_threshold,
  codes = codes)

mobility_summary <- table5_transitions %>%
 filter(transition_type %in% c("Entry", "Exit", "Poverty persistence")) %>%
 select(transition, transition_type, from_status, to_status, weighted_n, row_probability)

table6_income_composition <- make_income_composition_table(
  panel = panel_poverty,
  classified = classified_main,
  vars = vars,
  income_components = income_components,
  reference_year = project$reference_year)

table7_income_categories <- make_income_category_table(
  panel = panel_poverty,
  classified = classified_main,
  vars = vars,
  income_categories = income_categories,
  reference_year = project$reference_year
)


cat("Step 8: exporting tables, model outputs, and figure\n")
write_csv(panel_balance_check, file.path(project$out_dir, "panel_audit.csv"))
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
write_csv(table6_income_composition, file.path(project$table_dir, "table6_income_composition_by_poverty_profile.csv"))
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
        column_labels.border.bottom.width = gt::px(1)) }

  save_gt <- function(gt_tbl, stem) {
    gt::gtsave(gt_tbl, file.path(project$table_dir, paste0(stem, ".html")))
    gt::gtsave(gt_tbl, file.path(project$table_dir, paste0(stem, ".tex")))}

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
    select(table_family, category, unweighted_n, unweighted_share) %>%
      gt::gt(groupname_col = "table_family") %>%
      gt::tab_header(
        title = gt::md("**Table 2. Poverty group distribution**"),
        subtitle = "Unweighted counts and shares. Persistent poor is shown both in the mutually exclusive typology and as a subset of current poverty."
      ) %>%
    gt::cols_label(
      category = "Poverty profile",
      unweighted_n = "Unweighted N",
      unweighted_share = "Share"
    ) %>%
    gt::fmt_number(columns = unweighted_n, decimals = 0, use_seps = TRUE) %>%
    gt::fmt_percent(columns = unweighted_share, decimals = 2) %>%
      gt_theme()
  save_gt(table2_gt, "table2_poverty_group_distribution")

  table3_gt <- table3_duration %>%
    select(duration_years, unweighted_n, unweighted_share) %>%
      gt::gt() %>%
      gt::tab_header(
        title = gt::md("**Table 3. Poverty duration**"),
        subtitle = "Number of years below the 60% poverty threshold."
      ) %>%
    gt::cols_label(
      duration_years = "Years poor",
      unweighted_n = "Unweighted N",
      unweighted_share = "Share"
    ) %>%
    gt::fmt_number(columns = unweighted_n, decimals = 0, use_seps = TRUE) %>%
    gt::fmt_percent(columns = unweighted_share, decimals = 2) %>%
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
        "*Note:* Unweighted N and shares shown as N (%). Continuous cells show weighted mean (weighted SD); unweighted min-max. p-values: weighted chi-squared (categorical); weighted Kruskal-Wallis rank test (continuous). Currently poor (2019) = poor in 2019 regardless of spell history. Persistent poor is a strict subset of currently poor (2019). Source: TR-SILC 2016--2019."      )
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        "^1^ For this profile table, transient poor means poor in one or two years but not persistently poor. P-values are computed over the mutually exclusive four-year profile typology: never poor, transient poor, frequently poor, persistent poor."
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
  
  table5_display <- table5_transitions %>%
    mutate(
      table_group = paste0(formal_status, ": ", transition),
      from_status = factor(from_status, levels = c("Poor", "Non-poor")),
      to_status = factor(to_status, levels = c("Poor", "Non-poor"))
    ) %>%
    select(table_group, from_status, to_status, row_probability) %>%
    tidyr::pivot_wider(names_from = to_status, values_from = row_probability) %>%
    arrange(table_group, from_status)
  
  table5_gt <- table5_display %>%
    gt::gt(groupname_col = "table_group") %>%
      gt::tab_header(
        title = gt::md("**Table 5. Weighted poverty transition matrices by formality status**"),
        subtitle = "Row-conditional probabilities; 60% poverty threshold."
      ) %>%
    gt::cols_label(
      from_status = "Status at t",
      Poor = "Poor (t+1)",
      `Non-poor` = "Non-poor (t+1)"
      ) %>%
    gt::tab_spanner(
      label = "Status at t+1",
      columns = c(Poor, `Non-poor`)
    ) %>%
    gt::fmt_percent(columns = c(Poor, `Non-poor`), decimals = 1) %>%
    gt::tab_source_note(
      source_note = gt::md(
        "Row-conditional transition probabilities weighted by 4-year longitudinal weights. Formality status is measured at the origin year t using fi190; informal is defined by the configured informal social-security code. Persistence rate = P(poor t+1 | poor t); entry rate = P(poor t+1 | non-poor t). Jenkins (2000); Cappellari & Jenkins (2004)."
      )
    ) %>%
      gt_theme()
  
  save_gt(table5_gt, "table5_transition_matrices")

  table6_gt <- table6_income_composition %>%
    select(poverty_profile, income_variable, income_type, income_role, unweighted_n,
           unweighted_positive_n, mean_income, income_share) %>%
    gt::gt(groupname_col = "poverty_profile") %>%
    gt::tab_header(
      title = gt::md("**Table 6. Income composition by poverty profile, 2019**"),
      subtitle = "Income-component means and shares use longitudinal person weights. Missing component values are treated as zero."
    ) %>%
    gt::cols_label(
      income_variable = "Variable",
      income_type = "Income type",
      income_role = "Role",
      unweighted_n = "Unweighted N",
      unweighted_positive_n = "N > 0",
      mean_income = "Mean income",
      income_share = "Share of total income"
    ) %>%
    gt::fmt_number(columns = c(unweighted_n, unweighted_positive_n), decimals = 0, use_seps = TRUE) %>%
    gt::fmt_number(columns = mean_income, decimals = 2, use_seps = TRUE) %>%
    gt::fmt_percent(columns = income_share, decimals = 2) %>%
    gt_theme()
  save_gt(table6_gt, "table6_income_composition_by_poverty_profile")
  
  table7_display <- table7_income_categories %>%
    mutate(cell = paste0(
      formatC(mean_income, format = "f", digits = 2, big.mark = ","),
      " (",
      sprintf("%.2f%%", 100 * income_share),
      ")"
    )) %>%
    select(income_category, poverty_profile, cell) %>%
    tidyr::pivot_wider(names_from = poverty_profile, values_from = cell)
  
  table7_gt <- table7_display %>%
    gt::gt() %>%
    gt::tab_header(
      title = gt::md("**Table 7. Income categories by poverty profile, 2019**"),
      subtitle = "Cells show weighted mean income and share of the five-category income total: mean (share)."
    ) %>%
    gt::cols_label(
      income_category = "Income category",
      Overall = "Overall",
      `Never poor` = "Never poor",
      `Transient poor` = "Transient poor",
      `Persistent poor` = "Persistent poor"
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        "Wage income, entrepreneurial income, capital income, social transfers, and private transfers are summed from the configured component variables. Totals, taxes paid, and transfers given to other households are excluded from the five-category income total. Missing component values are treated as zero."
      )
    ) %>%
    gt_theme()
  save_gt(table7_gt, "table7_income_categories_by_poverty_profile")

}

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
  "\\input{tables/table6_income_composition_by_poverty_profile.tex}",
  "\\input{tables/table7_income_categories_by_poverty_profile.tex}",
  "",
  "\\input{figures/figure1_poverty_trends_ci.tex}" )
writeLines(latex_snippets, file.path(project_root, "latex_inputs.tex"))

cat("\nCompleted SILC poverty dynamics workflow.\n")
cat("Individuals in balanced panel:", panel_balance_check$individuals, "\n")
cat("Main threshold:", project$main_threshold, "% of annual median equivalised income\n")
cat("Outputs written to:", normalizePath(project$out_dir), normalizePath(project$table_dir),
    normalizePath(project$figure_dir), normalizePath(project$model_dir), "\n")

