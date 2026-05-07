# =============================================================================
# Descriptive profile tables for poverty typologies
#
# This script produces two publication-ready profile tables:
#
#   Table 4a  (Section 1 – Individual characteristics):
#     Age group, Sex, Education, Employment status, Social security
#     registration, and Equivalised income — by poverty typology.
#
#   Table 4b  (Section 2 – Household & regional characteristics):
#     Household size, No. of children under 18, No. of elderly (65+),
#     No. of earners, and NUTS-2 region (with labels) — by poverty
#     typology.
#
# Both tables report weighted column percentages N (%) for categorical
# variables, weighted Mean (SD) for continuous variables, and p-values
# from weighted chi-squared (categorical) or weighted Kruskal–Wallis
# (continuous) tests across the mutually exclusive four-year typology.
#
# Outputs are saved as CSV files and, if the gt package is available, as
# HTML and LaTeX files.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(haven)
  library(purrr)
  library(readr)
  library(sandwich)
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

source(file.path(project_root, "R", "00_config.R"))
source(file.path(project_root, "R", "01_functions.R"))

project$data_path  <- file.path(project_root, project$data_path)
project$out_dir    <- file.path(project_root, project$out_dir)
project$table_dir  <- file.path(project_root, project$table_dir)
project$figure_dir <- file.path(project_root, project$figure_dir)
project$model_dir  <- file.path(project_root, project$model_dir)

dir.create(project$out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(project$table_dir, showWarnings = FALSE, recursive = TRUE)

# ---- NUTS-2 region labels (26 İBBS Düzey-2 regions) ------------------------

nuts2_labels <- c(
  "TR10" = "İstanbul",
  "TR21" = "Tekirdağ, Edirne, Kırklareli",
  "TR22" = "Balıkesir, Çanakkale",
  "TR31" = "İzmir",
  "TR32" = "Aydın, Denizli, Muğla",
  "TR33" = "Manisa, Afyonkarahisar, Kütahya, Uşak",
  "TR41" = "Bursa, Eskişehir, Bilecik",
  "TR42" = "Kocaeli, Sakarya, Düzce, Bolu, Yalova",
  "TR51" = "Ankara",
  "TR52" = "Konya, Karaman",
  "TR61" = "Antalya, Isparta, Burdur",
  "TR62" = "Adana, Mersin",
  "TR63" = "Hatay, Kahramanmaraş, Osmaniye",
  "TR71" = "Kırıkkale, Aksaray, Niğde, Nevşehir, Kırşehir",
  "TR72" = "Kayseri, Sivas, Yozgat",
  "TR81" = "Zonguldak, Karabük, Bartın",
  "TR82" = "Kastamonu, Çankırı, Sinop",
  "TR83" = "Samsun, Tokat, Çorum, Amasya",
  "TR90" = "Trabzon, Ordu, Giresun, Rize, Artvin, Gümüşhane",
  "TRA1" = "Erzurum, Erzincan, Bayburt",
  "TRA2" = "Ağrı, Kars, Iğdır, Ardahan",
  "TRB1" = "Malatya, Elazığ, Bingöl, Tunceli",
  "TRB2" = "Van, Muş, Bitlis, Hakkari",
  "TRC1" = "Gaziantep, Adıyaman, Kilis",
  "TRC2" = "Şanlıurfa, Diyarbakır",
  "TRC3" = "Mardin, Batman, Şırnak, Siirt"
)

# ---- Read data and run the main pipeline ------------------------------------

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
panel_poverty  <- add_poverty_status(panel_income, poverty_lines, vars, project$thresholds)

cat("Classifying poverty spells ...\n")
classified_main <- classify_poverty_spells(
  panel       = panel_poverty,
  vars        = vars,
  panel_years = project$panel_years,
  reference_year = project$reference_year,
  threshold   = project$main_threshold
)

# ---- Helper utilities -------------------------------------------------------

fmt_n   <- function(x) formatC(round(x), format = "f", digits = 0, big.mark = ",")
fmt_pct <- function(x) ifelse(is.na(x), "", sprintf("%.1f%%", 100 * x))
fmt_p   <- function(p) {
  p <- as.numeric(p)[1]
  case_when(is.na(p) ~ "-", p < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p))
}

w_sd <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  mu <- weighted_mean_safe(x[ok], w[ok])
  sqrt(sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok]))
}

weight <- vars$longitudinal_weight

cell_n_pct <- function(data, var, level) {
  x <- data[[var]]; w <- data[[weight]]
  ok <- !is.na(x) & !is.na(w) & w > 0
  denom <- sum(w[ok], na.rm = TRUE)
  if (denom == 0) return("")
  n   <- sum(x[ok] == level, na.rm = TRUE)
  pct <- sum(w[ok][x[ok] == level], na.rm = TRUE) / denom
  paste0(fmt_n(n), " (", fmt_pct(pct), ")")
}

cell_mean_sd <- function(data, var) {
  x <- data[[var]]; w <- data[[weight]]
  paste0(fmt_n(weighted_mean_safe(x, w)), " (", fmt_n(w_sd(x, w)), ")")
}

# ---- Build the profile data frame (reference year only) ---------------------

id   <- vars$person_id
year <- vars$year

profile <- panel_poverty %>%
  filter(.data[[year]] == project$reference_year) %>%
  left_join(
    classified_main %>%
      select(all_of(id), poverty_group, n_poor_4yr, current_poor, persistent_poor),
    by = id
  ) %>%
  mutate(
    profile_spell_group = case_when(
      n_poor_4yr == 0                             ~ "Never poor",
      n_poor_4yr %in% c(1, 2) & !persistent_poor  ~ "Transient poor",
      persistent_poor                              ~ "Persistent poor",
      TRUE                                         ~ "Frequently poor"
    ),
    poverty_group_for_test = factor(
      profile_spell_group,
      levels = c("Never poor", "Transient poor", "Frequently poor", "Persistent poor")
    ),
    # --- Individual-level recodings ---
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
      case_when(
        .data[[vars$social_security]] %in% codes$likely_informal_social_security_values
          ~ "Not registered",
        !is.na(.data[[vars$social_security]]) ~ "Registered",
        TRUE ~ NA_character_
      ),
      levels = c("Registered", "Not registered")
    ),
    # --- NUTS-2 region with labels ---
    nuts2_code = as.character(.data[[vars$nuts2]]),
    nuts2_labelled = factor(
      ifelse(
        nuts2_code %in% names(nuts2_labels),
        paste0(nuts2_code, " - ", nuts2_labels[nuts2_code]),
        nuts2_code
      )
    )
  )

# ---- Define sub-populations -------------------------------------------------

groups <- list(
  overall        = profile,
  never_poor     = profile %>% filter(poverty_group == "Never poor"),
  transient_poor = profile %>% filter(profile_spell_group == "Transient poor"),
  currently_poor = profile %>% filter(current_poor),
  persistent_poor = profile %>% filter(persistent_poor)
)

headers <- tibble(
  col   = names(groups),
  label = c(
    paste0("**Overall**  \nN = ", fmt_n(nrow(groups$overall))),
    paste0("**Never poor**  \nN = ", fmt_n(nrow(groups$never_poor))),
    paste0("**Transient poor**  \nN = ", fmt_n(nrow(groups$transient_poor))),
    paste0("**Currently poor**  \n(2019)  \nN = ", fmt_n(nrow(groups$currently_poor))),
    paste0("**Persistent poor**  \n*(subset of current)*  \nN = ", fmt_n(nrow(groups$persistent_poor)))
  )
)

# ---- p-value functions (survey-weighted) ------------------------------------

weighted_chisq_p <- function(var) {
  df <- profile %>% filter(!is.na(.data[[var]]), !is.na(poverty_group_for_test))
  if (n_distinct(df[[var]]) < 2 || n_distinct(df$poverty_group_for_test) < 2) return(NA_real_)
  des <- survey::svydesign(ids = ~1, weights = as.formula(paste0("~", weight)), data = df)
  out <- tryCatch(
    survey::svychisq(as.formula(paste0("~ poverty_group_for_test + ", var)),
                     design = des, statistic = "F"),
    error = function(e) NULL
  )
  if (is.null(out)) NA_real_ else out$p.value
}

weighted_kw_p <- function(var) {
  df <- profile %>% filter(!is.na(.data[[var]]), !is.na(poverty_group_for_test))
  if (n_distinct(df$poverty_group_for_test) < 2) return(NA_real_)
  des <- survey::svydesign(ids = ~1, weights = as.formula(paste0("~", weight)), data = df)
  out <- tryCatch(
    survey::svyranktest(as.formula(paste0(var, " ~ poverty_group_for_test")),
                        design = des, test = "KruskalWallis"),
    error = function(e) NULL
  )
  if (is.null(out)) NA_real_ else out$p.value
}

# ---- Row builders -----------------------------------------------------------

section_row <- function(label) {
  tibble(
    row_type = "section", variable = label, p_value = "",
    overall = "", never_poor = "", transient_poor = "",
    currently_poor = "", persistent_poor = ""
  )
}

categorical_rows <- function(section, var, levels) {
  p <- fmt_p(weighted_chisq_p(var))
  bind_rows(
    section_row(section),
    purrr::map_dfr(seq_along(levels), function(i) {
      level <- levels[[i]]
      tibble(
        row_type = "category",
        variable = paste0("  ", level),
        p_value  = if_else(i == 1L, p, ""),
        overall        = cell_n_pct(groups$overall,        var, level),
        never_poor     = cell_n_pct(groups$never_poor,     var, level),
        transient_poor = cell_n_pct(groups$transient_poor, var, level),
        currently_poor = cell_n_pct(groups$currently_poor, var, level),
        persistent_poor = cell_n_pct(groups$persistent_poor, var, level)
      )
    })
  )
}

continuous_row <- function(label, var) {
  tibble(
    row_type = "continuous", variable = label,
    p_value        = fmt_p(weighted_kw_p(var)),
    overall        = cell_mean_sd(groups$overall,        var),
    never_poor     = cell_mean_sd(groups$never_poor,     var),
    transient_poor = cell_mean_sd(groups$transient_poor, var),
    currently_poor = cell_mean_sd(groups$currently_poor, var),
    persistent_poor = cell_mean_sd(groups$persistent_poor, var)
  )
}

# =============================================================================
# TABLE 4a — Section 1: Individual socio-economic characteristics
# =============================================================================

cat("Building Table 4a (individual characteristics) ...\n")

table4a <- bind_rows(
  categorical_rows("Age group",             "age_group",                levels(profile$age_group)),
  categorical_rows("Sex",                   "sex_recoded",              levels(profile$sex_recoded)),
  categorical_rows("Education (ISCED)",     "education_recoded",        levels(profile$education_recoded)),
  categorical_rows("Employment status",     "labour_recoded",           levels(profile$labour_recoded)),
  categorical_rows("Social security",       "social_security_recoded",  levels(profile$social_security_recoded)),
  continuous_row("Eq. income (TL), mean (SD)", "eq_income")
)

attr(table4a, "headers") <- headers

# =============================================================================
# TABLE 4b — Section 2: Household & regional characteristics
# =============================================================================

cat("Building Table 4b (household & regional characteristics) ...\n")

table4b <- bind_rows(
  continuous_row("Household size, mean (SD)",            "hh_size"),
  continuous_row("No. children under 18, mean (SD)",     "hh_children_u18"),
  continuous_row("No. elderly (65+), mean (SD)",         "hh_elderly_65plus"),
  continuous_row("No. earners (proxy), mean (SD)",       "hh_earners_proxy"),
  categorical_rows(
    "NUTS-2 region",
    "nuts2_labelled",
    sort(levels(profile$nuts2_labelled))
  )
)

attr(table4b, "headers") <- headers

# ---- Save CSV outputs -------------------------------------------------------

write_csv(table4a, file.path(project$table_dir, "table4a_individual_profile.csv"))
write_csv(table4b, file.path(project$table_dir, "table4b_household_regional_profile.csv"))

cat("CSV tables saved to:", normalizePath(project$table_dir), "\n")

# ---- Optional: gt formatted tables -----------------------------------------

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

  format_profile_gt <- function(tbl, title, subtitle, headers) {
    tbl %>%
      gt::gt() %>%
      gt::cols_hide(columns = row_type) %>%
      gt::cols_label(
        variable        = "",
        p_value         = gt::md("p-value^1^"),
        overall         = gt::md(headers$label[headers$col == "overall"]),
        never_poor      = gt::md(headers$label[headers$col == "never_poor"]),
        transient_poor  = gt::md(headers$label[headers$col == "transient_poor"]),
        currently_poor  = gt::md(headers$label[headers$col == "currently_poor"]),
        persistent_poor = gt::md(headers$label[headers$col == "persistent_poor"])
      ) %>%
      gt::tab_header(title = gt::md(title), subtitle = gt::md(subtitle)) %>%
      gt::tab_source_note(
        source_note = gt::md(
          paste0(
            "^1^ N (%) = unweighted count (weighted column percentage). ",
            "Mean (SD) = weighted mean (weighted standard deviation). ",
            "Transient poor = poor in one or two years but not persistently poor. ",
            "P-values are computed over the mutually exclusive four-year typology: ",
            "never poor, transient poor, frequently poor, persistent poor."
          )
        )
      ) %>%
      gt::tab_style(
        style     = list(gt::cell_text(weight = "bold")),
        locations = gt::cells_body(rows = row_type == "section")
      ) %>%
      gt::tab_style(
        style     = gt::cell_borders(sides = "left", color = "#808080", weight = gt::px(1)),
        locations = list(
          gt::cells_body(columns = p_value),
          gt::cells_column_labels(columns = p_value)
        )
      ) %>%
      gt::tab_style(
        style     = gt::cell_borders(sides = "left", color = "#d9d9d9", weight = gt::px(1)),
        locations = list(
          gt::cells_body(columns = currently_poor),
          gt::cells_column_labels(columns = currently_poor)
        )
      ) %>%
      gt::cols_align(align = "left", columns = everything()) %>%
      gt::cols_width(
        variable        ~ gt::px(190),
        p_value         ~ gt::px(90),
        overall         ~ gt::px(145),
        never_poor      ~ gt::px(145),
        transient_poor  ~ gt::px(145),
        currently_poor  ~ gt::px(145),
        persistent_poor ~ gt::px(145)
      ) %>%
      gt_theme()
  }

  # Table 4a – individual characteristics
  table4a_gt <- format_profile_gt(
    table4a,
    title    = "**Table 4a. Individual socio-economic profile by poverty type, 2019**",
    subtitle = paste0(
      "*Note:* Weighted column percentages shown as N (%). ",
      "Currently poor (2019) = poor in 2019 regardless of spell history. ",
      "Persistent poor is a strict subset of currently poor. ",
      "Source: TR-SILC 2016–2019."
    ),
    headers  = headers
  )
  save_gt(table4a_gt, "table4a_individual_profile")

  # Table 4b – household & regional characteristics (wider for NUTS-2 labels)
  table4b_gt <- table4b %>%
    gt::gt() %>%
    gt::cols_hide(columns = row_type) %>%
    gt::cols_label(
      variable        = "",
      p_value         = gt::md("p-value^1^"),
      overall         = gt::md(headers$label[headers$col == "overall"]),
      never_poor      = gt::md(headers$label[headers$col == "never_poor"]),
      transient_poor  = gt::md(headers$label[headers$col == "transient_poor"]),
      currently_poor  = gt::md(headers$label[headers$col == "currently_poor"]),
      persistent_poor = gt::md(headers$label[headers$col == "persistent_poor"])
    ) %>%
    gt::tab_header(
      title    = gt::md("**Table 4b. Household & regional profile by poverty type, 2019**"),
      subtitle = gt::md(
        paste0(
          "*Note:* Weighted column percentages shown as N (%). ",
          "Mean (SD) = weighted mean (weighted standard deviation). ",
          "NUTS-2 region labels follow TurkStat İBBS Düzey-2 classification. ",
          "Source: TR-SILC 2016–2019."
        )
      )
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(
        paste0(
          "^1^ P-values are computed over the mutually exclusive four-year typology: ",
          "never poor, transient poor, frequently poor, persistent poor."
        )
      )
    ) %>%
    gt::tab_style(
      style     = list(gt::cell_text(weight = "bold")),
      locations = gt::cells_body(rows = row_type == "section")
    ) %>%
    gt::tab_style(
      style     = gt::cell_borders(sides = "left", color = "#808080", weight = gt::px(1)),
      locations = list(
        gt::cells_body(columns = p_value),
        gt::cells_column_labels(columns = p_value)
      )
    ) %>%
    gt::tab_style(
      style     = gt::cell_borders(sides = "left", color = "#d9d9d9", weight = gt::px(1)),
      locations = list(
        gt::cells_body(columns = currently_poor),
        gt::cells_column_labels(columns = currently_poor)
      )
    ) %>%
    gt::cols_align(align = "left", columns = everything()) %>%
    gt::cols_width(
      variable        ~ gt::px(280),
      p_value         ~ gt::px(90),
      overall         ~ gt::px(145),
      never_poor      ~ gt::px(145),
      transient_poor  ~ gt::px(145),
      currently_poor  ~ gt::px(145),
      persistent_poor ~ gt::px(145)
    ) %>%
    gt_theme()
  save_gt(table4b_gt, "table4b_household_regional_profile")

  cat("gt tables (HTML + LaTeX) saved.\n")
}

cat("\nDescriptive profile tables completed.\n")
