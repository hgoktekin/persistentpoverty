# =============================================================================
# Reusable functions for poverty dynamics analysis
#
# Statistical conventions used here:
# - Modified OECD equivalence scale: 1.0 for household reference person,
#   0.5 for other members aged 14+, and 0.3 for children under 14.
# - Poverty thresholds: 50%, 60%, and 70% of the within-sample weighted median
#   of equivalised disposable household income in each year.
# - Persistent at-risk-of-poverty: poor in 2019 and poor in at least two of
#   2016, 2017, and 2018, following the Eurostat longitudinal definition.
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

require_vars <- function(data, required) {
  missing <- setdiff(unlist(required, use.names = FALSE), names(data))
  if (length(missing) > 0) {
    stop("Missing required variable(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_var_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- weighted_mean_safe(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

weighted_median <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

construct_balanced_panel <- function(raw, vars, panel_years, reference_year) {
  require_vars(raw, vars[c(
    "person_id", "household_id", "year", "household_income",
    "longitudinal_weight", "age", "household_reference_person"
  )])

  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight

  # TurkStat's 4-year longitudinal weight is attached to the reference wave.
  # We identify the longitudinal sample from non-missing reference-year weights,
  # then propagate the same panel weight back to all four person-years.
  ids_4year <- raw %>%
    filter(.data[[year]] == reference_year, !is.na(.data[[weight]]), .data[[weight]] > 0) %>%
    distinct(.data[[id]]) %>%
    pull(.data[[id]])

  weights_reference <- raw %>%
    filter(.data[[year]] == reference_year, .data[[id]] %in% ids_4year,
           !is.na(.data[[weight]]), .data[[weight]] > 0) %>%
    arrange(.data[[id]]) %>%
    distinct(.data[[id]], .keep_all = TRUE) %>%
    transmute(!!id := .data[[id]], panel_weight_reference = .data[[weight]])

  panel <- raw %>%
    filter(.data[[id]] %in% ids_4year, .data[[year]] %in% panel_years) %>%
    arrange(.data[[id]], .data[[year]]) %>%
    distinct(.data[[id]], .data[[year]], .keep_all = TRUE) %>%
    select(-all_of(weight)) %>%
    left_join(weights_reference, by = id) %>%
    rename(!!weight := panel_weight_reference)

  panel_check <- panel %>% count(.data[[id]], name = "n_years")

  if (!all(panel_check$n_years == length(panel_years))) {
    stop(
      "Panel is not balanced after deduplication. ",
      "Use the saved audit object or inspect duplicate person-year records.",
      call. = FALSE
    )
  }

  if (any(is.na(panel[[weight]]) | panel[[weight]] <= 0)) {
    stop("Missing or non-positive longitudinal weights after propagation.", call. = FALSE)
  }

  panel
}

add_equivalised_income <- function(panel, vars) {
  hh <- vars$household_id
  year <- vars$year
  age <- vars$age
  ref <- vars$household_reference_person
  income <- vars$household_income

  out <- panel %>%
    mutate(
      oecd_member_weight = case_when(
        .data[[ref]] == 1 ~ 1.0,
        .data[[ref]] != 1 & .data[[age]] >= 14 ~ 0.5,
        .data[[ref]] != 1 & .data[[age]] < 14 ~ 0.3,
        TRUE ~ NA_real_
      )
    ) %>%
    group_by(.data[[hh]], .data[[year]]) %>%
    mutate(
      hh_eq_size = sum(oecd_member_weight, na.rm = TRUE),
      hh_n_children_u14 = sum(.data[[age]] < 14, na.rm = TRUE),
      hh_n_adults_14plus = sum(.data[[age]] >= 14, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(eq_income = .data[[income]] / hh_eq_size)

  if (any(is.na(out$eq_income) | out$hh_eq_size <= 0)) {
    stop("Invalid equivalence scale or equivalised income. Check income, age, and household composition.", call. = FALSE)
  }

  out
}

add_household_context <- function(panel, vars, codes) {
  id <- vars$person_id
  hh <- vars$household_id
  year <- vars$year
  age <- vars$age
  labour <- vars$labour_status

  panel %>%
    mutate(
      is_child_u14 = .data[[age]] < 14,
      is_child_u18 = .data[[age]] < 18,
      is_elderly_65plus = .data[[age]] >= 65,
      is_earner_proxy = if_else(.data[[labour]] %in% codes$employed_labour_status, 1L, 0L, missing = 0L)
    ) %>%
    group_by(.data[[hh]], .data[[year]]) %>%
    mutate(
      hh_size = n_distinct(.data[[id]]),
      hh_children_u14 = sum(is_child_u14, na.rm = TRUE),
      hh_children_u18 = sum(is_child_u18, na.rm = TRUE),
      hh_elderly_65plus = sum(is_elderly_65plus, na.rm = TRUE),
      hh_earners_proxy = sum(is_earner_proxy, na.rm = TRUE),
      household_type_derived = case_when(
        hh_size == 1 ~ "Single person",
        hh_children_u18 > 0 & hh_elderly_65plus > 0 ~ "Children and elderly present",
        hh_children_u18 > 0 ~ "Children present",
        hh_elderly_65plus > 0 ~ "Elderly present",
        TRUE ~ "Adults only"
      )
    ) %>%
    ungroup()
}

compute_poverty_lines <- function(panel, vars, thresholds) {
  year <- vars$year
  weight <- vars$longitudinal_weight

  panel %>%
    group_by(.data[[year]]) %>%
    summarise(
      median_eq_income = weighted_median(eq_income, .data[[weight]]),
      .groups = "drop"
    ) %>%
    mutate(
      across(
        median_eq_income,
        list(
          poverty_line_50 = ~ .x * 0.50,
          poverty_line_60 = ~ .x * 0.60,
          poverty_line_70 = ~ .x * 0.70
        ),
        .names = "{.fn}"
      )
    ) %>%
    select(all_of(year), median_eq_income, paste0("poverty_line_", thresholds))
}

add_poverty_status <- function(panel, poverty_lines, vars, thresholds) {
  year <- vars$year

  out <- panel %>% left_join(poverty_lines, by = year)

  for (threshold in thresholds) {
    line_var <- paste0("poverty_line_", threshold)
    poor_var <- paste0("poor_", threshold)
    gap_var <- paste0("poverty_gap_", threshold)

    out <- out %>%
      mutate(
        "{poor_var}" := as.integer(eq_income < .data[[line_var]]),
        "{gap_var}" := pmax((.data[[line_var]] - eq_income) / .data[[line_var]], 0)
      )
  }

  out
}

classify_poverty_spells <- function(panel, vars, panel_years, reference_year, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  poor <- paste0("poor_", threshold)
  previous_years <- setdiff(panel_years, reference_year)

  wide <- panel %>%
    select(all_of(c(id, year, weight)), all_of(poor)) %>%
    pivot_wider(
      names_from = all_of(year),
      values_from = all_of(poor),
      names_prefix = "poor_"
    )

  poor_year_vars <- paste0("poor_", panel_years)
  previous_vars <- paste0("poor_", previous_years)
  current_var <- paste0("poor_", reference_year)

  wide %>%
    rowwise() %>%
    mutate(
      n_poor_4yr = sum(c_across(all_of(poor_year_vars)), na.rm = TRUE),
      n_poor_previous3 = sum(c_across(all_of(previous_vars)), na.rm = TRUE),
      current_poor = .data[[current_var]] == 1,
      never_poor = n_poor_4yr == 0,
      poor_once = n_poor_4yr == 1,
      poor_multiple = n_poor_4yr >= 2,
      persistent_poor = current_poor & n_poor_previous3 >= 2,
      current_poor_not_persistent = current_poor & !persistent_poor,
      # Mutually exclusive four-year typology:
      # never poor = never below the threshold;
      # transient poor = one poor year only;
      # persistent poor = Eurostat current-year anchored persistent poverty;
      # frequently poor = multiple poor years but not persistent by Eurostat rule.
      poverty_group = case_when(
        never_poor ~ "Never poor",
        poor_once ~ "Transient poor",
        persistent_poor ~ "Persistent poor",
        poor_multiple & !persistent_poor ~ "Frequently poor",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    mutate(
      threshold = threshold,
      poverty_group = factor(
        poverty_group,
        levels = c("Never poor", "Transient poor", "Frequently poor", "Persistent poor")
      )
    )
}

fgt_index <- function(income, poverty_line, weights, alpha = 0) {
  ok <- !is.na(income) & !is.na(poverty_line) & !is.na(weights) & weights > 0
  if (!any(ok)) return(NA_real_)
  gap <- pmax((poverty_line[ok] - income[ok]) / poverty_line[ok], 0)
  if (alpha == 0) {
    return(weighted_mean_safe(as.integer(gap > 0), weights[ok]))
  }
  weighted_mean_safe(gap^alpha, weights[ok])
}

make_table1_poverty_rates <- function(panel, vars, thresholds) {
  year <- vars$year
  weight <- vars$longitudinal_weight
  des <- svydesign(ids = ~1, weights = as.formula(paste0("~", weight)), data = panel)

  purrr::map_dfr(thresholds, function(threshold) {
    poor <- paste0("poor_", threshold)
    line <- paste0("poverty_line_", threshold)
    gap <- paste0("poverty_gap_", threshold)

    purrr::map_dfr(sort(unique(panel[[year]])), function(y) {
      d_y <- subset(des, get(year) == y)
      est <- svymean(as.formula(paste0("~", poor)), d_y, na.rm = TRUE)
      ci <- confint(est)
      data_y <- panel[panel[[year]] == y, ]

      tibble(
        year = y,
        threshold = threshold,
        poverty_line = unique(data_y[[line]])[1],
        poverty_rate = as.numeric(coef(est)[1]),
        poverty_rate_se = as.numeric(SE(est)[1]),
        poverty_rate_ci_low = as.numeric(ci[1, 1]),
        poverty_rate_ci_high = as.numeric(ci[1, 2]),
        FGT0_headcount = fgt_index(data_y$eq_income, data_y[[line]], data_y[[weight]], alpha = 0),
        FGT1_gap = fgt_index(data_y$eq_income, data_y[[line]], data_y[[weight]], alpha = 1),
        FGT2_severity = fgt_index(data_y$eq_income, data_y[[line]], data_y[[weight]], alpha = 2),
        unweighted_n = nrow(data_y),
        weighted_n = sum(data_y[[weight]], na.rm = TRUE)
      )
    })
  })
}

make_poverty_group_distribution <- function(classified, vars) {
  weight <- vars$longitudinal_weight
  total_w <- sum(classified[[weight]], na.rm = TRUE)

  mutual <- classified %>%
    count(poverty_group, wt = .data[[weight]], name = "weighted_n") %>%
    mutate(
      table_family = "Mutually exclusive four-year typology",
      population_share = weighted_n / total_w
    ) %>%
    rename(category = poverty_group)

  overlapping <- tibble(
    table_family = "Overlapping current-year flags",
    category = c("Current poor (2019)", "Persistent poor"),
    weighted_n = c(
      sum(classified[[weight]][classified$current_poor], na.rm = TRUE),
      sum(classified[[weight]][classified$persistent_poor], na.rm = TRUE)
    )
  ) %>%
    mutate(population_share = weighted_n / total_w)

  bind_rows(mutual, overlapping)
}

poverty_episode_summary <- function(status) {
  status <- as.integer(status)
  if (all(is.na(status))) {
    return(tibble(n_episodes = NA_integer_, max_spell_duration = NA_integer_, mean_spell_duration = NA_real_))
  }
  status[is.na(status)] <- 0L
  starts <- which(status == 1L & dplyr::lag(status, default = 0L) == 0L)
  ends <- which(status == 1L & dplyr::lead(status, default = 0L) == 0L)
  durations <- ends - starts + 1L
  tibble(
    n_episodes = length(durations),
    max_spell_duration = ifelse(length(durations) == 0, 0L, max(durations)),
    mean_spell_duration = ifelse(length(durations) == 0, 0, mean(durations))
  )
}

make_poverty_duration_table <- function(panel, classified, vars, panel_years, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  poor <- paste0("poor_", threshold)
  total_w <- sum(classified[[weight]], na.rm = TRUE)

  episodes <- panel %>%
    arrange(.data[[id]], .data[[year]]) %>%
    group_by(.data[[id]]) %>%
    summarise(
      poverty_episode_summary(.data[[poor]]),
      panel_weight = first(.data[[weight]]),
      .groups = "drop"
    )

  duration <- classified %>%
    mutate(duration_years = n_poor_4yr) %>%
    count(duration_years, wt = .data[[weight]], name = "weighted_n") %>%
    mutate(population_share = weighted_n / total_w)

  list(duration = duration, episodes = episodes)
}

make_transition_matrices <- function(panel, vars, panel_years, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  poor <- paste0("poor_", threshold)

  wide <- panel %>%
    select(all_of(c(id, year, weight)), all_of(poor)) %>%
    pivot_wider(names_from = all_of(year), values_from = all_of(poor), names_prefix = "poor_")

  purrr::map_dfr(seq_along(panel_years[-length(panel_years)]), function(i) {
    y0 <- panel_years[i]
    y1 <- panel_years[i + 1]
    from <- paste0("poor_", y0)
    to <- paste0("poor_", y1)

    wide %>%
      filter(!is.na(.data[[from]]), !is.na(.data[[to]])) %>%
      mutate(
        from_status = if_else(.data[[from]] == 1, "Poor", "Non-poor"),
        to_status = if_else(.data[[to]] == 1, "Poor", "Non-poor")
      ) %>%
      count(from_status, to_status, wt = .data[[weight]], name = "weighted_n") %>%
      group_by(from_status) %>%
      mutate(row_probability = weighted_n / sum(weighted_n)) %>%
      ungroup() %>%
      mutate(
        transition = paste(y0, y1, sep = "-"),
        transition_type = case_when(
          from_status == "Non-poor" & to_status == "Poor" ~ "Entry",
          from_status == "Poor" & to_status == "Non-poor" ~ "Exit",
          from_status == "Poor" & to_status == "Poor" ~ "Poverty persistence",
          TRUE ~ "Non-poor persistence"
        )
      )
  })
}

make_profile_table <- function(panel, classified, vars, codes, reference_year, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight

  fmt_n <- function(x) {
    formatC(round(x), format = "f", digits = 0, big.mark = ",")
  }

  fmt_pct <- function(x) {
    ifelse(is.na(x), "", sprintf("%.1f%%", 100 * x))
  }

  fmt_p <- function(p) {
    p <- as.numeric(p)[1]
    case_when(
      is.na(p) ~ "-",
      p < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p)
    )
  }

  w_sd <- function(x, w) {
    ok <- !is.na(x) & !is.na(w) & w > 0
    if (sum(ok) < 2) return(NA_real_)
    mu <- weighted_mean_safe(x[ok], w[ok])
    sqrt(sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok]))
  }

  cell_n_pct <- function(data, var, level) {
    x <- data[[var]]
    w <- data[[weight]]
    ok <- !is.na(x) & !is.na(w) & w > 0
    denom <- sum(w[ok], na.rm = TRUE)
    if (denom == 0) return("")
    n <- sum(x[ok] == level, na.rm = TRUE)
    pct <- sum(w[ok][x[ok] == level], na.rm = TRUE) / denom
    paste0(fmt_n(n), " (", fmt_pct(pct), ")")
  }

  cell_mean_sd <- function(data, var) {
    x <- data[[var]]
    w <- data[[weight]]
    paste0(fmt_n(weighted_mean_safe(x, w)), " (", fmt_n(w_sd(x, w)), ")")
  }

  profile <- panel %>%
    filter(.data[[year]] == reference_year) %>%
    left_join(
      classified %>%
        select(all_of(id), poverty_group, n_poor_4yr, current_poor, persistent_poor),
      by = id
    ) %>%
    mutate(
      # Used only for p-values. This is mutually exclusive, unlike the
      # descriptive "Currently poor (2019)" column.
      profile_spell_group = case_when(
        n_poor_4yr == 0 ~ "Never poor",
        n_poor_4yr %in% c(1, 2) & !persistent_poor ~ "Transient poor",
        persistent_poor ~ "Persistent poor",
        TRUE ~ "Frequently poor"
      ),
      poverty_group_for_test = factor(
        profile_spell_group,
        levels = c("Never poor", "Transient poor", "Frequently poor", "Persistent poor")
      ),
      sex_recoded = factor(
        case_when(
          .data[[vars$sex]] == codes$sex[["female"]] ~ "Female",
          .data[[vars$sex]] == codes$sex[["male"]] ~ "Male",
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
          .data[[vars$labour_status]] %in% c(1, 2) ~ "Employee",
          .data[[vars$labour_status]] %in% c(3, 4) ~ "Self-employed",
          .data[[vars$labour_status]] == 7 ~ "Retired",
          .data[[vars$labour_status]] %in% c(5, 6, 8, 9, 10) ~ "Inactive",
          TRUE ~ NA_character_
        ),
        levels = c("Employee", "Self-employed", "Retired", "Inactive")
      )
    )

  weighted_chisq_p <- function(var) {
    df <- profile %>%
      filter(!is.na(.data[[var]]), !is.na(poverty_group_for_test))
    if (n_distinct(df[[var]]) < 2 || n_distinct(df$poverty_group_for_test) < 2) {
      return(NA_real_)
    }
    des <- survey::svydesign(
      ids = ~1,
      weights = as.formula(paste0("~", weight)),
      data = df
    )
    out <- tryCatch(
      survey::svychisq(as.formula(paste0("~ poverty_group_for_test + ", var)),
                       design = des, statistic = "F"),
      error = function(e) NULL
    )
    if (is.null(out)) NA_real_ else out$p.value
  }

  weighted_kw_p <- function(var) {
    df <- profile %>%
      filter(!is.na(.data[[var]]), !is.na(poverty_group_for_test))
    if (n_distinct(df$poverty_group_for_test) < 2) {
      return(NA_real_)
    }
    des <- survey::svydesign(
      ids = ~1,
      weights = as.formula(paste0("~", weight)),
      data = df
    )
    out <- tryCatch(
      survey::svyranktest(as.formula(paste0(var, " ~ poverty_group_for_test")),
                          design = des, test = "KruskalWallis"),
      error = function(e) NULL
    )
    if (is.null(out)) NA_real_ else out$p.value
  }

  groups <- list(
    overall = profile,
    never_poor = profile %>% filter(poverty_group == "Never poor"),
    transient_poor = profile %>% filter(profile_spell_group == "Transient poor"),
    currently_poor = profile %>% filter(current_poor),
    persistent_poor = profile %>% filter(persistent_poor)
  )

  headers <- tibble(
    col = names(groups),
    label = c(
      paste0("**Overall**  \nN = ", fmt_n(nrow(groups$overall))),
      paste0("**Never poor**  \nN = ", fmt_n(nrow(groups$never_poor))),
      paste0("**Transient poor**  \nN = ", fmt_n(nrow(groups$transient_poor))),
      paste0("**Currently poor**  \n(2019)  \nN = ", fmt_n(nrow(groups$currently_poor))),
      paste0("**Persistent poor**  \n*(subset of current)*  \nN = ", fmt_n(nrow(groups$persistent_poor)))
    )
  )

  section_row <- function(label) {
    tibble(
      row_type = "section",
      variable = label,
      p_value = "",
      overall = "",
      never_poor = "",
      transient_poor = "",
      currently_poor = "",
      persistent_poor = ""
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
          p_value = if_else(i == 1L, p, ""),
          overall = cell_n_pct(groups$overall, var, level),
          never_poor = cell_n_pct(groups$never_poor, var, level),
          transient_poor = cell_n_pct(groups$transient_poor, var, level),
          currently_poor = cell_n_pct(groups$currently_poor, var, level),
          persistent_poor = cell_n_pct(groups$persistent_poor, var, level)
        )
      })
    )
  }

  continuous_row <- function(label, var) {
    tibble(
      row_type = "continuous",
      variable = label,
      p_value = fmt_p(weighted_kw_p(var)),
      overall = cell_mean_sd(groups$overall, var),
      never_poor = cell_mean_sd(groups$never_poor, var),
      transient_poor = cell_mean_sd(groups$transient_poor, var),
      currently_poor = cell_mean_sd(groups$currently_poor, var),
      persistent_poor = cell_mean_sd(groups$persistent_poor, var)
    )
  }

  out <- bind_rows(
    categorical_rows("Age group", "age_group", levels(profile$age_group)),
    categorical_rows("Sex", "sex_recoded", levels(profile$sex_recoded)),
    categorical_rows("Education (ISCED)", "education_recoded", levels(profile$education_recoded)),
    categorical_rows("Employment status", "labour_recoded", levels(profile$labour_recoded)),
    continuous_row("Eq. income (TL), mean (SD)", "eq_income")
  )

  attr(out, "headers") <- headers
  out
}

theme_thesis <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "serif", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.35),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3, linetype = "dashed"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, colour = "grey35", size = base_size - 2)
    )
}

make_poverty_trend_figure <- function(table1) {
  ggplot(table1, aes(x = year, y = poverty_rate, colour = factor(threshold), group = threshold)) +
    geom_ribbon(
      aes(ymin = poverty_rate_ci_low, ymax = poverty_rate_ci_high, fill = factor(threshold)),
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.9) +
    scale_x_continuous(breaks = sort(unique(table1$year))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_colour_manual(values = c("50" = "#1b1b1b", "60" = "#4d4d4d", "70" = "#8a8a8a")) +
    scale_fill_manual(values = c("50" = "#1b1b1b", "60" = "#4d4d4d", "70" = "#8a8a8a")) +
    labs(
      title = "At-risk-of-poverty rate by threshold, Türkiye SILC panel",
      x = "Survey year",
      y = "Weighted poverty rate",
      caption = "Notes: Thresholds are 50%, 60%, and 70% of the within-sample annual weighted median equivalised disposable income. Shaded bands show 95% confidence intervals."
    ) +
    theme_thesis()
}

fit_probit_model <- function(panel, classified, vars, codes, reference_year, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight

  model_data <- panel %>%
    filter(.data[[year]] == reference_year) %>%
    left_join(classified %>% select(all_of(id), current_poor), by = id) %>%
    mutate(
      poor_current = as.integer(current_poor),
      # Normalising weights leaves weighted point estimates unchanged but makes
      # iterative maximum-likelihood fitting numerically more stable.
      model_weight = .data[[weight]] / mean(.data[[weight]], na.rm = TRUE),
      female = as.integer(.data[[vars$sex]] == codes$sex[["female"]]),
      age = .data[[vars$age]],
      age_sq = age^2,
      education = factor(case_when(
        .data[[vars$education]] %in% c(0, 1, 2) ~ "Primary or below",
        .data[[vars$education]] %in% c(3, 4, 5) ~ "Secondary",
        .data[[vars$education]] %in% c(6, 7, 8) ~ "Tertiary",
        TRUE ~ NA_character_
      )),
      labour_status = factor(case_when(
        .data[[vars$labour_status]] %in% c(1, 2) ~ "Employee",
        .data[[vars$labour_status]] %in% c(3, 4) ~ "Self-employed",
        .data[[vars$labour_status]] == 7 ~ "Retired",
        .data[[vars$labour_status]] %in% c(5, 6, 8, 9, 10) ~ "Inactive/other",
        TRUE ~ NA_character_
      )),
      informal_proxy = factor(case_when(
        .data[[vars$social_security]] %in% codes$likely_informal_social_security_values ~ "Not registered",
        !is.na(.data[[vars$social_security]]) ~ "Registered/other",
        TRUE ~ NA_character_
      )),
      region = factor(.data[[vars$nuts2]]),
      household_type = factor(household_type_derived)
    ) %>%
    select(
      poor_current, female, age, age_sq, education, labour_status, informal_proxy,
      household_type, hh_children_u14, hh_earners_proxy, region, model_weight
    ) %>%
    drop_na()

  # Probit is appropriate because the dependent variable is binary and the model
  # constrains predicted probabilities to [0, 1]. Coefficients are latent-index
  # effects; average marginal effects are the main probability-scale quantities.
  full_formula <- poor_current ~ female + age + age_sq + education + labour_status +
    informal_proxy + household_type + hh_children_u14 + hh_earners_proxy + region

  fallback_formula <- poor_current ~ female + age + age_sq + education + labour_status +
    informal_proxy + household_type + hh_children_u14 + hh_earners_proxy

  fit_once <- function(formula) {
    suppressWarnings(glm(
      formula,
      data = model_data,
      family = quasibinomial(link = "probit"),
      weights = model_weight,
      control = glm.control(maxit = 100)
    ))
  }

  fit <- fit_once(full_formula)
  used_formula <- full_formula
  model_note <- "Full specification includes NUTS-2 region fixed effects."

  # Full regional specifications can be fragile in short panels when some
  # regions have sparse cells or near-perfect prediction. If that happens, the
  # script keeps the analysis reproducible by falling back to the transparent
  # parsimonious specification and records the reason in the model notes.
  unstable <- !isTRUE(fit$converged) || any(!is.finite(coef(fit))) ||
    any(abs(coef(fit)) > 20, na.rm = TRUE)

  if (unstable) {
    fit <- fit_once(fallback_formula)
    used_formula <- fallback_formula
    model_note <- paste(
      "The full NUTS-2 specification showed non-convergence or separation.",
      "Reported model omits region fixed effects; estimate regional models",
      "separately or collapse regions as a robustness check."
    )
  }

  if (!isTRUE(fit$converged)) {
    warning(
      "Fallback probit still did not converge. Inspect model_data for separation ",
      "or simplify categorical predictors further.",
      call. = FALSE
    )
  }

  robust_vcov <- sandwich::vcovHC(fit, type = "HC1")
  coef_table <- broom::tidy(fit) %>%
    mutate(
      robust_se = sqrt(diag(robust_vcov)),
      robust_z = estimate / robust_se,
      robust_p = 2 * pnorm(abs(robust_z), lower.tail = FALSE),
      conf_low = estimate - 1.96 * robust_se,
      conf_high = estimate + 1.96 * robust_se,
      model_note = model_note
    )

  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    marginal_effects <- marginaleffects::avg_slopes(fit, vcov = robust_vcov) %>%
      as_tibble() %>%
      mutate(model_note = model_note)
  } else {
    marginal_effects <- tibble(
      note = "Install the marginaleffects package to compute average marginal effects.",
      model_note = model_note
    )
  }

  model_diagnostics <- tibble(
    n_model = nrow(model_data),
    weighted_n_normalised = sum(model_data$model_weight, na.rm = TRUE),
    converged = isTRUE(fit$converged),
    max_abs_coefficient = max(abs(coef(fit)), na.rm = TRUE),
    formula = paste(deparse(used_formula), collapse = " "),
    note = model_note
  )

  list(model_data = model_data, fit = fit, robust_vcov = robust_vcov,
       coefficients = coef_table, marginal_effects = marginal_effects,
       diagnostics = model_diagnostics)
}
