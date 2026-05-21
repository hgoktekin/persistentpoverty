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
#setwd("/Users/haticegoktekin/Desktop/phd application/lisans tez/son denemeler/R")
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

# --- NEW HELPER: Significance Stars (Task 7) ---
add_significance_stars <- function(p) {
  p <- as.numeric(p)
  case_when(
    p < 0.01  ~ "***",
    p < 0.05  ~ "**",
    p < 0.1   ~ "*",
    TRUE      ~ ""
  )
}

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

# =============================================================================
#  Calculate attrition rate and compare baseline characteristics
# =============================================================================
calculate_attrition_rate <- function(raw, vars, codes, panel_years, reference_year) {
  year_var <- vars$year
  id_var   <- vars$person_id
  weight_var <- vars$longitudinal_weight

  by_year <- raw %>%
    group_by(.data[[year_var]]) %>%
    summarise(
      n_individuals_unweighted = n_distinct(.data[[id_var]]),
      weighted_population = sum(.data[[weight_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(.data[[year_var]])

  balanced_ids <- raw %>%
    group_by(.data[[id_var]]) %>%
    filter(n_distinct(.data[[year_var]]) == length(panel_years)) %>%
    distinct(.data[[id_var]]) %>%
    pull(.data[[id_var]])

  n_balanced <- length(balanced_ids)
  n_first    <- by_year %>%
    filter(.data[[year_var]] == min(panel_years)) %>%
    pull(n_individuals_unweighted)
  attrition_rate <- 1 - (n_balanced / n_first)

  # --- Expanded baseline comparison (first wave) ---
  baseline <- raw %>%
    filter(.data[[year_var]] == min(panel_years)) %>%
    mutate(
      stayer = .data[[id_var]] %in% balanced_ids,
      stayer_label = if_else(stayer, "Stayer", "Attritor"),
      female = if_else(.data[[vars$sex]] == codes$sex[["female"]], 1L, 0L, missing = NA_integer_),
      informal = if_else(
        .data[[vars$social_security]] %in% codes$likely_informal_social_security_values,
        1L, 0L, missing = 0L),
      edu_primary = if_else(.data[[vars$education]] %in% c(0, 1, 2), 1L, 0L, missing = NA_integer_),
      edu_secondary = if_else(.data[[vars$education]] %in% c(3, 4, 5), 1L, 0L, missing = NA_integer_),
      edu_tertiary = if_else(.data[[vars$education]] %in% c(6, 7, 8), 1L, 0L, missing = NA_integer_),
      employed = if_else(.data[[vars$labour_status]] %in% codes$employed_labour_status, 1L, 0L, missing = NA_integer_)
    )

  w <- baseline[[weight_var]]
  fmt_v <- function(x) sprintf("%.3f", x)
  fmt_p <- function(p) {
    p <- as.numeric(p)
    stars <- add_significance_stars(p)
    paste0(ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)), stars)
  }

  # Helper: weighted mean by group
  wmean_by <- function(var) {
    baseline %>%
      filter(!is.na(.data[[var]])) %>%
      group_by(stayer_label) %>%
      summarise(m = weighted_mean_safe(.data[[var]], .data[[weight_var]]), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = stayer_label, values_from = m)
  }

  # Helper: weighted proportion by group
  wprop_by <- function(var) {
    baseline %>%
      filter(!is.na(.data[[var]])) %>%
      group_by(stayer_label) %>%
      summarise(m = weighted_mean_safe(.data[[var]], .data[[weight_var]]), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = stayer_label, values_from = m)
  }

  # t-test p-value (survey-weighted)
  wt_p <- function(var) {
    df <- baseline %>% filter(!is.na(.data[[var]]))
    des <- survey::svydesign(ids = ~1, weights = as.formula(paste0("~", weight_var)), data = df)
    tt <- tryCatch(
      survey::svyttest(as.formula(paste0(var, " ~ stayer")), design = des),
      error = function(e) NULL)
    if (is.null(tt)) return(NA_real_)
    tt$p.value
  }

  build_row <- function(label, var) {
    vals <- wmean_by(var)
    p <- wt_p(var)
    tibble(
      Characteristic = label,
      Stayer = fmt_v(vals$Stayer),
      Attritor = fmt_v(vals$Attritor),
      Difference = fmt_v(vals$Stayer - vals$Attritor),
      p_value = fmt_p(p)
    )
  }

  n_row <- tibble(
    Characteristic = "N (unweighted)",
    Stayer = as.character(sum(baseline$stayer)),
    Attritor = as.character(sum(!baseline$stayer)),
    Difference = "",
    p_value = ""
  )

  attrition_table <- bind_rows(
    n_row,
    build_row("Age, mean", vars$age),
    build_row("Female, proportion", "female"),
    build_row("Household income, mean", vars$household_income),
    build_row("Primary or below, proportion", "edu_primary"),
    build_row("Secondary, proportion", "edu_secondary"),
    build_row("Tertiary, proportion", "edu_tertiary"),
    build_row("Employed, proportion", "employed"),
    build_row("Informal (not registered), proportion", "informal")
  )

  list(
    by_year = by_year,
    balanced_panel_n = n_balanced,
    first_year_n = n_first,
    attrition_rate_percent = attrition_rate * 100,
    retention_rate_percent = (1 - attrition_rate) * 100,
    attrition_table = attrition_table)
}

construct_balanced_panel <- function(raw, vars, panel_years, reference_year) {
  require_vars(raw, vars[c(
    "person_id", "household_id", "year", "household_income",
    "longitudinal_weight", "age", "household_reference_person")])

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
      call. = FALSE)
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
        TRUE ~ NA_real_)
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
  id     <- vars$person_id
  hh     <- vars$household_id
  year   <- vars$year
  age    <- vars$age
  labour <- vars$labour_status
  ss     <- vars$social_security
  sex_v  <- vars$sex

  panel %>%
    mutate(
      is_child_u14      = .data[[age]] < 14,
      is_child_u18      = .data[[age]] < 18,
      is_elderly_65plus = .data[[age]] >= 65,
      is_adult_14_65    = .data[[age]] >= 14 & .data[[age]] < 65,
      is_earner_proxy   = if_else(
        .data[[labour]] %in% codes$employed_labour_status, 1L, 0L, missing = 0L),
      # Binary person-level variables (0/1)
      informal_status = if_else(
        .data[[ss]] %in% codes$likely_informal_social_security_values,
        1L, 0L, missing = 0L),
      female = if_else(
        .data[[sex_v]] == codes$sex[["female"]], 1L, 0L, missing = NA_integer_)
    ) %>%
    group_by(.data[[hh]], .data[[year]]) %>%
    mutate(
      hh_size            = n_distinct(.data[[id]]),
      hh_children_u14    = sum(is_child_u14, na.rm = TRUE),
      hh_children_u18    = sum(is_child_u18, na.rm = TRUE),
      hh_elderly_65plus  = sum(is_elderly_65plus, na.rm = TRUE),
      hh_earners_proxy   = sum(is_earner_proxy, na.rm = TRUE),
      hh_formal_earners  = sum(is_formal_earner, na.rm = TRUE),
      hh_informal_earners = sum(is_informal_earner, na.rm = TRUE),
      hh_adults_14_65             = sum(is_adult_14_65, na.rm = TRUE),
      hh_adults_14_65_working     = sum(is_earner_proxy == 1L & is_adult_14_65, na.rm = TRUE),
      hh_adults_14_65_not_working = hh_adults_14_65 - hh_adults_14_65_working,
      # Standard dependency ratio
      dependency_ratio = if_else(
        hh_earners_proxy > 0,
        (hh_size - hh_earners_proxy) / hh_earners_proxy, NA_real_),
      # Modified (OECD-style) dependency ratio — FIXED
      # 0.3 × children <14 + 1.0 × adults 14-65 not working + 0.5 × elderly 65+
      # denominator: adults 14-65 working (not all earners)
      dependency_ratio_oecd = if_else(
        hh_adults_14_65_working > 0,
        (0.3 * hh_children_u14 +
         1.0 * hh_adults_14_65_not_working +
         0.5 * hh_elderly_65plus) / hh_adults_14_65_working,
        NA_real_),
      # Log of modified dependency ratio
      log_dependency_ratio_oecd = log(dependency_ratio_oecd + 0.001),
      household_type_derived = case_when(
        hh_size == 1 ~ "Single person",
        hh_children_u18 > 0 & hh_elderly_65plus > 0 ~ "Children and elderly present",
        hh_children_u18 > 0 ~ "Children present",
        hh_elderly_65plus > 0 ~ "Elderly present",
        TRUE ~ "Adults only")
    ) %>%
    ungroup()
}

add_employment_stability <- function(panel, vars) {
  id   <- vars$person_id
  year <- vars$year
  panel %>%
    arrange(.data[[id]], .data[[year]]) %>%
    group_by(.data[[id]]) %>%
    mutate(
      hh_earners_lag = lag(hh_earners_proxy, order_by = .data[[year]]),
      earner_loss = if_else(
        !is.na(hh_earners_lag) & hh_earners_proxy < hh_earners_lag, 1L, 0L)
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
          poverty_line_70 = ~ .x * 0.70),
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
      poor_once = n_poor_4yr %in% c(1,2),
      poor_multiple = n_poor_4yr == 3 & current_poor == 0 , #cant be current poor 
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
        TRUE ~ NA_character_)
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
  total_n <- nrow(classified)
  total_w <- sum(classified[[weight]], na.rm = TRUE)

  mutual <- classified %>%
    count(poverty_group, name = "unweighted_n") %>%
    left_join(
      classified %>%
        count(poverty_group, wt = .data[[weight]], name = "weighted_n"),
      by = "poverty_group") %>%
    mutate(
      table_family = "Mutually exclusive four-year typology",
      unweighted_share = unweighted_n / total_n,
      weighted_share = weighted_n / total_w
    ) %>%
    rename(category = poverty_group)

  overlapping <- tibble(
    table_family = "Overlapping current-year flags",
    category = c("Current poor (2019)", "Persistent poor"),
    unweighted_n = c(
      sum(classified$current_poor, na.rm = TRUE),
      sum(classified$persistent_poor, na.rm = TRUE)),
    weighted_n = c(
      sum(classified[[weight]][classified$current_poor], na.rm = TRUE),
      sum(classified[[weight]][classified$persistent_poor], na.rm = TRUE)
    )
  ) %>%
    mutate(
      unweighted_share = unweighted_n / total_n,
      weighted_share = weighted_n / total_w)
  bind_rows(mutual, overlapping)
}

poverty_episode_summary <- function(status) {
  status <- as.integer(status)
  if (all(is.na(status))) {
    return(tibble(n_episodes = NA_integer_, max_spell_duration = NA_integer_, 
                  mean_spell_duration = NA_real_))}
  status[is.na(status)] <- 0L
  starts <- which(status == 1L & dplyr::lag(status, default = 0L) == 0L)
  ends <- which(status == 1L & dplyr::lead(status, default = 0L) == 0L)
  durations <- ends - starts + 1L
  tibble(
    n_episodes = length(durations),
    max_spell_duration = ifelse(length(durations) == 0, 0L, max(durations)),
    mean_spell_duration = ifelse(length(durations) == 0, 0, mean(durations)))
}

make_poverty_duration_table <- function(panel, classified, vars, panel_years, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  poor <- paste0("poor_", threshold)
  total_n <- nrow(classified)
  total_w <- sum(classified[[weight]], na.rm = TRUE)

  episodes <- panel %>%
    arrange(.data[[id]], .data[[year]]) %>%
    group_by(.data[[id]]) %>%
    summarise(
      poverty_episode_summary(.data[[poor]]),
      panel_weight = first(.data[[weight]]),
      .groups = "drop")

  duration <- classified %>%
    mutate(duration_years = n_poor_4yr) %>%
    count(duration_years, name = "unweighted_n") %>%
    left_join(
      classified %>%
        mutate(duration_years = n_poor_4yr) %>%
        count(duration_years, wt = .data[[weight]], name = "weighted_n"),
      by = "duration_years"
    ) %>%
    mutate(
      unweighted_share = unweighted_n / total_n,
      weighted_share = weighted_n / total_w
    )
  list(duration = duration, episodes = episodes)}

make_transition_matrices <- function(panel, vars, panel_years, 
                                     threshold, codes = NULL) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  poor <- paste0("poor_", threshold)
  social_security <- vars$social_security
  
  wide <- panel %>%
    mutate(
      formal_status = case_when(
        is.null(codes) | is.null(codes$likely_informal_social_security_values) ~ NA_character_,
        is.na(.data[[social_security]]) ~ NA_character_,
        .data[[social_security]] %in% codes$likely_informal_social_security_values ~ "Informal",
        TRUE ~ "Formal"
      )
    ) %>%
    transmute(
      !!id := .data[[id]],
      !!year := .data[[year]],
      !!weight := .data[[weight]],
      poor_status = .data[[poor]],
      formal_status = formal_status
    ) %>%
    pivot_wider(
      names_from = all_of(year),
      values_from = c(poor_status, formal_status),
      names_sep = "_"
    )
  
  purrr::map_dfr(seq_along(panel_years[-length(panel_years)]), function(i) {
    y0 <- panel_years[i]
    y1 <- panel_years[i + 1]
    from <- paste0("poor_status_", y0)
    to <- paste0("poor_status_", y1)
    status_at_t <- paste0("formal_status_", y0)
    
    
    wide %>%
      filter(!is.na(.data[[from]]), !is.na(.data[[to]]), !is.na(.data[[status_at_t]])) %>%
      mutate(
        from_status = if_else(.data[[from]] == 1, "Poor", "Non-poor"),
        to_status = if_else(.data[[to]] == 1, "Poor", "Non-poor"),
        formal_status = .data[[status_at_t]] 
        ) %>%
      count(formal_status, from_status, to_status, wt = .data[[weight]], name = "weighted_n") %>%
      group_by(formal_status, from_status) %>%
      mutate(row_probability = weighted_n / sum(weighted_n)) %>%
      ungroup() %>%
      mutate(
        transition = paste(y0, y1, sep = "-"),
        transition_type = case_when(
          from_status == "Non-poor" & to_status == "Poor" ~ "Entry",
          from_status == "Poor" & to_status == "Non-poor" ~ "Exit",
          from_status == "Poor" & to_status == "Poor" ~ "Poverty persistence",
          TRUE ~ "Non-poor persistence")) } ) }

make_profile_table <- function(panel, classified, vars, codes, reference_year, threshold) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  fmt_n <- function(x) {formatC(round(x), format = "f", digits = 0, big.mark = ",")}
  fmt_num <- function(x) {ifelse(is.na(x), "", 
                          formatC(x, format = "f", digits = 2, big.mark = ","))}
  
  fmt_pct <- function(x) {ifelse(is.na(x), "", sprintf("%.2f%%", 100 * x))}
  fmt_p <- function(p) {
    p <- as.numeric(p)[1]
    stars <- add_significance_stars(p)
    case_when(is.na(p) ~ "-", p < 0.001 ~ paste0("<0.001", stars),
             TRUE ~ paste0(sprintf("%.3f", p), stars))}

  w_sd <- function(x, w) { ok <- !is.na(x) & !is.na(w) & w > 0
    if (sum(ok) < 2) return(NA_real_)
    mu <- weighted_mean_safe(x[ok], w[ok])
    sqrt(sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])) }

  cell_n_pct <- function(data, var, level) {
    x <- data[[var]]
    ok <- !is.na(x)
    denom <- sum(ok, na.rm = TRUE)
    if (denom == 0) return("")
    n <- sum(x[ok] == level, na.rm = TRUE)
    pct <- n / denom
    paste0(fmt_n(n), " (", fmt_pct(pct), ")")}

  cell_mean_sd <- function(data, var) {
    x <- data[[var]]
    w <- data[[weight]]
    paste0(fmt_n(weighted_mean_safe(x, w)), " (", fmt_n(w_sd(x, w)), ")") }

  profile <- panel %>%
    filter(.data[[year]] == reference_year) %>%
    left_join( classified %>%
        select(all_of(id), poverty_group, n_poor_4yr, current_poor, persistent_poor),
      by = id ) %>%
    mutate(
      # Used only for p-values. This is mutually exclusive, unlike the
      # descriptive "Currently poor (2019)" column.
      profile_spell_group = case_when(
        n_poor_4yr == 0 ~ "Never poor",
        n_poor_4yr %in% c(1, 2) & !persistent_poor ~ "Transient poor",
        persistent_poor ~ "Persistent poor",
        n_poor_4yr == 3 & !persistent_poor ~ "Frequently poor"),
      poverty_group_for_test = factor(
        profile_spell_group,
        levels = c("Never poor", "Transient poor", "Frequently poor", "Persistent poor")),
      sex_recoded = factor(
        case_when(
          .data[[vars$sex]] == codes$sex[["female"]] ~ "Female",
          .data[[vars$sex]] == codes$sex[["male"]] ~ "Male",
          TRUE ~ NA_character_),
        levels = c("Female", "Male")),
      age_group = cut(
        .data[[vars$age]],
        breaks = c(-Inf, 17, 24, 34, 44, 54, 64, Inf),
        labels = c("0-17", "18-24", "25-34", "35-44", "45-54", "55-64", "65+")),
      education_recoded = factor(
        case_when(
          .data[[vars$education]] %in% c(0, 1, 2) ~ "Primary or below",
          .data[[vars$education]] %in% c(3, 4, 5) ~ "Secondary",
          .data[[vars$education]] == 6 ~ "Tertiary",
          TRUE ~ NA_character_),
        levels = c("Primary or below", "Secondary", "Tertiary")),
      labour_recoded = factor(
        case_when(
          .data[[vars$labour_status]] %in% c(1, 2)          ~ "Employee",
          .data[[vars$labour_status]] %in% c(3, 4)          ~ "Self-employed",
          .data[[vars$labour_status]] == 5                   ~ "Unemployed",
          .data[[vars$labour_status]] == 7                   ~ "Retired",
          .data[[vars$labour_status]] %in% c(6, 8, 9, 10) ~ "Inactive",
          TRUE ~ NA_character_),
        levels = c("Employee", "Self-employed", "Unemplyoed", "Retired", "Inactive")
      ),
      informal_status_label = factor(
        if_else(informal_status == 1L, "Informal", "Formal"),
        levels = c("Formal", "Informal"))
    )

  weighted_chisq_p <- function(var) {
    df <- profile %>%
      filter(!is.na(.data[[var]]), !is.na(poverty_group_for_test))
    if (n_distinct(df[[var]]) < 2 || n_distinct(df$poverty_group_for_test) < 2) {
      return(NA_real_) }
    des <- survey::svydesign(
      ids = ~1, weights = as.formula(paste0("~", weight)), data = df )
    out <- tryCatch(
      survey::svychisq(as.formula(paste0("~ poverty_group_for_test + ", var)),
                       design = des, statistic = "F"),
      error = function(e) NULL )
    if (is.null(out)) NA_real_ else out$p.value }

  weighted_kw_p <- function(var) {
    df <- profile %>%
      filter(!is.na(.data[[var]]), !is.na(poverty_group_for_test))
    if (n_distinct(df$poverty_group_for_test) < 2) {
      return(NA_real_) }
    des <- survey::svydesign(
      ids = ~1, weights = as.formula(paste0("~", weight)), data = df )
    out <- tryCatch(
      survey::svyranktest(as.formula(paste0(var, " ~ poverty_group_for_test")),
                          design = des, test = "KruskalWallis"),
      error = function(e) NULL )
    if (is.null(out)) NA_real_ else out$p.value }

  groups <- list(
    overall = profile,
    never_poor = profile %>% filter(poverty_group == "Never poor"),
    transient_poor = profile %>% filter(profile_spell_group == "Transient poor"),
    currently_poor = profile %>% filter(current_poor),
    persistent_poor = profile %>% filter(persistent_poor) )

  headers <- tibble(
    col = names(groups),
    label = c(
      paste0("**Overall**  \nN = ", fmt_n(nrow(groups$overall))),
      paste0("**Never poor**  \nN = ", fmt_n(nrow(groups$never_poor))),
      paste0("**Transient poor**  \nN = ", fmt_n(nrow(groups$transient_poor))),
      paste0("**Currently poor**  \n(2019)  \nN = ", fmt_n(nrow(groups$currently_poor))),
      paste0("**Persistent poor**  \n*(subset of current)*  \nN = ", fmt_n(nrow(groups$persistent_poor)))
    ) )

  section_row <- function(label) {
    tibble(
      row_type = "section",
      variable = label,
      p_value = "",
      overall = "",
      never_poor = "",
      transient_poor = "",
      currently_poor = "",
      persistent_poor = "" ) }

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
        )  } ) ) }

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
    )  }

  out <- bind_rows(
    categorical_rows("Age group", "age_group", levels(profile$age_group)),
    categorical_rows("Sex", "sex_recoded", levels(profile$sex_recoded)),
    categorical_rows("Education", "education_recoded", levels(profile$education_recoded)),
    categorical_rows("Employment status", "labour_recoded", levels(profile$labour_recoded)),
    categorical_rows("Social Security Status", "informal_status_label", c("Formal", "Informal")),
  
    categorical_rows("Household type", "household_type_derived", levels(factor(profile$household_type_derived))),
    continuous_row("Household size, mean (SD); min-max", "hh_size"),
    continuous_row("Children under 14, mean (SD); min-max", "hh_children_u14"),
    continuous_row("Earners in household, mean (SD); min-max", "hh_earners_proxy"),
    continuous_row("Dependency ratio, mean (SD)", "dependency_ratio"),
    continuous_row("Modified dep. ratio (OECD), mean (SD)", "dependency_ratio_oecd"),
  
    continuous_row("Eq. income (TL), mean (SD)", "eq_income"),
    continuous_row("Eq. income (TL), mean (SD); min-max", "eq_income"))
    

  attr(out, "headers") <- headers
  out
}

make_income_composition_table <- function(panel, classified, vars, income_components, reference_year) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  available_components <- intersect(names(income_components), names(panel))
  
  if (length(available_components) == 0) {
    return(tibble(
      poverty_profile = character(),
      income_variable = character(),
      income_type = character(),
      unweighted_n = integer(),
      unweighted_positive_n = integer(),
      mean_income = numeric(),
      income_share = numeric()
    ))
  }
  
  component_vars <- setdiff(available_components, "fg140")
  denominator_var <- if ("fg140" %in% available_components) "fg140" else NA_character_
  
  profile <- panel %>%
    filter(.data[[year]] == reference_year) %>%
    left_join(
      classified %>%
        select(all_of(id), poverty_group, n_poor_4yr, current_poor, persistent_poor),
      by = id
    ) %>%
    mutate(
      profile_spell_group = case_when(
        n_poor_4yr == 0 ~ "Never poor",
        n_poor_4yr %in% c(1, 2) & !persistent_poor ~ "Transient poor",
        persistent_poor ~ "Persistent poor",
        TRUE ~ "Frequently poor"
      )
    )
  
  groups <- list(
    "Overall" = profile,
    "Never poor" = profile %>% filter(poverty_group == "Never poor"),
    "Transient poor" = profile %>% filter(profile_spell_group == "Transient poor"),
    "Currently poor (2019)" = profile %>% filter(current_poor),
    "Persistent poor" = profile %>% filter(persistent_poor)
  )
  
  purrr::imap_dfr(groups, function(data, group_name) {
    w <- data[[weight]]
    denom_total <- if (!is.na(denominator_var)) {
      sum(replace_na(as.numeric(data[[denominator_var]]), 0) * w, na.rm = TRUE)
    } else if (length(component_vars) > 0) {
      component_matrix <- as.data.frame(data[component_vars])
      component_total <- rowSums(dplyr::mutate(component_matrix, across(everything(), ~replace_na(as.numeric(.x), 0))))
      sum(component_total * w, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    purrr::map_dfr(available_components, function(v) {
      x <- replace_na(as.numeric(data[[v]]), 0)
      weighted_total <- sum(x * w, na.rm = TRUE)
      tibble(
        poverty_profile = group_name,
        income_variable = v,
        income_type = unname(income_components[[v]]),
        unweighted_n = nrow(data),
        unweighted_positive_n = sum(x > 0, na.rm = TRUE),
        mean_income = weighted_mean_safe(x, w),
        weighted_income_total = weighted_total,
        income_share = if_else(!is.na(denom_total) & denom_total != 0, weighted_total / denom_total, NA_real_)
      )
    })
  }) %>%
    mutate(
      income_role = if_else(income_variable == denominator_var, "Total income", "Income component"),
      income_type = factor(income_type, levels = unname(income_components[available_components])),
      poverty_profile = factor(
        poverty_profile,
        levels = c("Overall", "Never poor", "Transient poor", "Currently poor (2019)", "Persistent poor")
      )
    ) %>%
    arrange(poverty_profile, income_role == "Total income", income_variable)
}

make_income_category_table <- function(panel, classified, vars, income_categories, reference_year) {
  id <- vars$person_id
  year <- vars$year
  weight <- vars$longitudinal_weight
  
  profile <- panel %>%
    filter(.data[[year]] == reference_year) %>%
    left_join(
      classified %>%
        select(all_of(id), poverty_group, n_poor_4yr, persistent_poor),
      by = id
    ) %>%
    mutate(
      profile_spell_group = case_when(
        n_poor_4yr == 0 ~ "Never poor",
        n_poor_4yr %in% c(1, 2) & !persistent_poor ~ "Transient poor",
        persistent_poor ~ "Persistent poor",
        TRUE ~ "Frequently poor"
      )
    )
  
  for (category in names(income_categories)) {
    available_vars <- intersect(income_categories[[category]], names(profile))
    profile[[category]] <- if (length(available_vars) == 0) {
      0
    } else {
      rowSums(as.data.frame(profile[available_vars]) %>%
                mutate(across(everything(), ~replace_na(as.numeric(.x), 0))))
    }
  }
  
  category_names <- names(income_categories)
  profile <- profile %>%
    mutate(total_categorised_income = rowSums(across(all_of(category_names)), na.rm = TRUE))
  
  groups <- list(
    "Overall" = profile,
    "Never poor" = profile %>% filter(poverty_group == "Never poor"),
    "Transient poor" = profile %>% filter(profile_spell_group == "Transient poor"),
    "Persistent poor" = profile %>% filter(persistent_poor)
  )
  
  purrr::imap_dfr(groups, function(data, group_name) {
    w <- data[[weight]]
    denominator <- sum(data$total_categorised_income * w, na.rm = TRUE)
    
    purrr::map_dfr(category_names, function(category) {
      x <- data[[category]]
      weighted_total <- sum(x * w, na.rm = TRUE)
      tibble(
        poverty_profile = group_name,
        income_category = category,
        included_variables = paste(intersect(income_categories[[category]], names(panel)), collapse = ", "),
        unweighted_n = nrow(data),
        unweighted_positive_n = sum(x > 0, na.rm = TRUE),
        mean_income = weighted_mean_safe(x, w),
        weighted_income_total = weighted_total,
        income_share = if_else(denominator != 0, weighted_total / denominator, NA_real_)
      )
    })
  }) %>%
    mutate(
      poverty_profile = factor(
        poverty_profile,
        levels = c("Overall", "Never poor", "Transient poor", "Persistent poor")
      ),
      income_category = factor(income_category, levels = category_names)
    ) %>%
    arrange(income_category, poverty_profile)
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



