# =============================================================================
# Configuration for SILC Türkiye poverty dynamics analysis, 2016-2019
#
# Edit only this file when variable names differ from the defaults below.
# The empirical scripts source this configuration so that the analysis remains
# reproducible and auditable.
# =============================================================================
# To run all the project, run the code below : 
# source("/Users/haticegoktekin/Documents/Codex/WP-PERPOV/R/02_run_silc_poverty_analysis.R")

project <- list(
  data_path = "panel_16_19.dta",
  out_dir = "outputs",
  table_dir = "tables",
  figure_dir = "figures",
  model_dir = "models",
  panel_years = 2016:2019,
  reference_year = 2019,
  thresholds = c(50, 60, 70),
  main_threshold = 60
)

# Core SILC variable map.
# These names match the current Stata file in this workspace. If your cleaned
# extract uses English names, change the right-hand side only.
vars <- list(
  person_id = "fkimlik",
  household_id = "hkimlik",
  year = "hk010",
  household_income = "hg110",
  longitudinal_weight = "fk060_4",
  age = "fk070",
  household_reference_person = "fk095",
  sex = "fk090",
  education = "fe030",
  labour_status = "fi010",
  social_security = "fi020")

# Optional recoding choices. These are deliberately conservative and transparent.
# Verify the value labels in the TurkStat codebook before final publication.
codes <- list(
  sex = c(male = 1, female = 2),
  employed_labour_status = c(1, 2, 3, 4),
  likely_informal_social_security_values = c(2)
)

