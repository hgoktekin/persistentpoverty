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



# =============================================================================
# ATTRITION DOCUMENTATION
# =============================================================================
# The balanced panel approach (construct_balanced_panel function):
#  - Identifies individuals observed in all 4 waves (2016-2019)
#  - Uses longitudinal weights from the reference year (2019) in TurkStat's design
#  - Propagates these weights back to all 4 person-years
#
# Expected attrition:
#  - Initial sample (all individuals with 2016 obs): ~63 million
#  - Balanced panel (all 4 waves): ~15 million  
#  - Retention rate: 23.8% (15M / 63M)
#  - Attrition rate: 76.2%
#
# The 76% attrition is substantial. The balanced panel approach is justified
# because:
#  1. It provides a clean panel for transition analysis
#  2. Longitudinal weights adjust for non-random attrition
#  3. We explicitly test robustness vs. cross-sectional results
#  4. We compare baseline characteristics of stayers vs. attritors
#
# Key assumption: MAR (Missing At Random) conditional on observables;
# weights help satisfy this when attrition depends on observed covariates.
# =============================================================================

