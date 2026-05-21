# =============================================================================
# Configuration for SILC Türkiye poverty dynamics analysis, 2016-2019
#
# Edit only this file when variable names differ from the defaults below.
# The empirical scripts source this configuration so that the analysis remains
# reproducible and auditable.
# =============================================================================
# To run all the project, run the code below : 
# source("/Users/haticegoktekin/Documents/Codex/WP-PERPOV/R/02_run_silc_poverty_analysis.R")
setwd("/Users/haticegoktekin/Desktop/phd application/lisans tez/son denemeler/R")

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
  social_security = "fi190")
# household_typo = "hb050" # misclassified as

# ANNUAL INCOME
# Household income components available in the TR-SILC panel extract.
income_components <- c(
  fg010 = "Employee cash income",
  fg020 = "Employee non-cash income",
  fg030 = "Entrepreneur income",
  fg040 = "Non-cash entrepreneur income",
  
  fg070 = "Unemployment benefit",
  fg080 = "Retirement income",
  fg085 = "Retirement benefits",
  fg090 = "Pension for widow and orphan",
  fg100 = "Sickness/disability benefits",
  fg110 = "Disability, veteran, and disability retirement benefits",
  fg120 = "Education allowances",
  
  fg140 = "Total disposable income of the individual",
 
  hg010 = "Annual rent of the house",
  
  hg020 = "Income earned by children (under 15 years old)",
  
  hg030n = "Cash child support",
  hg030a = "Non-cash child support",
  hg040 = "Housing allowances",
  hg050n = "Cash social assistance",
  hg050a = "Non-cash social assistance",
  
  hg060n = "Cash assistance received from another person or household",
  hg060a = "Non-cash assistance received from another person or household",
  hg065n = "Alimony (nafaka)",
  
  hg070 = "Rental income (gayrimenkul)",
  hg080 = "Capital income (menkul)",
  
  hg090n = "Cash assistance given to another person or household",
  hg090a = "Non-cash assistance to another person or household",
  hg095 = "Paid alimony",
  
  hg100 = "Taxes paid",
  hg110 = "Total household disposable income"
)

# of household 
income_categories <- list(
  "Wage income" = c("fg010", "fg020", "hg020"),
  "Entrepreneurial income" = c("fg030", "fg040"),
  "Capital income" = c("hg070", "hg080"), # hg010 is paid rent (not included)
  "Social transfers" = c("fg070", "fg080", "fg085", "fg090", "fg100", "fg110",
                         "fg120","hg030n", "hg030a","hg040", "hg050n", "hg050a"),
  "Private transfers" = c("hg060n", "hg060a", "hg065n", 
                          "hg090n", "hg090a","hg095")
)

# Optional recoding choices. These are deliberately conservative and transparent.
# Verify the value labels in the TurkStat codebook before final publication.
codes <- list(
  sex = c(male = 1, female = 2),
  employed_labour_status = c(1, 2, 3, 4, 5),
  likely_informal_social_security_values = c(2) )

