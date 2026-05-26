/*=============================================================================
  00_master.do — Stata replication of the R poverty dynamics pipeline
  
  Persistent Poverty Analysis using TR-SILC 2016-2019
  
  This do-file replicates the full data preparation, variable construction,
  and descriptive analysis implemented across the R scripts:
    00_config.R, 01_functions.R, 02_run_silc_poverty_analysis.R,
    03_run_descriptive.R, 05_descriptives.R
  
  Sections:
    0. Configuration (paths, globals, thresholds)
    1. Data loading
    2. Balanced panel construction and weight propagation
    3. Modified OECD equivalence scale and equivalised income
    4. Household context and composition variables
    5. Employment stability (lagged earners, earner loss)
    6. Annual poverty thresholds and poverty status
    7. Poverty spell classification and Eurostat persistent poverty
    8. Variable recoding for descriptive tables and models
    9. Descriptive statistics
   10. Export analysis panel

  Input:  panel_16_19.dta
  Output: outputs/analysis_panel_2016_2019.dta
          outputs/poverty_lines.dta
          tables/*.csv
=============================================================================*/

clear all
set more off
set matsize 10000
cap log close
version 16

* ============================================================================
* SECTION 0: CONFIGURATION
* ============================================================================

* --- Paths (edit these to match your local directory structure) ---------------

global projdir   "."                          // project root
global datafile   "${projdir}/panel_16_19.dta"   // raw SILC panel data
global outdir     "${projdir}/outputs"
global tabledir   "${projdir}/tables"
global figuredir  "${projdir}/figures"
global modeldir   "${projdir}/models"

cap mkdir "${outdir}"
cap mkdir "${tabledir}"
cap mkdir "${figuredir}"
cap mkdir "${modeldir}"

* --- Panel parameters --------------------------------------------------------

global panel_years     "2016 2017 2018 2019"
global reference_year  2019
global first_year      2016
global thresholds      "50 60 70"
global main_threshold  60

* --- Core SILC variable map --------------------------------------------------
*     These names match the current Stata variable names in the data file.
*     If your extract uses different names, change the right-hand side only.

global v_person_id       "fkimlik"
global v_household_id    "hkimlik"
global v_year            "hk010"
global v_household_income "hg110"
global v_longitudinal_wt "fk060_4"
global v_age             "fk070"
global v_hh_ref_person   "fk095"
global v_sex             "fk090"
global v_education       "fe030"
global v_labour_status   "fi010"
global v_social_security "fi190"

* --- Code values -------------------------------------------------------------

global code_male     1
global code_female   2
global codes_employed "1 2 3 4 5"
global code_informal_ss 2    // social security value indicating informal

* --- Income components -------------------------------------------------------
*     (kept as globals for reference; used in income composition tables)

global inc_wage       "fg010 fg020 hg020"
global inc_entrepren  "fg030 fg040"
global inc_capital    "hg070 hg080"
global inc_social_tr  "fg070 fg080 fg085 fg090 fg100 fg110 fg120 hg030n hg030a hg040 hg050n hg050a"
global inc_private_tr "hg060n hg060a hg065n hg090n hg090a hg095"

display as text "=== Configuration loaded ==="
display as text "  Data:       ${datafile}"
display as text "  Panel:      ${first_year}-${reference_year}"
display as text "  Thresholds: ${thresholds}"

* ============================================================================
* SECTION 1: DATA LOADING
* ============================================================================

display as text _n "Step 1: Reading raw SILC panel data"

use "${datafile}", clear
describe, short
summarize ${v_person_id} ${v_year} ${v_household_income} ${v_longitudinal_wt} ${v_age}

* ============================================================================
* SECTION 2: BALANCED PANEL CONSTRUCTION AND WEIGHT PROPAGATION
* ============================================================================

display as text _n "Step 2: Constructing balanced panel and propagating longitudinal weights"

* --- 2a. Identify individuals with valid reference-year longitudinal weight --
*     TurkStat's 4-year longitudinal weight is attached to the reference wave.

preserve
  keep if ${v_year} == ${reference_year} & !missing(${v_longitudinal_wt}) & ${v_longitudinal_wt} > 0
  keep ${v_person_id} ${v_longitudinal_wt}
  rename ${v_longitudinal_wt} panel_weight_reference
  * De-duplicate (keep first occurrence per person)
  bysort ${v_person_id}: keep if _n == 1
  tempfile ref_weights
  save `ref_weights'
restore

* --- 2b. Keep only individuals present in the reference-year weight file -----

merge m:1 ${v_person_id} using `ref_weights', keep(match) nogenerate

* --- 2c. Restrict to panel years ---------------------------------------------

keep if inlist(${v_year}, 2016, 2017, 2018, 2019)

* --- 2d. De-duplicate person-year records ------------------------------------

bysort ${v_person_id} ${v_year}: keep if _n == 1

* --- 2e. Replace original weight with propagated reference-year weight -------

drop ${v_longitudinal_wt}
rename panel_weight_reference ${v_longitudinal_wt}

* --- 2f. Verify panel balance ------------------------------------------------

bysort ${v_person_id}: gen _n_years = _N
tab _n_years

* Assert fully balanced (4 years per person)
assert _n_years == 4
drop _n_years

* Check no missing/non-positive weights
assert !missing(${v_longitudinal_wt}) & ${v_longitudinal_wt} > 0

quietly distinct ${v_person_id}
local n_individuals = r(ndistinct)
local n_person_years = _N
display as text "  Balanced panel: `n_individuals' individuals, `n_person_years' person-years"

* ============================================================================
* SECTION 3: MODIFIED OECD EQUIVALENCE SCALE AND EQUIVALISED INCOME
* ============================================================================

display as text _n "Step 3: Calculating modified OECD equivalence scale and equivalised income"

* --- 3a. Assign OECD member weights ------------------------------------------
*     Reference person = 1.0; other adults (14+) = 0.5; children (<14) = 0.3

gen oecd_member_weight = .
replace oecd_member_weight = 1.0 if ${v_hh_ref_person} == 1
replace oecd_member_weight = 0.5 if ${v_hh_ref_person} != 1 & ${v_age} >= 14
replace oecd_member_weight = 0.3 if ${v_hh_ref_person} != 1 & ${v_age} < 14

* --- 3b. Household equivalence size and child/adult counts -------------------

bysort ${v_household_id} ${v_year}: egen hh_eq_size = total(oecd_member_weight)
bysort ${v_household_id} ${v_year}: egen hh_n_children_u14 = total(${v_age} < 14)
bysort ${v_household_id} ${v_year}: egen hh_n_adults_14plus = total(${v_age} >= 14)

* --- 3c. Equivalised income --------------------------------------------------

gen eq_income = ${v_household_income} / hh_eq_size

* Verification
assert !missing(eq_income) & hh_eq_size > 0

display as text "  Equivalised income computed. Summary:"
summarize eq_income, detail

* ============================================================================
* SECTION 4: HOUSEHOLD CONTEXT AND COMPOSITION VARIABLES
* ============================================================================

display as text _n "Step 4: Constructing household context variables"

* --- 4a. Person-level binary indicators --------------------------------------

gen is_child_u14      = (${v_age} < 14)
gen is_child_u18      = (${v_age} < 18)
gen is_elderly_65plus = (${v_age} >= 65)
gen is_adult_14_65    = (${v_age} >= 14 & ${v_age} < 65)

* Earner proxy: 1 if in employment (codes 1-5)
gen is_earner_proxy = 0
foreach code of numlist $codes_employed {
  replace is_earner_proxy = 1 if ${v_labour_status} == `code'
}

* Formal/informal earner classification
gen is_formal_earner = (is_earner_proxy == 1 & ${v_social_security} != ${code_informal_ss}) ///
                       if !missing(${v_social_security})
replace is_formal_earner = 0 if missing(is_formal_earner)

gen is_informal_earner = (is_earner_proxy == 1 & ${v_social_security} == ${code_informal_ss}) ///
                         if !missing(${v_social_security})
replace is_informal_earner = 0 if missing(is_informal_earner)

* Informal status (binary 0/1 for any person, not just earners)
gen informal_status = (${v_social_security} == ${code_informal_ss}) ///
                      if !missing(${v_social_security})
replace informal_status = 0 if missing(informal_status)

* Female indicator
gen female = (${v_sex} == ${code_female}) if !missing(${v_sex})

* --- 4b. Household-level aggregates ------------------------------------------

bysort ${v_household_id} ${v_year}: egen hh_size            = count(${v_person_id})
bysort ${v_household_id} ${v_year}: egen hh_children_u14    = total(is_child_u14)
bysort ${v_household_id} ${v_year}: egen hh_children_u18    = total(is_child_u18)
bysort ${v_household_id} ${v_year}: egen hh_elderly_65plus  = total(is_elderly_65plus)
bysort ${v_household_id} ${v_year}: egen hh_earners_proxy   = total(is_earner_proxy)
bysort ${v_household_id} ${v_year}: egen hh_formal_earners  = total(is_formal_earner)
bysort ${v_household_id} ${v_year}: egen hh_informal_earners = total(is_informal_earner)
bysort ${v_household_id} ${v_year}: egen hh_adults_14_65    = total(is_adult_14_65)

gen hh_adults_14_65_working = .
bysort ${v_household_id} ${v_year}: egen hh_adults_14_65_working_t = total(is_earner_proxy == 1 & is_adult_14_65 == 1)
replace hh_adults_14_65_working = hh_adults_14_65_working_t
drop hh_adults_14_65_working_t

gen hh_adults_14_65_not_working = hh_adults_14_65 - hh_adults_14_65_working

* --- 4c. Dependency ratios ---------------------------------------------------

* Standard dependency ratio
gen dependency_ratio = (hh_size - hh_earners_proxy) / hh_earners_proxy ///
                       if hh_earners_proxy > 0

* Modified (OECD-style) dependency ratio
* Numerator: 0.3 * children <14 + 1.0 * adults 14-65 not working + 0.5 * elderly 65+
* Denominator: adults 14-65 who are working
gen dependency_ratio_oecd = (0.3 * hh_children_u14 ///
                           + 1.0 * hh_adults_14_65_not_working ///
                           + 0.5 * hh_elderly_65plus) ///
                           / hh_adults_14_65_working ///
                           if hh_adults_14_65_working > 0

* Log of modified dependency ratio
gen log_dependency_ratio_oecd = ln(dependency_ratio_oecd + 0.001) ///
                                if !missing(dependency_ratio_oecd)

* --- 4d. Household type classification ---------------------------------------

gen household_type_derived = ""
replace household_type_derived = "Single person"                if hh_size == 1
replace household_type_derived = "Children and elderly present" if hh_size > 1 & hh_children_u18 > 0 & hh_elderly_65plus > 0
replace household_type_derived = "Children present"             if hh_size > 1 & hh_children_u18 > 0 & hh_elderly_65plus == 0 & household_type_derived == ""
replace household_type_derived = "Elderly present"              if hh_size > 1 & hh_children_u18 == 0 & hh_elderly_65plus > 0 & household_type_derived == ""
replace household_type_derived = "Adults only"                  if household_type_derived == ""

display as text "  Household types:"
tab household_type_derived

* ============================================================================
* SECTION 5: EMPLOYMENT STABILITY (LAGGED EARNERS, EARNER LOSS)
* ============================================================================

display as text _n "Step 5: Computing employment stability indicators"

sort ${v_person_id} ${v_year}

* Lagged household earner count
by ${v_person_id}: gen hh_earners_lag = hh_earners_proxy[_n-1] if _n > 1

* Earner loss indicator (household lost earners relative to previous year)
gen earner_loss = (hh_earners_proxy < hh_earners_lag) if !missing(hh_earners_lag)
replace earner_loss = 0 if missing(earner_loss)

* Other earners (earners in household excluding reference person's own earner status)
gen other_earners = max(hh_earners_proxy - is_earner_proxy, 0)

* ============================================================================
* SECTION 6: ANNUAL POVERTY THRESHOLDS AND POVERTY STATUS
* ============================================================================

display as text _n "Step 6: Computing annual within-sample poverty thresholds"

* --- 6a. Weighted median equivalised income by year --------------------------

foreach yr of numlist 2016 2017 2018 2019 {
  quietly summarize eq_income [aweight = ${v_longitudinal_wt}] if ${v_year} == `yr', detail
  local median_`yr' = r(p50)
  display as text "  `yr': weighted median eq. income = " %12.2f `median_`yr''
}

* --- 6b. Poverty lines at each threshold -------------------------------------

foreach thr of numlist $thresholds {
  gen poverty_line_`thr' = .
  foreach yr of numlist 2016 2017 2018 2019 {
    local line = `median_`yr'' * (`thr' / 100)
    replace poverty_line_`thr' = `line' if ${v_year} == `yr'
  }
}

* Store median equivalised income
gen median_eq_income = .
foreach yr of numlist 2016 2017 2018 2019 {
  replace median_eq_income = `median_`yr'' if ${v_year} == `yr'
}

* --- 6c. Poverty status indicators -------------------------------------------

foreach thr of numlist $thresholds {
  * Binary poverty indicator
  gen poor_`thr' = (eq_income < poverty_line_`thr')
  
  * Poverty gap ratio
  gen poverty_gap_`thr' = max((poverty_line_`thr' - eq_income) / poverty_line_`thr', 0)
}

display as text _n "  Poverty rates (60% threshold):"
tab ${v_year} poor_60 [aweight = ${v_longitudinal_wt}], row nofreq

* --- 6d. Save poverty lines --------------------------------------------------

preserve
  keep ${v_year} median_eq_income poverty_line_50 poverty_line_60 poverty_line_70
  bysort ${v_year}: keep if _n == 1
  sort ${v_year}
  export delimited using "${outdir}/poverty_lines.csv", replace
  save "${outdir}/poverty_lines.dta", replace
restore

* ============================================================================
* SECTION 7: POVERTY SPELL CLASSIFICATION AND EUROSTAT PERSISTENT POVERTY
* ============================================================================

display as text _n "Step 7: Classifying poverty spells (${main_threshold}% threshold)"

* --- 7a. Reshape to wide for spell counting ----------------------------------

preserve

  keep ${v_person_id} ${v_year} ${v_longitudinal_wt} poor_${main_threshold}
  
  * Create year suffix for reshape
  gen year_suffix = ${v_year}
  
  reshape wide poor_${main_threshold}, i(${v_person_id} ${v_longitudinal_wt}) j(year_suffix)

  * --- 7b. Count poor years ----------------------------------------------------
  
  egen n_poor_4yr = rowtotal(poor_${main_threshold}2016 poor_${main_threshold}2017 ///
                             poor_${main_threshold}2018 poor_${main_threshold}2019)
  
  * Count poor years in 2016-2018 (previous 3 years)
  egen n_poor_previous3 = rowtotal(poor_${main_threshold}2016 poor_${main_threshold}2017 ///
                                   poor_${main_threshold}2018)
  
  * Current-year poverty (reference year)
  gen current_poor = (poor_${main_threshold}${reference_year} == 1)
  
  * --- 7c. Eurostat persistent poverty -----------------------------------------
  *     Poor in reference year AND poor in at least 2 of the previous 3 years
  
  gen persistent_poor = (current_poor == 1 & n_poor_previous3 >= 2)
  
  * --- 7d. Mutually exclusive four-year typology --------------------------------
  
  gen poverty_group = ""
  replace poverty_group = "Never poor"      if n_poor_4yr == 0
  replace poverty_group = "Transient poor"  if inlist(n_poor_4yr, 1, 2) & poverty_group == ""
  replace poverty_group = "Persistent poor" if persistent_poor == 1
  replace poverty_group = "Frequently poor" if n_poor_4yr == 3 & current_poor == 0 & poverty_group == ""
  
  * Additional flags
  gen never_poor   = (n_poor_4yr == 0)
  gen poor_once    = inlist(n_poor_4yr, 1, 2)
  gen poor_multiple = (n_poor_4yr == 3 & current_poor == 0)
  gen current_poor_not_persistent = (current_poor == 1 & persistent_poor == 0)
  
  display as text "  Poverty group distribution (unweighted):"
  tab poverty_group
  
  display as text "  Poverty group distribution (weighted):"
  tab poverty_group [aweight = ${v_longitudinal_wt}]
  
  * --- 7e. Export poverty typology -----------------------------------------------
  
  export delimited using "${outdir}/poverty_typology_${main_threshold}_${reference_year}_anchor.csv", replace
  
  * Save for merging back
  keep ${v_person_id} n_poor_4yr n_poor_previous3 current_poor persistent_poor ///
       poverty_group never_poor poor_once poor_multiple current_poor_not_persistent
  tempfile spell_data
  save `spell_data'

restore

* --- 7f. Merge spell classification back to panel ----------------------------

merge m:1 ${v_person_id} using `spell_data', assert(match) nogenerate

* ============================================================================
* SECTION 8: VARIABLE RECODING FOR DESCRIPTIVE TABLES AND MODELS
* ============================================================================

display as text _n "Step 8: Recoding variables for analysis"

* --- 8a. Sex -----------------------------------------------------------------

gen sex_recoded = ""
replace sex_recoded = "Female" if ${v_sex} == ${code_female}
replace sex_recoded = "Male"   if ${v_sex} == ${code_male}

encode sex_recoded, gen(sex_factor)
label define sex_lbl 1 "Female" 2 "Male", replace
label values sex_factor sex_lbl

* --- 8b. Age groups ----------------------------------------------------------

gen age_group = ""
replace age_group = "0-17"  if ${v_age} >= 0  & ${v_age} <= 17
replace age_group = "18-24" if ${v_age} >= 18 & ${v_age} <= 24
replace age_group = "25-34" if ${v_age} >= 25 & ${v_age} <= 34
replace age_group = "35-44" if ${v_age} >= 35 & ${v_age} <= 44
replace age_group = "45-54" if ${v_age} >= 45 & ${v_age} <= 54
replace age_group = "55-64" if ${v_age} >= 55 & ${v_age} <= 64
replace age_group = "65+"   if ${v_age} >= 65 & !missing(${v_age})

encode age_group, gen(age_group_factor)

* Age binary indicators (for 06_model CRE probit)
gen age_young   = (${v_age} < 30)
gen age_elderly_flag = (${v_age} > 64)

* --- 8c. Education -----------------------------------------------------------

gen education_recoded = ""
replace education_recoded = "Primary or below" if inlist(${v_education}, 0, 1, 2)
replace education_recoded = "Secondary"        if inlist(${v_education}, 3, 4, 5)
replace education_recoded = "Tertiary"         if inlist(${v_education}, 6, 7, 8)

encode education_recoded, gen(education_factor)
label define edu_lbl 1 "Primary or below" 2 "Secondary" 3 "Tertiary", replace
label values education_factor edu_lbl

* --- 8d. Employment / labour status ------------------------------------------

gen labour_recoded = ""
replace labour_recoded = "Employee"      if inlist(${v_labour_status}, 1, 2)
replace labour_recoded = "Self-employed" if inlist(${v_labour_status}, 3, 4)
replace labour_recoded = "Unemployed"    if ${v_labour_status} == 5
replace labour_recoded = "Retired"       if ${v_labour_status} == 7
replace labour_recoded = "Inactive"      if inlist(${v_labour_status}, 6, 8, 9, 10)

encode labour_recoded, gen(labour_factor)

* --- 8e. Social security registration ----------------------------------------

gen social_security_recoded = ""
replace social_security_recoded = "Registered"     if informal_status == 0
replace social_security_recoded = "Not registered" if informal_status == 1

encode social_security_recoded, gen(ss_factor)

* --- 8f. Household type factor -----------------------------------------------

encode household_type_derived, gen(hh_type_factor)

* --- 8g. Poverty status factor (60%) -----------------------------------------

gen poverty_status_60 = ""
replace poverty_status_60 = "Non-poor" if poor_60 == 0
replace poverty_status_60 = "Poor"     if poor_60 == 1

encode poverty_status_60, gen(poverty_60_factor)

* --- 8h. Profile spell group (for model tables) ------------------------------
*     Mutually exclusive typology used for p-values in profile tables

gen profile_spell_group = ""
replace profile_spell_group = "Never poor"      if n_poor_4yr == 0
replace profile_spell_group = "Transient poor"  if inlist(n_poor_4yr, 1, 2) & persistent_poor == 0
replace profile_spell_group = "Persistent poor" if persistent_poor == 1
replace profile_spell_group = "Frequently poor" if n_poor_4yr == 3 & persistent_poor == 0

encode profile_spell_group, gen(spell_group_factor)

* ============================================================================
* SECTION 9: DESCRIPTIVE STATISTICS
* ============================================================================

display as text _n "Step 9: Descriptive statistics"

* --- 9a. Continuous variables ------------------------------------------------

display as text _n "--- Table A: Continuous Variables (all person-years) ---"

foreach var in ${v_age} ${v_household_income} eq_income hh_size ///
               hh_children_u14 hh_children_u18 hh_elderly_65plus ///
               hh_earners_proxy hh_formal_earners hh_informal_earners ///
               dependency_ratio dependency_ratio_oecd {
  display as text _n "  Variable: `var'"
  summarize `var', detail
}

* --- 9b. Categorical variables -----------------------------------------------

display as text _n "--- Table B: Categorical Variables (all person-years) ---"

tab sex_recoded
tab age_group
tab education_recoded
tab labour_recoded
tab social_security_recoded
tab poverty_status_60

* --- 9c. Poverty rates and FGT indices ---------------------------------------

display as text _n "--- Poverty rates (weighted) ---"

foreach thr of numlist $thresholds {
  display as text _n "  Threshold: `thr'%"
  foreach yr of numlist 2016 2017 2018 2019 {
    quietly summarize poor_`thr' [aweight = ${v_longitudinal_wt}] if ${v_year} == `yr'
    local rate = r(mean) * 100
    
    * FGT1: poverty gap
    quietly summarize poverty_gap_`thr' [aweight = ${v_longitudinal_wt}] if ${v_year} == `yr'
    local fgt1 = r(mean)
    
    * FGT2: poverty severity (gap squared)
    quietly gen _gap_sq = poverty_gap_`thr'^2 if ${v_year} == `yr'
    quietly summarize _gap_sq [aweight = ${v_longitudinal_wt}] if ${v_year} == `yr'
    local fgt2 = r(mean)
    quietly drop _gap_sq
    
    display as text "    `yr': rate = " %5.1f `rate' "%  FGT1 = " %6.4f `fgt1' "  FGT2 = " %6.4f `fgt2'
  }
}

* --- 9d. Poverty group distribution ------------------------------------------

display as text _n "--- Poverty group distribution ---"
tab poverty_group
tab poverty_group [aweight = ${v_longitudinal_wt}]

* --- 9e. Poverty duration (number of years poor) ----------------------------

display as text _n "--- Poverty duration ---"
tab n_poor_4yr
tab n_poor_4yr [aweight = ${v_longitudinal_wt}]

* --- 9f. Profile table: means by poverty group (reference year) --------------

display as text _n "--- Socio-demographic profile by poverty group (${reference_year}) ---"

preserve
  keep if ${v_year} == ${reference_year}
  
  foreach grp in "Never poor" "Transient poor" "Persistent poor" {
    display as text _n "  Group: `grp'"
    display as text "  Age:"
    summarize ${v_age} [aweight = ${v_longitudinal_wt}] if poverty_group == "`grp'"
    display as text "  Female share:"
    summarize female [aweight = ${v_longitudinal_wt}] if poverty_group == "`grp'"
    display as text "  Informal share:"
    summarize informal_status [aweight = ${v_longitudinal_wt}] if poverty_group == "`grp'"
    display as text "  Eq. income:"
    summarize eq_income [aweight = ${v_longitudinal_wt}] if poverty_group == "`grp'"
    display as text "  Household size:"
    summarize hh_size [aweight = ${v_longitudinal_wt}] if poverty_group == "`grp'"
    display as text "  Dependency ratio (OECD):"
    summarize dependency_ratio_oecd [aweight = ${v_longitudinal_wt}] if poverty_group == "`grp'"
  }
restore

* --- 9g. Transition matrices -------------------------------------------------

display as text _n "--- Poverty transitions (60% threshold) ---"

sort ${v_person_id} ${v_year}

forvalues t = 2016/2018 {
  local t1 = `t' + 1
  display as text _n "  Transition `t' -> `t1':"
  
  preserve
    * Keep only the two adjacent years
    keep if inlist(${v_year}, `t', `t1')
    
    * Reshape to wide
    gen yr = ${v_year}
    reshape wide poor_${main_threshold}, i(${v_person_id}) j(yr)
    
    * Transition table (weighted)
    tab poor_${main_threshold}`t' poor_${main_threshold}`t1' [aweight = ${v_longitudinal_wt}], ///
        row nofreq
  restore
}

* ============================================================================
* SECTION 10: EXPORT ANALYSIS PANEL
* ============================================================================

display as text _n "Step 10: Exporting analysis panel"

* Label key variables
label variable eq_income            "Equivalised disposable income"
label variable hh_eq_size           "Household equivalence size (mod. OECD)"
label variable hh_size              "Household size"
label variable hh_children_u14      "No. children under 14 in household"
label variable hh_children_u18      "No. children under 18 in household"
label variable hh_elderly_65plus    "No. elderly 65+ in household"
label variable hh_earners_proxy     "No. earners in household (proxy)"
label variable hh_formal_earners    "No. formal earners in household"
label variable hh_informal_earners  "No. informal earners in household"
label variable dependency_ratio      "Standard dependency ratio"
label variable dependency_ratio_oecd "Modified OECD dependency ratio"
label variable informal_status       "Informal (not registered) = 1"
label variable female                "Female = 1"
label variable earner_loss           "Household lost earner(s) vs prior year = 1"
label variable other_earners         "No. other earners excl. own earner status"
label variable poor_50               "Poor at 50% threshold"
label variable poor_60               "Poor at 60% threshold"
label variable poor_70               "Poor at 70% threshold"
label variable poverty_gap_50        "Poverty gap at 50%"
label variable poverty_gap_60        "Poverty gap at 60%"
label variable poverty_gap_70        "Poverty gap at 70%"
label variable n_poor_4yr            "No. years poor (out of 4)"
label variable n_poor_previous3      "No. years poor in 2016-2018"
label variable current_poor          "Poor in reference year (2019)"
label variable persistent_poor       "Eurostat persistent poor"
label variable poverty_group         "Poverty typology (4 groups)"
label variable age_young             "Age < 30"
label variable age_elderly_flag      "Age > 64"

* Save final panel
save "${outdir}/analysis_panel_2016_2019.dta", replace
export delimited using "${outdir}/analysis_panel_2016_2019.csv", replace

display as text _n "=== Pipeline complete ==="
display as text "  Panel saved to: ${outdir}/analysis_panel_2016_2019.dta"
display as text "  Poverty lines:  ${outdir}/poverty_lines.dta"
