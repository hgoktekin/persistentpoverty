# persistentpoverty
Persistent Poverty analysis using TR-SILC data between 2016-2019

# SILC Türkiye Poverty Dynamics Workflow
This workflow implements a reproducible empirical analysis of poverty dynamics in the balanced 2016-2019 Turkish SILC panel.

Run:
```
source("R/02_run_silc_poverty_analysis.R")
```

The pipeline proceeds in the order required for an article:

1. Data preparation and balanced panel construction.
2. Modified OECD equivalence scale and equivalised disposable income.
3. Annual within-sample poverty thresholds at 50%, 60%, and 70% of the weighted median.
4. Poverty classifications, including Eurostat persistent poverty anchored on 2019.
5. Weighted FGT poverty indices.
6. Descriptive tables, duration/spell tables, and transition matrices.
7. Weighted probit model with robust standard errors and marginal effects.
8. Publication-ready poverty trend figure with 95% confidence intervals.

Main scripts:
- `R/00_config.R`: paths, years, thresholds, and variable-name mapping.
- `R/01_functions.R`: reusable functions with methodological comments.
- `R/02_run_silc_poverty_analysis.R`: end-to-end analysis runner.

Key methodological choices are documented directly in the code. The balanced panel reduces attrition-related missingness in transition estimates, but it may introduce survivorship bias because the sample excludes individuals not observed in all four waves. Recommended robustness checks include comparing baseline characteristics of stayers and attritors, using cross-sectional weights for year-specific descriptive checks where available, and repeating the analysis at 50% and 70% thresholds.
