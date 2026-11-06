# SCP and TL-SCP

SCP jointly estimates predictor effects and groups predictors by fusing adjacent Cox coefficients in an estimated ordering. TL-SCP combines source and target ranks before applying the same fusion procedure to the target outcomes.

The code is organized by purpose:

```text
method/
├── methodology/
│   ├── scp_tlscp.R          SCP and TL-SCP estimation
│   ├── scp_penalties.R      Lasso, SCAD, MCP, and all-pair comparisons
│   └── fusion_solver.cpp   Compiled fusion solver for the simulations
├── simulation/
│   ├── helpers.R           Data generation and performance measures
│   └── rank_discordance.R  Source–target ranking-disagreement designs
├── real-data-analysis/
│   ├── analysis.R          Cohort preparation, fitting, and survival summaries
│   ├── fit_groups.R        Sex-adjusted grouped Cox fits and BIC selection
│   ├── fusion_solver.R     Fusion solver used in the clinical analyses
│   ├── ridge_cox.py        Preliminary Cox fits for the family bootstrap
│   ├── family_bootstrap.py Family-bootstrap standard errors
│   └── dcsm.py             DCSM comparison
└── visualization/
    ├── figures.R           Kaplan–Meier, alluvial, and simulation figures
    ├── alluvials.py        TikZ alluvial diagrams with mutation names
    ├── ranks.R             Weighted-rank and coefficient-comparison panels
    └── tables.R            LaTeX tables from saved mutation summaries
```

## Fitting SCP and TL-SCP

The main functions are `SCP()` and `trans_SCP()` in `methodology/scp_tlscp.R`. From the repository root:

```r
source("method/simulation/helpers.R")
source("method/methodology/scp_tlscp.R")

fit <- trans_SCP(Y_target, X_target, Y_source, X_source,
                 penalty = "distance", nclusters = 1:5)
selected <- which.min(fit$IC_path)
beta <- fit$Beta_path[, selected]
```

The `Y` matrices contain `time` and `status`, with status 1 for an event and 0 for censoring. Source and target design columns must identify the same predictors in the same order. `SCP()` takes `y`, `status`, and `X` separately. Both functions return the coefficient path, iteration counts, and BIC values in `IC_path`. Supplementary Figure A.5 provides a complete simulated example.

`methodology/fusion_solver.cpp` accelerates the simulation fits through RcppEigen. The simulation scripts call `use_fast_distance(root)` to load it. Figure A.5 uses `r_sort = TRUE` to retain the original R solver’s ordering of tied coefficient differences; Figures A.6–A.7 retain the stable ordering used in their simulations. The penalty comparisons in Supplementary Figures A.1–A.4 use `methodology/scp_penalties.R` and require the Gurobi R interface.

## Real-data analyses

The figure and table scripts call `fit_clinical()` in `real-data-analysis/analysis.R`. This prepares the cohorts, estimates preliminary mutation effects, combines ranks using mutation-specific carrier counts, and fits the target model with sex excluded from fusion. The clinical solver retains the tuning settings used for the manuscript analyses.

The reported MD Anderson labels include the manuscript's rank-based reassignment of R337H to Mid onset after fitting four groups. This step is explicit in `real-data-analysis/analysis.R`. It is not part of the general estimator in `methodology/`.

The scripts in `visualization/` read the saved summaries. They do not fit models or require individual clinical records.
