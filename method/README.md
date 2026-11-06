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

The figure and table scripts call `fit_clinical()` in `real-data-analysis/analysis.R` to prepare cohorts, estimate preliminary mutation effects, combine ranks using carrier counts for each mutation, and fit the target model. The model includes one indicator per retained mutation and a sex covariate, which is excluded from fusion. The clinical solver uses the tuning settings from the manuscript analyses.

The MD Anderson analysis fits four TL-SCP groups, then assigns the R337H singleton to Mid onset using its weighted rank. This reporting step is implemented in `real-data-analysis/analysis.R`, separately from the general estimator in `methodology/`. On a fresh run, `raw_TLSCP_group` records the label before the R337H reassignment; the supplied summaries contain the reported labels. The SCP comparison retains its four fitted groups.

The scripts in `visualization/` regenerate figures and tables from saved summaries without fitting models or reading individual clinical records.

## Data and variables

[`real-data-analysis/analysis.R`](real-data-analysis/analysis.R) prepares the clinical inputs described below. Figure and table READMEs describe their saved results.

### Public TP53 and YTA inputs

The files are in [`tables/Main-Table-1/raw-data/`](../tables/Main-Table-1/raw-data/). The analyses use these columns from `GermlineDownload_r21.csv`:

| Column | Meaning and use |
| --- | --- |
| `Individual_code` | Individual identifier used to reduce multiple records to one observation per carrier. |
| `Family_ID` | Family identifier used to keep families together in DCSM data splits and to resample families in the Figure 7 bootstrap. |
| `Class` | Clinical family classification. The target uses `LFS` and `TP53_Chompret`; the source uses `LFL`, `FH`, `noFH`, and `Other`. |
| `Germline_carrier` | Carrier status; the analysis retains `Confirmed` and `Obligatory` carriers. |
| `Sex` | `M` for male and `F` for female. The model codes the `gender` column as 1 for `M` and 0 otherwise. |
| `ProtDescription` | Protein-level TP53 mutation identifier, including the `p.` prefix; used to match mutations across cohorts and to the YTA catalog. |
| `Age_at_diagnosis` | Age in years at the recorded cancer diagnosis. |
| `Age` | Age in years at follow-up, used when diagnosis age is missing. |

The code sets `time` to `Age_at_diagnosis` when present and to `Age` otherwise. It sets `status` to 1 for a recorded diagnosis and 0 for censoring. After removing missing sex or mutation values, `p.?`, and missing or nonpositive times, it keeps the earliest recorded diagnosis for each `Individual_code`, or one follow-up record if no diagnosis is recorded. The outcome is therefore age at first recorded cancer or censoring.

The analysis retains mutations with at least two target carriers and at least one source carrier, excluding those with only two or three target carriers and one source carrier. The [TP53 Database user manual](https://tp53.cancer.gov/help) describes the remaining database fields.

`yta-catalog.csv` is the analysis extract of `Table_Classes_TP53variants.xlsx`. It contains `ProtDescription` and `cluster_rec`, the categorical YTA class (`A`, `B`, `C`, or `D`). The derived `YTAClass` matches these labels by mutation. The code assigns class `0` to protein descriptions containing `*`, `fs`, `del`, or `insIQ`, and labels remaining unmatched mutations `Unclassified`.

### MD Anderson inputs

Set `MDA_DATA` to an authorized analysis CSV with one row per carrier. It requires `Individual_code`, `Sex`, and `ProtDescription` as defined above, plus `time` (positive age in years at first cancer or censoring) and `status` (1 = cancer, 0 = censored). Figure 7 also requires a nonmissing `Family_ID` for every carrier. The CSV must already be prepared for analysis; the reader checks the required column names and duplicate individual identifiers.

The alternative `MDAprospectiveMUT_clean.RData` input must contain exactly one data frame with `pid`, `sex`, `aminoacidchange`, and `diagAge`. These map to `Individual_code`, `Sex`, `ProtDescription`, and diagnosis age in years. Optional columns `age` and `ageDeath` give follow-up age and age at death in years. When `diagAge` is missing, the loader uses `age`, then `ageDeath` if `age` is missing. `Family_ID` is carried through when present. The loader derives `status` from the presence of `diagAge`, removes incomplete or nonpositive observations, and keeps the earliest diagnosis per individual. Individual-level MD Anderson records are not distributed.

### Shared clinical results

The alluvial figures and mutation tables use annotation files with one row per `ProtDescription`. `YTAClass` is the functional class; `SCP` and `TL-SCP` are the reported onset-group labels. `target_n` and `source_n` give the number of carriers of that mutation in each analysis cohort.

Survival summaries contain a `curve` data frame with age in years (`time`), group label (`Group`), estimated cancer-free survival (`surv`), and pointwise confidence limits (`lower`, `upper`). The `counts` data frame gives `MutationN` and `PatientN` for each group; `p` is the log-rank p-value and `confidence` is the confidence level. The figure and table READMEs identify the groups in each file.

## Simulation results

The simulation scripts generate their own data. Saved result rows identify the replicate (`rep`), method (`method`), and experimental setting (`setting`, defined in each figure README). `n_S` is the source sample size and `K` the selected number of coefficient groups, where present. The main measures are:

| Variable | Definition |
| --- | --- |
| `relErr` | Squared relative coefficient error: `sum((beta_hat - beta)^2) / sum(beta^2)`. Lower values indicate smaller estimation error. |
| `NMI` | Mutual information between true and estimated groups divided by the larger group-label entropy, as computed in `simulation/helpers.R`. Higher values indicate better agreement. |
| `Runtime` | Elapsed fitting time in seconds, where recorded. |

The figure READMEs describe their saved summaries and target-only control files. Simulated survival times follow the scale of the data-generating model; event status is 1 for an event and 0 for censoring.
