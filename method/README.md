# SCP and TL-SCP

SCP groups predictors with similar Cox regression coefficients by fusing adjacent coefficients after ranking their estimated effects. TL-SCP incorporates information from a source cohort through a weighted combination of source and target rankings, while estimating the grouped effects using the target data.

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

For `trans_SCP()`, `Y_target` and `Y_source` are two-column matrices containing survival time and event status, coded as 1 for an event and 0 for censoring. The columns of `X_target` and `X_source` must correspond to the same predictors in the same order. The `SCP()` function takes survival time, event status, and the design matrix as separate arguments.

Both functions return the estimated coefficient path and BIC values (`IC_path`), along with iteration information. Supplementary Figure A.5 gives a complete example using simulated data.

The distance-to-set computations in Supplementary Figures A.5–A.7 use the C++ implementation in `methodology/fusion_solver.cpp`, loaded through RcppEigen. The penalty comparisons in Supplementary Figures A.1–A.4 use `methodology/scp_penalties.R` and require Gurobi.

## Real-data analyses

The real-data analyses use sex-adjusted Cox models to estimate mutation effects and identify groups with similar effects. For TL-SCP, the source and target mutation rankings are combined using mutation-specific carrier counts, while the grouped effects are estimated from the target cohort. The analyses are implemented in `real-data-analysis/analysis.R`.

For the MD Anderson analysis, TL-SCP initially produces four groups. The R337H singleton is then assigned to the Mid-onset group using its weighted rank. This adjustment is specific to the reported MD Anderson analysis and is not part of the general SCP or TL-SCP estimator. The SCP comparison retains the four fitted groups.

The plotting functions in `visualization/` use saved numerical results to regenerate the figures and tables without refitting the models.

## Data and variables

[`real-data-analysis/analysis.R`](real-data-analysis/analysis.R) prepares the clinical inputs described below. Figure and table READMEs describe their saved results.

### Public *TP53* and YTA inputs

The files are in [`tables/Main-Table-1/raw-data/`](../tables/Main-Table-1/raw-data/). The analyses use these columns from `GermlineDownload_r21.csv`:

| Column | Meaning and use |
| --- | --- |
| `Individual_code` | Individual identifier used to reduce multiple records to one observation per carrier. |
| `Family_ID` | Family identifier used to keep families together in DCSM data splits and to resample families in the Figure 7 bootstrap. |
| `Class` | Clinical family classification. The target uses `LFS` and `TP53_Chompret`; the source uses `LFL`, `FH`, `noFH`, and `Other`. |
| `Germline_carrier` | Carrier status; the analysis retains `Confirmed` and `Obligatory` carriers. |
| `Sex` | `M` for male and `F` for female. The model codes the `gender` column as 1 for `M` and 0 otherwise. |
| `ProtDescription` | Protein-level *TP53* mutation identifier, including the `p.` prefix; used to match mutations across cohorts and to the YTA catalog. |
| `Age_at_diagnosis` | Age in years at the recorded cancer diagnosis. |
| `Age` | Age in years at follow-up, used when diagnosis age is missing. |

The survival outcome is age at first recorded cancer diagnosis or censoring. When `Age_at_diagnosis` is recorded, we use it as `time` and set `status` to 1; otherwise, we use the follow-up age in `Age` and set `status` to 0. Survival times are in years. When an individual has multiple records, we retain the earliest cancer diagnosis or one follow-up record if no cancer is recorded.

We exclude observations with missing sex or mutation information, unknown mutation labels (`p.?`), or missing or nonpositive survival times.

The analysis retains mutations with at least two target carriers and at least one source carrier, excluding those with only two or three target carriers and one source carrier. The [*TP53* Database user manual](https://tp53.cancer.gov/help) describes the remaining database fields.

`yta-catalog.csv` is the analysis extract of `Table_Classes_TP53variants.xlsx`. It contains `ProtDescription` and `cluster_rec`, the categorical YTA class (`A`, `B`, `C`, or `D`). The derived `YTAClass` matches these labels by mutation. The code assigns class `0` to protein descriptions containing `*`, `fs`, `del`, or `insIQ`, and labels remaining unmatched mutations `Unclassified`.

### MD Anderson inputs

Refitting the original MD Anderson models requires individual-level data that are not included in this repository. Researchers with authorized access can specify the input file through `MDA_DATA`. The source cohort is prepared from the included public *TP53* Database records.

For a prepared CSV file, each row represents one carrier. The required variables are:

| Variable | Definition |
| --- | --- |
| `Individual_code` | Individual identifier, unique to each carrier. |
| `Sex` | `M` for male and `F` for female, with the model coding described above. |
| `ProtDescription` | Protein-level *TP53* mutation identifier, including the `p.` prefix. |
| `time` | Positive age in years at first cancer or censoring. |
| `status` | 1 for cancer and 0 for censoring. |

Figure 7 also requires a nonmissing `Family_ID` for every carrier to identify families for the bootstrap. The CSV must already be prepared for analysis; the reader checks the required column names and duplicate individual identifiers.

The alternative `MDAprospectiveMUT_clean.RData` input must contain exactly one data frame with `pid`, `sex`, `aminoacidchange`, and `diagAge`. These map to `Individual_code`, `Sex`, `ProtDescription`, and diagnosis age in years. Optional columns `age` and `ageDeath` give follow-up age and age at death in years. When `diagAge` is missing, the loader uses `age`, then `ageDeath` if `age` is missing. `Family_ID` is carried through when present. The loader derives `status` from the presence of `diagAge`, removes incomplete or nonpositive observations, and keeps the earliest diagnosis per individual.

### Synthetic MD Anderson data

The file `real-data-analysis/MDA-synthetic.csv` contains 202 simulated carriers and 28 mutation categories, matching the sample size and number of mutations in the MD Anderson target analysis. The TL-SCP analysis uses 329 carriers from the public *TP53* Database as the source cohort. All individual and family identifiers in the synthetic dataset are artificial.

Survival times are generated from a Cox proportional hazards model with three mutation-effect groups and independent right censoring. The generation script, `real-data-analysis/generate_MDA_synthetic.R`, uses a fixed random seed (`2024`). No individual-level MD Anderson records or fitted survival estimates are used in the simulation.

Generation uses base R, with a constant baseline hazard of 0.05 per year and no subject-specific noise. Mutation coefficients are 2, 0.5, and -1 for 10, 9, and 9 categories, respectively. Each carrier has one binary mutation indicator; sex is sampled with equal probabilities for M and F and has coefficient 0.5 (M = 1, F = 0). Mutations have 7 or 8 carriers, paired within families with a singleton where needed. Independent exponential censoring uses a common rate of approximately 0.0373 per year, calculated before drawing outcomes for 30% expected censoring. Observed time is the minimum of the event and censoring times, with `status = 1` when the event occurs first and 0 otherwise.

From the repository root, regenerate the synthetic target CSV with:

```sh
Rscript method/real-data-analysis/generate_MDA_synthetic.R
```

Set `MDA_DATA` to the synthetic target CSV to test the clinical analysis; the source data are read from the public *TP53* Database files. Run in a temporary copy because the analysis scripts overwrite their saved results and figures or tables:

```sh
workdir=$(mktemp -d)
cp -R method figures tables "$workdir/"
(
  cd "$workdir" || exit
  export MDA_DATA="$PWD/method/real-data-analysis/MDA-synthetic.csv"
  Rscript figures/Main-Figure-5/Main-Figure-5.R
  Rscript figures/Main-Figure-6/Main-Figure-6.R
  Rscript figures/Main-Figure-7/Main-Figure-7.R
  Rscript tables/Main-Table-1/Main-Table-1.R
  Rscript tables/Main-Table-2/Main-Table-2.R
  Rscript tables/Supplement-Table-F2/Supplement-Table-F2.R
)
```

The mutation names allow the scripts to use the existing YTA catalog and compare mutations shared with the public cohort. The MD Anderson reporting step, including the R337H reassignment, is applied unchanged to the synthetic fits. Main Figure 7 uses the existing family bootstrap with R175H as the reference mutation and requires Python with `numpy`, `pandas`, and `scipy`; set `PYTHON` if these are installed in a separate environment. Use the supplied `.rds` files to regenerate all published MD Anderson figures and tables.

<a id="shared-clinical-results"></a>

### Saved clinical results

The mutation annotation files contain one row per `ProtDescription`, with the YTA class (`YTAClass`), the onset groups assigned by SCP and TL-SCP, and the numbers of target and source carriers (`target_n` and `source_n`).

The saved survival summaries contain the estimated cancer-free survival curves and their pointwise confidence intervals. The `curve` data frame records age (`time`), group (`Group`), survival probability (`surv`), and confidence limits (`lower` and `upper`). The `counts` data frame contains the numbers of mutations (`MutationN`) and carriers (`PatientN`) in each group. The object also records the log-rank p-value (`p`) and confidence level (`confidence`).

On a fresh MD Anderson run, `raw_TLSCP_group` records the label before the R337H reassignment; the supplied summaries contain the reported labels.

## Simulation results

Each simulation generates its own data. The saved results identify the replicate (`rep`), method (`method`), and simulation setting (`setting`, or `D_target` in Supplementary Figure A.7). Depending on the experiment, they also record the source sample size (`n_S`) and number of selected groups (`K`). Performance is evaluated using the following measures:

| Variable | Definition |
| --- | --- |
| `relErr` | Squared relative coefficient error: `sum((beta_hat - beta)^2) / sum(beta^2)`. Lower values indicate smaller estimation error. |
| `NMI` | Mutual information between true and estimated groups divided by the larger group-label entropy, as computed in `simulation/helpers.R`. Higher values indicate better agreement. |
| `Runtime` | Elapsed fitting time in seconds, where recorded. |

The individual figure READMEs describe the simulation settings and saved results. Survival times are generated according to the model specified in each simulation, with event status coded as 1 for an event and 0 for censoring.

The simulation scripts load the C++ solver with `use_fast_distance(root)`. Figure A.5 uses `r_sort = TRUE` to retain the original R solver’s ordering of tied coefficient differences; Figures A.6–A.7 retain the stable ordering used in their simulations.
