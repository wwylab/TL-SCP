# Main Figure 7: Comparing mutation rankings across cohorts

The first two panels compare weighted mutation rankings from the *TP53* Database and MD Anderson analyses. The third compares preliminary Cox effects, relative to R175H, for mutations assigned to different onset groups. Standard errors are based on 1,000 family-bootstrap samples per cohort, with the two variances combined for the intervals.

## Data and variables

`Main-Figure-7.rds` contains two data frames:

- `ranks` has one row per `Mutation` and `cohort` (`NCI` for the public *TP53* Database, `MDA` for MD Anderson). `weighted_rank` averages the ascending preliminary Cox ranks using carrier counts: `target_weight = target_n / (target_n + source_n)` and source weight `1 - target_weight`. `group` gives the reported onset label; `highlight` marks discordant annotations.
- `forest` contains discordant mutations. `beta_NCI` and `beta_MDA` are preliminary target Cox log-hazard effects relative to R175H; `SE_NCI` and `SE_MDA` are family-bootstrap standard errors. `delta` is `beta_MDA - beta_NCI`; `SE_delta = sqrt(SE_NCI^2 + SE_MDA^2)`, and `lower`/`upper` give the approximate 95% interval.

Clinical input variables are defined in [Data and variables](../../method/README.md#data-and-variables).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-7/Main-Figure-7-from-rds.R
```

This reads `Main-Figure-7.rds` and recreates `Main-Figure-7.pdf` and its TikZ source (`Main-Figure-7.tex`). No model fitting is required.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-7/Main-Figure-7.R
```

The public *TP53* Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

To refit the MD Anderson analysis, set `MDA_DATA` to an authorized CSV or the reviewed `MDAprospectiveMUT_clean.RData` file; see [MD Anderson inputs](../../method/README.md#md-anderson-inputs). The saved `.rds` results regenerate this figure without individual patient data.

To test this script with the synthetic target cohort and public source data, follow the [synthetic-data instructions](../../method/README.md#synthetic-md-anderson-data), which use a temporary copy to preserve the published results.

Reported TL-SCP labels include the rank-based reassignment of R337H to Mid onset. The [clinical analysis documentation](../../method/README.md#real-data-analyses) gives the fitting and reporting rules.

The family bootstrap also requires `Family_ID` in the MD Anderson analysis data and Python with `numpy`, `pandas`, and `scipy`. Set `PYTHON` if these are installed in a separate environment.
