# Main Table 1: Distribution of YTA classes

This table summarizes mutation and carrier counts by YTA class for the two TP53 Database cohorts and the MD Anderson cohort. It describes the functional-class composition of the cohorts used in the survival analyses.

## Reproduce from saved results

From the repository root:

```sh
Rscript tables/Main-Table-1/Main-Table-1-from-rds.R
```

This reads `Main-Table-1.rds` and recreates `Main-Table-1.tex` as a LaTeX table (`booktabs` required). No model fitting is required.

## Rerun the analysis

```sh
Rscript tables/Main-Table-1/Main-Table-1.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

**MD Anderson data:** Individual-level clinical and family records are not distributed. To rerun the clinical analysis, set `MDA_DATA` to an authorized analysis CSV with `Individual_code`, `Sex` (`M`/`F`), `ProtDescription`, `time`, and `status` (1 = cancer, 0 = censored). The preparation code also accepts the reviewed `MDAprospectiveMUT_clean.RData` file. The saved `.rds` results can be used without access to these records.

The `raw-data/` folder contains the TP53 Database release 21 and YTA catalog; `yta-catalog.csv` holds the two catalog columns used by the R scripts. External datasets retain their original attribution and terms.
