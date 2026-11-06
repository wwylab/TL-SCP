# Main Table 1: Distribution of YTA classes

This table summarizes mutation and carrier counts by YTA class for the two TP53 Database cohorts and the MD Anderson cohort.

## Data and variables

`Main-Table-1.rds` contains `database` and `mda` data frames with one row per analysis mutation. For each `YTAClass`, the script counts mutations and sums `target_n` and `source_n` for the database cohorts, and `target_n` from `mda` for MD Anderson. `SCP` and `TL-SCP` are retained in the data but are not used in this table. See [Data and variables](../../method/README.md#data-and-variables) for field definitions and the public input schema.

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

The public TP53 Database release 21 and YTA catalog are in `tables/Main-Table-1/raw-data/` (`raw-data/` in this folder). The file `yta-catalog.csv` contains the two catalog columns used by the R scripts. Cohort definitions and mutation filters follow Section 5.1 of the manuscript. External datasets retain their original attribution and terms.

To refit the MD Anderson analysis, set `MDA_DATA` to an authorized CSV or the reviewed `MDAprospectiveMUT_clean.RData` file; see [MD Anderson inputs](../../method/README.md#md-anderson-inputs). The saved `.rds` results regenerate this table without individual patient data.
