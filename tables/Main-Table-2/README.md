# Main Table 2: Agreement between cohort-specific TL-SCP annotations

This table compares TL-SCP annotations for the 17 mutations shared by the primary TP53 Database and MD Anderson analyses. Nine mutations have the same onset assignment; the other eight differ by one adjacent group. Carrier counts for discordant mutations appear in parentheses.

## Data and variables

`Main-Table-2.rds` is a data frame joining the two mutation-annotation summaries by `ProtDescription`. Column suffixes `_Database` and `_MDA` identify the analysis. `TL-SCP_Database` and `TL-SCP_MDA` determine the table cells; `target_n_Database` and `target_n_MDA` supply the carrier counts shown for discordant mutations. Other annotation and source-count columns are retained for reference. Common field definitions are in [Data and variables](../../method/README.md#data-and-variables).

## Reproduce from saved results

From the repository root:

```sh
Rscript tables/Main-Table-2/Main-Table-2-from-rds.R
```

This reads `Main-Table-2.rds` and recreates `Main-Table-2.tex` as a LaTeX table (`booktabs` required). No model fitting is required.

## Rerun the analysis

```sh
Rscript tables/Main-Table-2/Main-Table-2.R
```

To refit the MD Anderson analysis, set `MDA_DATA` to an authorized CSV or the reviewed `MDAprospectiveMUT_clean.RData` file; see [MD Anderson inputs](../../method/README.md#md-anderson-inputs). The saved `.rds` results regenerate this table without individual patient data.
