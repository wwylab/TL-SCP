# Supplementary Table F.2: Mutation annotations at MD Anderson

This table lists 28 MD Anderson mutations, source-cohort carrier counts, and the corresponding onset annotations. The reported TL-SCP groups include the rank-based reassignment of R337H.

## Data and variables

`Supplement-Table-F2.rds` contains one row per MD Anderson analysis mutation, with `ProtDescription`, `target_n`, `source_n`, `TL-SCP`, `SCP`, and `YTAClass`. The table script orders mutations by the reported TL-SCP onset group, sums carrier counts, and displays `Unclassified` YTA annotations as `Unknown`. See [Data and variables](../../method/README.md#data-and-variables) for shared definitions.

## Reproduce from saved results

From the repository root:

```sh
Rscript tables/Supplement-Table-F2/Supplement-Table-F2-from-rds.R
```

This reads `Supplement-Table-F2.rds` and recreates `Supplement-Table-F2.tex` as a LaTeX table (`booktabs` required). No model fitting is required.

## Rerun the analysis

```sh
Rscript tables/Supplement-Table-F2/Supplement-Table-F2.R
```

To refit the MD Anderson analysis, set `MDA_DATA` to an authorized CSV or the reviewed `MDAprospectiveMUT_clean.RData` file; see [MD Anderson inputs](../../method/README.md#md-anderson-inputs). The saved `.rds` results regenerate this table without individual patient data.

Reported TL-SCP labels include the rank-based reassignment of R337H to Mid onset. The [clinical analysis documentation](../../method/README.md#real-data-analyses) gives the fitting and reporting rules.
