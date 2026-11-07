# Supplementary Table F.1: Mutation annotations in the primary analysis

This table reports carrier counts and TL-SCP, SCP, and YTA annotations for 31 mutations in the primary *TP53* Database analysis. Mutations are ordered by TL-SCP onset group.

## Data and variables

`Supplement-Table-F1.rds` is a mutation-level data frame with the six fields used in the table: `ProtDescription`, `target_n`, `source_n`, `TL-SCP`, `SCP`, and `YTAClass`. Rows are sorted by TL-SCP onset group and mutation; carrier counts are summed for the total row. The table displays `Unclassified` YTA annotations as `Unknown`. Common annotation and cohort definitions are in [Data and variables](../../method/README.md#data-and-variables).

## Reproduce from saved results

From the repository root:

```sh
Rscript tables/Supplement-Table-F1/Supplement-Table-F1-from-rds.R
```

This reads `Supplement-Table-F1.rds` and recreates `Supplement-Table-F1.tex` as a LaTeX table (`booktabs` required). No model fitting is required.

## Rerun the analysis

```sh
Rscript tables/Supplement-Table-F1/Supplement-Table-F1.R
```
