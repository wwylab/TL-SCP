# Supplementary Table F.1: Mutation annotations in the primary analysis

This table reports carrier counts and TL-SCP, SCP, and YTA annotations for 31 mutations in the primary TP53 Database analysis. Mutations are ordered by TL-SCP onset group.

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
