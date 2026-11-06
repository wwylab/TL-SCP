# Supplementary Figure A.1: Ranked-adjacent versus all-pair fusion

This simulation compares ranked-adjacent fusion with fusion across all predictor pairs under the Lasso penalty. It uses 300 observations, 30 predictors, four noise levels, and 50 replicates per setting.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A1/Supplement-Figure-A1-from-rds.R
```

This reads `Supplement-Figure-A1.rds` and recreates `Supplement-Figure-A1.pdf` and its TikZ source (`Supplement-Figure-A1.tex`). No model fitting is required.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A1/Supplement-Figure-A1.R
```

The simulation generates its data and updates `Supplement-Figure-A1.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A1-paired-comparisons.csv`.
