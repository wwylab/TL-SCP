# Supplementary Figure A.3: Fusion penalties under increasing censoring

This simulation examines estimation error and cluster recovery as the censoring proportion increases from 0.1 to 0.7. It uses 300 observations, 30 predictors, and 50 replicates per setting.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A3/Supplement-Figure-A3-from-rds.R
```

This reads `Supplement-Figure-A3.rds` and recreates `Supplement-Figure-A3.pdf` and its TikZ source (`Supplement-Figure-A3.tex`). No model fitting is required.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A3/Supplement-Figure-A3.R
```

The simulation generates its data and updates `Supplement-Figure-A3.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A3-paired-comparisons.csv`.
