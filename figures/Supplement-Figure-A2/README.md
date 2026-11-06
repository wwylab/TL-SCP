# Supplementary Figure A.2: Fusion penalties under increasing noise

This simulation compares Lasso, SCAD, MCP, and distance-to-set penalties across noise levels, using 300 observations and 30 predictors in three coefficient groups. Each setting has 50 replicates.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A2/Supplement-Figure-A2-from-rds.R
```

This reads `Supplement-Figure-A2.rds` and recreates `Supplement-Figure-A2.pdf` and its TikZ source (`Supplement-Figure-A2.tex`). No model fitting is required.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A2/Supplement-Figure-A2.R
```

The simulation generates its data and updates `Supplement-Figure-A2.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A2-paired-comparisons.csv`.
