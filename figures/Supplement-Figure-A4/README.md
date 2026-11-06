# Supplementary Figure A.4: Fusion penalties as dimension increases

This simulation compares four SCP penalties with 300 observations and 30, 60, 90, or 120 predictors. There are 50 replicates at the first three dimensions and 48 completed replicates at 120 predictors, matching the reported experiment.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A4/Supplement-Figure-A4-from-rds.R
```

This reads `Supplement-Figure-A4.rds` and recreates `Supplement-Figure-A4.pdf` and its TikZ source (`Supplement-Figure-A4.tex`). No model fitting is required.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A4/Supplement-Figure-A4.R
```

The simulation generates its data and updates `Supplement-Figure-A4.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A4-paired-comparisons.csv`.

The final two replicates at 120 predictors were not completed in the reported experiment; the script reproduces the same 48-replicate setting.
