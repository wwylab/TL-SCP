# Supplementary Figure A.3: Fusion penalties under increasing censoring

This simulation examines estimation error and cluster recovery as the censoring proportion increases from 0.1 to 0.7. It uses 300 observations, 30 predictors, and 50 replicates per setting.

## Data and variables

`Supplement-Figure-A3.rds` is a data frame with one row for each combination of `setting`, `rep`, and `method`. `setting` is the censoring proportion (0.1, 0.3, 0.5, or 0.7), and `method` identifies the Lasso, SCAD, MCP, or distance-to-set SCP penalty. The panels plot `relErr` and `NMI`; the saved `Runtime` column is not plotted. See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A3/Supplement-Figure-A3-from-rds.R
```

The script reads `Supplement-Figure-A3.rds` and regenerates `Supplement-Figure-A3.pdf` and its TikZ source (`Supplement-Figure-A3.tex`) without refitting the models.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A3/Supplement-Figure-A3.R
```

The script generates simulated data and overwrites `Supplement-Figure-A3.rds`. Copy the saved file before rerunning to retain the reported results.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A3-paired-comparisons.csv`.
