# Supplementary Figure A.2: Fusion penalties under increasing noise

This simulation compares Lasso, SCAD, MCP, and distance-to-set penalties across noise levels, using 300 observations and 30 predictors in three coefficient groups. Each setting has 50 replicates.

## Data and variables

`Supplement-Figure-A2.rds` is a data frame with one row for each combination of `setting`, `rep`, and `method`. `setting` is the noise standard deviation (0.1, 0.3, 0.7, or 1), and `method` identifies the Lasso, SCAD, MCP, or distance-to-set SCP penalty. The panels plot `relErr` and `NMI`; the saved `Runtime` column is not plotted. See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A2/Supplement-Figure-A2-from-rds.R
```

The script reads `Supplement-Figure-A2.rds` and regenerates `Supplement-Figure-A2.pdf` and its TikZ source (`Supplement-Figure-A2.tex`) without refitting the models.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A2/Supplement-Figure-A2.R
```

The script generates simulated data and overwrites `Supplement-Figure-A2.rds`. Copy the saved file before rerunning to retain the reported results.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A2-paired-comparisons.csv`.
