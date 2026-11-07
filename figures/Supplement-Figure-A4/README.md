# Supplementary Figure A.4: Fusion penalties as dimension increases

This simulation compares four SCP penalties with 300 observations and 30, 60, 90, or 120 predictors. There are 50 replicates at the first three dimensions and 48 completed replicates at 120 predictors.

## Data and variables

`Supplement-Figure-A4.rds` is a data frame with one row for each combination of `setting`, `rep`, and `method`. `setting` is the number of predictors (30, 60, 90, or 120), and `method` identifies the Lasso, SCAD, MCP, or distance-to-set SCP penalty. The panels plot `relErr` and `NMI`; the saved `Runtime` column is not plotted. See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A4/Supplement-Figure-A4-from-rds.R
```

The script reads `Supplement-Figure-A4.rds` and regenerates `Supplement-Figure-A4.pdf` and its TikZ source (`Supplement-Figure-A4.tex`) without refitting the models.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A4/Supplement-Figure-A4.R
```

The script generates simulated data and overwrites `Supplement-Figure-A4.rds`. Copy the saved file before rerunning to retain the reported results.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A4-paired-comparisons.csv`.

The final two replicates at 120 predictors were not completed in the reported experiment; the script reproduces the same 48-replicate setting.
