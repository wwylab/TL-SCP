# Supplementary Figure A.1: Ranked-adjacent versus all-pair fusion

This simulation compares ranked-adjacent fusion with fusion across all predictor pairs under the Lasso penalty. It uses 300 observations, 30 predictors, four noise levels, and 50 replicates per setting.

## Data and variables

`Supplement-Figure-A1.rds` is a data frame with one row for each combination of `setting`, `rep`, and `method`. `setting` is the noise standard deviation (0.1, 0.3, 0.7, or 1), and `method` distinguishes ranked-adjacent `SCP-Lasso` from `SCP-Allpair`. The three panels plot `relErr`, `NMI`, and `Runtime` (seconds). See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A1/Supplement-Figure-A1-from-rds.R
```

The script reads `Supplement-Figure-A1.rds` and regenerates `Supplement-Figure-A1.pdf` and its TikZ source (`Supplement-Figure-A1.tex`) without refitting the models.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A1/Supplement-Figure-A1.R
```

The script generates simulated data and overwrites `Supplement-Figure-A1.rds`. Copy the saved file before rerunning to retain the reported results.

The penalty comparisons require the Gurobi R interface and a valid Gurobi license. Gurobi is not needed to redraw the figure from saved results.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A1-paired-comparisons.csv`.
