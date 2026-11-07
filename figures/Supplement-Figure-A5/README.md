# Supplementary Figure A.5: Transfer with consistent effect rankings

Source and target coefficient rankings agree, although their effect magnitudes differ. With 200 target observations, this simulation examines how source sample size (500, 1,000, 1,500, or 2,000) affects estimation error and cluster recovery. Each setting has 50 replicates.

## Data and variables

`Supplement-Figure-A5.rds` is a data frame with one row for each combination of `setting`, `rep`, and `method`. `setting` is the source sample size (500, 1,000, 1,500, or 2,000), also stored as `n_S` in the supplied results. The panels plot `relErr` and `NMI` for `Cox-K-means`, `SCP`, and `TL-SCP`. See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A5/Supplement-Figure-A5-from-rds.R
```

The script reads `Supplement-Figure-A5.rds` and regenerates `Supplement-Figure-A5.pdf` and its TikZ source (`Supplement-Figure-A5.tex`) without refitting the models.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A5/Supplement-Figure-A5.R
```

The script generates simulated data and overwrites `Supplement-Figure-A5.rds`. Copy the saved file before rerunning to retain the reported results.

The distance-to-set fits use the RcppEigen solver provided in `method/`; compiling it requires a working C++ toolchain.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A5-paired-comparisons.csv`.
