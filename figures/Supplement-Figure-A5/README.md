# Supplementary Figure A.5: Transfer with consistent effect rankings

Source and target coefficient rankings agree, although their effect magnitudes differ. With 200 target observations, this simulation examines how source sample size (500, 1,000, 1,500, or 2,000) affects estimation error and cluster recovery. Each setting has 50 replicates.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A5/Supplement-Figure-A5-from-rds.R
```

This reads `Supplement-Figure-A5.rds` and recreates `Supplement-Figure-A5.pdf` and its TikZ source (`Supplement-Figure-A5.tex`). No model fitting is required.

## Rerun the simulation

```sh
Rscript figures/Supplement-Figure-A5/Supplement-Figure-A5.R
```

The simulation generates its data and updates `Supplement-Figure-A5.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The distance-to-set fits use the RcppEigen solver provided in `method/`; compiling it requires a working C++ toolchain.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A5-paired-comparisons.csv`.
