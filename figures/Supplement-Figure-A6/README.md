# Supplementary Figure A.6: Transfer under covariate shift

This simulation introduces Gaussian noise into the source covariates while preserving the ordering of mutation effects. It compares TL-SCP with target-only SCP and Cox K-means for source sample sizes of 500 and 1,000, using 50 replicates per setting.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A6/Supplement-Figure-A6-from-rds.R
```

This reads `Supplement-Figure-A6.rds` and recreates `Supplement-Figure-A6.pdf` and its TikZ source (`Supplement-Figure-A6.tex`). No model fitting is required.

## Rerun the simulation

First regenerate the shared target-only results:

```sh
Rscript figures/Supplement-Figure-A6/Supplement-Figure-A6-target-only.R
```

Then run:

```sh
Rscript figures/Supplement-Figure-A6/Supplement-Figure-A6.R
```

The simulation generates its data and updates `Supplement-Figure-A6.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The distance-to-set fits use the RcppEigen solver provided in `method/`; compiling it requires a working C++ toolchain.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A6-paired-comparisons.csv`.
