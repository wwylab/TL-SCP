# Supplementary Figure A.6: Transfer under covariate shift

This simulation introduces Gaussian noise into the source covariates while preserving the ordering of mutation effects. It compares TL-SCP with target-only SCP and Cox K-means for source sample sizes of 500 and 1,000, using 50 replicates per setting.

## Data and variables

`Supplement-Figure-A6.rds` is a data frame with simulation results by source sample size `n_S`, perturbation standard deviation `setting` (0.05, 0.10, 0.25, 0.35, or 0.50), `rep`, and `method`. The panels plot `relErr` and `NMI` for `Cox-K-means`, `SCP`, and `TL-SCP` separately for each source size.

`Supplement-Figure-A6-target-only.rds` contains the target-only results reused across source settings. Its `metrics` data frame holds the SCP and Cox K-means summaries; `beta_T` and `X_T` store the simulated target coefficients and design matrix. See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A6/Supplement-Figure-A6-from-rds.R
```

The script reads `Supplement-Figure-A6.rds` and regenerates `Supplement-Figure-A6.pdf` and its TikZ source (`Supplement-Figure-A6.tex`) without refitting the models.

## Rerun the simulation

First regenerate the shared target-only results:

```sh
Rscript figures/Supplement-Figure-A6/Supplement-Figure-A6-target-only.R
```

Then run:

```sh
Rscript figures/Supplement-Figure-A6/Supplement-Figure-A6.R
```

The script generates simulated data and overwrites `Supplement-Figure-A6.rds`. Copy the saved file before rerunning to retain the reported results.

The distance-to-set fits use the RcppEigen solver provided in `method/`; compiling it requires a working C++ toolchain.

The figure uses paired, two-sided Wilcoxon tests within each setting and metric, with Holm adjustment where multiple comparisons are made. Test results are saved in `Supplement-Figure-A6-paired-comparisons.csv`.
