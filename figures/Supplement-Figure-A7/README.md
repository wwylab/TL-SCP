# Supplementary Figure A.7: Transfer under ranking disagreement

This simulation evaluates TL-SCP as source and target effect rankings diverge. It uses six levels of population discordance (0 to 0.40), two source sample sizes (500 and 1,000), and 50 replicates per setting. The two columns show estimation error and clustering accuracy against empirical ranking disagreement. Vertical reference lines mark where the fitted population weighted ranking margin crosses zero.

## Data and variables

`Supplement-Figure-A7.rds` contains two data frames joined by source sample size `n_S`, designated population discordance `D_target` (0, 0.05, 0.10, 0.20, 0.30, or 0.40), and `rep`. The `metrics` frame holds `relErr` and `NMI` by `method`; `diagnostics` holds realized population discordance `D_realized`, empirical discordance `D_hat`, and the population weighted ranking margin `m_TL`. The panels plot the two performance measures against `D_hat`. The fitted `m_TL` curve determines the vertical zero-crossing reference lines.

`Supplement-Figure-A7-target-only.rds` contains the SCP and Cox K-means results reused across source settings. Its `metrics` frame holds those summaries, while `beta_T` and `X_T` store the simulated target coefficients and design matrix. See the [shared simulation variable definitions](../../method/README.md#simulation-results).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A7/Supplement-Figure-A7-from-rds.R
```

The script reads `Supplement-Figure-A7.rds` and regenerates `Supplement-Figure-A7.pdf` and its TikZ source (`Supplement-Figure-A7.tex`) without refitting the models.

## Rerun the simulation

First regenerate the shared target-only results:

```sh
Rscript figures/Supplement-Figure-A7/Supplement-Figure-A7-target-only.R
```

Then run:

```sh
Rscript figures/Supplement-Figure-A7/Supplement-Figure-A7.R
```

The script generates simulated data and overwrites `Supplement-Figure-A7.rds`. Copy the saved file before rerunning to retain the reported results.

The distance-to-set fits use the RcppEigen solver provided in `method/`; compiling it requires a working C++ toolchain.

Population discordance counts cross-group rank reversals; empirical discordance is computed from preliminary Cox estimates. Within each replicate, source ranks are perturbed to the designated population discordance. The plot uses LOESS (span 0.8); its vertical reference lines mark fitted zero crossings of the population weighted ranking margin (approximately 0.249 and 0.227). The summaries use all 600 transfer fits, including those that reached the iteration limit.

For a one-fit installation check, use `TLSCP_LIMIT=1 TLSCP_REPS=1` with the simulation command. This overwrites the saved results, so run the check in a separate copy of the repository.
