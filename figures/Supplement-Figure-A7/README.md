# Supplementary Figure A.7: Transfer under ranking disagreement

This simulation examines how TL-SCP performs as source and target effect rankings become less consistent. It considers six levels of population discordance (0–0.40), two source sample sizes (500 and 1,000), and 50 replicates per setting. The panels relate empirical ranking disagreement to estimation error, clustering accuracy, and the population weighted ranking margin.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Supplement-Figure-A7/Supplement-Figure-A7-from-rds.R
```

This reads `Supplement-Figure-A7.rds` and recreates `Supplement-Figure-A7.pdf` and its TikZ source (`Supplement-Figure-A7.tex`). No model fitting is required.

## Rerun the simulation

First regenerate the shared target-only results:

```sh
Rscript figures/Supplement-Figure-A7/Supplement-Figure-A7-target-only.R
```

Then run:

```sh
Rscript figures/Supplement-Figure-A7/Supplement-Figure-A7.R
```

The simulation generates its data and updates `Supplement-Figure-A7.rds`. To preserve the published numerical results, make a copy of the saved file before rerunning.

The distance-to-set fits use the RcppEigen solver provided in `method/`; compiling it requires a working C++ toolchain.

Population discordance counts cross-group rank reversals; empirical discordance is computed from preliminary Cox estimates. Within each replicate, source ranks are perturbed to the designated population discordance. The plot uses LOESS (span 0.8); its vertical reference lines mark fitted zero crossings of the population weighted ranking margin (approximately 0.249 and 0.227). All 600 transfer fits enter the reported summaries, including those that reached the iteration limit.

To check the installation with one fit before running the full simulation, use `TLSCP_LIMIT=1 TLSCP_REPS=1` with the simulation command. This overwrites the saved results, so run the check in a separate copy of the repository.
