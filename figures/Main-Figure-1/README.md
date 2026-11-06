# Main Figure 1: Cancer-free survival by YTA class

This figure compares cancer-free survival between YTA classes 0 and A among 303 LFS/Chompret carriers. Despite their distinct functional annotations, the two classes show limited separation in age at first cancer. The curves include 95% confidence intervals.

## Data and variables

`Main-Figure-1.rds` contains the Kaplan–Meier curves, confidence limits, group counts, and log-rank result for target carriers in YTA classes 0 and A. The shared [survival-summary definitions](../../method/README.md#shared-clinical-results) describe the fields used to draw this figure.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-1/Main-Figure-1-from-rds.R
```

This reads `Main-Figure-1.rds` and recreates `Main-Figure-1.pdf` and its TikZ source (`Main-Figure-1.tex`). No model fitting is required.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-1/Main-Figure-1.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.
