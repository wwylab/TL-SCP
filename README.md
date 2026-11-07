# TL-SCP: Reproducibility Materials

This repository contains the code and data for the simulation studies and real-data analyses in *Transfer Learning for Survival-based Clustering of Predictors with an Application to TP53 Mutation Annotation*. Saved R results and scripts are also included to reproduce the figures and tables in the manuscript and Supplementary Materials.

## Repository structure

The files are organized into three directories:

- `figures/`: Main Figures 1–7 and Supplementary Figures A.1–A.7.
- `tables/`: Main Tables 1–2 and Supplementary Tables F.1–F.2.
- `method/`: SCP and TL-SCP implementations, simulation code, real-data analyses, and plotting functions.

Each figure and table has its own folder containing the analysis script, saved results, reproduction script, and README. For example, `figures/Main-Figure-3/` contains:

- `Main-Figure-3.R` — Runs the analysis and generates the figure.
- `Main-Figure-3.rds` — Saved numerical results.
- `Main-Figure-3-from-rds.R` — Recreates the figure from the saved results.
- `Main-Figure-3.pdf` and `Main-Figure-3.tex` — Final figure and TikZ source.
- `README.md` — Details of the analysis and instructions for reproduction.

The table folders follow the same structure, with LaTeX files for the table output. The shared code is in `method/`.

## Data and variables

The [method README](method/README.md#data-and-variables) describes the clinical variables, cohort preparation, and saved R results. It also defines the [simulation measures](method/README.md#simulation-results). The README in each figure or table folder gives the additional details needed for that analysis.

## Getting started

### 1. Install dependencies

Install the required R packages:

```r
install.packages(c(
  "Matrix", "glmnet", "survival", "dplyr", "ggplot2",
  "patchwork", "ggrepel", "tikzDevice", "Rcpp", "RcppEigen"
))
```

Figures that use TikZ require a TeX distribution with `pdflatex` and the `mathpazo`, `tikz`, and `preview` packages.

### 2. Reproduce figures and tables from saved results

The supplied `.rds` files allow the figures and tables to be regenerated without refitting the models. This also applies to the MD Anderson results, for which individual-level patient data cannot be shared.

From the repository root, run:

```sh
Rscript figures/Main-Figure-3/Main-Figure-3-from-rds.R
Rscript tables/Supplement-Table-F1/Supplement-Table-F1-from-rds.R
```

### 3. Rerun the analyses

To rerun an analysis, use the script without the `-from-rds` suffix in the corresponding figure or table folder. The README in that folder gives the required inputs and commands.

The simulation scripts generate their own data, fit the competing methods, and save the results. Most simulation settings use 50 replicates, as specified in the manuscript and Supplementary Materials.

The original model fits used R 4.5.3 and `glmnet` 4.1-10. To reproduce these estimates, install the same versions:

```r
install.packages("remotes")
remotes::install_version(
  "glmnet", version = "4.1-10", upgrade = "never"
)
```

Some analyses require additional software:

- **Gurobi:** Penalty comparisons in Supplementary Figures A.1–A.4.
- **C++ compiler:** Distance-to-set optimization implemented using RcppEigen.
- **Python:** DCSM comparisons, family bootstrap analyses, and alluvial figures. The alluvial renderer uses only the Python standard library.

The corresponding folder READMEs list the additional dependencies and instructions.

## Simulation studies

The simulations examine the performance of SCP and TL-SCP under different noise levels, censoring rates, problem dimensions, source sample sizes, covariate shifts, and source–target ranking differences. We also compare different fusion penalties and evaluate estimation accuracy, group recovery, and computation time.

Supplementary Figure A.7 examines how TL-SCP performs as the source and target rankings become increasingly different. The simulations consider six levels of ranking discordance and two source sample sizes. The reference lines indicate where the fitted population weighted-ranking margin changes sign.

## Real-data analyses

We apply SCP and TL-SCP to group germline *TP53* mutations according to their associations with age at first cancer.

### *TP53* Database

The primary analysis uses release 21 of the public *TP53* Database and a yeast transactivation (YTA) mutation catalog. Both are included in `tables/Main-Table-1/raw-data/`.

The target cohort includes carriers classified as LFS/Chompret, while the source cohort includes carriers with LFL, FH, noFH, and Other family classifications. The mutation selection and filtering rules are described in the manuscript and implemented in the analysis scripts.

### MD Anderson cohort

We also analyze data from a prospective family-based Li-Fraumeni syndrome cohort at MD Anderson Cancer Center. Individual-level clinical and family data cannot be publicly shared.

The repository includes saved R results containing mutation annotations, carrier counts, survival summaries, and coefficient estimates. These files can be used to reproduce the reported figures and tables without access to the original patient records.

We also include a synthetic dataset with the same sample size and variable structure as the MD Anderson target cohort. Together with the public *TP53* Database source cohort, it allows readers to run the SCP and TL-SCP analyses without access to the restricted patient records. The resulting estimates are based on simulated data and will differ from those reported in the manuscript.

See the [synthetic-data instructions](method/README.md#synthetic-md-anderson-data) for the files and commands.

Researchers with authorized access can rerun the analyses by setting `MDA_DATA` to the appropriate input file, as described in the [method README](method/README.md#md-anderson-inputs).

For this analysis, the fitted TL-SCP model has four groups. The R337H singleton is subsequently assigned to the Mid-onset group according to its weighted rank, as reported in the manuscript.

## Code availability

SCP and TL-SCP were developed by Wenyi Wang's Statistical Bioinformatics Lab and Xiaoqian Liu.

The original TL-SCP implementation is available at:

https://github.com/Xiaoqian-Liu/TL-SCP

The code is distributed under the MIT license. External software and public datasets remain subject to their respective licenses and terms.

The DCSM comparison uses the implementation available at:

https://github.com/BojianHou/DCSM
