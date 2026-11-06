# TL-SCP: Reproducible Code and Analyses

This repository provides the code and results for the manuscript *Transfer Learning for Survival-based Clustering of Predictors with an Application to TP53 Mutation Annotation*.

Germline *TP53* mutations vary in their associations with age at first cancer, but many are observed in only a small number of carriers, making their individual effects difficult to estimate. Survival-based Clustering of Predictors (SCP) addresses this problem by grouping mutations with similar survival effects. Transfer learning-enhanced SCP (TL-SCP) further improves mutation grouping by incorporating ranking information from a related source cohort, while estimating the grouped effects using target-cohort data.

The repository includes the simulation studies, real-data analyses, and code used to produce the figures and tables in the manuscript and Supplementary Materials.

## Repository structure

The repository is organized into three directories:

- **`figures/`** — Main Figures 1–7 and Supplementary Figures A.1–A.7.
- **`tables/`** — Main Tables 1–2 and Supplementary Tables F.1–F.2.
- **`method/`** — Methodology, simulation, real-data analysis, and visualization code.

Each figure and table has its own folder, named according to its numbering in the manuscript. For example, `figures/Main-Figure-3/` contains:

- `Main-Figure-3.R` — Runs the analysis and generates the figure.
- `Main-Figure-3.rds` — Saved numerical results.
- `Main-Figure-3-from-rds.R` — Recreates the figure from the saved results.
- `Main-Figure-3.pdf` and `Main-Figure-3.tex` — Final figure and TikZ source.
- `README.md` — Details of the analysis and instructions for reproduction.

Tables follow the same organization, with LaTeX files containing the table output. The `method/` directory contains the core implementations of SCP and TL-SCP, along with supporting functions.

## Getting started

### 1. Install dependencies

The analyses were conducted in R. Install the required packages using:

```r
install.packages(c(
  "Matrix", "glmnet", "survival", "dplyr", "ggplot2",
  "patchwork", "ggrepel", "tikzDevice", "Rcpp", "RcppEigen"
))
```

Figures that use TikZ require a TeX distribution with `pdflatex` and the `mathpazo`, `tikz`, and `preview` packages.

### 2. Reproduce figures and tables from saved results

The repository includes saved numerical results (`.rds` files), allowing figures and tables to be reproduced without rerunning the analyses.

For example, from the repository root:

```sh
Rscript figures/Main-Figure-3/Main-Figure-3-from-rds.R
Rscript tables/Supplement-Table-F1/Supplement-Table-F1-from-rds.R
```

These scripts regenerate the corresponding outputs directly from the saved results. No model fitting or access to restricted clinical data is required.

### 3. Rerun the analyses

Each figure and table folder contains an analysis script and a README with instructions for running it.

The simulation scripts generate the simulated datasets, fit the competing methods, and produce the numerical results. The experiments generally use 50 simulation replicates per setting, as described in the manuscript and Supplementary Materials.

The original analyses used R 4.5.3 and `glmnet` 4.1-10. To reproduce the original model estimates, we recommend installing the same `glmnet` version:

```r
install.packages("remotes")
remotes::install_version(
  "glmnet", version = "4.1-10", upgrade = "never"
)
```

Some analyses require additional software:

- **Gurobi:** Penalty comparisons in Supplementary Figures A.1–A.4.
- **C++ compiler:** Distance-to-set optimization implemented using RcppEigen.
- **Python:** DCSM comparisons and family bootstrap analyses.

Instructions and additional dependencies are provided in the corresponding folder READMEs.

## Simulation studies

The simulation studies examine the performance of SCP and TL-SCP under different settings, including varying noise levels, censoring rates, problem dimensions, source sample sizes, covariate shifts, and inconsistencies in mutation-effect rankings.

The studies also compare alternative fusion penalties and assess when transfer learning improves estimation and mutation grouping.

In particular, Supplementary Figure A.7 examines the effect of increasing disagreement between source and target rankings. The experiment considers six levels of population ranking discordance and two source sample sizes, relating empirical ranking disagreement to estimation error, clustering accuracy, and the population weighted ranking margin.

All simulation scripts and saved results are organized by figure or table number.

## Real Data Analysis

The manuscript applies SCP and TL-SCP to annotate germline *TP53* mutations according to their associations with age at first cancer.

### *TP53* Database

The public *TP53* Database (release 21) and YTA mutation catalog are provided in `tables/Main-Table-1/raw-data/`.

The target cohort consists of carriers classified as LFS/Chompret, while the source cohort includes carriers with LFL, FH, noFH, and Other family classifications. The analysis scripts implement the mutation selection and filtering procedures described in the manuscript.

### MD Anderson cohort

The MD Anderson analyses use individual-level clinical and family data that cannot be publicly distributed.

To support reproducibility, the repository provides saved aggregate results, including mutation annotations, carrier counts, survival summaries, and estimated coefficients. These files do not contain patient identifiers, family identifiers, or individual-level survival records.

Researchers with authorized access to the clinical data can rerun the analyses by specifying the `MDA_DATA` environment variable, following the instructions in the relevant folders.

The MD Anderson analysis also includes a rank-based reassignment of R337H to the Mid-onset group, as described in the manuscript. This application-specific step is implemented explicitly in the corresponding scripts.

## Code availability

SCP and TL-SCP were developed by Wenyi Wang's Statistical Bioinformatics Lab and Xiaoqian Liu.

The original TL-SCP implementation is available at:

https://github.com/Xiaoqian-Liu/TL-SCP

The code is distributed under the MIT license. External software and public datasets remain subject to their respective licenses and terms.

The DCSM comparison uses the implementation available at:

https://github.com/BojianHou/DCSM