# TL-SCP: Reproducible Code and Analyses

This repository provides the code and results for the manuscript *Transfer Learning for Survival-based Clustering of Predictors with an Application to TP53 Mutation Annotation*.

Survival-based Clustering of Predictors (SCP) groups germline *TP53* mutations with similar effects on age at first cancer. Transfer learning-enhanced SCP (TL-SCP) uses mutation rankings from a related source cohort to form these groups and estimates their effects in the target cohort.

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

## Data and variables

The [method README](method/README.md#data-and-variables) describes the clinical inputs, cohort preparation, and shared result fields, along with the [simulation measures](method/README.md#simulation-results). Each figure and table README describes its saved `.rds` files and gives the commands for regenerating the output or rerunning the analysis.

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

The saved numerical results (`.rds` files) can be used to regenerate the figures and tables without refitting models or accessing restricted clinical data.

For example, from the repository root:

```sh
Rscript figures/Main-Figure-3/Main-Figure-3-from-rds.R
Rscript tables/Supplement-Table-F1/Supplement-Table-F1-from-rds.R
```

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
- **Python:** DCSM comparisons, family bootstrap analyses, and alluvial figures. The alluvial renderer uses only the Python standard library.

The corresponding folder READMEs list the additional dependencies and instructions.

## Simulation studies

The simulations compare SCP and TL-SCP across noise levels, censoring rates, problem dimensions, source sample sizes, covariate shifts, and disagreement in mutation-effect rankings. They also compare fusion penalties and measure estimation error and mutation-group recovery.

Supplementary Figure A.7 uses six levels of population ranking discordance and two source sample sizes to examine estimation error and clustering accuracy as disagreement between source and target rankings increases. Vertical reference lines mark where the fitted population weighted ranking margin crosses zero.

## Real data analysis

The manuscript applies SCP and TL-SCP to annotate germline *TP53* mutations according to their associations with age at first cancer.

### *TP53* Database

The public *TP53* Database (release 21) and YTA mutation catalog are provided in `tables/Main-Table-1/raw-data/`.

The target cohort consists of carriers classified as LFS/Chompret; the source cohort includes carriers with LFL, FH, noFH, and Other family classifications. The analysis scripts apply the mutation selection and filtering rules described in the manuscript.

### MD Anderson cohort

The MD Anderson analyses use individual-level clinical and family data that cannot be publicly distributed.

The repository includes aggregate mutation annotations, carrier counts, survival summaries, and estimated coefficients. These files do not contain patient identifiers, family identifiers, or individual-level survival records.

Researchers with authorized access can rerun the analyses by setting `MDA_DATA` to a file in one of the [documented input formats](method/README.md#md-anderson-inputs). The figure and table READMEs give the commands.

The MD Anderson analysis reassigns R337H to the Mid-onset group using its weighted rank, as described in the manuscript. The [method README](method/README.md#real-data-analyses) explains how the scripts apply this reporting step.

## Code availability

SCP and TL-SCP were developed by Wenyi Wang's Statistical Bioinformatics Lab and Xiaoqian Liu.

The original TL-SCP implementation is available at:

https://github.com/Xiaoqian-Liu/TL-SCP

The code is distributed under the MIT license. External software and public datasets remain subject to their respective licenses and terms.

The DCSM comparison uses the implementation available at:

https://github.com/BojianHou/DCSM
