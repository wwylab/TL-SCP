# Main Figure 2: Mutation annotations from DCSM and Cox K-means

Panel (a) compares YTA classes with DCSM patient clusters, assigning each mutation to its most frequent cluster (or marking ties). Panel (b) compares sex-adjusted Cox K-means annotations between the source and target cohorts; 16 of the 31 mutations change groups.

## Data and variables

`Main-Figure-2.rds` contains two mutation-level data frames:

- `dcsm` holds carrier counts in the `Early-onset`, `Mid-onset`, and `Late-onset` patient clusters. Panel (a) compares `YTAClass` with `majority_DCSM`, the most frequent cluster or `Tied`.
- `cox` pairs source and target results by `Mutation`. `beta_source` and `beta_target` are preliminary Cox log-hazard coefficients; `group_source` and `group_target` give the assignments used in panel (b).

Each ribbon counts mutations. Shared input and annotation definitions are in [Data and variables](../../method/README.md#data-and-variables).

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-2/Main-Figure-2-from-rds.R
```

This reads `Main-Figure-2.rds` and recreates `Main-Figure-2.pdf` and its TikZ source (`Main-Figure-2.tex`). No model fitting is required.

The alluvial figure uses Python 3 and LaTeX; no Python packages are needed for this step. Set `PYTHON` if your Python executable has a different name.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-2/Main-Figure-2.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

The DCSM comparison also requires Python (`numpy`, `pandas`, `scipy`, `scikit-learn`, and `torch`) and a checkout of [DCSM](https://github.com/BojianHou/DCSM). Set `DCSM_DIR` to the checkout and, if necessary, `PYTHON` to the intended Python executable. The original reproduction check used commit `6d67dc7ae17be224b8ac3fb7526ba213b3d74f7b`. The model uses three Weibull experts, a 16-unit embedding, and five initializations, with model selection separated by family. Only mutation-level cluster counts are saved.
