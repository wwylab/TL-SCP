# Main Figure 6: MD Anderson mutation annotations

The alluvial plot compares YTA classes with TL-SCP onset groups for 28 mutations in the MD Anderson cohort. The reported assignments contain 7 Early-, 16 Mid-, and 5 Late-onset mutations, with each mutation counted once.

## Data and variables

`Main-Figure-6.rds` has one row per MD Anderson analysis mutation. The plot compares `YTAClass` with the reported `TL-SCP` label for each `ProtDescription`. The file also contains `SCP` and carrier counts (`target_n`, `source_n`); ribbon widths count mutations. See [Data and variables](../../method/README.md#data-and-variables) for shared definitions.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-6/Main-Figure-6-from-rds.R
```

This reads `Main-Figure-6.rds` and recreates `Main-Figure-6.pdf` and its TikZ source (`Main-Figure-6.tex`). No model fitting is required.

The alluvial figure uses Python 3 and LaTeX; no Python packages are needed for this step. Set `PYTHON` if your Python executable has a different name.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-6/Main-Figure-6.R
```

The public *TP53* Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

To refit the MD Anderson analysis, set `MDA_DATA` to an authorized CSV or the reviewed `MDAprospectiveMUT_clean.RData` file; see [MD Anderson inputs](../../method/README.md#md-anderson-inputs). The saved `.rds` results regenerate this figure without individual patient data.

To test this script with the synthetic target cohort and public source data, follow the [synthetic-data instructions](../../method/README.md#synthetic-md-anderson-data), which use a temporary copy to preserve the published results.

Reported TL-SCP labels include the rank-based reassignment of R337H to Mid onset. The [clinical analysis documentation](../../method/README.md#real-data-analyses) gives the fitting and reporting rules.
