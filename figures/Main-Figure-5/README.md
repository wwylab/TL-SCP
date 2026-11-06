# Main Figure 5: TL-SCP onset groups at MD Anderson

This figure compares cancer-free survival across YTA classes and TL-SCP onset groups in the MD Anderson cohort (202 carriers with 28 recurrent mutations). TL-SCP uses ranking information from 329 TP53 Database source carriers. Shaded bands show pointwise 90% confidence intervals.

## Data and variables

`Main-Figure-5.rds` contains `yta` and `tlscp` Kaplan–Meier summaries for the MD Anderson target cohort. The YTA panel uses classes 0 and A–D, excluding unclassified mutations; the TL-SCP panel uses the reported onset groups. The [shared survival-summary definitions](../../method/README.md#shared-clinical-results) describe the curves, confidence limits, counts, and log-rank results.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-5/Main-Figure-5-from-rds.R
```

This reads `Main-Figure-5.rds` and recreates `Main-Figure-5.pdf` and its TikZ source (`Main-Figure-5.tex`). No model fitting is required.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-5/Main-Figure-5.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

To refit the MD Anderson analysis, set `MDA_DATA` to an authorized CSV or the reviewed `MDAprospectiveMUT_clean.RData` file; see [MD Anderson inputs](../../method/README.md#md-anderson-inputs). The saved `.rds` results regenerate this figure without individual patient data.

Reported TL-SCP labels include the rank-based reassignment of R337H to Mid onset. The [clinical analysis documentation](../../method/README.md#real-data-analyses) gives the fitting and reporting rules.
