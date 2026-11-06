# Main Figure 7: Comparing mutation rankings across cohorts

The first two panels compare weighted mutation rankings from the TP53 Database and MD Anderson analyses. The third compares preliminary Cox effects, relative to R175H, for mutations assigned to different onset groups. Standard errors are based on 1,000 family-bootstrap samples per cohort, with the two variances combined for the intervals.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-7/Main-Figure-7-from-rds.R
```

This reads `Main-Figure-7.rds` and recreates `Main-Figure-7.pdf` and its TikZ source (`Main-Figure-7.tex`). No model fitting is required.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-7/Main-Figure-7.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

**MD Anderson data:** Individual-level clinical and family records are not distributed. To rerun the clinical analysis, set `MDA_DATA` to an authorized analysis CSV with `Individual_code`, `Sex` (`M`/`F`), `ProtDescription`, `time`, and `status` (1 = cancer, 0 = censored). The preparation code also accepts the reviewed `MDAprospectiveMUT_clean.RData` file. The saved `.rds` results can be used without access to these records.

**Reported onset groups:** The MD Anderson analysis fits four TL-SCP groups and then assigns the R337H singleton to Mid onset using its weighted rank. This reporting step is explicit in the analysis code; `raw_TLSCP_group` preserves the group assigned by the fitted model. The SCP comparison retains its four fitted groups.

The family bootstrap additionally requires `Family_ID` in the MD Anderson analysis data and Python with `numpy`, `pandas`, and `scipy`. Set `PYTHON` if these are installed in a separate environment.
