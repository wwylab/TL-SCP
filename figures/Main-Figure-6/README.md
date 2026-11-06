# Main Figure 6: MD Anderson mutation annotations

The alluvial plot compares YTA classes with TL-SCP onset groups for 28 mutations in the MD Anderson cohort. The reported assignments contain 7 Early-, 16 Mid-, and 5 Late-onset mutations, with each mutation counted once.

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

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.

**MD Anderson data:** Individual-level clinical and family records are not distributed. To rerun the clinical analysis, set `MDA_DATA` to an authorized analysis CSV with `Individual_code`, `Sex` (`M`/`F`), `ProtDescription`, `time`, and `status` (1 = cancer, 0 = censored). The preparation code also accepts the reviewed `MDAprospectiveMUT_clean.RData` file. The saved `.rds` results can be used without access to these records.

**Reported onset groups:** The MD Anderson analysis fits four TL-SCP groups and then assigns the R337H singleton to Mid onset using its weighted rank. This reporting step is explicit in the analysis code; `raw_TLSCP_group` preserves the group assigned by the fitted model. The SCP comparison retains its four fitted groups.
