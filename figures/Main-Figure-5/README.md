# Main Figure 5: TL-SCP onset groups at MD Anderson

This figure compares cancer-free survival across YTA classes and TL-SCP onset groups in the MD Anderson cohort (202 carriers with 28 recurrent mutations). TL-SCP uses ranking information from 329 TP53 Database source carriers. Shaded bands show pointwise 90% confidence intervals.

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

**MD Anderson data:** Individual-level clinical and family records are not distributed. To rerun the clinical analysis, set `MDA_DATA` to an authorized analysis CSV with `Individual_code`, `Sex` (`M`/`F`), `ProtDescription`, `time`, and `status` (1 = cancer, 0 = censored). The preparation code also accepts the reviewed `MDAprospectiveMUT_clean.RData` file. The saved `.rds` results can be used without access to these records.

**Reported onset groups:** The MD Anderson analysis fits four TL-SCP groups and then assigns the R337H singleton to Mid onset using its weighted rank. This reporting step is explicit in the analysis code; `raw_TLSCP_group` preserves the group assigned by the fitted model. The SCP comparison retains its four fitted groups.
