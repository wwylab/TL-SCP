# Supplementary Table F.2: Mutation annotations at MD Anderson

This table lists 28 MD Anderson mutations, source-cohort carrier counts, and the corresponding onset annotations. The reported TL-SCP groups include the rank-based reassignment of R337H.

## Reproduce from saved results

From the repository root:

```sh
Rscript tables/Supplement-Table-F2/Supplement-Table-F2-from-rds.R
```

This reads `Supplement-Table-F2.rds` and recreates `Supplement-Table-F2.tex` as a LaTeX table (`booktabs` required). No model fitting is required.

## Rerun the analysis

```sh
Rscript tables/Supplement-Table-F2/Supplement-Table-F2.R
```

**MD Anderson data:** Individual-level clinical and family records are not distributed. To rerun the clinical analysis, set `MDA_DATA` to an authorized analysis CSV with `Individual_code`, `Sex` (`M`/`F`), `ProtDescription`, `time`, and `status` (1 = cancer, 0 = censored). The preparation code also accepts the reviewed `MDAprospectiveMUT_clean.RData` file. The saved `.rds` results can be used without access to these records.

**Reported onset groups:** The MD Anderson analysis fits four TL-SCP groups and then assigns the R337H singleton to Mid onset using its weighted rank. This reporting step is explicit in the analysis code; `raw_TLSCP_group` preserves the group assigned by the fitted model. The SCP comparison retains its four fitted groups.
