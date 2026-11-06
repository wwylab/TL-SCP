# Main Table 2: Agreement between cohort-specific TL-SCP annotations

This table compares TL-SCP annotations for the 17 mutations shared by the primary TP53 Database and MD Anderson analyses. Nine mutations have the same onset assignment; the other eight differ by one adjacent group. Carrier counts for discordant mutations appear in parentheses.

## Reproduce from saved results

From the repository root:

```sh
Rscript tables/Main-Table-2/Main-Table-2-from-rds.R
```

This reads `Main-Table-2.rds` and recreates `Main-Table-2.tex` as a LaTeX table (`booktabs` required). No model fitting is required.

## Rerun the analysis

```sh
Rscript tables/Main-Table-2/Main-Table-2.R
```

**MD Anderson data:** Individual-level clinical and family records are not distributed. To rerun the clinical analysis, set `MDA_DATA` to an authorized analysis CSV with `Individual_code`, `Sex` (`M`/`F`), `ProtDescription`, `time`, and `status` (1 = cancer, 0 = censored). The preparation code also accepts the reviewed `MDAprospectiveMUT_clean.RData` file. The saved `.rds` results can be used without access to these records.
