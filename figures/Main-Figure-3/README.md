# Main Figure 3: TL-SCP onset groups in the TP53 Database

TL-SCP groups 31 recurrent mutations among 393 LFS/Chompret carriers using ranking information from 356 source-cohort carriers. The figure compares cancer-free survival across the three TL-SCP onset groups and YTA classes. Shaded bands are pointwise 90% confidence intervals, and the legends report mutation and carrier counts.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-3/Main-Figure-3-from-rds.R
```

This reads `Main-Figure-3.rds` and recreates `Main-Figure-3.pdf` and its TikZ source (`Main-Figure-3.tex`). No model fitting is required.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-3/Main-Figure-3.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.
