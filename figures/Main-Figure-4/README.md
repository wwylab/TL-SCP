# Main Figure 4: Comparing TL-SCP, SCP, and YTA annotations

These alluvial plots compare the TL-SCP onset groups with YTA classes and target-only SCP for 31 recurrent mutations. Each ribbon represents a mutation, showing how source-informed ranking changes its annotation.

## Reproduce from saved results

From the repository root:

```sh
Rscript figures/Main-Figure-4/Main-Figure-4-from-rds.R
```

This reads `Main-Figure-4.rds` and recreates `Main-Figure-4.pdf` and its TikZ source (`Main-Figure-4.tex`). No model fitting is required.

The alluvial figure uses Python 3 and LaTeX; no Python packages are needed for this step. Set `PYTHON` if your Python executable has a different name.

## Rerun the analysis

```sh
Rscript figures/Main-Figure-4/Main-Figure-4.R
```

The public TP53 Database release 21 and YTA catalog are provided in `tables/Main-Table-1/raw-data/`. The cohort definitions and mutation filters follow Section 5.1 of the manuscript.
