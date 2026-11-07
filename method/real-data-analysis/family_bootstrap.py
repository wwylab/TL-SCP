"""Family-bootstrap standard errors for mutation effects relative to R175H."""

import sys
import numpy as np
import pandas as pd
from ridge_cox import fit

rng = np.random.default_rng(20260915)
rows = []
for cohort, file in zip(["NCI", "MDA"], sys.argv[1:3]):
    d = pd.read_csv(file)
    if d.Family_ID.isna().any():
        raise ValueError("Family_ID is required for every carrier.")
    levels = sorted(d.ProtDescription.unique())
    reference = levels.index("p.R175H")
    coefficient, result, objective = fit(d, levels)
    point = coefficient[1:] - coefficient[reference + 1]
    families = sorted(d.Family_ID.astype(str).unique())
    codes = pd.Categorical(d.Family_ID.astype(str), categories=families).codes
    mutation_codes = pd.Categorical(d.ProtDescription, categories=levels).codes
    boot = np.full((1000, len(levels)), np.nan)
    for b in range(1000):
        family_counts = np.bincount(
            rng.integers(len(families), size=len(families)), minlength=len(families)
        )
        weights = family_counts[codes]
        counts = np.bincount(mutation_codes, weights=weights, minlength=len(levels))
        coef, result, objective = fit(d, levels, weights)
        if not np.isfinite(coef).all() or max(abs(objective(result.x)[1])) > 1e-5:
            continue
        if counts[reference] > 0:
            values = coef[1:] - coef[reference + 1]
            values[counts == 0] = np.nan
            boot[b] = values
        if (b + 1) % 100 == 0:
            print(cohort, b + 1, flush=True)
    for j, mutation in enumerate(levels):
        values = boot[np.isfinite(boot[:, j]), j]
        rows.append(
            dict(
                cohort=cohort,
                Mutation=mutation,
                stable_beta=point[j],
                bootstrap_SE=np.std(values, ddof=1),
                valid_replicates=len(values),
            )
        )
pd.DataFrame(rows).to_csv(sys.argv[3], index=False)
