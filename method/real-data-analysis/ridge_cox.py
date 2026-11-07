"""Standardized ridge Cox fit with Breslow ties and optional case weights."""

import numpy as np, pandas as pd
from scipy.optimize import minimize


def fit(d, levels, weights=None):
    w = np.ones(len(d)) if weights is None else weights.astype(float)
    take = w > 0
    w = w[take]
    d = d.loc[take].reset_index(drop=True)
    x = np.column_stack(
        [(d.Sex == "M").astype(float)]
        + [(d.ProtDescription == m).astype(float) for m in levels]
    )
    mean = np.average(x, axis=0, weights=w)
    sd = np.sqrt(np.average((x - mean) ** 2, axis=0, weights=w))
    sd[sd < 1e-10] = 1
    x = (x - mean) / sd
    t = d.time.to_numpy()
    e = d.status.to_numpy() * w
    order = np.argsort(-t, kind="stable")
    x = x[order]
    t = t[order]
    e = e[order]
    w = w[order]
    ends = np.r_[np.flatnonzero(t[1:] != t[:-1]), len(t) - 1]
    starts = np.r_[0, ends[:-1] + 1]
    ev = np.add.reduceat(e, starts)
    keep = ev > 0
    ends = ends[keep]
    ev = ev[keep]
    eventx = e @ x
    N = w.sum()

    def fun(b):
        eta = x @ b
        shift = eta.max()
        z = w * np.exp(eta - shift)
        s0 = np.cumsum(z)[ends]
        s1 = np.cumsum(z[:, None] * x, axis=0)[ends]
        loss = (-eventx @ b + ev @ (np.log(s0) + shift)) / N + 1e-4 * (b @ b) / 2
        grad = (-eventx + (ev[:, None] * s1 / s0[:, None]).sum(axis=0)) / N + 1e-4 * b
        return loss, grad

    res = minimize(
        fun,
        np.zeros(x.shape[1]),
        jac=True,
        method="L-BFGS-B",
        options={"gtol": 1e-8, "ftol": 1e-12, "maxiter": 1000},
    )
    return res.x / sd, res, fun
