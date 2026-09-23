"""State-specific node-wise SEM:  X = B X + Gamma Z + eps,  diag(B)=0.
B is returned as B[target, source]  (= effect of source on target, i.e. source -> target).

Solver "lw" (default, fast, exact for the model with unpenalised context Z):
  1. residualise every protein on the context (intercept + non-constant Z + dataset dummies).
     By Frisch-Waugh-Lovell, the protein->protein coefficients of the joint regression
     equal the node-wise regressions among these residuals.
  2. node-wise regressions among residuals via the precision matrix of the Ledoit-Wolf
     shrunk covariance:  B_ij = -Theta_ij / Theta_ii.
  Refitting after removing proteins re-estimates the shrinkage AND re-inverts the smaller
  covariance -> a genuine new fit of the reduced system (not a row/col deletion of the old B).
  Speed trick: S = R'R/n and M = (R^2)'(R^2)/n are computed once per state; every refit only
  needs their sub-blocks + one small matrix inverse.

Solver "enet": cross-validated elastic net per node on [other proteins, context] — the solver
  used for the original 80x80 general graph. Exact same model, but ~10^4x slower.
"""
import numpy as np
import pandas as pd
from config import MIN_RESID_DF


class StateSystem:
    def __init__(self, X, Znum, datasets, solver="lw"):
        """X: n x P DataFrame (expr_z, NaN allowed); Znum: n x q numeric context; datasets: n labels."""
        self.proteins = list(X.columns)
        self.n = len(X)
        self.solver = solver
        Xf = X.fillna(0.0).values.astype(float)
        # ---- context: intercept + non-constant Z (NaN -> state mean) + non-constant dataset dummies
        import config as cfg
        self.shrink_scale = cfg.SHRINK_SCALE
        C = []
        self.context_cols = []
        zcols = [] if cfg.CONTEXT == "none" else [c for c in Znum.columns
                                                   if not (cfg.CONTEXT == "no_pmi" and c == "pmi")]
        for c in zcols:
            z = Znum[c].astype(float)
            if z.notna().sum() < 2:
                continue
            z = z.fillna(z.mean())
            if z.std() > 1e-12:
                C.append(z.values); self.context_cols.append(c)
        for d in sorted(set(datasets))[1:]:
            dd = (np.asarray(datasets) == d).astype(float)
            if 0 < dd.sum() < len(dd):
                C.append(dd); self.context_cols.append("ds_" + d)
        self.z_dropped = False
        if self.n - len(C) - 1 < MIN_RESID_DF:
            C, self.context_cols, self.z_dropped = [], [], True
        Cm = np.column_stack([np.ones(self.n)] + C) if C else np.ones((self.n, 1))
        self.Cm, self.Xf = Cm, Xf
        coef, *_ = np.linalg.lstsq(Cm, Xf, rcond=None)
        R = Xf - Cm @ coef
        R -= R.mean(0)
        self.S = R.T @ R / self.n
        R2 = R * R
        self.M = R2.T @ R2 / self.n

    # ------------------------------------------------------------------
    def fit(self, idx):
        """B over the protein subset idx (list of column indices), shape k x k, [target, source]."""
        k = len(idx)
        if k == 1:
            return np.zeros((1, 1))
        if self.solver == "enet":
            return self._fit_enet(idx)
        ix = np.ix_(idx, idx)
        S, M, n = self.S[ix], self.M[ix], self.n
        mu = np.trace(S) / k
        S2 = S * S
        delta = (S2.sum() - 2 * mu * np.trace(S) + k * mu * mu) / k
        beta = min((M.sum() - S2.sum()) / (k * n), delta)
        shrink = 0.0 if beta == 0 else beta / delta
        shrink = min(max(shrink * self.shrink_scale, 0.0), 1.0)
        Sig = (1 - shrink) * S
        Sig[np.diag_indices(k)] += shrink * mu
        Th = np.linalg.inv(Sig)
        B = -Th / np.diag(Th)[:, None]
        np.fill_diagonal(B, 0.0)
        return B

    def _fit_enet(self, idx):
        from sklearn.linear_model import ElasticNetCV
        k = len(idx)
        B = np.zeros((k, k))
        Xs = self.Xf[:, idx]
        for a in range(k):
            others = [b for b in range(k) if b != a]
            D = np.hstack([Xs[:, others], self.Cm[:, 1:]])
            m = ElasticNetCV(l1_ratio=[.5, .8], n_alphas=40, cv=3, max_iter=5000).fit(D, Xs[:, a])
            B[a, others] = m.coef_[:len(others)]
        return B
