"""Helpers to read back what the DAG1 pipeline saved.

  from load import refit, dag
  proteins, B = refit("DAG1/runs/results_thr_l0.1_tA0.025_t0.1/amyA/1", ordering=0, step=5)   # B after 5 sinks removed (75x75)
  G = dag("DAG1/results/amyA/1", 17)                                # DAG #17 as a labelled DataFrame
"""
import numpy as np
import pandas as pd
from pathlib import Path


def refit(state_dir, ordering, step):
    """Refit matrix used at removal step `step` (0 = full 80x80) of ordering `ordering`.
    Returns (protein list, B[target, source] DataFrame)."""
    z = np.load(Path(state_dir) / "refits.npz")
    rid = z["step_refit"][ordering, step]
    sizes = z["sizes"]
    off = int((sizes[:rid] ** 2).sum())
    k = int(sizes[rid])
    B = z["flat"][off:off + k * k].reshape(k, k)
    prots = list(z["proteins"][z["masks"][rid]])
    return prots, pd.DataFrame(B, index=prots, columns=prots)


def dag(state_dir, m):
    z = np.load(Path(state_dir) / "dags.npz")
    P = list(z["proteins"])
    return pd.DataFrame(z["G"][m].astype(int), index=P, columns=P)   # [source, target]
