"""Sampler noise vs data noise for the sweep winner (reliable states, n>=150).

  seed agreement  : consensus edge frequency / mean ordering position, winner (seed A) vs seed 111 vs 222
                    -> same people, same B; only the random sampling differs
  N=300 vs N=100  : does more orderings change the consensus?
  best DAG        : Jaccard of the selected best DAG across seeds
Compare with split-half agreement from sweep/stage2 (different people) = data noise.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

H = Path(__file__).resolve().parent            # DAG1/parameter_sweep
RUNS = H.parent / "runs"                       # where run_all.py / stage3.sh write
DET = H / "results" / "details"
runs = {"seedA": "results_sweep_winner", "seed111": "results_sweep_winner_seed111",
        "seed222": "results_sweep_winner_seed222", "N300": "results_sweep_winner_N300"}
iu = np.triu_indices(80, 1)


def load(run, v, s):
    d = RUNS / runs[run] / v / str(s)
    z = np.load(d / "dags.npz"); G = z["G"]; F = G.mean(0)
    O = pd.read_csv(d / "orderings.csv", index_col=0).values
    P = list(z["proteins"]); ix = {p: i for i, p in enumerate(P)}
    pos = np.zeros(80)
    for o in O:
        for k, p in enumerate(o): pos[ix[p]] += k
    best = G[int(np.argmax(z["scores"]))]
    return (F + F.T)[iu], pos / len(O), best


S = pd.read_csv(RUNS / runs["seed111"] / "summary.csv", dtype={"state": str})
rows = []
for _, r in S.iterrows():
    L = {k: load(k, r.variable, r.state) for k in runs}
    def rep(a, b, i): return spearmanr(L[a][i], L[b][i])[0]
    def jac(a, b):
        x, y = L[a][2], L[b][2]; return (x & y).sum() / (x | y).sum()
    rows.append(dict(variable=r.variable, state=r.state, n=r.n,
                     seed_edges=np.mean([rep("seedA", "seed111", 0), rep("seedA", "seed222", 0), rep("seed111", "seed222", 0)]),
                     seed_pos=np.mean([rep("seedA", "seed111", 1), rep("seedA", "seed222", 1), rep("seed111", "seed222", 1)]),
                     N300_vs_100_edges=rep("seedA", "N300", 0), N300_vs_100_pos=rep("seedA", "N300", 1),
                     bestDAG_jaccard_seeds=np.mean([jac("seedA", "seed111"), jac("seedA", "seed222"), jac("seed111", "seed222")])))
D = pd.DataFrame(rows)
D.to_csv(DET / "stability_winner.csv", index=False)
sp = pd.read_csv(DET / "stage2_per_state.csv")
w = sp[sp.config.str.contains("LAMBDA=0.12") & sp.config.str.contains("TAU_A_RATIO=0.5") &
       sp.config.str.contains("T_SOFTMAX=0.05") & sp.config.str.contains("SHRINK_SCALE=1.0") &
       sp.config.str.contains("CONTEXT=none") & sp.config.str.contains("LAMBDA_MODE=abs")]
print("winner, reliable states (median):")
print(f"  sampler noise  (same people, different seeds): edges {D.seed_edges.median():.3f}   ordering {D.seed_pos.median():.3f}")
print(f"  more orderings (N=300 vs 100):                 edges {D.N300_vs_100_edges.median():.3f}   ordering {D.N300_vs_100_pos.median():.3f}")
print(f"  data noise     (different halves of people):   edges {w.rep_edges.median():.3f}   ordering {w.rep_pos.median():.3f}")
print(f"  best-DAG edge overlap across seeds (Jaccard):  {D.bestDAG_jaccard_seeds.median():.3f}")
