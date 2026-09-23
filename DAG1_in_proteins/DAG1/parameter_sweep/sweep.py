"""Parameter sweep over the DAG1 pipeline (architecture unchanged) with a C-FREE selection criterion.

Cleanliness = reproducibility. Each reliable state (n>=150) is split once into two random halves
(same split for every configuration). The full pipeline runs on each half; we measure
  rep_edges : Spearman between the two halves' consensus edge frequencies (over 3160 pairs)
  rep_pos   : Spearman between the two halves' mean ordering positions (80 proteins)
  rep_top   : Jaccard overlap of the two halves' top-5% most frequent edges
Pre-registered selection rule (fixed before any results):
  score = mean(rep_edges, rep_pos); keep only configs with median sink entropy in [0.2, 0.95]
  (not uniform, not a greedy/deterministic ordering); pick the max.
Agreement with the literature matrix C is computed and reported but NEVER used for selection.

  cd DAG1/parameter_sweep
  python sweep.py --stage 1          # sampler grid: lambda x tau_A ratio x t   (72 configs)
  python sweep.py --stage 2          # estimator/threshold variants around stage-1 top 3
"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import argparse, itertools, json, time, zlib
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "methodology"))   # DAG1/methodology
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
from scipy.stats import spearmanr, kendalltau
import config as cfg
from data import build_aligned, z_numeric, list_states
from sem import StateSystem
from probs import edge_probs, resolve_lambda, sink_probs
from ordering import sample_orderings
from dags import build_and_select

OUT = Path(__file__).resolve().parent / "results" / "details"   # stage CSVs + logs live here
OUT.mkdir(parents=True, exist_ok=True)
CM = pd.read_excel(cfg.PRIOR_V3, sheet_name="Matrix_80x80", index_col=0)
KNOBS = ["LAMBDA_MODE", "LAMBDA", "LAMBDA_Q", "TAU_A_RATIO", "T_SOFTMAX", "SHRINK_SCALE", "CONTEXT"]
DEFAULT = dict(LAMBDA_MODE="abs", LAMBDA=0.10, LAMBDA_Q=0.85, TAU_A_RATIO=0.25, T_SOFTMAX=0.1,
               SHRINK_SCALE=1.0, CONTEXT="all")


def cid(c):
    return "|".join(f"{k}={c[k]}" for k in KNOBS)


def apply(c):
    for k in KNOBS:
        setattr(cfg, k, c[k])
    cfg.TAU_A = c["TAU_A_RATIO"] * c["LAMBDA"]


def c_metrics(F, P, rng):
    """Agreement of consensus edge frequency F[i,j] (i->j) with C. Reported only."""
    C = CM.reindex(index=P, columns=P).values.astype(float); np.fill_diagonal(C, 0)
    iu = np.triu_indices(len(P), 1)
    Cu = np.where(np.abs(C[iu]) >= np.abs(C.T[iu]), C[iu], C.T[iu])
    U = (F + F.T)[iu]
    out = {}
    null = np.array([[U[rng.permutation(len(U))][Cu == c].mean() for c in (5, 3, -1, -2)] for _ in range(100)])
    for k, c in enumerate((5, 3, -1, -2)):
        v = U[Cu == c].mean()
        out[f"C_enrich_{c:+d}"] = v / U.mean()
        out[f"C_z_{c:+d}"] = (v - null[:, k].mean()) / (null[:, k].std() + 1e-12)
    asym = C > C.T
    a, d = (F * asym).sum(), (F.T * asym).sum()
    out["C_orient"] = a / (a + d) if a + d else np.nan
    return out


def run_half(keys, X, Znum, meta, c, seed, systems):
    apply(c)
    skey = (c["CONTEXT"], c["SHRINK_SCALE"])
    if skey not in systems:
        systems[skey] = StateSystem(X.loc[keys], Znum.loc[keys], meta.loc[keys, "dataset"].values)
    system = systems[skey]
    rng = np.random.default_rng(seed)
    P = system.proteins
    B0 = system.fit(list(range(len(P))))
    resolve_lambda(B0)
    A0 = edge_probs(B0)
    orderings, *_ = sample_orderings(system, cfg.N_ORDERINGS, rng)
    G, scores, best, ok = build_and_select(orderings, A0, rng)
    F = G.mean(0)
    pos = np.empty(orderings.shape, float)
    for m, o in enumerate(orderings):
        pos[m, o] = np.arange(len(P))
    p = sink_probs(A0)
    prs = rng.choice(len(orderings), size=(150, 2))
    kt = np.nanmean([kendalltau(pos[a], pos[b])[0] for a, b in prs if a != b])
    diag = dict(lambda_used=cfg.LAMBDA, density=float(G.sum((1, 2)).mean() / (len(P) * (len(P) - 1) / 2)),
                sink_entropy=float(-(p * np.log(p + 1e-300)).sum() / np.log(len(p))),
                sink_pmax=float(p.max()), kendall=float(kt),
                best_minus_mean=float(scores[best] - scores.mean()), acyclic=bool(ok))
    diag.update(c_metrics(F, P, rng))
    return F.astype(np.float32), pos.mean(0).astype(np.float32), diag


_G = {}


def _init():
    X, Z, meta = build_aligned()
    _G.update(X=X, Znum=z_numeric(Z), meta=meta)


def job(var, label, keys, half, configs):
    seed = (cfg.SEED + zlib.crc32(f"{var}|{label}|half{half}".encode())) % (2**32)
    systems, res = {}, []
    for c in configs:
        F, pos, diag = run_half(keys, _G["X"], _G["Znum"], _G["meta"], c, seed, systems)
        res.append((cid(c), F, pos, diag))
    return var, label, half, res


def stage_configs(stage):
    if stage == 1:
        grid = itertools.product([0.06, 0.08, 0.10, 0.12, 0.15, 0.20], [0.5, 0.25, 0.125], [0.02, 0.05, 0.1, 0.2])
        return [dict(DEFAULT, LAMBDA=l, TAU_A_RATIO=r, T_SOFTMAX=t) for l, r, t in grid]
    top = pd.read_csv(OUT / "stage1_ranking.csv").query("eligible").head(3)
    cs = []
    for _, r in top.iterrows():
        base = dict(DEFAULT, LAMBDA=r.LAMBDA, TAU_A_RATIO=r.TAU_A_RATIO, T_SOFTMAX=r.T_SOFTMAX)
        for ctx, sh in itertools.product(["all", "no_pmi", "none"], [1.0, 1.5, 2.0, 3.0]):
            cs.append(dict(base, CONTEXT=ctx, SHRINK_SCALE=sh))
        for q in [0.75, 0.80, 0.85, 0.90, 0.95]:
            cs.append(dict(base, LAMBDA_MODE="quantile", LAMBDA_Q=q))
    uniq = {cid(c): c for c in cs}
    return list(uniq.values())


def summarise(stage, store):
    rows = []
    iu = np.triu_indices(80, 1)
    for (var, label, key), halves in store.items():
        if len(halves) < 2:
            continue
        (F1, p1, d1), (F2, p2, d2) = halves[0], halves[1]
        U1, U2 = (F1 + F1.T)[iu], (F2 + F2.T)[iu]
        k = int(0.05 * len(U1))
        t1, t2 = set(np.argsort(-U1)[:k]), set(np.argsort(-U2)[:k])
        row = dict(variable=var, state=label, config=key,
                   rep_edges=spearmanr(U1, U2)[0], rep_pos=spearmanr(p1, p2)[0],
                   rep_top=len(t1 & t2) / len(t1 | t2))
        for kk in d1:
            if kk != "acyclic":
                row[kk] = np.nanmean([d1[kk], d2[kk]])
        row["acyclic"] = d1["acyclic"] and d2["acyclic"]
        rows.append(row)
    D = pd.DataFrame(rows)
    D.to_csv(OUT / f"stage{stage}_per_state.csv", index=False)
    agg = D.groupby("config").median(numeric_only=True).reset_index()
    agg["acyclic_all"] = D.groupby("config").acyclic.all().values
    for k in KNOBS:
        agg[k] = agg.config.str.extract(rf"{k}=([^|]+)")[0]
    for k in ["LAMBDA", "LAMBDA_Q", "TAU_A_RATIO", "T_SOFTMAX", "SHRINK_SCALE"]:
        agg[k] = agg[k].astype(float)
    agg["score"] = (agg.rep_edges + agg.rep_pos) / 2
    agg["eligible"] = agg.sink_entropy.between(0.2, 0.95)
    agg = agg.sort_values(["eligible", "score"], ascending=[False, False])
    agg.to_csv(OUT / f"stage{stage}_ranking.csv", index=False)
    return agg


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--min_n", type=int, default=150)
    ap.add_argument("--limit_states", type=int, default=None)   # for quick tests
    ap.add_argument("--limit_configs", type=int, default=None)
    a = ap.parse_args()
    X, Z, meta = build_aligned()
    states = [(v, L, k) for v, L, k in list_states(Z) if len(k) >= a.min_n][:a.limit_states]
    configs = stage_configs(a.stage)[:a.limit_configs]
    tasks = []
    for v, L, keys in states:
        rng = np.random.default_rng(zlib.crc32(f"split|{v}|{L}".encode()))
        perm = rng.permutation(len(keys)); h = len(keys) // 2
        for half, idx in enumerate([perm[:h], perm[h:2 * h]]):
            tasks.append((v, L, keys[np.sort(idx)], half, configs))
    print(f"stage {a.stage}: {len(configs)} configs x {len(states)} states x 2 halves = "
          f"{len(configs)*len(tasks)} pipeline runs, {a.jobs} workers", flush=True)
    t0, store, done = time.time(), {}, 0
    with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init) as ex:
        for f in as_completed([ex.submit(job, *t) for t in tasks]):
            var, label, half, res = f.result(); done += 1
            for key, F, pos, diag in res:
                store.setdefault((var, label, key), []).append((F, pos, diag))
            print(f"  [{done}/{len(tasks)}] {var}={label} half{half}  {time.time()-t0:.0f}s", flush=True)
    agg = summarise(a.stage, store)
    cols = ["config", "eligible", "score", "rep_edges", "rep_pos", "rep_top", "density", "sink_entropy",
            "kendall", "C_z_+5", "C_z_-1", "C_orient"]
    print(agg[cols].head(15).round(3).to_string(index=False))
    print(f"done in {time.time()-t0:.0f}s")
