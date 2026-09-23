"""Run the full pipeline for ONE trajectory state (e.g. amyA = 1).

  python run_state.py --var amyA --state 1
"""
import argparse, json, time, zlib
import numpy as np
import pandas as pd
import config as cfg
from data import build_aligned, z_numeric, list_states
from sem import StateSystem
from probs import edge_probs, resolve_lambda
from ordering import sample_orderings
from dags import build_and_select
from diagnostics import first_sink_entropy, ordering_kendall, dag_density


def reliability(n):
    return ("p>n (overfit)" if n < 80 else "exploratory" if n < 150
            else "workable" if n < 300 else "comfortable")


def run_state(var, label, keys, X, Znum, meta, solver=cfg.SOLVER, n_orderings=cfg.N_ORDERINGS,
              outroot=cfg.RESULTS):
    t0 = time.time()
    out = outroot / var / str(label)
    out.mkdir(parents=True, exist_ok=True)
    n = len(keys)
    ds = meta.loc[keys, "dataset"]
    row = dict(variable=var, state=label, n=n, reliability=reliability(n),
               **{f"n_{d}": int((ds == d).sum()) for d in ["ROSMAP", "Diverse", "Banner"]})
    if n < cfg.MIN_N:
        row.update(status="skipped (n<MIN_N)")
        json.dump(row, open(out / "meta.json", "w"), indent=2)
        return row
    seed = (cfg.SEED + zlib.crc32(f"{var}|{label}".encode())) % (2**32)
    rng = np.random.default_rng(seed)
    P = list(X.columns)

    system = StateSystem(X.loc[keys], Znum.loc[keys], ds.values, solver=solver)
    B0 = system.fit(list(range(len(P))))
    lam = resolve_lambda(B0)
    A0 = edge_probs(B0)
    pd.DataFrame(B0, index=P, columns=P).to_csv(out / "B_init.csv")            # [target, source]
    pd.DataFrame(A0, index=P, columns=P).to_csv(out / "A_init.csv")            # A_ij = P(i->j)

    orderings, step_refit, masks, refitB = sample_orderings(system, n_orderings, rng)
    pd.DataFrame([[P[i] for i in o] for o in orderings],
                 columns=[f"pos{k+1}" for k in range(len(P))]).to_csv(out / "orderings.csv", index_label="ordering")

    G, scores, best, acyclic_ok = build_and_select(orderings, A0, rng)
    np.savez_compressed(out / "dags.npz", G=G, scores=scores, proteins=np.array(P))
    pd.DataFrame(G[best].astype(int), index=P, columns=P).to_csv(out / "best_dag.csv")  # [source, target]
    bi, bj = np.nonzero(G[best])
    pd.DataFrame({"source": [P[i] for i in bi], "target": [P[j] for j in bj],
                  "A": A0[bi, bj], "B_coef": B0[bj, bi]}).sort_values("A", ascending=False) \
      .to_csv(out / "best_dag_edges.csv", index=False)

    if cfg.SAVE_REFITS:
        sizes = np.array([m.sum() for m in masks])
        flat = np.concatenate([b.ravel() for b in refitB]).astype(np.float32)
        np.savez_compressed(out / "refits.npz", flat=flat, sizes=sizes,
                            masks=np.array(masks), step_refit=step_refit,
                            orderings=orderings, proteins=np.array(P))

    sink_H, sink_pmax = first_sink_entropy(A0)
    row.update(a_map=cfg.A_MAP, lambda_used=round(lam, 4), mean_A=round(float(A0[~np.eye(len(P), dtype=bool)].mean()), 4),
               dag_density=round(dag_density(G), 4), sink_entropy=round(sink_H, 4),
               sink_pmax=round(sink_pmax, 4), ordering_kendall=round(ordering_kendall(orderings), 4),
               best_minus_mean=round(float(scores[best] - scores.mean()), 4))
    row.update(status="ok", solver=solver, context_cols=len(system.context_cols),
               z_context_dropped=system.z_dropped, acyclic_ok=acyclic_ok,
               best_dag=int(best), best_score=round(float(scores[best]), 4),
               score_mean=round(float(scores.mean()), 4), score_sd=round(float(scores.std()), 4),
               best_edges=int(G[best].sum()), mean_edges=round(float(G.sum((1, 2)).mean()), 1),
               unique_refits=len(refitB), seed=int(seed), seconds=round(time.time() - t0, 1))
    json.dump(row, open(out / "meta.json", "w"), indent=2)
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--var", required=True)
    ap.add_argument("--state", required=True)
    ap.add_argument("--solver", default=cfg.SOLVER)
    A = ap.parse_args()
    X, Z, meta = build_aligned()
    Znum = z_numeric(Z)
    for v, L, keys in list_states(Z):
        if v == A.var and str(L) == A.state:
            print(run_state(v, L, keys, X, Znum, meta, solver=A.solver))
            break
    else:
        print("state not found; available:",
              [(v, L) for v, L, _ in list_states(Z) if v == A.var])
