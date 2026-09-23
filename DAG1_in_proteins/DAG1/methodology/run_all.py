"""Run every state of every trajectory variable, in parallel across states.

  python run_all.py                 # all cores; default settings reproduce DAG1/results (run 3)
  python run_all.py --jobs 4
  python run_all.py --vars amyA braaksc
  python run_all.py --set LAMBDA=0.12 T_SOFTMAX=0.05 --out my_run   # any config override

Output goes to DAG1/runs/<name>/ (default name from the settings, e.g. results_thr_l0.1_tA0.025_t0.1).
Add --no_refits to skip the ~3 GB of refit matrices.
"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")          # one BLAS thread per worker -> no oversubscription
import argparse, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import config as cfg
from data import build_aligned, z_numeric, list_states
from run_state import run_state

_G = {}


def _apply_overrides(a):
    for kv in a.set:
        k, v = kv.split("=", 1)
        old = getattr(cfg, k)
        setattr(cfg, k, type(old)(v) if not isinstance(old, str) else v)
    if cfg.LAMBDA_MODE == "abs" and not any(kv.startswith("TAU_A=") for kv in a.set):
        cfg.TAU_A = cfg.TAU_A_RATIO * cfg.LAMBDA
    if a.seed is not None: cfg.SEED = a.seed
    if a.n_orderings is not None: cfg.N_ORDERINGS = a.n_orderings
    if a.no_refits: cfg.SAVE_REFITS = False
    cfg.RESULTS = cfg.RUNS / (a.out or f"results_{cfg.run_name()}")   # DAG1/runs/<name>/
    cfg.RESULTS.mkdir(parents=True, exist_ok=True)


def _init(overrides=None):
    if overrides is not None:
        _apply_overrides(overrides)
    X, Z, meta = build_aligned()
    _G.update(X=X, Znum=z_numeric(Z), meta=meta)


def _job(args):
    var, label, keys, solver = args
    return run_state(var, label, keys, _G["X"], _G["Znum"], _G["meta"], solver=solver,
                     n_orderings=cfg.N_ORDERINGS, outroot=cfg.RESULTS)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--vars", nargs="*", default=None)
    ap.add_argument("--solver", default=cfg.SOLVER)
    ap.add_argument("--set", nargs="*", default=[], help="override config, e.g. --set LAMBDA=0.12 T_SOFTMAX=0.05")
    ap.add_argument("--out", default=None, help="results folder name (default from RUN_NAME)")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--n_orderings", type=int, default=None)
    ap.add_argument("--min_n", type=int, default=0, help="only states with n >= this")
    ap.add_argument("--no_refits", action="store_true")
    a = ap.parse_args()
    _apply_overrides(a)
    X, Z, meta = build_aligned()
    states = [(v, L, k, a.solver) for v, L, k in list_states(Z)
              if (not a.vars or v in a.vars) and len(k) >= a.min_n]
    states.sort(key=lambda s: -len(s[2]))            # big states first -> better load balance
    print(f"{len(states)} states, {a.jobs} workers, solver={a.solver}", flush=True)
    t0, rows = time.time(), []
    with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init, initargs=(a,)) as ex:
        futs = [ex.submit(_job, s) for s in states]
        for f in as_completed(futs):
            r = f.result(); rows.append(r)
            print(f"  [{len(rows)}/{len(states)}] {r['variable']}={r['state']} n={r['n']} "
                  f"{r.get('status')} best={r.get('best_score')} {r.get('seconds','')}s", flush=True)
    order = {v: i for i, v in enumerate(cfg.STATE_VARS)}
    summ = pd.DataFrame(rows)
    summ["_o"] = summ.variable.map(order)
    summ.sort_values(["_o", "state"]).drop(columns="_o").to_csv(cfg.RESULTS / "summary.csv", index=False)
    print(f"done in {time.time()-t0:.0f}s -> {cfg.RESULTS/'summary.csv'}")
