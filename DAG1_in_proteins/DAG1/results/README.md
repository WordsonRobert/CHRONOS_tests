# DAG1/results — the committed pipeline run (run 3)

Settings (the defaults in `../methodology/config.py`):
`A_ij = σ((|B_ji| − 0.10) / 0.025)`, softmax temperature t = 0.1, 100 orderings per state, Ledoit–Wolf
state SEM with all 14 trajectory variables + dataset as context, fixed seed per state.
Regenerate with `cd ../methodology && python run_all.py` (output to `DAG1/runs/`, identical files).

## Contents
- `summary.csv` — one row per state (57): n per cohort, reliability, diagnostics, best-DAG score, timing.
- `<variable>/<state>/` — per state:
  | file | contents |
  |---|---|
  | `B_init.csv` | the state's 80×80 SEM, `[target, source]` (keeps the sign of each effect) |
  | `A_init.csv` | edge probabilities, `A[i,j] = P(i→j)` |
  | `orderings.csv` | the 100 sampled orderings (pos1 = source-most, pos80 = final sink) |
  | `dags.npz` | all 100 DAGs (`G[m, source, target]`) and their scores — read with `../methodology/load.py` |
  | `best_dag.csv` / `best_dag_edges.csv` | the highest-scoring DAG (matrix / edge list with A and B) |
  | `meta.json` | same fields as the summary row |
- `analysis_vs_C/dag_vs_C_per_state.csv` — agreement with the literature matrix C per state
  (made by `../parameter_sweep/compare_to_C.py`).

The refit matrices (`refits.npz`, B at every removal step of every ordering, ~3 GB) are not included; a rerun
without `--no_refits` writes them.

## Summary of this run
| | value |
|---|---|
| states | 57 (13 variables; MMSE in boxes of 6, plaque/tangle in boxes of 3) |
| reliability | 10 comfortable (n ≥ 300), 11 workable, 5 exploratory, **31 with n < 80 (fitted to noise)** |
| DAG density | median 14% of possible edges (5–28%) |
| first-sink entropy / ordering Kendall τ | 0.66 / 0.25 (concentrated but stochastic) |
| all DAGs acyclic | yes |
| +5 edges vs a typical pair | 1.75× (above chance in 51/57 states) |
| −1 edges vs a typical pair | 0.92× |
| direction agreement with C | 0.49 (chance) |

## How to read it
- Prefer the **edge frequency across the 100 DAGs** (`dags.npz`) over the single best DAG: the best DAG shares
  only ~19% of its edges between reruns with different seeds, while the consensus is stable (0.86).
- Compare an edge's frequency with the **state's base rate** (mean frequency of all pairs in that state).
- **Arrow direction is not identified** — the state SEMs are nearly symmetric, and TMED-family proteins sit early
  in orderings because each has a strong partner, not because of biology.
- States with n < 80 (the `reliability` column) should not be interpreted.
