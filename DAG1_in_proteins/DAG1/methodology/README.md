# DAG1 methodology — state-conditioned DAGs (one DAG set per trajectory state)

For every value (or box) of 13 trajectory variables (all 14 except pmi):

1. take only the subjects in that state (from the 14 CSVs in `DATA/`, one row per subject)
2. fit the 80-protein SEM  `X = B X + Γ Z + ε`  → `B` (80×80, `[target, source]`)
3. edge probabilities `A_ij = P(i→j) = σ((|B_ji| − λ) / τ_A)`   (λ = 0.10, τ_A = 0.025)
4. `p_sink[i] = 1 − max_j A_ij`,  `p = softmax(p_sink / t)` (t = 0.1), sample a sink
5. remove it and **refit** the SEM on the remaining proteins (same subjects), repeat to 1 → one ordering
6. 100 orderings → 100 DAGs (`G_ij ~ Bernoulli(A_ij)` only if i precedes j; always acyclic)
7. best DAG = argmax `S(G;A) = mean_{(i,j)∈E_G} A_ij`, with A from step 3 of that state

The literature matrix C is not used anywhere in building the DAGs.

## Run (from this folder, `DAG1/methodology/`)
```
python run_all.py                     # every state, all CPU cores  -> DAG1/runs/results_thr_l0.1_tA0.025_t0.1/
python run_all.py --jobs 8 --vars amyA braaksc
python run_all.py --no_refits         # skip the ~3 GB of refit matrices
python run_all.py --set LAMBDA=0.12 T_SOFTMAX=0.05 --out my_run      # any config override -> DAG1/runs/my_run/
python run_state.py --var amyA --state 1
```
Inputs are found automatically: `DATA/` (the 14 CSVs) two levels up. New runs always go to `DAG1/runs/`,
never into `DAG1/results/` (the committed run). With the default settings a rerun reproduces
`DAG1/results/` exactly (orderings, all 100 DAGs, B) — seeds are fixed per state.

Everything is set in `config.py`: variables and box widths (MMSE boxes of 6, Plaque/Tangle boxes of 3),
λ, τ_A, t, number of orderings, seed, solver, shrinkage, context Z, and whether to save refits.

## Modules
| file | job |
|---|---|
| `config.py` | all settings and paths |
| `data.py` | builds the aligned subject table from the 14 CSVs; defines states/boxes |
| `sem.py` | state SEM fit (`lw` fast closed form; `enet` = original elastic-net solver) |
| `probs.py` | B→A, sink scores, softmax |
| `ordering.py` | recursive sink sampling with refit at every step |
| `dags.py` | Bernoulli DAGs, score, best-DAG selection, acyclicity check |
| `diagnostics.py` | sink entropy, ordering Kendall τ, DAG density (reported in `summary.csv`) |
| `run_state.py` / `run_all.py` | one state / all states in parallel |
| `load.py` | read back any saved refit matrix or DAG |

## Output — `<results folder>/<variable>/<state>/`
| file | contents |
|---|---|
| `B_init.csv` | the state's 80×80 SEM, `[target, source]` |
| `A_init.csv` | edge probabilities, `A[i,j] = P(i→j)` |
| `orderings.csv` | 100 orderings (pos1 = source-most, pos80 = final sink) |
| `dags.npz` | all 100 DAGs (`G[m, source, target]`) + their scores |
| `best_dag.csv` / `best_dag_edges.csv` | the selected DAG (matrix / edge list with A and B) |
| `refits.npz` | every refit matrix B^(80…1) of every ordering (read with `load.refit`; not in the repo, ~3 GB) |
| `meta.json` | n per cohort, reliability, diagnostics, scores, timing |

`summary.csv` has one row per state.

## Solver note
The recursion needs ~3,240 node regressions per ordering × 100 orderings × 57 states ≈ 18.5 M fits.
With the cross-validated elastic net that's ~77 h on 2 cores. The default `lw` solver fits the same model:
proteins are residualised on the context Z (exact by Frisch–Waugh–Lovell), then node-wise regressions come
from the inverse of the Ledoit–Wolf shrunk covariance. Each refit re-estimates the shrinkage and re-inverts
the smaller covariance, so it is a real refit. Checked on the full population: identical to sklearn's
LedoitWolf (1e-16), and corr 0.887 with the original elastic-net 80×80 (`DAG1/general_graph`) with 100% sign
agreement on every edge with |B| > 0.02. The whole run takes ~110 s on 2 cores.

## How the B→A map was settled (three runs)
| run | B→A map | DAG density | sink entropy | ordering Kendall τ | verdict |
|---|---|---|---|---|---|
| 1 | σ(B/0.1) | 52% | 0.97 | 0.11 | B≈0 → A=0.5 coin flips; near-random |
| 2 | σ((\|B\|−0.02)/0.01) | 68% | 1.00 | 0.00 | λ below the median \|B\| (0.038) of the state SEMs → every A≈1, sinks uniform |
| **3** | σ((\|B\|−0.10)/0.025) | **14%** | **0.66** | **0.25** | concentrated but stochastic — **this is `DAG1/results/`** |

Only run 3 is committed. Runs 1–2 can be regenerated with
`python run_all.py --no_refits --set A_MAP=signed --out run1` and
`python run_all.py --no_refits --set LAMBDA=0.02 TAU_A=0.01 --out run2`.

λ=0.10 was chosen from the state-SEM scale before looking at any biological result: median |B| = 0.038,
each protein's strongest outgoing |B| = 0.11–0.40, and a λ scan on the saved B matrices showed 0.08–0.12 is
where DAG density lands in 10–30% and the sink score spreads. The sign of B is kept separately in
`B_init.csv` / `best_dag_edges.csv`.

## Parameter sweep
123 configurations, selected by a C-free reproducibility criterion — see `../parameter_sweep/README.md`.
Conclusions do not change across the grid, and the committed run-3 settings are representative.

## Sample sizes
`summary.csv → reliability`: n < 80 = p>n (overfit), < 150 exploratory, < 300 workable, ≥ 300 comfortable.
Many boxes are small (braak 0: 6, MMSE-at-dx 0–6: 6, apoe 22: 8), and those DAGs are fitted to noise.
