# DAG1

Causal graphs over the 80 CHRONOS proteins, built from the protein data only (the literature matrix C is
used only to compare against afterwards).

| folder | contents |
|---|---|
| `general_graph/` | **one pooled 80×80 graph** B from `X = BX + ΓZ + ε` over all 1,570 subjects (elastic net, node-wise) |
| `methodology/` | **the state-conditioned pipeline**: for each trajectory state, fit the state SEM, sample 100 orderings by recursive sink removal with a refit at every step, draw 100 DAGs, pick the best |
| `results/` | the committed pipeline run (run 3): 57 states, all DAGs, best DAG, diagnostics, C comparison |
| `parameter_sweep/` | 123-config sweep (C-free reproducibility criterion), stability checks, C / region robustness across the top configs |
| `runs/` | created by reruns; everything new goes here |

Order to read: `methodology/README.md` → `results/README.md` → `parameter_sweep/README.md`.
`general_graph/` is the unconditioned reference graph that the state SEMs were checked against.
