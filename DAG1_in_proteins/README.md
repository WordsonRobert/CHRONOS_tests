# DAG1 in proteins — state-conditioned causal DAGs over the 80 CHRONOS proteins

Proteomics from three AMP-AD cohorts (ROSMAP, Diverse Cohorts, Banner; 1,570 subjects) are used to build
causal DAGs over the 80 CHRONOS proteins, separately for each value of 13 trajectory variables
(Braak stage, amyloid staging, CERAD, MMSE, plaque/tangle load, APOE, sex, …). The literature
confidence matrix C is never used to build any graph; it is only compared against afterwards.

## Layout
| path | what it is | README |
|---|---|---|
| `DATA/` | the 14 per-variable CSVs (one row per protein × subject) — the input to everything | `DATA/README.md` |
| `CHRONOS_80x80_confidence.xlsx` | literature matrix C (v3: scores −2…+5, `Matrix_80x80`, `Sign_80x80`, `Edge_List`, `Rubric`) — used only in comparisons | — |
| `ground_truth/` | where along each trajectory the literature-scored protein pairs co-move / stay flat (exploratory) | its `README.md` |
| `DAG1/general_graph/` | one pooled 80×80 causal graph B (no conditioning) | its `README.md` |
| `DAG1/methodology/` | the state-conditioned DAG pipeline (code) | its `README.md` |
| `DAG1/results/` | the committed run of the pipeline (run 3): 57 states × 100 DAGs + best DAG | its `README.md` |
| `DAG1/parameter_sweep/` | 123-configuration sweep, stability checks, C / region comparisons | its `README.md` |
| `DAG1/runs/` | created when you rerun anything; new outputs go here, never over committed results | — |

## Quick start
```
pip install numpy pandas scipy scikit-learn openpyxl
cd DAG1/methodology && python run_all.py --no_refits      # reruns all 57 states (~2 min on 2 cores)
```
Every script uses paths relative to its own location and reproduces the committed outputs exactly
(fixed seeds). Tested with Python 3.11.

## Main results (details in the folder READMEs)
- **+5 literature edges** (structurally resolved complexes) are over-represented among DAG edges in essentially
  every state (median 1.75×, above chance in 51/57 states); −1 edges are mildly under-represented.
- **Edge direction** agrees with C's direction only at chance level (~0.49) — the pipeline does not identify direction.
- **States flagged by the trajectory regions** show no more literature agreement than other states.
- **Reproducibility:** fitting a state on two halves of its subjects gives edge agreement ≈ 0.12, about the same
  as between two different states; the sampler itself is stable (seed agreement 0.86). Use the consensus edge
  frequency over the 100 DAGs per state, not the single best DAG. States with n < 80 are fitted to noise.

## Data
The CSVs in `DATA/` and `DAG1/general_graph/aligned_XZ.csv` are individual-level data derived from AMP-AD
(Synapse) and are shared here only with project members under the AMP-AD data-use terms.
