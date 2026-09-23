# DAG1 parameter sweep

Sweeps every tunable parameter of the pipeline in `../methodology` without changing its structure, and checks
whether any setting gives cleaner (more reproducible) DAGs and whether the conclusions depend on the settings.

**Selection criterion is C-free: reproducibility.** Each state with n ≥ 150 (21 states) is split once into two
random halves; the full pipeline runs on each half, and we measure how well the halves agree on
consensus edges (`rep_edges`), mean ordering positions (`rep_pos`) and top-5% edges (`rep_top`).
Rule fixed before results: maximise mean(rep_edges, rep_pos) among configs with median sink entropy in
[0.2, 0.95]. Agreement with the literature matrix C is reported for every config but never used to choose.

## Stages
| stage | what | configs |
|---|---|---|
| 1 | λ {0.06, 0.08, 0.10, 0.12, 0.15, 0.20} × τ_A/λ {0.125, 0.25, 0.5} × t {0.02, 0.05, 0.1, 0.2} | 72 |
| 2 | around the stage-1 top 3: context Z {all, no_pmi, none} × shrinkage ×{1, 1.5, 2, 3}, and per-state quantile λ {0.75…0.95} | 51 |
| 3 | full 57-state runs of the top configs, 2 extra seeds, 300 vs 100 orderings, C / region tests | — |

## Run (from this folder)
```
python sweep.py --stage 1 --jobs 8            # ~30 min on 2 cores -> results/details/stage1_*.csv
python sweep.py --stage 2 --jobs 8            # reads stage1_ranking.csv -> results/details/stage2_*.csv
bash results/details/stage3.sh                 # ~20 min; full runs go to DAG1/runs/
python compare_to_C.py [results folder]        # default DAG1/results -> <folder>/analysis_vs_C/
python region_test.py ../results ../runs/...   # -> results/region_test_summary.csv
python stability_check.py                      # needs the stage-3 runs -> results/details/stability_winner.csv
```
All scripts import the pipeline from `../methodology` and read C from `CHRONOS_80x80_confidence.xlsx`
at the top level; they run from any working directory. Reruns reproduce the committed tables exactly.

## Files in `results/`
| file | contents |
|---|---|
| `details/stage1_ranking.csv`, `details/stage2_ranking.csv` | one row per config: reproducibility, diagnostics, C agreement (medians over states), `eligible`, `score` |
| `details/stage1_per_state.csv`, `details/stage2_per_state.csv` | the same per config × state |
| `details/stability_winner.csv` | winner: seed-to-seed and 300-vs-100-ordering agreement per state |
| `details/stage*.log`, `details/stage3.sh` | run logs and the stage-3 script |
| `region_test_summary.csv` | C agreement and region test for run 3 and the 6 stage-3 configs |

## Results
**Winner:** λ = 0.12, τ_A = 0.5·λ, t = 0.05, shrinkage ×1, context = none (dataset only).
Rerun: `cd ../methodology && python run_all.py --set LAMBDA=0.12 TAU_A_RATIO=0.5 T_SOFTMAX=0.05 CONTEXT=none --out results_sweep_winner`.
It mainly improves ordering reproducibility (0.65 vs ~0.50 for run 3); edge reproducibility does not move.

| agreement (median, states with n ≥ 150) | edges | ordering |
|---|---|---|
| same subjects, different seeds (sampler noise) | 0.86 | 0.99 |
| 300 vs 100 orderings | 0.90 | 1.00 |
| **different halves of the subjects (data noise)** | **0.125** | 0.65 |
| best single DAG across seeds (Jaccard) | 0.19 | |

- Edge reproducibility between halves stayed at 0.10–0.15 in **all 123 configs**. The raw state SEM |B| reproduces
  0.12 between halves vs 0.11 between different states (plain correlations 0.57 vs 0.56): the limit is the
  80-protein partial regression at n = 150–900, not the sampler.
- The committed run 3 (`DAG1/results`) was chosen before the sweep; the sweep shows the conclusions below hold
  across the whole grid, so run 3 is representative.

| run | +5 enrichment | states with +5 above chance | −1 enrichment | direction agreement with C | flagged-vs-other p (+5 / +4 / −1) |
|---|---|---|---|---|---|
| run 3 (`DAG1/results`) | 1.75 | 51/57 | 0.92 | 0.49 | 0.37 / 0.05 / 0.31 |
| winner | 1.35 | 50 | 0.97 | 0.47 | 0.40 / 0.08 / 0.40 |
| top 2 | 1.32 | 50 | 0.98 | 0.48 | 0.48 / 0.11 / 0.24 |
| top 3 | 1.33 | 42 | 0.98 | 0.49 | 0.63 / 0.18 / 0.35 |
| top 4 | 1.30 | 50 | 0.97 | 0.48 | 0.58 / 0.03 / 0.31 |
| top 5 | 1.36 | 48 | 0.97 | 0.49 | 0.38 / 0.07 / 0.35 |
| best with context | 1.30 | 47 | 0.98 | 0.50 | 0.48 / 0.03 / 0.44 |

"Flagged" = states picked out by the trajectory regions; p is one-sided Mann-Whitney vs all other states.
Only +4 dips below 0.05 (2 of 7 runs), which is what chance gives across 42 tests.
