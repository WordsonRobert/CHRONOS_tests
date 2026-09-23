#!/bin/bash
# Stage 3 of the parameter sweep: full runs of the top configs + stability checks + C / S-region robustness.
# Run from anywhere:  bash DAG1/parameter_sweep/results/details/stage3.sh      (~20 min on 2 cores)
# Full runs are written to DAG1/runs/<name>/ ; run 3 itself is the committed DAG1/results.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"          # DAG1/parameter_sweep/results/details
SWEEP="$(cd "$HERE/../.." && pwd)"             # DAG1/parameter_sweep
METH="$(cd "$SWEEP/../methodology" && pwd)"   # DAG1/methodology
RUNS="$SWEEP/../runs"                          # DAG1/runs
RUN3="$SWEEP/../results"                       # DAG1/results (run 3)

cd "$METH"
W="--set LAMBDA=0.12 TAU_A_RATIO=0.5 T_SOFTMAX=0.05 SHRINK_SCALE=1.0 CONTEXT=none"
python3 run_all.py --jobs 2 --no_refits $W --out results_sweep_winner      # winner, all 57 states
python3 run_all.py --jobs 2 --no_refits --set LAMBDA=0.15 TAU_A_RATIO=0.5 T_SOFTMAX=0.1 CONTEXT=none  --out results_sweep_top2
python3 run_all.py --jobs 2 --no_refits --set LAMBDA=0.12 TAU_A_RATIO=0.5 T_SOFTMAX=0.05 SHRINK_SCALE=1.5 CONTEXT=none --out results_sweep_top3
python3 run_all.py --jobs 2 --no_refits --set LAMBDA_MODE=quantile LAMBDA_Q=0.75 TAU_A_RATIO=0.5 T_SOFTMAX=0.05 --out results_sweep_top4
python3 run_all.py --jobs 2 --no_refits --set LAMBDA_MODE=quantile LAMBDA_Q=0.9 TAU_A_RATIO=0.5 T_SOFTMAX=0.1 --out results_sweep_top5
python3 run_all.py --jobs 2 --no_refits --set LAMBDA=0.15 TAU_A_RATIO=0.5 T_SOFTMAX=0.1 CONTEXT=all --out results_sweep_best_with_context
# stability of the winner: 2 more seeds, and 300 orderings (reliable states only)
python3 run_all.py --jobs 2 --no_refits $W --min_n 150 --seed 111 --out results_sweep_winner_seed111
python3 run_all.py --jobs 2 --no_refits $W --min_n 150 --seed 222 --out results_sweep_winner_seed222
python3 run_all.py --jobs 2 --no_refits $W --min_n 150 --n_orderings 300 --out results_sweep_winner_N300

cd "$SWEEP"
DIRS="$RUN3 $RUNS/results_sweep_winner $RUNS/results_sweep_top2 $RUNS/results_sweep_top3 $RUNS/results_sweep_top4 $RUNS/results_sweep_top5 $RUNS/results_sweep_best_with_context"
for d in $DIRS; do python3 compare_to_C.py "$d"; done
python3 region_test.py $DIRS            # -> DAG1/parameter_sweep/results/region_test_summary.csv
python3 stability_check.py              # -> DAG1/parameter_sweep/results/details/stability_winner.csv
echo STAGE3_DONE
