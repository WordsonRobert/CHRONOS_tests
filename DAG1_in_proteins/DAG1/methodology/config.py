"""DAG1 configuration — every knob lives here.

Paths are relative to this file:  <repo folder>/DAG1/methodology/config.py
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent            # DAG1/methodology
DAG1 = HERE.parent                                # DAG1
ROOT = DAG1.parent                                # top-level folder (DATA/, confidence matrices)
CSV_DIR = ROOT / "DATA"                           # the 14 per-variable CSVs
PRIOR_V3 = ROOT / "CHRONOS_80x80_confidence.xlsx"      # literature matrix C, v3 contents (comparisons only, never used to build DAGs)
RUNS = DAG1 / "runs"                              # new runs are written here (DAG1/results = committed run 3)

# ---- trajectory variables to condition on (13 = all 14 minus pmi) ----
# None  -> each distinct value is its own state
# (w,M) -> continuous: boxes of width w on the scale [0, M] (top value goes in the last box)
STATE_VARS = {
    "sex": None,
    "apoe_genotype": None,
    "cts_mmse30_lv": (6, 30),            # MMSE boxes of 6
    "cts_mmse30_first_ad_dx": (6, 30),   # MMSE boxes of 6
    "braaksc": None,
    "ceradsc": None,
    "amyThal": None,
    "amyA": None,
    "amyCerad": None,
    "amyAny": None,
    "reag": None,
    "PlaqueTotal": (3, 15),              # boxes of 3 (same 0-15 scale as TangleTotal)
    "TangleTotal": (3, 15),              # boxes of 3
}
# all 14 variables are available as context Z inside each state's SEM
# (the conditioning variable is constant inside a discrete state and drops out automatically)
Z_VARS = ["sex", "apoe_genotype", "cts_mmse30_lv", "cts_mmse30_first_ad_dx", "pmi",
          "braaksc", "ceradsc", "amyThal", "amyA", "amyCerad", "amyAny", "reag",
          "PlaqueTotal", "TangleTotal"]

# ---- algorithm parameters ----
SOLVER = "lw"          # "lw"   = closed-form node-wise SEM (Ledoit-Wolf shrinkage), fast, default
                       # "enet" = cross-validated elastic net per node (the original 80x80 solver; very slow)
A_MAP = "threshold"    # how B becomes edge probability A_ij = P(i->j)   (B is [target, source])
                       #  "threshold": A_ij = sigmoid((|B_ji| - LAMBDA) / TAU_A)   -> B~0 gives A~0
                       #  "signed":    A_ij = sigmoid(B_ji / TAU_SIGMOID)          (run 1; B=0 gives 0.5)
LAMBDA = 0.10          # effect-size threshold, calibrated to the state-SEM |B| scale (median |B| ~0.04)
TAU_A = 0.025
TAU_SIGMOID = 0.1
T_SOFTMAX = 0.1        # p   = softmax(p_sink / T_SOFTMAX)       ("small temperature")
LAMBDA_MODE = "abs"    # "abs": LAMBDA as given | "quantile": LAMBDA = LAMBDA_Q-quantile of |B_init| per state
LAMBDA_Q = 0.85        #   (used only in quantile mode; TAU_A then = TAU_A_RATIO * LAMBDA)
TAU_A_RATIO = 0.25
SHRINK_SCALE = 1.0     # multiplies the Ledoit-Wolf shrinkage intensity (1 = pure Ledoit-Wolf)
CONTEXT = "all"        # context Z adjusted for: "all" (14 Z + dataset) | "no_pmi" | "none" (dataset only)
def run_name():
    return f"thr_l{LAMBDA}_tA{TAU_A}_t{T_SOFTMAX}" if A_MAP == "threshold" else "signed"


RUN_NAME = run_name()
RESULTS = RUNS / f"results_{RUN_NAME}"
N_ORDERINGS = 100      # one Bernoulli DAG is drawn per ordering -> 100 DAGs per state
MIN_N = 5              # states with fewer subjects are skipped
MIN_RESID_DF = 10      # if n - (#context columns) - 1 < this, context Z is dropped for that state
SAVE_REFITS = True     # save every refit matrix B^(80..1) of every ordering (big: ~15-35 MB/state)
SEED = 20260923
