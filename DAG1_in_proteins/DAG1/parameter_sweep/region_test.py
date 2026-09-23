"""Do the S-region-flagged states show more of the literature pattern than other states?
Run on one or more results folders (each must already have analysis_vs_C/ from compare_to_C.py).

  python region_test.py ../results ../runs/results_sweep_winner ...   (from DAG1/parameter_sweep)
For each folder and C class: enrichment in flagged vs other states (Mann-Whitney, one-sided in the
predicted direction), raw and after regressing out log(n); plus global +5 / -1 enrichment.
"""
import sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr

rows = []
for d in sys.argv[1:]:
    D = pd.read_csv(Path(d) / "analysis_vs_C" / "dag_vs_C_per_state.csv", dtype={"state": str})
    D = D[np.isfinite(D["all100_enrich_+5"])]
    ln = np.log(D.n)
    name = Path(d).resolve().name
    r = dict(results="run3 (DAG1/results)" if name == "results" else name,
             enrich5_all=D["all100_enrich_+5"].median(),
             z5_states_gt2=int((D["all100_z_+5"] > 2).sum()),
             enrich_m1_all=D["all100_enrich_-1"].median(),
             orient=D["all100_orient_agree"].median(),
             rho_enrich5_n=spearmanr(D["all100_enrich_+5"], D.n)[0])
    for c in ["+5", "+4", "+3", "+1", "-1", "-2"]:
        y = D[f"all100_enrich_{c}"]; f = D[f"flag_{c}"]
        res = y - np.polyval(np.polyfit(ln, y, 1), ln)
        alt = "greater" if c.startswith("+") else "less"
        r[f"p_{c}"] = mannwhitneyu(y[f], y[~f], alternative=alt).pvalue
        r[f"p_{c}_nadj"] = mannwhitneyu(res[f], res[~f], alternative=alt).pvalue
    b = D[D.variable == "braaksc"].set_index("state")["all100_enrich_+5"]
    r["braak_+5_by_stage"] = " ".join(f"{s}:{b[s]:.2f}" for s in b.index)
    rows.append(r)
T = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(T.round(3).to_string(index=False))
out = Path(__file__).resolve().parent / "results" / "region_test_summary.csv"
T.to_csv(out, index=False); print("saved", out)
