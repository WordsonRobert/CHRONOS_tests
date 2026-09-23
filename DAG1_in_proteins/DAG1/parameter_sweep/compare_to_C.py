"""After-the-fact comparison of DAG1 DAGs with the literature matrix C (never used to build them).

Per state, over all 100 DAGs (and the best DAG):
  presence_c   = fraction of C-class-c directed cells (C row=source, col=target) that are DAG edges
  undirected_c = fraction of class-c pairs with an edge in either direction
  enrich_c     = undirected_c / undirected rate of ALL pairs   (>1 = over-represented)
  z_c          = enrichment vs 200 label-shuffles of C over the same pairs
  orient_agree = among pairs where C_ij > C_ji and the DAG has an edge between them,
                 fraction oriented i->j (literature direction); 0.5 = chance
"""
import numpy as np, pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "methodology"))   # DAG1/methodology
import config as cfg

# usage: python compare_to_C.py [results folder]   (default: DAG1/results = run 3)
RES = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else cfg.DAG1 / "results"
Cm = pd.read_excel(cfg.PRIOR_V3, sheet_name="Matrix_80x80", index_col=0)
CLASSES = [5, 4, 3, 2, 1, -1, -2]
rng = np.random.default_rng(0)


def state_stats(d):
    z = np.load(d / "dags.npz"); G = z["G"]; P = list(z["proteins"]); best = int(np.argmax(z["scores"]))
    C = Cm.reindex(index=P, columns=P).values.astype(float); np.fill_diagonal(C, 0)
    iu = np.triu_indices(len(P), 1)
    # undirected class per pair = the stronger-|C| direction (keeps sign)
    Cu = np.where(np.abs(C[iu]) >= np.abs(C.T[iu]), C[iu], C.T[iu])
    out = {}
    for tag, Gs in (("all100", G), ("best", G[best:best + 1])):
        U = (Gs | Gs.transpose(0, 2, 1))[:, iu[0], iu[1]].mean(0)      # P(edge between pair)
        base = U.mean()
        null = np.array([[U[rng.permutation(len(U))][Cu == c].mean() for c in CLASSES] for _ in range(200)])
        for k, c in enumerate(CLASSES):
            m = Cu == c
            val = U[m].mean()
            out[f"{tag}_enrich_{c:+d}"] = val / base
            out[f"{tag}_z_{c:+d}"] = (val - null[:, k].mean()) / null[:, k].std()
        # orientation agreement with literature asymmetry
        asym = C > C.T
        E = Gs.sum(0) if tag == "all100" else Gs[0].astype(int)
        agree = (E * asym).sum(); disagree = (E.T * asym).sum()
        out[f"{tag}_orient_agree"] = agree / (agree + disagree) if agree + disagree else np.nan
    return out


# --- states flagged by the PMI-free S-regions (mapped onto DAG1 state labels / boxes) ---
FLAG = {
 "+5": [("TangleTotal","9-12"),("cts_mmse30_lv","0-6"),("cts_mmse30_lv","6-12"),("amyA","3"),("amyThal","3"),
        ("braaksc","3"),("braaksc","4"),("braaksc","5"),("ceradsc","1"),("cts_mmse30_first_ad_dx","12-18"),("cts_mmse30_lv","24-30")],
 "+4": [("amyCerad","1"),("amyCerad","2"),("reag","2"),("braaksc","4"),("braaksc","5"),("cts_mmse30_lv","18-24")],
 "+3": [("amyA","2"),("amyA","3"),("amyCerad","1"),("amyThal","1"),("amyThal","2"),("braaksc","4"),("braaksc","5"),("ceradsc","1")],
 "+1": [("amyA","2"),("amyCerad","1"),("amyCerad","3"),("amyThal","1"),("amyThal","2"),("braaksc","1"),("braaksc","2"),("reag","2"),("ceradsc","1")],
 "-1": [("amyCerad","2"),("braaksc","3"),("braaksc","4"),("braaksc","5"),("ceradsc","1")],
 "-2": [("amyCerad","2"),("braaksc","3"),("braaksc","4"),("braaksc","5"),("cts_mmse30_lv","24-30")],
}

if __name__ == "__main__":
    S = pd.read_csv(RES / "summary.csv", dtype={"state": str})
    rows = []
    for _, r in S.iterrows():
        st = state_stats(RES / r.variable / str(r.state))
        rows.append(dict(variable=r.variable, state=r.state, n=r.n, reliability=r.reliability, **st,
                         **{f"flag_{k}": (r.variable, str(r.state)) in v for k, v in FLAG.items()}))
    D = pd.DataFrame(rows)
    out = RES / "analysis_vs_C"; out.mkdir(exist_ok=True)
    D.to_csv(out / "dag_vs_C_per_state.csv", index=False)
    print("saved", out / "dag_vs_C_per_state.csv")
