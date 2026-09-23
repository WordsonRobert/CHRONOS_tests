"""Build expression_explorer.html: 13 interactive plots (one per clinical variable), all 80 proteins in each.

Per variable and bucket (same buckets as the DAG states):
  levels   = mean protein level (within-cohort z-score) of everyone in the bucket, +/- SE, n
  held     = change vs the first bucket with every other clinical variable, sex, APOE, PMI and cell types held
             constant (cumulative sum of the pooled adjacent-bucket steps in bucket_contrasts.csv; APOE: vs 33)
  effect   = per-step effect with everything else held (effects_per_step.csv, cell version): beta, p, q
Buckets with fewer than MIN_N people are left out of the plots.
  python build_explorer.py        (run one_at_a_time.py first)
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
CW = HERE.parent
for d in ["DAG2-with_cell_composition", "DAG2"]:
    if (CW / d).exists(): sys.path.insert(0, str(CW / d / "methodology")); break
import config as cfg
cfg.COMPOSITION = "cell"
from data import build_aligned

MIN_N = 10
VARS = [("sex", "Sex", "female → male", "Sex"),
        ("apoe_genotype", "APOE genotype", "vs 33 when held", "Risk gene"),
        ("cts_mmse30_lv", "MMSE, last visit", "0–30 · higher = better", "Cognition"),
        ("cts_mmse30_first_ad_dx", "MMSE at AD diagnosis", "0–30 · higher = better · ROSMAP only", "Cognition"),
        ("braaksc", "Braak stage", "0–6 · tangle spread", "Tau"),
        ("TangleTotal", "Tangle total", "0–15 · Banner only", "Tau"),
        ("ceradsc", "CERAD score", "0–3 · ROSMAP + Banner", "Amyloid"),
        ("amyCerad", "CERAD score (Diverse)", "0–3 · Diverse only", "Amyloid"),
        ("amyThal", "Thal phase", "1–5 · amyloid spread · Diverse", "Amyloid"),
        ("amyA", "A-score", "1–3 · binned Thal · Diverse", "Amyloid"),
        ("amyAny", "Any amyloid", "no → yes · Diverse", "Amyloid"),
        ("PlaqueTotal", "Plaque total", "0–15 · Banner only", "Amyloid"),
        ("reag", "NIA-Reagan", "0–3 · Braak + CERAD combined · Diverse", "Combined")]

X, Z, meta = build_aligned()
P = list(X.columns)
eff = pd.read_csv(HERE / "effects_per_step.csv")
con = pd.read_csv(HERE / "bucket_contrasts.csv", dtype={"from_bucket": str, "to_bucket": str})
eff_c = eff[eff.version == "cell"]; con_c = con[con.version == "cell"]


def buckets(v):
    s = Z[v]
    if v == "sex":
        return s.map({"female": "female", "male": "male"}), ["female", "male"]
    if v == "apoe_genotype":
        lab = s.map(lambda g: str(int(g)) if pd.notna(g) else np.nan)
        return lab, sorted(lab.dropna().unique())
    spec = cfg.STATE_VARS.get(v)
    x = pd.to_numeric(s, errors="coerce")
    if spec is None:
        lab = x.map(lambda a: str(int(a)) if pd.notna(a) else np.nan)
        return lab, [str(int(a)) for a in sorted(x.dropna().unique())]
    w, M = spec
    b = np.minimum(np.floor(x / w), M // w - 1)
    lab = b.map(lambda a: f"{int(a * w)}-{int(a * w + w)}" if pd.notna(a) else np.nan)
    return lab, [f"{int(a * w)}-{int(a * w + w)}" for a in sorted(b.dropna().unique())]


out = {"proteins": P, "vars": []}
for v, name, note, family in VARS:
    lab, order = buckets(v)
    keep = [b for b in order if (lab == b).sum() >= MIN_N]
    hidden = [f"{b} (n={(lab == b).sum()})" for b in order if b not in keep]
    counts = {b: {ds: int(((lab == b) & (meta.dataset == ds)).sum()) for ds in ["ROSMAP", "Diverse", "Banner"]} for b in keep}
    mean, se, nn = [], [], []
    for p in P:
        m, s_, n_ = [], [], []
        for b in keep:
            y = X.loc[(lab == b).values, p].dropna()
            m.append(round(float(y.mean()), 4)); s_.append(round(float(y.std(ddof=1) / np.sqrt(len(y))), 4) if len(y) > 1 else None)
            n_.append(int(len(y)))
        mean.append(m); se.append(s_); nn.append(n_)
    # held-constant path
    cc = con_c[con_c.variable == v]
    held, held_p = [], []
    for p in P:
        cp = cc[cc.protein == p].set_index(["from_bucket", "to_bucket"])
        path, pv = [], []
        if v == "apoe_genotype":
            for b in keep:
                if b == "33": path.append(0.0); pv.append(None)
                elif ("33", b) in cp.index:
                    r = cp.loc[("33", b)]; path.append(round(float(r.beta), 4)); pv.append(round(float(r.p), 5))
                else: path.append(None); pv.append(None)
        else:
            acc = 0.0; broken = False
            for i, b in enumerate(keep):
                if i == 0: path.append(0.0); pv.append(None); continue
                key = (keep[i - 1], b)
                if broken or key not in cp.index: broken = True; path.append(None); pv.append(None); continue
                r = cp.loc[key]; acc += float(r.beta); path.append(round(acc, 4)); pv.append(round(float(r.p), 5))
        held.append(path); held_p.append(pv)
    e = eff_c[eff_c.variable == v].set_index("protein").reindex(P)
    effect = [None if pd.isna(r.beta) else [round(float(r.beta), 4), float(f"{r.p:.3g}"), float(f"{r.q:.3g}")] for r in e.itertuples()]
    ref = "33" if v == "apoe_genotype" else keep[0]
    out["vars"].append(dict(key=v, name=name, note=note, family=family, buckets=keep, counts=counts, hidden=hidden,
                            ref=ref, mean=mean, se=se, n=nn, held=held, held_p=held_p, effect=effect,
                            n_total=int(lab.isin(keep).sum())))

tpl = (HERE / "explorer_template.html").read_text()
html = tpl.replace("/*__DATA__*/null", json.dumps(out, separators=(",", ":")))
(HERE / "expression_explorer.html").write_text(html)
print("wrote expression_explorer.html", f"{len(html) / 1e6:.2f} MB")
