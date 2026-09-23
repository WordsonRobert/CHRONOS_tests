"""Change ONE clinical variable, hold everything else the same, and see how each of the 80 proteins moves.

This is the mirror image of the DAG state pipeline (there: fix one variable, let the rest vary).

For each focal variable v (the 13 trajectory variables, bucketed exactly as in the DAG states):
  protein  =  a  +  effect of v's bucket  +  [every OTHER clinical variable measured for that person, raw values]
                                           +  sex + APOE e4 dosage + log PMI  (+ 5 cell-type scores)
  fitted separately in every cohort (and every group of people with the same set of measured variables =
  a "stratum"), then pooled across strata by inverse-variance weighting.

Rules that keep "hold everything else the same" meaningful:
  * derived variables are never held constant, because they are re-codings of other variables:
        amyA   = binned amyThal         amyAny = (amyCerad > 0)         reag = f(Braak, amyCerad)
    and when a derived variable is the focal one, its parents are not held either (they are the same information)
  * a stratum is used only if n >= 40 and v still moves once the others are fixed
    (remaining variation = 1 - R2 of v on the others; reported)
  * two versions: cell composition held constant ("cell", primary, like DAG2) and not held ("nocell")

Outputs (this folder):
  effects_per_step.csv     v x protein: change in protein (z units) per bucket step of v, SE, p, BH q (per v), n
  bucket_contrasts.csv     v x protein x adjacent buckets (k vs k-1): the step-by-step changes
  effects_matrix_z.csv     80 x 13 z-scores (cell version), for heat maps
  bucket_contrasts_per_stratum.csv   the same, per cohort / stratum before pooling (replication check)
  design_check.csv         per v: strata, people, which variables were held, remaining variation
  marginal_vs_held.csv     how much holding everything else changes the answer vs holding only sex/APOE/PMI
  heatmap_cell.png
  python one_at_a_time.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm, spearmanr

HERE = Path(__file__).resolve().parent
CW = HERE.parent
sys.path.insert(0, str(CW / "DAG2-with_cell_composition" / "methodology"))
if not (CW / "DAG2-with_cell_composition").exists():                       # cloud copy is called DAG2
    sys.path.insert(0, str(CW / "DAG2" / "methodology"))
import config as cfg
cfg.COMPOSITION = "cell"
from data import build_aligned, z_numeric
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 300)

X, Z, meta = build_aligned(); Zn = z_numeric(Z)
PROT = list(X.columns)
CELL = ["cell_neuron", "cell_astrocyte", "cell_microglia", "cell_oligodendrocyte", "cell_endothelial"]
Zn["log_pmi"] = np.log(Zn["pmi"])
Zn["apoe_e4"] = Zn["apoe_genotype"]                       # z_numeric already codes APOE as e4 dosage 0/1/2
apoe_raw = Z["apoe_genotype"]

CLIN = ["cts_mmse30_lv", "cts_mmse30_first_ad_dx", "braaksc", "ceradsc", "amyThal", "amyA", "amyCerad",
        "amyAny", "reag", "PlaqueTotal", "TangleTotal"]
FOCAL = ["sex", "apoe_genotype"] + CLIN
PARENTS = {"amyA": ["amyThal"], "amyAny": ["amyCerad"], "reag": ["braaksc", "amyCerad"]}
DERIVED = set(PARENTS)
NEVER_HELD = DERIVED | {"cts_mmse30_first_ad_dx"}         # MMSE at first AD diagnosis exists for 126 people only
BOXES = cfg.STATE_VARS                                      # same buckets as the DAG states


def bucket(v, x):
    """Ordered bucket index (0, 1, 2, ...) and label, same buckets as the DAG states."""
    spec = BOXES.get(v)
    if spec is None:
        return x
    w, M = spec
    return np.minimum(np.floor(x / w), M // w - 1)          # top value goes into the last box


def labels(v, b):
    spec = BOXES.get(v)
    if spec is None: return str(int(b)) if v != "sex" else ("male" if b == 1 else "female")
    w, _ = spec; lo = int(b * w); return f"{lo}-{lo + w}"


def focal_values(v):
    if v == "apoe_genotype":            # per-step = per e4 allele; buckets = genotypes (reference 33)
        return Zn["apoe_e4"], apoe_raw.map(lambda g: str(int(g)) if pd.notna(g) else np.nan)
    if v == "sex":
        return Zn["sex"], Zn["sex"]
    b = bucket(v, Zn[v].astype(float))
    return b, b


def fit(y, x_step, dummies, C):
    """OLS; returns (step beta, se), dummy betas + covariance (for contrasts)."""
    out = {}
    for name, D in [("step", x_step[:, None]), ("dum", dummies)]:
        if D is None or D.shape[1] == 0: continue
        M = np.column_stack([np.ones(len(y)), D, C])
        b, *_ = np.linalg.lstsq(M, y, rcond=None); r = y - M @ b
        dof = len(y) - np.linalg.matrix_rank(M)
        V = (r @ r / dof) * np.linalg.pinv(M.T @ M)
        k = D.shape[1]
        out[name] = (b[1:1 + k], V[1:1 + k, 1:1 + k])
    return out


GLOBAL_ORDER = {}


def run(version):
    rows_step, rows_con, rows_des = [], [], []
    for v in FOCAL:
        step_all, bval_all = focal_values(v)
        GLOBAL_ORDER[v] = sorted(bval_all.dropna().unique(), key=lambda s: str(s) if v == "apoe_genotype" else float(s))
        for ds in ["ROSMAP", "Diverse", "Banner"]:
            m0 = (meta.dataset == ds).values & step_all.notna().values
            if m0.sum() < 40: continue
            # everything else measured in this cohort that may be held constant
            others = [c for c in CLIN if c != v and c not in NEVER_HELD and c not in PARENTS.get(v, [])
                      and Zn.loc[m0, c].notna().any()]
            base = ["sex", "log_pmi"] + (["apoe_e4"] if v != "apoe_genotype" else []) + (CELL if version == "cell" else [])
            base = [c for c in base if c != v]
            # strata = people sharing the same set of measured "others"
            pattern = Zn.loc[m0, others].notna().apply(lambda r: "".join("1" if x else "0" for x in r), axis=1)
            for pat, idx in pattern.groupby(pattern).groups.items():
                held = [c for c, ok in zip(others, pat) if ok == "1"]
                keys = [k for k in idx if Zn.loc[k, base].notna().all()]
                if len(keys) < 40: continue
                xs = step_all.loc[keys].values.astype(float)
                bv = bval_all.loc[keys]
                if len(np.unique(xs)) < 2: continue
                C = Zn.loc[keys, held + base].values.astype(float)
                # remaining variation of the focal variable after holding the others
                Mc = np.column_stack([np.ones(len(keys)), C])
                res = xs - Mc @ np.linalg.lstsq(Mc, xs, rcond=None)[0]
                remain = float(res.var() / xs.var())
                levels = sorted(bv.unique(), key=lambda s: (float(s) if v != "apoe_genotype" else {"33": -1}.get(s, float(s))))
                ref = "33" if v == "apoe_genotype" and "33" in levels else levels[0]
                lev = [l for l in levels if l != ref]
                cnt = bv.value_counts()
                lev = [l for l in lev if cnt[l] >= 5]                                 # skip buckets with < 5 people
                keep = bv.isin([ref] + lev).values                                   # drop people in those tiny buckets
                keys = [k for k, kk in zip(keys, keep) if kk]
                xs, bv, C = xs[keep], bv[keep], C[keep]
                if len(keys) < 40 or len(np.unique(xs)) < 2: continue
                Mc = np.column_stack([np.ones(len(keys)), C])
                res = xs - Mc @ np.linalg.lstsq(Mc, xs, rcond=None)[0]
                remain = float(res.var() / xs.var())
                D = np.column_stack([(bv.values == l).astype(float) for l in lev]) if lev else np.zeros((len(keys), 0))
                rows_des.append(dict(version=version, variable=v, cohort=ds, stratum_n=len(keys),
                                     held=", ".join(held) if held else "(none measured)",
                                     remaining_variation=round(remain, 3),
                                     buckets=" ".join(f"{labels(v, float(l)) if v not in ('apoe_genotype',) else l}:{cnt[l]}" for l in [ref] + lev),
                                     used=remain >= 0.05))
                if remain < 0.05: continue            # v can't move once everything else is fixed
                for p in PROT:
                    y = X.loc[keys, p].values.astype(float)
                    ok = np.isfinite(y)
                    if ok.sum() < 40: continue
                    f = fit(y[ok], xs[ok], D[ok][:, (D[ok].sum(0) >= 5)] if D.shape[1] else None, C[ok])
                    bs, Vs = f["step"]
                    rows_step.append(dict(variable=v, cohort=ds, protein=p, beta=bs[0], se=np.sqrt(Vs[0, 0]), n=int(ok.sum())))
                    if "dum" in f:
                        lv = [l for l, s in zip(lev, D[ok].sum(0) >= 5) if s]
                        bd, Vd = f["dum"]
                        chain = [ref] + lv
                        for a in range(1, len(chain)):
                            # contrast bucket a vs a-1  (reference has coefficient 0)
                            w = np.zeros(len(lv)); w[a - 1] = 1
                            if a >= 2: w[a - 2] = -1
                            if v == "apoe_genotype":            # contrasts vs 33 for every genotype
                                w = np.zeros(len(lv)); w[a - 1] = 1; prev = ref
                            else:
                                prev = chain[a - 1]
                                g = GLOBAL_ORDER[v]
                                if g.index(chain[a]) - g.index(prev) != 1: continue   # never jump over a missing bucket
                            rows_con.append(dict(variable=v, cohort=ds, protein=p,
                                                 from_bucket=labels(v, float(prev)) if v != "apoe_genotype" else prev,
                                                 to_bucket=labels(v, float(chain[a])) if v != "apoe_genotype" else chain[a],
                                                 beta=w @ bd, se=np.sqrt(w @ Vd @ w), n=int(ok.sum())))
    return pd.DataFrame(rows_step), pd.DataFrame(rows_con), pd.DataFrame(rows_des)


def pool(df, by):
    def agg(g):
        w = 1 / g.se ** 2; b = (w * g.beta).sum() / w.sum(); se = 1 / np.sqrt(w.sum())
        return pd.Series(dict(beta=b, se=se, z=b / se, p=2 * norm.sf(abs(b / se)), n=g.n.sum(), strata=len(g),
                              strata_same_sign=int((np.sign(g.beta) == np.sign(b)).sum())))
    return df.groupby(by, sort=False).apply(agg, include_groups=False).reset_index()


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


allS, allC, allD, allRaw = [], [], [], []
for version in ["cell", "nocell"]:
    S, Cn, D = run(version)
    PS = pool(S, ["variable", "protein"]); PS["q"] = PS.groupby("variable").p.transform(bh); PS.insert(0, "version", version)
    PC = pool(Cn, ["variable", "from_bucket", "to_bucket", "protein"]); PC["q"] = PC.groupby(["variable", "from_bucket", "to_bucket"]).p.transform(bh)
    PC.insert(0, "version", version)
    allS.append(PS); allC.append(PC); allD.append(D)
    Cn.insert(0, "version", version); allRaw.append(Cn)
PS = pd.concat(allS); PC = pd.concat(allC); D = pd.concat(allD)
PS.to_csv(HERE / "effects_per_step.csv", index=False)
PC.to_csv(HERE / "bucket_contrasts.csv", index=False)
D.to_csv(HERE / "design_check.csv", index=False)
R = pd.concat(allRaw); R["p"] = 2 * norm.sf(abs(R.beta / R.se))
R.to_csv(HERE / "bucket_contrasts_per_stratum.csv", index=False)       # same contrasts, before pooling (per cohort/stratum)
Mz = PS[PS.version == "cell"].pivot(index="protein", columns="variable", values="z").reindex(columns=[v for v in FOCAL if v in set(PS.variable)])
Mz.to_csv(HERE / "effects_matrix_z.csv")

# ---------------- marginal (hold only sex / APOE / PMI) vs everything held ----------------
def marginal():
    rows = []
    for v in FOCAL:
        step_all, _ = focal_values(v)
        for ds in ["ROSMAP", "Diverse", "Banner"]:
            base = ["sex", "log_pmi"] + (["apoe_e4"] if v != "apoe_genotype" else []) + CELL
            base = [c for c in base if c != v]
            keys = X.index[(meta.dataset == ds).values & step_all.notna().values & Zn[base].notna().all(axis=1).values]
            if len(keys) < 40 or step_all.loc[keys].nunique() < 2: continue
            C = Zn.loc[keys, base].values.astype(float); xs = step_all.loc[keys].values.astype(float)
            for p in PROT:
                y = X.loc[keys, p].values.astype(float); ok = np.isfinite(y)
                b, V = fit(y[ok], xs[ok], None, C[ok])["step"]
                rows.append(dict(variable=v, cohort=ds, protein=p, beta=b[0], se=np.sqrt(V[0, 0]), n=int(ok.sum())))
    M = pool(pd.DataFrame(rows), ["variable", "protein"]); M["q"] = M.groupby("variable").p.transform(bh)
    return M
MG = marginal()
H = PS[PS.version == "cell"]
cmp = []
for v in Mz.columns:
    a = MG[MG.variable == v].set_index("protein"); b = H[H.variable == v].set_index("protein").reindex(a.index)
    cmp.append(dict(variable=v, sig_q05_hold_only_sex_apoe_pmi_cell=int((a.q < .05).sum()),
                    sig_q05_hold_everything=int((b.q < .05).sum()),
                    corr_of_80_effects=round(float(np.corrcoef(a.beta, b.beta)[0, 1]), 2),
                    median_abs_effect_ratio=round(float(np.median(np.abs(b.beta)) / np.median(np.abs(a.beta))), 2)))
CMP = pd.DataFrame(cmp); CMP.to_csv(HERE / "marginal_vs_held.csv", index=False)

# ---------------- printout ----------------
NAMES = {"sex": "sex (male vs female)", "apoe_genotype": "APOE (per e4 allele)", "cts_mmse30_lv": "MMSE last visit (per box of 6, higher = better)",
         "cts_mmse30_first_ad_dx": "MMSE at AD dx (per box of 6, higher = better)", "braaksc": "Braak (per stage)", "ceradsc": "CERAD R/B (per step)",
         "amyThal": "Thal phase (per phase)", "amyA": "A-score (per step)", "amyCerad": "CERAD Diverse (per step)",
         "amyAny": "any amyloid (yes vs no)", "reag": "NIA-Reagan (per step)", "PlaqueTotal": "plaque total (per box of 3)",
         "TangleTotal": "tangle total (per box of 3)"}
print("=== DESIGN: what was held constant, and can the variable still move? (cell version) ===")
dd = D[D.version == "cell"]
print(dd[["variable", "cohort", "stratum_n", "remaining_variation", "used", "held"]].to_string(index=False))

print("\n=== EFFECTS (cell composition held too). z-units of protein per step of the variable ===")
for v in Mz.columns:
    s = H[H.variable == v].sort_values("beta")
    sig = s[s.q < 0.05]
    t = s.set_index("protein").loc[["TMED2", "TMED9", "TMED10"]]
    print(f"\n{NAMES[v]}: n={int(s.n.median())}  proteins changed at FDR 5%: {len(sig)}/80 "
          f"(down {int((sig.beta < 0).sum())}, up {int((sig.beta > 0).sum())})")
    if len(sig):
        dn = sig[sig.beta < 0].head(6); up = sig[sig.beta > 0].tail(6)[::-1]
        if len(dn): print("   down:", ", ".join(f"{r.protein} {r.beta:+.2f}" for r in dn.itertuples()))
        if len(up): print("   up:  ", ", ".join(f"{r.protein} {r.beta:+.2f}" for r in up.itertuples()))
    print("   TMED:", ", ".join(f"{p} {r.beta:+.3f} (p={r.p:.2g}, q={r.q:.2g})" for p, r in t.iterrows()))

print("\n=== does holding everything else matter? (vs holding only sex, APOE, PMI, cell types) ===")
print(CMP.to_string(index=False))
nc = PS[PS.version == "nocell"].set_index(["variable", "protein"])
print("\n=== cell composition NOT held: proteins at FDR 5% per variable ===")
print((nc.q < .05).groupby(level=0).sum().reindex(Mz.columns).to_string())

# ---------------- heat map ----------------
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    Q = H.pivot(index="protein", columns="variable", values="q").reindex(index=Mz.index, columns=Mz.columns)
    order = Mz.fillna(0).abs().sum(axis=1).sort_values(ascending=False).index
    Mp, Qp = Mz.loc[order], Q.loc[order]
    fig, ax = plt.subplots(figsize=(9, 20))
    lim = np.nanpercentile(np.abs(Mp.values), 98)
    im = ax.imshow(Mp.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(Mp.shape[0]):
        for j in range(Mp.shape[1]):
            if Qp.iloc[i, j] < 0.05: ax.text(j, i, "*", ha="center", va="center", fontsize=8)
    ax.set_xticks(range(Mp.shape[1])); ax.set_xticklabels([NAMES[c] for c in Mp.columns], rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(Mp.shape[0])); ax.set_yticklabels(Mp.index, fontsize=7)
    for lab in ax.get_yticklabels():
        if lab.get_text().startswith("TMED"): lab.set_fontweight("bold"); lab.set_color("darkgreen")
    ax.set_title("Change ONE variable, hold everything else\n(z-score of protein change per step; * = FDR < 5%; red = up, blue = down)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="z")
    fig.tight_layout(); fig.savefig(HERE / "heatmap_cell.png", dpi=150)
    print("\nsaved heatmap_cell.png")
except ImportError:
    print("matplotlib not installed - skipped heat map")
