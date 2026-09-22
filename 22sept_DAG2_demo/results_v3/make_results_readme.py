"""Builds a results README from a results folder (tables are generated, not hand-copied).

python make_results_readme.py                                   # old prior  -> README_RESULTS.md
python make_results_readme.py --results results_v3 --prior CHRONOS_80x80_confidence_v3.xlsx \
       --out README_RESULTS_v3.md --label v3 --compare_with results
Needs {results}/extra_stats.json (from extra_stats.py).
"""
import argparse, json, os, sys, numpy as np, pandas as pd
from pathlib import Path
_ap = argparse.ArgumentParser()
_ap.add_argument("--results", default="results"); _ap.add_argument("--prior", default="../CHRONOS_80x80_confidence.xlsx")
_ap.add_argument("--out", default="README_RESULTS.md"); _ap.add_argument("--label", default="")
_ap.add_argument("--compare_with", default=None, help="another results folder to compare against (e.g. the old prior)")
A = _ap.parse_args()
sys.argv = [sys.argv[0]]; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chronos_regions as cr
R = Path(A.results); DS = ["ROSMAP", "Diverse", "Banner", "BannerLFQ"]; MAIN = DS[:3]
X = json.loads((R / "extra_stats.json").read_text())
PN, PC = cr.load_prior(A.prior)
LABEL = f" ({A.label} prior)" if A.label else ""
S = {d: json.loads((R / f"{d}_summary.json").read_text()) for d in DS}
V = {d: pd.read_csv(R / f"{d}_variables.csv") for d in DS}
G = {d: pd.read_csv(R / f"{d}_regions.csv") for d in DS}
E = {d: pd.read_csv(R / f"{d}_causal_edges.csv") for d in DS}
FIX = {d: pd.read_csv(R / f"{d}_regions_usable_null_check.csv") for d in DS if (R / f"{d}_regions_usable_null_check.csv").exists()}
out = []
w = out.append

def md(df, fl=4):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f":
            isz = c in ("z", "z_fixed", "z_orig")
            ispq = (c.startswith("q") or c.startswith("p_") or c == "p" or c == "q")
            def f(x, isz=isz, ispq=ispq):
                if pd.isna(x): return "—"
                if isz: return f"{x:.2f}"
                if ispq and x != 0 and abs(x) < 1e-3: return f"{x:.2e}"
                return f"{x:.{fl}f}"
            df[c] = df[c].map(f)
        else:
            df[c] = df[c].map(lambda x: "—" if (x is None or (isinstance(x, float) and np.isnan(x)) or str(x) in ("nan", "")) else x)
    cols = list(df.columns)
    s = "| " + " | ".join(cols) + " |\n|" + "|".join("---" for _ in cols) + "|\n"
    for _, r in df.iterrows():
        s += "| " + " | ".join(str(r[c]) for c in cols) + " |\n"
    return s

w(f"""# CHRONOS causal graphs — all results{LABEL}

Everything the pipeline produced, dataset by dataset, with as little interpretation as possible.
How each number is computed is in `README_HOW_IT_WORKS.md`; section numbers like (HIW 4.7) point there.
All tables below were generated from the files in `{R}/` by `make_results_readme.py`,
using the confidence matrix `{os.path.basename(A.prior)}`.

Region results are given twice: the pipeline's original scoring, and a corrected re-scoring that
draws the random comparison groups only from people who have adjusted data (HIW 9.1). Section 5.1
has both counts.

Datasets: **ROSMAP**, **Diverse** (AMP-AD Diverse Cohorts), **Banner** (TMT) are the three
independent graphs. **BannerLFQ** is the same Banner people measured with a second method
(label-free), included as a technical check only.

---
""")

# 1 data as loaded
w("\n## 1. Data as loaded\n\n")
w(f"""| | ROSMAP | Diverse | Banner (TMT) | BannerLFQ |
|---|---|---|---|---|
| proteomics samples in file | 400 | 1,086 | 220 | 225 |
| proteins in file (unique gene symbols) | 8,252 | 9,152 | 9,728 | 5,041 |
| samples matched to a person | 400 | 1,086 | 198 | 198 |
| people after averaging replicates | 400 | 980 | 198 | 198 |
| people excluded as shared with another dataset | 0 | 16 (ROSMAP donors) | 8 (possible Diverse donors) | 8 |
| **people analysed** | **400** | **964** | **190** | **190** |
| people with usable data after adjustment for Z* | {X['ROSMAP']['usable_after_adjustment']} | **{X['Diverse']['usable_after_adjustment']}** | {X['Banner']['usable_after_adjustment']} | {X['BannerLFQ']['usable_after_adjustment']} |
| graph nodes found (of 80) | 80 | 80 | 80 | 58 |
| missing fraction, 80 nodes × people (raw) | 5.6% | 1.0% | 4.9% | 10.9% |
| nodes with no missing values | 58 | 69 | 66 | 42 |
| people with all nodes measured | 40 | 525 | 0 | 0 |
| median people per protein pair (raw) | 400 | 964 | 190 | 190 |
| smallest people per protein pair (raw) | 120 | 593 | 8 | 0 |
| "rest-of-proteome" proteins used for PCs (≥95% observed) | 5,683 | 7,286 | 6,184 | 2,649 |
| candidate trajectory variables | 29 | 33 | 25 | 25 |

Most-missing graph nodes (fraction of people missing):

| ROSMAP | Diverse | Banner | BannerLFQ |
|---|---|---|---|
| SIGMAR1 0.46 | SREBF2 0.32 | GORASP1 0.73 | SEC23B 0.95 |
| AP4B1 0.42 | AP4M1 0.18 | SREBF2 0.64 | MAPK14 0.94 |
| AP4M1 0.42 | GORASP1 0.17 | AP4M1 0.58 | SEC16A 0.85 |
| AP3B2 0.40 | ABCA7 0.05 | ABCA7 0.40 | IDE 0.73 |
| SREBF2 0.36 | SIGMAR1 0.03 | AP4E1 0.36 | GGA1 0.65 |
| GORASP1 0.34 | AP4E1 0.01 | AP4B1 0.36 | COPG2 0.52 |
| ABCA7 0.30 | MAP2K3 0.01 | MAPKAPK2 0.22 | PREB 0.35 |
| RTN3 0.28 | PLCG2 0.01 | SIGMAR1 0.18 | ITPR2 0.34 |

Nodes absent from BannerLFQ (22): ABCA7, ADAM17, AP4B1, AP4E1, AP4M1, AP4S1, BACE1, ECE1, GGA2,
GORASP1, MAP2K3, MAPKAPK2, NCSTN, PLCG2, PLG, PSEN1, SEC24A, SEC24D, SIGMAR1, SREBF2, SYK, VLDLR.

Candidate trajectory variables (besides 5 cell-type scores and 10 proteome PCs):
- **ROSMAP (14):** Study, msex, educ, apoe_genotype, age_at_visit_max, age_first_ad_dx, age_death, cts_mmse30_first_ad_dx, cts_mmse30_lv, pmi, braaksc, ceradsc, cogdx, dcfdx_lv
- **Diverse (18):** dataContributionGroup, cohort, sex, race, isHispanic, ageDeath, PMI, apoeGenotype, amyThal, amyA, amyCerad, Braak, mayoDx, amyAny, bScore, reag, ADoutcome, derivedOutcomeBasedOnMayoDx
- **Banner, BannerLFQ (10):** sex, ageDeath, apoeGenotype, pmi, diagnosis, CERAD, Braak, PlaqueTotal, TangleTotal, lastMMSE

Cell-type score markers missing from a dataset (all other listed markers were used): ROSMAP microglia
lacked CD68, CSF1R, CX3CR1 (used AIF1, ITGAM, P2RY12, HLA-DRA); Diverse microglia lacked CD68; Banner
neuron lacked RBFOX3 and Banner microglia lacked CD68. In BannerLFQ the microglia score exists for only 32 people and the
endothelial score for 101.

---
""")

# 2 prior
d_ = PC >= cr.MIN_PRIOR; U_ = np.triu(np.maximum(PC, PC.T) >= cr.MIN_PRIOR, 1); B_ = np.triu(d_ & d_.T, 1)
EQ_ = np.triu(d_ & d_.T & (np.abs(PC - PC.T) < 1e-9), 1)
cd_ = pd.Series(PC[d_]).round(2).value_counts().sort_index()
od_ = sorted(zip(d_.sum(1), PN), reverse=True)[:5]; id_ = sorted(zip(d_.sum(0), PN), reverse=True)[:5]
w(f"""## 2. The confidence graph C

| | value |
|---|---|
| file | `{os.path.basename(A.prior)}` (first sheet), raw value / 5, diagonal set to 0 |
| directed prior edges (C ≥ 0.05) | {int(d_.sum())} |
| unordered pairs with a prior edge | {int(U_.sum())} of 3,160 |
| one-way pairs | {int(U_.sum() - B_.sum())} |
| two-way pairs | {int(B_.sum())} ({int(EQ_.sum())} with equal C both ways) |
| C values among edges | {", ".join(f"{k:g}: {v}" for k, v in cd_.items())} |
| largest out-degree | {", ".join(f"{p} {n}" for n, p in od_)} |
| largest in-degree | {", ".join(f"{p} {n}" for n, p in id_)} |

""")
_xl = pd.ExcelFile(A.prior)
if "Edge_List" in _xl.sheet_names and "Provenance" in pd.read_excel(A.prior, sheet_name="Edge_List", nrows=1).columns:
    EL = pd.read_excel(A.prior, sheet_name="Edge_List"); EL = EL[EL.Source != EL.Target]
    ELe = EL[EL.Confidence >= 1]
    w("Provenance of the directed prior edges (from the `Edge_List` sheet):\n\n" +
      "| provenance | edges |\n|---|---|\n" + "".join(f"| {k} | {v} |\n" for k, v in ELe.Provenance.value_counts().items()) + "\n")
    if "Sign_80x80" in _xl.sheet_names:
        SG = pd.read_excel(A.prior, sheet_name="Sign_80x80", index_col=0).values.astype(float); np.fill_diagonal(SG, 0)
        sv = pd.Series(SG[d_]).value_counts()
        w(f"`Sign_80x80` among directed prior edges: +1 {int(sv.get(1, 0))}, 0 {int(sv.get(0, 0))}, −1 {int(sv.get(-1, 0))}. "
          "The pipeline does not read this sheet.\n\n")
    neg = pd.read_excel(A.prior, index_col=0).values; negc = pd.Series(neg[~np.eye(80, dtype=bool)]).value_counts().sort_index()
    w("Raw off-diagonal values in the matrix: " + ", ".join(f"{int(k)}: {v}" for k, v in negc.items()) +
      ". Values ≤ 0 are non-edges for the pipeline.\n\n")
w("---\n\n")

# 3 signal
w("## 3. Signal: how well each dataset recovers C (AUROC, HIW 4.4)\n\n0.5 = no relation to C. \"Signed\" (used by the pipeline) ranks pairs by ρ; \"unsigned\" by abs(ρ).\n\n")
w("| | " + " | ".join(DS) + " |\n|---|" + "---|" * len(DS) + "\n")
for lab, key, kk in [("signed AUROC, raw", "raw", "signed_auroc"), ("signed AUROC, adjusted for Z*", "adjusted", "signed_auroc"),
                     ("unsigned AUROC, raw", "raw", "unsigned_auroc"), ("unsigned AUROC, adjusted", "adjusted", "unsigned_auroc"),
                     ("pairs scored (n ≥ 20)", "raw", "pairs_scored"), ("prior pairs among them", "raw", "prior_pairs_scored")]:
    vals = [X[d][key][kk] for d in DS]
    fmt = (lambda v: f"**{v:.4f}**") if lab.startswith("signed AUROC, adj") else (lambda v: f"{v:.4f}" if isinstance(v, float) else f"{v:,}")
    w(f"| {lab} | " + " | ".join(fmt(v) for v in vals) + " |\n")
w("\nCorrelation between prior-linked pairs vs other pairs:\n\n")
hdr = [f"{d} {k}" for d in DS for k in ("raw", "adj")]
w("| | " + " | ".join(hdr) + " |\n|---|" + "---|" * len(hdr) + "\n")
for lab, kk, sg in [("prior pairs: mean ρ", "prior_mean_rho", True), ("prior pairs: fraction ρ > 0", "prior_frac_pos", False), ("prior pairs: mean abs(ρ)", "prior_mean_abs", False),
                    ("other pairs: mean ρ", "other_mean_rho", True), ("other pairs: fraction ρ > 0", "other_frac_pos", False), ("other pairs: mean abs(ρ)", "other_mean_abs", False)]:
    vals = [X[d][k][kk] for d in DS for k in ("raw", "adjusted")]
    w(f"| {lab} | " + " | ".join((f"{v:+.3f}" if sg else f"{v:.3f}") for v in vals) + " |\n")
w("\n---\n\n")

# 4 variables
w("## 4. Trajectory variables (HIW 4.6)\n\n### 4.1 Selected Z*, in order\n\n")
for d in DS:
    v = V[d]; sel = S[d]["Z_star"]
    rows = []
    for k, name in enumerate(sel, 1):
        r = v[(v.step == k) & (v.variable == name)].iloc[0]
        rows.append(dict(step=k, variable=name, kind=r.kind, AUROC_after=r.auroc, shuffled_mean=r.null_mean,
                         gain=r.gain_vs_current, z=r.z_vs_shuffled, people=int(r.n_subjects)))
    w(f"**{d}** — start AUROC {S[d]['auroc_raw']:.4f}, end {S[d]['auroc_adjusted']:.4f}\n\n")
    w(md(pd.DataFrame(rows)) + "\n")
w("""Proteome PC descriptions (share of rest-of-proteome variance; cell-type score it tracks most, Spearman ρ):

| PC | ROSMAP | Diverse | Banner |
|---|---|---|---|
| 1 | 9.9%, neuron +0.59 | 23.3%, neuron −0.89 | 12.2%, oligodendrocyte −0.59 |
| 2 | 6.1%, neuron +0.27 | 8.8%, astrocyte +0.67 | 6.9%, oligodendrocyte −0.35 |
| 3 | 4.9%, oligodendrocyte −0.47 | 4.0%, oligodendrocyte +0.25 | 5.6%, neuron +0.23 |
| 4 | 4.6%, oligodendrocyte −0.48 | 3.2%, endothelial −0.17 | 3.6%, endothelial −0.38 |
| 5 | 3.2%, astrocyte −0.47 | 3.0%, microglia −0.24 | 3.1%, neuron −0.31 |
| 6 | 2.9%, oligodendrocyte −0.13 | 2.7%, microglia +0.21 | 2.6%, astrocyte −0.20 |
| 7 | 2.4%, microglia −0.23 | 2.0%, oligodendrocyte −0.22 | 2.2%, microglia −0.30 |
| 8 | 2.3%, oligodendrocyte +0.29 | 1.9%, astrocyte +0.14 | 2.1%, oligodendrocyte −0.08 |
| 9 | 2.1%, astrocyte +0.27 | 1.8%, astrocyte −0.37 | 2.0%, endothelial −0.27 |
| 10 | 1.7%, astrocyte +0.32 | 1.5%, astrocyte +0.24 | 1.8%, astrocyte +0.39 |

BannerLFQ PC1: 24.1% of variance, tracks microglia (ρ = −0.42); PC9: 1.2%, microglia (ρ = −0.20).

""")
w("### 4.2 Every candidate at step 1 (each variable alone vs no adjustment)\n\n"
  "`gain` = AUROC change vs unadjusted; `z` = vs 10 shuffled copies; eligible = z > 2 and gain > 0.0001.\n\n")
for d in DS:
    x = V[d][V[d].step == 1].sort_values("auroc", ascending=False)
    x = x[["variable", "kind", "auroc", "null_mean", "gain_vs_current", "z_vs_shuffled", "n_subjects"]].rename(
        columns={"null_mean": "shuffled_mean", "gain_vs_current": "gain", "z_vs_shuffled": "z", "n_subjects": "people"})
    x["z"] = x["z"].round(2)
    w(f"**{d}** (unadjusted AUROC {S[d]['auroc_raw']:.4f})\n\n" + md(x) + "\n")
w("### 4.3 Clinical/pathology variables at every step\n\nz vs shuffled for each non-proteome variable at each selection step (the full per-step lists are in `{ds}_variables.csv`). Bold = selected at that step.\n\n")
for d in DS:
    v = V[d]; cl = v[v.kind != "proteome-derived"]
    piv = cl.pivot_table(index="variable", columns="step", values="z_vs_shuffled", aggfunc="first").round(2)
    piv.columns = [f"step {c}" for c in piv.columns]
    sel = S[d]["Z_star"]
    piv = piv.reset_index()
    for c in piv.columns[1:]:
        piv[c] = piv[c].map(lambda x: "" if pd.isna(x) else f"{x:.2f}").astype(object)
    for k, name in enumerate(sel, 1):
        if name in set(piv.variable) and f"step {k}" in piv.columns:
            i = piv.index[piv.variable == name][0]; piv.loc[i, f"step {k}"] = f"**{piv.loc[i, f'step {k}']}**"
    w(f"**{d}**\n\n" + md(piv.fillna("")) + "\n")
w("---\n\n")

# 5 regions
w("## 5. Regions (HIW 4.7)\n\n### 5.1 Overview\n\n")
rows = []
for d in DS:
    g = G[d]
    row = dict(dataset=d, regions_tested=len(g), depth1=int((g.depth == 1).sum()), depth2=int((g.depth == 2).sum()),
               z_above_2=int((g.z > 2).sum()), z_below_minus2=int((g.z < -2).sum()), STRONG_original=int(g.strong.sum()),
               min_q_original=round(g.q_fdr.min(), 3), whole_cohort_AUROC=S[d]["auroc_adjusted"])
    if d in FIX:
        f = FIX[d]; row.update(STRONG_corrected=str(int(f.strong_fixed.sum())), min_q_corrected=round(f.q_fixed.min(), 3),
                               z_above_2_corrected=str(int((f.z_fixed > 2).sum())))
    rows.append(row)
w(md(pd.DataFrame(rows).fillna("—"), 3) + "\n")
_miss = [d for d in DS if d not in FIX]
_us = ", ".join(f"{d} {X[d]['usable_after_adjustment']}/{X[d]['people']}" for d in DS)
w(f"""`original` = the pipeline's output (`{{ds}}_regions.csv`). `corrected` = same regions re-scored
with random subsets drawn only from people who have adjusted data and matched to the region's
usable size (`{{ds}}_regions_usable_null_check.csv`, HIW 9.1). People with adjusted data / people analysed:
{_us}. Where the two numbers are equal, the two versions differ only by random draws.""" +
  (f" The corrected check was not run for: {', '.join(_miss)}." if _miss else "") + "\n\n")

cols = ["region", "depth", "n", "auroc", "null_mean", "z", "q_fdr"]
fcols = ["region", "n_listed", "n_usable", "auroc", "null_usable", "z_fixed", "q_fixed"]
for d in DS:
    g = G[d]
    w(f"### 5.{DS.index(d)+2} {d}\n\nWhole-cohort adjusted AUROC: {S[d]['auroc_adjusted']:.4f}\n\n")
    w("**Top 15 by z (original)**\n\n" + md(g.sort_values("z", ascending=False).head(15)[cols]) + "\n")
    if d in FIX:
        w("**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data\n\n"
          + md(FIX[d].sort_values("z_fixed", ascending=False).head(15)[fcols]) + "\n")
    w("**Bottom 8 by z (C least visible, original)**\n\n" + md(g.sort_values("z").head(8)[cols]) + "\n")
    d1 = g[g.depth == 1].copy(); d1["var"] = d1.region.str.split(" ").str[0]
    cl = d1[~d1["var"].str.startswith("celltype")].sort_values(["var", "region"])
    w("**All single-variable clinical/pathology regions (original)**\n\n" + md(cl[["region", "n", "auroc", "null_mean", "z", "q_fdr"]]) + "\n")
    ct = d1[d1["var"].str.startswith("celltype")].sort_values(["var", "region"])
    w("**All single-variable cell-type regions (original)**\n\n" + md(ct[["region", "n", "auroc", "null_mean", "z", "q_fdr"]]) + "\n")
    if d == "Diverse":
        st = g[g.strong][cols]
        if len(st):
            w(f"**The {len(st)} regions flagged STRONG in the original run** ({int(FIX[d].strong_fixed.sum())} remain STRONG with the corrected null)\n\n" + md(st) + "\n")
        nanr = g[g.auroc.isna()][["region", "n", "auroc", "null_mean", "p_vs_random_same_size", "q_fdr"]]
        if len(nanr): w("**Regions with no usable people after adjustment** (AUROC NaN, p set to the 1/201 floor; HIW 9.2)\n\n" + md(nanr) + "\n")
    if S[d].get("best_region"):
        w(f"`{d}_best_region_W.csv` was refitted in `{S[d]['best_region']}`, the original top STRONG region.\n\n")
w("---\n\n")

# 6 graphs
w("## 6. The graphs (HIW 4.8–4.10)\n\n### 6.1 Edge counts and sizes\n\n")
dirc = {d: tuple(X[d]["arrow_setting"][k] for k in ("one_way_agrees", "one_way_disagrees_downweighted", "two_way_picked", "low_n_placeholder")) for d in DS}
rows = []
for d in DS:
    e = E[d]; s = e[e.data_supported == 1]; wp = e[e.W > 0]
    rows.append(dict(dataset=d, prior_edges=len(e), testable_n20=int((e.n_pairwise >= 20).sum()), W_positive=len(wp),
                     q_below_05_any_direction=int((e.q_fdr < 0.05).sum()), data_supported=len(s),
                     pct_of_prior_edges=round(100 * len(s) / len(e), 1), median_n_W_positive=int(wp.n_pairwise.median()),
                     min_n_W_positive=int(wp.n_pairwise.min())))
w(md(pd.DataFrame(rows)) + "\n")
rows = []
for d in DS:
    e = E[d]; s = e[e.data_supported == 1]; wp = e[e.W > 0]
    rows.append(dict(dataset=d, W_median_all=wp.W.median(), W_max=e.W.max(), W_median_supported=s.W.median(), W_min_supported=s.W.min(),
                     beta_median_supported=s.beta.median(), beta_min=s.beta.min(), beta_max=s.beta.max(),
                     frac_beta_positive=round((s.beta > 0).mean(), 3), n_abs_beta_above_2=int((s.beta.abs() > 2).sum())))
w(md(pd.DataFrame(rows)) + "\n")
w("### 6.2 How each prior pair's arrow and weight were set (HIW 4.8)\n\n"
  "| | one-way pair, direction test agrees with C (full weight) | one-way pair, direction test disagrees (kept in C's direction at 0.3 weight) | two-way pair (direction test picks the arrow) | n < 20 (placeholder 0.1·C) |\n|---|---|---|---|---|\n")
for d in DS:
    a, b, c, z = dirc[d]; w(f"| {d} | {a} | {b} | {c} | {z} |\n")
w("\nShare of one-way pairs where the direction test agrees with C: " + ", ".join(
    (f"{d} {100*dirc[d][0]/(dirc[d][0]+dirc[d][1]):.0f}%" if (dirc[d][0] + dirc[d][1]) else f"{d} — (no one-way pairs)") for d in DS) + ".\n\n")
w("### 6.3 Supported edges by prior confidence\n\n")
rows = []
allc = E["ROSMAP"].prior_C.value_counts().sort_index()
for c in allc.index:
    r = {"prior C": c, "prior edges": int(allc[c])}
    for d in DS:
        s = E[d][E[d].data_supported == 1]; k = int((s.prior_C == c).sum()); r[d] = f"{k} ({100*k/allc[c]:.0f}%)"
    rows.append(r)
w(md(pd.DataFrame(rows), 1) + "\n")
w("### 6.4 Hubs of the supported graph (number of supported edges touching each protein)\n\n")
for d in DS:
    s = E[d][E[d].data_supported == 1]; deg = pd.concat([s.source, s.target]).value_counts()
    allp = sorted(set(E[d].source) | set(E[d].target)); zero = [p for p in allp if p not in deg.index]
    w(f"- **{d}** top 15: " + ", ".join(f"{k} {v}" for k, v in deg.head(15).items()) + "\n")
    w(f"  - proteins with no supported edge ({len(zero)}): " + ", ".join(zero) + "\n")
w("\n### 6.5 Supported edges, strongest first (top 40 per dataset; full lists in `{ds}_causal_edges.csv`)\n\n")
for d in DS:
    s = E[d][E[d].data_supported == 1].sort_values("W", ascending=False).head(40)
    w(f"**{d}** ({int(E[d].data_supported.sum())} supported)\n\n" + md(s[["source", "target", "prior_C", "W", "beta", "n_pairwise", "q_fdr"]]) + "\n")
w("### 6.6 Selected edges across datasets\n\n")
pairs = [("TMED2", "TMED10"), ("TMED9", "TMED2"), ("TMED9", "TMED10"), ("APP", "BACE1"), ("LRP1", "APP"), ("APOE", "APP"),
         ("CLU", "APOE"), ("SORL1", "APP"), ("PSEN1", "NCSTN"), ("BIN1", "BACE1"), ("PICALM", "APP")]
rows = []
for u, v in pairs:
    for d in DS:
        e = E[d]; x = e[((e.source == u) & (e.target == v)) | ((e.source == v) & (e.target == u))]
        for _, r in x.iterrows():
            rows.append(dict(pair=f"{u}–{v}", dataset=d, edge=f"{r.source}→{r.target}", prior_C=r.prior_C, W=r.W, beta=r.beta,
                             n=int(r.n_pairwise), q=r.q_fdr, supported=int(r.data_supported)))
w(md(pd.DataFrame(rows).drop_duplicates()) + "\n---\n\n")

# 7 comparison
cp = pd.read_csv(R / "comparison_pairs.csv"); cons = pd.read_csv(R / "consensus_edges.csv")
w("## 7. Agreement between the three graphs (HIW 5)\n\n### 7.1 Pairwise overlap of supported pairs vs chance\n\n")
w(md(cp.rename(columns={"direction_agree": "same_direction"})) + "\n")
cons["oneway"] = (cons.prior_C_a_to_b < cr.MIN_PRIOR) | (cons.prior_C_b_to_a < cr.MIN_PRIOR)
two_ = cons[cons.n_datasets == 2]; tri_ = cons[cons.n_datasets == 3]
pc_ = {f"{a_}+{b_}": int(((two_[f"W_{a_}"] > 0) & (two_[f"W_{b_}"] > 0)).sum()) for a_, b_ in [("ROSMAP", "Diverse"), ("Diverse", "Banner"), ("ROSMAP", "Banner")]}
ow_, tw_ = cons[cons.oneway], cons[~cons.oneway]; tow_, ttw_ = tri_[tri_.oneway], tri_[~tri_.oneway]
tcons_ = ttw_[ttw_.direction_consistent == 1]
w(f"""### 7.2 Consensus pairs

- pairs supported in ≥ 2 datasets: **{len(cons)}** (in exactly 2: {len(two_)}; in all 3: **{len(tri_)}**)
- in exactly 2: """ + ", ".join(f"{k} {v}" for k, v in pc_.items()) + f"""
- arrow identical in every dataset that supports the pair: {int(cons.direction_consistent.sum())} of {len(cons)}; among the {len(tri_)} triple pairs: {int(tri_.direction_consistent.sum())} of {len(tri_)}
- split by prior type — one-way C pairs (arrow fixed by C): {len(ow_)} of {len(cons)}, {int(ow_.direction_consistent.sum())} consistent; two-way C pairs (arrow picked by the direction test): {len(tw_)}, of which {int(tw_.direction_consistent.sum())} consistent and {int((tw_.direction_consistent == 0).sum())} not
- among the {len(tri_)} triple pairs: {len(tow_)} one-way ({int(tow_.direction_consistent.sum())} consistent); {len(ttw_)} two-way, of which {len(tcons_)} consistent""" + (" (" + ", ".join(f"{a_}–{b_}" for a_, b_ in zip(tcons_.protein_a, tcons_.protein_b)) + ")" if 0 < len(tcons_) <= 15 else "") + f""" and {int((ttw_.direction_consistent == 0).sum())} not

**All {len(tri_)} pairs supported in all three datasets**

""")
t = cons[cons.n_datasets == 3][["protein_a", "protein_b", "prior_C_a_to_b", "prior_C_b_to_a", "dir_ROSMAP", "W_ROSMAP", "dir_Diverse", "W_Diverse", "dir_Banner", "W_Banner", "direction_consistent"]]
w(md(t, 3) + "\n")
two = cons[cons.n_datasets == 2].copy()
for d in MAIN:
    two[f"W_{d}"] = [("—" if (pd.isna(a) or a == "") else f"{b:.3f}") for a, b in zip(two[f"dir_{d}"], two[f"W_{d}"])]
w(f"**Pairs supported in exactly two datasets** ({len(two)})\n\n" + md(two[["protein_a", "protein_b", "dir_ROSMAP", "W_ROSMAP", "dir_Diverse", "W_Diverse", "dir_Banner", "W_Banner"]], 3) + "\n---\n\n")

# 7b comparison with another prior's run
if A.compare_with:
    R0 = Path(A.compare_with)
    S0 = {d: json.loads((R0 / f"{d}_summary.json").read_text()) for d in DS}
    E0 = {d: pd.read_csv(R0 / f"{d}_causal_edges.csv") for d in DS}
    X0 = json.loads((R0 / "extra_stats.json").read_text())
    w(f"## 7b. This run vs the run in `{R0}/` (previous prior)\n\n")
    rows = []
    for d in DS:
        rows.append(dict(dataset=d, AUROC_raw_prev=S0[d]["auroc_raw"], AUROC_raw_this=S[d]["auroc_raw"],
                         AUROC_adj_prev=S0[d]["auroc_adjusted"], AUROC_adj_this=S[d]["auroc_adjusted"],
                         prior_edges_prev=S0[d]["prior_edges"], prior_edges_this=S[d]["prior_edges"],
                         supported_prev=S0[d]["data_supported_edges"], supported_this=S[d]["data_supported_edges"],
                         STRONG_orig_prev=S0[d]["n_strong_regions"], STRONG_orig_this=S[d]["n_strong_regions"]))
    w(md(pd.DataFrame(rows)) + "\n")
    w("| dataset | Z* previous | Z* this run |\n|---|---|---|\n" + "".join(f"| {d} | {', '.join(S0[d]['Z_star'])} | {', '.join(S[d]['Z_star'])} |\n" for d in DS) + "\n")
    und = lambda df: {tuple(sorted(x)) for x in zip(df.source, df.target)}
    rows = []
    for d in DS:
        a0, a1 = und(E0[d][E0[d].data_supported == 1]), und(E[d][E[d].data_supported == 1])
        p0, p1 = und(E0[d]), und(E[d])
        rows.append(dict(dataset=d, supported_pairs_prev=len(a0), supported_pairs_this=len(a1), supported_in_both=len(a0 & a1),
                         prev_supported_not_prior_now=len(a0 - p1), this_supported_not_prior_before=len(a1 - p0),
                         jaccard=round(len(a0 & a1) / max(1, len(a0 | a1)), 3)))
    w("Supported protein pairs (unordered) in the two runs:\n\n" + md(pd.DataFrame(rows)) + "\n")
    w("`prev_supported_not_prior_now` = pairs supported before whose link is no longer a prior edge in this matrix; "
      "`this_supported_not_prior_before` = supported pairs that were not prior edges before.\n\n")
    if "Edge_List" in _xl.sheet_names and "Provenance" in pd.read_excel(A.prior, sheet_name="Edge_List", nrows=1).columns:
        prov = {(r.Source, r.Target): r.Provenance for r in EL.itertuples()}
        rows = []
        for d in DS:
            s_ = E[d][E[d].data_supported == 1]
            vc = pd.Series([prov.get((a_, b_), "?") for a_, b_ in zip(s_.source, s_.target)]).value_counts()
            r = {"dataset": d, "supported": len(s_)}; r.update({k: int(v) for k, v in vc.items()}); rows.append(r)
        w("Provenance of each run's supported edges (the arrow chosen, looked up in `Edge_List`):\n\n" + md(pd.DataFrame(rows).fillna(0)) + "\n")
        ch = ELe[ELe.Provenance == "CHRONOS"]
        w(f"The matrix has {len(ch)} off-diagonal edges tagged `CHRONOS`, i.e. set using the earlier 3-cohort run on these same datasets.\n\n")
    w("---\n\n")

# 8 banner technical
w("""## 8. Banner-specific checks

### 8.1 TMT batch mapping (HIW 3.1)

Mean of Y-chromosome proteins (RPS4Y1, DDX3Y, EIF1AY, KDM5D, USP9Y, NLGN4Y) per sample, 198 samples with recorded sex:

| | n | mean | min | max |
|---|---|---|---|---|
| female | 86 | −0.541 | −1.817 | −0.072 |
| male | 112 | +0.243 | −0.220 | +0.483 |

- best single cut-off (0.022): 99.0% correctly classified (3 females above the lowest male, 2 males below the highest female)
- with batch labels randomly shuffled (200 times): best accuracy mean 57.8%, maximum 64.1%
- BannerLFQ detects only one Y protein (USP9Y, in 20 samples), so this check is not possible there

### 8.2 Banner TMT vs Banner LFQ (same 190 people)

""")
_rt = (R / "comparison_report.txt").read_text(); _sec6 = _rt[_rt.index("6. TECHNICAL"):] if "6. TECHNICAL" in _rt else ""
w(f"""| | TMT | LFQ |
|---|---|---|
| graph nodes present | {S['Banner']['proteins_found']} | {S['BannerLFQ']['proteins_found']} |
| signed AUROC raw → adjusted | {S['Banner']['auroc_raw']:.4f} → {S['Banner']['auroc_adjusted']:.4f} | {S['BannerLFQ']['auroc_raw']:.4f} → {S['BannerLFQ']['auroc_adjusted']:.4f} |
| Z* | {", ".join(S['Banner']['Z_star'])} | {", ".join(S['BannerLFQ']['Z_star'])} |
| supported edges | {S['Banner']['data_supported_edges']} | {S['BannerLFQ']['data_supported_edges']} |
| STRONG regions (original) | {S['Banner']['n_strong_regions']} | {S['BannerLFQ']['n_strong_regions']} |

From `comparison_report.txt`:

```
{_sec6.strip()}
```

""")
w("""### 8.3 Possible overlap with Diverse

- Diverse contains 43 donors from the Banner cohort; 17 of them are among the 964 Diverse people analysed.
- No ID links the two systems. Matching on sex + APOE + age at death + Braak gave 8 Banner donors
  (09-17, 11-14, 01-16, 06-09, 05-57, 05-35, 13-46, 13-66), matched to 7 Diverse donors; all 8 were removed from Banner and BannerLFQ.

---

## 9. File index (`results/`)

| file | content |
|---|---|
| `{ROSMAP,Diverse,Banner,BannerLFQ}_causal_W.csv` | 80 × 80 weights |
| `…_causal_beta.csv` | 80 × 80 effect sizes |
| `…_causal_N.csv` | 80 × 80 people per pair |
| `…_causal_edges.csv` | 688 prior edges with W, β, n, p, q, data_supported |
| `…_variables.csv` | every variable at every selection step |
| `…_regions.csv` | every region (original scoring) |
| `{ROSMAP,Diverse,Banner}_regions_usable_null_check.csv` | regions re-scored with the corrected null |
| `Diverse_best_region_W.csv` | graph refitted in the original top Diverse region |
| `…_summary.json`, `…_report.txt`, `….log`, `…_individuals.txt` | run summaries, logs, people used |
| `comparison_report.txt`, `comparison_pairs.csv`, `consensus_edges.csv` | cross-dataset comparison |
""")
_txt = "".join(out)
_txt = _txt.replace("| `…_causal_edges.csv` | 688 prior edges with", f"| `…_causal_edges.csv` | {S['ROSMAP']['prior_edges']} prior edges with")
_txt = _txt.replace("## 9. File index (`results/`)", f"## 9. File index (`{R}/`)")
Path(A.out).write_text(_txt)
print("written", sum(len(x) for x in out), "chars")
