"""
CHRONOS - Hypothesis 1: TMED2/TMED10/TMED9 Braak Trajectory
=============================================================
Script 1: Tests A through H
  A - Raw expression trajectory across Braak (descriptive)
  B - Unadjusted continuous Braak regression
  C - Jonckheere-Terpstra monotonic trend test
  D - Covariate-adjusted Braak regression (age, sex, PMI, batch)
  E - + neuronal composition adjustment (marker-gene based)
  F - + full cell-type composition adjustment
  G - Nonlinear Braak trajectory (linear vs quadratic)
  H - Pre-symptomatic contrast (Braak 0-2 vs 3-4 vs 5-6)

Input:
  ROSMAP_DLPFC_logCPM.tsv
  ROSMAP_clinical.csv
  ROSMAP_biospecimen_metadata.csv
  ROSMAP_assay_rnaSeq_metadata.csv (for batch)

Output:
  tmed_testA_trajectory.png
  tmed_testB_regression.csv
  tmed_testC_jt.csv
  tmed_testD_adjusted.csv
  tmed_testEF_composition.csv
  tmed_testG_nonlinear.csv
  tmed_testH_presymptomatic.csv
  tmed_script1_summary.csv
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

# ── GENE IDs ───────────────────────────────────────────────────────────────────

TMED_GENES = {
    "TMED2":  "ENSG00000086598",
    "TMED10": "ENSG00000170348",
    "TMED9":  "ENSG00000184840",
}

# Marker genes for cell-type deconvolution (well-established markers)
CELL_MARKERS = {
    "neuron":     ["ENSG00000132639",  # SNAP25
                   "ENSG00000197971",  # MBP — wait, that's oligo; use SYP
                   "ENSG00000102003",  # SYP
                   "ENSG00000067715",  # SYT1
                   "ENSG00000108556"], # CHRM1 — skip, use RBFOX3
    "neuron_clean": ["ENSG00000102003",  # SYP
                     "ENSG00000067715",  # SYT1
                     "ENSG00000132639",  # SNAP25
                     "ENSG00000008056",  # SYN1
                     "ENSG00000157542"], # SYN2
    "astrocyte":  ["ENSG00000131095",  # GFAP
                   "ENSG00000026025",  # VIM
                   "ENSG00000116016",  # EPAS1 — skip; use AQP4
                   "ENSG00000171885"], # AQP4
    "astrocyte_clean": ["ENSG00000131095",  # GFAP
                        "ENSG00000171885",  # AQP4
                        "ENSG00000135926"], # TTYH1
    "microglia":  ["ENSG00000197249",  # AIF1 (IBA1)
                   "ENSG00000101439",  # CSF1R
                   "ENSG00000168685"], # IL7R — skip; use TMEM119
    "microglia_clean": ["ENSG00000197249",  # AIF1
                        "ENSG00000101439"], # CSF1R
    "oligodendrocyte": ["ENSG00000197971",  # MBP
                        "ENSG00000123560",  # PLP1
                        "ENSG00000197891"], # MOG — check
    "oligodendrocyte_clean": ["ENSG00000197971",  # MBP
                               "ENSG00000123560"], # PLP1
}

COLORS = {"TMED2": "steelblue", "TMED10": "seagreen", "TMED9": "firebrick"}

# ── 1. LOAD DATA ───────────────────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

# fix age
clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+", "90")
clinical["age_at_visit_max"] = pd.to_numeric(clinical["age_at_visit_max"], errors="coerce")
for col in ["braaksc","msex","pmi","cogdx","ceradsc"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

print(f"  Expr matrix: {expr.shape}")

# ── 2. LINK SAMPLES ────────────────────────────────────────────────────────────

print("Linking samples to metadata...")
sample_ids = expr.columns.tolist()
rna_bio = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()

clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi","cogdx","ceradsc","apoe_genotype"]
meta = rna_bio.merge(clinical[clin_cols], on="individualID", how="left")

# add batch from rna_meta
rna_meta_sub = rna_meta[["specimenID","rnaBatch","sequencingBatch"]].drop_duplicates()
meta = meta.merge(rna_meta_sub, on="specimenID", how="left")
meta = meta.set_index("specimenID")

print(f"  Samples linked: {len(meta)}")
print(f"  Braak available: {meta['braaksc'].notna().sum()}")
print(f"  Braak distribution:\n{meta['braaksc'].value_counts().sort_index()}")

# ── 3. EXTRACT TMED GENES ──────────────────────────────────────────────────────

tmed_df = pd.DataFrame({
    name: expr.loc[eid] for name, eid in TMED_GENES.items() if eid in expr.index
})
print(f"\n  TMED genes found: {list(tmed_df.columns)}")

# merge with metadata
df = tmed_df.join(meta)
df_braak = df.dropna(subset=["braaksc"]).copy()
df_braak["braaksc"] = df_braak["braaksc"].astype(int)
print(f"  Samples with Braak: {len(df_braak)}")

# ── 4. CELL TYPE SCORES (marker-gene mean) ─────────────────────────────────────

print("\nEstimating cell-type composition from marker genes...")
ct_scores = {}
ct_map = {
    "neuron_score":    CELL_MARKERS["neuron_clean"],
    "astrocyte_score": CELL_MARKERS["astrocyte_clean"],
    "microglia_score": CELL_MARKERS["microglia_clean"],
    "oligo_score":     CELL_MARKERS["oligodendrocyte_clean"],
}
for ct, markers in ct_map.items():
    present = [m for m in markers if m in expr.index]
    if present:
        vals = expr.loc[present]
        z = (vals - vals.mean(axis=1).values[:,None]) / (vals.std(axis=1).values[:,None] + 1e-10)
        ct_scores[ct] = z.mean(axis=0)
        print(f"  {ct}: {len(present)}/{len(markers)} markers found")

ct_df = pd.DataFrame(ct_scores)
df_braak = df_braak.join(ct_df, how="left")

# ── 5. OLS HELPER ──────────────────────────────────────────────────────────────

def ols_result(y, X_df, term):
    """OLS regression, returns stats for one term."""
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y.values
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs, _, _, _ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm @ coeffs
    n, p = len(y_v), Xm.shape[1]
    mse = np.sum(resid**2) / (n - p)
    try:
        cov = mse * np.linalg.inv(Xm.T @ Xm)
    except np.linalg.LinAlgError:
        return {"beta": np.nan, "se": np.nan, "ci_lo": np.nan,
                "ci_hi": np.nan, "p": np.nan, "r2": np.nan, "n": n}
    se = np.sqrt(np.diag(cov))
    cols = ["intercept"] + list(X_df.columns)
    idx = cols.index(term)
    beta, se_b = coeffs[idx], se[idx]
    t = beta / se_b
    pv = 2 * stats.t.sf(abs(t), df=n-p)
    # partial R²
    ss_res_full = np.sum(resid**2)
    Xr = np.delete(Xm, idx, axis=1)
    cr, _, _, _ = np.linalg.lstsq(Xr, y_v, rcond=None)
    ss_res_red = np.sum((y_v - Xr @ cr)**2)
    partial_r2 = (ss_res_red - ss_res_full) / ss_res_red
    return {"beta": beta, "se": se_b,
            "ci_lo": beta - 1.96*se_b, "ci_hi": beta + 1.96*se_b,
            "p": pv, "r2": partial_r2, "n": int(mask.sum())}

# ── TEST A: DESCRIPTIVE TRAJECTORY ────────────────────────────────────────────

print("\n=== TEST A: Descriptive trajectory ===")
braak_stats = []
for braak in sorted(df_braak["braaksc"].unique()):
    grp = df_braak[df_braak["braaksc"] == braak]
    row = {"braak": braak, "n": len(grp)}
    for gene in TMED_GENES:
        v = grp[gene].dropna()
        row[f"{gene}_mean"] = v.mean()
        row[f"{gene}_sd"]   = v.std()
        row[f"{gene}_se"]   = v.sem()
        row[f"{gene}_med"]  = v.median()
        row[f"{gene}_q25"]  = v.quantile(0.25)
        row[f"{gene}_q75"]  = v.quantile(0.75)
    braak_stats.append(row)
    print(f"  Braak {braak}: n={len(grp)}, "
          + ", ".join(f"{g}={row[f'{g}_mean']:.3f}" for g in TMED_GENES))

braak_stats_df = pd.DataFrame(braak_stats)

# ── TEST B: UNADJUSTED REGRESSION ─────────────────────────────────────────────

print("\n=== TEST B: Unadjusted regression ===")
testB = []
for gene in TMED_GENES:
    X = df_braak[["braaksc"]]
    res = ols_result(df_braak[gene], X, "braaksc")
    print(f"  {gene}: β={res['beta']:.4f} [{res['ci_lo']:.4f},{res['ci_hi']:.4f}], "
          f"p={res['p']:.3e}, R²={res['r2']:.4f}, n={res['n']}")
    testB.append({"test":"B_unadjusted","gene":gene,**res})

# ── TEST C: JONCKHEERE-TERPSTRA ────────────────────────────────────────────────

print("\n=== TEST C: Jonckheere-Terpstra monotonic trend ===")

def jonckheere_terpstra(groups_data):
    """JT test statistic and p-value (normal approximation)."""
    k = len(groups_data)
    U = 0
    for i in range(k-1):
        for j in range(i+1, k):
            xi, xj = groups_data[i], groups_data[j]
            for a in xi:
                for b in xj:
                    if b > a: U += 1
                    elif b == a: U += 0.5
    # expected value and variance under H0
    N = sum(len(g) for g in groups_data)
    ns = [len(g) for g in groups_data]
    E = (N**2 - sum(n**2 for n in ns)) / 4
    n2 = sum(n**2 for n in ns)
    n3 = sum(n**3 for n in ns)
    V = (N**2*(2*N+3) - sum(n**2*(2*n+3) for n in ns)) / 72
    z = (U - E) / np.sqrt(V)
    p = 2 * stats.norm.sf(abs(z))
    return z, p

testC = []
for gene in TMED_GENES:
    groups = [df_braak[df_braak["braaksc"]==b][gene].dropna().values
              for b in sorted(df_braak["braaksc"].unique())]
    z, p = jonckheere_terpstra(groups)
    print(f"  {gene}: JT z={z:.3f}, p={p:.3e}")
    testC.append({"test":"C_JT","gene":gene,"JT_z":z,"p":p})

# ── TEST D: COVARIATE-ADJUSTED ─────────────────────────────────────────────────

print("\n=== TEST D: Covariate-adjusted (age, sex, PMI, batch) ===")
df_D = df_braak.dropna(subset=["age_at_visit_max","msex","pmi"]).copy()

# encode rnaBatch as dummies if available
batch_cols = []
if "rnaBatch" in df_D.columns and df_D["rnaBatch"].notna().sum() > 10:
    batch_dummies = pd.get_dummies(df_D["rnaBatch"], prefix="batch", drop_first=True)
    df_D = pd.concat([df_D, batch_dummies], axis=1)
    batch_cols = list(batch_dummies.columns)

cov_D = ["age_at_visit_max","msex","pmi"] + batch_cols
testD = []
for gene in TMED_GENES:
    X = df_D[["braaksc"] + cov_D]
    res = ols_result(df_D[gene], X, "braaksc")
    print(f"  {gene}: β={res['beta']:.4f} [{res['ci_lo']:.4f},{res['ci_hi']:.4f}], "
          f"p={res['p']:.3e}, R²={res['r2']:.4f}, n={res['n']}")
    testD.append({"test":"D_adjusted","gene":gene,**res})

# ── TEST E: + NEURON PROPORTION ────────────────────────────────────────────────

print("\n=== TEST E: + neuronal composition ===")
testE = []
if "neuron_score" in df_D.columns:
    for gene in TMED_GENES:
        X = df_D[["braaksc","neuron_score"] + cov_D]
        res = ols_result(df_D[gene], X, "braaksc")
        print(f"  {gene}: β={res['beta']:.4f}, p={res['p']:.3e} (after neuron adjustment)")
        testE.append({"test":"E_neuron_adj","gene":gene,**res})
else:
    print("  neuron_score not available")

# ── TEST F: FULL CELL COMPOSITION ──────────────────────────────────────────────

print("\n=== TEST F: Full cell-type composition ===")
testF = []
ct_cols = [c for c in ["neuron_score","astrocyte_score","microglia_score","oligo_score"]
           if c in df_D.columns]
if ct_cols:
    for gene in TMED_GENES:
        X = df_D[["braaksc"] + ct_cols + cov_D]
        res = ols_result(df_D[gene], X, "braaksc")
        print(f"  {gene}: β={res['beta']:.4f}, p={res['p']:.3e} (full cell-type adj)")
        testF.append({"test":"F_fullcell_adj","gene":gene,**res})

# ── TEST G: NONLINEAR (LINEAR vs QUADRATIC) ────────────────────────────────────

print("\n=== TEST G: Nonlinear trajectory ===")
testG = []
df_G = df_D.copy()
df_G["braaksc2"] = df_G["braaksc"] ** 2

for gene in TMED_GENES:
    # linear
    Xl = df_G[["braaksc"] + cov_D]
    resl = ols_result(df_G[gene], Xl, "braaksc")
    # quadratic
    Xq = df_G[["braaksc","braaksc2"] + cov_D]
    resq = ols_result(df_G[gene], Xq, "braaksc2")

    # compare R² improvement
    mask = df_G[["braaksc","braaksc2"] + cov_D + [gene]].dropna().index
    y_v = df_G.loc[mask, gene].values
    Xl_m = np.column_stack([np.ones(len(mask))] +
                            [df_G.loc[mask,c].values for c in ["braaksc"]+cov_D])
    Xq_m = np.column_stack([np.ones(len(mask))] +
                            [df_G.loc[mask,c].values for c in ["braaksc","braaksc2"]+cov_D])
    cl, _, _, _ = np.linalg.lstsq(Xl_m, y_v, rcond=None)
    cq, _, _, _ = np.linalg.lstsq(Xq_m, y_v, rcond=None)
    ss_lin = np.sum((y_v - Xl_m@cl)**2)
    ss_qua = np.sum((y_v - Xq_m@cq)**2)
    # F-test for quadratic term
    n, p = len(y_v), Xq_m.shape[1]
    F = ((ss_lin - ss_qua)/1) / (ss_qua/(n-p))
    pF = stats.f.sf(F, 1, n-p)

    print(f"  {gene}: linear β={resl['beta']:.4f}  | quadratic term β={resq['beta']:.4f}, p={resq['p']:.3e} | F-test p={pF:.3e}")
    testG.append({"gene":gene,
                  "beta_linear":resl["beta"], "p_linear":resl["p"],
                  "beta_quad":resq["beta"], "p_quad":resq["p"],
                  "F_quad":F, "p_F":pF})

# ── TEST H: PRE-SYMPTOMATIC CONTRAST ──────────────────────────────────────────

print("\n=== TEST H: Pre-symptomatic contrast (0-2 vs 3-4 vs 5-6) ===")
def braak_group(b):
    if b <= 2: return "early (0-2)"
    elif b <= 4: return "mid (3-4)"
    else: return "late (5-6)"

df_braak["braak_group"] = df_braak["braaksc"].apply(braak_group)
testH = []
for gene in TMED_GENES:
    early = df_braak[df_braak["braak_group"]=="early (0-2)"][gene].dropna()
    mid   = df_braak[df_braak["braak_group"]=="mid (3-4)"][gene].dropna()
    late  = df_braak[df_braak["braak_group"]=="late (5-6)"][gene].dropna()
    # ANOVA
    F, pF = stats.f_oneway(early, mid, late)
    # pairwise t-tests
    t_em, p_em = stats.ttest_ind(early, mid)
    t_ml, p_ml = stats.ttest_ind(mid, late)
    t_el, p_el = stats.ttest_ind(early, late)
    print(f"  {gene}: early={early.mean():.3f} mid={mid.mean():.3f} late={late.mean():.3f} | ANOVA p={pF:.3e}")
    print(f"    early→mid: p={p_em:.3e}  mid→late: p={p_ml:.3e}  early→late: p={p_el:.3e}")
    testH.append({"gene":gene,
                  "mean_early":early.mean(),"mean_mid":mid.mean(),"mean_late":late.mean(),
                  "n_early":len(early),"n_mid":len(mid),"n_late":len(late),
                  "ANOVA_F":F,"ANOVA_p":pF,
                  "p_early_mid":p_em,"p_mid_late":p_ml,"p_early_late":p_el})

# ── FDR CORRECTION ON PRIMARY TESTS (B, D) ────────────────────────────────────

print("\n=== FDR correction (BH) on primary Braak tests ===")
for label, tests in [("Test B (unadjusted)", testB), ("Test D (adjusted)", testD)]:
    ps = [t["p"] for t in tests]
    _, fdrs, _, _ = multipletests(ps, method="fdr_bh")
    for t, fdr in zip(tests, fdrs):
        t["fdr"] = fdr
    print(f"  {label}: " + "  ".join(f"{t['gene']} fdr={fdr:.3f}" for t,fdr in zip(tests,fdrs)))

# ── SAVE RESULTS ──────────────────────────────────────────────────────────────

all_results = testB + testC + testD + testE + testF
summary_df = pd.DataFrame(all_results)
summary_df.to_csv("tmed_script1_summary.csv", index=False)
pd.DataFrame(testG).to_csv("tmed_testG_nonlinear.csv", index=False)
pd.DataFrame(testH).to_csv("tmed_testH_presymptomatic.csv", index=False)
braak_stats_df.to_csv("tmed_testA_stats.csv", index=False)
print("\nSaved CSVs.")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

# Plot A: trajectory with SE bands
ax = fig.add_subplot(gs[0, :2])
for gene in TMED_GENES:
    means = braak_stats_df[f"{gene}_mean"]
    ses   = braak_stats_df[f"{gene}_se"]
    braaks = braak_stats_df["braak"]
    ax.plot(braaks, means, "o-", color=COLORS[gene], label=gene, linewidth=2, markersize=6)
    ax.fill_between(braaks, means-ses, means+ses, alpha=0.2, color=COLORS[gene])
ax.set_xlabel("Braak stage", fontsize=11)
ax.set_ylabel("Mean logCPM ± SE", fontsize=11)
ax.set_title("Test A: Raw expression trajectory across Braak", fontsize=11)
ax.legend(fontsize=10)
ax.set_xticks(range(7))

# Plot A2: sample sizes per Braak
ax2 = fig.add_subplot(gs[0, 2])
ax2.bar(braak_stats_df["braak"], braak_stats_df["n"], color="grey", alpha=0.7)
ax2.set_xlabel("Braak stage"); ax2.set_ylabel("n samples")
ax2.set_title("Sample sizes per Braak", fontsize=10)
ax2.set_xticks(range(7))

# Plot B+D: beta comparison (unadjusted vs adjusted)
ax3 = fig.add_subplot(gs[1, :2])
genes = list(TMED_GENES.keys())
x = np.arange(len(genes)); w = 0.35
bB = [t["beta"] for t in testB]
bD = [t["beta"] for t in testD]
eB = [1.96*t["se"] for t in testB]
eD = [1.96*t["se"] for t in testD]
bars1 = ax3.bar(x-w/2, bB, w, yerr=eB, label="Unadjusted", color="steelblue", alpha=0.8, capsize=4)
bars2 = ax3.bar(x+w/2, bD, w, yerr=eD, label="Adjusted", color="firebrick", alpha=0.8, capsize=4)
ax3.axhline(0, color="black", linewidth=0.8)
ax3.set_xticks(x); ax3.set_xticklabels(genes, fontsize=11)
ax3.set_ylabel("β per Braak unit (95% CI)", fontsize=10)
ax3.set_title("Tests B & D: Braak regression coefficients", fontsize=11)
ax3.legend(fontsize=9)

# mark significance
for i, (tB, tD) in enumerate(zip(testB, testD)):
    for xi, t, offset in [(i-w/2, tB, 0), (i+w/2, tD, 0)]:
        sig = "***" if t["p"]<0.001 else "**" if t["p"]<0.01 else "*" if t["p"]<0.05 else ""
        if sig:
            ypos = t["beta"] + 1.96*t["se"] + 0.005
            ax3.text(xi, ypos, sig, ha="center", fontsize=10)

# Plot E+F: Braak beta after cell-type adjustment
ax4 = fig.add_subplot(gs[1, 2])
if testE and testF:
    bE = [t["beta"] for t in testE]
    bF = [t["beta"] for t in testF]
    ax4.bar(x-w/2, bE, w, label="+neuron", color="seagreen", alpha=0.8)
    ax4.bar(x+w/2, bF, w, label="+all cells", color="orange", alpha=0.8)
    ax4.axhline(0, color="black", linewidth=0.8)
    ax4.set_xticks(x); ax4.set_xticklabels(genes, fontsize=9)
    ax4.set_ylabel("β (Braak)", fontsize=10)
    ax4.set_title("Tests E&F: Cell-type adjustment", fontsize=10)
    ax4.legend(fontsize=8)

# Plot G: linear vs quadratic fits for each gene
for i, gene in enumerate(TMED_GENES):
    ax = fig.add_subplot(gs[2, i])
    x_v = df_braak["braaksc"].values
    y_v = df_braak[gene].values
    mask = ~(np.isnan(x_v)|np.isnan(y_v))
    x_v, y_v = x_v[mask], y_v[mask]
    ax.scatter(x_v + np.random.normal(0, 0.1, len(x_v)), y_v,
               s=4, alpha=0.3, color=COLORS[gene])
    xl = np.linspace(0, 6, 100)
    # linear fit
    cl = np.polyfit(x_v, y_v, 1)
    ax.plot(xl, np.polyval(cl, xl), color="black", linewidth=1.5, label="linear")
    # quadratic fit
    cq = np.polyfit(x_v, y_v, 2)
    ax.plot(xl, np.polyval(cq, xl), color="firebrick", linewidth=1.5,
            linestyle="--", label="quadratic")
    gG = next((t for t in testG if t["gene"]==gene), None)
    pF_str = f"F-test p={gG['p_F']:.2e}" if gG else ""
    ax.set_title(f"{gene}\n{pF_str}", fontsize=10)
    ax.set_xlabel("Braak"); ax.set_ylabel("logCPM", fontsize=9)
    ax.legend(fontsize=7)

plt.suptitle("CHRONOS — TMED2/TMED10/TMED9 Braak Trajectory (Tests A–H)", fontsize=13, y=1.01)
plt.savefig("tmed_trajectory_script1.png", dpi=200, bbox_inches="tight")
print("Saved: tmed_trajectory_script1.png")

# ── PLOT H: Pre-symptomatic ────────────────────────────────────────────────────

fig2, axes = plt.subplots(1, 3, figsize=(13, 5))
fig2.suptitle("Test H: Pre-symptomatic contrast", fontsize=12)
groups_order = ["early (0-2)", "mid (3-4)", "late (5-6)"]
for ax, gene in zip(axes, TMED_GENES):
    data = [df_braak[df_braak["braak_group"]==g][gene].dropna().values
            for g in groups_order]
    ax.boxplot(data, labels=["0-2", "3-4", "5-6"], patch_artist=True,
               boxprops=dict(facecolor=COLORS[gene], alpha=0.6))
    ax.set_xlabel("Braak group"); ax.set_ylabel("logCPM")
    ax.set_title(gene, fontsize=11)

plt.tight_layout()
plt.savefig("tmed_testH_presymptomatic.png", dpi=200, bbox_inches="tight")
print("Saved: tmed_testH_presymptomatic.png")
print("\nAll done.")
