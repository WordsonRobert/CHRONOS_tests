"""
Hypothesis 9: Cell-type composition analysis
Tests A-K: composition shift, nested TMED models, neuron-marker controls, mediation decomposition
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ── paths ──────────────────────────────────────────────────────────────────────
DATA = "/home/wordson22/projects/CHRONOS/DATA/"
OUT  = "./"

# ── Ensembl IDs ────────────────────────────────────────────────────────────────
TMED2_ID  = "ENSG00000086598"
TMED10_ID = "ENSG00000170348"
TMED9_ID  = "ENSG00000184840"

# neuronal markers (Test E/K)
NEURON_CTRL = {
    "RBFOX3":  "ENSG00000167281",
    "SYT1":    "ENSG00000067715",
    "SNAP25":  "ENSG00000132639",
    "NEFL":    "ENSG00000277586",
    "TUBB3":   "ENSG00000258947",
    "MAP2":    "ENSG00000078018",
    "SYN1":    "ENSG00000008056",
    "CAMK2A":  "ENSG00000070808",
    "GRIN1":   "ENSG00000176884",
    "SLC17A7": "ENSG00000104888",
}

# cell-type markers for composition scores
NEURON_M = ["ENSG00000102003","ENSG00000067715","ENSG00000132639",
            "ENSG00000008056","ENSG00000157542"]
ASTRO_M  = ["ENSG00000131095","ENSG00000171885"]
MICRO_M  = ["ENSG00000204472","ENSG00000138185"]
OLIGO_M  = ["ENSG00000197971","ENSG00000123560"]

# non-neuronal markers (Test F)
ASTRO_CTRL = {"GFAP":"ENSG00000131095","AQP4":"ENSG00000171885",
              "S100B":"ENSG00000160307","VIM":"ENSG00000026025"}
MICRO_CTRL = {"AIF1":"ENSG00000204472","CSF1R":"ENSG00000138185",
              "TMEM119":"ENSG00000183160","P2RY12":"ENSG00000169313"}

# ── load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = pd.to_numeric(
    clinical["age_at_visit_max"].replace("90+","90"), errors="coerce")
for col in ["braaksc","msex","pmi","cogdx"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

clin_cols = ["individualID","braaksc","cogdx","msex","pmi","age_at_visit_max"]
sample_ids = expr.columns.tolist()
rna_bio    = biospec[biospec["specimenID"].isin(sample_ids)][
    ["specimenID","individualID"]].drop_duplicates()
rna_batch  = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
meta = (rna_bio
        .merge(clinical[clin_cols], on="individualID", how="left")
        .merge(rna_batch, on="specimenID", how="left")
        .set_index("specimenID"))

# ── cell type scores ───────────────────────────────────────────────────────────
def ct_score(markers, name):
    present = [m for m in markers if m in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns, name=name)
    vals = expr.loc[present]
    z = (vals - vals.mean(axis=1).values[:,None])/(vals.std(axis=1).values[:,None]+1e-10)
    s = z.mean(axis=0); s.name = name; return s

meta["neuron_score"] = ct_score(NEURON_M, "neuron_score")
meta["astro_score"]  = ct_score(ASTRO_M,  "astro_score")
meta["micro_score"]  = ct_score(MICRO_M,  "micro_score")
meta["oligo_score"]  = ct_score(OLIGO_M,  "oligo_score")

# batch dummies
batch_dum  = pd.get_dummies(meta["rnaBatch"], prefix="batch", drop_first=True)
meta       = pd.concat([meta, batch_dum], axis=1)
batch_cols = list(batch_dum.columns)

# extract TMED genes
for eid, name in [(TMED2_ID,"TMED2"),(TMED10_ID,"TMED10"),(TMED9_ID,"TMED9")]:
    meta[name] = expr.loc[eid] if eid in expr.index else np.nan

# AD binary
meta["AD"] = np.where(meta["cogdx"]==1, 0, np.where(meta["cogdx"]>=4, 1, np.nan))

df = meta.dropna(subset=["braaksc","TMED2","TMED10","TMED9",
                           "neuron_score","astro_score","micro_score","oligo_score"])
print(f"Working n = {len(df)}")

# ── OLS helper ─────────────────────────────────────────────────────────────────
def ols_term(y_s, X_df, term):
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1)|np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm@coeffs
    n, p  = len(y_v), Xm.shape[1]
    mse   = np.sum(resid**2)/(n-p)
    cov   = mse*np.linalg.inv(Xm.T@Xm)
    se    = np.sqrt(np.diag(cov))
    cols  = ["intercept"]+list(X_df.columns)
    idx   = cols.index(term)
    beta, se_b = coeffs[idx], se[idx]
    t  = beta/se_b
    pv = 2*stats.t.sf(abs(t), df=n-p)
    ss_res = np.sum(resid**2); ss_tot = np.sum((y_v-y_v.mean())**2)
    r2 = 1 - ss_res/ss_tot
    return dict(beta=beta, se=se_b, ci_lo=beta-1.96*se_b, ci_hi=beta+1.96*se_b,
                p=pv, n=int(mask.sum()), r2=r2)

def get_resid(y_s, X_df):
    """Residuals of y regressed on X."""
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1)|np.isnan(y_v))
    idx_valid = np.where(mask)[0]
    Xm2, y2 = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm2, y2, rcond=None)
    resid = np.full(len(y_v), np.nan)
    resid[idx_valid] = y2 - Xm2@coeffs
    return pd.Series(resid, index=y_s.index)

cov_base  = ["msex","pmi","age_at_visit_max"] + batch_cols
cov_cells = ["neuron_score","astro_score","micro_score","oligo_score"]

# ══════════════════════════════════════════════════════════════════════════════
# TEST A: Composition shift across Braak and AD
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST A: Composition shift ──")
testA_rows = []
for ct in ["neuron_score","astro_score","micro_score","oligo_score"]:
    sub = df[[ct,"braaksc","AD"]+cov_base].dropna()
    # Braak model
    r1 = ols_term(sub[ct], sub[["braaksc"]+cov_base], "braaksc")
    # AD model (subset with AD defined)
    sub_ad = sub.dropna(subset=["AD"])
    r2 = ols_term(sub_ad[ct], sub_ad[["AD"]+cov_base], "AD") if len(sub_ad)>30 else {}
    testA_rows.append(dict(cell_type=ct,
                           beta_braak=r1["beta"], p_braak=r1["p"], n_braak=r1["n"],
                           beta_AD=r2.get("beta",np.nan), p_AD=r2.get("p",np.nan)))
    print(f"  {ct}: β_Braak={r1['beta']:.4f} p={r1['p']:.4f}  |  "
          f"β_AD={r2.get('beta',np.nan):.4f} p={r2.get('p',np.nan):.4f}")

testA_df = pd.DataFrame(testA_rows)
# FDR
ps = testA_df["p_braak"].values
_, fdr, _, _ = multipletests(ps, method="fdr_bh")
testA_df["fdr_braak"] = fdr

# ══════════════════════════════════════════════════════════════════════════════
# TEST B: Nested TMED models (M1–M4)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST B: Nested TMED models ──")
testB_rows = []
for gene in ["TMED2","TMED10","TMED9"]:
    sub = df[[gene,"braaksc"]+cov_base+cov_cells].dropna()
    # M1: raw Braak only
    r1 = ols_term(sub[gene], sub[["braaksc"]], "braaksc")
    # M2: + demographics/batch
    r2 = ols_term(sub[gene], sub[["braaksc"]+cov_base], "braaksc")
    # M3: + neuron only
    r3 = ols_term(sub[gene], sub[["braaksc"]+cov_base+["neuron_score"]], "braaksc")
    # M4: + all cells
    r4 = ols_term(sub[gene], sub[["braaksc"]+cov_base+cov_cells], "braaksc")

    pct_change_neuron = 100*(r3["beta"]-r1["beta"])/r1["beta"] if r1["beta"]!=0 else np.nan
    pct_change_all    = 100*(r4["beta"]-r1["beta"])/r1["beta"] if r1["beta"]!=0 else np.nan

    testB_rows.append(dict(
        gene=gene,
        M1_beta=r1["beta"], M1_p=r1["p"],
        M2_beta=r2["beta"], M2_p=r2["p"],
        M3_beta=r3["beta"], M3_p=r3["p"],
        M4_beta=r4["beta"], M4_p=r4["p"],
        pct_change_neuron=pct_change_neuron,
        pct_change_all=pct_change_all,
    ))
    print(f"  {gene}:")
    print(f"    M1 raw:         β={r1['beta']:.4f}  p={r1['p']:.4f}")
    print(f"    M2 +covariates: β={r2['beta']:.4f}  p={r2['p']:.4f}")
    print(f"    M3 +neuron:     β={r3['beta']:.4f}  p={r3['p']:.4f}  ({pct_change_neuron:.1f}% change)")
    print(f"    M4 +all cells:  β={r4['beta']:.4f}  p={r4['p']:.4f}  ({pct_change_all:.1f}% change)")

testB_df = pd.DataFrame(testB_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST C: Correlation matrix — TMED genes vs cell fractions
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST C: TMED ↔ cell-fraction correlations ──")
corr_rows = []
for gene in ["TMED2","TMED10","TMED9"]:
    for ct in cov_cells:
        sub = df[[gene,ct]].dropna()
        r_p, p_p = stats.pearsonr(sub[gene], sub[ct])
        r_s, p_s = stats.spearmanr(sub[gene], sub[ct])
        corr_rows.append(dict(gene=gene, cell=ct,
                              pearson_r=r_p, pearson_p=p_p,
                              spearman_r=r_s, spearman_p=p_s))
        print(f"  {gene} ~ {ct}: r={r_p:.3f} p={p_p:.4f}  ρ={r_s:.3f} p={p_s:.4f}")

corr_df = pd.DataFrame(corr_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST D: Residualized TMED expression vs Braak
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST D: Residualized TMED vs Braak ──")
resid_data = {}
for gene in ["TMED2","TMED10","TMED9"]:
    sub = df[[gene,"braaksc"]+cov_base+cov_cells].dropna()
    # regress out everything except braaksc
    X_no_braak = sub[cov_base+cov_cells]
    resid = get_resid(sub[gene], X_no_braak)
    resid_data[gene] = pd.DataFrame({"resid":resid.values, "braaksc":sub["braaksc"].values})
    r, p = stats.pearsonr(resid_data[gene]["braaksc"], resid_data[gene]["resid"])
    print(f"  {gene} residual ~ Braak: r={r:.3f}  p={p:.4f}  n={len(sub)}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST E: Neuronal marker control panel
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST E: Neuronal marker control panel ──")
testE_rows = []
for gname, eid in NEURON_CTRL.items():
    if eid not in expr.index:
        print(f"  {gname} not found, skipping"); continue
    df[gname] = expr.loc[eid]
    sub = df[[gname,"braaksc"]+cov_base].dropna()
    r_no  = ols_term(sub[gname], sub[["braaksc"]+cov_base], "braaksc")
    sub2  = df[[gname,"braaksc"]+cov_base+cov_cells].dropna()
    r_yes = ols_term(sub2[gname], sub2[["braaksc"]+cov_base+cov_cells], "braaksc")
    pct   = 100*(r_yes["beta"]-r_no["beta"])/r_no["beta"] if r_no["beta"]!=0 else np.nan
    testE_rows.append(dict(gene=gname,
                           beta_no_cells=r_no["beta"], p_no_cells=r_no["p"],
                           beta_with_cells=r_yes["beta"], p_with_cells=r_yes["p"],
                           pct_change=pct))
    print(f"  {gname}: β_no_cells={r_no['beta']:.4f} p={r_no['p']:.4f}  "
          f"β_with_cells={r_yes['beta']:.4f} p={r_yes['p']:.4f}  ({pct:.1f}%)")

testE_df = pd.DataFrame(testE_rows)
if len(testE_df) > 0:
    _, fdr_e, _, _ = multipletests(testE_df["p_no_cells"].values, method="fdr_bh")
    testE_df["fdr_no_cells"] = fdr_e

# ══════════════════════════════════════════════════════════════════════════════
# TEST F: Non-neuronal marker controls
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST F: Non-neuronal marker controls ──")
testF_rows = []
for gname, eid in {**ASTRO_CTRL, **MICRO_CTRL}.items():
    if eid not in expr.index:
        continue
    df[gname] = expr.loc[eid]
    sub  = df[[gname,"braaksc"]+cov_base].dropna()
    r_no = ols_term(sub[gname], sub[["braaksc"]+cov_base], "braaksc")
    sub2 = df[[gname,"braaksc"]+cov_base+cov_cells].dropna()
    r_yes= ols_term(sub2[gname], sub2[["braaksc"]+cov_base+cov_cells], "braaksc")
    pct  = 100*(r_yes["beta"]-r_no["beta"])/r_no["beta"] if r_no["beta"]!=0 else np.nan
    cell_class = "Astrocyte" if gname in ASTRO_CTRL else "Microglia"
    testF_rows.append(dict(gene=gname, cell_class=cell_class,
                           beta_no_cells=r_no["beta"], p_no_cells=r_no["p"],
                           beta_with_cells=r_yes["beta"], p_with_cells=r_yes["p"],
                           pct_change=pct))
    print(f"  [{cell_class}] {gname}: β={r_no['beta']:.4f} → {r_yes['beta']:.4f}  ({pct:.1f}%)")

testF_df = pd.DataFrame(testF_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST G: Composition decomposition (mediation-style)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST G: Composition decomposition ──")
testG_rows = []
for gene in ["TMED2","TMED10","TMED9"]:
    sub = df[[gene,"braaksc","neuron_score"]+cov_base].dropna()
    # Path A: Braak → neuron_score
    rA = ols_term(sub["neuron_score"], sub[["braaksc"]+cov_base], "braaksc")
    # Path B: neuron_score → gene (adjusted for Braak)
    rB = ols_term(sub[gene], sub[["neuron_score","braaksc"]+cov_base], "neuron_score")
    # Total effect: gene ~ Braak + covariates
    rT = ols_term(sub[gene], sub[["braaksc"]+cov_base], "braaksc")
    # Direct effect: gene ~ Braak + neuron + covariates
    rD = ols_term(sub[gene], sub[["braaksc","neuron_score"]+cov_base], "braaksc")
    # Indirect effect ≈ path A × path B (product-of-coefficients)
    indirect = rA["beta"] * rB["beta"]
    pct_mediated = 100*indirect/rT["beta"] if rT["beta"]!=0 else np.nan
    testG_rows.append(dict(gene=gene,
                           path_A=rA["beta"], path_B=rB["beta"],
                           total=rT["beta"], direct=rD["beta"],
                           indirect=indirect, pct_mediated=pct_mediated))
    print(f"  {gene}: total={rT['beta']:.4f}  direct={rD['beta']:.4f}  "
          f"indirect≈{indirect:.4f}  ({pct_mediated:.1f}% via neuron)")

testG_df = pd.DataFrame(testG_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST I: Pre-symptomatic subset (Braak ≤ 2)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST I: Pre-symptomatic (Braak ≤ 2) ──")
df_early = df[df["braaksc"] <= 2].copy()
print(f"  Braak ≤ 2: n={len(df_early)}")
testI_rows = []
for gene in ["TMED2","TMED10","TMED9"]:
    sub = df_early[[gene,"braaksc"]+cov_base+cov_cells].dropna()
    if len(sub) < 20: continue
    r_no  = ols_term(sub[gene], sub[["braaksc"]+cov_base], "braaksc")
    r_yes = ols_term(sub[gene], sub[["braaksc"]+cov_base+cov_cells], "braaksc")
    pct   = 100*(r_yes["beta"]-r_no["beta"])/r_no["beta"] if r_no["beta"]!=0 else np.nan
    testI_rows.append(dict(gene=gene,
                           beta_no=r_no["beta"], p_no=r_no["p"],
                           beta_yes=r_yes["beta"], p_yes=r_yes["p"],
                           pct_change=pct, n=r_yes["n"]))
    print(f"  {gene}: β={r_no['beta']:.4f}→{r_yes['beta']:.4f}  ({pct:.1f}%)  n={r_yes['n']}")

testI_df = pd.DataFrame(testI_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST J: Composition-adjusted composite (H8)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST J: Composition-adjusted H8 composite ──")
df["H8_ratio"] = df["TMED9"] - df["TMED2"] - df["TMED10"]
for label, term, outcome in [("Braak","braaksc","H8_ratio"),("AD","AD","H8_ratio")]:
    sub = df[["H8_ratio",term]+cov_base+cov_cells].dropna()
    if len(sub) < 20: continue
    r_no  = ols_term(sub["H8_ratio"], sub[[term]+cov_base], term)
    r_yes = ols_term(sub["H8_ratio"], sub[[term]+cov_base+cov_cells], term)
    print(f"  H8_ratio ~ {label}: β={r_no['beta']:.4f} p={r_no['p']:.4f}  "
          f"→ +cells: β={r_yes['beta']:.4f} p={r_yes['p']:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST K: Killer control — TMED vs NeuronMarkerScore
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST K: TMED vs NeuronMarkerScore ──")
# Build NeuronMarkerScore from control panel genes
ctrl_genes_present = [g for g,e in NEURON_CTRL.items() if e in expr.index]
if ctrl_genes_present:
    ctrl_mat = df[[g for g in ctrl_genes_present if g in df.columns]].copy()
    ctrl_z   = (ctrl_mat - ctrl_mat.mean()) / (ctrl_mat.std() + 1e-10)
    df["NeuronMarkerScore"] = ctrl_z.mean(axis=1)
    print(f"  NeuronMarkerScore from {len(ctrl_genes_present)} genes")

    testK_rows = []
    for gene in ["TMED2","TMED10","TMED9"]:
        sub = df[[gene,"braaksc","NeuronMarkerScore"]+cov_base].dropna()
        # Braak ~ NeuronMarkerScore
        rN = ols_term(sub["braaksc"], sub[["NeuronMarkerScore"]+cov_base], "NeuronMarkerScore")
        # Braak ~ gene
        rG = ols_term(sub["braaksc"], sub[[gene]+cov_base], gene)
        # Braak ~ gene + NeuronMarkerScore
        rGN = ols_term(sub["braaksc"], sub[[gene,"NeuronMarkerScore"]+cov_base], gene)
        testK_rows.append(dict(gene=gene,
                               beta_alone=rG["beta"], p_alone=rG["p"],
                               beta_joint=rGN["beta"], p_joint=rGN["p"],
                               beta_neuron=rN["beta"], p_neuron=rN["p"]))
        print(f"  {gene}: β_alone={rG['beta']:.4f} p={rG['p']:.4f}  "
              f"β_joint(+NeuronScore)={rGN['beta']:.4f} p={rGN['p']:.4f}")
    testK_df = pd.DataFrame(testK_rows)
else:
    print("  No NeuronMarkerScore genes found")
    testK_df = pd.DataFrame()

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(3, 3, figsize=(16, 14))
fig.suptitle("H9: Cell-type composition analysis", fontsize=14, fontweight="bold")

# 1. Test A: composition shift vs Braak
ax = axes[0,0]
for ct, col in zip(["neuron_score","astro_score","micro_score","oligo_score"],
                   ["steelblue","darkorange","green","purple"]):
    sub = df[["braaksc",ct]].dropna()
    means = sub.groupby("braaksc")[ct].mean()
    ses   = sub.groupby("braaksc")[ct].sem()
    ax.errorbar(means.index, means.values, yerr=ses.values, label=ct.replace("_score",""),
                marker="o", capsize=3, color=col)
ax.set_xlabel("Braak stage"); ax.set_ylabel("z-score")
ax.set_title("Cell-type scores across Braak")
ax.legend(fontsize=7)

# 2. Test B: β_Braak across nested models
ax = axes[0,1]
genes = ["TMED2","TMED10","TMED9"]
x = np.arange(len(genes)); w = 0.2
models = ["M1_beta","M2_beta","M3_beta","M4_beta"]
labels = ["M1 raw","M2 +cov","M3 +neuron","M4 +all"]
cols   = ["steelblue","darkorange","green","red"]
for i,(m,lab,c) in enumerate(zip(models,labels,cols)):
    betas = [testB_df.loc[testB_df["gene"]==g,m].values[0] for g in genes]
    ax.bar(x + i*w, betas, w, label=lab, color=c, alpha=0.8)
ax.axhline(0, color="k", linewidth=0.8)
ax.set_xticks(x+w*1.5); ax.set_xticklabels(genes)
ax.set_ylabel("β Braak"); ax.set_title("Nested model comparison")
ax.legend(fontsize=7)

# 3. Test C: correlation heatmap
ax = axes[0,2]
ct_labels = ["neuron","astro","micro","oligo"]
gene_labels = ["TMED2","TMED10","TMED9"]
mat = np.zeros((3,4))
for i,g in enumerate(gene_labels):
    for j,ct in enumerate(cov_cells):
        row = corr_df[(corr_df["gene"]==g)&(corr_df["cell"]==ct)]
        mat[i,j] = row["pearson_r"].values[0] if len(row)>0 else 0
im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(ct_labels)
ax.set_yticks(range(3)); ax.set_yticklabels(gene_labels)
for i in range(3):
    for j in range(4):
        ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=8)
plt.colorbar(im, ax=ax)
ax.set_title("TMED ↔ cell-fraction Pearson r")

# 4. Test D: residualized scatter (3 panels combined)
ax = axes[1,0]
colors = {"TMED2":"steelblue","TMED10":"darkorange","TMED9":"green"}
for gene in ["TMED2","TMED10","TMED9"]:
    d = resid_data[gene]
    means = d.groupby("braaksc")["resid"].mean()
    ax.plot(means.index, means.values, marker="o", label=gene, color=colors[gene])
ax.axhline(0, color="k", linestyle="--", alpha=0.5)
ax.set_xlabel("Braak stage"); ax.set_ylabel("Residual expression")
ax.set_title("Residualized TMED vs Braak\n(composition + demographics removed)")
ax.legend(fontsize=8)

# 5. Test E: neuronal marker β_Braak before/after cells
ax = axes[1,1]
if len(testE_df) > 0:
    x = np.arange(len(testE_df)); w = 0.35
    ax.bar(x-w/2, testE_df["beta_no_cells"], w, label="No cells", color="steelblue", alpha=0.8)
    ax.bar(x+w/2, testE_df["beta_with_cells"], w, label="+cells", color="darkorange", alpha=0.8)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(testE_df["gene"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("β Braak"); ax.set_title("Neuronal marker controls\n(β before/after cell adjustment)")
    ax.legend(fontsize=7)

# 6. Test G: mediation decomposition
ax = axes[1,2]
if len(testG_df) > 0:
    x = np.arange(len(testG_df)); w = 0.25
    ax.bar(x-w, testG_df["total"], w, label="Total", color="steelblue")
    ax.bar(x, testG_df["direct"], w, label="Direct", color="darkorange")
    ax.bar(x+w, testG_df["indirect"], w, label="Indirect (via neuron)", color="green")
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(testG_df["gene"])
    ax.set_ylabel("β Braak"); ax.set_title("Decomposition: total / direct / indirect")
    ax.legend(fontsize=7)

# 7. Test I: pre-symptomatic nested
ax = axes[2,0]
if len(testI_df) > 0:
    x = np.arange(len(testI_df)); w = 0.35
    ax.bar(x-w/2, testI_df["beta_no"], w, label="No cells", color="steelblue", alpha=0.8)
    ax.bar(x+w/2, testI_df["beta_yes"], w, label="+cells", color="darkorange", alpha=0.8)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(testI_df["gene"])
    ax.set_ylabel("β Braak"); ax.set_title("Pre-symptomatic (Braak ≤ 2)")
    ax.legend(fontsize=7)

# 8. Test K: TMED vs NeuronMarkerScore
ax = axes[2,1]
if len(testK_df) > 0:
    x = np.arange(len(testK_df)); w = 0.35
    ax.bar(x-w/2, testK_df["beta_alone"], w, label="Gene alone", color="steelblue", alpha=0.8)
    ax.bar(x+w/2, testK_df["beta_joint"], w, label="+NeuronScore", color="darkorange", alpha=0.8)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(testK_df["gene"])
    ax.set_ylabel("β Braak (gene term)"); ax.set_title("Test K: gene vs NeuronMarkerScore")
    ax.legend(fontsize=7)

# 9. Test F: non-neuronal markers
ax = axes[2,2]
if len(testF_df) > 0:
    x = np.arange(len(testF_df)); w = 0.35
    ax.bar(x-w/2, testF_df["beta_no_cells"], w, label="No cells", color="steelblue", alpha=0.8)
    ax.bar(x+w/2, testF_df["beta_with_cells"], w, label="+cells", color="darkorange", alpha=0.8)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(testF_df["gene"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("β Braak"); ax.set_title("Non-neuronal marker controls")
    ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUT + "h9_plots.png", dpi=150)
plt.close()
print(f"\nPlot saved: {OUT}h9_plots.png")

# ── save tables ───────────────────────────────────────────────────────────────
testA_df.to_csv(OUT+"h9_testA_composition_shift.csv", index=False)
testB_df.to_csv(OUT+"h9_testB_nested_models.csv", index=False)
corr_df.to_csv(OUT+"h9_testC_correlations.csv", index=False)
testE_df.to_csv(OUT+"h9_testE_neuronal_controls.csv", index=False)
testF_df.to_csv(OUT+"h9_testF_nonneuronal_controls.csv", index=False)
testG_df.to_csv(OUT+"h9_testG_decomposition.csv", index=False)
testI_df.to_csv(OUT+"h9_testI_presymptomatic.csv", index=False)
if len(testK_df)>0:
    testK_df.to_csv(OUT+"h9_testK_killer_control.csv", index=False)

print("\n══ H9 COMPLETE ══")
