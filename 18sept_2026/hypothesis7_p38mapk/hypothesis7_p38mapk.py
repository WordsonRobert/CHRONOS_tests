"""
CHRONOS - Hypothesis 7: p38-MAPK (Route 4)
===========================================
Tests A-H: MAPK14 and p38 pathway ~ Braak
Key challenge: severe microglial confounding (DAM activation in AD).

Tests:
  A - MAPK14 ~ Braak (raw + adjusted + cell composition)
  B - Nested models highlighting microglial confounding
  C - MAPK14 ~ microglia_score (direct confound test)
  D - MAPK14 ~ neuron_score (neuronal signal test)
  E - AD vs NCI (secondary)
  F - p38 pathway panel (MAP2K3, MAP2K6, MAPKAPK2)
  G - (blocked: requires snRNA)
  H - MAPK14 × microglia interaction (exploratory)

Input : ROSMAP_DLPFC_logCPM.tsv, ROSMAP_clinical.csv, biospecimen, rna_meta
Output: h7_results.csv, h7_plots.png
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

PATHWAY_GENES = {
    "MAPK14":   "ENSG00000112062",  # p38α — primary
    "MAP2K3":   "ENSG00000034152",  # upstream kinase
    "MAP2K6":   "ENSG00000108984",  # upstream kinase
    "MAPKAPK2": "ENSG00000162889",  # downstream substrate
}

NEURON_MARKERS = ["ENSG00000102003","ENSG00000067715",
                  "ENSG00000132639","ENSG00000008056","ENSG00000157542"]
ASTRO_MARKERS  = ["ENSG00000131095","ENSG00000171885"]
MICRO_MARKERS  = ["ENSG00000197249","ENSG00000101439"]
OLIGO_MARKERS  = ["ENSG00000197971","ENSG00000123560"]

# ── 1. LOAD ────────────────────────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+","90")
for col in ["age_at_visit_max","braaksc","msex","pmi","cogdx"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

sample_ids = expr.columns.tolist()
rna_bio   = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()
rna_batch = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi","cogdx"]
meta = rna_bio.merge(clinical[clin_cols], on="individualID", how="left")\
              .merge(rna_batch, on="specimenID", how="left")\
              .set_index("specimenID")

# ── 2. EXTRACT GENES + CELL SCORES ────────────────────────────────────────────

print("\nChecking genes...")
gene_df = {}
for name, eid in PATHWAY_GENES.items():
    if eid in expr.index:
        gene_df[name] = expr.loc[eid]
        print(f"  {name}: found")
    else:
        print(f"  {name}: NOT FOUND")
gene_df = pd.DataFrame(gene_df)

def ct_score(markers, name):
    present = [m for m in markers if m in expr.index]
    if not present: return None
    vals = expr.loc[present]
    z = (vals - vals.mean(axis=1).values[:,None]) / (vals.std(axis=1).values[:,None] + 1e-10)
    s = z.mean(axis=0); s.name = name; return s

neuron_s = ct_score(NEURON_MARKERS,"neuron_score")
astro_s  = ct_score(ASTRO_MARKERS, "astro_score")
micro_s  = ct_score(MICRO_MARKERS, "micro_score")
oligo_s  = ct_score(OLIGO_MARKERS, "oligo_score")

df = gene_df.join(meta).join(neuron_s).join(astro_s).join(micro_s).join(oligo_s)
df = df.dropna(subset=["braaksc","age_at_visit_max","msex","pmi"]).copy()
df["braaksc"] = df["braaksc"].astype(float)

batch_dum  = pd.get_dummies(df["rnaBatch"], prefix="batch", drop_first=True)
df = pd.concat([df, batch_dum], axis=1)
batch_cols = list(batch_dum.columns)
cell_cols  = [c for c in ["neuron_score","astro_score","micro_score","oligo_score"] if c in df.columns]
base_covs  = ["age_at_visit_max","msex","pmi"] + batch_cols

print(f"\n  Analysis samples: {len(df)}")

# ── OLS HELPER ─────────────────────────────────────────────────────────────────

def ols_term(y_s, X_df, term):
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1)|np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    if len(y_v) < 10:
        return dict(beta=np.nan,se=np.nan,p=np.nan,n=len(y_v),ci_lo=np.nan,ci_hi=np.nan)
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm@coeffs
    n,p = len(y_v), Xm.shape[1]
    mse = np.sum(resid**2)/(n-p)
    try: cov = mse*np.linalg.inv(Xm.T@Xm)
    except: return dict(beta=np.nan,se=np.nan,p=np.nan,n=n,ci_lo=np.nan,ci_hi=np.nan)
    se = np.sqrt(np.diag(cov))
    cols = ["intercept"]+list(X_df.columns)
    idx = cols.index(term)
    beta,se_b = coeffs[idx],se[idx]
    t = beta/se_b
    pv = 2*stats.t.sf(abs(t),df=n-p)
    return dict(beta=beta,se=se_b,ci_lo=beta-1.96*se_b,ci_hi=beta+1.96*se_b,p=pv,n=int(mask.sum()))

# ── TEST A: MAPK14 ~ BRAAK ────────────────────────────────────────────────────

print("\n=== TEST A: MAPK14 ~ Braak ===")
r_raw  = ols_term(df["MAPK14"], df[["braaksc"]], "braaksc")
r_adj  = ols_term(df["MAPK14"], df[["braaksc"]+base_covs], "braaksc")
r_full = ols_term(df["MAPK14"], df[["braaksc"]+base_covs+cell_cols], "braaksc")
print(f"  Raw  : β={r_raw['beta']:.4f}, p={r_raw['p']:.3e}")
print(f"  Adj  : β={r_adj['beta']:.4f}, p={r_adj['p']:.3e}")
print(f"  Full : β={r_full['beta']:.4f}, p={r_full['p']:.3e}")

# ── TEST B: NESTED MODELS — MICROGLIA KEY ──────────────────────────────────────

print("\n=== TEST B: Nested models — microglial confounding ===")
no_micro  = [c for c in cell_cols if c != "micro_score"]
m1 = ols_term(df["MAPK14"], df[["braaksc"]], "braaksc")
m2 = ols_term(df["MAPK14"], df[["braaksc"]+base_covs], "braaksc")
m3 = ols_term(df["MAPK14"], df[["braaksc"]+base_covs+no_micro], "braaksc")
m4 = ols_term(df["MAPK14"], df[["braaksc"]+base_covs+cell_cols], "braaksc")

nested = [("M1 (Braak only)",m1),("M2 (+age/sex/PMI/batch)",m2),
          ("M3 (+neuron/astro/oligo)",m3),("M4 (+microglia)",m4)]
for label, res in nested:
    print(f"  {label}: β={res['beta']:.4f}, p={res['p']:.3e}")

print(f"\n  β_Braak change when adding microglia: {m3['beta']:.4f} → {m4['beta']:.4f}")
print(f"  (attenuation = {(m3['beta']-m4['beta'])/abs(m3['beta'])*100:.1f}% if same sign)")

# ── TEST C: MAPK14 ~ MICROGLIA ────────────────────────────────────────────────

print("\n=== TEST C: MAPK14 ~ microglia score ===")
r_mic_raw = ols_term(df["MAPK14"], df[["micro_score"]], "micro_score")
r_mic_adj = ols_term(df["MAPK14"], df[["micro_score","braaksc"]+base_covs], "micro_score")
rho_mic, p_mic = stats.pearsonr(df["MAPK14"].dropna(),
                                df.loc[df["MAPK14"].notna(),"micro_score"].dropna())
print(f"  r(MAPK14, micro_score) = {rho_mic:.3f}, p={p_mic:.3e}")
print(f"  MAPK14 ~ micro_score (raw): β={r_mic_raw['beta']:.4f}, p={r_mic_raw['p']:.3e}")
print(f"  MAPK14 ~ micro_score (+Braak): β={r_mic_adj['beta']:.4f}, p={r_mic_adj['p']:.3e}")

# ── TEST D: MAPK14 ~ NEURON ────────────────────────────────────────────────────

print("\n=== TEST D: MAPK14 ~ neuron score ===")
r_neu_raw = ols_term(df["MAPK14"], df[["neuron_score"]], "neuron_score")
r_neu_adj = ols_term(df["MAPK14"], df[["neuron_score","braaksc"]+base_covs], "neuron_score")
rho_neu, p_neu = stats.pearsonr(df["MAPK14"].dropna(),
                                df.loc[df["MAPK14"].notna(),"neuron_score"].dropna())
print(f"  r(MAPK14, neuron_score) = {rho_neu:.3f}, p={p_neu:.3e}")
print(f"  MAPK14 ~ neuron_score (raw): β={r_neu_raw['beta']:.4f}, p={r_neu_raw['p']:.3e}")
print(f"  MAPK14 ~ neuron_score (+Braak): β={r_neu_adj['beta']:.4f}, p={r_neu_adj['p']:.3e}")

# ── TEST E: AD VS NCI ─────────────────────────────────────────────────────────

print("\n=== TEST E: AD vs NCI ===")
df["AD"] = (df["cogdx"] >= 4).astype(float)
r_ad = ols_term(df["MAPK14"], df[["AD","braaksc"]+base_covs+cell_cols], "AD")
print(f"  MAPK14 ~ AD (adj Braak+cells): β={r_ad['beta']:.4f}, p={r_ad['p']:.3e}")

# ── TEST F: PATHWAY PANEL ─────────────────────────────────────────────────────

print("\n=== TEST F: p38 pathway panel ~ Braak ===")
panel_results = []
for gene in ["MAPK14","MAP2K3","MAP2K6","MAPKAPK2"]:
    if gene not in df.columns: continue
    r = ols_term(df[gene], df[["braaksc"]+base_covs+cell_cols], "braaksc")
    print(f"  {gene}: β={r['beta']:.4f}, p={r['p']:.3e}")
    panel_results.append({"gene":gene,**r})

# FDR
ps = [r["p"] for r in panel_results]
_, fdrs, _, _ = multipletests(ps, method="fdr_bh")
for r,fdr in zip(panel_results,fdrs):
    r["fdr"] = fdr
    print(f"    FDR: {fdr:.3f}")

# ── TEST H: MAPK14 × MICROGLIA INTERACTION ────────────────────────────────────

print("\n=== TEST H: MAPK14 × microglia interaction (exploratory) ===")
df["braak_x_micro"] = df["braaksc"] * df["micro_score"]
r_int = ols_term(df["MAPK14"],
                 df[["braaksc","micro_score","braak_x_micro"]+base_covs],
                 "braak_x_micro")
print(f"  Braak × micro interaction: β={r_int['beta']:.4f}, p={r_int['p']:.3e}")

# ── SAVE ───────────────────────────────────────────────────────────────────────

all_res = [
    {"test":"A_raw","gene":"MAPK14","term":"braaksc",**r_raw},
    {"test":"A_adj","gene":"MAPK14","term":"braaksc",**r_adj},
    {"test":"A_full","gene":"MAPK14","term":"braaksc",**r_full},
    {"test":"B_M1","gene":"MAPK14","term":"braaksc",**m1},
    {"test":"B_M2","gene":"MAPK14","term":"braaksc",**m2},
    {"test":"B_M3","gene":"MAPK14","term":"braaksc",**m3},
    {"test":"B_M4","gene":"MAPK14","term":"braaksc",**m4},
    {"test":"C_micro_raw","gene":"MAPK14","term":"micro_score",**r_mic_raw},
    {"test":"D_neuron_raw","gene":"MAPK14","term":"neuron_score",**r_neu_raw},
    {"test":"E_AD","gene":"MAPK14","term":"AD",**r_ad},
]
pd.DataFrame(all_res).to_csv("h7_results.csv", index=False)
pd.DataFrame(panel_results).to_csv("h7_pathway_panel.csv", index=False)
print("\nSaved CSVs.")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("CHRONOS — p38-MAPK Route (Hypothesis 7)", fontsize=13)

# Plot 1: MAPK14 trajectory
ax = axes[0,0]
bs = df.groupby("braaksc")["MAPK14"].agg(["mean","sem"]).reset_index()
ax.errorbar(bs["braaksc"], bs["mean"], yerr=bs["sem"],
            fmt="o-", color="firebrick", linewidth=2, capsize=4)
ax.set_xlabel("Braak stage"); ax.set_ylabel("MAPK14 (logCPM)")
ax.set_title(f"MAPK14 ~ Braak\nraw β={r_raw['beta']:.4f}, p={r_raw['p']:.2e}", fontsize=10)

# Plot 2: Nested model betas
ax = axes[0,1]
labels = ["M1\n(Braak)","M2\n(+cov)","M3\n(+cells\nno micro)","M4\n(+micro)"]
betas  = [m1["beta"],m2["beta"],m3["beta"],m4["beta"]]
colors = ["steelblue","steelblue","orange","firebrick"]
ax.bar(labels, betas, color=colors, alpha=0.8)
ax.axhline(0,color="black",linewidth=0.8)
ax.set_ylabel("β (Braak term)")
ax.set_title("Test B: Nested models\nMicroglial confounding", fontsize=10)

# Plot 3: MAPK14 vs microglia scatter
ax = axes[0,2]
x = df["micro_score"].values; y = df["MAPK14"].values
mask = ~(np.isnan(x)|np.isnan(y))
ax.scatter(x[mask], y[mask], s=5, alpha=0.3, color="purple")
m,b = np.polyfit(x[mask], y[mask], 1)
xl = np.linspace(x[mask].min(), x[mask].max(), 100)
ax.plot(xl, m*xl+b, color="black", linewidth=1.5)
ax.set_xlabel("Microglia score"); ax.set_ylabel("MAPK14 (logCPM)")
ax.set_title(f"MAPK14 vs Microglia\nr={rho_mic:.3f}, p={p_mic:.1e}", fontsize=10)

# Plot 4: MAPK14 vs neuron scatter
ax = axes[1,0]
x = df["neuron_score"].values
mask = ~(np.isnan(x)|np.isnan(y))
ax.scatter(x[mask], y[mask], s=5, alpha=0.3, color="steelblue")
m,b = np.polyfit(x[mask], y[mask], 1)
xl = np.linspace(x[mask].min(), x[mask].max(), 100)
ax.plot(xl, m*xl+b, color="black", linewidth=1.5)
ax.set_xlabel("Neuron score"); ax.set_ylabel("MAPK14 (logCPM)")
ax.set_title(f"MAPK14 vs Neuron\nr={rho_neu:.3f}, p={p_neu:.1e}", fontsize=10)

# Plot 5: Pathway panel betas
ax = axes[1,1]
pan_df = pd.DataFrame(panel_results)
if not pan_df.empty:
    cols_p = ["steelblue" if g=="MAPK14" else "grey" for g in pan_df["gene"]]
    bars = ax.bar(pan_df["gene"], pan_df["beta"], color=cols_p, alpha=0.8)
    ax.axhline(0,color="black",linewidth=0.8)
    for bar, pv in zip(bars, pan_df["p"]):
        sig = "***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else ""
        if sig: ax.text(bar.get_x()+bar.get_width()/2,
                        bar.get_height()+0.001, sig, ha="center", fontsize=10)
    ax.set_ylabel("β (Braak, full model)")
    ax.set_title("Test F: p38 pathway panel\n~ Braak (full model)", fontsize=10)

# Plot 6: Microglia score ~ Braak
ax = axes[1,2]
bs2 = df.groupby("braaksc")["micro_score"].agg(["mean","sem"]).reset_index()
ax.errorbar(bs2["braaksc"], bs2["mean"], yerr=bs2["sem"],
            fmt="o-", color="purple", linewidth=2, capsize=4)
r_mb, p_mb = stats.pearsonr(df["braaksc"].values, df["micro_score"].values)
ax.set_xlabel("Braak stage"); ax.set_ylabel("Microglia score")
ax.set_title(f"Microglia ~ Braak\nr={r_mb:.3f}, p={p_mb:.1e}", fontsize=10)

plt.tight_layout()
plt.savefig("h7_plots.png", dpi=200, bbox_inches="tight")
print("Saved: h7_plots.png")
print("\nDone.")
