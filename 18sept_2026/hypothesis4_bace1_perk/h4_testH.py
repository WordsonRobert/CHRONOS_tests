"""
CHRONOS - Hypothesis 4, Test H: BACE1 RNA vs Protein
======================================================
Compare BACE1 RNA and protein trajectories across Braak.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"
PROT_FILE = DATA + "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"
BACE1_RNA_ID = "ENSG00000186318"

# ── 1. LOAD ────────────────────────────────────────────────────────────────────

print("Loading RNA...")
expr = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
bace1_rna = expr.loc[BACE1_RNA_ID]; bace1_rna.name = "BACE1_RNA"

print("Loading proteomics...")
prot = pd.read_csv(PROT_FILE, index_col=0)
bace1_prot_row = next((i for i in prot.index if "BACE1" in str(i).upper()), None)
print(f"  BACE1 protein row: {bace1_prot_row}")
bace1_prot = prot.loc[bace1_prot_row]; bace1_prot.name = "BACE1_prot"

print("Loading metadata...")
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")
prot_meta = pd.read_csv(DATA + "ROSMAP_assay_proteomics_TMTquantitation_metadata.csv")

clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+","90")
for col in ["age_at_visit_max","braaksc","msex","pmi","cogdx"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi","cogdx"]

# ── 2. LINK RNA ────────────────────────────────────────────────────────────────

rna_bio = biospec[biospec["specimenID"].isin(bace1_rna.index)][["specimenID","individualID"]].drop_duplicates()
rna_batch = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
rna_linked = rna_bio.merge(clinical[clin_cols], on="individualID", how="left").merge(rna_batch, on="specimenID", how="left").set_index("specimenID")
rna_df = bace1_rna.to_frame().join(rna_linked).dropna(subset=["braaksc","msex","age_at_visit_max","pmi"])
rna_df["braaksc"] = rna_df["braaksc"].astype(float)
print(f"\n  RNA linked: {len(rna_df)} samples")

# ── 3. LINK PROTEIN ────────────────────────────────────────────────────────────

# protein columns = batchChannel (e.g. b01.127C)
# prot_meta: batchChannel → specimenID → individualID
prot_meta_clean = prot_meta[~prot_meta["isAssayControl"].astype(bool)].copy()
prot_meta_clean = prot_meta_clean[["specimenID","batchChannel"]].drop_duplicates()
prot_bio = biospec[biospec["specimenID"].isin(prot_meta_clean["specimenID"])][["specimenID","individualID"]].drop_duplicates()
prot_link = prot_meta_clean.merge(prot_bio, on="specimenID", how="left")

# bace1_prot index = batchChannel
prot_series = bace1_prot.reset_index()
prot_series.columns = ["batchChannel","BACE1_prot"]

prot_full = prot_series.merge(prot_link, on="batchChannel", how="left")
prot_full = prot_full.merge(clinical[clin_cols], on="individualID", how="left")
prot_df = prot_full.dropna(subset=["braaksc","msex","age_at_visit_max","pmi","BACE1_prot"]).copy()
prot_df["braaksc"] = prot_df["braaksc"].astype(float)
print(f"  Protein linked: {len(prot_df)} samples")

# ── 4. PAIRED DONORS ──────────────────────────────────────────────────────────

rna_ind_map = rna_linked[["individualID"]].join(rna_df[["BACE1_RNA"]], how="inner")
prt_ind_map = prot_df[["individualID","BACE1_prot"]].dropna()
paired = rna_ind_map.reset_index().merge(prt_ind_map, on="individualID", how="inner")
print(f"  Paired donors (RNA+protein): {len(paired)}")

# ── 5. OLS ────────────────────────────────────────────────────────────────────

def ols_term(y_s, X_df, term):
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1)|np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    if len(y_v) < 10: return dict(beta=np.nan,se=np.nan,p=np.nan,n=len(y_v),ci_lo=np.nan,ci_hi=np.nan)
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm@coeffs
    n,p = len(y_v), Xm.shape[1]
    mse = np.sum(resid**2)/(n-p)
    try: cov = mse*np.linalg.inv(Xm.T@Xm)
    except: return dict(beta=np.nan,se=np.nan,p=np.nan,n=n,ci_lo=np.nan,ci_hi=np.nan)
    se = np.sqrt(np.diag(cov))
    cols = ["intercept"]+list(X_df.columns)
    idx2 = cols.index(term)
    beta,se_b = coeffs[idx2],se[idx2]
    t = beta/se_b
    pv = 2*stats.t.sf(abs(t),df=n-p)
    return dict(beta=beta,se=se_b,ci_lo=beta-1.96*se_b,ci_hi=beta+1.96*se_b,p=pv,n=int(mask.sum()))

# RNA models
batch_dum = pd.get_dummies(rna_df["rnaBatch"], prefix="batch", drop_first=True)
rna_df2 = pd.concat([rna_df, batch_dum], axis=1)
batch_cols = list(batch_dum.columns)
base_rna = ["age_at_visit_max","msex","pmi"] + batch_cols
base_prot = ["age_at_visit_max","msex","pmi"]

print("\n=== BACE1 RNA ~ Braak ===")
res_rna = ols_term(rna_df2["BACE1_RNA"], rna_df2[["braaksc"]+base_rna], "braaksc")
print(f"  β={res_rna['beta']:.4f}, p={res_rna['p']:.3e}, n={res_rna['n']}")

print("\n=== BACE1 Protein ~ Braak ===")
res_prot = ols_term(prot_df["BACE1_prot"], prot_df[["braaksc"]+base_prot], "braaksc")
print(f"  β={res_prot['beta']:.4f}, p={res_prot['p']:.3e}, n={res_prot['n']}")

# AD models
rna_df2["AD"] = (rna_df2["cogdx"] >= 4).astype(float)
prot_df["AD"] = (prot_df["cogdx"] >= 4).astype(float)
res_rna_ad  = ols_term(rna_df2["BACE1_RNA"],  rna_df2[["braaksc","AD"]+base_rna],  "AD")
res_prot_ad = ols_term(prot_df["BACE1_prot"], prot_df[["braaksc","AD"]+base_prot], "AD")
print(f"\n  RNA  ~ AD (adj Braak): β={res_rna_ad['beta']:.4f}, p={res_rna_ad['p']:.3e}")
print(f"  Prot ~ AD (adj Braak): β={res_prot_ad['beta']:.4f}, p={res_prot_ad['p']:.3e}")

# Paired correlation
if len(paired) > 10:
    r_rp, p_rp = stats.pearsonr(paired["BACE1_RNA"], paired["BACE1_prot"])
    print(f"\n  RNA vs Protein (n={len(paired)}): r={r_rp:.3f}, p={p_rp:.3e}")
else:
    r_rp, p_rp = np.nan, np.nan
    print(f"\n  Paired donors: {len(paired)} — insufficient for correlation")

# Pattern
rna_dir  = "DOWN" if res_rna["beta"]<0  and res_rna["p"]<0.05  else \
           "UP"   if res_rna["beta"]>0  and res_rna["p"]<0.05  else "FLAT"
prot_dir = "DOWN" if res_prot["beta"]<0 and res_prot["p"]<0.05 else \
           "UP"   if res_prot["beta"]>0 and res_prot["p"]<0.05 else "FLAT"
interp = {
    ("DOWN","DOWN"): "Both reduced — broadly decreased BACE1. Does not support uORF upregulation.",
    ("DOWN","FLAT"): "RNA down, protein stable — possible post-transcriptional compensation. Weakly compatible with uORF.",
    ("DOWN","UP"):   "RNA down, protein UP — strong post-transcriptional regulation. Most compatible with uORF mechanism.",
    ("FLAT","UP"):   "RNA stable, protein up — post-transcriptional upregulation.",
    ("FLAT","FLAT"): "Both flat — no bulk evidence for BACE1 change.",
    ("UP","UP"):     "Both up — transcriptional + protein increase.",
    ("UP","FLAT"):   "RNA up, protein flat — transcript change doesn't translate to protein.",
}
print(f"\n=== PATTERN: RNA={rna_dir}, Protein={prot_dir} ===")
print(f"  {interp.get((rna_dir,prot_dir),'See individual results.')}")

# ── 6. SAVE ───────────────────────────────────────────────────────────────────

pd.DataFrame([
    {"measure":"RNA","model":"~Braak",**res_rna},
    {"measure":"RNA","model":"~AD",**res_rna_ad},
    {"measure":"Protein","model":"~Braak",**res_prot},
    {"measure":"Protein","model":"~AD",**res_prot_ad},
]).to_csv("h4_testH_results.csv", index=False)

# ── 7. PLOT ───────────────────────────────────────────────────────────────────

print("\nPlotting...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("CHRONOS — BACE1 RNA vs Protein ~ Braak (Test H)", fontsize=13)

# RNA trajectory
ax = axes[0]
bs = rna_df.groupby("braaksc")["BACE1_RNA"].agg(["mean","sem"]).reset_index()
ax.errorbar(bs["braaksc"], bs["mean"], yerr=bs["sem"], fmt="o-", color="firebrick", linewidth=2, capsize=4)
xl = np.linspace(0,6,100)
m,b = np.polyfit(rna_df["braaksc"].values, rna_df["BACE1_RNA"].values, 1)
ax.plot(xl, m*xl+b, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.set_xlabel("Braak stage"); ax.set_ylabel("BACE1 RNA (logCPM)")
ax.set_title(f"BACE1 RNA ~ Braak\nβ={res_rna['beta']:.4f}, p={res_rna['p']:.2e} → {rna_dir}", fontsize=10)

# Protein trajectory
ax = axes[1]
if len(prot_df) > 0:
    bs2 = prot_df.groupby("braaksc")["BACE1_prot"].agg(["mean","sem"]).reset_index()
    ax.errorbar(bs2["braaksc"], bs2["mean"], yerr=bs2["sem"], fmt="o-", color="steelblue", linewidth=2, capsize=4)
    if len(prot_df) > 5:
        m2,b2 = np.polyfit(prot_df["braaksc"].values, prot_df["BACE1_prot"].values, 1)
        ax.plot(xl, m2*xl+b2, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.set_xlabel("Braak stage"); ax.set_ylabel("BACE1 Protein (log2 ratio)")
    ax.set_title(f"BACE1 Protein ~ Braak\nβ={res_prot['beta']:.4f}, p={res_prot['p']:.2e} → {prot_dir}", fontsize=10)
else:
    ax.text(0.5, 0.5, "Protein linking failed\ncheck metadata", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("BACE1 Protein — no data", fontsize=10)

# RNA vs protein scatter
ax = axes[2]
if len(paired) > 10:
    ax.scatter(paired["BACE1_RNA"], paired["BACE1_prot"], s=10, alpha=0.5, color="purple")
    m3,b3 = np.polyfit(paired["BACE1_RNA"].values, paired["BACE1_prot"].values, 1)
    xl3 = np.linspace(paired["BACE1_RNA"].min(), paired["BACE1_RNA"].max(), 100)
    ax.plot(xl3, m3*xl3+b3, color="black", linewidth=1.5)
    ax.set_title(f"RNA vs Protein (n={len(paired)})\nr={r_rp:.3f}, p={p_rp:.2e}", fontsize=10)
    ax.set_xlabel("BACE1 RNA (logCPM)"); ax.set_ylabel("BACE1 Protein")
else:
    ax.text(0.5, 0.5, f"Only {len(paired)} paired donors\nInsufficient for correlation",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_title("RNA vs Protein — insufficient pairs", fontsize=10)

plt.tight_layout()
plt.savefig("h4_testH.png", dpi=200, bbox_inches="tight")
print("Saved: h4_testH.png")
print("\nDone.")
