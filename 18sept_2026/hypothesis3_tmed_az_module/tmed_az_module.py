"""
CHRONOS - Hypothesis 3: TMED9/TMED10 vs Active Zone Module
===========================================================
Test A: Do TMED9/TMED10 correlate with the presynaptic AZ module?
Test B: Is the correlation specific (AZ vs control module)?
Test C: Does TMED10 expression become Braak-independent after conditioning on AZ score?
        (Replicates the APP result from hypothesis 2)

Input : ROSMAP_DLPFC_logCPM.tsv, ROSMAP_clinical.csv, ROSMAP_biospecimen_metadata.csv
Output: tmed_az_results.csv, tmed_az_plots.png
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

TARGET_GENES = {
    "TMED9":  "ENSG00000184840",
    "TMED10": "ENSG00000170348",
    "TMED2":  "ENSG00000086598",   # include for completeness
}

AZ_MODULE = {
    "CASK":   "ENSG00000147044",
    "STX1A":  "ENSG00000106089",
    "SNAP25": "ENSG00000132639",
    "VAMP2":  "ENSG00000220205",
    "SYN1":   "ENSG00000008056",
    "SYN2":   "ENSG00000157542",
    "RIMS1":  "ENSG00000079841",
    "UNC13A": "ENSG00000130477",
    "DNM1":   "ENSG00000106976",
    "SV2A":   "ENSG00000197912",
    "SYP":    "ENSG00000102003",
    "SYT1":   "ENSG00000067715",
}

CTRL_MODULE = {
    "ACTB":  "ENSG00000075624",
    "GAPDH": "ENSG00000111640",
    "LDHA":  "ENSG00000134333",
    "VIM":   "ENSG00000026025",
    "GFAP":  "ENSG00000131095",
    "AIF1":  "ENSG00000197249",
    "MBP":   "ENSG00000197971",
    "PLP1":  "ENSG00000123560",
}

# Reference: APP results from hypothesis 2 (for comparison)
APP_AZ_r   = 0.871
APP_CTRL_r = 0.475

COLORS = {"TMED9": "firebrick", "TMED10": "seagreen", "TMED2": "steelblue"}

# ── 1. LOAD ────────────────────────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")

clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+","90")
clinical["age_at_visit_max"] = pd.to_numeric(clinical["age_at_visit_max"], errors="coerce")
for col in ["braaksc","msex","pmi"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

sample_ids = expr.columns.tolist()
rna_bio = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()
clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi"]
meta = rna_bio.merge(clinical[clin_cols], on="individualID", how="left").set_index("specimenID")

# ── 2. EXTRACT GENES ───────────────────────────────────────────────────────────

def get_genes(gene_dict):
    rows = {n: expr.loc[e] for n,e in gene_dict.items() if e in expr.index}
    return pd.DataFrame(rows)

target_df = get_genes(TARGET_GENES)
az_df     = get_genes(AZ_MODULE)
ctrl_df   = get_genes(CTRL_MODULE)

print(f"  Target genes   : {list(target_df.columns)}")
print(f"  AZ module genes: {az_df.shape[1]}/{len(AZ_MODULE)}")
print(f"  Ctrl genes     : {ctrl_df.shape[1]}/{len(CTRL_MODULE)}")

# ── 3. MODULE SCORES ───────────────────────────────────────────────────────────

def module_score(df):
    z = (df - df.mean()) / (df.std() + 1e-10)
    return z.mean(axis=1)

az_score   = module_score(az_df);   az_score.name   = "AZscore"
ctrl_score = module_score(ctrl_df); ctrl_score.name = "CTRLscore"

# ── TEST A+B: CORRELATIONS ─────────────────────────────────────────────────────

print("\n=== TEST A+B: Correlations with AZ and control modules ===")
print(f"  (Reference: APP r_AZ=0.871, r_ctrl=0.475)\n")

results = []
for gene in TARGET_GENES:
    if gene not in target_df.columns:
        continue
    y = target_df[gene]
    common = y.index.intersection(az_score.index)
    y_v    = y.loc[common].values
    az_v   = az_score.loc[common].values
    ct_v   = ctrl_score.loc[common].values

    r_az,  p_az  = stats.pearsonr(y_v, az_v)
    r_ct,  p_ct  = stats.pearsonr(y_v, ct_v)
    specificity  = r_az - r_ct   # how much stronger is AZ vs control

    print(f"  {gene}:")
    print(f"    vs AZ module  : r={r_az:.3f}, p={p_az:.2e}")
    print(f"    vs Ctrl module: r={r_ct:.3f}, p={p_ct:.2e}")
    print(f"    Specificity (AZ - ctrl): {specificity:.3f}")
    print()

    results.append({"gene": gene, "r_az": r_az, "p_az": p_az,
                    "r_ctrl": r_ct, "p_ctrl": p_ct,
                    "specificity": specificity})

res_df = pd.DataFrame(results)

# ── TEST C: BRAAK INDEPENDENCE AFTER AZ CONDITIONING ─────────────────────────

print("=== TEST C: Does Braak effect survive AZ conditioning? ===")
all_df = target_df.join(az_score).join(meta).dropna(subset=["braaksc","msex","age_at_visit_max","pmi"])
base_covs = ["age_at_visit_max","msex","pmi"]

def ols_term(y_s, X_df, term):
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1)|np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm@coeffs
    n,p = len(y_v), Xm.shape[1]
    mse = np.sum(resid**2)/(n-p)
    try:
        cov = mse*np.linalg.inv(Xm.T@Xm)
    except:
        return dict(beta=np.nan,se=np.nan,p=np.nan,n=n)
    se = np.sqrt(np.diag(cov))
    cols = ["intercept"]+list(X_df.columns)
    idx = cols.index(term)
    beta,se_b = coeffs[idx],se[idx]
    t = beta/se_b
    pv = 2*stats.t.sf(abs(t),df=n-p)
    return dict(beta=beta,se=se_b,
                ci_lo=beta-1.96*se_b,ci_hi=beta+1.96*se_b,p=pv,n=int(mask.sum()))

testC = []
for gene in ["TMED9","TMED10","TMED2"]:
    if gene not in all_df.columns: continue

    # M1: gene ~ Braak + covariates (no AZ)
    m1 = ols_term(all_df[gene], all_df[["braaksc"]+base_covs], "braaksc")

    # M2: gene ~ AZscore + Braak + covariates
    m2_az    = ols_term(all_df[gene], all_df[["AZscore","braaksc"]+base_covs], "AZscore")
    m2_braak = ols_term(all_df[gene], all_df[["AZscore","braaksc"]+base_covs], "braaksc")

    # M3: gene ~ AZscore + covariates (no Braak)
    m3 = ols_term(all_df[gene], all_df[["AZscore"]+base_covs], "AZscore")

    print(f"  {gene}:")
    print(f"    M1 Braak only : β_Braak={m1['beta']:.4f}, p={m1['p']:.3e}")
    print(f"    M2 AZ+Braak   : β_AZ={m2_az['beta']:.4f} p={m2_az['p']:.3e} | β_Braak={m2_braak['beta']:.4f} p={m2_braak['p']:.3e}")
    print(f"    M3 AZ only    : β_AZ={m3['beta']:.4f}, p={m3['p']:.3e}")
    print()

    testC.append({"gene":gene,
                  "M1_beta_braak":m1["beta"], "M1_p_braak":m1["p"],
                  "M2_beta_az":m2_az["beta"], "M2_p_az":m2_az["p"],
                  "M2_beta_braak":m2_braak["beta"], "M2_p_braak":m2_braak["p"],
                  "M3_beta_az":m3["beta"], "M3_p_az":m3["p"]})

# ── TEST C extra: AZscore ~ Braak ─────────────────────────────────────────────

az_braak = ols_term(all_df["AZscore"], all_df[["braaksc"]+base_covs], "braaksc")
print(f"  AZscore ~ Braak: β={az_braak['beta']:.4f}, p={az_braak['p']:.3e}")

# ── SAVE ───────────────────────────────────────────────────────────────────────

res_df.to_csv("tmed_az_correlations.csv", index=False)
pd.DataFrame(testC).to_csv("tmed_az_braak_models.csv", index=False)
print("\nSaved CSVs.")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle("CHRONOS — TMED9/TMED10 vs Active Zone Module", fontsize=13)

genes = [g for g in ["TMED9","TMED10","TMED2"] if g in target_df.columns]

# Row 1: scatter vs AZ module
for i, gene in enumerate(genes):
    ax = axes[0, i]
    y  = target_df[gene]
    common = y.index.intersection(az_score.index)
    xv = az_score.loc[common].values
    yv = y.loc[common].values
    r, p = stats.pearsonr(xv, yv)
    ax.scatter(xv, yv, s=6, alpha=0.4, color=COLORS[gene])
    m, b = np.polyfit(xv, yv, 1)
    xl = np.linspace(xv.min(), xv.max(), 100)
    ax.plot(xl, m*xl+b, color="black", linewidth=1.5)
    ax.set_xlabel("Presynaptic AZ module score", fontsize=9)
    ax.set_ylabel(f"{gene} (logCPM)", fontsize=9)
    ax.set_title(f"{gene} vs AZ module\nr={r:.3f}, p={p:.1e}", fontsize=10)
    # add APP reference line as annotation
    ax.annotate(f"APP ref: r=0.871", xy=(0.05,0.92), xycoords="axes fraction",
                fontsize=7, color="grey")

# Row 2 left: specificity bar chart (AZ vs ctrl r, all genes + APP reference)
ax = axes[1, 0]
all_genes = genes + ["APP"]
r_az_vals   = [res_df[res_df.gene==g]["r_az"].values[0]   if g in res_df.gene.values else APP_AZ_r   for g in all_genes]
r_ctrl_vals = [res_df[res_df.gene==g]["r_ctrl"].values[0] if g in res_df.gene.values else APP_CTRL_r for g in all_genes]
x = np.arange(len(all_genes)); w = 0.35
ax.bar(x-w/2, r_az_vals,   w, label="AZ module",   color="steelblue", alpha=0.8)
ax.bar(x+w/2, r_ctrl_vals, w, label="Ctrl module", color="grey",      alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(all_genes, fontsize=10)
ax.axhline(0, color="black", linewidth=0.5)
ax.set_ylabel("Pearson r"); ax.set_title("AZ vs Ctrl specificity\n(APP = reference)", fontsize=10)
ax.legend(fontsize=8)

# Row 2 middle: Braak beta before/after AZ conditioning
ax = axes[1, 1]
tc_df = pd.DataFrame(testC)
x = np.arange(len(tc_df)); w = 0.35
ax.bar(x-w/2, tc_df["M1_beta_braak"], w, label="Braak only (M1)",   color="steelblue", alpha=0.8)
ax.bar(x+w/2, tc_df["M2_beta_braak"], w, label="Braak + AZ (M2)",   color="firebrick", alpha=0.8)
ax.axhline(0, color="black", linewidth=0.5)
ax.set_xticks(x); ax.set_xticklabels(tc_df["gene"], fontsize=10)
ax.set_ylabel("β (Braak term)")
ax.set_title("Test C: Braak effect before/after\nAZ score conditioning", fontsize=10)
ax.legend(fontsize=8)

# Row 2 right: AZ score ~ Braak
ax = axes[1, 2]
xv = all_df["braaksc"].values
yv = all_df["AZscore"].values
ax.scatter(xv + np.random.normal(0,0.1,len(xv)), yv, s=5, alpha=0.3, color="grey")
m,b = np.polyfit(xv, yv, 1)
xl = np.linspace(0, 6, 100)
ax.plot(xl, m*xl+b, color="firebrick", linewidth=2)
r, p = stats.pearsonr(xv, yv)
ax.set_xlabel("Braak stage"); ax.set_ylabel("AZ module score")
ax.set_title(f"AZscore ~ Braak\nr={r:.3f}, p={p:.1e}\nβ={az_braak['beta']:.4f}", fontsize=10)

plt.tight_layout()
plt.savefig("tmed_az_plots.png", dpi=200, bbox_inches="tight")
print("Saved: tmed_az_plots.png")
print("\nDone.")
