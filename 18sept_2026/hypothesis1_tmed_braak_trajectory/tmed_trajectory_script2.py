"""
CHRONOS - Hypothesis 1: TMED2/TMED10/TMED9
Script 2: Tests I through L
  I  - TMED2 ↔ TMED10 relationship + Braak interaction
  J  - TMED9 ↔ TMED2, TMED9 ↔ TMED10 + Braak interactions
  K  - TMED2/TMED10 log-ratio ~ Braak
  L  - Joint p24 score ~ Braak (let data decide TMED9 direction)

Input : lipids_with_apoe4.csv already in folder (uses ROSMAP_DLPFC_logCPM.tsv)
Output: tmed_script2_summary.csv, tmed_script2.png
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

TMED_GENES = {
    "TMED2":  "ENSG00000086598",
    "TMED10": "ENSG00000170348",
    "TMED9":  "ENSG00000184840",
}
COLORS = {"TMED2": "steelblue", "TMED10": "seagreen", "TMED9": "firebrick"}

# ── 1. LOAD (same as script 1) ─────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+","90")
clinical["age_at_visit_max"] = pd.to_numeric(clinical["age_at_visit_max"], errors="coerce")
for col in ["braaksc","msex","pmi"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

sample_ids = expr.columns.tolist()
rna_bio = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()
clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi"]
meta = rna_bio.merge(clinical[clin_cols], on="individualID", how="left")
rna_meta_sub = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
meta = meta.merge(rna_meta_sub, on="specimenID", how="left").set_index("specimenID")

tmed_df = pd.DataFrame({n: expr.loc[e] for n,e in TMED_GENES.items() if e in expr.index})
df = tmed_df.join(meta).dropna(subset=["braaksc","age_at_visit_max","msex","pmi"]).copy()
df["braaksc"] = df["braaksc"].astype(float)

# encode batch
batch_dummies = pd.get_dummies(df["rnaBatch"], prefix="batch", drop_first=True)
df = pd.concat([df, batch_dummies], axis=1)
batch_cols = list(batch_dummies.columns)
base_covs = ["age_at_visit_max","msex","pmi"] + batch_cols

print(f"  Samples: {len(df)}")

# ── OLS HELPER ─────────────────────────────────────────────────────────────────

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
    return dict(beta=beta,se=se_b,ci_lo=beta-1.96*se_b,ci_hi=beta+1.96*se_b,p=pv,n=int(mask.sum()))

# ── TEST I: TMED2 ↔ TMED10 + BRAAK INTERACTION ────────────────────────────────

print("\n=== TEST I: TMED2 ↔ TMED10 relationship ===")

# I-1: TMED10 ~ TMED2 + covariates (baseline relationship)
res_i1 = ols_term(df["TMED10"], df[["TMED2"]+base_covs], "TMED2")
print(f"  TMED10 ~ TMED2: β={res_i1['beta']:.4f}, p={res_i1['p']:.3e}")

# I-2: TMED2 ~ TMED10 + covariates
res_i2 = ols_term(df["TMED2"], df[["TMED10"]+base_covs], "TMED10")
print(f"  TMED2 ~ TMED10: β={res_i2['beta']:.4f}, p={res_i2['p']:.3e}")

# I-3: TMED10 ~ TMED2 + Braak + TMED2*Braak (interaction)
df["TMED2_x_Braak"] = df["TMED2"] * df["braaksc"]
res_i3_main = ols_term(df["TMED10"], df[["TMED2","braaksc","TMED2_x_Braak"]+base_covs], "TMED2")
res_i3_int  = ols_term(df["TMED10"], df[["TMED2","braaksc","TMED2_x_Braak"]+base_covs], "TMED2_x_Braak")
print(f"  Interaction TMED2:Braak → TMED10: β={res_i3_int['beta']:.4f}, p={res_i3_int['p']:.3e}")
print(f"  (Does TMED2→TMED10 coupling change with pathology?)")

testI = [
    {"test":"I1","predictor":"TMED2","outcome":"TMED10",**res_i1},
    {"test":"I2","predictor":"TMED10","outcome":"TMED2",**res_i2},
    {"test":"I3_main","predictor":"TMED2","outcome":"TMED10","term":"TMED2",**res_i3_main},
    {"test":"I3_int","predictor":"TMED2:Braak","outcome":"TMED10","term":"interaction",**res_i3_int},
]

# ── TEST J: TMED9 RELATIONSHIPS ────────────────────────────────────────────────

print("\n=== TEST J: TMED9 relationships ===")

# J-1: TMED9 ~ TMED2
res_j1 = ols_term(df["TMED9"], df[["TMED2"]+base_covs], "TMED2")
print(f"  TMED9 ~ TMED2: β={res_j1['beta']:.4f}, p={res_j1['p']:.3e}")

# J-2: TMED9 ~ TMED10
res_j2 = ols_term(df["TMED9"], df[["TMED10"]+base_covs], "TMED10")
print(f"  TMED9 ~ TMED10: β={res_j2['beta']:.4f}, p={res_j2['p']:.3e}")

# J-3: TMED9 ~ TMED2 + Braak + TMED2*Braak
res_j3_int = ols_term(df["TMED9"], df[["TMED2","braaksc","TMED2_x_Braak"]+base_covs], "TMED2_x_Braak")
print(f"  Interaction TMED2:Braak → TMED9: β={res_j3_int['beta']:.4f}, p={res_j3_int['p']:.3e}")

# J-4: TMED9 ~ TMED10 + Braak + TMED10*Braak
df["TMED10_x_Braak"] = df["TMED10"] * df["braaksc"]
res_j4_int = ols_term(df["TMED9"], df[["TMED10","braaksc","TMED10_x_Braak"]+base_covs], "TMED10_x_Braak")
print(f"  Interaction TMED10:Braak → TMED9: β={res_j4_int['beta']:.4f}, p={res_j4_int['p']:.3e}")

testJ = [
    {"test":"J1","predictor":"TMED2","outcome":"TMED9",**res_j1},
    {"test":"J2","predictor":"TMED10","outcome":"TMED9",**res_j2},
    {"test":"J3_int","predictor":"TMED2:Braak","outcome":"TMED9",**res_j3_int},
    {"test":"J4_int","predictor":"TMED10:Braak","outcome":"TMED9",**res_j4_int},
]

# ── TEST K: LOG-RATIO ──────────────────────────────────────────────────────────

print("\n=== TEST K: TMED2/TMED10 log-ratio ~ Braak ===")
df["ratio_T2_T10"] = df["TMED2"] - df["TMED10"]  # on log scale = log ratio

res_k = ols_term(df["ratio_T2_T10"], df[["braaksc"]+base_covs], "braaksc")
print(f"  log(TMED2/TMED10) ~ Braak: β={res_k['beta']:.4f}, p={res_k['p']:.3e}")
print(f"  (positive β = TMED2 gains relative to TMED10 with pathology; negative = TMED2 loses)")

# also TMED9/TMED2
df["ratio_T9_T2"] = df["TMED9"] - df["TMED2"]
res_k2 = ols_term(df["ratio_T9_T2"], df[["braaksc"]+base_covs], "braaksc")
print(f"  log(TMED9/TMED2) ~ Braak: β={res_k2['beta']:.4f}, p={res_k2['p']:.3e}")

testK = [
    {"test":"K1","ratio":"log(TMED2/TMED10)",**res_k},
    {"test":"K2","ratio":"log(TMED9/TMED2)",**res_k2},
]

# ── TEST L: JOINT p24 SCORE ────────────────────────────────────────────────────

print("\n=== TEST L: Joint p24 score ~ Braak ===")

# z-score each gene
for g in ["TMED2","TMED10","TMED9"]:
    df[f"z_{g}"] = (df[g] - df[g].mean()) / df[g].std()

# check direction from Test B data: TMED2 β<0, TMED10 β>0 (not sig), TMED9 β<0
# Don't force — use actual signs from unadjusted regression
# TMED2: negative, TMED10: positive (slight), TMED9: negative
# Score = mean of z-scores with sign flipped if going opposite direction
# Since hypothesis is TMED2/TMED10 down and TMED9 up:
# but data shows TMED9 also goes down → don't flip TMED9

# Version 1: hypothesis-driven (TMED9 flipped as if it should go up)
df["p24_score_hyp"] = (df["z_TMED2"] + df["z_TMED10"] - df["z_TMED9"]) / 3

# Version 2: data-driven (all same sign direction — all slightly negative)
df["p24_score_data"] = (df["z_TMED2"] + df["z_TMED10"] + df["z_TMED9"]) / 3

res_l1 = ols_term(df["p24_score_hyp"],  df[["braaksc"]+base_covs], "braaksc")
res_l2 = ols_term(df["p24_score_data"], df[["braaksc"]+base_covs], "braaksc")
print(f"  p24_score (hypothesis: TMED9 flipped) ~ Braak: β={res_l1['beta']:.4f}, p={res_l1['p']:.3e}")
print(f"  p24_score (data-driven: all same dir) ~ Braak: β={res_l2['beta']:.4f}, p={res_l2['p']:.3e}")

testL = [
    {"test":"L1_hyp","score":"p24_hypothesis_driven",**res_l1},
    {"test":"L2_data","score":"p24_data_driven",**res_l2},
]

# ── SAVE ───────────────────────────────────────────────────────────────────────

all_res = testI + testJ + testK + testL
pd.DataFrame(all_res).to_csv("tmed_script2_summary.csv", index=False)
print("\nSaved: tmed_script2_summary.csv")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("CHRONOS — TMED Inter-gene Relationships & p24 Score (Tests I–L)", fontsize=12)

# I: TMED2 vs TMED10 scatter
ax = axes[0,0]
ax.scatter(df["TMED2"], df["TMED10"], s=5, alpha=0.4, color="steelblue")
m,b = np.polyfit(df["TMED2"].dropna(), df.loc[df["TMED2"].notna(),"TMED10"].dropna(), 1)
xl = np.linspace(df["TMED2"].min(), df["TMED2"].max(), 100)
ax.plot(xl, m*xl+b, color="black", linewidth=1.5)
ax.set_xlabel("TMED2 (logCPM)"); ax.set_ylabel("TMED10 (logCPM)")
ax.set_title(f"Test I: TMED2 ↔ TMED10\nβ={res_i1['beta']:.3f}, p={res_i1['p']:.1e}", fontsize=10)

# J: TMED9 vs TMED2 scatter
ax = axes[0,1]
ax.scatter(df["TMED2"], df["TMED9"], s=5, alpha=0.4, color="firebrick")
m,b = np.polyfit(df["TMED2"].dropna(), df.loc[df["TMED2"].notna(),"TMED9"].dropna(), 1)
ax.plot(xl, m*xl+b, color="black", linewidth=1.5)
ax.set_xlabel("TMED2 (logCPM)"); ax.set_ylabel("TMED9 (logCPM)")
ax.set_title(f"Test J: TMED2 ↔ TMED9\nβ={res_j1['beta']:.3f}, p={res_j1['p']:.1e}", fontsize=10)

# J: TMED9 vs TMED10 scatter
ax = axes[0,2]
ax.scatter(df["TMED10"], df["TMED9"], s=5, alpha=0.4, color="seagreen")
m,b = np.polyfit(df["TMED10"].dropna(), df.loc[df["TMED10"].notna(),"TMED9"].dropna(), 1)
xl2 = np.linspace(df["TMED10"].min(), df["TMED10"].max(), 100)
ax.plot(xl2, m*xl2+b, color="black", linewidth=1.5)
ax.set_xlabel("TMED10 (logCPM)"); ax.set_ylabel("TMED9 (logCPM)")
ax.set_title(f"Test J: TMED10 ↔ TMED9\nβ={res_j2['beta']:.3f}, p={res_j2['p']:.1e}", fontsize=10)

# K: log-ratio across Braak
ax = axes[1,0]
braak_means = df.groupby("braaksc")["ratio_T2_T10"].mean()
braak_ses   = df.groupby("braaksc")["ratio_T2_T10"].sem()
ax.errorbar(braak_means.index, braak_means.values, yerr=braak_ses.values,
            fmt="o-", color="purple", linewidth=2, capsize=4)
ax.axhline(braak_means.values.mean(), color="grey", linestyle="--", linewidth=0.8)
ax.set_xlabel("Braak stage"); ax.set_ylabel("log(TMED2/TMED10)")
ax.set_title(f"Test K: TMED2/TMED10 ratio ~ Braak\nβ={res_k['beta']:.4f}, p={res_k['p']:.3e}", fontsize=10)

# L: p24 score across Braak (both versions)
ax = axes[1,1]
for score, color, label in [
    ("p24_score_hyp",  "steelblue", "Hypothesis-driven"),
    ("p24_score_data", "firebrick", "Data-driven"),
]:
    bm = df.groupby("braaksc")[score].mean()
    bs = df.groupby("braaksc")[score].sem()
    ax.errorbar(bm.index, bm.values, yerr=bs.values,
                fmt="o-", color=color, label=label, linewidth=2, capsize=3)
ax.axhline(0, color="black", linewidth=0.5)
ax.set_xlabel("Braak stage"); ax.set_ylabel("p24 score (z-scored)")
ax.set_title("Test L: p24 composite score ~ Braak", fontsize=10)
ax.legend(fontsize=8)

# Interaction: does TMED2→TMED10 coupling change with Braak?
ax = axes[1,2]
braak_groups = [(0,2,"0-2","steelblue"), (3,4,"3-4","orange"), (5,6,"5-6","firebrick")]
for lo,hi,label,color in braak_groups:
    sub = df[(df["braaksc"]>=lo)&(df["braaksc"]<=hi)]
    ax.scatter(sub["TMED2"], sub["TMED10"], s=5, alpha=0.4, color=color, label=label)
    m,b = np.polyfit(sub["TMED2"].values, sub["TMED10"].values, 1)
    xl_ = np.linspace(sub["TMED2"].min(), sub["TMED2"].max(), 50)
    ax.plot(xl_, m*xl_+b, color=color, linewidth=1.5)
ax.set_xlabel("TMED2"); ax.set_ylabel("TMED10")
ax.set_title(f"Test I: TMED2→TMED10 by Braak group\ninteraction p={res_i3_int['p']:.3e}", fontsize=10)
ax.legend(fontsize=8, title="Braak")

plt.tight_layout()
plt.savefig("tmed_script2.png", dpi=200, bbox_inches="tight")
print("Saved: tmed_script2.png")
print("\nDone.")
