"""
CHRONOS - Hypothesis 6: APOE4 Stratification
=============================================
Tests A-I: Does APOE4 status associate with TMED2/TMED10/TMED9/AZ module?

Primary:
  A - TMED2 ~ APOE4 + covariates
  B - TMED10 ~ APOE4 + covariates
  C - TMED9 ~ APOE4 + covariates
  D - AZscore ~ APOE4 + covariates

Secondary mechanistic:
  E - log(TMED2/TMED10) ~ APOE4
  F - log(TMED9/TMED2) ~ APOE4
  G - log(TMED9/TMED10) ~ APOE4

Interaction:
  H - APOE4 x Braak for all four primary outcomes

Pre-symptomatic:
  I - Repeat A-D in Braak<=2 (primary + strict NCI)

Input : ROSMAP_DLPFC_logCPM.tsv, ROSMAP_clinical.csv,
        ROSMAP_biospecimen_metadata.csv, ROSMAP_assay_rnaSeq_metadata.csv
Output: h6_primary_results.csv, h6_interactions.csv,
        h6_presymptomatic.csv, h6_plots.png
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
AZ_MODULE = {
    "CASK":   "ENSG00000147044", "STX1A":  "ENSG00000106089",
    "SNAP25": "ENSG00000132639", "VAMP2":  "ENSG00000220205",
    "SYN1":   "ENSG00000008056", "SYN2":   "ENSG00000157542",
    "RIMS1":  "ENSG00000079841", "UNC13A": "ENSG00000130477",
    "DNM1":   "ENSG00000106976", "SV2A":   "ENSG00000197912",
    "SYP":    "ENSG00000102003", "SYT1":   "ENSG00000067715",
}
NEURON_MARKERS = ["ENSG00000102003","ENSG00000067715",
                  "ENSG00000132639","ENSG00000008056","ENSG00000157542"]
ASTRO_MARKERS  = ["ENSG00000131095","ENSG00000171885"]
MICRO_MARKERS  = ["ENSG00000197249","ENSG00000101439"]
OLIGO_MARKERS  = ["ENSG00000197971","ENSG00000123560"]

COLORS = {"TMED2":"steelblue","TMED10":"seagreen","TMED9":"firebrick","AZscore":"purple"}

# ── 1. LOAD ────────────────────────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+","90")
for col in ["age_at_visit_max","braaksc","msex","pmi","cogdx","apoe_genotype"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

# APOE4 carrier status
def apoe4_carrier(g):
    if g in [34,44]: return 1
    elif g in [22,23,33]: return 0
    return np.nan

clinical["apoe4"] = clinical["apoe_genotype"].apply(apoe4_carrier)

sample_ids = expr.columns.tolist()
rna_bio   = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()
rna_batch = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi","cogdx","apoe4"]
meta = rna_bio.merge(clinical[clin_cols], on="individualID", how="left")\
              .merge(rna_batch, on="specimenID", how="left")\
              .set_index("specimenID")

# ── 2. EXTRACT GENES + SCORES ──────────────────────────────────────────────────

tmed_df = pd.DataFrame({n: expr.loc[e] for n,e in TMED_GENES.items() if e in expr.index})

az_df = pd.DataFrame({n: expr.loc[e] for n,e in AZ_MODULE.items() if e in expr.index})
az_z  = (az_df - az_df.mean()) / (az_df.std() + 1e-10)
az_score = az_z.mean(axis=1); az_score.name = "AZscore"

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

df = tmed_df.join(az_score).join(meta).join(neuron_s).join(astro_s).join(micro_s).join(oligo_s)
df = df.dropna(subset=["apoe4","braaksc","age_at_visit_max","msex","pmi"]).copy()

batch_dum  = pd.get_dummies(df["rnaBatch"], prefix="batch", drop_first=True)
df = pd.concat([df, batch_dum], axis=1)
batch_cols = list(batch_dum.columns)
cell_cols  = [c for c in ["neuron_score","astro_score","micro_score","oligo_score"] if c in df.columns]
base_covs  = ["age_at_visit_max","msex","pmi","braaksc"] + batch_cols
full_covs  = base_covs + cell_cols

print(f"  Total samples with APOE4 status: {len(df)}")
print(f"  APOE4+: {int(df['apoe4'].sum())}  APOE4-: {int((df['apoe4']==0).sum())}")

# ── TEST A: COHORT SUMMARY ─────────────────────────────────────────────────────

print("\n=== TEST A: Cohort summary ===")
for grp, label in [(1,"APOE4+"),(0,"APOE4-")]:
    sub = df[df["apoe4"]==grp]
    print(f"  {label}: n={len(sub)}, age={sub['age_at_visit_max'].mean():.1f}±{sub['age_at_visit_max'].std():.1f}, "
          f"sex(M%)={sub['msex'].mean():.0%}, Braak={sub['braaksc'].mean():.2f}")

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

# ── TESTS B/C/D/E_primary: PRIMARY OUTCOMES ~ APOE4 ───────────────────────────

print("\n=== PRIMARY TESTS: Outcome ~ APOE4 (full cohort) ===")
outcomes = ["TMED2","TMED10","TMED9","AZscore"]
primary_results = []

for outcome in outcomes:
    if outcome not in df.columns: continue
    r_raw  = ols_term(df[outcome], df[["apoe4"]], "apoe4")
    r_adj  = ols_term(df[outcome], df[["apoe4"]+base_covs], "apoe4")
    r_full = ols_term(df[outcome], df[["apoe4"]+full_covs], "apoe4")
    print(f"  {outcome}:")
    print(f"    Raw  : β={r_raw['beta']:.4f}, p={r_raw['p']:.3e}")
    print(f"    Adj  : β={r_adj['beta']:.4f}, p={r_adj['p']:.3e}")
    print(f"    Full : β={r_full['beta']:.4f}, p={r_full['p']:.3e}, n={r_full['n']}")
    primary_results.append({"outcome":outcome,
        "raw_beta":r_raw["beta"],"raw_p":r_raw["p"],
        "adj_beta":r_adj["beta"],"adj_p":r_adj["p"],
        "full_beta":r_full["beta"],"full_p":r_full["p"],
        "n":r_full["n"]})

# FDR on adjusted p-values
adj_ps = [r["adj_p"] for r in primary_results]
_, fdrs, _, _ = multipletests(adj_ps, method="fdr_bh")
for r, fdr in zip(primary_results, fdrs):
    r["fdr"] = fdr
    print(f"  {r['outcome']} FDR: {fdr:.3f}")

# ── TESTS E/F/G: RATIOS ~ APOE4 ───────────────────────────────────────────────

print("\n=== SECONDARY: Log-ratios ~ APOE4 ===")
df["ratio_T2_T10"]  = df["TMED2"]  - df["TMED10"]
df["ratio_T9_T2"]   = df["TMED9"]  - df["TMED2"]
df["ratio_T9_T10"]  = df["TMED9"]  - df["TMED10"]

ratio_results = []
for ratio, label in [("ratio_T2_T10","log(TMED2/TMED10)"),
                      ("ratio_T9_T2", "log(TMED9/TMED2)"),
                      ("ratio_T9_T10","log(TMED9/TMED10)")]:
    r = ols_term(df[ratio], df[["apoe4"]+full_covs], "apoe4")
    print(f"  {label}: β={r['beta']:.4f}, p={r['p']:.3e}")
    ratio_results.append({"ratio":label,**r})

# ── TEST H: APOE4 × BRAAK INTERACTION ─────────────────────────────────────────

print("\n=== TEST H: APOE4 × Braak interaction ===")
df["apoe4_x_braak"] = df["apoe4"] * df["braaksc"]
interaction_results = []

for outcome in outcomes:
    if outcome not in df.columns: continue
    r_int = ols_term(df[outcome],
                     df[["apoe4","braaksc","apoe4_x_braak"]+
                        ["age_at_visit_max","msex","pmi"]+batch_cols+cell_cols],
                     "apoe4_x_braak")
    r_main = ols_term(df[outcome],
                      df[["apoe4","braaksc","apoe4_x_braak"]+
                         ["age_at_visit_max","msex","pmi"]+batch_cols+cell_cols],
                      "apoe4")
    print(f"  {outcome}: interaction β={r_int['beta']:.4f}, p={r_int['p']:.3e} | main β={r_main['beta']:.4f}")
    interaction_results.append({"outcome":outcome,
        "main_beta":r_main["beta"],"main_p":r_main["p"],
        "int_beta":r_int["beta"],"int_p":r_int["p"],"n":r_int["n"]})

# ── TEST I: PRE-SYMPTOMATIC (BRAAK<=2) ────────────────────────────────────────

print("\n=== TEST I: Pre-symptomatic (Braak<=2) ===")
df_pre    = df[df["braaksc"]<=2].copy()
df_strict = df[(df["braaksc"]<=2) & (df["cogdx"]==1)].copy()
presymp_results = []

for cohort_name, df_c in [("Braak<=2", df_pre), ("Braak<=2+NCI", df_strict)]:
    print(f"\n  {cohort_name} (n={len(df_c)}, APOE4+={int(df_c['apoe4'].sum())}):")
    for outcome in outcomes:
        if outcome not in df_c.columns: continue
        r = ols_term(df_c[outcome], df_c[["apoe4"]+full_covs], "apoe4")
        print(f"    {outcome}: β={r['beta']:.4f}, p={r['p']:.3e}, n={r['n']}")
        presymp_results.append({"cohort":cohort_name,"outcome":outcome,
            "beta":r["beta"],"p":r["p"],"n":r["n"]})

# ── SAVE ───────────────────────────────────────────────────────────────────────

pd.DataFrame(primary_results).to_csv("h6_primary_results.csv", index=False)
pd.DataFrame(interaction_results).to_csv("h6_interactions.csv", index=False)
pd.DataFrame(presymp_results).to_csv("h6_presymptomatic.csv", index=False)
pd.DataFrame(ratio_results).to_csv("h6_ratios.csv", index=False)
print("\nSaved CSVs.")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("CHRONOS — APOE4 Stratification (Hypothesis 6)", fontsize=13)

# Plot 1: Primary betas (adj vs full)
ax = axes[0,0]
pr_df = pd.DataFrame(primary_results)
x = np.arange(len(pr_df)); w=0.35
ax.bar(x-w/2, pr_df["adj_beta"],  w, label="Adjusted",   color="steelblue", alpha=0.8)
ax.bar(x+w/2, pr_df["full_beta"], w, label="+Cell comp", color="firebrick", alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(pr_df["outcome"], fontsize=10)
ax.axhline(0,color="black",linewidth=0.8)
# significance markers
for i, row in pr_df.iterrows():
    for xi, pv in [(i-w/2, row["adj_p"]), (i+w/2, row["full_p"])]:
        sig = "***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else ""
        if sig:
            ypos = max(row["adj_beta"], row["full_beta"]) + 0.005
            ax.text(xi, ypos, sig, ha="center", fontsize=10)
ax.set_ylabel("β (APOE4 term)"); ax.set_title("Primary outcomes ~ APOE4\nadj vs +cell comp", fontsize=10)
ax.legend(fontsize=8)

# Plot 2: Boxplot of AZscore by APOE4
ax = axes[0,1]
pos = df[df["apoe4"]==0]["AZscore"].dropna()
car = df[df["apoe4"]==1]["AZscore"].dropna()
ax.boxplot([pos.values, car.values], labels=["APOE4−","APOE4+"],
           patch_artist=True,
           boxprops=dict(facecolor="steelblue", alpha=0.6))
t,p = stats.ttest_ind(pos, car)
ax.set_ylabel("AZ module score")
ax.set_title(f"AZscore by APOE4\nt-test p={p:.3e}", fontsize=10)

# Plot 3: Ratios
ax = axes[0,2]
rd = pd.DataFrame(ratio_results)
colors_r = ["steelblue","firebrick","seagreen"]
bars = ax.bar(range(len(rd)), rd["beta"], color=colors_r, alpha=0.8)
ax.set_xticks(range(len(rd))); ax.set_xticklabels(rd["ratio"], rotation=20, fontsize=8)
ax.axhline(0,color="black",linewidth=0.8)
for bar, pv in zip(bars, rd["p"]):
    sig = "***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else ""
    if sig: ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.001, sig, ha="center")
ax.set_ylabel("β (APOE4)"); ax.set_title("Log-ratios ~ APOE4\n(full model)", fontsize=10)

# Plot 4: Interactions
ax = axes[1,0]
ir_df = pd.DataFrame(interaction_results)
ax.bar(range(len(ir_df)), ir_df["int_beta"], color="purple", alpha=0.8)
ax.set_xticks(range(len(ir_df))); ax.set_xticklabels(ir_df["outcome"])
ax.axhline(0,color="black",linewidth=0.8)
ax.set_ylabel("β (APOE4×Braak)")
ax.set_title("Test H: APOE4×Braak interaction", fontsize=10)

# Plot 5: Pre-symptomatic
ax = axes[1,1]
ps_df = pd.DataFrame(presymp_results)
x = np.arange(len(outcomes)); w=0.35
b_pre = [ps_df[(ps_df.cohort=="Braak<=2")&(ps_df.outcome==o)]["beta"].values[0]
         if len(ps_df[(ps_df.cohort=="Braak<=2")&(ps_df.outcome==o)])>0 else np.nan
         for o in outcomes]
b_str = [ps_df[(ps_df.cohort=="Braak<=2+NCI")&(ps_df.outcome==o)]["beta"].values[0]
         if len(ps_df[(ps_df.cohort=="Braak<=2+NCI")&(ps_df.outcome==o)])>0 else np.nan
         for o in outcomes]
ax.bar(x-w/2, b_pre, w, label="Braak≤2",     color="steelblue", alpha=0.8)
ax.bar(x+w/2, b_str, w, label="Braak≤2+NCI", color="firebrick", alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(outcomes)
ax.axhline(0,color="black",linewidth=0.8)
ax.set_ylabel("β (APOE4)"); ax.set_title("Test I: Pre-symptomatic\nAPOE4 effect", fontsize=10)
ax.legend(fontsize=8)

# Plot 6: Braak distribution by APOE4 group
ax = axes[1,2]
for grp, label, color in [(0,"APOE4−","steelblue"),(1,"APOE4+","firebrick")]:
    sub = df[df["apoe4"]==grp]["braaksc"].dropna()
    ax.hist(sub, bins=range(8), alpha=0.6, color=color, label=f"{label} (n={len(sub)})")
ax.set_xlabel("Braak stage"); ax.set_ylabel("Count")
ax.set_title("Braak distribution by APOE4", fontsize=10)
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("h6_plots.png", dpi=200, bbox_inches="tight")
print("Saved: h6_plots.png")
print("\nDone.")
