"""
CHRONOS - Hypothesis 5: Age-dependent decline in pre-symptomatic window
========================================================================
Question: Do TMED9/TMED10 decline with age in people with minimal AD pathology?
This tests whether normal aging drives p24 system decline before AD develops.

Tests A-I as specified:
  A - Define cohorts (Braak<=2, strict NCI)
  B - TMED9 ~ age (primary + adjusted)
  C - TMED10 ~ age
  D - TMED2 ~ age (control gene)
  E - Full pathology adjustment (braak + ceradsc + cell composition)
  F - Strict NCI-only sensitivity analysis
  G - Nested models showing cell composition effect on age slope
  H - Nonlinear age trajectory
  I - Joint TMED9/TMED10 aging signature

Input : ROSMAP_DLPFC_logCPM.tsv, ROSMAP_clinical.csv,
        ROSMAP_biospecimen_metadata.csv, ROSMAP_assay_rnaSeq_metadata.csv
Output: h5_results.csv, h5_plots.png, h5_nested_models.csv
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
    "TMED9":  "ENSG00000184840",
    "TMED10": "ENSG00000170348",
    "TMED2":  "ENSG00000086598",
}
COLORS = {"TMED9":"firebrick","TMED10":"seagreen","TMED2":"steelblue"}

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
for col in ["age_at_visit_max","braaksc","msex","pmi","cogdx","ceradsc"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

sample_ids = expr.columns.tolist()
rna_bio = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()
rna_batch = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
clin_cols = ["individualID","braaksc","msex","age_at_visit_max","pmi","cogdx","ceradsc"]
meta = rna_bio.merge(clinical[clin_cols], on="individualID", how="left")\
              .merge(rna_batch, on="specimenID", how="left")\
              .set_index("specimenID")

# TMED genes
tmed_df = pd.DataFrame({n: expr.loc[e] for n,e in TMED_GENES.items() if e in expr.index})

# cell type scores
def ct_score(markers, name):
    present = [m for m in markers if m in expr.index]
    if not present: return None
    vals = expr.loc[present]
    z = (vals - vals.mean(axis=1).values[:,None]) / (vals.std(axis=1).values[:,None] + 1e-10)
    s = z.mean(axis=0); s.name = name; return s

neuron_s = ct_score(NEURON_MARKERS, "neuron_score")
astro_s  = ct_score(ASTRO_MARKERS,  "astro_score")
micro_s  = ct_score(MICRO_MARKERS,  "micro_score")
oligo_s  = ct_score(OLIGO_MARKERS,  "oligo_score")

df_full = tmed_df.join(meta).join(neuron_s).join(astro_s).join(micro_s).join(oligo_s)
df_full = df_full.dropna(subset=["braaksc","age_at_visit_max","msex","pmi"]).copy()
df_full["braaksc"] = df_full["braaksc"].astype(float)

batch_dum = pd.get_dummies(df_full["rnaBatch"], prefix="batch", drop_first=True)
df_full = pd.concat([df_full, batch_dum], axis=1)
batch_cols = list(batch_dum.columns)
cell_cols  = [c for c in ["neuron_score","astro_score","micro_score","oligo_score"] if c in df_full.columns]

print(f"  Full cohort: {len(df_full)} samples")

# ── TEST A: DEFINE COHORTS ─────────────────────────────────────────────────────

print("\n=== TEST A: Cohort definitions ===")

# Primary: Braak 0-2
df_primary = df_full[df_full["braaksc"] <= 2].copy()

# Strict: Braak 0-2 + NCI (cogdx == 1)
df_strict = df_full[(df_full["braaksc"] <= 2) & (df_full["cogdx"] == 1)].copy()

for name, df in [("Primary (Braak<=2)", df_primary), ("Strict (Braak<=2 + NCI)", df_strict)]:
    print(f"\n  {name}: n={len(df)}")
    print(f"    Braak dist: {dict(df['braaksc'].value_counts().sort_index())}")
    print(f"    cogdx dist: {dict(df['cogdx'].value_counts().sort_index())}")
    print(f"    Age: mean={df['age_at_visit_max'].mean():.1f}, "
          f"SD={df['age_at_visit_max'].std():.1f}, "
          f"range={df['age_at_visit_max'].min():.0f}-{df['age_at_visit_max'].max():.0f}")

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

# ── TESTS B/C/D: GENE ~ AGE IN PRIMARY COHORT ─────────────────────────────────

print("\n=== TESTS B/C/D: Gene ~ Age (primary cohort, Braak<=2) ===")
base_covs   = ["msex","pmi"] + batch_cols
path_covs   = ["braaksc","ceradsc"]

all_results = []
for gene in ["TMED9","TMED10","TMED2"]:
    if gene not in df_primary.columns: continue

    r_raw  = ols_term(df_primary[gene], df_primary[["age_at_visit_max"]], "age_at_visit_max")
    r_adj  = ols_term(df_primary[gene], df_primary[["age_at_visit_max"]+base_covs], "age_at_visit_max")
    r_path = ols_term(df_primary[gene], df_primary[["age_at_visit_max"]+base_covs+path_covs].dropna(subset=path_covs, how="any")
                      if True else df_primary[["age_at_visit_max"]+base_covs], "age_at_visit_max")

    # with pathology — handle NaNs in ceradsc
    df_p = df_primary.dropna(subset=["ceradsc"])
    r_path = ols_term(df_p[gene], df_p[["age_at_visit_max"]+base_covs+path_covs], "age_at_visit_max")
    r_full = ols_term(df_p[gene], df_p[["age_at_visit_max"]+base_covs+path_covs+cell_cols], "age_at_visit_max")

    print(f"\n  {gene} (n_primary={len(df_primary)}):")
    print(f"    Raw           : β={r_raw['beta']:.4f} [{r_raw['ci_lo']:.4f},{r_raw['ci_hi']:.4f}], p={r_raw['p']:.3e}")
    print(f"    + sex/PMI/batch: β={r_adj['beta']:.4f}, p={r_adj['p']:.3e}")
    print(f"    + pathology   : β={r_path['beta']:.4f}, p={r_path['p']:.3e}, n={r_path['n']}")
    print(f"    + cells       : β={r_full['beta']:.4f}, p={r_full['p']:.3e}")

    all_results.append({"cohort":"primary","gene":gene,
        "raw_beta":r_raw["beta"],"raw_p":r_raw["p"],
        "adj_beta":r_adj["beta"],"adj_p":r_adj["p"],
        "path_beta":r_path["beta"],"path_p":r_path["p"],
        "full_beta":r_full["beta"],"full_p":r_full["p"],
        "n_raw":r_raw["n"],"n_full":r_full["n"]})

# FDR on primary adjusted
primary_ps = [r["adj_p"] for r in all_results if r["cohort"]=="primary"]
_, fdrs, _, _ = multipletests(primary_ps, method="fdr_bh")
for r, fdr in zip([r for r in all_results if r["cohort"]=="primary"], fdrs):
    r["fdr"] = fdr
    print(f"  {r['gene']} FDR (adj): {fdr:.3f}")

# ── TEST F: STRICT NCI-ONLY ────────────────────────────────────────────────────

print("\n=== TEST F: Strict NCI-only (Braak<=2 + cogdx==1) ===")
df_s = df_strict.dropna(subset=["ceradsc"])
for gene in ["TMED9","TMED10","TMED2"]:
    if gene not in df_s.columns: continue
    r = ols_term(df_s[gene], df_s[["age_at_visit_max"]+base_covs+path_covs+cell_cols], "age_at_visit_max")
    print(f"  {gene}: β={r['beta']:.4f}, p={r['p']:.3e}, n={r['n']}")
    all_results.append({"cohort":"strict_NCI","gene":gene,
        "full_beta":r["beta"],"full_p":r["p"],"n_full":r["n"],
        "raw_beta":np.nan,"raw_p":np.nan,"adj_beta":np.nan,"adj_p":np.nan,
        "path_beta":np.nan,"path_p":np.nan,"fdr":np.nan})

# ── TEST G: NESTED MODELS ──────────────────────────────────────────────────────

print("\n=== TEST G: Nested models — cell composition effect on β_age ===")
nested = []
df_g = df_primary.dropna(subset=["ceradsc"]).copy()
for gene in ["TMED9","TMED10","TMED2"]:
    if gene not in df_g.columns: continue
    m1 = ols_term(df_g[gene], df_g[["age_at_visit_max"]], "age_at_visit_max")
    m2 = ols_term(df_g[gene], df_g[["age_at_visit_max"]+base_covs+path_covs], "age_at_visit_max")
    m3 = ols_term(df_g[gene], df_g[["age_at_visit_max"]+base_covs+path_covs+cell_cols], "age_at_visit_max")
    print(f"  {gene}:")
    print(f"    M1 (age only)       : β={m1['beta']:.4f}, p={m1['p']:.3e}")
    print(f"    M2 (+pathology)     : β={m2['beta']:.4f}, p={m2['p']:.3e}")
    print(f"    M3 (+cells)         : β={m3['beta']:.4f}, p={m3['p']:.3e}")
    nested.append({"gene":gene,
        "M1_beta":m1["beta"],"M1_p":m1["p"],
        "M2_beta":m2["beta"],"M2_p":m2["p"],
        "M3_beta":m3["beta"],"M3_p":m3["p"]})

# ── TEST H: NONLINEAR ─────────────────────────────────────────────────────────

print("\n=== TEST H: Nonlinear age trajectory ===")
df_h = df_primary.dropna(subset=["ceradsc"]).copy()
df_h["age2"] = df_h["age_at_visit_max"] ** 2
for gene in ["TMED9","TMED10","TMED2"]:
    if gene not in df_h.columns: continue
    r_lin  = ols_term(df_h[gene], df_h[["age_at_visit_max"]+base_covs+path_covs], "age_at_visit_max")
    r_quad = ols_term(df_h[gene], df_h[["age_at_visit_max","age2"]+base_covs+path_covs], "age2")

    # F-test
    cov_lin = ["age_at_visit_max"]+base_covs+path_covs
    cov_qua = ["age_at_visit_max","age2"]+base_covs+path_covs
    mask = df_h[cov_qua+[gene]].dropna().index
    y_v = df_h.loc[mask,gene].values
    Xl = np.column_stack([np.ones(len(mask))]+[df_h.loc[mask,c].values for c in cov_lin])
    Xq = np.column_stack([np.ones(len(mask))]+[df_h.loc[mask,c].values for c in cov_qua])
    cl,_,_,_ = np.linalg.lstsq(Xl,y_v,rcond=None)
    cq,_,_,_ = np.linalg.lstsq(Xq,y_v,rcond=None)
    ssl = np.sum((y_v-Xl@cl)**2); ssq = np.sum((y_v-Xq@cq)**2)
    n,p = len(y_v),Xq.shape[1]
    F = ((ssl-ssq)/1)/(ssq/(n-p))
    pF = stats.f.sf(F,1,n-p)
    print(f"  {gene}: linear β={r_lin['beta']:.4f} | quadratic term β={r_quad['beta']:.6f}, F-test p={pF:.3e}")

# ── TEST I: JOINT AGING SIGNATURE ─────────────────────────────────────────────

print("\n=== TEST I: Joint TMED9/TMED10 aging signature ===")
for cohort_name, df_c in [("Primary", df_primary), ("Strict NCI", df_strict)]:
    df_c2 = df_c.dropna(subset=["ceradsc"]).copy()
    if len(df_c2) < 20: continue
    for g in ["TMED9","TMED10"]:
        if g in df_c2.columns:
            df_c2[f"z_{g}"] = (df_c2[g]-df_c2[g].mean())/df_c2[g].std()
    if "z_TMED9" in df_c2.columns and "z_TMED10" in df_c2.columns:
        df_c2["aging_p24"] = (df_c2["z_TMED9"] + df_c2["z_TMED10"]) / 2
        r = ols_term(df_c2["aging_p24"],
                     df_c2[["age_at_visit_max"]+base_covs+path_covs+cell_cols],
                     "age_at_visit_max")
        print(f"  {cohort_name}: aging_p24 ~ age: β={r['beta']:.4f}, p={r['p']:.3e}, n={r['n']}")

# ── SAVE ───────────────────────────────────────────────────────────────────────

pd.DataFrame(all_results).to_csv("h5_results.csv", index=False)
pd.DataFrame(nested).to_csv("h5_nested_models.csv", index=False)
print("\nSaved CSVs.")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("CHRONOS — Age-dependent Decline in Pre-symptomatic Window (Hypothesis 5)", fontsize=12)

genes = ["TMED9","TMED10","TMED2"]

# Row 1: age scatter in primary cohort
for i, gene in enumerate(genes):
    ax = axes[0,i]
    if gene not in df_primary.columns: continue
    x = df_primary["age_at_visit_max"].values
    y = df_primary[gene].values
    mask = ~(np.isnan(x)|np.isnan(y))
    x,y = x[mask],y[mask]
    ax.scatter(x,y,s=8,alpha=0.4,color=COLORS[gene])
    m,b = np.polyfit(x,y,1)
    xl = np.linspace(x.min(),x.max(),100)
    ax.plot(xl,m*xl+b,color="black",linewidth=1.5)
    r = next((res for res in all_results if res["gene"]==gene and res["cohort"]=="primary"),{})
    ax.set_xlabel("Age at visit max"); ax.set_ylabel(f"{gene} (logCPM)")
    ax.set_title(f"{gene} ~ Age (Braak≤2)\nβ={r.get('adj_beta',np.nan):.4f}, p={r.get('adj_p',np.nan):.2e}", fontsize=10)

# Row 2 left: nested model betas (M1, M2, M3)
ax = axes[1,0]
nested_df = pd.DataFrame(nested)
if not nested_df.empty:
    x = np.arange(len(nested_df)); w=0.25
    ax.bar(x-w,   nested_df["M1_beta"], w, label="M1 age only",    color="steelblue",alpha=0.8)
    ax.bar(x,     nested_df["M2_beta"], w, label="M2 +pathology",  color="orange",   alpha=0.8)
    ax.bar(x+w,   nested_df["M3_beta"], w, label="M3 +cells",      color="firebrick",alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(nested_df["gene"])
    ax.axhline(0,color="black",linewidth=0.8)
    ax.set_ylabel("β (age term)"); ax.set_title("Test G: Nested models\nβ_age attenuation", fontsize=10)
    ax.legend(fontsize=8)

# Row 2 middle: primary vs strict NCI comparison
ax = axes[1,1]
primary_betas = {r["gene"]:r["full_beta"] for r in all_results if r["cohort"]=="primary"}
strict_betas  = {r["gene"]:r["full_beta"] for r in all_results if r["cohort"]=="strict_NCI"}
x = np.arange(len(genes)); w=0.35
b1 = [primary_betas.get(g,np.nan) for g in genes]
b2 = [strict_betas.get(g,np.nan)  for g in genes]
ax.bar(x-w/2, b1, w, label="Primary (Braak≤2)",    color="steelblue", alpha=0.8)
ax.bar(x+w/2, b2, w, label="Strict (Braak≤2+NCI)", color="firebrick", alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(genes)
ax.axhline(0,color="black",linewidth=0.8)
ax.set_ylabel("β (age, full model)")
ax.set_title("Primary vs Strict NCI\nfull model β_age", fontsize=10)
ax.legend(fontsize=8)

# Row 2 right: age distribution in cohorts
ax = axes[1,2]
ax.hist(df_primary["age_at_visit_max"].dropna(), bins=15, alpha=0.6,
        color="steelblue", label=f"Primary (n={len(df_primary)})")
ax.hist(df_strict["age_at_visit_max"].dropna(),  bins=15, alpha=0.6,
        color="firebrick", label=f"Strict NCI (n={len(df_strict)})")
ax.set_xlabel("Age at visit max"); ax.set_ylabel("Count")
ax.set_title("Age distribution\nin pre-symptomatic cohorts", fontsize=10)
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("h5_plots.png", dpi=200, bbox_inches="tight")
print("Saved: h5_plots.png")
print("\nDone.")
