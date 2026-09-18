"""
CHRONOS - APP/APLP Presynaptic Localization — Test C
======================================================
Does the APP/AZ-module relationship survive AD pathology adjustment?
And does APP/APLP expression itself change with AD pathology?

Models:
  1. APP/APLP1/APLP2 ~ Braak + age + sex + PMI
  2. APP/APLP1/APLP2 ~ AZscore + Braak + age + sex + PMI
  3. AZscore ~ Braak + age + sex + PMI
  4. APP/APLP1/APLP2 ~ AZscore + age + sex + PMI (no pathology — baseline)

Input : ROSMAP_DLPFC_logCPM.tsv, ROSMAP_clinical.csv, ROSMAP_biospecimen_metadata.csv
Output: app_aplp_testC_results.csv, app_aplp_testC.png
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

TARGET_GENES = {
    "APP":   "ENSG00000142192",
    "APLP1": "ENSG00000105290",
    "APLP2": "ENSG00000084234",
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

# ── 1. LOAD ────────────────────────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")

# fix age
clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+", "90")
clinical["age_at_visit_max"] = pd.to_numeric(clinical["age_at_visit_max"], errors="coerce")

# ── 2. LINK SAMPLES TO CLINICAL ────────────────────────────────────────────────

print("Linking samples to clinical...")
sample_ids = expr.columns.tolist()
rna = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()

clin_cols = ["individualID", "braaksc", "cogdx", "ceradsc",
             "apoe_genotype", "msex", "age_at_visit_max", "pmi"]
for col in ["braaksc","cogdx","ceradsc","msex","pmi"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

meta = rna.merge(clinical[clin_cols], on="individualID", how="left")
meta = meta.set_index("specimenID")
print(f"  Samples with full metadata: {meta.dropna().shape[0]} / {len(meta)}")

# ── 3. EXTRACT GENES & COMPUTE AZ SCORE ───────────────────────────────────────

def get_expr(gene_dict, expr):
    rows = {name: expr.loc[eid] for name, eid in gene_dict.items() if eid in expr.index}
    return pd.DataFrame(rows)  # samples x genes

target_df = get_expr(TARGET_GENES, expr)
az_df     = get_expr(AZ_MODULE, expr)

az_z = (az_df - az_df.mean()) / (az_df.std() + 1e-10)
az_score = az_z.mean(axis=1)
az_score.name = "AZscore"

# combine everything into one dataframe
all_df = target_df.join(az_score).join(meta)
all_df = all_df.dropna(subset=["braaksc","msex","age_at_visit_max","pmi"])
print(f"  Samples for regression: {len(all_df)}")

# ── 4. OLS HELPER ─────────────────────────────────────────────────────────────

def ols(y, X_df):
    """Simple OLS. Returns dict with coeff, se, t, p for each predictor."""
    X = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y.values
    coeffs, _, _, _ = np.linalg.lstsq(X, y_v, rcond=None)
    resid = y_v - X @ coeffs
    n, p = len(y_v), X.shape[1]
    mse = np.sum(resid**2) / (n - p)
    cov = mse * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = coeffs / se
    pv = 2 * stats.t.sf(np.abs(t), df=n-p)
    cols = ["intercept"] + list(X_df.columns)
    return pd.DataFrame({"coef": coeffs, "se": se, "t": t, "p": pv}, index=cols)

# ── 5. RUN MODELS ─────────────────────────────────────────────────────────────

base_covs = ["age_at_visit_max", "msex", "pmi"]
results = []

print("\n=== MODEL 1: Gene ~ Braak + covariates ===")
for gene in TARGET_GENES:
    X = all_df[base_covs + ["braaksc"]]
    res = ols(all_df[gene], X)
    row = res.loc["braaksc"]
    print(f"  {gene}: beta={row['coef']:.4f}, p={row['p']:.3e}")
    results.append({"model": "M1_braak_only", "gene": gene,
                    "term": "braaksc", **row.to_dict()})

print("\n=== MODEL 2: Gene ~ AZscore + Braak + covariates ===")
for gene in TARGET_GENES:
    X = all_df[base_covs + ["AZscore", "braaksc"]]
    res = ols(all_df[gene], X)
    for term in ["AZscore", "braaksc"]:
        row = res.loc[term]
        print(f"  {gene} | {term}: beta={row['coef']:.4f}, p={row['p']:.3e}")
        results.append({"model": "M2_az_plus_braak", "gene": gene,
                        "term": term, **row.to_dict()})

print("\n=== MODEL 3: AZscore ~ Braak + covariates ===")
X = all_df[base_covs + ["braaksc"]]
res = ols(all_df["AZscore"], X)
row = res.loc["braaksc"]
print(f"  AZscore: beta={row['coef']:.4f}, p={row['p']:.3e}")
results.append({"model": "M3_az_braak", "gene": "AZscore",
                "term": "braaksc", **row.to_dict()})

print("\n=== MODEL 4: Gene ~ AZscore only (baseline, no pathology) ===")
for gene in TARGET_GENES:
    X = all_df[base_covs + ["AZscore"]]
    res = ols(all_df[gene], X)
    row = res.loc["AZscore"]
    print(f"  {gene}: beta={row['coef']:.4f}, p={row['p']:.3e}")
    results.append({"model": "M4_az_only", "gene": gene,
                    "term": "AZscore", **row.to_dict()})

res_df = pd.DataFrame(results)
res_df.to_csv("app_aplp_testC_results.csv", index=False)

# ── 6. PLOT ────────────────────────────────────────────────────────────────────

print("\nPlotting...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Test C: APP/APLP — AD Pathology & AZ Module Regression", fontsize=12)

genes = list(TARGET_GENES.keys())
colors = {"APP": "firebrick", "APLP1": "steelblue", "APLP2": "seagreen"}

# Plot 1: Gene ~ Braak (M1 beta)
ax = axes[0]
m1 = res_df[res_df.model == "M1_braak_only"]
betas = [m1[m1.gene==g]["coef"].values[0] for g in genes]
pvals = [m1[m1.gene==g]["p"].values[0] for g in genes]
bars = ax.bar(genes, betas, color=[colors[g] for g in genes], alpha=0.8)
for bar, p in zip(bars, pvals):
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.001,
            sig, ha="center", va="bottom", fontsize=11)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_ylabel("Beta (per Braak unit)", fontsize=10)
ax.set_title("M1: Gene ~ Braak\n(+ age, sex, PMI)", fontsize=10)

# Plot 2: AZscore beta in M2 (does AZ association survive pathology adjustment?)
ax = axes[1]
m2_az = res_df[(res_df.model=="M2_az_plus_braak") & (res_df.term=="AZscore")]
m4_az = res_df[(res_df.model=="M4_az_only") & (res_df.term=="AZscore")]
x = np.arange(len(genes))
w = 0.35
b1 = [m4_az[m4_az.gene==g]["coef"].values[0] for g in genes]
b2 = [m2_az[m2_az.gene==g]["coef"].values[0] for g in genes]
ax.bar(x-w/2, b1, w, label="M4: AZ only",       color="steelblue", alpha=0.8)
ax.bar(x+w/2, b2, w, label="M2: AZ + Braak",    color="firebrick", alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(genes)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_ylabel("Beta (AZscore term)", fontsize=10)
ax.set_title("AZ→Gene: before vs after\nBraak adjustment", fontsize=10)
ax.legend(fontsize=8)

# Plot 3: AZscore ~ Braak scatter
ax = axes[2]
x_v = all_df["braaksc"].values
y_v = all_df["AZscore"].values
ax.scatter(x_v, y_v, s=6, alpha=0.4, color="grey")
m, b = np.polyfit(x_v, y_v, 1)
xl = np.linspace(x_v.min(), x_v.max(), 100)
ax.plot(xl, m*xl+b, color="firebrick", linewidth=1.5)
r, p = stats.pearsonr(x_v, y_v)
ax.set_xlabel("Braak stage", fontsize=10)
ax.set_ylabel("AZ module score", fontsize=10)
ax.set_title(f"M3: AZscore ~ Braak\nr={r:.3f}, p={p:.1e}", fontsize=10)

plt.tight_layout()
plt.savefig("app_aplp_testC.png", dpi=200, bbox_inches="tight")
print("Saved: app_aplp_testC.png")
print("\nDone.")
