"""
CHRONOS - APP/APLP Presynaptic Localization
============================================
Test A: Are APP/APLP1/APLP2 expressed in bulk DLPFC?
Test B: Do they correlate with the presynaptic active-zone module?
Test B+: Is this specific to the presynaptic module vs a generic neuronal signal?

Input : ROSMAP_DLPFC_logCPM.tsv
        ROSMAP_clinical.csv
        ROSMAP_biospecimen_metadata.csv
Output: app_aplp_bulk_results.csv
        app_aplp_correlations.png
        app_aplp_module_specificity.png
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

# ── GENE SETS ──────────────────────────────────────────────────────────────────

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

# Control module: generic housekeeping / non-neuronal genes
CTRL_MODULE = {
    "ACTB":   "ENSG00000075624",  # actin
    "GAPDH":  "ENSG00000111640",  # housekeeping
    "LDHA":   "ENSG00000134333",  # metabolic
    "VIM":    "ENSG00000026025",  # astrocyte/mesenchymal
    "GFAP":   "ENSG00000131095",  # astrocyte
    "AIF1":   "ENSG00000197249",  # microglia
    "MBP":    "ENSG00000197971",  # oligodendrocyte
    "PLP1":   "ENSG00000123560",  # oligodendrocyte
}

# ── 1. LOAD DATA ───────────────────────────────────────────────────────────────

print("Loading bulk RNA...")
expr = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
print(f"  Matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")

# ── 2. EXTRACT GENES ───────────────────────────────────────────────────────────

def get_genes(gene_dict, expr):
    found = {}
    for name, eid in gene_dict.items():
        if eid in expr.index:
            found[name] = expr.loc[eid]
        else:
            print(f"  WARNING: {name} ({eid}) not found")
    return pd.DataFrame(found)

print("\nExtracting gene sets...")
target_df = get_genes(TARGET_GENES, expr)   # samples x 3
az_df     = get_genes(AZ_MODULE,    expr)   # samples x 12
ctrl_df   = get_genes(CTRL_MODULE,  expr)   # samples x 8

print(f"  Target genes found : {target_df.shape[1]}/3")
print(f"  AZ module genes    : {az_df.shape[1]}/{len(AZ_MODULE)}")
print(f"  Control genes found: {ctrl_df.shape[1]}/{len(CTRL_MODULE)}")

# ── 3. TEST A: EXPRESSION SUMMARY ─────────────────────────────────────────────

print("\n=== TEST A: Expression levels ===")
for gene in TARGET_GENES:
    if gene in target_df.columns:
        vals = target_df[gene].dropna()
        print(f"  {gene}: mean={vals.mean():.2f}, median={vals.median():.2f}, "
              f"min={vals.min():.2f}, max={vals.max():.2f}, n={len(vals)}")

# ── 4. COMPUTE MODULE SCORES ───────────────────────────────────────────────────

# z-score each gene within module, then average
def module_score(df):
    z = (df - df.mean()) / (df.std() + 1e-10)
    return z.mean(axis=1)

az_score   = module_score(az_df)
ctrl_score = module_score(ctrl_df)

print(f"\n  AZ module score: mean={az_score.mean():.3f}, std={az_score.std():.3f}")
print(f"  Ctrl score     : mean={ctrl_score.mean():.3f}, std={ctrl_score.std():.3f}")

# ── 5. TEST B: CORRELATIONS ────────────────────────────────────────────────────

print("\n=== TEST B: APP/APLP correlation with AZ module ===")
results = []
for gene in TARGET_GENES:
    if gene not in target_df.columns:
        continue
    y = target_df[gene]
    common = y.index.intersection(az_score.index)
    y_v = y.loc[common].values
    az_v = az_score.loc[common].values
    ct_v = ctrl_score.loc[common].values

    r_az, p_az = stats.pearsonr(y_v, az_v)
    r_ct, p_ct = stats.pearsonr(y_v, ct_v)

    print(f"  {gene} vs AZ module  : r={r_az:.3f}, p={p_az:.2e}")
    print(f"  {gene} vs Ctrl module: r={r_ct:.3f}, p={p_ct:.2e}")
    print()

    results.append({"gene": gene,
                    "r_az": r_az, "p_az": p_az,
                    "r_ctrl": r_ct, "p_ctrl": p_ct})

res_df = pd.DataFrame(results)
res_df.to_csv("app_aplp_bulk_results.csv", index=False)

# ── 6. PLOTS ───────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.suptitle("APP/APLP1/APLP2 — Bulk DLPFC Expression & Presynaptic Module", fontsize=13)

genes = [g for g in TARGET_GENES if g in target_df.columns]

# Row 1: correlation with AZ module
for i, gene in enumerate(genes):
    ax = axes[0, i]
    y = target_df[gene]
    common = y.index.intersection(az_score.index)
    x_v = az_score.loc[common].values
    y_v = y.loc[common].values
    r, p = stats.pearsonr(x_v, y_v)

    ax.scatter(x_v, y_v, s=6, alpha=0.4, color="steelblue")
    m, b = np.polyfit(x_v, y_v, 1)
    xl = np.linspace(x_v.min(), x_v.max(), 100)
    ax.plot(xl, m*xl+b, color="firebrick", linewidth=1.5)
    ax.set_xlabel("Presynaptic AZ module score", fontsize=9)
    ax.set_ylabel(f"{gene} (logCPM)", fontsize=9)
    ax.set_title(f"{gene} vs AZ module\nr={r:.3f}, p={p:.1e}", fontsize=10)

# Row 2: AZ vs Ctrl correlation comparison (bar plot per gene)
ax = axes[1, 0]
x = np.arange(len(genes))
w = 0.35
r_az_vals   = [res_df[res_df.gene==g]["r_az"].values[0]   for g in genes]
r_ctrl_vals = [res_df[res_df.gene==g]["r_ctrl"].values[0] for g in genes]
ax.bar(x - w/2, r_az_vals,   w, label="AZ module",   color="steelblue", alpha=0.8)
ax.bar(x + w/2, r_ctrl_vals, w, label="Ctrl module", color="grey",      alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(genes)
ax.set_ylabel("Pearson r", fontsize=9)
ax.set_title("Module specificity", fontsize=10)
ax.axhline(0, color="black", linewidth=0.5)
ax.legend(fontsize=8)

# Row 2 middle: AZ module gene expression heatmap (mean)
ax = axes[1, 1]
az_means = az_df.mean().sort_values(ascending=False)
ax.barh(az_means.index, az_means.values, color="steelblue", alpha=0.8)
ax.set_xlabel("Mean logCPM", fontsize=9)
ax.set_title("AZ module gene expression", fontsize=10)

# Row 2 right: target gene expression bar
ax = axes[1, 2]
tgt_means = target_df.mean().sort_values(ascending=False)
ax.bar(tgt_means.index, tgt_means.values, color="firebrick", alpha=0.8)
ax.set_ylabel("Mean logCPM", fontsize=9)
ax.set_title("APP/APLP1/APLP2 expression", fontsize=10)

plt.tight_layout()
plt.savefig("app_aplp_bulk_results.png", dpi=200, bbox_inches="tight")
print("Saved: app_aplp_bulk_results.png")
print("\nDone.")
