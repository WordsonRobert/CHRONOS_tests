"""
CHRONOS - Hypothesis 4: BACE1 and PERK Pathway
================================================
Tests A-H:
  A - BACE1 ~ Braak (unadjusted + adjusted + cell composition)
  B - EIF2AK3 (PERK) ~ Braak
  C - ATF4/CHOP individual target genes ~ Braak
  D - ATF4_score and CHOP_score ~ Braak
  E - Pattern classification (BACE1 vs ATF4/CHOP)
  F - GSK3B ~ Braak
  G - AD vs NCI secondary (cogdx)
  H - BACE1 RNA vs protein (if available in proteomics)

Input : ROSMAP_DLPFC_logCPM.tsv, ROSMAP_clinical.csv,
        ROSMAP_biospecimen_metadata.csv, ROSMAP_assay_rnaSeq_metadata.csv
        proteomics file (for Test H)
Output: h4_results.csv, h4_plots.png, h4_pattern.txt
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

# ── GENE SETS ──────────────────────────────────────────────────────────────────

GENES = {
    # Primary
    "BACE1":    "ENSG00000186318",
    "EIF2AK3":  "ENSG00000172071",  # PERK
    "GSK3B":    "ENSG00000082701",
    # ATF4 branch targets
    "ATF4":     "ENSG00000128272",
    "ASNS":     "ENSG00000070669",
    "TRIB3":    "ENSG00000101255",
    "SLC7A5":   "ENSG00000103257",
    # CHOP branch targets
    "DDIT3":    "ENSG00000175197",   # CHOP
    "PPP1R15A": "ENSG00000087074",   # GADD34
}

# Cell type markers for composition estimate
NEURON_MARKERS = ["ENSG00000102003","ENSG00000067715",
                  "ENSG00000132639","ENSG00000008056","ENSG00000157542"]
ASTRO_MARKERS  = ["ENSG00000131095","ENSG00000171885"]
MICRO_MARKERS  = ["ENSG00000197249","ENSG00000101439"]
OLIGO_MARKERS  = ["ENSG00000197971","ENSG00000123560"]

COLORS = {
    "BACE1":"firebrick","EIF2AK3":"steelblue","GSK3B":"purple",
    "ATF4":"orange","ASNS":"seagreen","TRIB3":"darkcyan",
    "SLC7A5":"olive","DDIT3":"brown","PPP1R15A":"grey"
}

# ── 1. LOAD ────────────────────────────────────────────────────────────────────

print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = clinical["age_at_visit_max"].replace("90+","90")
clinical["age_at_visit_max"] = pd.to_numeric(clinical["age_at_visit_max"], errors="coerce")
for col in ["braaksc","msex","pmi","cogdx"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

sample_ids = expr.columns.tolist()
rna_bio = biospec[biospec["specimenID"].isin(sample_ids)][["specimenID","individualID"]].drop_duplicates()
meta = rna_bio.merge(
    clinical[["individualID","braaksc","msex","age_at_visit_max","pmi","cogdx"]],
    on="individualID", how="left"
)
rna_meta_sub = rna_meta[["specimenID","rnaBatch"]].drop_duplicates()
meta = meta.merge(rna_meta_sub, on="specimenID", how="left").set_index("specimenID")
print(f"  Samples: {len(meta)}, Braak available: {meta['braaksc'].notna().sum()}")

# ── 2. EXTRACT GENES ───────────────────────────────────────────────────────────

print("\nChecking gene availability...")
gene_df = {}
for name, eid in GENES.items():
    if eid in expr.index:
        gene_df[name] = expr.loc[eid]
        print(f"  {name}: found")
    else:
        print(f"  {name}: NOT FOUND")
gene_df = pd.DataFrame(gene_df)

# ── 3. CELL TYPE SCORES ────────────────────────────────────────────────────────

def ct_score(markers):
    present = [m for m in markers if m in expr.index]
    if not present: return None
    vals = expr.loc[present]
    z = (vals - vals.mean(axis=1).values[:,None]) / (vals.std(axis=1).values[:,None] + 1e-10)
    return z.mean(axis=0)

neuron_score = ct_score(NEURON_MARKERS); neuron_score.name = "neuron_score"
astro_score  = ct_score(ASTRO_MARKERS);  astro_score.name  = "astro_score"
micro_score  = ct_score(MICRO_MARKERS);  micro_score.name  = "micro_score"
oligo_score  = ct_score(OLIGO_MARKERS);  oligo_score.name  = "oligo_score"

df = gene_df.join(meta).join(neuron_score).join(astro_score).join(micro_score).join(oligo_score)
df_braak = df.dropna(subset=["braaksc","age_at_visit_max","msex","pmi"]).copy()
df_braak["braaksc"] = df_braak["braaksc"].astype(float)

batch_dummies = pd.get_dummies(df_braak["rnaBatch"], prefix="batch", drop_first=True)
df_braak = pd.concat([df_braak, batch_dummies], axis=1)
batch_cols = list(batch_dummies.columns)
base_covs  = ["age_at_visit_max","msex","pmi"] + batch_cols
cell_covs  = ["neuron_score","astro_score","micro_score","oligo_score"]
print(f"\n  Analysis samples: {len(df_braak)}")

# ── OLS HELPER ─────────────────────────────────────────────────────────────────

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
    idx = cols.index(term)
    beta,se_b = coeffs[idx],se[idx]
    t = beta/se_b
    pv = 2*stats.t.sf(abs(t),df=n-p)
    return dict(beta=beta,se=se_b,ci_lo=beta-1.96*se_b,ci_hi=beta+1.96*se_b,p=pv,n=int(mask.sum()))

# ── BRAAK TRAJECTORY STATS ────────────────────────────────────────────────────

def braak_stats(df, gene):
    rows = []
    for b in sorted(df["braaksc"].unique()):
        v = df[df["braaksc"]==b][gene].dropna()
        rows.append({"braak":b,"n":len(v),"mean":v.mean(),"se":v.sem(),"median":v.median()})
    return pd.DataFrame(rows)

# ── TESTS A, B, F: PRIMARY GENES ~ BRAAK ──────────────────────────────────────

print("\n=== TESTS A/B/F: Primary genes ~ Braak ===")
primary_genes = ["BACE1","EIF2AK3","GSK3B"]
primary_results = []

for gene in primary_genes:
    if gene not in df_braak.columns: continue

    r_raw  = ols_term(df_braak[gene], df_braak[["braaksc"]], "braaksc")
    r_adj  = ols_term(df_braak[gene], df_braak[["braaksc"]+base_covs], "braaksc")
    r_cell = ols_term(df_braak[gene], df_braak[["braaksc"]+base_covs+cell_covs], "braaksc")

    print(f"  {gene}:")
    print(f"    Raw      : β={r_raw['beta']:.4f}, p={r_raw['p']:.3e}")
    print(f"    Adjusted : β={r_adj['beta']:.4f}, p={r_adj['p']:.3e}")
    print(f"    +Cell    : β={r_cell['beta']:.4f}, p={r_cell['p']:.3e}")

    primary_results.append({"gene":gene,
        "raw_beta":r_raw["beta"],"raw_p":r_raw["p"],
        "adj_beta":r_adj["beta"],"adj_p":r_adj["p"],
        "cell_beta":r_cell["beta"],"cell_p":r_cell["p"]})

# ── TEST C: ATF4/CHOP TARGETS ~ BRAAK ─────────────────────────────────────────

print("\n=== TEST C: ATF4/CHOP individual targets ~ Braak ===")
pathway_genes = ["ATF4","ASNS","TRIB3","SLC7A5","DDIT3","PPP1R15A"]
pathway_results = []

for gene in pathway_genes:
    if gene not in df_braak.columns: continue
    r_adj  = ols_term(df_braak[gene], df_braak[["braaksc"]+base_covs], "braaksc")
    r_cell = ols_term(df_braak[gene], df_braak[["braaksc"]+base_covs+cell_covs], "braaksc")
    print(f"  {gene}: β_adj={r_adj['beta']:.4f} p={r_adj['p']:.3e} | β_cell={r_cell['beta']:.4f} p={r_cell['p']:.3e}")
    pathway_results.append({"gene":gene,
        "adj_beta":r_adj["beta"],"adj_p":r_adj["p"],
        "cell_beta":r_cell["beta"],"cell_p":r_cell["p"]})

# FDR across all genes
all_res = primary_results + pathway_results
all_ps  = [r["adj_p"] for r in all_res]
_, fdrs, _, _ = multipletests(all_ps, method="fdr_bh")
for r, fdr in zip(all_res, fdrs):
    r["fdr"] = fdr

print("\n  FDR-corrected summary:")
for r in all_res:
    print(f"    {r['gene']}: β={r['adj_beta']:.4f}, p={r['adj_p']:.3e}, FDR={r['fdr']:.3f}")

# ── TEST D: ATF4 AND CHOP MODULE SCORES ───────────────────────────────────────

print("\n=== TEST D: ATF4 and CHOP module scores ~ Braak ===")
atf4_genes = [g for g in ["ASNS","TRIB3","SLC7A5"] if g in df_braak.columns]
chop_genes = [g for g in ["DDIT3","PPP1R15A"] if g in df_braak.columns]

if atf4_genes:
    atf4_z = df_braak[atf4_genes].apply(lambda x: (x-x.mean())/x.std())
    df_braak["ATF4_score"] = atf4_z.mean(axis=1)
    r_atf4 = ols_term(df_braak["ATF4_score"], df_braak[["braaksc"]+base_covs+cell_covs], "braaksc")
    print(f"  ATF4_score ~ Braak: β={r_atf4['beta']:.4f}, p={r_atf4['p']:.3e}")

if chop_genes:
    chop_z = df_braak[chop_genes].apply(lambda x: (x-x.mean())/x.std())
    df_braak["CHOP_score"] = chop_z.mean(axis=1)
    r_chop = ols_term(df_braak["CHOP_score"], df_braak[["braaksc"]+base_covs+cell_covs], "braaksc")
    print(f"  CHOP_score ~ Braak: β={r_chop['beta']:.4f}, p={r_chop['p']:.3e}")

# ── TEST E: PATTERN CLASSIFICATION ────────────────────────────────────────────

print("\n=== TEST E: Pattern classification ===")
bace1_up  = next((r for r in all_res if r["gene"]=="BACE1"), None)
atf4_sig  = any(r["adj_p"]<0.05 for r in all_res if r["gene"] in ["ASNS","TRIB3","SLC7A5"])
chop_sig  = any(r["adj_p"]<0.05 for r in all_res if r["gene"] in ["DDIT3","PPP1R15A"])

if bace1_up:
    bace1_direction = "UP" if bace1_up["adj_beta"]>0 and bace1_up["adj_p"]<0.05 else \
                      "DOWN" if bace1_up["adj_beta"]<0 and bace1_up["adj_p"]<0.05 else "FLAT"
else:
    bace1_direction = "NOT FOUND"

pattern_str = f"""
PATTERN CLASSIFICATION (Test E)
================================
BACE1 direction  : {bace1_direction}
ATF4 targets sig : {'YES' if atf4_sig else 'NO'}
CHOP targets sig : {'YES' if chop_sig else 'NO'}

Interpretation:
"""
if bace1_direction=="UP" and not atf4_sig and not chop_sig:
    pattern_str += "BACE1 UP + ATF4/CHOP FLAT → consistent with attenuated PERK-like pattern (Hou 2017)"
elif bace1_direction=="UP" and (atf4_sig or chop_sig):
    pattern_str += "BACE1 UP + ATF4/CHOP UP → broad ISR activation (less specific)"
elif bace1_direction=="FLAT" and (atf4_sig or chop_sig):
    pattern_str += "BACE1 FLAT + ATF4/CHOP UP → stress response without BACE1 transcript change"
else:
    pattern_str += f"BACE1 {bace1_direction} + ATF4/CHOP mixed → see individual results"

print(pattern_str)

# ── TEST G: AD VS NCI ─────────────────────────────────────────────────────────

print("\n=== TEST G: AD vs NCI (cogdx) ===")
df_braak["AD"] = (df_braak["cogdx"] >= 4).astype(float)
for gene in ["BACE1","EIF2AK3","GSK3B"]:
    if gene not in df_braak.columns: continue
    r = ols_term(df_braak[gene], df_braak[["AD","braaksc"]+base_covs+cell_covs], "AD")
    print(f"  {gene} ~ AD (adj for Braak+cells): β={r['beta']:.4f}, p={r['p']:.3e}")

# ── TEST H: BACE1 RNA VS PROTEIN ─────────────────────────────────────────────

print("\n=== TEST H: BACE1 RNA vs Protein ===")
prot_file = DATA + "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"
try:
    prot = pd.read_csv(prot_file, index_col=0)
    bace1_prot = None
    for idx in prot.index:
        if "BACE1" in str(idx) or "BACE" in str(idx):
            bace1_prot = prot.loc[idx]
            print(f"  Found BACE1 in proteomics: {idx}")
            break
    if bace1_prot is None:
        # try searching columns
        bace1_cols = [c for c in prot.index if "BACE" in str(c).upper()]
        if bace1_cols:
            bace1_prot = prot.loc[bace1_cols[0]]
            print(f"  Found: {bace1_cols[0]}")
        else:
            print("  BACE1 not found in proteomics index — skipping Test H")
except Exception as e:
    print(f"  Could not load proteomics: {e}")
    bace1_prot = None

# ── SAVE RESULTS ──────────────────────────────────────────────────────────────

pd.DataFrame(all_res).to_csv("h4_gene_results.csv", index=False)
with open("h4_pattern.txt","w") as f:
    f.write(pattern_str)
print("\nSaved CSVs.")

# ── PLOTS ──────────────────────────────────────────────────────────────────────

print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("CHRONOS — BACE1/PERK Pathway ~ Braak (Hypothesis 4)", fontsize=13)

# Plot 1: BACE1 trajectory
ax = axes[0,0]
if "BACE1" in df_braak.columns:
    bs = braak_stats(df_braak, "BACE1")
    ax.errorbar(bs["braak"], bs["mean"], yerr=bs["se"], fmt="o-",
                color="firebrick", linewidth=2, capsize=4)
    ax.set_xlabel("Braak stage"); ax.set_ylabel("BACE1 (logCPM)")
    b1 = next((r for r in all_res if r["gene"]=="BACE1"),{})
    ax.set_title(f"BACE1 ~ Braak\nβ_adj={b1.get('adj_beta',np.nan):.4f}, p={b1.get('adj_p',np.nan):.2e}", fontsize=10)

# Plot 2: EIF2AK3 trajectory
ax = axes[0,1]
if "EIF2AK3" in df_braak.columns:
    bs = braak_stats(df_braak, "EIF2AK3")
    ax.errorbar(bs["braak"], bs["mean"], yerr=bs["se"], fmt="o-",
                color="steelblue", linewidth=2, capsize=4)
    ax.set_xlabel("Braak stage"); ax.set_ylabel("EIF2AK3/PERK (logCPM)")
    b1 = next((r for r in all_res if r["gene"]=="EIF2AK3"),{})
    ax.set_title(f"EIF2AK3 (PERK) ~ Braak\nβ_adj={b1.get('adj_beta',np.nan):.4f}, p={b1.get('adj_p',np.nan):.2e}", fontsize=10)

# Plot 3: GSK3B trajectory
ax = axes[0,2]
if "GSK3B" in df_braak.columns:
    bs = braak_stats(df_braak, "GSK3B")
    ax.errorbar(bs["braak"], bs["mean"], yerr=bs["se"], fmt="o-",
                color="purple", linewidth=2, capsize=4)
    ax.set_xlabel("Braak stage"); ax.set_ylabel("GSK3B (logCPM)")
    b1 = next((r for r in all_res if r["gene"]=="GSK3B"),{})
    ax.set_title(f"GSK3B ~ Braak\nβ_adj={b1.get('adj_beta',np.nan):.4f}, p={b1.get('adj_p',np.nan):.2e}", fontsize=10)

# Plot 4: ATF4/CHOP pathway betas
ax = axes[1,0]
pw = pd.DataFrame(pathway_results)
if not pw.empty:
    colors_bar = [COLORS.get(g,"grey") for g in pw["gene"]]
    bars = ax.bar(pw["gene"], pw["adj_beta"], color=colors_bar, alpha=0.8)
    ax.axhline(0, color="black", linewidth=0.8)
    for bar, p in zip(bars, pw["adj_p"]):
        sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        if sig:
            ax.text(bar.get_x()+bar.get_width()/2,
                    bar.get_height()+0.001, sig, ha="center", fontsize=10)
    ax.set_ylabel("β (Braak, adjusted)"); ax.set_title("ATF4/CHOP targets ~ Braak", fontsize=10)
    ax.tick_params(axis="x", rotation=30)

# Plot 5: ATF4 and CHOP scores
ax = axes[1,1]
for score, color, label in [("ATF4_score","orange","ATF4 targets"),
                              ("CHOP_score","brown","CHOP targets")]:
    if score in df_braak.columns:
        bs = df_braak.groupby("braaksc")[score].agg(["mean","sem"]).reset_index()
        ax.errorbar(bs["braaksc"], bs["mean"], yerr=bs["sem"], fmt="o-",
                    color=color, label=label, linewidth=2, capsize=3)
ax.axhline(0, color="black", linewidth=0.5)
ax.set_xlabel("Braak stage"); ax.set_ylabel("Module score (z-scored)")
ax.set_title("ATF4 and CHOP module scores ~ Braak", fontsize=10)
ax.legend(fontsize=8)

# Plot 6: summary beta comparison (adjusted vs cell-adjusted)
ax = axes[1,2]
summary_df = pd.DataFrame(all_res)
if not summary_df.empty:
    x = np.arange(len(summary_df)); w=0.35
    ax.bar(x-w/2, summary_df["adj_beta"],  w, label="Adjusted",      color="steelblue", alpha=0.8)
    ax.bar(x+w/2, summary_df["cell_beta"], w, label="+Cell comp",     color="firebrick", alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(summary_df["gene"], rotation=35, fontsize=8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("β (Braak term)")
    ax.set_title("All genes: adj vs cell-composition adj", fontsize=10)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("h4_plots.png", dpi=200, bbox_inches="tight")
print("Saved: h4_plots.png")
print("\nDone.")
