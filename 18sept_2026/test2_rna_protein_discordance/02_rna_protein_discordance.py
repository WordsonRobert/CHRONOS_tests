"""
02_rna_protein_discordance.py
Test 2: RNA vs Protein Discordance for TMED2, TMED10, TMED9
CHRONOS / ROSMAP
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, warnings
warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"
OUT  = "./"

GENES = {
    "TMED2":  "ENSG00000086598",
    "TMED10": "ENSG00000170348",
    "TMED9":  "ENSG00000184840",
}

CELL_MARKERS = {
    "neuron_score": ["ENSG00000102003","ENSG00000067715","ENSG00000132639",
                     "ENSG00000008056","ENSG00000157542"],
    "astro_score":  ["ENSG00000131095","ENSG00000171885"],
    "micro_score":  ["ENSG00000204472","ENSG00000138185"],
    "oligo_score":  ["ENSG00000197971","ENSG00000123560"],
}
CELL_COLS = list(CELL_MARKERS.keys())

# ── OLS helpers ─────────────────────────────────────────────────────────
def ols_term(y_s, X_df, term):
    cols = list(X_df.columns)
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    y_v = y_s.values if hasattr(y_s,"values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(y_v))
    if mask.sum() < len(cols)+3:
        return dict(beta=np.nan,se=np.nan,p=np.nan,n=int(mask.sum()))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm@coeffs
    n, p = len(y_v), Xm.shape[1]
    mse = np.sum(resid**2)/(n-p)
    try:
        cov = mse*np.linalg.inv(Xm.T@Xm)
    except np.linalg.LinAlgError:
        return dict(beta=np.nan,se=np.nan,p=np.nan,n=int(mask.sum()))
    se = np.sqrt(np.diag(cov))
    idx = (["intercept"]+cols).index(term)
    beta, se_b = coeffs[idx], se[idx]
    t = beta/se_b
    pv = 2*stats.t.sf(abs(t), df=n-p)
    return dict(beta=beta,se=se_b,p=pv,n=int(mask.sum()))

def ols_r2(y_s, X_df):
    cols = list(X_df.columns)
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    y_v = y_s.values if hasattr(y_s,"values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    ss_res = np.sum((y_v - Xm@coeffs)**2)
    ss_tot = np.sum((y_v - y_v.mean())**2)
    return (1-ss_res/ss_tot if ss_tot>0 else np.nan), int(mask.sum())

def ols_resid(y_s, X_df):
    cols = list(X_df.columns)
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    y_v = y_s.values if hasattr(y_s,"values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(y_v))
    out = np.full(len(y_s), np.nan)
    Xm_m, y_m = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm_m, y_m, rcond=None)
    out[mask] = y_m - Xm_m@coeffs
    return pd.Series(out, index=y_s.index if hasattr(y_s,"index") else None)

# ─────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────────────────
print("Loading data...")
expr     = pd.read_csv(DATA+"ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA+"ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA+"ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA+"ROSMAP_assay_rnaSeq_metadata.csv")
prot_raw = pd.read_csv(DATA+"C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv",
                       index_col=0)
prot_meta = pd.read_csv(DATA+"ROSMAP_assay_proteomics_TMTquantitation_metadata.csv")

# clean clinical
clinical["age_at_visit_max"] = pd.to_numeric(
    clinical["age_at_visit_max"].replace("90+","90"), errors="coerce")
for col in ["braaksc","msex","pmi","cogdx"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")
clinical["AD"] = np.nan
clinical.loc[clinical["cogdx"]==1, "AD"] = 0
clinical.loc[clinical["cogdx"]>=4, "AD"] = 1

# ─────────────────────────────────────────────────────────────────────────
# BUILD RNA COHORT
# ─────────────────────────────────────────────────────────────────────────
print("Building RNA cohort...")

# rna_meta specimenID → biospec → individualID
rna_bio = biospec[["individualID","specimenID"]].drop_duplicates()
rna_link = rna_meta[["specimenID","rnaBatch"]].merge(rna_bio, on="specimenID", how="left")
rna_link = rna_link[rna_link["specimenID"].isin(expr.columns)]
rna_link = rna_link.merge(
    clinical[["individualID","braaksc","msex","pmi","age_at_visit_max","cogdx","AD"]],
    on="individualID", how="left")
rna_link = rna_link.dropna(subset=["individualID"])

# rnaBatch: encode only if it has >1 unique non-null value; otherwise drop it
_rna_batch_vals = rna_link["rnaBatch"].dropna().unique()
if len(_rna_batch_vals) > 1:
    rna_link["rnaBatch_cat"] = pd.Categorical(rna_link["rnaBatch"]).codes.astype(float)
    rna_link.loc[rna_link["rnaBatch"].isna(), "rnaBatch_cat"] = np.nan
    _USE_RNA_BATCH = True
else:
    rna_link["rnaBatch_cat"] = 0.0   # constant — will be excluded from models
    _USE_RNA_BATCH = False
print(f"  rnaBatch unique values: {len(_rna_batch_vals)}  — using in model: {_USE_RNA_BATCH}")

print(f"  RNA samples linked: {len(rna_link)}")

# cell scores
def cell_score(markers):
    avail = [m for m in markers if m in expr.index]
    if not avail: return pd.Series(np.nan, index=expr.columns)
    z = expr.loc[avail].apply(lambda r: (r-r.mean())/r.std(), axis=1)
    return z.mean(axis=0)

for cname, markers in CELL_MARKERS.items():
    scores = cell_score(markers)
    rna_link[cname] = rna_link["specimenID"].map(scores)

# pull TMED RNA values
for name, eid in GENES.items():
    if eid in expr.index:
        rna_link[f"rna_{name}"] = rna_link["specimenID"].map(expr.loc[eid])

# ─────────────────────────────────────────────────────────────────────────
# BUILD PROTEIN COHORT
# ─────────────────────────────────────────────────────────────────────────
print("Building protein cohort...")

# prot_meta: specimenID format is like "ROSMAP.DLPFC.b01.127C.R8316516"
# last part after final "." is the individualID (R-number)
# batchChannel is like "b01.127C"
pm = prot_meta[prot_meta["isAssayControl"]==False].copy()
pm = pm[["specimenID","batchChannel","batch"]].drop_duplicates()

# extract individualID from specimenID (last segment after last ".")
pm["individualID"] = pm["specimenID"].str.split(".").str[-1]

# keep only batchChannels that exist in prot_raw columns
pm = pm[pm["batchChannel"].isin(prot_raw.columns)]
print(f"  Protein samples (non-control, in data): {len(pm)}")

# map each batchChannel → its protein values, then group by individualID
# prot_raw rows = proteins, cols = batchChannels
# we need: for each protein gene, batchChannel → value → individualID
# Build a per-protein series: individualID → value (mean if duplicated)

# find TMED proteins
tmed_prot_map = {}
for name in GENES:
    matches = [g for g in prot_raw.index if g.startswith(name+"|")]
    if matches:
        tmed_prot_map[name] = matches[0]
print(f"  TMED proteins found: {tmed_prot_map}")

# build a small df: individualID × [TMED proteins + batch]
bc_to_id   = pm.set_index("batchChannel")["individualID"].to_dict()
bc_to_batch = pm.set_index("batchChannel")["batch"].to_dict()

prot_cols = list(tmed_prot_map.values())
sub_prot  = prot_raw.loc[prot_cols, pm["batchChannel"].values]  # proteins × batchChannels

# transpose → batchChannels × proteins, add individualID
sub_t = sub_prot.T.copy()
sub_t.index.name = "batchChannel"
sub_t = sub_t.reset_index()
sub_t["individualID"] = sub_t["batchChannel"].map(bc_to_id)
sub_t["prot_batch"]   = sub_t["batchChannel"].map(bc_to_batch).astype(str)

# average per individual (rare duplicates)
sub_t = sub_t.drop(columns=["batchChannel"])
sub_t[prot_cols] = sub_t[prot_cols].apply(pd.to_numeric, errors="coerce")
# batch is a string — extract before numeric mean
prot_batch_mode = sub_t.groupby("individualID")["prot_batch"].first().reset_index()
prot_df = sub_t.drop(columns=["prot_batch"]).groupby("individualID").mean().reset_index()
prot_df = prot_df.merge(prot_batch_mode, on="individualID")

# rename columns to short names
for name, pname in tmed_prot_map.items():
    prot_df.rename(columns={pname: f"prot_{name}"}, inplace=True)

# merge with clinical
prot_clin = prot_df.merge(
    clinical[["individualID","braaksc","msex","pmi","age_at_visit_max","cogdx","AD"]],
    on="individualID", how="left")
prot_clin["prot_batch_cat"] = pd.Categorical(prot_clin["prot_batch"]).codes

# add cell scores to prot_clin via rna_link
cell_from_rna = rna_link[["individualID"] + CELL_COLS].drop_duplicates("individualID")
prot_clin = prot_clin.merge(cell_from_rna, on="individualID", how="left")

print(f"  Protein individuals (after clinical merge): {len(prot_clin)}")

# ─────────────────────────────────────────────────────────────────────────
# TEST A: AVAILABILITY
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST A: AVAILABILITY ===")
print("="*60)

rna_ids  = set(rna_link["individualID"].dropna())
prot_ids = set(prot_clin["individualID"].dropna())
meta_ids = set(clinical["individualID"].dropna())
paired   = rna_ids & prot_ids & meta_ids

print(f"\nRNA samples:              {len(rna_ids)}")
print(f"Protein individuals:      {len(prot_ids)}")
print(f"RNA ∩ protein ∩ metadata: {len(paired)}")

print("\nGene availability in paired:")
for name, eid in GENES.items():
    rna_ok  = eid in expr.index
    prot_ok = name in tmed_prot_map
    n_paired = 0
    if rna_ok and prot_ok:
        pr = rna_link[rna_link["individualID"].isin(paired) & rna_link[f"rna_{name}"].notna()]
        pp = prot_clin[prot_clin["individualID"].isin(paired) & prot_clin[f"prot_{name}"].notna()]
        n_paired = len(set(pr["individualID"]) & set(pp["individualID"]))
    print(f"  {name}: RNA={'YES' if rna_ok else 'NO'}, "
          f"Protein={'YES ('+tmed_prot_map.get(name,'')+')' if prot_ok else 'NOT FOUND'}, "
          f"Paired n={n_paired}")

print(f"\nProtein batch col: prot_batch_cat (from 'batch' in metadata)")

if len(paired)==0:
    print("ERROR: 0 paired subjects — check individualID format mismatch")
    print("  rna_link individualID sample:", list(rna_link["individualID"].dropna()[:3]))
    print("  prot_clin individualID sample:", list(prot_clin["individualID"].dropna()[:3]))
    import sys; sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────
# COVARIATE COLUMN SETS
# ─────────────────────────────────────────────────────────────────────────
_RNA_BATCH_COLS = ["rnaBatch_cat"] if _USE_RNA_BATCH else []

RNA_COVS_SIMPLE = ["AD","age_at_visit_max","msex","pmi"] + _RNA_BATCH_COLS
RNA_COVS_CELLS  = RNA_COVS_SIMPLE + CELL_COLS

PROT_COVS_SIMPLE = ["AD","age_at_visit_max","msex","pmi","prot_batch_cat"]
PROT_COVS_CELLS  = PROT_COVS_SIMPLE + CELL_COLS

RNA_BRAAK_COVS   = ["braaksc","age_at_visit_max","msex","pmi"] + _RNA_BATCH_COLS + CELL_COLS
PROT_BRAAK_COVS  = ["braaksc","age_at_visit_max","msex","pmi","prot_batch_cat"] + CELL_COLS

# ─────────────────────────────────────────────────────────────────────────
# TEST B: RNA ~ AD
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST B: RNA ~ AD ===")
print("="*60)

rna_ad = rna_link.dropna(subset=RNA_COVS_CELLS)
rna_ad = rna_ad[rna_ad["AD"].isin([0,1])]
print(f"\nn = {len(rna_ad)}")
print(f"{'Gene':<10} {'β_AD(simple)':>14} {'p(simple)':>10} {'β_AD(+cells)':>14} {'p(+cells)':>10}")
print("-"*55)

testB = {}
for name in GENES:
    col = f"rna_{name}"
    if col not in rna_ad.columns: print(f"  {name}: missing"); continue
    r_s = ols_term(rna_ad[col], rna_ad[RNA_COVS_SIMPLE], "AD")
    r_c = ols_term(rna_ad[col], rna_ad[RNA_COVS_CELLS],  "AD")
    testB[name] = {"simple":r_s,"cells":r_c}
    print(f"  {name:<8} {r_s['beta']:>14.4f} {r_s['p']:>10.4f} {r_c['beta']:>14.4f} {r_c['p']:>10.4f}  n={r_c['n']}")

# ─────────────────────────────────────────────────────────────────────────
# TEST C: Protein ~ AD
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST C: PROTEIN ~ AD ===")
print("="*60)

prot_ad = prot_clin.dropna(subset=PROT_COVS_CELLS)
prot_ad = prot_ad[prot_ad["AD"].isin([0,1])]
print(f"\nn = {len(prot_ad)}")
print(f"{'Gene':<10} {'β_AD(simple)':>14} {'p(simple)':>10} {'β_AD(+cells)':>14} {'p(+cells)':>10}")
print("-"*55)

testC = {}
for name in GENES:
    col = f"prot_{name}"
    if col not in prot_ad.columns: print(f"  {name}: missing"); continue
    r_s = ols_term(prot_ad[col], prot_ad[PROT_COVS_SIMPLE], "AD")
    r_c = ols_term(prot_ad[col], prot_ad[PROT_COVS_CELLS],  "AD")
    testC[name] = {"simple":r_s,"cells":r_c}
    print(f"  {name:<8} {r_s['beta']:>14.4f} {r_s['p']:>10.4f} {r_c['beta']:>14.4f} {r_c['p']:>10.4f}  n={r_c['n']}")

# ─────────────────────────────────────────────────────────────────────────
# TEST D: Braak — RNA vs Protein
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST D: BRAAK — RNA vs PROTEIN ===")
print("="*60)

rna_braak  = rna_link.dropna(subset=RNA_BRAAK_COVS)
prot_braak = prot_clin.dropna(subset=PROT_BRAAK_COVS)
print(f"\nRNA n={len(rna_braak)}  Protein n={len(prot_braak)}")
print(f"\n{'Gene':<10} {'RNA β_Braak':>12} {'RNA p':>8} {'Prot β_Braak':>13} {'Prot p':>8}")
print("-"*55)

testD = {}
for name in GENES:
    rc = f"rna_{name}";  pc = f"prot_{name}"
    r = ols_term(rna_braak[rc],  rna_braak[RNA_BRAAK_COVS],  "braaksc") if rc in rna_braak.columns  else None
    p = ols_term(prot_braak[pc], prot_braak[PROT_BRAAK_COVS],"braaksc") if pc in prot_braak.columns else None
    testD[name] = {"rna":r,"prot":p}
    rb = f"{r['beta']:12.4f}" if r else "         N/A"
    rp = f"{r['p']:8.4f}"     if r else "     N/A"
    pb = f"{p['beta']:13.4f}" if p else "          N/A"
    pp = f"{p['p']:8.4f}"     if p else "     N/A"
    print(f"  {name:<8} {rb} {rp} {pb} {pp}")

# ─────────────────────────────────────────────────────────────────────────
# PAIRED COHORT (for E/F/G/H/I)
# ─────────────────────────────────────────────────────────────────────────
# merge RNA and protein on individualID (paired subjects only)
rna_p = rna_link[rna_link["individualID"].isin(paired)].copy()
prot_p = prot_clin[prot_clin["individualID"].isin(paired)].copy()

# one RNA row per individual (take first if duplicates)
rna_p = rna_p.drop_duplicates("individualID").set_index("individualID")
prot_p = prot_p.drop_duplicates("individualID").set_index("individualID")

paired_idx = rna_p.index.intersection(prot_p.index)

rna_cols  = ["braaksc","msex","pmi","age_at_visit_max","AD","rnaBatch_cat"] + CELL_COLS + \
            [f"rna_{n}" for n in GENES if f"rna_{n}" in rna_p.columns]
prot_cols2 = [f"prot_{n}" for n in GENES if f"prot_{n}" in prot_p.columns] + ["prot_batch_cat"]

paired_df = rna_p.loc[paired_idx, [c for c in rna_cols if c in rna_p.columns]].copy()
for c in [c for c in prot_cols2 if c in prot_p.columns]:
    paired_df[c] = prot_p.loc[paired_idx, c]

print(f"\nPaired subjects (merged df): {len(paired_df)}")

PAIRED_R_COVS = ["age_at_visit_max","msex","pmi"] + _RNA_BATCH_COLS + CELL_COLS
PAIRED_P_COVS = ["age_at_visit_max","msex","pmi","prot_batch_cat"] + CELL_COLS

# ─────────────────────────────────────────────────────────────────────────
# TEST E: STANDARDIZED DISCORDANCE
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST E: STANDARDIZED DISCORDANCE (PAIRED) ===")
print("="*60)

print(f"\n{'Gene':<10} {'β_RNA_z':>10} {'p_RNA':>8} {'β_Prot_z':>10} {'p_Prot':>8} {'Δβ':>8} {'D':>8} {'n':>5}")
print("-"*65)

testE = {}
for name in GENES:
    rc = f"rna_{name}"; pc = f"prot_{name}"
    if rc not in paired_df.columns or pc not in paired_df.columns:
        print(f"  {name}: missing"); continue
    sub = paired_df[[rc,pc,"AD"] + PAIRED_R_COVS + ["prot_batch_cat"]].dropna()
    if len(sub) < 20: print(f"  {name}: n={len(sub)} too small"); continue
    rna_z  = pd.Series((sub[rc]-sub[rc].mean())/sub[rc].std(), index=sub.index)
    prot_z = pd.Series((sub[pc]-sub[pc].mean())/sub[pc].std(), index=sub.index)
    cov_r  = sub[["AD"] + PAIRED_R_COVS]
    cov_p  = sub[["AD"] + PAIRED_P_COVS]
    r_r = ols_term(rna_z,  cov_r, "AD")
    r_p = ols_term(prot_z, cov_p, "AD")
    db = r_p["beta"] - r_r["beta"]
    D  = abs(r_p["beta"])/abs(r_r["beta"]) if abs(r_r["beta"])>1e-10 else np.nan
    testE[name] = {"rna":r_r,"prot":r_p,"delta_beta":db,"D":D}
    print(f"  {name:<8} {r_r['beta']:>10.4f} {r_r['p']:>8.4f} {r_p['beta']:>10.4f} {r_p['p']:>8.4f} "
          f"{db:>8.4f} {D:>8.4f} {r_r['n']:>5}")

# ─────────────────────────────────────────────────────────────────────────
# TEST F: RNA-PROTEIN SPEARMAN
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST F: RNA-PROTEIN SPEARMAN (PAIRED) ===")
print("="*60)

for name in GENES:
    rc = f"rna_{name}"; pc = f"prot_{name}"
    if rc not in paired_df.columns or pc not in paired_df.columns: continue
    print(f"\n  {name}")
    sub_all = paired_df[[rc,pc,"AD"]].dropna()
    sub_nci = sub_all[sub_all["AD"]==0]
    sub_ad  = sub_all[sub_all["AD"]==1]
    for label, sub in [("all",sub_all),("NCI",sub_nci),("AD",sub_ad)]:
        if len(sub)<5: print(f"    {label}: n too small"); continue
        rho,pv = spearmanr(sub[rc], sub[pc])
        print(f"    {label:3}: rho={rho:7.4f}  p={pv:.4f}  n={len(sub)}")

# ─────────────────────────────────────────────────────────────────────────
# TEST G: PROTEIN ~ RNA + BRAAK + covariates (PAIRED)
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST G: PROTEIN ~ RNA + BRAAK (conditional, PAIRED) ===")
print("="*60)

for name in GENES:
    rc = f"rna_{name}"; pc = f"prot_{name}"
    if rc not in paired_df.columns or pc not in paired_df.columns: continue
    sub = paired_df[[rc,pc,"braaksc","AD"] + PAIRED_P_COVS].dropna()
    if len(sub)<20: print(f"\n  {name}: n={len(sub)} too small"); continue
    print(f"\n  {name}  (n={len(sub)})")

    # Protein ~ Braak + covariates
    r0 = ols_term(sub[pc], sub[["braaksc"]+PAIRED_P_COVS], "braaksc")
    print(f"    Protein ~ Braak + cov:         β_Braak={r0['beta']:7.4f}  p={r0['p']:.4f}")

    # Protein ~ Braak + RNA + covariates
    r1 = ols_term(sub[pc], sub[[rc,"braaksc"]+PAIRED_P_COVS], "braaksc")
    r1r = ols_term(sub[pc], sub[[rc,"braaksc"]+PAIRED_P_COVS], rc)
    print(f"    Protein ~ Braak+RNA+cov:       β_Braak={r1['beta']:7.4f}  p={r1['p']:.4f}  β_RNA={r1r['beta']:7.4f}  p_RNA={r1r['p']:.4f}")

    # with AD
    ra0 = ols_term(sub[pc], sub[["AD"]+PAIRED_P_COVS], "AD")
    ra1 = ols_term(sub[pc], sub[[rc,"AD"]+PAIRED_P_COVS], "AD")
    print(f"    Protein ~ AD + cov:            β_AD={ra0['beta']:7.4f}  p={ra0['p']:.4f}")
    print(f"    Protein ~ AD+RNA+cov:          β_AD={ra1['beta']:7.4f}  p={ra1['p']:.4f}")

# ─────────────────────────────────────────────────────────────────────────
# TEST H: ΔR²
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST H: INCREMENTAL R² (PAIRED) ===")
print("="*60)
print(f"\n{'Gene':<10} {'R²(M0)':>10} {'R²(M1+RNA)':>12} {'ΔR²':>8} {'β_Braak M0':>12} {'β_Braak M1':>12}")
print("-"*60)

for name in GENES:
    rc = f"rna_{name}"; pc = f"prot_{name}"
    if rc not in paired_df.columns or pc not in paired_df.columns: continue
    sub = paired_df[[rc,pc,"braaksc"]+PAIRED_P_COVS].dropna()
    if len(sub)<20: continue
    r2_m0,_ = ols_r2(sub[pc], sub[["braaksc"]+PAIRED_P_COVS])
    r2_m1,_ = ols_r2(sub[pc], sub[[rc,"braaksc"]+PAIRED_P_COVS])
    b0 = ols_term(sub[pc], sub[["braaksc"]+PAIRED_P_COVS], "braaksc")["beta"]
    b1 = ols_term(sub[pc], sub[[rc,"braaksc"]+PAIRED_P_COVS], "braaksc")["beta"]
    print(f"  {name:<8} {r2_m0:>10.4f} {r2_m1:>12.4f} {r2_m1-r2_m0:>8.4f} {b0:>12.4f} {b1:>12.4f}")

# ─────────────────────────────────────────────────────────────────────────
# TEST I: PROTEIN RESIDUAL ~ BRAAK/AD
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST I: PROTEIN RESIDUAL ~ BRAAK/AD (PAIRED) ===")
print("="*60)
print(f"\n{'Gene':<10} {'β_resid~Braak':>15} {'p_Braak':>9} {'β_resid~AD':>12} {'p_AD':>7}")
print("-"*55)

for name in GENES:
    rc = f"rna_{name}"; pc = f"prot_{name}"
    if rc not in paired_df.columns or pc not in paired_df.columns: continue
    resid_base = [rc,"age_at_visit_max","msex","pmi","prot_batch_cat"] + CELL_COLS
    sub = paired_df[[rc,pc,"braaksc","AD"]+PAIRED_P_COVS].dropna()
    if len(sub)<20: continue
    prot_resid = ols_resid(sub[pc], sub[[c for c in resid_base if c in sub.columns]])
    prot_resid = prot_resid.dropna()
    sub2 = sub.loc[prot_resid.index]
    rb = ols_term(prot_resid, sub2[["braaksc"]], "braaksc")
    ra = ols_term(prot_resid, sub2[["AD"]],      "AD")
    print(f"  {name:<8} {rb['beta']:>15.4f} {rb['p']:>9.4f} {ra['beta']:>12.4f} {ra['p']:>7.4f}")

# ─────────────────────────────────────────────────────────────────────────
# TEST J: PLOTS
# ─────────────────────────────────────────────────────────────────────────
print("\n=== TEST J: PLOTS ===")
gene_list = [n for n in GENES if f"rna_{n}" in paired_df.columns and f"prot_{n}" in paired_df.columns]

if gene_list:
    fig, axes = plt.subplots(len(gene_list), 3, figsize=(14, 4*len(gene_list)))
    if len(gene_list)==1: axes = axes[None,:]
    for i, name in enumerate(gene_list):
        rc = f"rna_{name}"; pc = f"prot_{name}"
        sub = paired_df[[rc,pc,"braaksc","AD"]].dropna()
        ad_lab = sub["AD"].map({0:"NCI",1:"AD"}).fillna("Unknown")

        ax = axes[i,0]
        for grp,col in [("NCI","#2196F3"),("AD","#F44336"),("Unknown","#9E9E9E")]:
            m = ad_lab==grp
            if m.sum(): ax.scatter(sub.loc[m,rc], sub.loc[m,pc], alpha=0.5, s=18, c=col, label=grp)
        rho,pv = spearmanr(sub[rc],sub[pc])
        ax.set_xlabel(f"{name} RNA (logCPM)")
        ax.set_ylabel(f"{name} Protein (log2)")
        ax.set_title(f"{name}: RNA vs Protein")
        ax.legend(fontsize=7)
        ax.text(0.05,0.95,f"ρ={rho:.3f} p={pv:.3f}", transform=ax.transAxes, fontsize=8, va="top")

        axes[i,1].scatter(sub["braaksc"], sub[rc], alpha=0.4, s=18, c="#4CAF50")
        axes[i,1].set_xlabel("Braak"); axes[i,1].set_ylabel(f"{name} RNA"); axes[i,1].set_title(f"{name}: Braak vs RNA")

        axes[i,2].scatter(sub["braaksc"], sub[pc], alpha=0.4, s=18, c="#FF9800")
        axes[i,2].set_xlabel("Braak"); axes[i,2].set_ylabel(f"{name} Protein"); axes[i,2].set_title(f"{name}: Braak vs Protein")

    plt.tight_layout()
    plt.savefig(OUT+"test_J_paired_scatterplots.png", dpi=120)
    plt.close()
    print("  Saved: test_J_paired_scatterplots.png")

# ─────────────────────────────────────────────────────────────────────────
# TEST K: SUMMARY TABLES
# ─────────────────────────────────────────────────────────────────────────
print("\n"+"="*60)
print("=== TEST K: SUMMARY TABLES ===")
print("="*60)

print("\n--- AD effect (composition-adjusted) ---")
print(f"{'Gene':<10} {'RNA β_AD':>10} {'RNA p':>8} {'Prot β_AD':>11} {'Prot p':>8}")
print("-"*45)
for name in GENES:
    r = testB.get(name,{}).get("cells"); p = testC.get(name,{}).get("cells")
    rb = f"{r['beta']:10.4f}" if r else "       N/A"
    rp = f"{r['p']:8.4f}"     if r else "     N/A"
    pb = f"{p['beta']:11.4f}" if p else "        N/A"
    pp = f"{p['p']:8.4f}"     if p else "     N/A"
    print(f"  {name:<8} {rb} {rp} {pb} {pp}")

print("\n--- Braak effect (composition-adjusted) ---")
print(f"{'Gene':<10} {'RNA β_Braak':>12} {'RNA p':>8} {'Prot β_Braak':>13} {'Prot p':>8}")
print("-"*48)
for name in GENES:
    r = testD.get(name,{}).get("rna"); p = testD.get(name,{}).get("prot")
    rb = f"{r['beta']:12.4f}" if r else "         N/A"
    rp = f"{r['p']:8.4f}"     if r else "     N/A"
    pb = f"{p['beta']:13.4f}" if p else "          N/A"
    pp = f"{p['p']:8.4f}"     if p else "     N/A"
    print(f"  {name:<8} {rb} {rp} {pb} {pp}")

print("\n--- Standardized discordance (paired) ---")
print(f"{'Gene':<10} {'β_RNA_z':>9} {'β_Prot_z':>10} {'Δβ':>8} {'D':>8}")
print("-"*42)
for name in GENES:
    e = testE.get(name)
    if not e: print(f"  {name:<8}  N/A"); continue
    print(f"  {name:<8} {e['rna']['beta']:>9.4f} {e['prot']['beta']:>10.4f} {e['delta_beta']:>8.4f} {e['D']:>8.4f}")

print("\n"+"="*60)
print("=== DONE ===")
print("="*60)
