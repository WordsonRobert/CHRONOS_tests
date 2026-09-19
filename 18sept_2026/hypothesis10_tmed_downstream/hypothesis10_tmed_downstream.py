"""
hypothesis10_tmed_downstream.py
H10 / Test 3: TMED → APP / BACE1 / Amyloid
CHRONOS / ROSMAP  —  Sep 2026

Tests
  A  Raw TMED ↔ downstream RNA (Spearman + Pearson, n=15 pairs)
  B  Composition-adjusted TMED → downstream RNA (OLS β, FDR)
  C  Within-diagnosis stratification for key pairs
  D  TMED RNA ↔ pathology (Spearman, primary + secondary + exploratory)
  E  Composition-adjusted amyloid models (RNA → amyloid_load)
  F  TMED protein ↔ amyloid (Spearman + OLS)
  G  Protein ↔ amyloid within NCI / AD
  H  Conditional protein model (protein + RNA → amyloid)
  I  Incremental R² / AIC  (TMED protein beyond covariates)
  J  BACE1 protein ↔ TMED protein  (if BACE1 in proteomics)
  K  Directionality matrix
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────────
DATA    = "/home/wordson22/projects/CHRONOS/DATA/"
OUT     = "/home/wordson22/projects/CHRONOS/analysis/hypothesis10_tmed_downstream/"

import os; os.makedirs(OUT, exist_ok=True)

# ── gene IDs ───────────────────────────────────────────────────────────────
TMED_IDS = {
    "TMED2" : "ENSG00000086598",
    "TMED10": "ENSG00000170348",
    "TMED9" : "ENSG00000184840",
}
TARGET_IDS = {
    "BACE1" : "ENSG00000186318",
    "APP"   : "ENSG00000142192",
    "PSEN1" : "ENSG00000080815",
    "PSEN2" : "ENSG00000164000",
    "ADAM10": "ENSG00000137845",
}
ALL_IDS = {**TMED_IDS, **TARGET_IDS}

CELL_MARKERS = {
    "neuron_score": ["ENSG00000102003","ENSG00000067715","ENSG00000132639","ENSG00000008056","ENSG00000157542"],
    "astro_score" : ["ENSG00000131095","ENSG00000171885"],
    "micro_score" : ["ENSG00000204472","ENSG00000138185"],
    "oligo_score" : ["ENSG00000197971","ENSG00000123560"],
}
CELL_COLS = list(CELL_MARKERS.keys())

PATHOLOGY_OUTCOMES = ["ceradsc","braaksc","cogdx"]
PRIMARY_PATH   = ["ceradsc"]      # CERAD = amyloid burden proxy (1=definite, 4=no AD)
SECONDARY_PATH = ["braaksc"]
EXPLORATORY_PATH = ["cogdx"]

# ── helpers ────────────────────────────────────────────────────────────────

def ols_term(y_s, X_df, term):
    """OLS: return β, SE, p, n for one predictor 'term' in X_df."""
    cols = list(X_df.columns)
    Xm   = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    yv   = y_s.values if hasattr(y_s, "values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(yv))
    n_ok = int(mask.sum())
    if n_ok < len(cols) + 3:
        return dict(beta=np.nan, se=np.nan, p=np.nan, n=n_ok)
    Xm, yv = Xm[mask], yv[mask]
    coeffs, _, _, _ = np.linalg.lstsq(Xm, yv, rcond=None)
    resid  = yv - Xm @ coeffs
    n, p   = len(yv), Xm.shape[1]
    mse    = np.sum(resid**2) / (n - p)
    try:
        cov = mse * np.linalg.inv(Xm.T @ Xm)
    except np.linalg.LinAlgError:
        return dict(beta=np.nan, se=np.nan, p=np.nan, n=n_ok)
    se_all = np.sqrt(np.diag(cov))
    idx    = (["intercept"] + cols).index(term)
    beta, se_b = coeffs[idx], se_all[idx]
    t  = beta / se_b
    pv = 2 * stats.t.sf(abs(t), df=n - p)
    return dict(beta=beta, se=se_b, p=pv, n=n_ok)


def r2_model(y_s, X_df):
    """Return R² for OLS y ~ X_df (complete cases)."""
    cols = list(X_df.columns)
    Xm   = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    yv   = y_s.values if hasattr(y_s, "values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(yv))
    if mask.sum() < len(cols) + 3:
        return np.nan, np.nan
    Xm, yv = Xm[mask], yv[mask]
    coeffs, _, _, _ = np.linalg.lstsq(Xm, yv, rcond=None)
    yhat = Xm @ coeffs
    ss_res = np.sum((yv - yhat)**2)
    ss_tot = np.sum((yv - yv.mean())**2)
    n, p = len(yv), Xm.shape[1]
    r2 = 1 - ss_res / ss_tot
    aic = n * np.log(ss_res / n) + 2 * p
    return r2, aic


def spearman_sub(x, y):
    """Spearman on non-NaN pairs."""
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 5:
        return np.nan, np.nan, int(mask.sum())
    rho, p = stats.spearmanr(x[mask], y[mask])
    return rho, p, int(mask.sum())


def pearson_sub(x, y):
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 5:
        return np.nan, np.nan, int(mask.sum())
    r, p = stats.pearsonr(x[mask], y[mask])
    return r, p, int(mask.sum())


def bh_fdr(pvals):
    """BH FDR; NaN p-values get NaN FDR."""
    out = np.full(len(pvals), np.nan)
    mask = ~np.isnan(pvals)
    if mask.sum() == 0:
        return out
    _, adj, _, _ = multipletests(np.array(pvals)[mask], method="fdr_bh")
    out[mask] = adj
    return out


def section(title):
    print("\n" + "="*60)
    print(f"=== {title} ===")
    print("="*60)


# ══════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════
print("Loading data...")

expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv",   sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")
prot_raw = pd.read_csv(DATA + "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv", index_col=0)
prot_meta= pd.read_csv(DATA + "ROSMAP_assay_proteomics_TMTquantitation_metadata.csv")

# ── clinical numeric ────────────────────────────────────────────────────────
clinical["age_at_visit_max"] = pd.to_numeric(
    clinical["age_at_visit_max"].replace("90+", "90"), errors="coerce")
for col in ["braaksc","msex","pmi","cogdx","ceradsc"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")
clinical["AD"] = np.nan
clinical.loc[clinical["cogdx"] == 1, "AD"] = 0
clinical.loc[clinical["cogdx"] >= 4, "AD"] = 1

COV_BASE = ["age_at_visit_max","msex","pmi"]

# ── RNA cohort ──────────────────────────────────────────────────────────────
print("Building RNA cohort...")
rna_bio  = biospec[["individualID","specimenID"]].drop_duplicates()
rna_link = rna_meta[["specimenID","rnaBatch"]].merge(rna_bio, on="specimenID", how="left")
rna_link = rna_link[rna_link["specimenID"].isin(expr.columns)]
rna_link = rna_link.merge(
    clinical[["individualID","braaksc","msex","pmi","age_at_visit_max",
              "cogdx","AD","ceradsc"]],
    on="individualID", how="left")
rna_link = rna_link.dropna(subset=["individualID"])

# rnaBatch usable?
_rb_vals = rna_link["rnaBatch"].dropna().unique()
if len(_rb_vals) > 1:
    rna_link["rnaBatch_cat"] = pd.Categorical(rna_link["rnaBatch"]).codes.astype(float)
    rna_link.loc[rna_link["rnaBatch"].isna(), "rnaBatch_cat"] = np.nan
    _USE_RNA_BATCH = True
else:
    rna_link["rnaBatch_cat"] = 0.0
    _USE_RNA_BATCH = False
_RNA_BATCH_COLS = ["rnaBatch_cat"] if _USE_RNA_BATCH else []
print(f"  rnaBatch usable: {_USE_RNA_BATCH}")

# cell scores
def cell_score(markers):
    avail = [m for m in markers if m in expr.index]
    if not avail:
        return pd.Series(np.nan, index=expr.columns)
    z = expr.loc[avail].apply(lambda r: (r - r.mean()) / r.std(), axis=1)
    return z.mean(axis=0)

for cname, markers in CELL_MARKERS.items():
    scores = cell_score(markers)
    rna_link[cname] = rna_link["specimenID"].map(scores)

# gene expression columns
for name, eid in ALL_IDS.items():
    if eid in expr.index:
        rna_link[f"rna_{name}"] = rna_link["specimenID"].map(expr.loc[eid])
    else:
        print(f"  WARNING: {name} ({eid}) not found in expression matrix")
        rna_link[f"rna_{name}"] = np.nan

print(f"  RNA samples linked: {len(rna_link)}")

# ── protein cohort ──────────────────────────────────────────────────────────
print("Building protein cohort...")
pm = prot_meta[prot_meta["isAssayControl"] == False].copy()
pm = pm[["specimenID","batchChannel","batch"]].drop_duplicates()
pm["individualID"] = pm["specimenID"].str.split(".").str[-1]
pm = pm[pm["batchChannel"].isin(prot_raw.columns)]
pm = pm.merge(clinical[["individualID","braaksc","msex","pmi","age_at_visit_max",
                         "cogdx","AD","ceradsc"]],
              on="individualID", how="left")

# TMED proteins
tmed_prot_map = {}
for name in TMED_IDS:
    matches = [g for g in prot_raw.index if g.startswith(name + "|")]
    if matches:
        tmed_prot_map[name] = matches[0]
print(f"  TMED proteins: {tmed_prot_map}")

# BACE1 protein (Test J)
bace1_prot_col = None
bace1_matches = [g for g in prot_raw.index if g.startswith("BACE1|") or g.startswith("BACE|")]
if bace1_matches:
    bace1_prot_col = bace1_matches[0]
    print(f"  BACE1 protein found: {bace1_prot_col}")
else:
    print("  BACE1 protein: NOT FOUND in proteomics")

# add protein values to pm
for name, col in tmed_prot_map.items():
    pm[f"prot_{name}"] = pm["batchChannel"].map(prot_raw.loc[col])
if bace1_prot_col:
    pm["prot_BACE1"] = pm["batchChannel"].map(prot_raw.loc[bace1_prot_col])

# aggregate to individual (mean across channels)
prot_batch_mode = pm.groupby("individualID")["batch"].first().reset_index()
prot_df = pm.drop(columns=["batch","batchChannel","specimenID"], errors="ignore").groupby("individualID").mean().reset_index()
prot_df = prot_df.merge(prot_batch_mode, on="individualID")
prot_df["prot_batch_cat"] = pd.Categorical(prot_df["batch"]).codes.astype(float)
PROT_BATCH_COLS = ["prot_batch_cat"]

print(f"  Protein individuals: {len(prot_df)}")

# ── paired cohort ───────────────────────────────────────────────────────────
RNA_COV_CELLS = COV_BASE + _RNA_BATCH_COLS + CELL_COLS

rna_paired = rna_link.dropna(subset=RNA_COV_CELLS).copy()
rna_paired = rna_paired[rna_paired["AD"].isin([0, 1])]

prot_paired = prot_df.dropna(subset=COV_BASE + PROT_BATCH_COLS).copy()

_prot_cols_for_paired = [f"prot_{n}" for n in tmed_prot_map]
if bace1_prot_col:
    _prot_cols_for_paired += ["prot_BACE1"]
paired = rna_paired[["individualID"] + [f"rna_{n}" for n in list(TMED_IDS)+list(TARGET_IDS)]
                    + RNA_COV_CELLS + ["AD","braaksc","ceradsc"]].merge(
    prot_df[["individualID"] + _prot_cols_for_paired + PROT_BATCH_COLS],
    on="individualID", how="inner")
print(f"  Paired (RNA+Prot+cells): {len(paired)}")


# ══════════════════════════════════════════════════════════════
# TEST A: Raw TMED ↔ downstream RNA (Spearman + Pearson)
# ══════════════════════════════════════════════════════════════
section("TEST A: RAW TMED ↔ DOWNSTREAM RNA")

rows_a = []
for tmed in TMED_IDS:
    for tgt in TARGET_IDS:
        x = rna_link[f"rna_{tmed}"].values.astype(float)
        y = rna_link[f"rna_{tgt}"].values.astype(float)
        rho, p_sp, n = spearman_sub(x, y)
        r,   p_pe, _ = pearson_sub(x, y)
        rows_a.append(dict(TMED=tmed, Target=tgt,
                           rho=rho, p_spearman=p_sp,
                           r_pearson=r, p_pearson=p_pe, n=n))

df_a = pd.DataFrame(rows_a)
df_a["FDR"] = bh_fdr(df_a["p_spearman"].values)

print(f"\n{'TMED':<8} {'Target':<8} {'rho':>7} {'p_sp':>9} {'FDR':>9} {'r':>7} {'p_pe':>9} {'n':>5}")
print("-"*65)
for _, r in df_a.iterrows():
    print(f"{r['TMED']:<8} {r['Target']:<8} {r['rho']:>7.4f} {r['p_spearman']:>9.4f} "
          f"{r['FDR']:>9.4f} {r['r_pearson']:>7.4f} {r['p_pearson']:>9.4f} {int(r['n']):>5}")

df_a.to_csv(OUT + "h10_A_raw_rna_correlations.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST B: Composition-adjusted TMED → downstream RNA
# ══════════════════════════════════════════════════════════════
section("TEST B: COMPOSITION-ADJUSTED TMED → DOWNSTREAM RNA")

rows_b = []
sub_b = rna_link.dropna(subset=RNA_COV_CELLS + [f"rna_{t}" for t in TMED_IDS]).copy()

for tmed in TMED_IDS:
    for tgt in TARGET_IDS:
        y_col = f"rna_{tgt}"
        x_col = f"rna_{tmed}"
        sub   = sub_b.dropna(subset=[y_col, x_col])
        X_df  = sub[[x_col] + RNA_COV_CELLS].rename(columns={x_col: tmed})
        res   = ols_term(sub[y_col], X_df, tmed)
        rows_b.append(dict(TMED=tmed, Target=tgt, **res))

df_b = pd.DataFrame(rows_b)
df_b["FDR"] = bh_fdr(df_b["p"].values)

print(f"\n{'TMED':<8} {'Target':<8} {'β':>9} {'SE':>8} {'p':>9} {'FDR':>9} {'n':>5}")
print("-"*57)
for _, r in df_b.iterrows():
    print(f"{r['TMED']:<8} {r['Target']:<8} {r['beta']:>9.4f} {r['se']:>8.4f} "
          f"{r['p']:>9.4f} {r['FDR']:>9.4f} {int(r['n']):>5}")

df_b.to_csv(OUT + "h10_B_adjusted_rna_models.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST C: Within-diagnosis stratification (key pairs)
# ══════════════════════════════════════════════════════════════
section("TEST C: WITHIN-DIAGNOSIS STRATIFICATION (KEY PAIRS)")

KEY_PAIRS = [("TMED2","BACE1"), ("TMED2","APP"), ("TMED9","BACE1")]

rows_c = []
for tmed, tgt in KEY_PAIRS:
    for group, label in [(None,"ALL"), (0,"NCI"), (1,"AD")]:
        if group is None:
            sub = rna_link
        else:
            sub = rna_link[rna_link["AD"] == group]
        x = sub[f"rna_{tmed}"].values.astype(float)
        y = sub[f"rna_{tgt}"].values.astype(float)
        rho, p, n = spearman_sub(x, y)
        rows_c.append(dict(TMED=tmed, Target=tgt, Group=label, rho=rho, p=p, n=n))

df_c = pd.DataFrame(rows_c)
print(f"\n{'TMED':<8} {'Target':<8} {'Group':<6} {'rho':>8} {'p':>9} {'n':>5}")
print("-"*47)
for _, r in df_c.iterrows():
    print(f"{r['TMED']:<8} {r['Target']:<8} {r['Group']:<6} {r['rho']:>8.4f} {r['p']:>9.4f} {int(r['n']):>5}")

df_c.to_csv(OUT + "h10_C_stratified.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST D: TMED RNA ↔ Pathology (Spearman)
# ══════════════════════════════════════════════════════════════
section("TEST D: TMED RNA ↔ PATHOLOGY (Spearman)")

rows_d = []
for tmed in TMED_IDS:
    for path in PATHOLOGY_OUTCOMES:
        priority = ("primary" if path in PRIMARY_PATH
                    else "secondary" if path in SECONDARY_PATH
                    else "exploratory")
        x = rna_link[f"rna_{tmed}"].values.astype(float)
        y = rna_link[path].values.astype(float)
        rho, p, n = spearman_sub(x, y)
        rows_d.append(dict(TMED=tmed, Outcome=path, Priority=priority, rho=rho, p=p, n=n))

df_d = pd.DataFrame(rows_d)
df_d["FDR"] = bh_fdr(df_d["p"].values)

print(f"\n{'TMED':<8} {'Outcome':<14} {'Priority':<12} {'rho':>8} {'p':>9} {'FDR':>9} {'n':>5}")
print("-"*68)
for _, r in df_d.iterrows():
    print(f"{r['TMED']:<8} {r['Outcome']:<14} {r['Priority']:<12} "
          f"{r['rho']:>8.4f} {r['p']:>9.4f} {r['FDR']:>9.4f} {int(r['n']):>5}")

df_d.to_csv(OUT + "h10_D_rna_pathology.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST E: Composition-adjusted amyloid models (RNA)
# ══════════════════════════════════════════════════════════════
section("TEST E: COMPOSITION-ADJUSTED AMYLOID MODELS (RNA)")

sub_e = rna_link.dropna(subset=RNA_COV_CELLS + ["ceradsc"]).copy()

rows_e = []
for tmed in TMED_IDS:
    x_col = f"rna_{tmed}"
    sub   = sub_e.dropna(subset=[x_col])
    X_df  = sub[[x_col] + RNA_COV_CELLS].rename(columns={x_col: tmed})
    res   = ols_term(sub["ceradsc"], X_df, tmed)
    r2, _ = r2_model(sub["ceradsc"], X_df)
    rows_e.append(dict(TMED=tmed, Outcome="ceradsc", r2=r2, **res))

df_e = pd.DataFrame(rows_e)
df_e["FDR"] = bh_fdr(df_e["p"].values)

print(f"\n{'TMED':<8} {'β':>9} {'SE':>8} {'p':>9} {'FDR':>9} {'R²':>7} {'n':>5}")
print("-"*57)
for _, r in df_e.iterrows():
    print(f"{r['TMED']:<8} {r['beta']:>9.4f} {r['se']:>8.4f} {r['p']:>9.4f} "
          f"{r['FDR']:>9.4f} {r['r2']:>7.4f} {int(r['n']):>5}")

df_e.to_csv(OUT + "h10_E_rna_amyloid_adjusted.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST F: TMED protein ↔ amyloid (Spearman + OLS)
# ══════════════════════════════════════════════════════════════
section("TEST F: TMED PROTEIN ↔ AMYLOID")

sub_f = prot_df.dropna(subset=COV_BASE + PROT_BATCH_COLS + ["ceradsc"]).copy()

rows_f_sp = []
rows_f_ols = []

for name in tmed_prot_map:
    pc = f"prot_{name}"
    # Spearman raw
    x = sub_f[pc].values.astype(float)
    y = sub_f["ceradsc"].values.astype(float)
    rho, p, n = spearman_sub(x, y)
    rows_f_sp.append(dict(TMED=name, rho=rho, p=p, n=n, model="raw"))

    # OLS simple (demographics + protein batch)
    sub   = sub_f.dropna(subset=[pc])
    X_df  = sub[[pc] + COV_BASE + PROT_BATCH_COLS].rename(columns={pc: name})
    res   = ols_term(sub["ceradsc"], X_df, name)
    r2, _ = r2_model(sub["ceradsc"], X_df)
    rows_f_ols.append(dict(TMED=name, model="demo+batch", r2=r2, **res))

    # OLS + cells (if cells available in paired)
    sub_fc = paired.dropna(subset=[pc, "ceradsc"] + RNA_COV_CELLS).copy()
    if len(sub_fc) >= 10:
        X_df2 = sub_fc[[pc] + COV_BASE + PROT_BATCH_COLS + CELL_COLS].rename(columns={pc: name})
        res2  = ols_term(sub_fc["ceradsc"], X_df2, name)
        r2b, _= r2_model(sub_fc["ceradsc"], X_df2)
        rows_f_ols.append(dict(TMED=name, model="demo+batch+cells", r2=r2b, **res2))

df_f_sp  = pd.DataFrame(rows_f_sp)
df_f_ols = pd.DataFrame(rows_f_ols)

print("\n  Spearman (raw):")
print(f"  {'TMED':<8} {'rho':>8} {'p':>9} {'n':>5}")
for _, r in df_f_sp.iterrows():
    print(f"  {r['TMED']:<8} {r['rho']:>8.4f} {r['p']:>9.4f} {int(r['n']):>5}")

print("\n  OLS (adjusted):")
print(f"  {'TMED':<8} {'Model':<22} {'β':>9} {'p':>9} {'R²':>7} {'n':>5}")
for _, r in df_f_ols.iterrows():
    print(f"  {r['TMED']:<8} {r['model']:<22} {r['beta']:>9.4f} {r['p']:>9.4f} {r['r2']:>7.4f} {int(r['n']):>5}")

df_f_sp.to_csv(OUT  + "h10_F_protein_amyloid_spearman.csv", index=False)
df_f_ols.to_csv(OUT + "h10_F_protein_amyloid_ols.csv",     index=False)


# ══════════════════════════════════════════════════════════════
# TEST G: Protein ↔ amyloid within NCI / AD
# ══════════════════════════════════════════════════════════════
section("TEST G: PROTEIN ↔ AMYLOID WITHIN NCI / AD")

rows_g = []
for name in tmed_prot_map:
    pc = f"prot_{name}"
    for group, label in [(None,"ALL"), (0,"NCI"), (1,"AD")]:
        if group is None:
            sub = prot_df
        else:
            sub = prot_df[prot_df["AD"] == group]
        x = sub[pc].values.astype(float)
        y = sub["ceradsc"].values.astype(float)
        rho, p, n = spearman_sub(x, y)
        rows_g.append(dict(TMED=name, Group=label, rho=rho, p=p, n=n))

df_g = pd.DataFrame(rows_g)
print(f"\n{'TMED':<8} {'Group':<6} {'rho':>8} {'p':>9} {'n':>5}")
print("-"*40)
for _, r in df_g.iterrows():
    print(f"{r['TMED']:<8} {r['Group']:<6} {r['rho']:>8.4f} {r['p']:>9.4f} {int(r['n']):>5}")

df_g.to_csv(OUT + "h10_G_protein_amyloid_stratified.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST H: Conditional protein model (protein + RNA → amyloid)
# ══════════════════════════════════════════════════════════════
section("TEST H: CONDITIONAL PROTEIN MODEL (protein + RNA → amyloid)")

rows_h = []
sub_h = paired.dropna(subset=["ceradsc"] + COV_BASE + PROT_BATCH_COLS + CELL_COLS).copy()

for name in tmed_prot_map:
    pc   = f"prot_{name}"
    rc   = f"rna_{name}"
    sub  = sub_h.dropna(subset=[pc, rc])

    # protein alone (+ demo + batch + cells)
    X_pa = sub[[pc] + COV_BASE + PROT_BATCH_COLS + CELL_COLS].rename(columns={pc: name})
    res_p= ols_term(sub["ceradsc"], X_pa, name)
    r2_p, aic_p = r2_model(sub["ceradsc"], X_pa)

    # protein + RNA
    X_pr  = sub[[pc, rc] + COV_BASE + PROT_BATCH_COLS + CELL_COLS].rename(columns={pc: name, rc: f"{name}_rna"})
    res_pc= ols_term(sub["ceradsc"], X_pr, name)
    r2_pr, aic_pr = r2_model(sub["ceradsc"], X_pr)

    rows_h.append(dict(TMED=name, model="protein+cov",
                       beta_prot=res_p["beta"], p_prot=res_p["p"],
                       r2=r2_p, aic=aic_p, n=res_p["n"]))
    rows_h.append(dict(TMED=name, model="protein+RNA+cov",
                       beta_prot=res_pc["beta"], p_prot=res_pc["p"],
                       r2=r2_pr, aic=aic_pr, n=res_pc["n"]))

df_h = pd.DataFrame(rows_h)
print(f"\n{'TMED':<8} {'Model':<22} {'β_prot':>9} {'p_prot':>9} {'R²':>7} {'n':>5}")
print("-"*60)
for _, r in df_h.iterrows():
    print(f"{r['TMED']:<8} {r['model']:<22} {r['beta_prot']:>9.4f} "
          f"{r['p_prot']:>9.4f} {r['r2']:>7.4f} {int(r['n']):>5}")

df_h.to_csv(OUT + "h10_H_conditional_protein.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST I: Incremental R² / ΔAIC (TMED protein beyond covariates)
# ══════════════════════════════════════════════════════════════
section("TEST I: INCREMENTAL R² / ΔAIC")

rows_i = []
sub_i = prot_df.dropna(subset=COV_BASE + PROT_BATCH_COLS + ["ceradsc"]).copy()

# null model
X_null = sub_i[COV_BASE + PROT_BATCH_COLS]
r2_null, aic_null = r2_model(sub_i["ceradsc"], X_null)

print(f"\nNull (demo+batch): R²={r2_null:.4f}, AIC={aic_null:.2f}, n={len(sub_i)}")
print(f"\n{'TMED':<8} {'R²_null':>9} {'R²_TMED':>9} {'ΔR²':>8} {'AIC_null':>10} {'AIC_TMED':>10} {'ΔAIC':>8} {'β':>9} {'p':>9} {'n':>5}")
print("-"*90)

for name in tmed_prot_map:
    pc  = f"prot_{name}"
    sub = sub_i.dropna(subset=[pc])
    X0  = sub[COV_BASE + PROT_BATCH_COLS]
    Xt  = sub[[pc] + COV_BASE + PROT_BATCH_COLS].rename(columns={pc: name})
    r2_0, aic_0 = r2_model(sub["ceradsc"], X0)
    r2_t, aic_t = r2_model(sub["ceradsc"], Xt)
    res = ols_term(sub["ceradsc"], Xt, name)
    rows_i.append(dict(TMED=name, r2_null=r2_0, r2_tmed=r2_t,
                       delta_r2=r2_t - r2_0, aic_null=aic_0, aic_tmed=aic_t,
                       delta_aic=aic_t - aic_0, beta=res["beta"], p=res["p"], n=res["n"]))
    print(f"{name:<8} {r2_0:>9.4f} {r2_t:>9.4f} {r2_t-r2_0:>8.4f} "
          f"{aic_0:>10.2f} {aic_t:>10.2f} {aic_t-aic_0:>8.2f} "
          f"{res['beta']:>9.4f} {res['p']:>9.4f} {int(res['n']):>5}")

df_i = pd.DataFrame(rows_i)
df_i.to_csv(OUT + "h10_I_incremental_r2.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST J: BACE1 protein ↔ TMED protein
# ══════════════════════════════════════════════════════════════
section("TEST J: BACE1 PROTEIN ↔ TMED PROTEIN")

if bace1_prot_col:
    sub_j = prot_df.dropna(subset=COV_BASE + PROT_BATCH_COLS + ["prot_BACE1"]).copy()
    rows_j = []
    for name in tmed_prot_map:
        pc = f"prot_{name}"
        sub = sub_j.dropna(subset=[pc])

        # Spearman
        rho, p, n = spearman_sub(sub[pc].values.astype(float),
                                 sub["prot_BACE1"].values.astype(float))
        # OLS: BACE1_prot ~ TMED_prot + covariates
        X_df = sub[[pc] + COV_BASE + PROT_BATCH_COLS].rename(columns={pc: name})
        res  = ols_term(sub["prot_BACE1"], X_df, name)
        rows_j.append(dict(TMED=name, rho=rho, p_sp=p, n_sp=n,
                           beta=res["beta"], p_ols=res["p"], n_ols=res["n"]))

    df_j = pd.DataFrame(rows_j)
    print(f"\n{'TMED':<8} {'rho':>8} {'p_sp':>9} {'β':>9} {'p_ols':>9} {'n':>5}")
    print("-"*48)
    for _, r in df_j.iterrows():
        print(f"{r['TMED']:<8} {r['rho']:>8.4f} {r['p_sp']:>9.4f} "
              f"{r['beta']:>9.4f} {r['p_ols']:>9.4f} {int(r['n_ols']):>5}")

    # conditional: BACE1_prot ~ TMED_prot + TMED_RNA + cov (paired)
    print("\n  Conditional (protein + RNA) in paired:")
    for name in tmed_prot_map:
        pc = f"prot_{name}"
        rc = f"rna_{name}"
        sub = paired.dropna(subset=[pc, rc, "prot_BACE1"] + COV_BASE + PROT_BATCH_COLS)
        if len(sub) < 10:
            print(f"  {name}: insufficient n={len(sub)}")
            continue
        X_df = sub[[pc, rc] + COV_BASE + PROT_BATCH_COLS].rename(
            columns={pc: name, rc: f"{name}_rna"})
        res_p = ols_term(sub["prot_BACE1"], X_df, name)
        res_r = ols_term(sub["prot_BACE1"], X_df, f"{name}_rna")
        print(f"  {name}: β_prot={res_p['beta']:.4f} p={res_p['p']:.4f} | "
              f"β_rna={res_r['beta']:.4f} p={res_r['p']:.4f}  n={res_p['n']}")

    df_j.to_csv(OUT + "h10_J_bace1_protein.csv", index=False)
else:
    print("  BACE1 protein not found — Test J skipped.")
    df_j = pd.DataFrame()


# ══════════════════════════════════════════════════════════════
# TEST K: Directionality matrix
# ══════════════════════════════════════════════════════════════
section("TEST K: DIRECTIONALITY MATRIX")

# Helper: sign + significance label
def sig_label(beta, p, fdr=None):
    if np.isnan(beta) or np.isnan(p):
        return "n/a"
    sign  = "−" if beta < 0 else "+"
    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    return f"{sign}({stars})"

relationships = [
    ("TMED2",  "BACE1",  "−"),
    ("TMED2",  "APP",    "−"),
    ("TMED9",  "BACE1",  "+"),
    ("TMED2",  "ceradsc","−"),
    ("TMED10", "ceradsc","−"),
    ("TMED9",  "ceradsc","?"),
]

# pull from earlier results
rna_raw_lookup   = df_a.set_index(["TMED","Target"])
rna_adj_lookup   = df_b.set_index(["TMED","Target"])
path_lookup      = df_d.set_index(["TMED","Outcome"])
prot_ols_lookup  = df_f_ols[df_f_ols["model"]=="demo+batch"].set_index("TMED")

print(f"\n{'Relationship':<26} {'Expected':<10} {'RNA raw':>10} {'RNA adj':>10} {'Prot adj':>10}")
print("-"*68)

rows_k = []
for (tmed, target, expected) in relationships:
    # RNA raw
    try:
        row = rna_raw_lookup.loc[(tmed, target)]
        rna_raw = sig_label(row["rho"], row["p_spearman"])
    except:
        rna_raw = "n/a"
    # RNA adj
    try:
        row = rna_adj_lookup.loc[(tmed, target)]
        rna_adj = sig_label(row["beta"], row["p"])
    except:
        # target is a pathology outcome, not a gene
        try:
            row = path_lookup.loc[(tmed, target)]
            rna_adj = sig_label(row["rho"], row["p"])
        except:
            rna_adj = "n/a"
    # Prot adj
    if target == "ceradsc":
        try:
            row = prot_ols_lookup.loc[tmed]
            prot_adj = sig_label(row["beta"], row["p"])
        except:
            prot_adj = "n/a"
    else:
        prot_adj = "n/a"

    rel_str = f"{tmed} → {target}"
    print(f"{rel_str:<26} {expected:<10} {rna_raw:>10} {rna_adj:>10} {prot_adj:>10}")
    rows_k.append(dict(TMED=tmed, Target=target, Expected=expected,
                       RNA_raw=rna_raw, RNA_adj=rna_adj, Prot_adj=prot_adj))

df_k = pd.DataFrame(rows_k)
df_k.to_csv(OUT + "h10_K_directionality_matrix.csv", index=False)


# ══════════════════════════════════════════════════════════════
# SAVE COMBINED SUMMARY
# ══════════════════════════════════════════════════════════════
section("DONE")
print(f"\nOutputs written to: {OUT}")
print("  h10_A_raw_rna_correlations.csv")
print("  h10_B_adjusted_rna_models.csv")
print("  h10_C_stratified.csv")
print("  h10_D_rna_pathology.csv")
print("  h10_E_rna_amyloid_adjusted.csv")
print("  h10_F_protein_amyloid_spearman.csv")
print("  h10_F_protein_amyloid_ols.csv")
print("  h10_G_protein_amyloid_stratified.csv")
print("  h10_H_conditional_protein.csv")
print("  h10_I_incremental_r2.csv")
print("  h10_J_bace1_protein.csv  (if BACE1 protein found)")
print("  h10_K_directionality_matrix.csv")
