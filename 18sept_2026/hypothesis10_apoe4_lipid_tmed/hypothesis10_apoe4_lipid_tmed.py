"""
Test 6: APOE4 → Lipid → TMED bridge
CHRONOS / ROSMAP

Tests:
  A. APOE4 → lipid genes (RNA: SGMS1/2, SMPD1/3, CERT1)
  B. APOE4 → TMED genes (TMED2/10/9)
  C. Lipid gene ↔ TMED RNA (15 pairwise Spearman)
  D. APOE4 × lipid-gene interaction → TMED (15 pairs)
  E. Lipidomics SM species ↔ TMED RNA (all SM × 3 TMED, BH-FDR)
  F. SM-C18 species ↔ TMED (if annotated)
  G. Nested models: APOE4 / SM-C18 / both → TMED2

Outputs: ~/projects/CHRONOS/analysis/hypothesis10_apoe4_lipid_tmed/
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import os, warnings
warnings.filterwarnings("ignore")

DATA    = "/home/wordson22/projects/CHRONOS/DATA/"
RESULTS = "/home/wordson22/projects/CHRONOS/results/"
OUT     = "/home/wordson22/projects/CHRONOS/analysis/hypothesis10_apoe4_lipid_tmed/"
os.makedirs(OUT, exist_ok=True)

def section(title):
    print("\n" + "="*60)
    print(f"=== {title} ===")
    print("="*60)

# ── helpers ───────────────────────────────────────────────────
def spearman_safe(x, y):
    mask = ~(np.isnan(x.astype(float)) | np.isnan(y.astype(float)))
    if mask.sum() < 10:
        return np.nan, np.nan, int(mask.sum())
    r, p = stats.spearmanr(x[mask].astype(float), y[mask].astype(float))
    return float(r), float(p), int(mask.sum())

def ols_fit(y, X):
    """OLS returning coeffs, se, p per term. X includes intercept col."""
    y = np.array(y, dtype=float)
    X = np.array(X, dtype=float)
    mask = ~(np.isnan(y) | np.isnan(X).any(axis=1))
    y, X = y[mask], X[mask]
    n, k = X.shape
    if n < k + 3:
        return None, int(n)
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coeffs
    mse   = np.sum(resid**2) / (n - k)
    cov   = mse * np.linalg.pinv(X.T @ X)
    se    = np.sqrt(np.diag(cov))
    t     = coeffs / np.where(se > 1e-12, se, np.nan)
    p     = [2*(1 - stats.t.cdf(abs(ti), df=n-k)) if not np.isnan(ti) else np.nan for ti in t]
    return dict(coeffs=coeffs, se=se, p=p, n=int(n)), int(n)

# ── 1. LOAD EXPRESSION ────────────────────────────────────────
print("Loading expression...")
expr = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
import json
with open(DATA + "ensg_to_symbol.json") as f:
    ensg_to_symbol = json.load(f)
expr.index = [ensg_to_symbol.get(i, i) for i in expr.index]
expr = expr[~expr.index.duplicated(keep="first")]

# map specimens → individuals
bio = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_bio = bio[bio["assay"]=="rnaSeq"].dropna(subset=["individualID"])
spec2ind = dict(zip(rna_bio["specimenID"], rna_bio["individualID"]))
valid = {c: spec2ind[c] for c in expr.columns if c in spec2ind}
expr = expr[[c for c in expr.columns if c in valid]].copy()
expr.columns = [valid[c] for c in expr.columns]
print(f"Expression: {expr.shape}")

# ── 2. LOAD CLINICAL ──────────────────────────────────────────
print("Loading clinical...")
clin = pd.read_csv(RESULTS + "master_covariates.csv")
shared = list(set(expr.columns) & set(clin["individualID"]))
expr   = expr[shared]
clin_s = clin[clin["individualID"].isin(shared)].set_index("individualID")
# clean censored age values
for c in clin_s.columns:
    if "age" in c.lower():
        clin_s[c] = clin_s[c].astype(str).str.replace(r'\+.*', '', regex=True)
        clin_s[c] = pd.to_numeric(clin_s[c], errors="coerce")
print(f"Shared RNA individuals: {len(shared)}")

# age column
age_col = next((c for c in clin_s.columns if "age" in c.lower()), None)
COV = [c for c in [age_col, "msex", "pmi", "braaksc", "cogdx"] if c and c in clin_s.columns]
print(f"Covariates: {COV}")

# APOE4 carrier flag
apoe4_col = next((c for c in clin_s.columns if "apoe4" in c.lower()), None)
if apoe4_col is None:
    # try to derive from apoe_genotype
    gt_col = next((c for c in clin_s.columns if "apoe" in c.lower() and "genotype" in c.lower()), None)
    if gt_col:
        clin_s["apoe4_carrier"] = clin_s[gt_col].astype(str).str.contains("4").astype(float)
        apoe4_col = "apoe4_carrier"
print(f"APOE4 column: {apoe4_col}")

# gene lists
LIPID_GENES = ["SGMS1","SGMS2","SMPD1","SMPD3","CERT1"]
TMED_GENES  = ["TMED2","TMED10","TMED9"]
ALL_GENES   = LIPID_GENES + TMED_GENES

# check availability
for g in ALL_GENES:
    print(f"  {g}: {'found' if g in expr.index else 'MISSING'}")

# pull expression into clin_s
for g in ALL_GENES:
    if g in expr.index:
        clin_s[f"rna_{g}"] = expr.loc[g].reindex(clin_s.index)

# ── 3. LOAD LIPIDOMICS ────────────────────────────────────────
print("\nLoading brain lipidomics...")
lip_raw = pd.read_csv(DATA + "ROSMAP_Brain_Lipidomic data_batch correction_combined Pos and Neg.csv",
                      index_col=0)
print(f"  Lipidomics raw: {lip_raw.shape}  (rows=biospecimen, cols=lipids+batch)")

# biospecimen → individual mapping (lipidomics-specific)
lip_bio = pd.read_csv(DATA + "ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv")
print(f"  Lipidomics biospec cols: {list(lip_bio.columns)}")

# the specimenID in lipidomics metadata
spec_col = "specimenID" if "specimenID" in lip_bio.columns else lip_bio.columns[1]
ind_col  = "individualID" if "individualID" in lip_bio.columns else lip_bio.columns[0]
lip_spec2ind = dict(zip(lip_bio[spec_col], lip_bio[ind_col]))

# rows of lip_raw are biospecimen IDs — map to individual
lip_raw["individualID"] = lip_raw.index.map(lip_spec2ind)
lip_raw = lip_raw.dropna(subset=["individualID"])
lip_raw = lip_raw.drop_duplicates(subset=["individualID"])
lip_raw = lip_raw.set_index("individualID")

# drop non-lipid columns
drop_cols = [c for c in ["sequence order","batch"] if c in lip_raw.columns]
lip_raw   = lip_raw.drop(columns=drop_cols, errors="ignore")
lip_raw   = lip_raw.apply(pd.to_numeric, errors="coerce")
print(f"  Lipidomics individuals: {len(lip_raw)}, lipid species: {lip_raw.shape[1]}")

# SM species
sm_cols = [c for c in lip_raw.columns if c.upper().startswith("SM")]
print(f"  SM species: {len(sm_cols)}")

# overlap with RNA
shared_lip = list(set(lip_raw.index) & set(clin_s.index))
print(f"  RNA + lipidomics overlap: {len(shared_lip)}")

lip_s   = lip_raw.loc[shared_lip]
clin_sl = clin_s.loc[shared_lip]   # clinical subset with lipidomics

# ── TEST A: APOE4 → LIPID GENES ──────────────────────────────
section("TEST A: APOE4 → LIPID GENES (RNA)")

rows_a = []
if apoe4_col and apoe4_col in clin_s.columns:
    for g in LIPID_GENES:
        col = f"rna_{g}"
        if col not in clin_s.columns: continue
        x = clin_s[apoe4_col].values.astype(float)
        y = clin_s[col].values.astype(float)
        r, p, n = spearman_safe(x, y)
        rows_a.append(dict(gene=g, rho=r, p=p, n=n))
    df_a = pd.DataFrame(rows_a)
    if len(df_a):
        _, fdr, _, _ = multipletests(df_a["p"].fillna(1), method="fdr_bh")
        df_a["FDR"] = fdr
        print(f"\n{'Gene':<12} {'rho':>8} {'p':>10} {'FDR':>10} {'n':>6}")
        print("-"*50)
        for _, r in df_a.iterrows():
            print(f"{r['gene']:<12} {r['rho']:>8.4f} {r['p']:>10.4f} {r['FDR']:>10.4f} {int(r['n']):>6}")
        df_a.to_csv(OUT + "h6_A_apoe4_lipid_genes.csv", index=False)
else:
    print("  APOE4 column not found — skipping Test A")
    df_a = pd.DataFrame()

# ── TEST B: APOE4 → TMED GENES ───────────────────────────────
section("TEST B: APOE4 → TMED GENES (RNA)")

rows_b = []
if apoe4_col and apoe4_col in clin_s.columns:
    for g in TMED_GENES:
        col = f"rna_{g}"
        if col not in clin_s.columns: continue
        x = clin_s[apoe4_col].values.astype(float)
        y = clin_s[col].values.astype(float)
        r, p, n = spearman_safe(x, y)
        rows_b.append(dict(gene=g, rho=r, p=p, n=n))
    df_b = pd.DataFrame(rows_b)
    if len(df_b):
        _, fdr, _, _ = multipletests(df_b["p"].fillna(1), method="fdr_bh")
        df_b["FDR"] = fdr
        print(f"\n{'Gene':<12} {'rho':>8} {'p':>10} {'FDR':>10} {'n':>6}")
        print("-"*50)
        for _, r in df_b.iterrows():
            print(f"{r['gene']:<12} {r['rho']:>8.4f} {r['p']:>10.4f} {r['FDR']:>10.4f} {int(r['n']):>6}")
        df_b.to_csv(OUT + "h6_B_apoe4_tmed_genes.csv", index=False)
else:
    print("  APOE4 column not found — skipping Test B")
    df_b = pd.DataFrame()

# ── TEST C: LIPID GENE ↔ TMED RNA (15 pairs) ─────────────────
section("TEST C: LIPID GENE ↔ TMED RNA (15 pairwise Spearman)")

rows_c = []
for lg in LIPID_GENES:
    for tg in TMED_GENES:
        lc, tc = f"rna_{lg}", f"rna_{tg}"
        if lc not in clin_s.columns or tc not in clin_s.columns: continue
        x = clin_s[lc].values.astype(float)
        y = clin_s[tc].values.astype(float)
        r, p, n = spearman_safe(x, y)
        rows_c.append(dict(lipid_gene=lg, tmed_gene=tg, rho=r, p=p, n=n))

df_c = pd.DataFrame(rows_c)
if len(df_c):
    _, fdr, _, _ = multipletests(df_c["p"].fillna(1), method="fdr_bh")
    df_c["FDR"] = fdr
    print(f"\n{'Lipid gene':<12} {'TMED gene':<12} {'rho':>8} {'p':>10} {'FDR':>10} {'n':>6}")
    print("-"*58)
    for _, r in df_c.iterrows():
        print(f"{r['lipid_gene']:<12} {r['tmed_gene']:<12} {r['rho']:>8.4f} {r['p']:>10.4f} {r['FDR']:>10.4f} {int(r['n']):>6}")
    df_c.to_csv(OUT + "h6_C_lipid_tmed_pairs.csv", index=False)

# ── TEST D: APOE4 × LIPID-GENE INTERACTION → TMED ────────────
section("TEST D: APOE4 × LIPID-GENE INTERACTION → TMED")

rows_d = []
if apoe4_col and apoe4_col in clin_s.columns:
    for lg in LIPID_GENES:
        for tg in TMED_GENES:
            lc, tc = f"rna_{lg}", f"rna_{tg}"
            if lc not in clin_s.columns or tc not in clin_s.columns: continue

            sub = clin_s[[tc, lc, apoe4_col] + COV].dropna()
            if len(sub) < 20: continue

            y      = sub[tc].values.astype(float)
            lip_v  = sub[lc].values.astype(float)
            ap_v   = sub[apoe4_col].values.astype(float)
            inter  = lip_v * ap_v
            cov_v  = sub[COV].values.astype(float)

            # center lip and apoe4 for interaction
            lip_c = lip_v - lip_v.mean()
            ap_c  = ap_v  - ap_v.mean()
            inter_c = lip_c * ap_c

            X = np.column_stack([np.ones(len(y)), lip_c, ap_c, inter_c, cov_v])
            res, n = ols_fit(y, X)
            if res is None: continue

            # coefficient index: 0=intercept,1=lip,2=apoe4,3=interaction,4+=cov
            b_inter = res["coeffs"][3]
            p_inter = res["p"][3]
            rows_d.append(dict(lipid_gene=lg, tmed_gene=tg,
                               beta_interaction=b_inter, p_interaction=p_inter, n=n))

    df_d = pd.DataFrame(rows_d)
    if len(df_d):
        _, fdr, _, _ = multipletests(df_d["p_interaction"].fillna(1), method="fdr_bh")
        df_d["FDR"] = fdr
        print(f"\n{'Lipid gene':<12} {'TMED gene':<12} {'β_interact':>12} {'p':>10} {'FDR':>10} {'n':>6}")
        print("-"*66)
        for _, r in df_d.iterrows():
            print(f"{r['lipid_gene']:<12} {r['tmed_gene']:<12} {r['beta_interaction']:>12.4f} "
                  f"{r['p_interaction']:>10.4f} {r['FDR']:>10.4f} {int(r['n']):>6}")
        df_d.to_csv(OUT + "h6_D_interaction.csv", index=False)
else:
    print("  APOE4 column not found — skipping Test D")
    df_d = pd.DataFrame()

# composite scores
print("\nBuilding composite scores...")
def zscore_col(s):
    m, sd = s.mean(), s.std()
    return (s - m) / sd if sd > 1e-8 else s * 0

for g in ALL_GENES:
    col = f"rna_{g}"
    if col in clin_s.columns:
        clin_s[f"z_{g}"] = zscore_col(clin_s[col])

# SM synthesis, breakdown, TMED scores
avail_syn = [g for g in ["SGMS1","SGMS2"] if f"z_{g}" in clin_s.columns]
avail_brk = [g for g in ["SMPD1","SMPD3"] if f"z_{g}" in clin_s.columns]
avail_tm2 = [g for g in ["TMED2","TMED10"] if f"z_{g}" in clin_s.columns]
avail_tm3 = [g for g in ["TMED2","TMED10","TMED9"] if f"z_{g}" in clin_s.columns]

if avail_syn:
    clin_s["SM_synthesis"] = clin_s[[f"z_{g}" for g in avail_syn]].mean(axis=1)
if avail_brk:
    clin_s["SM_breakdown"] = clin_s[[f"z_{g}" for g in avail_brk]].mean(axis=1)
if "CERT1" in [g for g in LIPID_GENES if f"z_{g}" in clin_s.columns]:
    clin_s["CERT1_z"] = clin_s["z_CERT1"]
if avail_tm2:
    clin_s["TMED_score2"] = clin_s[[f"z_{g}" for g in avail_tm2]].mean(axis=1)
if avail_tm3:
    clin_s["TMED_score3"] = clin_s[[f"z_{g}" for g in avail_tm3]].mean(axis=1)

print(f"  SM_synthesis ({'+'.join(avail_syn)}), SM_breakdown ({'+'.join(avail_brk)})")
print(f"  TMED_score2 ({'+'.join(avail_tm2)}), TMED_score3 ({'+'.join(avail_tm3)})")

# composite ↔ TMED score
rows_comp = []
for lipsc, lipname in [("SM_synthesis","SM_synthesis"),("SM_breakdown","SM_breakdown"),("CERT1_z","CERT1")]:
    if lipsc not in clin_s.columns: continue
    for tmsc, tmname in [("TMED_score2","TMED_score2"),("TMED_score3","TMED_score3")]:
        if tmsc not in clin_s.columns: continue
        x = clin_s[lipsc].values.astype(float)
        y = clin_s[tmsc].values.astype(float)
        r, p, n = spearman_safe(x, y)
        rows_comp.append(dict(lipid_score=lipname, tmed_score=tmname, rho=r, p=p, n=n))

if rows_comp:
    df_comp = pd.DataFrame(rows_comp)
    _, fdr, _, _ = multipletests(df_comp["p"].fillna(1), method="fdr_bh")
    df_comp["FDR"] = fdr
    print(f"\nComposite score associations:")
    print(f"\n{'Lipid score':<16} {'TMED score':<14} {'rho':>8} {'p':>10} {'FDR':>10}")
    print("-"*60)
    for _, r in df_comp.iterrows():
        print(f"{r['lipid_score']:<16} {r['tmed_score']:<14} {r['rho']:>8.4f} {r['p']:>10.4f} {r['FDR']:>10.4f}")
    df_comp.to_csv(OUT + "h6_C_composite_scores.csv", index=False)

# ── TEST E: SM LIPIDOMICS ↔ TMED RNA ─────────────────────────
section("TEST E: SM LIPIDOMICS SPECIES ↔ TMED RNA")

print(f"  Testing {len(sm_cols)} SM species × {len(TMED_GENES)} TMED genes = {len(sm_cols)*len(TMED_GENES)} tests")
print(f"  Overlap sample n = {len(shared_lip)}")

rows_e = []
for tg in TMED_GENES:
    tc = f"rna_{tg}"
    if tc not in clin_sl.columns: continue
    tmed_v = clin_sl[tc].values.astype(float)
    for sm in sm_cols:
        sm_v = lip_s[sm].values.astype(float)
        r, p, n = spearman_safe(sm_v, tmed_v)
        rows_e.append(dict(sm_species=sm, tmed_gene=tg, rho=r, p=p, n=n))

df_e = pd.DataFrame(rows_e).dropna(subset=["p"])
_, fdr, _, _ = multipletests(df_e["p"], method="fdr_bh")
df_e["FDR"] = fdr
df_e = df_e.sort_values("p")

n_fdr = (df_e["FDR"] < 0.05).sum()
n_nom = (df_e["p"] < 0.05).sum()
print(f"\n  FDR < 0.05: {n_fdr}  |  nominal p < 0.05: {n_nom}  |  total tests: {len(df_e)}")
print(f"\n  Top 10 associations (by p-value):")
print(f"  {'SM species':<30} {'TMED':>8} {'rho':>8} {'p':>10} {'FDR':>10}")
print("  " + "-"*70)
for _, r in df_e.head(10).iterrows():
    print(f"  {r['sm_species']:<30} {r['tmed_gene']:>8} {r['rho']:>8.4f} {r['p']:>10.4f} {r['FDR']:>10.4f}")

df_e.to_csv(OUT + "h6_E_sm_lipidomics_tmed.csv", index=False)

# ── TEST F: SM-C18 SPECIES ────────────────────────────────────
section("TEST F: SM-C18 SPECIES ↔ TMED")

# look for C18 annotation in column names
c18_cols = [c for c in sm_cols if "18:0" in c or "C18" in c or "/18:" in c]
print(f"  SM-C18-annotated species found: {len(c18_cols)}")
for c in c18_cols:
    print(f"    {c}")

rows_f = []
if c18_cols:
    # use mean across C18 species if multiple
    lip_s["SM_C18_mean"] = lip_s[c18_cols].mean(axis=1)
    clin_sl = clin_sl.copy()
    clin_sl["SM_C18"] = lip_s["SM_C18_mean"].reindex(clin_sl.index)

    for tg in TMED_GENES:
        tc = f"rna_{tg}"
        if tc not in clin_sl.columns: continue
        x = clin_sl["SM_C18"].values.astype(float)
        y = clin_sl[tc].values.astype(float)
        r, p, n = spearman_safe(x, y)
        rows_f.append(dict(species="SM_C18_mean", tmed_gene=tg, rho=r, p=p, n=n))
        print(f"  SM-C18 ↔ {tg}: rho={r:.4f}, p={p:.4f}, n={n}")

    # APOE4 → SM-C18
    if apoe4_col and apoe4_col in clin_sl.columns:
        x = clin_sl[apoe4_col].values.astype(float)
        y = clin_sl["SM_C18"].values.astype(float)
        r, p, n = spearman_safe(x, y)
        print(f"\n  APOE4 ↔ SM-C18: rho={r:.4f}, p={p:.4f}, n={n}")

    pd.DataFrame(rows_f).to_csv(OUT + "h6_F_smc18_tmed.csv", index=False)
else:
    print("  No C18-specific SM columns found in annotation — Test F skipped")
    clin_sl["SM_C18"] = np.nan

# ── TEST G: NESTED MODELS ─────────────────────────────────────
section("TEST G: NESTED MODELS — APOE4 / SM / BOTH → TMED2")

# use SM-C18 if available, else top SM from Test E
if "SM_C18" in clin_sl.columns and clin_sl["SM_C18"].notna().sum() > 20:
    sm_pred_col = "SM_C18"
    sm_pred_name = "SM-C18"
else:
    # use top SM species for TMED2 from Test E
    top_e = df_e[df_e["tmed_gene"]=="TMED2"].sort_values("p").head(1)
    if len(top_e):
        best_sm = top_e.iloc[0]["sm_species"]
        clin_sl = clin_sl.copy()
        clin_sl["SM_best"] = lip_s[best_sm].reindex(clin_sl.index)
        sm_pred_col = "SM_best"
        sm_pred_name = best_sm
        print(f"  Using top SM for TMED2: {best_sm}")
    else:
        sm_pred_col = None
        sm_pred_name = None

rows_g = []
cov_g = [c for c in COV if c in clin_sl.columns]

for tg in TMED_GENES:
    tc = f"rna_{tg}"
    if tc not in clin_sl.columns: continue
    y_all = clin_sl[tc].values.astype(float)

    # base covariates only
    Xcov = np.column_stack([np.ones(len(clin_sl))] + [clin_sl[c].values.astype(float) for c in cov_g])
    res0, n0 = ols_fit(y_all, Xcov)
    r2_0 = np.nan
    if res0:
        yhat = Xcov @ res0["coeffs"]
        mask0 = ~(np.isnan(y_all)|np.isnan(Xcov).any(axis=1))
        ss_res = np.sum((y_all[mask0]-yhat[mask0])**2)
        ss_tot = np.sum((y_all[mask0]-y_all[mask0].mean())**2)
        r2_0 = 1 - ss_res/ss_tot if ss_tot>0 else np.nan

    # Model A: + APOE4
    if apoe4_col and apoe4_col in clin_sl.columns:
        ap_v = clin_sl[apoe4_col].values.astype(float)
        Xa = np.column_stack([Xcov, ap_v])
        resa, na = ols_fit(y_all, Xa)
        if resa:
            b_ap = resa["coeffs"][-1]; p_ap = resa["p"][-1]
            mask_a = ~(np.isnan(y_all)|np.isnan(Xa).any(axis=1))
            yhat_a = Xa[mask_a] @ resa["coeffs"]
            ss_res_a = np.sum((y_all[mask_a]-yhat_a)**2)
            ss_tot_a = np.sum((y_all[mask_a]-y_all[mask_a].mean())**2)
            r2_a = 1 - ss_res_a/ss_tot_a if ss_tot_a>0 else np.nan
            rows_g.append(dict(tmed=tg, model="MA: covariates+APOE4",
                               predictor="APOE4", beta=b_ap, p=p_ap, r2=r2_a, n=na))

    # Model B: + SM predictor
    if sm_pred_col and sm_pred_col in clin_sl.columns:
        sm_v = clin_sl[sm_pred_col].values.astype(float)
        Xb = np.column_stack([Xcov, sm_v])
        resb, nb = ols_fit(y_all, Xb)
        if resb:
            b_sm = resb["coeffs"][-1]; p_sm = resb["p"][-1]
            mask_b = ~(np.isnan(y_all)|np.isnan(Xb).any(axis=1))
            yhat_b = Xb[mask_b] @ resb["coeffs"]
            ss_res_b = np.sum((y_all[mask_b]-yhat_b)**2)
            ss_tot_b = np.sum((y_all[mask_b]-y_all[mask_b].mean())**2)
            r2_b = 1 - ss_res_b/ss_tot_b if ss_tot_b>0 else np.nan
            rows_g.append(dict(tmed=tg, model=f"MB: covariates+{sm_pred_name}",
                               predictor=sm_pred_name, beta=b_sm, p=p_sm, r2=r2_b, n=nb))

    # Model C: + APOE4 + SM
    if (apoe4_col and apoe4_col in clin_sl.columns
            and sm_pred_col and sm_pred_col in clin_sl.columns):
        ap_v  = clin_sl[apoe4_col].values.astype(float)
        sm_v  = clin_sl[sm_pred_col].values.astype(float)
        Xc = np.column_stack([Xcov, ap_v, sm_v])
        resc, nc = ols_fit(y_all, Xc)
        if resc:
            b_ap2 = resc["coeffs"][-2]; p_ap2 = resc["p"][-2]
            b_sm2 = resc["coeffs"][-1]; p_sm2 = resc["p"][-1]
            mask_c = ~(np.isnan(y_all)|np.isnan(Xc).any(axis=1))
            yhat_c = Xc[mask_c] @ resc["coeffs"]
            ss_res_c = np.sum((y_all[mask_c]-yhat_c)**2)
            ss_tot_c = np.sum((y_all[mask_c]-y_all[mask_c].mean())**2)
            r2_c = 1-ss_res_c/ss_tot_c if ss_tot_c>0 else np.nan
            rows_g.append(dict(tmed=tg, model=f"MC: covariates+APOE4+{sm_pred_name}",
                               predictor="APOE4", beta=b_ap2, p=p_ap2, r2=r2_c, n=nc))
            rows_g.append(dict(tmed=tg, model=f"MC: covariates+APOE4+{sm_pred_name}",
                               predictor=sm_pred_name, beta=b_sm2, p=p_sm2, r2=r2_c, n=nc))

df_g = pd.DataFrame(rows_g)
if len(df_g):
    print(f"\n{'TMED':<8} {'Model':<40} {'Predictor':<16} {'β':>10} {'p':>10} {'R²':>8} {'n':>6}")
    print("-"*100)
    for _, r in df_g.iterrows():
        print(f"{r['tmed']:<8} {r['model']:<40} {r['predictor']:<16} "
              f"{r['beta']:>10.4f} {r['p']:>10.4f} {r['r2']:>8.4f} {int(r['n']):>6}")
    df_g.to_csv(OUT + "h6_G_nested_models.csv", index=False)

print("\n" + "="*60)
print("=== DONE ===")
print("="*60)
print(f"\nOutputs written to: {OUT}")
for f in sorted(os.listdir(OUT)):
    print(f"  {f}")
