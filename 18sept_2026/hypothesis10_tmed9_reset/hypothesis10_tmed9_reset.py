"""
hypothesis10_tmed9_reset.py
H4 / Test 4: TMED9 → ER stress / RESET
CHRONOS / ROSMAP  —  Sep 2026

Tests
  A  TMED9 ↔ individual ER-stress markers (Spearman)
  B  TMED9 ↔ ER-stress composite score
  C  TMED9 ↔ disease severity (Braak, CERAD, cogdx)
  D  Composition-adjusted pathology models
  E  RESET vs compensation (nested models: TMED9 ~ TMED2 + ER_score + Braak)
  F  snRNA pseudobulk: TMED9 ~ pathology within each cell type
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────────
DATA = "/home/wordson22/projects/CHRONOS/DATA/"
OUT  = "/home/wordson22/projects/CHRONOS/analysis/hypothesis10_tmed9_reset/"

import os; os.makedirs(OUT, exist_ok=True)

# ── gene IDs ───────────────────────────────────────────────────────────────
TMED9_ID  = "ENSG00000184840"
TMED2_ID  = "ENSG00000086598"
TMED10_ID = "ENSG00000170348"

ER_MARKER_IDS = {
    "EIF2AK3" : "ENSG00000172071",   # PERK
    "ATF4"    : "ENSG00000128272",
    "DDIT3"   : "ENSG00000175197",   # CHOP
    "HSPA5"   : "ENSG00000044574",   # GRP78/BiP
    "PPP1R15A": "ENSG00000087074",   # GADD34
}

CELL_MARKERS = {
    "neuron_score": ["ENSG00000102003","ENSG00000067715","ENSG00000132639","ENSG00000008056","ENSG00000157542"],
    "astro_score" : ["ENSG00000131095","ENSG00000171885"],
    "micro_score" : ["ENSG00000204472","ENSG00000138185"],
    "oligo_score" : ["ENSG00000197971","ENSG00000123560"],
}
CELL_COLS = list(CELL_MARKERS.keys())
COV_BASE  = ["age_at_visit_max","msex","pmi"]

# ── helpers ────────────────────────────────────────────────────────────────

def ols_term(y_s, X_df, term):
    cols = list(X_df.columns)
    Xm   = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    yv   = y_s.values if hasattr(y_s,"values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(yv))
    n_ok = int(mask.sum())
    if n_ok < len(cols)+3:
        return dict(beta=np.nan,se=np.nan,p=np.nan,n=n_ok)
    Xm,yv = Xm[mask],yv[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm,yv,rcond=None)
    resid = yv - Xm@coeffs
    n,p   = len(yv),Xm.shape[1]
    mse   = np.sum(resid**2)/(n-p)
    try:
        cov = mse*np.linalg.inv(Xm.T@Xm)
    except np.linalg.LinAlgError:
        return dict(beta=np.nan,se=np.nan,p=np.nan,n=n_ok)
    se_all = np.sqrt(np.diag(cov))
    idx    = (["intercept"]+cols).index(term)
    beta,se_b = coeffs[idx],se_all[idx]
    t  = beta/se_b
    pv = 2*stats.t.sf(abs(t),df=n-p)
    return dict(beta=beta,se=se_b,p=pv,n=n_ok)


def r2_model(y_s, X_df):
    cols = list(X_df.columns)
    Xm   = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    yv   = y_s.values if hasattr(y_s,"values") else np.array(y_s)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(yv))
    if mask.sum() < len(cols)+3:
        return np.nan, np.nan
    Xm,yv = Xm[mask],yv[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm,yv,rcond=None)
    yhat = Xm@coeffs
    ss_res = np.sum((yv-yhat)**2)
    ss_tot = np.sum((yv-yv.mean())**2)
    n,p = len(yv),Xm.shape[1]
    r2  = 1 - ss_res/ss_tot
    aic = n*np.log(ss_res/n) + 2*p
    return r2, aic


def spearman_sub(x,y):
    mask = ~(np.isnan(x)|np.isnan(y))
    if mask.sum()<5:
        return np.nan,np.nan,int(mask.sum())
    rho,p = stats.spearmanr(x[mask],y[mask])
    return rho,p,int(mask.sum())


def bh_fdr(pvals):
    out  = np.full(len(pvals),np.nan)
    mask = ~np.isnan(pvals)
    if mask.sum()==0: return out
    _,adj,_,_ = multipletests(np.array(pvals)[mask],method="fdr_bh")
    out[mask]  = adj
    return out


def section(title):
    print("\n"+"="*60)
    print(f"=== {title} ===")
    print("="*60)


# ══════════════════════════════════════════════════════════════
# LOAD BULK RNA DATA
# ══════════════════════════════════════════════════════════════
print("Loading bulk RNA data...")
expr     = pd.read_csv(DATA+"ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA+"ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA+"ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA+"ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = pd.to_numeric(
    clinical["age_at_visit_max"].replace("90+","90"), errors="coerce")
for col in ["braaksc","msex","pmi","cogdx","ceradsc"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")
clinical["AD"] = np.nan
clinical.loc[clinical["cogdx"]==1,"AD"] = 0
clinical.loc[clinical["cogdx"]>=4,"AD"] = 1

# RNA linking
rna_bio  = biospec[["individualID","specimenID"]].drop_duplicates()
rna_link = rna_meta[["specimenID","rnaBatch"]].merge(rna_bio, on="specimenID", how="left")
rna_link = rna_link[rna_link["specimenID"].isin(expr.columns)]
rna_link = rna_link.merge(
    clinical[["individualID","braaksc","msex","pmi","age_at_visit_max","cogdx","AD","ceradsc"]],
    on="individualID", how="left").dropna(subset=["individualID"])

# rnaBatch
_rb = rna_link["rnaBatch"].dropna().unique()
if len(_rb)>1:
    rna_link["rnaBatch_cat"] = pd.Categorical(rna_link["rnaBatch"]).codes.astype(float)
    rna_link.loc[rna_link["rnaBatch"].isna(),"rnaBatch_cat"] = np.nan
    _USE_BATCH = True
else:
    rna_link["rnaBatch_cat"] = 0.0
    _USE_BATCH = False
BATCH_COLS = ["rnaBatch_cat"] if _USE_BATCH else []

# cell scores
def cell_score(markers):
    avail = [m for m in markers if m in expr.index]
    if not avail: return pd.Series(np.nan, index=expr.columns)
    z = expr.loc[avail].apply(lambda r:(r-r.mean())/r.std(), axis=1)
    return z.mean(axis=0)

for cname,markers in CELL_MARKERS.items():
    rna_link[cname] = rna_link["specimenID"].map(cell_score(markers))

# TMED9, TMED2, ER markers
ALL_GENE_IDS = {"TMED9":TMED9_ID,"TMED2":TMED2_ID,"TMED10":TMED10_ID,**ER_MARKER_IDS}
for name,eid in ALL_GENE_IDS.items():
    if eid in expr.index:
        rna_link[f"rna_{name}"] = rna_link["specimenID"].map(expr.loc[eid])
    else:
        print(f"  WARNING: {name} ({eid}) not in expression matrix")
        rna_link[f"rna_{name}"] = np.nan

# ER composite score (z-score each marker, average)
er_cols = [f"rna_{m}" for m in ER_MARKER_IDS]
er_avail = [c for c in er_cols if not rna_link[c].isna().all()]
if er_avail:
    z_er = rna_link[er_avail].apply(lambda col:(col-col.mean())/col.std(), axis=0)
    rna_link["ER_score"] = z_er.mean(axis=1)
else:
    rna_link["ER_score"] = np.nan

COV_CELLS = COV_BASE + BATCH_COLS + CELL_COLS
rna_sub = rna_link.dropna(subset=COV_CELLS+["rna_TMED9"]).copy()
print(f"  Bulk RNA n (complete covariates): {len(rna_sub)}")


# ══════════════════════════════════════════════════════════════
# TEST A: TMED9 ↔ ER-stress markers (Spearman)
# ══════════════════════════════════════════════════════════════
section("TEST A: TMED9 ↔ ER-STRESS MARKERS")

rows_a = []
for marker in ER_MARKER_IDS:
    x = rna_link["rna_TMED9"].values.astype(float)
    y = rna_link[f"rna_{marker}"].values.astype(float)
    rho,p,n = spearman_sub(x,y)
    rows_a.append(dict(Marker=marker, rho=rho, p=p, n=n))

df_a = pd.DataFrame(rows_a)
df_a["FDR"] = bh_fdr(df_a["p"].values)

print(f"\n{'Marker':<12} {'rho':>8} {'p':>9} {'FDR':>9} {'n':>5}")
print("-"*47)
for _,r in df_a.iterrows():
    print(f"{r['Marker']:<12} {r['rho']:>8.4f} {r['p']:>9.4f} {r['FDR']:>9.4f} {int(r['n']):>5}")

df_a.to_csv(OUT+"h4_A_er_markers.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST B: TMED9 ↔ ER-stress composite
# ══════════════════════════════════════════════════════════════
section("TEST B: TMED9 ↔ ER-STRESS COMPOSITE")

x = rna_link["rna_TMED9"].values.astype(float)
y = rna_link["ER_score"].values.astype(float)
rho_b, p_b, n_b = spearman_sub(x, y)
print(f"\n  TMED9 ↔ ER-stress composite:")
print(f"  rho = {rho_b:.4f}  p = {p_b:.4f}  n = {n_b}")

# also composition-adjusted
sub_b = rna_link.dropna(subset=COV_CELLS+["rna_TMED9","ER_score"])
X_b   = sub_b[["ER_score"]+COV_CELLS]
res_b = ols_term(sub_b["rna_TMED9"], X_b, "ER_score")
print(f"\n  Adjusted (ER_score → TMED9):")
print(f"  β = {res_b['beta']:.4f}  p = {res_b['p']:.4f}  n = {res_b['n']}")

pd.DataFrame([dict(rho=rho_b, p_spearman=p_b, n_sp=n_b,
                   beta_adj=res_b["beta"], p_adj=res_b["p"], n_adj=res_b["n"])]
             ).to_csv(OUT+"h4_B_er_composite.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST C: TMED9 ↔ disease severity (Spearman)
# ══════════════════════════════════════════════════════════════
section("TEST C: TMED9 ↔ DISEASE SEVERITY")

PATH_OUTCOMES = {"braaksc":"Braak","ceradsc":"CERAD","cogdx":"cogdx"}
rows_c = []
for col,label in PATH_OUTCOMES.items():
    x = rna_link["rna_TMED9"].values.astype(float)
    y = rna_link[col].values.astype(float)
    rho,p,n = spearman_sub(x,y)
    rows_c.append(dict(Outcome=label, rho=rho, p=p, n=n))

df_c = pd.DataFrame(rows_c)
df_c["FDR"] = bh_fdr(df_c["p"].values)

print(f"\n{'Outcome':<10} {'rho':>8} {'p':>9} {'FDR':>9} {'n':>5}")
print("-"*45)
for _,r in df_c.iterrows():
    print(f"{r['Outcome']:<10} {r['rho']:>8.4f} {r['p']:>9.4f} {r['FDR']:>9.4f} {int(r['n']):>5}")

df_c.to_csv(OUT+"h4_C_disease_severity.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST D: Composition-adjusted pathology models
# ══════════════════════════════════════════════════════════════
section("TEST D: COMPOSITION-ADJUSTED PATHOLOGY MODELS")

rows_d = []
for path_col, path_label in PATH_OUTCOMES.items():
    sub = rna_link.dropna(subset=["rna_TMED9", path_col])

    # Model 1: simple
    X1  = sub[[path_col]]
    r1  = ols_term(sub["rna_TMED9"], X1, path_col)
    rows_d.append(dict(Outcome=path_label, Model="simple",
                       beta=r1["beta"], p=r1["p"], n=r1["n"]))

    # Model 2: + demographics + batch
    sub2 = sub.dropna(subset=COV_BASE+BATCH_COLS)
    X2   = sub2[[path_col]+COV_BASE+BATCH_COLS]
    r2   = ols_term(sub2["rna_TMED9"], X2, path_col)
    rows_d.append(dict(Outcome=path_label, Model="+demo+batch",
                       beta=r2["beta"], p=r2["p"], n=r2["n"]))

    # Model 3: + cells
    sub3 = sub2.dropna(subset=CELL_COLS)
    X3   = sub3[[path_col]+COV_BASE+BATCH_COLS+CELL_COLS]
    r3   = ols_term(sub3["rna_TMED9"], X3, path_col)
    rows_d.append(dict(Outcome=path_label, Model="+demo+batch+cells",
                       beta=r3["beta"], p=r3["p"], n=r3["n"]))

df_d = pd.DataFrame(rows_d)
print(f"\n{'Outcome':<10} {'Model':<22} {'β':>9} {'p':>9} {'n':>5}")
print("-"*57)
for _,r in df_d.iterrows():
    print(f"{r['Outcome']:<10} {r['Model']:<22} {r['beta']:>9.4f} {r['p']:>9.4f} {int(r['n']):>5}")

df_d.to_csv(OUT+"h4_D_composition_adjusted.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST E: RESET vs compensation (nested models)
# ══════════════════════════════════════════════════════════════
section("TEST E: RESET vs COMPENSATION (NESTED MODELS)")

sub_e = rna_link.dropna(subset=["rna_TMED9","rna_TMED2","ER_score","braaksc"]+COV_CELLS).copy()
print(f"\n  n for nested models: {len(sub_e)}")

models = [
    ("M1: TMED2",                          ["rna_TMED2"]),
    ("M2: TMED2 + ER",                     ["rna_TMED2","ER_score"]),
    ("M3: TMED2 + ER + Braak",             ["rna_TMED2","ER_score","braaksc"]),
    ("M4: TMED2 + ER + Braak + cov+cells", ["rna_TMED2","ER_score","braaksc"]+COV_CELLS),
]

rows_e = []
print(f"\n{'Model':<38} {'term':<12} {'β':>9} {'p':>9} {'R²':>7} {'n':>5}")
print("-"*80)
for label, predictors in models:
    X_df = sub_e[predictors]
    r2,_ = r2_model(sub_e["rna_TMED9"], X_df)
    for term in predictors:
        res = ols_term(sub_e["rna_TMED9"], X_df, term)
        rows_e.append(dict(Model=label, Term=term,
                           beta=res["beta"], p=res["p"], r2=r2, n=res["n"]))
        print(f"{label:<38} {term:<12} {res['beta']:>9.4f} {res['p']:>9.4f} "
              f"{r2:>7.4f} {int(res['n']):>5}")

df_e = pd.DataFrame(rows_e)
df_e.to_csv(OUT+"h4_E_nested_models.csv", index=False)

# TMED9 ~ TMED2 alone: key compensation test
res_comp = ols_term(rna_link.dropna(subset=["rna_TMED9","rna_TMED2"])["rna_TMED9"],
                    rna_link.dropna(subset=["rna_TMED9","rna_TMED2"])[["rna_TMED2"]],
                    "rna_TMED2")
print(f"\n  TMED9 ~ TMED2 (simple): β={res_comp['beta']:.4f} p={res_comp['p']:.4f} "
      f"n={res_comp['n']}  (negative β = compensation)")


# ══════════════════════════════════════════════════════════════
# TEST F: snRNA pseudobulk — TMED9 within each cell type
# ══════════════════════════════════════════════════════════════
section("TEST F: snRNA PSEUDOBULK — TMED9 WITHIN CELL TYPES")

print("\nLoading snRNA data (this may take a minute)...")

try:
    # Load sparse matrix
    mat   = mmread(DATA+"filtered_count_matrix.mtx").tocsc()
    genes = pd.read_csv(DATA+"filtered_gene_row_names.txt", header=None, names=["gene"])
    meta  = pd.read_csv(DATA+"filtered_column_metadata.txt", sep="\t")
    id_map= pd.read_csv(DATA+"snRNAseqPFC_BA10_id_mapping.csv")

    print(f"  Matrix: {mat.shape[0]} genes × {mat.shape[1]} cells")
    print(f"  Metadata columns: {list(meta.columns)}")
    print(f"  ID map columns: {list(id_map.columns)}")

    # Find TMED9 row
    tmed9_rows = genes[genes["gene"].str.contains("TMED9", case=False, na=False)].index.tolist()
    if not tmed9_rows:
        # try ENSEMBL ID
        tmed9_rows = genes[genes["gene"]==TMED9_ID].index.tolist()
    if not tmed9_rows:
        print("  WARNING: TMED9 not found in snRNA gene list — checking first 20 gene names:")
        print(f"  {genes['gene'].head(20).tolist()}")
        raise ValueError("TMED9 not found")

    tmed9_idx = tmed9_rows[0]
    print(f"  TMED9 found at row {tmed9_idx}: {genes.loc[tmed9_idx,'gene']}")

    # Extract TMED9 counts per cell
    tmed9_vec = np.array(mat[tmed9_idx, :].todense()).flatten()

    # Detect cell-type column
    ct_col = None
    for cand in ["broad.cell.type","cellType","cell_type","Celltype","cluster"]:
        if cand in meta.columns:
            ct_col = cand; break
    if ct_col is None:
        print(f"  Could not find cell-type column. Available: {list(meta.columns)}")
        raise ValueError("No cell-type column")
    print(f"  Cell-type column: '{ct_col}'")
    print(f"  Cell types: {sorted(meta[ct_col].unique())}")

    # Donor column
    donor_col = None
    for cand in ["projid","Subject","donor","individualID","donorID","TAG"]:
        if cand in meta.columns:
            donor_col = cand; break
    if donor_col is None:
        # use id_map to find
        print(f"  Could not find donor column. Available: {list(meta.columns)}")
        raise ValueError("No donor column")
    print(f"  Donor column: '{donor_col}'")

    meta = meta.copy()
    meta["TMED9_count"] = tmed9_vec

    # Map donor → projid if needed (id_map links TAG/barcode to projid)
    if donor_col != "projid" and "projid" not in meta.columns:
        # id_map likely has the bridge
        id_map_cols = list(id_map.columns)
        print(f"  id_map columns: {id_map_cols}")
        shared = [c for c in id_map_cols if c in meta.columns]
        if shared:
            meta = meta.merge(id_map[id_map_cols], on=shared[0], how="left")
            if "projid" in meta.columns:
                donor_col = "projid"

    # Merge clinical
    clin_sn = clinical[["projid","braaksc","ceradsc","cogdx","AD"]].copy()
    clin_sn["_projid_key"] = clin_sn["projid"].astype(str).str.strip()
    meta["_projid_key"] = meta[donor_col].astype(str).str.strip()
    meta = meta.merge(clin_sn.drop(columns=["projid"]), on="_projid_key", how="left")
    print(f"  After clinical merge: {(~meta['braaksc'].isna()).sum()} cells with Braak data")

    # Library-size normalize per cell (log1p CPM)
    col_sums = np.array(mat.sum(axis=0)).flatten()
    col_sums[col_sums==0] = 1
    tmed9_norm = np.log1p(tmed9_vec / col_sums * 1e6)
    meta["TMED9_logCPM"] = tmed9_norm

    # Pseudobulk: mean logCPM per donor per cell type
    pb = meta.groupby([donor_col, ct_col]).agg(
        TMED9_mean=("TMED9_logCPM","mean"),
        n_cells=("TMED9_logCPM","count"),
        braaksc=("braaksc","first"),
        ceradsc=("ceradsc","first"),
        cogdx=("cogdx","first"),
        AD=("AD","first")
    ).reset_index()

    # Only keep donor×celltype with ≥5 cells
    pb = pb[pb["n_cells"] >= 5]
    print(f"\n  Pseudobulk entries (≥5 cells): {len(pb)}")

    cell_types = sorted(pb[ct_col].unique())
    print(f"  Cell types with data: {cell_types}")

    rows_f = []
    print(f"\n{'Cell type':<25} {'Braak ρ':>9} {'p':>9} {'CERAD ρ':>9} {'p':>9} {'n donors':>9}")
    print("-"*72)
    for ct in cell_types:
        sub_ct = pb[pb[ct_col]==ct].dropna(subset=["TMED9_mean"])
        n_don  = sub_ct[donor_col].nunique()
        # Braak
        rho_b,p_b,_ = spearman_sub(sub_ct["TMED9_mean"].values.astype(float),
                                    sub_ct["braaksc"].values.astype(float))
        # CERAD
        rho_c,p_c,_ = spearman_sub(sub_ct["TMED9_mean"].values.astype(float),
                                    sub_ct["ceradsc"].values.astype(float))
        rows_f.append(dict(cell_type=ct, rho_braak=rho_b, p_braak=p_b,
                           rho_cerad=rho_c, p_cerad=p_c, n_donors=n_don))
        print(f"{ct:<25} {rho_b:>9.4f} {p_b:>9.4f} {rho_c:>9.4f} {p_c:>9.4f} {n_don:>9}")

    df_f = pd.DataFrame(rows_f)
    df_f["FDR_braak"] = bh_fdr(df_f["p_braak"].values)
    df_f["FDR_cerad"] = bh_fdr(df_f["p_cerad"].values)
    df_f.to_csv(OUT+"h4_F_snrna_pseudobulk.csv", index=False)

except Exception as e:
    print(f"\n  snRNA test failed: {e}")
    print("  Continuing without Test F.")
    df_f = pd.DataFrame()


# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
section("DONE")
print(f"\nOutputs written to: {OUT}")
print("  h4_A_er_markers.csv")
print("  h4_B_er_composite.csv")
print("  h4_C_disease_severity.csv")
print("  h4_D_composition_adjusted.csv")
print("  h4_E_nested_models.csv")
print("  h4_F_snrna_pseudobulk.csv  (if snRNA loaded)")
