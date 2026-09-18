"""
Hypothesis 8: Three-protein composite biomarker score
Tests A-K: ratio composite, PCA, AUC, DeLong, Braak R², CV, permutation
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

# ── paths ──────────────────────────────────────────────────────────────────────
DATA = "/home/wordson22/projects/CHRONOS/DATA/"
OUT  = "./"

# ── Ensembl IDs ────────────────────────────────────────────────────────────────
TMED2_ID  = "ENSG00000086598"
TMED10_ID = "ENSG00000170348"
TMED9_ID  = "ENSG00000184840"

# ── load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = pd.to_numeric(
    clinical["age_at_visit_max"].replace("90+", "90"), errors="coerce")
for col in ["braaksc", "msex", "pmi", "cogdx"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

# ── link samples ───────────────────────────────────────────────────────────────
clin_cols = ["individualID", "braaksc", "cogdx", "msex", "pmi", "age_at_visit_max"]
sample_ids = expr.columns.tolist()
rna_bio    = biospec[biospec["specimenID"].isin(sample_ids)][
    ["specimenID", "individualID"]].drop_duplicates()
rna_batch  = rna_meta[["specimenID", "rnaBatch"]].drop_duplicates()
meta = (rna_bio
        .merge(clinical[clin_cols], on="individualID", how="left")
        .merge(rna_batch, on="specimenID", how="left")
        .set_index("specimenID"))

# ── cell type scores ───────────────────────────────────────────────────────────
NEURON_M = ["ENSG00000102003","ENSG00000067715","ENSG00000132639",
            "ENSG00000008056","ENSG00000157542"]
ASTRO_M  = ["ENSG00000131095","ENSG00000171885"]
MICRO_M  = ["ENSG00000204472","ENSG00000138185"]
OLIGO_M  = ["ENSG00000197971","ENSG00000123560"]

def ct_score(markers, name):
    present = [m for m in markers if m in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns, name=name)
    vals = expr.loc[present]
    z = (vals - vals.mean(axis=1).values[:,None]) / (vals.std(axis=1).values[:,None] + 1e-10)
    s = z.mean(axis=0); s.name = name; return s

meta["neuron_score"]  = ct_score(NEURON_M, "neuron_score")
meta["astro_score"]   = ct_score(ASTRO_M,  "astro_score")
meta["micro_score"]   = ct_score(MICRO_M,  "micro_score")
meta["oligo_score"]   = ct_score(OLIGO_M,  "oligo_score")

# ── batch dummies ──────────────────────────────────────────────────────────────
batch_dum  = pd.get_dummies(meta["rnaBatch"], prefix="batch", drop_first=True)
meta       = pd.concat([meta, batch_dum], axis=1)
batch_cols = list(batch_dum.columns)

# ── extract TMED genes ─────────────────────────────────────────────────────────
for eid, name in [(TMED2_ID,"TMED2"),(TMED10_ID,"TMED10"),(TMED9_ID,"TMED9")]:
    if eid in expr.index:
        meta[name] = expr.loc[eid]
    else:
        print(f"  WARNING: {name} ({eid}) not found"); meta[name] = np.nan

# ── working dataframe ──────────────────────────────────────────────────────────
needed = ["braaksc","cogdx","msex","pmi","age_at_visit_max",
          "neuron_score","astro_score","micro_score","oligo_score",
          "TMED2","TMED10","TMED9"] + batch_cols
df = meta[needed].dropna(subset=["braaksc","TMED2","TMED10","TMED9"])
print(f"Working n = {len(df)}")

# ── OLS helper ─────────────────────────────────────────────────────────────────
def ols_term(y_s, X_df, term):
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm @ coeffs
    n, p  = len(y_v), Xm.shape[1]
    mse   = np.sum(resid**2) / (n - p)
    cov   = mse * np.linalg.inv(Xm.T @ Xm)
    se    = np.sqrt(np.diag(cov))
    cols  = ["intercept"] + list(X_df.columns)
    idx   = cols.index(term)
    beta, se_b = coeffs[idx], se[idx]
    t  = beta / se_b
    pv = 2 * stats.t.sf(abs(t), df=n - p)
    # R²
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y_v - y_v.mean())**2)
    r2 = 1 - ss_res/ss_tot
    return dict(beta=beta, se=se_b, ci_lo=beta-1.96*se_b, ci_hi=beta+1.96*se_b,
                p=pv, n=int(mask.sum()), r2=r2)

def r2_model(y_s, X_df):
    """Returns R² for full model (no specific term)."""
    Xm = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    y_v = y_s.values
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(y_v))
    Xm, y_v = Xm[mask], y_v[mask]
    coeffs,_,_,_ = np.linalg.lstsq(Xm, y_v, rcond=None)
    resid = y_v - Xm @ coeffs
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y_v - y_v.mean())**2)
    return 1 - ss_res/ss_tot, int(mask.sum())

# ── AD binary ─────────────────────────────────────────────────────────────────
df2 = df.dropna(subset=["cogdx"])
ad_mask = df2["cogdx"].isin([1, 4, 5, 6])
df_ad   = df2[ad_mask].copy()
df_ad["AD"] = (df_ad["cogdx"] >= 4).astype(int)
print(f"AD classification n={len(df_ad)}  NCI={( df_ad['AD']==0).sum()}  AD={( df_ad['AD']==1).sum()}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST A: Build H8_ratio = TMED9 - TMED2 - TMED10  (log-scale ratio)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST A: H8 ratio composite ──")
df["H8_ratio"] = df["TMED9"] - df["TMED2"] - df["TMED10"]
df_ad["H8_ratio"] = df_ad["TMED9"] - df_ad["TMED2"] - df_ad["TMED10"]

r_tmed2,  p_tmed2  = stats.pearsonr(df["H8_ratio"].dropna(), df.loc[df["H8_ratio"].notna(),"TMED2"])
r_tmed10, p_tmed10 = stats.pearsonr(df["H8_ratio"].dropna(), df.loc[df["H8_ratio"].notna(),"TMED10"])
r_tmed9,  p_tmed9  = stats.pearsonr(df["H8_ratio"].dropna(), df.loc[df["H8_ratio"].notna(),"TMED9"])

print(f"  H8_ratio mean={df['H8_ratio'].mean():.3f}  SD={df['H8_ratio'].std():.3f}")
print(f"  corr(ratio, TMED2) = {r_tmed2:.3f}  p={p_tmed2:.4f}")
print(f"  corr(ratio, TMED10)= {r_tmed10:.3f}  p={p_tmed10:.4f}")
print(f"  corr(ratio, TMED9) = {r_tmed9:.3f}  p={p_tmed9:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST B: PCA of TMED2, TMED10, TMED9
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST B: PCA ──")
pca_df = df[["TMED2","TMED10","TMED9"]].dropna()
scaler = StandardScaler()
X_scaled = scaler.fit_transform(pca_df)
pca = PCA(n_components=3)
pca.fit(X_scaled)
pca_scores = pca.transform(X_scaled)
df.loc[pca_df.index, "PC1"] = pca_scores[:,0]
df.loc[pca_df.index, "PC2"] = pca_scores[:,1]
df.loc[pca_df.index, "PC3"] = pca_scores[:,2]
df_ad["PC1"] = df.loc[df_ad.index, "PC1"] if "PC1" in df.columns else np.nan
df_ad.loc[df_ad.index, "PC1"] = df.loc[df_ad.index.intersection(df.index), "PC1"]

print(f"  Explained variance: PC1={pca.explained_variance_ratio_[0]:.3f}  "
      f"PC2={pca.explained_variance_ratio_[1]:.3f}  PC3={pca.explained_variance_ratio_[2]:.3f}")
print(f"  PC1 loadings: TMED2={pca.components_[0,0]:.3f}  "
      f"TMED10={pca.components_[0,1]:.3f}  TMED9={pca.components_[0,2]:.3f}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST C: AUC (AD vs NCI) for each feature
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST C: AUC comparison ──")
features = ["TMED2","TMED10","TMED9","H8_ratio","PC1"]
auc_results = {}

for feat in features:
    sub = df_ad[["AD", feat]].dropna()
    if len(sub) < 20 or sub["AD"].nunique() < 2:
        auc_results[feat] = dict(auc=np.nan, n=len(sub))
        continue
    y_true = sub["AD"].values
    y_score = sub[feat].values
    auc = roc_auc_score(y_true, y_score)
    # flip if AUC < 0.5 (direction)
    if auc < 0.5:
        auc = 1 - auc
        y_score = -y_score
    auc_results[feat] = dict(auc=auc, n=len(sub), y_true=y_true, y_score=y_score)
    print(f"  {feat:12s}  AUC={auc:.3f}  n={len(sub)}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST D: DeLong AUC comparison (composite vs each individual gene)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST D: DeLong test ──")

def delong_auc_variance(y_true, y_score):
    """Returns AUC and its variance via DeLong method."""
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    m, n = len(pos), len(neg)
    # placement values
    vpos = np.array([np.mean(p > neg) + 0.5*np.mean(p == neg) for p in pos])
    vneg = np.array([np.mean(n_ < pos) + 0.5*np.mean(n_ == pos) for n_ in neg])
    auc  = vpos.mean()
    s10  = np.var(vpos, ddof=1) / m
    s01  = np.var(vneg, ddof=1) / n
    var  = s10 + s01
    return auc, var

def delong_compare(y_true, y_score_a, y_score_b):
    """Two-sided DeLong test: H0: AUC_a == AUC_b"""
    pos = y_true == 1
    neg = y_true == 0
    m, n = pos.sum(), neg.sum()
    # placement for a
    def placements(score, ytrue):
        p = score[ytrue==1]; g = score[ytrue==0]
        vp = np.array([np.mean(pi > g) + 0.5*np.mean(pi==g) for pi in p])
        vn = np.array([np.mean(gi < p) + 0.5*np.mean(gi==p) for gi in g])
        return vp, vn
    vpa, vna = placements(y_score_a, y_true)
    vpb, vnb = placements(y_score_b, y_true)
    auc_a = vpa.mean(); auc_b = vpb.mean()
    # covariance matrix
    s_aa = (np.cov(vpa,vpa)[0,1]/m + np.cov(vna,vna)[0,1]/n) if m>1 and n>1 else 0
    s_bb = (np.cov(vpb,vpb)[0,1]/m + np.cov(vnb,vnb)[0,1]/n) if m>1 and n>1 else 0
    s_ab = (np.cov(vpa,vpb)[0,1]/m + np.cov(vna,vnb)[0,1]/n) if m>1 and n>1 else 0
    var_diff = s_aa + s_bb - 2*s_ab
    if var_diff <= 0:
        return auc_a, auc_b, np.nan, np.nan
    z = (auc_a - auc_b) / np.sqrt(var_diff)
    p = 2 * stats.norm.sf(abs(z))
    return auc_a, auc_b, z, p

delong_rows = []
ref_feat = "H8_ratio"
ref = auc_results.get(ref_feat, {})
if "y_true" in ref:
    for comp in ["TMED2","TMED10","TMED9"]:
        c = auc_results.get(comp, {})
        if "y_true" in c:
            sub_idx = pd.Series(df_ad.index).isin(
                df_ad[[ref_feat, comp, "AD"]].dropna().index)
            sub = df_ad[[ref_feat, comp, "AD"]].dropna()
            auc_r, auc_c, z, p = delong_compare(
                sub["AD"].values, sub[ref_feat].values, sub[comp].values)
            delong_rows.append(dict(reference=ref_feat, comparator=comp,
                                    auc_ref=auc_r, auc_comp=auc_c, z=z, p=p))
            print(f"  {ref_feat} vs {comp}: AUC_ref={auc_r:.3f}  AUC_comp={auc_c:.3f}  z={z:.3f}  p={p:.4f}")

delong_df = pd.DataFrame(delong_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST E/F: Braak R² for individual genes vs composite
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST E/F: Braak R² ──")
cov_cols = ["msex","pmi","age_at_visit_max","neuron_score","astro_score",
            "micro_score","oligo_score"] + batch_cols
braak_rows = []

for feat in ["TMED2","TMED10","TMED9","H8_ratio","PC1"]:
    sub = df[[feat,"braaksc"] + cov_cols].dropna()
    if len(sub) < 30:
        braak_rows.append(dict(feature=feat, r2=np.nan, n=len(sub))); continue
    # null model R² (covariates only)
    X_null = sub[cov_cols]
    r2_null, _ = r2_model(sub["braaksc"], X_null)
    # full model
    X_full = sub[[feat] + cov_cols]
    r2_full, n = r2_model(sub["braaksc"], X_full)
    delta_r2 = r2_full - r2_null
    braak_rows.append(dict(feature=feat, r2_null=r2_null, r2_full=r2_full,
                           delta_r2=delta_r2, n=n))
    print(f"  {feat:12s}  R²_null={r2_null:.4f}  R²_full={r2_full:.4f}  ΔR²={delta_r2:.4f}  n={n}")

braak_df = pd.DataFrame(braak_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST G: Three-gene multivariable model
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST G: Three-gene multivariable Braak model ──")
sub = df[["TMED2","TMED10","TMED9","braaksc"] + cov_cols].dropna()
if len(sub) > 30:
    for gene in ["TMED2","TMED10","TMED9"]:
        res = ols_term(sub["braaksc"], sub[[gene] + cov_cols], gene)
        print(f"  {gene}: β={res['beta']:.4f}  p={res['p']:.4f}  n={res['n']}")
    # all three together
    res3_tmed2  = ols_term(sub["braaksc"], sub[["TMED2","TMED10","TMED9"] + cov_cols], "TMED2")
    res3_tmed10 = ols_term(sub["braaksc"], sub[["TMED2","TMED10","TMED9"] + cov_cols], "TMED10")
    res3_tmed9  = ols_term(sub["braaksc"], sub[["TMED2","TMED10","TMED9"] + cov_cols], "TMED9")
    print(f"\n  Three-gene joint model (n={res3_tmed2['n']}):")
    print(f"  TMED2  β={res3_tmed2['beta']:.4f}  p={res3_tmed2['p']:.4f}")
    print(f"  TMED10 β={res3_tmed10['beta']:.4f}  p={res3_tmed10['p']:.4f}")
    print(f"  TMED9  β={res3_tmed9['beta']:.4f}  p={res3_tmed9['p']:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST H: Repeated 5-fold CV (5 folds × 10 repeats)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST H: Cross-validated AUC ──")
cv_rows = []
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)

for feat in ["TMED2","TMED10","TMED9","H8_ratio","PC1"]:
    sub = df_ad[["AD", feat]].dropna()
    if len(sub) < 40 or sub["AD"].nunique() < 2:
        cv_rows.append(dict(feature=feat, mean_auc=np.nan, sd_auc=np.nan))
        continue
    X = sub[[feat]].values
    y = sub["AD"].values
    aucs = []
    for train_idx, test_idx in rskf.split(X, y):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        if len(np.unique(y_te)) < 2:
            continue
        scl = StandardScaler(); X_tr_s = scl.fit_transform(X_tr); X_te_s = scl.transform(X_te)
        clf = LogisticRegression(max_iter=1000, random_state=42)
        clf.fit(X_tr_s, y_tr)
        prob = clf.predict_proba(X_te_s)[:,1]
        aucs.append(roc_auc_score(y_te, prob))
    m_auc = np.mean(aucs); s_auc = np.std(aucs)
    cv_rows.append(dict(feature=feat, mean_auc=m_auc, sd_auc=s_auc, n_folds=len(aucs)))
    print(f"  {feat:12s}  CV-AUC = {m_auc:.3f} ± {s_auc:.3f}  ({len(aucs)} folds)")

cv_df = pd.DataFrame(cv_rows)

# ══════════════════════════════════════════════════════════════════════════════
# TEST I: Permutation test (500 permutations) for H8_ratio
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST I: Permutation test ──")
N_PERM = 500
sub_perm = df_ad[["AD","H8_ratio"]].dropna()
if len(sub_perm) > 20 and sub_perm["AD"].nunique() > 1:
    obs_auc = roc_auc_score(sub_perm["AD"], sub_perm["H8_ratio"])
    if obs_auc < 0.5: obs_auc = 1 - obs_auc
    null_aucs = []
    rng = np.random.default_rng(42)
    for _ in range(N_PERM):
        y_perm = rng.permutation(sub_perm["AD"].values)
        a = roc_auc_score(y_perm, sub_perm["H8_ratio"].values)
        null_aucs.append(max(a, 1-a))
    perm_p = np.mean(np.array(null_aucs) >= obs_auc)
    print(f"  H8_ratio: observed AUC={obs_auc:.3f}  permutation p={perm_p:.4f}  (n={N_PERM} perms)")
else:
    obs_auc, null_aucs, perm_p = np.nan, [], np.nan
    print("  Insufficient data for permutation test")

# ══════════════════════════════════════════════════════════════════════════════
# TEST K: corr(Ratio, PC1) and corr(Ratio, PC2)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── TEST K: Ratio ↔ PC correlation ──")
sub_k = df[["H8_ratio","PC1","PC2"]].dropna()
if len(sub_k) > 10:
    r1, p1 = stats.pearsonr(sub_k["H8_ratio"], sub_k["PC1"])
    r2, p2 = stats.pearsonr(sub_k["H8_ratio"], sub_k["PC2"])
    print(f"  corr(H8_ratio, PC1) = {r1:.3f}  p={p1:.4f}")
    print(f"  corr(H8_ratio, PC2) = {r2:.3f}  p={p2:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle("H8: Three-protein composite biomarker", fontsize=14, fontweight="bold")

# 1. ROC curves
ax = axes[0,0]
from sklearn.metrics import roc_curve
colors = {"TMED2":"steelblue","TMED10":"darkorange","TMED9":"green",
          "H8_ratio":"red","PC1":"purple"}
for feat in ["TMED2","TMED10","TMED9","H8_ratio","PC1"]:
    r = auc_results.get(feat, {})
    if "y_true" in r:
        fpr, tpr, _ = roc_curve(r["y_true"], r["y_score"])
        ax.plot(fpr, tpr, label=f"{feat} (AUC={r['auc']:.3f})", color=colors[feat])
ax.plot([0,1],[0,1],"k--", alpha=0.5)
ax.set_xlabel("FPR"); ax.set_ylabel("TPR"); ax.set_title("ROC curves (AD vs NCI)")
ax.legend(fontsize=7); ax.set_xlim(0,1); ax.set_ylim(0,1)

# 2. CV-AUC barplot
ax = axes[0,1]
cv_plot = cv_df.dropna(subset=["mean_auc"])
if len(cv_plot) > 0:
    ax.bar(cv_plot["feature"], cv_plot["mean_auc"],
           yerr=cv_plot["sd_auc"], capsize=4,
           color=[colors.get(f,"gray") for f in cv_plot["feature"]])
    ax.axhline(0.5, color="k", linestyle="--", alpha=0.5)
    ax.set_ylabel("CV-AUC (mean ± SD)"); ax.set_title("5-fold × 10-rep CV AUC")
    ax.set_ylim(0.3, 1.0)
    ax.tick_params(axis="x", rotation=15)

# 3. ΔR² plot (Braak)
ax = axes[0,2]
br_plot = braak_df.dropna(subset=["delta_r2"])
if len(br_plot) > 0:
    ax.bar(br_plot["feature"], br_plot["delta_r2"],
           color=[colors.get(f,"gray") for f in br_plot["feature"]])
    ax.set_ylabel("ΔR² (above covariates)"); ax.set_title("Braak ~ feature ΔR²")
    ax.tick_params(axis="x", rotation=15)

# 4. Permutation null distribution
ax = axes[1,0]
if null_aucs:
    ax.hist(null_aucs, bins=30, color="lightgray", edgecolor="k", label="Null AUC")
    ax.axvline(obs_auc, color="red", linestyle="--", label=f"Observed={obs_auc:.3f}")
    ax.set_xlabel("AUC"); ax.set_ylabel("Count")
    ax.set_title(f"Permutation test (H8_ratio)\np={perm_p:.4f}")
    ax.legend(fontsize=8)

# 5. H8_ratio distribution by AD
ax = axes[1,1]
sub_box = df_ad[["AD","H8_ratio"]].dropna()
if len(sub_box) > 0:
    g0 = sub_box.loc[sub_box["AD"]==0,"H8_ratio"]
    g1 = sub_box.loc[sub_box["AD"]==1,"H8_ratio"]
    ax.boxplot([g0, g1], tick_labels=["NCI","AD"])
    ax.set_ylabel("H8_ratio (TMED9 − TMED2 − TMED10)")
    ax.set_title(f"H8_ratio by diagnosis\nNCI n={len(g0)}  AD n={len(g1)}")

# 6. PCA biplot
ax = axes[1,2]
sub_pca = df[["PC1","PC2","braaksc"]].dropna()
if len(sub_pca) > 0:
    sc = ax.scatter(sub_pca["PC1"], sub_pca["PC2"],
                    c=sub_pca["braaksc"], cmap="RdYlGn_r", s=10, alpha=0.7)
    plt.colorbar(sc, ax=ax, label="Braak")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title("PCA score plot (coloured by Braak)")

plt.tight_layout()
plt.savefig(OUT + "h8_plots.png", dpi=150)
plt.close()
print(f"\nPlot saved: {OUT}h8_plots.png")

# ══════════════════════════════════════════════════════════════════════════════
# SAVE TABLES
# ══════════════════════════════════════════════════════════════════════════════
auc_df = pd.DataFrame([
    dict(feature=k, auc=v.get("auc"), n=v.get("n"))
    for k,v in auc_results.items()
])
auc_df.to_csv(OUT + "h8_aucs.csv", index=False)
braak_df.to_csv(OUT + "h8_braak.csv", index=False)
cv_df.to_csv(OUT + "h8_cv_aucs.csv", index=False)
if len(delong_df) > 0:
    delong_df.to_csv(OUT + "h8_delong.csv", index=False)

print("\n══ H8 COMPLETE ══")
print(f"  AUC table: {OUT}h8_aucs.csv")
print(f"  Braak R²:  {OUT}h8_braak.csv")
print(f"  CV AUC:    {OUT}h8_cv_aucs.csv")
print(f"  DeLong:    {OUT}h8_delong.csv")
