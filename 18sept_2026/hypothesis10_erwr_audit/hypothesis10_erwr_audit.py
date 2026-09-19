"""
hypothesis10_erwr_audit.py
H5 / Test 5: ROSMAP eRWR — audit, normalization fix, biological association
CHRONOS / ROSMAP  —  Sep 2026

Tests
  A  Implementation audit  (synthetic perturbation — does score respond?)
  B  Normalization comparison  (per-sample min-max vs global z vs global percentile)
  C  Biological association  (eRWR ~ Braak, CERAD, cogdx)
  D  Composition-adjusted models
  E  Residualized score sensitivity
  F  Network null  (500 matched random networks)
  G  Incremental validity  (eRWR vs TMED2 alone)
"""

import numpy as np
import pandas as pd
import networkx as nx
import pickle, json, os, warnings
import scipy.stats as stats
import torch
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

DATA    = "/home/wordson22/projects/CHRONOS/DATA/"
RESULTS = "/home/wordson22/projects/CHRONOS/results/"
OUT     = "/home/wordson22/projects/CHRONOS/analysis/hypothesis10_erwr_audit/"
os.makedirs(OUT, exist_ok=True)

# ── helpers ────────────────────────────────────────────────────────────────

def bh_fdr(pvals):
    out  = np.full(len(pvals), np.nan)
    mask = ~np.isnan(pvals)
    if mask.sum() == 0: return out
    _, adj, _, _ = multipletests(np.array(pvals)[mask], method="fdr_bh")
    out[mask] = adj
    return out

def ols_term(y_s, X_df, term):
    cols = list(X_df.columns)
    Xm   = np.column_stack([np.ones(len(X_df))] + [pd.to_numeric(X_df[c], errors="coerce").values for c in cols]).astype(float)
    yv   = pd.to_numeric(y_s, errors="coerce").values.astype(float)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(yv))
    n_ok = int(mask.sum())
    if n_ok < len(cols) + 3:
        return dict(beta=np.nan, se=np.nan, p=np.nan, n=n_ok, r2=np.nan)
    Xm, yv = Xm[mask], yv[mask]
    coeffs, _, _, _ = np.linalg.lstsq(Xm, yv, rcond=None)
    resid = yv - Xm @ coeffs
    n, p  = len(yv), Xm.shape[1]
    mse   = np.sum(resid**2) / (n - p)
    try:
        cov = mse * np.linalg.inv(Xm.T @ Xm)
    except np.linalg.LinAlgError:
        return dict(beta=np.nan, se=np.nan, p=np.nan, n=n_ok, r2=np.nan)
    se_all = np.sqrt(np.diag(cov))
    idx    = (["intercept"] + cols).index(term)
    beta, se_b = coeffs[idx], se_all[idx]
    t  = beta / se_b
    pv = 2 * stats.t.sf(abs(t), df=n - p)
    yhat = Xm @ coeffs
    r2   = 1 - np.sum(resid**2) / np.sum((yv - yv.mean())**2)
    return dict(beta=beta, se=se_b, p=pv, n=n_ok, r2=r2)

def r2_full(y_s, X_df):
    cols = list(X_df.columns)
    Xm   = np.column_stack([np.ones(len(X_df))] + [pd.to_numeric(X_df[c], errors="coerce").values for c in cols]).astype(float)
    yv   = pd.to_numeric(y_s, errors="coerce").values.astype(float)
    mask = ~(np.isnan(Xm).any(axis=1) | np.isnan(yv))
    if mask.sum() < len(cols) + 3: return np.nan, np.nan
    Xm, yv = Xm[mask], yv[mask]
    coeffs, _, _, _ = np.linalg.lstsq(Xm, yv, rcond=None)
    yhat = Xm @ coeffs
    r2   = 1 - np.sum((yv - yhat)**2) / np.sum((yv - yv.mean())**2)
    n, p = len(yv), Xm.shape[1]
    aic  = n * np.log(np.sum((yv - yhat)**2) / n) + 2 * p
    return r2, aic

def spearman_safe(x, y):
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 5: return np.nan, np.nan, int(mask.sum())
    rho, p = stats.spearmanr(x[mask], y[mask])
    return rho, p, int(mask.sum())

def section(title):
    print("\n" + "="*60)
    print(f"=== {title} ===")
    print("="*60)


# ══════════════════════════════════════════════════════════════
# LOAD SHARED DATA
# ══════════════════════════════════════════════════════════════
print("Loading data...")

# Network
with open(RESULTS + "checkpoints/network.pkl", "rb") as f:
    G = pickle.load(f)
TARGET = "TMED2"
nbrs   = list(G.neighbors(TARGET))
print(f"Network: {G.number_of_nodes()} nodes, TMED2 neighbours: {len(nbrs)}")

# Expression
expr = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
with open(DATA + "ensg_to_symbol.json") as f:
    ensg_to_symbol = json.load(f)
expr.index = [ensg_to_symbol.get(i, i) for i in expr.index]
expr = expr[~expr.index.duplicated(keep="first")]

# Sample → individual mapping
bio = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_bio = bio[bio["assay"] == "rnaSeq"].dropna(subset=["individualID"])
specimen_to_ind = dict(zip(rna_bio["specimenID"], rna_bio["individualID"]))
valid_cols = {col: specimen_to_ind[col] for col in expr.columns if col in specimen_to_ind}
expr_sub = expr[[c for c in expr.columns if c in valid_cols]].copy()
expr_sub.columns = [valid_cols[c] for c in expr_sub.columns]

# Clinical
clin = pd.read_csv(RESULTS + "master_covariates.csv")
shared_inds = list(set(expr_sub.columns) & set(clin["individualID"]))
expr_sub  = expr_sub[shared_inds]
clin_sub  = clin[clin["individualID"].isin(shared_inds)].copy()
print(f"Shared individuals: {len(shared_inds)}")

# Cell composition scores (same markers as other hypotheses)
CELL_MARKERS = {
    "neuron_score": ["SYP","SYN1","SNAP25","SLC17A7","GABRB2"],
    "astro_score" : ["GFAP","AQP4"],
    "micro_score" : ["C1QA","TMEM119"],
    "oligo_score" : ["MBP","PLP1"],
}
CELL_COLS = list(CELL_MARKERS.keys())
for cname, markers in CELL_MARKERS.items():
    avail = [m for m in markers if m in expr_sub.index]
    if avail:
        z = expr_sub.loc[avail].apply(lambda r: (r - r.mean()) / r.std(), axis=1)
        score = z.mean(axis=0)
        clin_sub[cname] = clin_sub["individualID"].map(
            dict(zip(expr_sub.columns, score.reindex(expr_sub.columns))))
    else:
        clin_sub[cname] = np.nan

_age_col = next((c for c in clin_sub.columns if "age" in c.lower()), None)
COV_BASE  = ([_age_col] if _age_col else []) + ["msex", "pmi"]
print(f"Age column: {_age_col}  |  COV_BASE: {COV_BASE}")
COV_CELLS = COV_BASE + CELL_COLS

# Build subgraph
sub_nodes = set([TARGET] + nbrs)
sub       = G.subgraph([n for n in sub_nodes if G.degree(n) > 0])
sub_node_list = list(sub.nodes)
Adj = nx.to_numpy_array(sub, nodelist=sub_node_list)
node_idx  = {n: i for i, n in enumerate(sub_node_list)}
secMs     = [n for n in nbrs if n in node_idx]
sub_genes_in_expr = [g for g in sub_node_list if g in expr_sub.index]
print(f"Subgraph: {len(sub_node_list)} nodes, {len(secMs)} secMs")
print(f"Subgraph genes with expression: {len(sub_genes_in_expr)}/{len(sub_node_list)}")

# Global statistics for normalization (computed once across ALL individuals)
global_mean = {}
global_std  = {}
global_rank = {}   # gene → array of per-individual percentile ranks
all_genes_sub = [g for g in sub_node_list if g in expr_sub.index]
for g in all_genes_sub:
    vals = expr_sub.loc[g].values.astype(float)
    global_mean[g] = np.nanmean(vals)
    global_std[g]  = max(np.nanstd(vals), 1e-8)
    ranks = stats.rankdata(np.nan_to_num(vals, nan=np.nanmedian(vals)))
    global_rank[g] = dict(zip(expr_sub.columns, ranks / (len(vals) - 1)))

print("Global normalization stats computed.")

# ── eRWR core ─────────────────────────────────────────────────────────────

def rwr_torch(p0_vec, adj, expr_vec, n_prop=20, a_prop=0.1):
    A   = torch.tensor(adj, dtype=torch.float32)
    Di  = torch.diag(1.0 / torch.clamp(A.sum(0), min=1e-8))
    W   = A @ Di
    e   = torch.tensor(expr_vec, dtype=torch.float32)
    p0t = torch.tensor(p0_vec.reshape(-1, 1), dtype=torch.float32)
    Wa  = torch.diag(e) @ W
    Wa  = Wa + torch.diag(1 - torch.clamp(Wa.sum(0), max=1.0))
    pk  = p0t.clone()
    for _ in range(n_prop):
        pk = (1 - a_prop) * (Wa @ pk) + a_prop * p0t
    return pk.detach().numpy().reshape(-1)

p0_vec = np.array([1.0 if n == TARGET else 0.0 for n in sub_node_list])

def get_score(ind, norm="persample"):
    """
    norm: 'persample' | 'globalz' | 'globalrank'
    Returns mean eRWR score over secMs.
    """
    filler = float(np.median(
        [expr_sub.loc[g, ind] for g in all_genes_sub
         if not np.isnan(expr_sub.loc[g, ind])]) if all_genes_sub else 1.0)

    vec = np.zeros(len(sub_node_list))
    for i, n in enumerate(sub_node_list):
        if n in expr_sub.index and not np.isnan(expr_sub.loc[n, ind]):
            raw = float(expr_sub.loc[n, ind])
        else:
            raw = filler

        if norm == "persample":
            vec[i] = raw   # will be normalized below after full vector built
        elif norm == "globalz":
            vec[i] = (raw - global_mean.get(n, raw)) / global_std.get(n, 1.0)
        elif norm == "globalrank":
            vec[i] = global_rank.get(n, {}).get(ind, 0.5)

    if norm == "persample":
        vmin, vmax = vec.min(), vec.max()
        if vmax > vmin:
            vec = (vec - vmin) / (vmax - vmin)
        # zero TARGET to median so it doesn't drive its own score
        vec[node_idx[TARGET]] = float(np.median(vec))
    elif norm == "globalz":
        # shift to non-negative, then scale to [0,1] to keep RWR weights bounded
        vec = vec - vec.min() + 1e-6
        vmax = vec.max()
        if vmax > 1e-8:
            vec = vec / vmax
        vec[node_idx[TARGET]] = float(np.median(vec))
    elif norm == "globalrank":
        vec[node_idx[TARGET]] = float(np.median(vec))

    pk = rwr_torch(p0_vec, Adj, vec)
    return float(np.mean([pk[node_idx[g]] for g in secMs if g in node_idx]))


# ══════════════════════════════════════════════════════════════
# TEST A: IMPLEMENTATION AUDIT (synthetic perturbation)
# ══════════════════════════════════════════════════════════════
section("TEST A: IMPLEMENTATION AUDIT")

# Use a single representative individual (median TMED2 expression)
tmed2_vals = expr_sub.loc["TMED2"].values.astype(float)
rep_ind = expr_sub.columns[np.argsort(np.abs(tmed2_vals - np.nanmedian(tmed2_vals)))[0]]
print(f"\n  Representative individual: {rep_ind}")

rows_a = []
for norm in ["persample", "globalz", "globalrank"]:
    # baseline
    s_base = get_score(rep_ind, norm=norm)

    # perturb: add +2 logCPM to all TMED2-subnetwork genes in expression
    expr_backup = expr_sub.copy()
    for g in sub_genes_in_expr:
        expr_sub.loc[g, rep_ind] += 2.0
    s_high = get_score(rep_ind, norm=norm)
    expr_sub = expr_backup.copy()

    # perturb: subtract 2
    for g in sub_genes_in_expr:
        expr_sub.loc[g, rep_ind] -= 2.0
    s_low = get_score(rep_ind, norm=norm)
    expr_sub = expr_backup.copy()

    responds = (s_high != s_base) or (s_low != s_base)
    rows_a.append(dict(norm=norm, s_low=s_low, s_base=s_base, s_high=s_high,
                       delta_high=s_high - s_base, delta_low=s_low - s_base,
                       responds=responds))
    print(f"  [{norm}]  low={s_low:.6f}  base={s_base:.6f}  high={s_high:.6f}  "
          f"Δhigh={s_high-s_base:+.6f}  responds={responds}")

df_a = pd.DataFrame(rows_a)
df_a.to_csv(OUT + "h5_A_audit.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST B: NORMALIZATION COMPARISON
# ══════════════════════════════════════════════════════════════
section("TEST B: NORMALIZATION COMPARISON")

scores_by_norm = {}
for norm in ["persample", "globalz", "globalrank"]:
    print(f"\n  Computing eRWR [{norm}] for {len(shared_inds)} samples...")
    sc = {}
    for i, ind in enumerate(shared_inds):
        if i % 100 == 0: print(f"    {i}/{len(shared_inds)}")
        try:
            sc[ind] = get_score(ind, norm=norm)
        except Exception as e:
            sc[ind] = np.nan
    scores_by_norm[norm] = sc
    vals = np.array(list(sc.values()))
    print(f"  [{norm}]  mean={np.nanmean(vals):.6f}  SD={np.nanstd(vals):.6f}  "
          f"min={np.nanmin(vals):.6f}  max={np.nanmax(vals):.6f}")

# Assemble into dataframe
scores_df = pd.DataFrame({"individualID": shared_inds})
for norm, sc in scores_by_norm.items():
    scores_df[f"erwr_{norm}"] = scores_df["individualID"].map(sc)

scores_df = scores_df.merge(
    clin_sub[["individualID","cogdx","braaksc","ceradsc","presymp_group"] + COV_CELLS],
    on="individualID", how="left")

print("\n  Score SD summary:")
rows_b = []
for norm in ["persample", "globalz", "globalrank"]:
    col = f"erwr_{norm}"
    v   = scores_df[col].dropna()
    nci = scores_df[scores_df["cogdx"]==1][col].dropna()
    ad  = scores_df[scores_df["cogdx"]==4][col].dropna()
    rho, p_rho, n_rho = spearman_safe(
        scores_df["braaksc"].values.astype(float),
        scores_df[col].values.astype(float))
    _, p_mw = stats.mannwhitneyu(nci, ad, alternative="two-sided") if (len(nci)>5 and len(ad)>5) else (np.nan, np.nan)
    rows_b.append(dict(norm=norm, mean=v.mean(), SD=v.std(),
                       min=v.min(), max=v.max(),
                       rho_braak=rho, p_braak=p_rho,
                       p_NCI_vs_AD=p_mw, n=len(v)))
    print(f"  [{norm}]  SD={v.std():.6f}  rho_Braak={rho:.4f} p={p_rho:.4f}  "
          f"p_NCI_vs_AD={p_mw:.4f}")

df_b = pd.DataFrame(rows_b)
df_b.to_csv(OUT + "h5_B_normalization.csv", index=False)

# Use globalz as the primary score going forward (best expected variance)
PRIMARY_NORM = "globalz"
PRIMARY_COL  = f"erwr_{PRIMARY_NORM}"
print(f"\n  Primary score for Tests C-G: {PRIMARY_COL}")


# ══════════════════════════════════════════════════════════════
# TEST C: BIOLOGICAL ASSOCIATION
# ══════════════════════════════════════════════════════════════
section("TEST C: BIOLOGICAL ASSOCIATION")

PATH_OUTCOMES = {"braaksc": "Braak", "ceradsc": "CERAD", "cogdx": "cogdx"}
rows_c = []
print(f"\n{'Outcome':<10} {'ρ':>8} {'p':>9} {'FDR':>9} {'n':>5}")
print("-"*45)
for col, label in PATH_OUTCOMES.items():
    x = scores_df[PRIMARY_COL].values.astype(float)
    y = scores_df[col].values.astype(float)
    rho, p, n = spearman_safe(x, y)
    rows_c.append(dict(Outcome=label, rho=rho, p=p, n=n))

df_c = pd.DataFrame(rows_c)
df_c["FDR"] = bh_fdr(df_c["p"].values)
for _, r in df_c.iterrows():
    print(f"{r['Outcome']:<10} {r['rho']:>8.4f} {r['p']:>9.4f} {r['FDR']:>9.4f} {int(r['n']):>5}")

# NCI vs AD Mann-Whitney
nci_s = scores_df[scores_df["cogdx"]==1][PRIMARY_COL].dropna()
ad_s  = scores_df[scores_df["cogdx"]==4][PRIMARY_COL].dropna()
if len(nci_s) > 5 and len(ad_s) > 5:
    _, p_mw = stats.mannwhitneyu(nci_s, ad_s, alternative="two-sided")
    d_cohen = (nci_s.mean() - ad_s.mean()) / np.sqrt((nci_s.std()**2 + ad_s.std()**2)/2)
    rb = 1 - 2*stats.mannwhitneyu(nci_s, ad_s).statistic / (len(nci_s)*len(ad_s))
    pct_change = (ad_s.mean() - nci_s.mean()) / abs(nci_s.mean()) * 100
    print(f"\n  NCI vs AD: n_NCI={len(nci_s)}, n_AD={len(ad_s)}")
    print(f"  Mann-Whitney p={p_mw:.4f}, Cohen's d={d_cohen:.4f}, rank-biserial r={rb:.4f}")
    print(f"  AD vs NCI change: {pct_change:+.2f}%")
    rows_c.append(dict(Outcome="NCI_vs_AD_MW", rho=rb, p=p_mw, n=len(nci_s)+len(ad_s)))

df_c.to_csv(OUT + "h5_C_biological.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST D: COMPOSITION-ADJUSTED MODELS
# ══════════════════════════════════════════════════════════════
section("TEST D: COMPOSITION-ADJUSTED MODELS")

rows_d = []
print(f"\n{'Outcome':<8} {'Model':<28} {'β':>9} {'p':>9} {'R²':>7} {'n':>5}")
print("-"*65)
for path_col, path_label in PATH_OUTCOMES.items():
    sub = scores_df.dropna(subset=[PRIMARY_COL, path_col])

    # simple
    r1 = ols_term(sub[PRIMARY_COL], sub[[path_col]], path_col)
    rows_d.append(dict(Outcome=path_label, Model="simple",
                       beta=r1["beta"], p=r1["p"], r2=r1["r2"], n=r1["n"]))
    print(f"{path_label:<8} {'simple':<28} {r1['beta']:>9.4f} {r1['p']:>9.4f} {r1['r2']:>7.4f} {r1['n']:>5}")

    # + demographics
    sub2 = sub.dropna(subset=COV_BASE)
    r2 = ols_term(sub2[PRIMARY_COL], sub2[[path_col]+COV_BASE], path_col)
    rows_d.append(dict(Outcome=path_label, Model="+demo",
                       beta=r2["beta"], p=r2["p"], r2=r2["r2"], n=r2["n"]))
    print(f"{path_label:<8} {'+demo':<28} {r2['beta']:>9.4f} {r2['p']:>9.4f} {r2['r2']:>7.4f} {r2['n']:>5}")

    # + demographics + cells
    sub3 = sub2.dropna(subset=CELL_COLS)
    r3 = ols_term(sub3[PRIMARY_COL], sub3[[path_col]+COV_BASE+CELL_COLS], path_col)
    rows_d.append(dict(Outcome=path_label, Model="+demo+cells",
                       beta=r3["beta"], p=r3["p"], r2=r3["r2"], n=r3["n"]))
    print(f"{path_label:<8} {'+demo+cells':<28} {r3['beta']:>9.4f} {r3['p']:>9.4f} {r3['r2']:>7.4f} {r3['n']:>5}")

df_d = pd.DataFrame(rows_d)
df_d.to_csv(OUT + "h5_D_composition_adjusted.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST E: RESIDUALIZED SCORE SENSITIVITY
# ══════════════════════════════════════════════════════════════
section("TEST E: RESIDUALIZED SCORE SENSITIVITY")

sub_e = scores_df.dropna(subset=[PRIMARY_COL] + CELL_COLS).copy()
# Regress composition out of score
Xc = np.column_stack([np.ones(len(sub_e))] + [sub_e[c].values for c in CELL_COLS])
yc = sub_e[PRIMARY_COL].values
coeffs_c, _, _, _ = np.linalg.lstsq(Xc, yc, rcond=None)
sub_e["erwr_resid"] = yc - Xc @ coeffs_c

rows_e = []
print(f"\n{'Outcome':<10} {'ρ (raw)':>10} {'ρ (resid)':>10} {'p (resid)':>10} {'n':>5}")
print("-"*50)
for path_col, path_label in PATH_OUTCOMES.items():
    rho_raw, p_raw, _ = spearman_safe(
        sub_e[PRIMARY_COL].values.astype(float),
        sub_e[path_col].values.astype(float))
    rho_res, p_res, n_res = spearman_safe(
        sub_e["erwr_resid"].values.astype(float),
        sub_e[path_col].values.astype(float))
    rows_e.append(dict(Outcome=path_label, rho_raw=rho_raw, p_raw=p_raw,
                       rho_resid=rho_res, p_resid=p_res, n=n_res))
    print(f"{path_label:<10} {rho_raw:>10.4f} {rho_res:>10.4f} {p_res:>10.4f} {n_res:>5}")

df_e = pd.DataFrame(rows_e)
df_e.to_csv(OUT + "h5_E_residualized.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST F: NETWORK NULL (500 matched random networks)
# ══════════════════════════════════════════════════════════════
section("TEST F: NETWORK NULL")

print(f"\n  Running 500 random-network controls (globalz norm)...")
np.random.seed(42)
all_genes_in_net = [g for g in G.nodes if g in expr_sub.index]

# TMED2 network mean score across all individuals
tmed2_scores = np.array([scores_by_norm[PRIMARY_NORM].get(ind, np.nan) for ind in shared_inds])
tmed2_mean   = np.nanmean(tmed2_scores)

null_means = []
n_iter = 0
for perm_i in range(600):
    if len(null_means) >= 500: break
    rand_gene = np.random.choice(all_genes_in_net)
    rand_nbrs = list(G.neighbors(rand_gene))
    if len(rand_nbrs) < 3: continue
    rand_sub_nodes = set([rand_gene] + rand_nbrs)
    rand_sub  = G.subgraph([n for n in rand_sub_nodes if G.degree(n) > 0])
    rand_node_list = list(rand_sub.nodes)
    rand_Adj  = nx.to_numpy_array(rand_sub, nodelist=rand_node_list)
    rand_idx  = {n: i for i, n in enumerate(rand_node_list)}
    rand_secMs = [n for n in rand_nbrs if n in rand_idx]
    rand_p0   = np.array([1.0 if n == rand_gene else 0.0 for n in rand_node_list])
    rand_genes_in_expr = [g for g in rand_node_list if g in expr_sub.index]
    if len(rand_genes_in_expr) < 3: continue

    # compute mean score over a random subset of 50 individuals (faster)
    sample_inds = np.random.choice(shared_inds, size=min(50, len(shared_inds)), replace=False)
    sc_list = []
    for ind in sample_inds:
        filler = float(np.nanmedian(
            [expr_sub.loc[g, ind] for g in rand_genes_in_expr
             if not np.isnan(expr_sub.loc[g, ind])]))
        vec = np.zeros(len(rand_node_list))
        for i, n in enumerate(rand_node_list):
            if n in expr_sub.index and not np.isnan(expr_sub.loc[n, ind]):
                raw = float(expr_sub.loc[n, ind])
            else:
                raw = filler
            # globalz
            gm = global_mean.get(n, raw)
            gs = global_std.get(n, 1.0)
            vec[i] = (raw - gm) / gs
        vec = vec - vec.min() + 1e-6
        vmax_null = vec.max()
        if vmax_null > 1e-8:
            vec = vec / vmax_null
        if rand_gene in rand_idx:
            vec[rand_idx[rand_gene]] = float(np.median(vec))
        pk = rwr_torch(rand_p0, rand_Adj, vec)
        sc_list.append(float(np.mean([pk[rand_idx[g]] for g in rand_secMs if g in rand_idx])))

    if sc_list:
        null_means.append(np.mean(sc_list))
    n_iter += 1

null_arr  = np.array(null_means)
null_mean = np.mean(null_arr)
null_std  = np.std(null_arr)
z_null    = (tmed2_mean - null_mean) / null_std
p_null    = stats.norm.sf(abs(z_null)) * 2

print(f"  Null: mean={null_mean:.6f}, SD={null_std:.6f}, n={len(null_arr)}")
print(f"  TMED2 mean score: {tmed2_mean:.6f}")
print(f"  Z = {z_null:.3f}, p = {p_null:.4f}")

df_f = pd.DataFrame(dict(null_mean=null_arr))
df_f.to_csv(OUT + "h5_F_null.csv", index=False)
pd.DataFrame([dict(tmed2_mean=tmed2_mean, null_mean=null_mean, null_sd=null_std,
                   z=z_null, p=p_null, n_null=len(null_arr))]
             ).to_csv(OUT + "h5_F_null_summary.csv", index=False)


# ══════════════════════════════════════════════════════════════
# TEST G: INCREMENTAL VALIDITY (eRWR vs TMED2 alone)
# ══════════════════════════════════════════════════════════════
section("TEST G: INCREMENTAL VALIDITY")

# Add TMED2 expression to scores_df
if "TMED2" in expr_sub.index:
    tmed2_expr = expr_sub.loc["TMED2"]
    scores_df["TMED2_expr"] = scores_df["individualID"].map(
        lambda ind: float(tmed2_expr.loc[ind]) if ind in tmed2_expr.index else np.nan)

rows_g = []
print(f"\n{'Outcome':<8} {'Model':<30} {'β_eRWR':>9} {'p_eRWR':>9} {'R²':>7} {'ΔAIC':>8} {'n':>5}")
print("-"*75)

for path_col, path_label in PATH_OUTCOMES.items():
    sub_g = scores_df.dropna(subset=[PRIMARY_COL, "TMED2_expr", path_col] + COV_BASE + CELL_COLS).copy()
    cov_cols = COV_BASE + CELL_COLS

    # M1: covariates only
    r2_m1, aic_m1 = r2_full(sub_g[path_col], sub_g[cov_cols])

    # M2: + TMED2
    r2_m2, aic_m2 = r2_full(sub_g[path_col], sub_g[cov_cols + ["TMED2_expr"]])
    rm2 = ols_term(sub_g[path_col], sub_g[cov_cols + ["TMED2_expr"]], "TMED2_expr")

    # M3: + eRWR (no TMED2)
    r2_m3, aic_m3 = r2_full(sub_g[path_col], sub_g[cov_cols + [PRIMARY_COL]])
    rm3 = ols_term(sub_g[path_col], sub_g[cov_cols + [PRIMARY_COL]], PRIMARY_COL)

    # M4: + TMED2 + eRWR
    r2_m4, aic_m4 = r2_full(sub_g[path_col], sub_g[cov_cols + ["TMED2_expr", PRIMARY_COL]])
    rm4 = ols_term(sub_g[path_col], sub_g[cov_cols + ["TMED2_expr", PRIMARY_COL]], PRIMARY_COL)

    for label, rm, r2, aic in [
        ("M2: +TMED2",          rm2, r2_m2, aic_m2),
        ("M3: +eRWR",           rm3, r2_m3, aic_m3),
        ("M4: +TMED2+eRWR",     rm4, r2_m4, aic_m4),
    ]:
        delta_aic = aic - aic_m1
        rows_g.append(dict(Outcome=path_label, Model=label,
                           beta=rm["beta"], p=rm["p"], r2=r2,
                           delta_aic=delta_aic, n=rm["n"]))
        b_str = f"{rm['beta']:>9.4f}" if not np.isnan(rm["beta"]) else "       nan"
        p_str = f"{rm['p']:>9.4f}"    if not np.isnan(rm["p"])    else "       nan"
        print(f"{path_label:<8} {label:<30} {b_str} {p_str} {r2:>7.4f} {delta_aic:>8.2f} {rm['n']:>5}")

df_g = pd.DataFrame(rows_g)
df_g.to_csv(OUT + "h5_G_incremental.csv", index=False)

scores_df.to_csv(OUT + "h5_scores_all_norms.csv", index=False)


# ══════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════
section("DONE")
print(f"\nOutputs written to: {OUT}")
for f in ["h5_A_audit.csv","h5_B_normalization.csv","h5_C_biological.csv",
          "h5_D_composition_adjusted.csv","h5_E_residualized.csv",
          "h5_F_null.csv","h5_F_null_summary.csv","h5_G_incremental.csv",
          "h5_scores_all_norms.csv"]:
    print(f"  {f}")
