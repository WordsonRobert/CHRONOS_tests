"""
chronos_regions.py  —  CHRONOS causal graphs that learn from the confidence graph
=================================================================================

One 80x80 directed causal graph per dataset (ROSMAP, Diverse, Banner), fitted
independently. The confidence graph C is the target: it decides which edges
may exist, and every search step below is scored by how well the data
recovers C.

For each dataset the script answers three questions:

 1. SIGNAL      How well does this dataset's protein co-variation recover C?
                Score = AUROC: rank all 3,160 protein pairs by |Spearman rho|;
                AUROC = chance that a prior edge (C>0) outranks a non-edge.
                0.5 = no signal, 1.0 = data reproduces C perfectly.

 2. VARIABLES   Which trajectory variables Z (anything in the dataset except
                the 80 proteins: clinical, pathology, genetics, technical
                batch) sharpen that signal when adjusted for?
                Each Z's gain is compared with the gain from adjusting for a
                shuffled copy of Z (same values, broken link to subjects), so
                a variable is selected only if it beats noise (z > 2).

 3. REGIONS     Which subgroups of people (e.g. Braak stage V-VI, APOE e4
                carriers, age >= 90, one cohort/race) show C most strongly?
                Each region's AUROC is compared with 200 random subsets of the
                SAME size, because smaller groups are noisier. A region is
                "strong" if it beats >= 95% of those random subsets.
                Best regions are then refined once with a second variable.

Outputs (per dataset, in --out):
  {ds}_causal_W.csv      80x80 weights  W_ij = C_ij x |rho_ij|, LiNGAM direction
  {ds}_causal_beta.csv   80x80 effect sizes (regression of target on parents)
  {ds}_causal_N.csv      80x80 pairwise sample size behind each edge
  {ds}_causal_edges.csv  prior edges with W, beta, n, p, FDR q, data_supported
  {ds}_variables.csv     every candidate Z with gain vs shuffled-Z null
  {ds}_regions.csv       every region tested with AUROC vs same-size null
  {ds}_best_region_W.csv graph refitted inside the strongest region (if any)
  {ds}_report.txt        plain-English summary

Missing data: pairwise-complete everywhere (no subject is dropped because
some other protein is missing); n_ij is stored for every edge.
"""

import argparse, json, re, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm, t as tdist

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(22)

PRIOR_SCALE = 5.0
MIN_PRIOR = 0.05
MIN_PAIR_N = 20        # pairs with fewer subjects are not scored
MIN_REGION_N = 40      # smallest subgroup worth fitting a graph in
N_NULL_REGION = 200    # random same-size subsets per region
N_NULL_VAR = 10        # shuffled-Z repeats per variable
MAX_Z = 5
SIGNED = True
DESCR = {}
N_PROTEOME_PCS = 10

# --------------------------------------------------------------------------
# DATASET RECIPES: where proteomics + metadata live and how IDs chain
# --------------------------------------------------------------------------
RECIPES = {
    "ROSMAP": dict(
        prot="C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv",
        fmt="gene_pipe",
        assay="ROSMAP_assay_proteomics_TMTquantitation_metadata.csv",
        biospec="ROSMAP_biospecimen_metadata.csv",
        clin="ROSMAP_clinical.csv",
    ),
    "Diverse": dict(
        prot="n1086_residual_log2_batch.csv",
        fmt="gene_pipe",
        assay="AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata_260622.csv",
        biospec="AMP-AD_DiverseCohorts_biospecimen_metadata.csv",
        clin="AMP-AD_DiverseCohorts_individual_metadata.csv",
    ),
    "Banner": dict(
        # TMT, made by prep_banner.py from bannertmt_22batchmulticoncensus_Proteins.txt
        prot="Banner_TMT_log2ratio_GIS.csv",
        fmt="gene_pipe",
        assay="Banner_TMTquantitation_assay_metadata.csv",
        biospec="Banner_biospecimen_metadata.csv",
        clin="Banner_individual_metadata.csv",
        batch_from_channel=True,   # TMT run = prefix of batchChannel (b01..b22)
    ),
    "BannerLFQ": dict(
        # label-free, made by prep_banner.py from banner_proteomics_pfc_proteinoutput.txt
        # (robustness check only: different technology, same donors as Banner TMT)
        prot="Banner_LFQ_log2.csv",
        fmt="gene_pipe",
        assay="Banner_proteomics_assay_metadata.csv",
        biospec="Banner_biospecimen_metadata.csv",
        clin="Banner_individual_metadata.csv",
    ),
}

# columns that are identifiers / bookkeeping, never trajectory variables
NOT_Z = re.compile(r"(^unnamed|individualid|specimenid|projid|sampleid|^id$|source|"
                   r"multiconsensus|lotnumber|^platform$|fragmentation|resolution|"
                   r"tmttype|batchchannel|iscontrol|isassaycontrol|controltype|"
                   r"submissionbatch|^assay$|^organ$|^tissue$|brodmann|^species$)", re.I)
TECHNICAL = re.compile(r"(batch|pmi|ph$|study|cohort|datacontribution)", re.I)
PROTEOME = re.compile(r"^(proteomePC|celltype_)")


# ==========================================================================
# 1. LOADING + ID CHAIN
# ==========================================================================

def load_prior(path):
    P = pd.read_excel(path, index_col=0)
    proteins = [str(p) for p in P.index]
    C = P.values.astype(float) / PRIOR_SCALE
    np.fill_diagonal(C, 0)
    return proteins, C


def load_proteomics(path, fmt, proteins):
    sep = "\t" if str(path).endswith((".txt", ".tsv")) else ","
    if fmt == "gene_pipe":
        raw = pd.read_csv(path, index_col=0, sep=sep)
        raw.index = [str(i).split("|")[0].strip() for i in raw.index]
    else:
        # "auto": any rows-are-proteins table. Gene symbol from a GENE|UniProt
        # index, a gene-symbol column, or GN=SYMBOL inside a description column;
        # sample columns = numeric columns (abundances), bookkeeping dropped.
        raw = pd.read_csv(path, sep=sep, low_memory=False)
        gcol = next((c for c in raw.columns if re.fullmatch(r"(gene.?symbol|gene.?name|genes?|symbol)", str(c), re.I)), None)
        if gcol is not None:
            genes = raw[gcol].astype(str).str.split(r"[;|]").str[0].str.strip()
        else:
            dcol = next((c for c in raw.columns if re.search(r"description", str(c), re.I)), None)
            if dcol is not None:
                genes = raw[dcol].astype(str).str.extract(r"GN=([A-Za-z0-9\-]+)")[0]
            else:
                genes = raw.iloc[:, 0].astype(str).str.split("|").str[0].str.strip()
        num = raw.select_dtypes("number")
        num = num.loc[:, [c for c in num.columns if not re.search(
            r"(#|peptide|psm|score|coverage|mw|kda|length|protein group id|q.?value|pep$|calc)", str(c), re.I)]]
        if num.shape[1] < 20:
            raise ValueError(f"{path.name}: only {num.shape[1]} numeric sample columns — this is not an "
                             "abundance matrix (the Banner ProteinGroups file only has found/not-found flags)")
        num.index = genes.values
        raw = num[num.index.notna()]
    raw = raw[~raw.index.duplicated(keep="first")]
    found = [p for p in proteins if p in raw.index]
    X = raw.loc[found].T.apply(pd.to_numeric, errors="coerce")
    X.index = X.index.astype(str)
    rest = raw.drop(index=found).T.apply(pd.to_numeric, errors="coerce")
    rest.index = rest.index.astype(str)
    rest = rest.loc[:, [c for c in rest.columns if isinstance(c, str) and c not in ("", "NA", "nan")]]
    return X, found, rest


def sample_to_individual(sample_ids, assay_path, biospec_path):
    bio = pd.read_csv(biospec_path, dtype=str)
    spec2ind = dict(zip(bio["specimenID"].str.strip(), bio["individualID"].str.strip()))
    bc2spec, spec2batch = {}, {}
    if assay_path and Path(assay_path).exists():
        a = pd.read_csv(assay_path, dtype=str)
        a["specimenID"] = a["specimenID"].str.strip()
        bcol = next((c for c in a.columns if c.lower() == "batchchannel"), None)
        if bcol:
            # a batchChannel can appear in several releases; keep the one whose
            # specimen is in the biospecimen file
            for _, r in a.iterrows():
                bc = str(r[bcol]).strip()
                if bc not in bc2spec or (r["specimenID"] in spec2ind and bc2spec[bc] not in spec2ind):
                    bc2spec[bc] = r["specimenID"]
        if "batch" in a.columns:
            spec2batch = dict(zip(a["specimenID"], a["batch"].astype(str)))
    out, batch = {}, {}
    for s in sample_ids:
        spec = s if s in spec2ind else bc2spec.get(s)
        if spec in spec2ind:
            out[s] = spec2ind[spec]
            if spec in spec2batch:
                batch[s] = spec2batch[spec]
    return out, batch


def collapse_to_individuals(X, s2i, batch):
    """One row per person: replicate specimens of the same person are averaged."""
    keep = [s for s in X.index if s in s2i]
    Xk = X.loc[keep].copy()
    Xk["__ind"] = [s2i[s] for s in keep]
    b = pd.Series({s2i[s]: batch.get(s) for s in keep})
    Xi = Xk.groupby("__ind").mean()
    Xi.index.name = None
    return Xi, b.groupby(level=0).first()


def clean_covariates(clin, proteins):
    """Every non-protein, non-ID column -> numeric (text top-codes like '90+'
    -> 90; text categories -> ordered codes). Returns numeric table + label maps."""
    prot = {p.lower() for p in proteins}
    num, labels = {}, {}
    for col in clin.columns:
        if NOT_Z.search(col) or col.lower() in prot:
            continue
        s = clin[col]
        if s.notna().sum() < MIN_REGION_N or s.nunique(dropna=True) < 2:
            continue
        if s.value_counts(normalize=True, dropna=True).iloc[0] > 0.95:
            continue
        txt = s.astype(str).str.strip().str.replace(r"^(>=|>|<=|<)", "", regex=True) \
                                        .str.replace(r"\+$", "", regex=True)
        txt = txt.where(s.notna() & ~s.astype(str).str.lower().isin(
            ["nan", "na", "missing or unknown", "unknown", "not applicable", ""]))
        v = pd.to_numeric(txt, errors="coerce")
        if v.notna().sum() >= 0.8 * txt.notna().sum() and v.nunique() >= 2:
            num[col] = v
        else:
            lv = sorted(txt.dropna().unique(), key=_natural_key)
            if len(lv) > 12:
                continue  # free text / near-unique labels, not a variable
            m = {l: i for i, l in enumerate(lv)}
            num[col] = txt.map(m).astype(float)
            labels[col] = lv
    return pd.DataFrame(num, index=clin.index), labels


def _natural_key(s):
    """Sort 'Stage II' before 'Stage IV', '2' before '10'."""
    roman = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}
    key = []
    for t in re.split(r"(\d+(?:\.\d+)?|\b[IV]+\b)", str(s)):
        if not t:
            continue
        if re.fullmatch(r"\d+(?:\.\d+)?", t):
            key.append((0, float(t), ""))
        elif t in roman:
            key.append((0, float(roman[t]), ""))
        else:
            key.append((1, 0.0, t.lower()))
    return key


CELL_MARKERS = {
    "neuron":          ["SNAP25", "SYT1", "STMN2", "GAP43", "SYN1", "RBFOX3", "NEFL", "CAMK2A"],
    "astrocyte":       ["GFAP", "ALDH1L1", "AQP4", "SLC1A2", "SLC1A3", "GJA1"],
    "microglia":       ["AIF1", "CD68", "ITGAM", "CSF1R", "P2RY12", "CX3CR1", "HLA-DRA"],
    "oligodendrocyte": ["MBP", "MOG", "PLP1", "MAG", "CNP", "MOBP"],
    "endothelial":     ["CLDN5", "VWF", "PECAM1", "FLT1", "ESAM"],
}


def proteome_covariates(rest_ind, log):
    """Trajectory variables built from everything in the dataset that is NOT one
    of the 80 nodes: (a) cell-type scores = mean z-score of marker proteins,
    (b) the top principal components of the rest of the proteome, each labelled
    by the cell-type score it tracks most closely."""
    Z = (rest_ind - rest_ind.mean()) / rest_ind.std()
    cov, info = {}, []
    for ct, genes in CELL_MARKERS.items():
        g = [x for x in genes if x in Z.columns]
        if len(g) >= 2:
            cov[f"celltype_{ct}"] = Z[g].mean(axis=1, skipna=True)
            info.append(f"celltype_{ct}: mean of {', '.join(g)}")
    dense = Z.loc[:, rest_ind.notna().mean() >= 0.95].fillna(0.0)
    if dense.shape[1] >= 50:
        U, S, _ = np.linalg.svd(dense.values - dense.values.mean(0), full_matrices=False)
        ve = S ** 2 / np.sum(S ** 2)
        for k in range(min(N_PROTEOME_PCS, U.shape[1])):
            pc = pd.Series(U[:, k] * S[k], index=dense.index)
            lab = ""
            best = 0
            for ct in CELL_MARKERS:
                c = f"celltype_{ct}"
                if c in cov:
                    r = pc.corr(cov[c], method="spearman")
                    if abs(r) > abs(best):
                        best, lab = r, ct
            name = f"proteomePC{k + 1}"
            cov[name] = pc
            info.append(f"{name}: {100 * ve[k]:.1f}% of rest-of-proteome variance; "
                        f"tracks {lab} (rho={best:+.2f})" if lab else
                        f"{name}: {100 * ve[k]:.1f}% of variance")
        log(f"  Rest-of-proteome: {dense.shape[1]} well-measured proteins -> {N_PROTEOME_PCS} PCs + "
            f"{sum(1 for k in cov if k.startswith('celltype'))} cell-type scores")
    return pd.DataFrame(cov), info


# ==========================================================================
# 2. FAST PAIRWISE-COMPLETE SPEARMAN + SIGNAL SCORE
# ==========================================================================

def spearman_pairwise(Xa):
    """Spearman rho for every protein pair using the subjects that have both.
    (Ranks computed per protein over its observed values, then pairwise
    Pearson on ranks with mask algebra — fast and very close to exact.)"""
    M = np.isfinite(Xa).astype(float)
    R = np.where(M > 0, 0.0, 0.0)
    for j in range(Xa.shape[1]):
        ok = M[:, j] > 0
        if ok.sum():
            R[ok, j] = rankdata(Xa[ok, j])
    N = M.T @ M
    S1 = R.T @ M
    S2 = (R ** 2).T @ M
    SXY = R.T @ R
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = SXY - S1 * S1.T / N
        vx = S2 - S1 ** 2 / N
        rho = cov / np.sqrt(vx * vx.T)
    rho[~np.isfinite(rho)] = 0.0
    np.fill_diagonal(rho, 0.0)
    return np.clip(rho, -1, 1), N.astype(int)


def auroc_vs_prior(rho, N, edge_mask, iu):
    """AUROC of rho separating prior edges from non-edges (pairs with n>=MIN_PAIR_N).
    Signed by default: confidence-graph edges are overwhelmingly 'rise together'
    links, and in both ROSMAP and Diverse signed rho recovers C better than |rho|."""
    a = rho[iu] if SIGNED else np.abs(rho[iu])
    e = edge_mask[iu]; ok = N[iu] >= MIN_PAIR_N
    a, e = a[ok], e[ok]
    if e.sum() < 5 or (~e).sum() < 5:
        return np.nan
    r = rankdata(a)
    n1, n0 = e.sum(), (~e).sum()
    return (r[e].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def residualise(Xa, Z):
    """Remove the linear effect of Z from every protein, pairwise-complete:
    each protein uses every subject where it and all of Z are observed.
    Subjects with incomplete Z become NaN (never mixed raw/adjusted)."""
    out = np.full_like(Xa, np.nan)
    zok = np.all(np.isfinite(Z), axis=1)
    if zok.sum() < MIN_REGION_N:
        return out
    D = np.column_stack([np.ones(zok.sum()), Z[zok]])
    Xz = Xa[zok]
    for j in range(Xa.shape[1]):
        ok = np.isfinite(Xz[:, j])
        if ok.sum() < MIN_PAIR_N:
            continue
        b, *_ = np.linalg.lstsq(D[ok], Xz[ok, j], rcond=None)
        col = np.full(Xz.shape[0], np.nan)
        col[ok] = Xz[ok, j] - D[ok] @ b
        out[np.where(zok)[0], j] = col
    return out


def design(cov, cols, labels):
    """Numeric design for adjustment: categorical -> one-hot (drop first)."""
    parts = []
    for c in cols:
        v = cov[c].values
        if c in labels:
            lv = np.unique(v[np.isfinite(v)])
            for l in lv[1:]:
                d = (v == l).astype(float); d[~np.isfinite(v)] = np.nan
                parts.append(d)
        else:
            parts.append(v.astype(float))
    return np.column_stack(parts) if parts else np.empty((len(cov), 0))


# ==========================================================================
# 3. VARIABLE (Z) SELECTION — gain vs shuffled-Z null
# ==========================================================================

def select_variables(Xa, cov, labels, edge_mask, iu, log):
    rho0, N0 = spearman_pairwise(Xa)
    base = auroc_vs_prior(rho0, N0, edge_mask, iu)
    log(f"  Signal with no adjustment: AUROC = {base:.4f}")
    selected, rows, cur = [], [], base
    remaining = list(cov.columns)
    for step in range(MAX_Z):
        best = None
        for c in remaining:
            Z = design(cov, selected + [c], labels)
            s = auroc_vs_prior(*spearman_pairwise(residualise(Xa, Z)), edge_mask, iu)
            nulls = []
            for _ in range(N_NULL_VAR):
                covp = cov.copy()
                v = covp[c].values.copy(); ok = np.isfinite(v)
                v[ok] = RNG.permutation(v[ok]); covp[c] = v
                Zp = design(covp, selected + [c], labels)
                nulls.append(auroc_vs_prior(*spearman_pairwise(residualise(Xa, Zp)), edge_mask, iu))
            nulls = np.array(nulls)
            z = (s - nulls.mean()) / (nulls.std() + 1e-9)
            rows.append(dict(step=step + 1, variable=c, adjusted_for=" + ".join(selected + [c]),
                             auroc=round(s, 5), null_mean=round(nulls.mean(), 5),
                             gain_vs_current=round(s - cur, 5), z_vs_shuffled=round(z, 2),
                             n_subjects=int(np.all(np.isfinite(Z), axis=1).sum()),
                             kind=("technical" if TECHNICAL.search(c) else
                                   "proteome-derived" if PROTEOME.search(c) else "clinical/pathology")))
            if z > 2 and s > cur + 1e-4 and (best is None or s > best[0]):
                best = (s, c, z)
        if best is None:
            break
        cur = best[0]; selected.append(best[1]); remaining.remove(best[1])
        d = DESCR.get(best[1], "")
        log(f"  + {best[1]:<28} AUROC -> {cur:.4f}   (z vs shuffled = {best[2]:.1f})" + (f"   [{d}]" if d else ""))
    if not selected:
        log("  No variable beat its shuffled null (z > 2)")
    return selected, pd.DataFrame(rows), base, cur


# ==========================================================================
# 4. REGION SEARCH — subgroup AUROC vs random same-size subsets
# ==========================================================================

def candidate_regions(cov, labels):
    regs = []
    for c in cov.columns:
        if re.search(r"batch|^proteomePC", c, re.I):
            continue  # batches aren't biology; PCs aren't interpretable subgroups
        v = cov[c].values
        ok = np.isfinite(v)
        if c in labels or len(np.unique(v[ok])) <= 6:
            for l in np.unique(v[ok]):
                name = labels[c][int(l)] if c in labels else _fmt(l)
                regs.append((f"{c} = {name}", v == l))
        else:
            q1, q2 = np.nanquantile(v, [1 / 3, 2 / 3])
            med = np.nanmedian(v)
            regs += [(f"{c} <= {_fmt(q1)} (low third)", ok & (v <= q1)),
                     (f"{c} {_fmt(q1)}-{_fmt(q2)} (middle third)", ok & (v > q1) & (v <= q2)),
                     (f"{c} > {_fmt(q2)} (high third)", ok & (v > q2)),
                     (f"{c} <= {_fmt(med)} (lower half)", ok & (v <= med)),
                     (f"{c} > {_fmt(med)} (upper half)", ok & (v > med))]
    return [(n, m) for n, m in regs if MIN_REGION_N <= m.sum() <= 0.9 * len(m)]


def _fmt(x):
    return str(int(x)) if float(x).is_integer() else f"{x:.3g}"


_null_cache = {}
def null_for_size(Xa, n, edge_mask, iu):
    if n not in _null_cache:
        vals = []
        for _ in range(N_NULL_REGION):
            idx = RNG.choice(Xa.shape[0], n, replace=False)
            vals.append(auroc_vs_prior(*spearman_pairwise(Xa[idx]), edge_mask, iu))
        _null_cache[n] = np.array(vals)
    return _null_cache[n]


def score_region(Xa, mask, edge_mask, iu):
    n = int(mask.sum())
    s = auroc_vs_prior(*spearman_pairwise(Xa[mask]), edge_mask, iu)
    # bucket sizes to the nearest 10 so the null can be reused
    null = null_for_size(Xa, max(MIN_REGION_N, int(round(n / 10.0)) * 10), edge_mask, iu)
    p = (1 + np.sum(null >= s)) / (1 + len(null))
    return n, s, null.mean(), null.std(), p


def search_regions(Xa, cov, labels, edge_mask, iu, log):
    rows = []
    regs = candidate_regions(cov, labels)
    log(f"  Testing {len(regs)} single-variable regions ...")
    for name, m in regs:
        n, s, mu, sd, p = score_region(Xa, m, edge_mask, iu)
        rows.append(dict(region=name, depth=1, n=n, auroc=s, null_mean=mu, null_sd=sd,
                         z=(s - mu) / (sd + 1e-9), p_vs_random_same_size=p))
    df = pd.DataFrame(rows)
    # refine the 5 best single regions with a second variable
    top = df.sort_values("z", ascending=False).head(5)
    masks = dict(regs)
    rows2, seen = [], set()
    for _, r in top.iterrows():
        m1 = masks[r.region]
        var1 = r.region.split(" ")[0]
        for name, m2 in regs:
            if name.split(" ")[0] == var1:
                continue
            key = frozenset([r.region, name])
            if key in seen:
                continue
            seen.add(key)
            m = m1 & m2
            if m.sum() < MIN_REGION_N:
                continue
            n, s, mu, sd, p = score_region(Xa, m, edge_mask, iu)
            rows2.append(dict(region=f"{r.region}  AND  {name}", depth=2, n=n, auroc=s,
                              null_mean=mu, null_sd=sd, z=(s - mu) / (sd + 1e-9),
                              p_vs_random_same_size=p))
            masks[f"{r.region}  AND  {name}"] = m
    df = pd.concat([df, pd.DataFrame(rows2)], ignore_index=True)
    # Benjamini-Hochberg over all regions tested
    df["q_fdr"] = bh(df["p_vs_random_same_size"].values)
    df["strong"] = (df["q_fdr"] < 0.10) & (df["z"] > 2)
    df = df.sort_values("z", ascending=False).reset_index(drop=True)
    return df, masks


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ==========================================================================
# 5. FINAL GRAPH (prior-constrained, LiNGAM direction, pairwise-complete)
# ==========================================================================

def lingam_direction(x, y):
    x = (x - x.mean()) / (x.std() + 1e-12); y = (y - y.mean()) / (y.std() + 1e-12)
    e_xy = y - np.mean(x * y) * x
    e_yx = x - np.mean(x * y) * y
    return abs(np.mean(e_yx * y ** 2) - np.mean(e_yx) * np.mean(y ** 2)) - \
           abs(np.mean(e_xy * x ** 2) - np.mean(e_xy) * np.mean(x ** 2))


def build_graph(Xa, C):
    rho, N = spearman_pairwise(Xa)
    m = Xa.shape[1]
    W = np.zeros((m, m)); P = np.ones((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            n = N[i, j]
            if n >= 4 and abs(rho[i, j]) < 1:
                tt = rho[i, j] * np.sqrt((n - 2) / (1 - rho[i, j] ** 2))
                P[i, j] = P[j, i] = 2 * tdist.sf(abs(tt), n - 2)
            if C[i, j] < MIN_PRIOR and C[j, i] < MIN_PRIOR:
                continue
            r = abs(rho[i, j])
            if n < MIN_PAIR_N:
                if C[i, j] >= C[j, i]: W[i, j] = C[i, j] * 0.1
                else: W[j, i] = C[j, i] * 0.1
                continue
            ok = np.isfinite(Xa[:, i]) & np.isfinite(Xa[:, j])
            s = lingam_direction(Xa[ok, i], Xa[ok, j])
            fwd = s > 0
            if fwd and C[i, j] >= MIN_PRIOR: W[i, j] = C[i, j] * r
            elif (not fwd) and C[j, i] >= MIN_PRIOR: W[j, i] = C[j, i] * r
            elif fwd: W[j, i] = C[j, i] * r * 0.3
            else: W[i, j] = C[i, j] * r * 0.3
    return W, rho, N, P


def compute_beta(Xa, W):
    m = Xa.shape[1]; B = np.zeros((m, m))
    for j in range(m):
        par = [i for i in range(m) if W[i, j] > 0]
        if not par: continue
        A = np.column_stack([Xa[:, j], Xa[:, par]])
        ok = np.all(np.isfinite(A), axis=1)
        if ok.sum() < len(par) + 10:
            # too few people have target + ALL parents: fall back to one
            # parent at a time (pairwise) so beta isn't silently zero
            for i in par:
                o = np.isfinite(Xa[:, i]) & np.isfinite(Xa[:, j])
                if o.sum() >= MIN_PAIR_N:
                    B[i, j] = np.polyfit(Xa[o, i], Xa[o, j], 1)[0]
            continue
        D = np.column_stack([np.ones(ok.sum()), Xa[ok][:, par]])
        coef, *_ = np.linalg.lstsq(D, Xa[ok, j], rcond=None)
        B[par, j] = coef[1:]
    return B


def expand(mat, found, proteins, fill=0.0):
    idx = {p: i for i, p in enumerate(proteins)}
    out = np.full((len(proteins), len(proteins)), fill, dtype=float)
    fi = [idx[p] for p in found]
    out[np.ix_(fi, fi)] = mat
    return out


def edge_table(proteins, C, W80, B80, N80, P80):
    rows = []
    for i in range(len(proteins)):
        for j in range(len(proteins)):
            if i != j and C[i, j] >= MIN_PRIOR:
                rows.append(dict(source=proteins[i], target=proteins[j], prior_C=round(C[i, j], 3),
                                 W=round(W80[i, j], 4), beta=round(B80[i, j], 4),
                                 n_pairwise=int(N80[i, j]), p_corr=P80[i, j]))
    df = pd.DataFrame(rows)
    df["q_fdr"] = bh(df["p_corr"].values)
    df["data_supported"] = ((df["q_fdr"] < 0.05) & (df["W"] > 0) & (df["n_pairwise"] >= MIN_PAIR_N)).astype(int)
    return df.sort_values("W", ascending=False)


# ==========================================================================
# 6. MAIN
# ==========================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(RECIPES))
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--prior", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--prot", default=None, help="override proteomics file")
    ap.add_argument("--fmt", default=None, choices=["gene_pipe", "auto"])
    ap.add_argument("--exclude_individuals", default=None,
                    help="text file of individualIDs to drop (e.g. ROSMAP donors when fitting Diverse)")
    ap.add_argument("--fast", action="store_true", help="fewer null draws (quick test)")
    a = ap.parse_args()
    global N_NULL_REGION, N_NULL_VAR
    if a.fast:
        N_NULL_REGION, N_NULL_VAR = 50, 4

    ds, D, R = a.dataset, Path(a.data_dir), dict(RECIPES[a.dataset])
    if a.prot: R["prot"] = a.prot
    if a.fmt: R["fmt"] = a.fmt
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 70); log(f"CHRONOS — {ds}"); log("=" * 70)
    proteins, C = load_prior(a.prior)
    if not R["prot"]:
        sys.exit(f"{ds}: no abundance matrix set. Download it and pass --prot.")
    X, found, rest = load_proteomics(D / R["prot"], R["fmt"], proteins)
    log(f"Proteins: {len(found)}/{len(proteins)} of the confidence-graph nodes found")

    s2i, batch = sample_to_individual(list(X.index), D / R["assay"], D / R["biospec"])
    if R.get("batch_from_channel"):
        batch = {s: s.split(".")[0] for s in s2i}
    log(f"Samples: {X.shape[0]}  |  matched to a person: {len(s2i)}")
    Xi, b = collapse_to_individuals(X, s2i, batch)
    rest_i, _ = collapse_to_individuals(rest, s2i, batch)
    if a.exclude_individuals:
        ex = set(Path(a.exclude_individuals).read_text().split())
        before = len(Xi); Xi = Xi.loc[~Xi.index.isin(ex)]; rest_i = rest_i.loc[Xi.index]
        log(f"Excluded {before - len(Xi)} people shared with another dataset")
    log(f"People (one row each, replicates averaged): {len(Xi)}")

    clin = pd.read_csv(D / R["clin"], dtype=str)
    idc = next(c for c in clin.columns if c.lower() == "individualid")
    clin[idc] = clin[idc].str.strip()
    clin = clin.drop_duplicates(idc).set_index(idc).reindex(Xi.index)
    clin["batch"] = b.reindex(Xi.index).values
    cov, labels = clean_covariates(clin, proteins)
    pcov, pinfo = proteome_covariates(rest_i.reindex(Xi.index), log)
    for line in pinfo:
        log("    " + line)
        k, v = line.split(": ", 1)
        DESCR[k] = v
    cov = pd.concat([cov, pcov.reindex(Xi.index)], axis=1)
    (out / f"{ds}_individuals.txt").write_text("\n".join(Xi.index))
    log(f"Trajectory-variable candidates (everything except the 80 proteins): {len(cov.columns)}")
    log("  " + ", ".join(cov.columns))

    Xa = Xi[found].values.astype(float)
    fi = [proteins.index(p) for p in found]
    Cs = C[np.ix_(fi, fi)]
    edge_mask = (np.maximum(Cs, Cs.T) >= MIN_PRIOR)
    iu = np.triu_indices(len(found), 1)
    rho_full, Nf = spearman_pairwise(Xa)
    log(f"Median people per protein pair: {int(np.median(Nf[iu]))}  "
        f"(people with ALL {len(found)} proteins: {int(np.all(np.isfinite(Xa), 1).sum())})")

    log("\n[1-2] SIGNAL + TRAJECTORY VARIABLES")
    Zstar, vardf, base, adj = select_variables(Xa, cov, labels, edge_mask, iu, log)
    vardf.to_csv(out / f"{ds}_variables.csv", index=False)

    Xfin = residualise(Xa, design(cov, Zstar, labels)) if Zstar else Xa
    log("\n[3] REGIONS (who shows the confidence graph most strongly, AFTER adjusting for Z*)")
    regdf, masks = search_regions(Xfin, cov, labels, edge_mask, iu, log)
    regdf.round(5).to_csv(out / f"{ds}_regions.csv", index=False)
    log(f"  Whole cohort AUROC after adjustment {adj:.4f}")
    for _, r in regdf.head(8).iterrows():
        flag = "STRONG" if r.strong else "      "
        log(f"  {flag} AUROC {r.auroc:.4f}  (random same-size {r.null_mean:.4f})  z={r.z:5.2f}  n={r.n:<4} {r.region}")
    weak = regdf.sort_values("z").head(3)
    log("  Weakest regions (C least visible):")
    for _, r in weak.iterrows():
        log(f"         AUROC {r.auroc:.4f}  (random same-size {r.null_mean:.4f})  z={r.z:5.2f}  n={r.n:<4} {r.region}")

    log("\n[4] FINAL GRAPH (whole cohort, adjusted for selected variables)")
    W, rho, N, P = build_graph(Xfin, C[np.ix_(fi, fi)])
    B = compute_beta(Xfin, W)
    W80, B80 = expand(W, found, proteins), expand(B, found, proteins)
    N80, P80 = expand(N, found, proteins), expand(P, found, proteins, fill=1.0)
    names = dict(index=proteins, columns=proteins)
    pd.DataFrame(W80, **names).to_csv(out / f"{ds}_causal_W.csv")
    pd.DataFrame(B80, **names).to_csv(out / f"{ds}_causal_beta.csv")
    pd.DataFrame(N80.astype(int), **names).to_csv(out / f"{ds}_causal_N.csv")
    et = edge_table(proteins, C, W80, B80, N80, P80)
    et.to_csv(out / f"{ds}_causal_edges.csv", index=False)
    testable = (et.n_pairwise >= MIN_PAIR_N).sum()
    log(f"  Prior edges: {len(et)} | testable (n>={MIN_PAIR_N}): {testable} | "
        f"non-zero W: {(et.W > 0).sum()} | data-supported (FDR<0.05): {et.data_supported.sum()}")

    strong = regdf[regdf.strong]
    best_region = None
    if len(strong):
        best_region = strong.iloc[0].region
        m = masks[best_region]
        Wr, *_ = build_graph(Xfin[m], C[np.ix_(fi, fi)])
        pd.DataFrame(expand(Wr, found, proteins), **names).to_csv(out / f"{ds}_best_region_W.csv")
        log(f"  Best region graph saved: {best_region}  (n={int(m.sum())})")

    summ = dict(dataset=ds, n_people=int(len(Xa)), proteins_found=len(found),
                auroc_raw=round(base, 4), auroc_adjusted=round(adj, 4), Z_star=Zstar,
                n_strong_regions=int(regdf.strong.sum()), best_region=best_region,
                top_regions=regdf.head(5)[["region", "n", "auroc", "z", "q_fdr"]].round(4).to_dict("records"),
                data_supported_edges=int(et.data_supported.sum()), prior_edges=int(len(et)))
    (out / f"{ds}_summary.json").write_text(json.dumps(summ, indent=2, default=str))
    (out / f"{ds}_report.txt").write_text("\n".join(lines))
    log("\nDone.")


if __name__ == "__main__":
    main()
