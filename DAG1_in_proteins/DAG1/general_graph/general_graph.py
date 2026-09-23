"""
ONE general 80x80 causal graph over the proteins (all subjects pooled, no trajectory conditioning).

Model (linear SEM):   X = B X + Gamma Z + eps ,  diag(B)=0
  X = 80 protein levels (expr_z, within-dataset z-score)
  Z = 14 trajectory/context variables (+ dataset dummies)  -> nuisance
  B = 80x80 weighted DIRECTED protein->protein graph  (B[i,j] = effect of j on i,
      i.e. j -> i). Continuous weights, asymmetry allowed, cycles allowed.

Estimated node-wise: for each protein i, regress X_i on the other 79 proteins + Z +
dataset dummies (elastic net, CV). Coefficients on proteins = row i of B; on Z = row i
of Gamma. The literature matrix C is NOT used here — compared afterwards
(general_graph_compare.py).

Input: aligned_XZ.csv in this folder — one row per subject (1,570 = ROSMAP 400 + Diverse 980
+ Banner 190): 80 protein columns (expr_z), 14 Z columns (harmonised as in DATA/), dataset.
It was originally assembled from the raw AMP-AD Synapse files; it carries the same values as
DATA/ (same subjects, same expr_z and harmonisation). Run from anywhere:
  python DAG1/general_graph/general_graph.py
Writes B_causal_graph.csv, Gamma_Z_to_X.csv, node_R2.csv, _fit_summary.json next to this file.
"""
import warnings, json, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from pathlib import Path
from sklearn.linear_model import ElasticNetCV

HERE = Path(__file__).resolve().parent
OUT = HERE
DATASETS = ["ROSMAP", "Diverse", "Banner"]
ZVARS = ["sex", "apoe_genotype", "cts_mmse30_lv", "cts_mmse30_first_ad_dx", "pmi", "braaksc",
         "ceradsc", "amyThal", "amyA", "amyCerad", "amyAny", "reag", "PlaqueTotal", "TangleTotal"]
# (apoe_genotype in aligned_XZ.csv is already the e4 dosage 0/1/2; sex is 1=male/0=female)

aligned = pd.read_csv(HERE / "aligned_XZ.csv", index_col=0)
proteins = [c for c in aligned.columns if c not in ZVARS + ["dataset"]]
X = aligned[proteins].astype(float)
Z = aligned[ZVARS].astype(float)
dsser = aligned["dataset"]
# standardise Z across pooled data (ignoring NaN), then NaN->0 (neutral)
Zs = (Z - Z.mean()) / Z.std(ddof=0); Zs = Zs.fillna(0.0)
# dataset dummies (Diverse, Banner; ROSMAP = reference)
D = pd.get_dummies(dsser).reindex(columns=DATASETS).astype(float)
Dext = D[["Diverse", "Banner"]]
Xf = X.fillna(0.0)                                           # predictors: NaN protein -> 0 (=mean)
print(f"pooled subjects: {X.shape[0]}  (ROSMAP {sum(dsser=='ROSMAP')}, Diverse {sum(dsser=='Diverse')}, Banner {sum(dsser=='Banner')})")
print(f"X: {X.shape}  Z(std): {Zs.shape}  dummies: {Dext.shape}")

# ---- node-wise elastic net -> B, Gamma ----
Pn = len(proteins)
B = np.zeros((Pn, Pn)); G = np.zeros((Pn, len(ZVARS))); r2 = np.zeros(Pn)
Zmat = Zs.values; Dmat = Dext.values; Xfm = Xf.values; Xtrue = X.values
t0 = time.time()
for i in range(Pn):
    yobs = ~np.isnan(Xtrue[:, i])                              # rows where target observed
    y = Xtrue[yobs, i]
    predcols = [j for j in range(Pn) if j != i]
    Xd = np.hstack([Xfm[yobs][:, predcols], Zmat[yobs], Dmat[yobs]])
    m = ElasticNetCV(l1_ratio=[.5, .8], n_alphas=40, cv=3, max_iter=5000, n_jobs=-1)
    m.fit(Xd, y)
    coef = m.coef_
    B[i, predcols] = coef[:len(predcols)]
    G[i, :] = coef[len(predcols):len(predcols) + len(ZVARS)]
    r2[i] = m.score(Xd, y)
    if i % 20 == 0: print(f"  node {i}/{Pn}  elapsed {time.time()-t0:.0f}s")
Bdf = pd.DataFrame(B, index=proteins, columns=proteins)        # B[target i, source j] = j->i
Gdf = pd.DataFrame(G, index=proteins, columns=ZVARS)
Bdf.to_csv(OUT / "B_causal_graph.csv"); Gdf.to_csv(OUT / "Gamma_Z_to_X.csv")
pd.Series(r2, index=proteins, name="node_R2").to_csv(OUT / "node_R2.csv")
nnz = (B != 0).sum()
print(f"\nB: {nnz} nonzero directed edges / {Pn*(Pn-1)} possible ({100*nnz/(Pn*(Pn-1)):.1f}%)")
print(f"median node R2: {np.median(r2):.3f}")
json.dump(dict(n_subjects=int(X.shape[0]), nonzero_edges=int(nnz), median_R2=float(np.median(r2))),
          open(OUT / "_fit_summary.json", "w"), indent=2)
