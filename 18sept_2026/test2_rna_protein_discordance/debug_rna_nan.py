"""
debug_rna_nan.py — figure out why RNA AD regression returns NaN
"""
import numpy as np
import pandas as pd
from scipy import stats
import warnings; warnings.filterwarnings("ignore")

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

GENES = {"TMED2":"ENSG00000086598","TMED10":"ENSG00000170348","TMED9":"ENSG00000184840"}
CELL_MARKERS = {
    "neuron_score": ["ENSG00000102003","ENSG00000067715","ENSG00000132639","ENSG00000008056","ENSG00000157542"],
    "astro_score":  ["ENSG00000131095","ENSG00000171885"],
    "micro_score":  ["ENSG00000204472","ENSG00000138185"],
    "oligo_score":  ["ENSG00000197971","ENSG00000123560"],
}
CELL_COLS = list(CELL_MARKERS.keys())

expr     = pd.read_csv(DATA+"ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA+"ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA+"ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA+"ROSMAP_assay_rnaSeq_metadata.csv")

clinical["age_at_visit_max"] = pd.to_numeric(clinical["age_at_visit_max"].replace("90+","90"), errors="coerce")
for col in ["braaksc","msex","pmi","cogdx"]: clinical[col] = pd.to_numeric(clinical[col], errors="coerce")
clinical["AD"] = np.nan
clinical.loc[clinical["cogdx"]==1,"AD"] = 0
clinical.loc[clinical["cogdx"]>=4,"AD"] = 1

rna_bio  = biospec[["individualID","specimenID"]].drop_duplicates()
rna_link = rna_meta[["specimenID","rnaBatch"]].merge(rna_bio, on="specimenID", how="left")
rna_link = rna_link[rna_link["specimenID"].isin(expr.columns)]
rna_link = rna_link.merge(clinical[["individualID","braaksc","msex","pmi","age_at_visit_max","cogdx","AD"]], on="individualID", how="left")
rna_link = rna_link.dropna(subset=["individualID"])
rna_link["rnaBatch_cat"] = pd.Categorical(rna_link["rnaBatch"]).codes

def cell_score(markers):
    avail = [m for m in markers if m in expr.index]
    if not avail: return pd.Series(np.nan, index=expr.columns)
    z = expr.loc[avail].apply(lambda r:(r-r.mean())/r.std(), axis=1)
    return z.mean(axis=0)

for cname, markers in CELL_MARKERS.items():
    scores = cell_score(markers)
    rna_link[cname] = rna_link["specimenID"].map(scores)

for name, eid in GENES.items():
    if eid in expr.index:
        rna_link[f"rna_{name}"] = rna_link["specimenID"].map(expr.loc[eid])

RNA_COVS_CELLS = ["AD","age_at_visit_max","msex","pmi","rnaBatch_cat"] + CELL_COLS
rna_ad = rna_link.dropna(subset=RNA_COVS_CELLS)
rna_ad = rna_ad[rna_ad["AD"].isin([0,1])]

print(f"n after dropna: {len(rna_ad)}")
print(f"rnaBatch_cat unique values: {sorted(rna_ad['rnaBatch_cat'].unique())}")
print(f"AD value counts: {rna_ad['AD'].value_counts().to_dict()}")
print(f"cell score NaN counts: {rna_ad[CELL_COLS].isna().sum().to_dict()}")

# check rank of design matrix
col = "rna_TMED2"
sub = rna_ad.dropna(subset=[col])
Xm = np.column_stack([np.ones(len(sub))] + [sub[c].values for c in RNA_COVS_CELLS])
print(f"\nDesign matrix shape: {Xm.shape}")
print(f"Matrix rank: {np.linalg.matrix_rank(Xm)}")
print(f"Expected rank (full): {Xm.shape[1]}")

# check for constant/duplicate columns
for i, cname in enumerate(["intercept"]+RNA_COVS_CELLS):
    col_vals = Xm[:,i]
    print(f"  col {cname}: min={col_vals.min():.3f} max={col_vals.max():.3f} nunique={len(np.unique(col_vals))}")
