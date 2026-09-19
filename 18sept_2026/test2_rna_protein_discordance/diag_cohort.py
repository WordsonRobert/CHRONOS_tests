"""
diag_cohort.py — diagnose why RNA and protein cohorts have 0 samples
"""
import pandas as pd
import numpy as np
import os

DATA = "/home/wordson22/projects/CHRONOS/DATA/"

print("=== FILES IN DATA DIR ===")
for f in sorted(os.listdir(DATA)):
    print(" ", f)

print("\n=== LOADING ===")
expr     = pd.read_csv(DATA + "ROSMAP_DLPFC_logCPM.tsv", sep="\t", index_col=0)
clinical = pd.read_csv(DATA + "ROSMAP_clinical.csv")
biospec  = pd.read_csv(DATA + "ROSMAP_biospecimen_metadata.csv")
rna_meta = pd.read_csv(DATA + "ROSMAP_assay_rnaSeq_metadata.csv")

print(f"expr shape: {expr.shape}  — first 3 col names: {list(expr.columns[:3])}")
print(f"clinical cols: {list(clinical.columns[:10])}")
print(f"biospec cols:  {list(biospec.columns)}")
print(f"rna_meta cols: {list(rna_meta.columns)}")

print("\n=== biospec assay values ===")
if "assay" in biospec.columns:
    print(biospec["assay"].value_counts().head(10))
else:
    print("NO 'assay' column in biospec")
    print("biospec columns:", list(biospec.columns))

print("\n=== rna_meta specimenID overlap with expr.columns ===")
if "specimenID" in rna_meta.columns:
    overlap = rna_meta["specimenID"].isin(expr.columns).sum()
    print(f"rna_meta specimenIDs in expr cols: {overlap} / {len(rna_meta)}")
    print("rna_meta first 3 specimenIDs:", list(rna_meta["specimenID"].head(3)))
    print("expr first 3 col names:       ", list(expr.columns[:3]))
else:
    print("NO specimenID in rna_meta")

print("\n=== biospec individualID linkage ===")
print("biospec first 5 rows:")
print(biospec.head())

print("\n=== clinical individualID ===")
print(f"clinical n rows: {len(clinical)}")
print(f"clinical 'individualID' col exists: {'individualID' in clinical.columns}")
if "individualID" in clinical.columns:
    print(f"clinical individualID sample: {list(clinical['individualID'].head(3))}")

print("\n=== PROTEOMICS ===")
prot_raw = pd.read_csv(
    DATA + "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv",
    index_col=0)
print(f"prot_raw shape: {prot_raw.shape}")
print(f"prot_raw first 3 col names: {list(prot_raw.columns[:3])}")
print(f"prot_raw first 3 row names: {list(prot_raw.index[:3])}")

# check if any proteomics metadata
print("\n=== Proteomics metadata candidates ===")
prot_meta_names = [f for f in os.listdir(DATA) if "proteom" in f.lower() or "TMT" in f or "tmt" in f.lower()]
print("candidates:", prot_meta_names)
for fn in prot_meta_names:
    pm = pd.read_csv(DATA + fn)
    print(f"\n  {fn}: shape={pm.shape}, cols={list(pm.columns)}")
    if "batchChannel" in pm.columns:
        print(f"    batchChannel sample: {list(pm['batchChannel'].head(3))}")
        overlap_bc = pm["batchChannel"].isin(prot_raw.columns).sum()
        print(f"    batchChannel in prot_raw cols: {overlap_bc}/{len(pm)}")
    if "specimenID" in pm.columns:
        print(f"    specimenID sample: {list(pm['specimenID'].head(3))}")
    if "isAssayControl" in pm.columns:
        print(f"    isAssayControl values: {pm['isAssayControl'].unique()[:5]}")
