"""
prep_banner.py — turn Banner raw search outputs into a clean protein x sample matrix
that chronos_regions.py reads with fmt=gene_pipe (rows GENE|Accession, cols = sample IDs).

TMT  (bannertmt_22batchmulticoncensus_Proteins.txt, Proteome Discoverer):
   value = log2( sample abundance / GIS-pool 126 abundance in the same batch )
   -> the pool channel cancels batch effects (same idea as ROSMAP's processed file)
   -> per-sample median centring (equal loading)
   columns renamed F{k}: {chan}  ->  b{k:02d}.{chan}  (= batchChannel in the assay metadata)
LFQ  (banner_proteomics_pfc_proteinoutput.txt, MaxQuant):
   value = log2(LFQ intensity), 0 -> missing, per-sample median centring
   columns renamed LFQ.intensity.b1_002_07 -> b1_002_07_lfq (= specimenID)
Isoforms: when a gene has several rows, keep the row with the most observed samples.
"""
import re, sys, numpy as np, pandas as pd
from pathlib import Path

def keep_best(df, genes, acc):
    df = df.copy(); df["_g"] = genes.values; df["_a"] = acc.values
    df["_n"] = df.drop(columns=["_g", "_a"]).notna().sum(axis=1)
    df = df[df["_g"].notna() & (df["_g"] != "nan")].sort_values("_n", ascending=False).drop_duplicates("_g")
    df.index = df["_g"] + "|" + df["_a"]
    return df.drop(columns=["_g", "_a", "_n"])

def center(m):
    return m - m.median(axis=0)

def tmt(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    ab = {c: re.match(r"Abundances \(Normalized\): F(\d+): (\w+),", c) for c in d.columns}
    ab = {c: m for c, m in ab.items() if m}
    out = {}
    for k in sorted({int(m.group(1)) for m in ab.values()}):
        cols = {m.group(2): c for c, m in ab.items() if int(m.group(1)) == k}
        ref = d[cols["126"]].where(d[cols["126"]] > 0)
        for ch, c in cols.items():
            if ch == "126":
                continue
            out[f"b{k:02d}.{ch}"] = np.log2(d[c].where(d[c] > 0) / ref)
    m = center(pd.DataFrame(out))
    genes = d["Gene Symbol"].astype(str).str.split(";").str[0].str.strip()
    return keep_best(m, genes, d["Accession"].astype(str))

def lfq(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    lc = [c for c in d.columns if c.startswith("LFQ.intensity.")]
    m = np.log2(d[lc].where(d[lc] > 0))
    m.columns = [c.replace("LFQ.intensity.", "") + "_lfq" for c in lc]
    m = center(m)
    genes = d["GN"].astype(str).str.split(";").str[0].str.strip()
    acc = d["Unique.ID"].astype(str).str.split("|").str[-1]
    return keep_best(m, genes, acc)

if __name__ == "__main__":
    D = Path(sys.argv[1]); out = Path(sys.argv[2]) if len(sys.argv) > 2 else D
    t = tmt(D / "bannertmt_22batchmulticoncensus_Proteins.txt")
    t.to_csv(out / "Banner_TMT_log2ratio_GIS.csv"); print("TMT", t.shape, "missing", round(t.isna().mean().mean(), 3))
    l = lfq(D / "banner_proteomics_pfc_proteinoutput.txt")
    l.to_csv(out / "Banner_LFQ_log2.csv"); print("LFQ", l.shape, "missing", round(l.isna().mean().mean(), 3))
