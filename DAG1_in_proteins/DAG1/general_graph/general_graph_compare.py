"""Compare the data-derived B to the literature matrices C (confidence) and Sign
(effect direction). B is NOT built from C; this is an after-the-fact agreement check.

  python DAG1/general_graph/general_graph_compare.py      (reads B from this folder, C = ../../CHRONOS_80x80_confidence.xlsx)
"""
import numpy as np, pandas as pd
from pathlib import Path
OUT=Path(__file__).resolve().parent                       # DAG1/general_graph
PRIOR=OUT.parent.parent/"CHRONOS_80x80_confidence.xlsx"   # top-level folder
B=pd.read_csv(OUT/"B_causal_graph.csv",index_col=0)        # B[target,source] = source->target
P=list(B.index); idx={p:i for i,p in enumerate(P)}
Cm=pd.read_excel(PRIOR,sheet_name="Matrix_80x80",index_col=0).reindex(index=P,columns=P)
Sg=pd.read_excel(PRIOR,sheet_name="Sign_80x80",index_col=0).reindex(index=P,columns=P)
Bv=B.values; Cv=Cm.values.astype(float); Sv=Sg.values.astype(float)
np.fill_diagonal(Cv,0); np.fill_diagonal(Sv,0)

# ---- directed edge table: s->t  (C row=source col=target; B[t,s]) ----
rows=[]
for s in range(len(P)):
    for t in range(len(P)):
        if s==t: continue
        c=Cv[s,t]; sg=Sv[s,t]; b=Bv[t,s]
        if c==0 and abs(b)<1e-9: continue
        rows.append((P[s],P[t],int(c),int(sg),round(float(b),4)))
E=pd.DataFrame(rows,columns=["source","target","C_conf","C_sign","B_data"])
E["absB"]=E.B_data.abs()
E.to_csv(OUT/"edge_table_B_vs_C.csv",index=False)

print("=== 1. Does |B| grow with literature confidence?  (existence check) ===")
g=E.groupby("C_conf").agg(n=("B_data","size"),frac_nonzero=("B_data",lambda x:(x.abs()>1e-9).mean()),
                          mean_absB=("absB","mean")).round(3)
print(g.to_string())

print("\n=== 2. Sign agreement: for confident edges (C>=+3) with a defined effect sign, does sign(B) match? ===")
for thr in [5,4,3]:
    sub=E[(E.C_conf>=thr)&(E.C_sign!=0)&(E.B_data.abs()>1e-9)]
    if len(sub):
        agree=(np.sign(sub.B_data)==np.sign(sub.C_sign)).mean()
        print(f"  C>=+{thr}: {len(sub)} edges w/ sign & nonzero B, sign-match={agree:.2f}")

print("\n=== 3. Negative-confidence edges (C<0): are they ~absent in B? ===")
for c in [-1,-2]:
    sub=E[E.C_conf==c]
    if len(sub): print(f"  C={c}: {len(sub)} edges, mean|B|={sub.absB.mean():.3f}, frac|B|<0.02={ (sub.absB<0.02).mean():.2f}")
pos=E[E.C_conf>=3]; print(f"  (for contrast) C>=+3: mean|B|={pos.absB.mean():.3f}")

print("\n=== 4. Directional asymmetry in B (is direction even identified?) ===")
Boff=Bv.copy(); np.fill_diagonal(Boff,0)
mask=(Boff!=0)|(Boff.T!=0)
sym=np.corrcoef(Boff[mask], Boff.T[mask])[0,1]
print(f"  corr(B, B^T) over edges = {sym:.3f}   (near 1 => direction NOT identified by node-wise regression)")

print("\n=== 5. Spotlight edges ===")
def show(s,t):
    r=E[(E.source==s)&(E.target==t)]
    rr=E[(E.source==t)&(E.target==s)]
    b1=r.B_data.values[0] if len(r) else 0; c1=r.C_conf.values[0] if len(r) else 0
    b2=rr.B_data.values[0] if len(rr) else 0; c2=rr.C_conf.values[0] if len(rr) else 0
    print(f"  {s}->{t}: B={b1:+.3f} (C={c1:+d})    {t}->{s}: B={b2:+.3f} (C={c2:+d})")
for a,b in [("APP","BACE1"),("TMED2","TMED10"),("TMED2","TMED9"),("TMED10","APP"),
            ("PSEN1","NCSTN"),("COPA","COPB1"),("APOE","TMED2")]:
    show(a,b)
print("\nwrote edge_table_B_vs_C.csv")
