"""
Adaptive trajectory-region discovery (23 Sept methodology).

For each ordered trajectory variable Z_k and each dataset:
  1. rank-normalise Z_k -> u in [0,1]  (subjects sorted along the axis)
  2. slide a window of `w` subjects; inside each window fit a LOCAL slope
        beta_i(u) = slope of protein i's (normalised) level vs Z_k
     for every protein i.
  3. for a confidence-c edge (i,j): local behaviour B_ij(u) = (sgn beta_i, sgn beta_j)
     first-pass rule: BOTH slopes same (non-zero) sign  ==  co-moving == satisfies rule.
  4. aggregate edges by literature score c in {-2,-1,+1,+2,+3,+4,+5}:
        p_c(u) = fraction of confidence-c edges satisfying the rule at u
        (denominator = edges whose BOTH slopes are defined in that window)
  5. discovered set  S_c = { u : p_c(u) > tau }  ->  unions of intervals.

Output: one master table of discovered regions + the full p_c(u) curves per
(dataset, variable) for plotting.  No value invented; edges/scores read straight
from CHRONOS_80x80_confidence.xlsx (Matrix_80x80), protein levels from the
per-variable CSVs (expr_z = within-dataset per-protein z-score).

Rule is deliberately the simple co-movement one (sign-agnostic). Making it
sign-aware (using the Sign_80x80 sheet, so inhibitory edges expect OPPOSITE
directions) is a one-line change flagged in the README.
"""
import sys, warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from pathlib import Path

# paths relative to this file (works whatever this folder is renamed to):
#   <top folder>/<this folder>/trajectory_regions.py ; data in <top folder>/DATA
OUT     = Path(__file__).resolve().parent                              # outputs written next to this script
CSV_DIR = OUT.parent / "DATA"                                          # the 14 per-variable CSVs
PRIORXL = OUT.parent / "CHRONOS_80x80_confidence.xlsx"             # literature matrix C (v3 contents)
(OUT/"curves").mkdir(parents=True, exist_ok=True)

# ---- parameters -----------------------------------------------------------
W        = 40      # window size in subjects
STEP     = 2       # slide step (evaluate every STEP-th window centre)
# per-score adaptive threshold tau_c (per dataset,variable,score):
#   tau_c = clip( mean_c + Z_ADAPT * std_c , floor=FLOOR , cap=0.97 )
# i.e. "where does confidence-c co-movement rise clearly above its own baseline".
Z_ADAPT  = 1.25
FLOOR    = 0.62
FLAT_R   = 0.15    # |local Pearson r| below this = "no change" (flat) along the axis.
                   #   positive edges: BOTH proteins must have |r|>=FLAT_R and same sign (co-move)
                   #   negative edges: BOTH proteins must have |r|<FLAT_R (inert / didn't change)
MIN_RUN  = 3       # a region must stay above tau_c for >=MIN_RUN consecutive windows
MIN_VALID= 15      # min non-missing subjects to fit a protein's local slope
MIN_PEOPLE_REGION = 15   # min subjects a discovered interval must span
SCORES   = [-2,-1,1,2,3,4,5]

# variables with a usable ordered axis (the rolling slope is meaningless for
# an unordered category or a 0/1 flag, so those are reported as skipped)
ORDERED = ["cts_mmse30_lv","cts_mmse30_first_ad_dx","pmi","braaksc","ceradsc",
           "amyThal","amyA","amyCerad","reag","PlaqueTotal","TangleTotal"]
SKIPPED = {"sex":"binary / unordered", "apoe_genotype":"unordered category",
           "amyAny":"binary 0/1 flag"}

# ---- confidence edges -----------------------------------------------------
Mc = pd.read_excel(PRIORXL, sheet_name="Matrix_80x80", index_col=0)
PROT = list(Mc.index)
idx  = {p:i for i,p in enumerate(PROT)}
Cmat = Mc.values.astype(float); np.fill_diagonal(Cmat, np.nan)
EDGES = {c: np.argwhere(Cmat == c) for c in SCORES}   # directed (i,j) cells

def local_corr(Zw, Xw):
    """vectorised local Pearson r of every protein (cols of Xw) vs Zw, within a
    window. NaN where <MIN_VALID valid subjects or no variance.
    r carries BOTH the trend direction (sign) and the strength (|r|), so one
    statistic serves co-movement (positive scores) and flatness (negative scores)."""
    valid = ~np.isnan(Xw)                      # (w, P)
    n = valid.sum(0)
    z = Zw[:, None]
    zbar = np.where(n>0, (valid*z).sum(0)/np.maximum(n,1), np.nan)
    xbar = np.where(n>0, np.nansum(Xw,0)/np.maximum(n,1), np.nan)
    dz = np.where(valid, z - zbar, 0.0)
    dx = np.where(valid, np.nan_to_num(Xw) - xbar, 0.0)
    num  = (dz*dx).sum(0)
    denz = (dz*dz).sum(0)
    denx = (dx*dx).sum(0)
    r = np.where((n>=MIN_VALID) & (denz>0) & (denx>0),
                 num/np.sqrt(np.where((denz*denx)==0,1,denz*denx)), np.nan)
    return r

def pc_curve(values, X):
    """values: (n,) sorted axis; X: (n,P) protein z-scores sorted the same way.
    returns u_centres and p_c for each score."""
    n = len(values)
    centres = list(range(0, n-W+1, STEP))
    u_all = (np.arange(n))/(n-1)
    out = {c: [] for c in SCORES}
    uc, valc = [], []
    for t in centres:
        r = local_corr(values[t:t+W], X[t:t+W])   # (P,) Pearson r vs axis, or nan
        s = np.sign(r); flat = np.abs(r) < FLAT_R  # flat = no local trend
        uc.append(u_all[t:t+W].mean())
        valc.append(np.nanmedian(values[t:t+W]))
        for c in SCORES:
            e = EDGES[c]
            if len(e)==0: out[c].append(np.nan); continue
            ri, rj = r[e[:,0]], r[e[:,1]]
            defined = ~np.isnan(ri) & ~np.isnan(rj)
            if defined.sum()==0: out[c].append(np.nan); continue
            if c > 0:
                # expected (positive edge): the two proteins move the SAME direction
                sat = defined & (np.sign(ri)==np.sign(rj)) & (ri!=0) & (rj!=0)
            else:
                # expected (negative edge): pair is inert -> BOTH flat / no change
                sat = defined & (np.abs(ri)<FLAT_R) & (np.abs(rj)<FLAT_R)
            out[c].append(sat.sum()/defined.sum())
    return np.array(uc), np.array(valc), {c:np.array(out[c]) for c in SCORES}

def intervals_above(u, val, p, tau):
    """contiguous runs where p>tau -> list of (u_lo,u_hi,val_lo,val_hi,support,npeople)."""
    mask = p > tau
    regs=[]; i=0
    while i < len(mask):
        if mask[i] and not np.isnan(p[i]):
            j=i
            while j+1<len(mask) and mask[j+1] and not np.isnan(p[j+1]): j+=1
            if (j-i+1) < MIN_RUN:      # not a sustained region, skip
                i=j+1; continue
            ulo,uhi = u[i],u[j]
            vlo,vhi = val[i],val[j]
            # region spans from the first window's start to the last window's end,
            # so credit the full window width W (centres are W/2 in from the edges)
            npeople = int(round((uhi-ulo)*(NglobalN-1))) + W
            support = float(np.nanmean(p[i:j+1]))
            regs.append((ulo,uhi,vlo,vhi,support,npeople))
            i=j+1
        else:
            i+=1
    return regs

rows=[]; NglobalN=0
for var in ORDERED:
    f = CSV_DIR/f"{var}.csv"
    if not f.exists(): continue
    d = pd.read_csv(f)
    for ds in d["dataset"].unique():
        sub = d[d.dataset==ds]
        piv = sub.pivot_table(index="individualID", columns="protein", values="expr_z")
        vals = sub.groupby("individualID")[var].first().reindex(piv.index)
        piv = piv.reindex(columns=PROT)                     # align to 80-protein order
        order = np.argsort(vals.values, kind="mergesort")   # sort subjects by axis
        v = vals.values[order].astype(float)
        X = piv.values[order]
        if len(v) < W+STEP:            # too few subjects to slide a window
            continue
        NglobalN = len(v)
        u, valc, pc = pc_curve(v, X)
        # save curves
        cur = pd.DataFrame({"u":u, "value_median":valc})
        for c in SCORES: cur[f"p_{c:+d}"] = pc[c]
        cur.to_csv(OUT/"curves"/f"{ds}_{var}.csv", index=False)
        # discovered regions per score, adaptive tau_c
        for c in SCORES:
            pcv = pc[c]; good = pcv[~np.isnan(pcv)]
            if len(good) < 3: continue
            tau_c = float(np.clip(good.mean() + Z_ADAPT*good.std(), FLOOR, 0.97))
            for (ulo,uhi,vlo,vhi,support,npeople) in intervals_above(u,valc,pcv,tau_c):
                if npeople < MIN_PEOPLE_REGION: continue
                rows.append(dict(dataset=ds, variable=var, confidence=c,
                                 u_lo=round(ulo,3), u_hi=round(uhi,3),
                                 value_lo=round(float(vlo),3), value_hi=round(float(vhi),3),
                                 support_fraction=round(support,3), tau_c=round(tau_c,3),
                                 n_people=npeople,
                                 n_edges_c=int(len(EDGES[c])), n_subjects=len(v)))
        # also record the overall mean p_c (baseline context)
COLS=["dataset","variable","confidence","u_lo","u_hi","value_lo","value_hi",
      "support_fraction","tau_c","n_people","n_edges_c","n_subjects"]
tab = pd.DataFrame(rows, columns=COLS)
if len(tab): tab = tab.sort_values(["dataset","variable","confidence","u_lo"])
tab.to_csv(OUT/"discovered_regions.csv", index=False)

# baseline table: mean p_c(u) per (dataset,variable,score) for context
base=[]
for cf in sorted((OUT/"curves").glob("*.csv")):
    name=cf.stem; d=pd.read_csv(cf)
    ds,var = name.split("_",1)
    for c in SCORES:
        col=f"p_{c:+d}"
        base.append(dict(dataset=ds,variable=var,confidence=c,
                         mean_pc=round(float(np.nanmean(d[col])),3),
                         max_pc=round(float(np.nanmax(d[col])),3)))
pd.DataFrame(base).to_csv(OUT/"baseline_pc.csv", index=False)

print(f"discovered {len(tab)} regions across {tab[['dataset','variable']].drop_duplicates().shape[0]} (dataset,variable) combos")
print("scores present per confidence:", {c:int(len(EDGES[c])) for c in SCORES})
print("\nregions per confidence level:")
print(tab.groupby('confidence').size().to_string() if len(tab) else "  (none above tau)")
print("\nsample of discovered regions:")
print(tab.head(20).to_string(index=False) if len(tab) else "none")
json.dump(dict(W=W,STEP=STEP,Z_ADAPT=Z_ADAPT,FLOOR=FLOOR,FLAT_R=FLAT_R,MIN_VALID=MIN_VALID,
               MIN_RUN=MIN_RUN,MIN_PEOPLE_REGION=MIN_PEOPLE_REGION,skipped=SKIPPED,
               rule="positive c: both |r|>=FLAT_R & same sign (co-move); "
                    "negative c: both |r|<FLAT_R (flat/no change)"),
          open(OUT/"_params.json","w"), indent=2)
