"""
compare_three.py — compare the per-dataset CHRONOS graphs made by chronos_regions.py

Compares DATA-SUPPORTED edges (prior edge + correlation significant at FDR<0.05),
not all prior edges: the prior alone fixes which 507 edges get a non-zero W, so
comparing W>0 sets is circular (that's why the old comparison showed
precision = 1.000 in every dataset).

Usage:  python compare_three.py --results results
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import hypergeom

ap = argparse.ArgumentParser(); ap.add_argument("--results", default="results"); a = ap.parse_args()
R = Path(a.results)
ds = [d for d in ["ROSMAP", "Diverse", "Banner"] if (R / f"{d}_causal_edges.csv").exists()]
E = {d: pd.read_csv(R / f"{d}_causal_edges.csv") for d in ds}
S = {d: json.loads((R / f"{d}_summary.json").read_text()) for d in ds}
lines = []
def log(s=""): print(s); lines.append(s)

log("=" * 72); log("CHRONOS — cross-dataset comparison  (" + ", ".join(ds) + ")"); log("=" * 72)

log("\n1. HOW WELL EACH DATASET RECOVERS THE CONFIDENCE GRAPH  (AUROC; 0.5 = chance)")
log(f"   {'dataset':<9}{'people':>7}{'raw':>8}{'adjusted':>10}   trajectory variables Z*")
for d in ds:
    s = S[d]
    log(f"   {d:<9}{s['n_people']:>7}{s['auroc_raw']:>8.3f}{s['auroc_adjusted']:>10.3f}   {s['Z_star']}")

log("\n2. REGIONS WHERE THE CONFIDENCE GRAPH IS STRONGEST (after adjusting for Z*)")
for d in ds:
    s = S[d]
    log(f"   {d}: {s['n_strong_regions']} regions beat random same-size subsets (FDR<0.10, z>2)")
    for r in s["top_regions"][:3]:
        log(f"      AUROC {r['auroc']:.3f}  z={r['z']:.2f}  q={r['q_fdr']:.3f}  n={r['n']}  {r['region']}")

key = lambda df: set(zip(df.source, df.target))
sup = {d: E[d][E[d].data_supported == 1] for d in ds}
und = lambda st: {tuple(sorted(e)) for e in st}
log("\n3. DATA-SUPPORTED EDGES  (prior edge AND correlation FDR < 0.05)")
for d in ds:
    log(f"   {d:<9} {len(sup[d]):>4} / {len(E[d])} prior edges "
        f"({100 * len(sup[d]) / len(E[d]):.0f}%)")

if len(ds) >= 2:
    log("\n4. AGREEMENT BETWEEN DATASETS")
    log("   'expected' = overlap if each dataset's supported edges were a random pick of")
    log("   the prior edges both datasets could test (hypergeometric); p = chance of >= shared")
    log(f"   {'pair':<20}{'Jaccard':>9}{'shared':>8}{'expected':>10}{'fold':>6}{'p':>10}{'same direction':>16}")
    rows = []
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            a_, b_ = und(key(sup[ds[i]])), und(key(sup[ds[j]]))
            sh = a_ & b_
            testable = lambda d: und(key(E[d][E[d].n_pairwise >= 20]))
            U = testable(ds[i]) & testable(ds[j]); a2, b2 = a_ & U, b_ & U
            exp = len(a2) * len(b2) / max(1, len(U))
            pv = hypergeom.sf(len(sh) - 1, len(U), len(a2), len(b2))
            jac = len(sh) / max(1, len(a_ | b_))
            # direction: which way W points in each dataset for shared pairs
            def dirn(d, u, v):
                w = E[d].set_index(["source", "target"])["W"]
                return np.sign(w.get((u, v), 0) - w.get((v, u), 0))
            agree = [dirn(ds[i], u, v) == dirn(ds[j], u, v) for u, v in sh]
            da = np.mean(agree) if agree else np.nan
            log(f"   {ds[i] + ' vs ' + ds[j]:<20}{jac:>9.3f}{len(sh):>8}{exp:>10.1f}{len(sh)/max(exp,1e-9):>6.2f}{pv:>10.1e}{100 * da:>15.1f}%")
            rows.append(dict(pair=f"{ds[i]}-{ds[j]}", jaccard=jac, shared=len(sh), expected_by_chance=exp,
                             fold=len(sh)/max(exp,1e-9), p_hypergeom=pv, direction_agree=da))
    pd.DataFrame(rows).to_csv(R / "comparison_pairs.csv", index=False)

    log("\n5. CONSENSUS PAIRS  (protein pair data-supported in >= 2 datasets)")
    log("   Counted on pairs, not arrows: direction is mostly set by the prior, and where the")
    log("   prior is two-way the data's direction choice is unstable (see 'same direction').")
    allk = {}
    for d in ds:
        for _, r in sup[d].iterrows():
            allk.setdefault(tuple(sorted((r.source, r.target))), {})[d] = r
    cons = []
    for (u, v), dd in allk.items():
        if len(dd) >= 2:
            row = dict(protein_a=u, protein_b=v, n_datasets=len(dd),
                       prior_C_a_to_b=float(E[ds[0]].set_index(["source", "target"]).prior_C.get((u, v), 0)),
                       prior_C_b_to_a=float(E[ds[0]].set_index(["source", "target"]).prior_C.get((v, u), 0)))
            for d in ds:
                row[f"dir_{d}"] = f"{dd[d].source}->{dd[d].target}" if d in dd else ""
                row[f"W_{d}"] = dd[d].W if d in dd else 0.0
                row[f"beta_{d}"] = dd[d].beta if d in dd else 0.0
            dirs = [row[f"dir_{d}"] for d in ds if row[f"dir_{d}"]]
            row["direction_consistent"] = int(len(set(dirs)) == 1)
            cons.append(row)
    cdf = pd.DataFrame(cons)
    if len(cdf):
        cdf["mean_W"] = cdf[[f"W_{d}" for d in ds]].mean(axis=1)
        cdf = cdf.sort_values(["n_datasets", "mean_W"], ascending=False)
        cdf.to_csv(R / "consensus_edges.csv", index=False)
        log(f"   {len(cdf)} consensus pairs ({(cdf.n_datasets == len(ds)).sum()} in all {len(ds)}); "
            f"direction consistent in {cdf.direction_consistent.sum()}; top 20:")
        for _, r in cdf.head(20).iterrows():
            log(f"   {r.protein_a:>8} - {r.protein_b:<8} " +
                "  ".join(f"{d}:{r[f'dir_{d}'] or '-':<16} W={r[f'W_{d}']:.3f}" for d in ds))
    else:
        log("   none")

if (R / "BannerLFQ_causal_edges.csv").exists() and "Banner" in ds:
    L = pd.read_csv(R / "BannerLFQ_causal_edges.csv")
    T = E["Banner"]
    m = T.merge(L, on=["source", "target"], suffixes=("_tmt", "_lfq"))
    m = m[(m.n_pairwise_tmt >= 20) & (m.n_pairwise_lfq >= 20)]
    a_, b_ = und(key(m[m.data_supported_tmt == 1])), und(key(m[m.data_supported_lfq == 1]))
    U = und(key(m)); sh = a_ & b_
    exp = len(a_) * len(b_) / max(1, len(U)); pv = hypergeom.sf(len(sh) - 1, len(U), len(a_), len(b_))
    log("\n6. TECHNICAL CHECK: Banner TMT vs Banner LFQ (same people, different mass-spec method)")
    log(f"   edges testable in both: {len(U)} pairs | supported TMT {len(a_)}, LFQ {len(b_)}, both {len(sh)} "
        f"(chance {exp:.1f}, p={pv:.1e})")
    log(f"   Spearman of |W| across shared testable edges: "
        f"{pd.Series(m.W_tmt.abs()).corr(pd.Series(m.W_lfq.abs()), method='spearman'):.3f}")

(R / "comparison_report.txt").write_text("\n".join(lines))
