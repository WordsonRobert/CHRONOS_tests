"""
extra_stats.py — per-dataset statistics that need the data itself (not just results/ files):
signed/unsigned AUROC raw and adjusted, correlation summaries for prior vs non-prior pairs,
how each prior pair's arrow was set, missingness, and the corrected region null (HIW 9.1).
Reads Z* from results/{ds}_summary.json, so it works for any prior. Read-only except for
{results}/extra_stats.json and {results}/{ds}_regions_usable_null_check.csv.

python extra_stats.py --results results    --prior ../CHRONOS_80x80_confidence.xlsx    --data DATA_DIR
python extra_stats.py --results results_v3 --prior CHRONOS_80x80_confidence_v3.xlsx   --data DATA_DIR
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from pathlib import Path
sys.argv, _argv = [sys.argv[0]], sys.argv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chronos_regions as cr
sys.argv = _argv
ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True); ap.add_argument("--prior", required=True); ap.add_argument("--data", required=True)
ap.add_argument("--banner_tmt", default="data/Banner_TMT_log2ratio_GIS.csv"); ap.add_argument("--banner_lfq", default="data/Banner_LFQ_log2.csv")
ap.add_argument("--rosmap_ids", default="results/ROSMAP_individuals.txt")
ap.add_argument("--banner_exclude", default="banner_exclude_possible_diverse_overlap.txt")
ap.add_argument("--skip_regions", action="store_true")
a = ap.parse_args(); R = Path(a.results)
CFG = {"ROSMAP": (None, None), "Diverse": (None, a.rosmap_ids),
       "Banner": (os.path.abspath(a.banner_tmt), a.banner_exclude), "BannerLFQ": (os.path.abspath(a.banner_lfq), a.banner_exclude)}

def prep(ds, prot, excl, Zstar):
    Rc = dict(cr.RECIPES[ds]); D = Path(a.data)
    if prot: Rc["prot"] = prot
    proteins, C = cr.load_prior(a.prior)
    X, found, rest = cr.load_proteomics(D / Rc["prot"], Rc["fmt"], proteins)
    s2i, batch = cr.sample_to_individual(list(X.index), D / Rc["assay"], D / Rc["biospec"])
    if Rc.get("batch_from_channel"): batch = {s: s.split(".")[0] for s in s2i}
    Xi, b = cr.collapse_to_individuals(X, s2i, batch); rest_i, _ = cr.collapse_to_individuals(rest, s2i, batch)
    if excl:
        ex = set(Path(excl).read_text().split()); Xi = Xi.loc[~Xi.index.isin(ex)]; rest_i = rest_i.loc[Xi.index]
    clin = pd.read_csv(D / Rc["clin"], dtype=str); idc = next(c for c in clin.columns if c.lower() == "individualid")
    clin[idc] = clin[idc].str.strip(); clin = clin.drop_duplicates(idc).set_index(idc).reindex(Xi.index); clin["batch"] = b.reindex(Xi.index).values
    cov, labels = cr.clean_covariates(clin, proteins); pcov, _ = cr.proteome_covariates(rest_i.reindex(Xi.index), lambda s: None)
    cov = pd.concat([cov, pcov.reindex(Xi.index)], axis=1)
    Xa = Xi[found].values.astype(float); fi = [proteins.index(p) for p in found]; Cs = C[np.ix_(fi, fi)]
    Xf = cr.residualise(Xa, cr.design(cov, Zstar, labels)) if Zstar else Xa
    return found, Cs, Xa, Xf, cov, labels

def unsigned_auroc(rho, N, em, iu):
    r = np.abs(rho[iu]); e = em[iu]; ok = N[iu] >= cr.MIN_PAIR_N; r, e = r[ok], e[ok]
    rk = cr.rankdata(r); n1, n0 = e.sum(), (~e).sum(); return float((rk[e].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

out = {}
for ds, (prot, excl) in CFG.items():
    if not (R / f"{ds}_summary.json").exists(): continue
    Z = json.loads((R / f"{ds}_summary.json").read_text())["Z_star"]
    found, Cs, Xa, Xf, cov, labels = prep(ds, prot, excl, Z)
    iu = np.triu_indices(len(found), 1); em = np.maximum(Cs, Cs.T) >= cr.MIN_PRIOR
    usable = np.isfinite(Xf).any(1)
    st = dict(people=int(len(Xa)), usable_after_adjustment=int(usable.sum()), nodes=len(found),
              missing_fraction=round(float(np.isnan(Xa).mean()), 4))
    for lab, M in [("raw", Xa), ("adjusted", Xf)]:
        rho, N = cr.spearman_pairwise(M); r = rho[iu]; n = N[iu]; ok = n >= cr.MIN_PAIR_N; e = em[iu]
        st[lab] = dict(pairs_scored=int(ok.sum()), prior_pairs_scored=int((e & ok).sum()),
                       signed_auroc=round(float(cr.auroc_vs_prior(rho, N, em, iu)), 4), unsigned_auroc=round(unsigned_auroc(rho, N, em, iu), 4),
                       prior_mean_rho=round(float(r[e & ok].mean()), 3), prior_frac_pos=round(float(np.mean(r[e & ok] > 0)), 3), prior_mean_abs=round(float(np.abs(r[e & ok]).mean()), 3),
                       other_mean_rho=round(float(r[~e & ok].mean()), 3), other_frac_pos=round(float(np.mean(r[~e & ok] > 0)), 3), other_mean_abs=round(float(np.abs(r[~e & ok]).mean()), 3))
    rho, N = cr.spearman_pairwise(Xf); m = len(found)
    cnt = dict(one_way_agrees=0, one_way_disagrees_downweighted=0, two_way_picked=0, low_n_placeholder=0)
    for i in range(m):
        for j in range(i + 1, m):
            if Cs[i, j] < cr.MIN_PRIOR and Cs[j, i] < cr.MIN_PRIOR: continue
            if N[i, j] < cr.MIN_PAIR_N: cnt["low_n_placeholder"] += 1; continue
            ok = np.isfinite(Xf[:, i]) & np.isfinite(Xf[:, j]); fwd = cr.lingam_direction(Xf[ok, i], Xf[ok, j]) > 0
            if Cs[i, j] >= cr.MIN_PRIOR and Cs[j, i] >= cr.MIN_PRIOR: cnt["two_way_picked"] += 1
            elif (fwd and Cs[i, j] >= cr.MIN_PRIOR) or ((not fwd) and Cs[j, i] >= cr.MIN_PRIOR): cnt["one_way_agrees"] += 1
            else: cnt["one_way_disagrees_downweighted"] += 1
    st["arrow_setting"] = cnt
    if not a.skip_regions:
        reg = pd.read_csv(R / f"{ds}_regions.csv"); masks = dict(cr.candidate_regions(cov, labels))
        rng = np.random.default_rng(5); cache = {}; rows = []; U = np.where(usable)[0]
        for _, r in reg.iterrows():
            parts = r.region.split("  AND  ")
            if not all(p in masks for p in parts): continue
            mk = np.logical_and.reduce([masks[p] for p in parts]) & usable; k = int(mk.sum())
            s = cr.auroc_vs_prior(*cr.spearman_pairwise(Xf[mk]), em, iu) if k >= cr.MIN_REGION_N else np.nan
            key = max(cr.MIN_REGION_N, int(round(k / 10)) * 10)
            if key not in cache:
                cache[key] = np.array([cr.auroc_vs_prior(*cr.spearman_pairwise(Xf[rng.choice(U, min(key, len(U)), replace=False)]), em, iu) for _ in range(cr.N_NULL_REGION)])
            nl = cache[key]; p = (1 + np.sum(nl >= s)) / (1 + len(nl)) if np.isfinite(s) else np.nan
            rows.append(dict(region=r.region, depth=r.depth, n_listed=r.n, n_usable=k, auroc=s, null_usable=nl.mean(),
                             z_fixed=(s - nl.mean()) / (nl.std() + 1e-9), p_fixed=p, z_orig=r.z, q_orig=r.q_fdr))
        o = pd.DataFrame(rows); okp = o.p_fixed.notna(); o.loc[okp, "q_fixed"] = cr.bh(o.loc[okp, "p_fixed"].values)
        o["strong_fixed"] = (o.q_fixed < 0.10) & (o.z_fixed > 2)
        o.round(5).to_csv(R / f"{ds}_regions_usable_null_check.csv", index=False)
        st["regions_corrected"] = dict(rescored=int(len(o)), strong=int(o.strong_fixed.sum()), z_above_2=int((o.z_fixed > 2).sum()), min_q=round(float(o.q_fixed.min()), 4))
    out[ds] = st; print(ds, json.dumps(st))
(R / "extra_stats.json").write_text(json.dumps(out, indent=2))
