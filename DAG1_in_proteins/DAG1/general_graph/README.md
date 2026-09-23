# The general 80×80 causal graph  B

`X = B X + Γ Z + ε`,  diag(B)=0. One graph over all 80 proteins, estimated across the whole
pooled population, with the 14 trajectory variables as context. **C (literature) is NOT used to
build B** — only compared afterward.

## How it was built
- **One aligned table, one row per subject** (1,570 people: ROSMAP 400 + Diverse 980 + Banner 190).
  Subjects are pooled, NOT stacked 14×. `aligned_XZ.csv` = the table (80 X columns = expr_z,
  within-dataset z-score; 14 Z columns; dataset label).
- **Node-wise elastic-net regression**: for each protein *i*, regress X_i on the other 79 proteins
  + 14 Z (standardised) + 2 dataset dummies. Coefficients on proteins → row *i* of **B**; on Z →
  row *i* of **Γ**. Elastic net (l1∈{.5,.8}, 3-fold CV) handles the heavy collinearity among
  complex members. Missing protein predictors filled at 0 (=mean); missing Z at 0 after standardising.
- **Convention**: `B_causal_graph.csv` is `B[target_row, source_col]` — i.e. cell (i,j) = effect of
  protein j on protein i = **j → i**. Weighted, continuous, asymmetry allowed, cycles allowed.

## Run
```
python DAG1/general_graph/general_graph.py           # refits B from aligned_XZ.csv (~5 s)
python DAG1/general_graph/general_graph_compare.py   # B vs literature C (v3) -> edge_table_B_vs_C.csv
```
Both run from anywhere (paths are relative to the scripts) and reproduce the committed files.
`aligned_XZ.csv` was assembled from the raw AMP-AD Synapse files and carries the same subjects and values as
`DATA/`; `apoe_genotype` in it is already the e4 dosage (0/1/2) and `sex` is 1 = male, 0 = female.

## Files
- `B_causal_graph.csv` — the 80×80 graph (rows=target/effect, cols=source/cause).
- `Gamma_Z_to_X.csv` — 80×14 context effects (nuisance: how each Z drives each protein).
- `aligned_XZ.csv` — the per-subject (X,Z) table B was fit from.
- `node_R2.csv` — variance of each protein explained (median 0.54).
- `edge_table_B_vs_C.csv` — every directed edge: source, target, C_conf, C_sign, B_data, |B|.
- `_fit_summary.json`.

## Fit
1,570 subjects; **3,417 nonzero directed edges (54%)**; median node R² **0.54**.

## Agreement with the literature (after the fact)
1. **|B| grows with confidence only at the top.** C=+5 edges: 73% nonzero, mean|B|=0.082 — far above
   every other tier (all ~0.02–0.04). So the data strongly corroborates the **+5 structural-complex**
   edges, but the mid tiers (+1…+4) are **not** clearly separated from noise by magnitude.
2. **Sign agreement** (confident edges with a defined effect sign): +4 → 71%, +3 → 65%, +5 → 58%.
   Better than chance, not stellar (collinear complex members can flip partial-coefficient sign).
3. **Negative-confidence edges are weak, as they should be**: C=−1 mean|B|=0.019 (68% below 0.02),
   C=−2 mean|B|=0.025 — vs C≥+3 mean|B|=0.045. The "shouldn't exist" edges are mostly near-zero.
4. **Spotlights**: TMED2–TMED10 strong (0.40/0.54, the core heterodimer holds); PSEN1–NCSTN
   (0.14/0.18, γ-secretase holds); COPI members modest; **APP–BACE1 ≈ 0 in both directions**
   despite C=+3/+5 — the canonical β-secretase edge washes out once you condition on all 78 other
   proteins + Z (mediated / not co-varying at the abundance level); TMED10–APP weak (~−0.03).

## Caveats — important
- **Direction is essentially NOT identified.** corr(B, Bᵀ) = 0.876: B comes out nearly symmetric,
  so `B_ij` ≈ `B_ji` for almost every pair (see the TMED/PSEN1 spotlights). This is expected —
  node-wise Gaussian-type regression recovers *partial associations*, which are symmetric; it cannot
  tell A→B from B→A. The directional-asymmetry comparison you were excited about (B_ij big, B_ji≈0
  agreeing with C) **cannot be read off this B**. Getting real direction needs a non-Gaussian method
  (LiNGAM / DirectLiNGAM, or its cyclic variant LiNG) or interventional/temporal data.
- **These are conditional associations, not guaranteed causal effects.** Observational cross-sectional
  data doesn't identify every entry of B causally (you flagged this yourself). B is the correct
  *mathematical* target of your SEM; the causal reading requires the usual SEM assumptions.
- **Partialling changes meaning.** Because every protein is conditioned on all others, a real but
  *mediated* edge (A→M→B) shrinks toward zero in B even if A and B are strongly marginally correlated.
  APP–BACE1 is likely an example. If you also want the marginal (unconditioned) graph, that's a
  separate, easy pass.
- **Pooling + expr_z**: cohort mean differences are removed by within-dataset z-scoring and absorbed
  by dataset dummies, but any cohort-specific *wiring* differences are averaged into one graph by design.

## Natural next step
Add a **DirectLiNGAM** (or cyclic-LiNGAM) layer on the same aligned (X,Z) to get an orientation for
each edge — that's what turns this symmetric weighted B into the directional graph whose asymmetry you
can actually compare to C's Source→Target asymmetry.
