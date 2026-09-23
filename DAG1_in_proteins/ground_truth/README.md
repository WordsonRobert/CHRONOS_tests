# Adaptive trajectory-region discovery — 23 Sept methodology

Data-driven version of the "confidence-score → trajectory regions" idea. Instead of
hand-partitioning MMSE/Braak/etc. into fixed bins, the molecular edges themselves
decide where along each trajectory the behaviour changes.

## What it does (per dataset, per ordered variable Z_k)
1. **Rank-normalise** Z_k → u ∈ [0,1] (subjects sorted along the axis).
2. **Rolling local slope**: slide a window of `W=40` subjects; inside each window fit
   `β_i(u)` = slope of protein i's normalised level (`expr_z`) vs Z_k, for all 80 proteins.
3. **Edge behaviour** (sign-of-confidence aware — the local statistic is the
   in-window Pearson r of each protein vs the axis):
   - **positive edge (+1…+5)**: expected behaviour = the two proteins **co-move**
     (`sgn r_i == sgn r_j`). Region = where confidence-c edges co-move.
   - **negative edge (−1,−2)**: the pair is believed *unrelated*, so the expected
     behaviour is that it **stays inert** — BOTH proteins **flat / no change**
     (`|r_i| < FLAT_R` and `|r_j| < FLAT_R`, FLAT_R=0.15 ≈ no local trend at n=40).
     Region = where confidence-c edges are flat. (This is the fix to the first pass,
     which wrongly asked negative edges to co-move and so returned empty S₋₁/S₋₂.)
4. **Aggregate by literature score** c ∈ {−2,−1,+1,+2,+3,+4,+5} (read from
   `CHRONOS_80x80_confidence.xlsx`, sheet `Matrix_80x80`):
   `p_c(u)` = fraction of confidence-c edges that co-move at u
   (denominator = edges whose both slopes are defined in that window).
5. **Discover regions**: `S_c = { u : p_c(u) > τ_c }`, merged into intervals.
   `τ_c` is **adaptive per (dataset,variable,score)**: `τ_c = clip(mean_c + 1.25·sd_c, 0.62, 0.97)`
   — "where does confidence-c co-movement rise clearly above its own baseline".
   A region must also stay above `τ_c` for **≥3 consecutive windows** (not a 1-window flicker)
   and cover **≥15 people**.

## Files
- `discovered_regions.csv` — the master table: `dataset, variable, confidence, u_lo, u_hi,
  value_lo, value_hi` (the u-range mapped back to real variable values), `support_fraction`
  (mean p_c in the region), `tau_c`, `n_people`, `n_edges_c`, `n_subjects`.
- `curves/<dataset>_<variable>.csv` — the full `p_c(u)` curve for every score (for plotting).
- `baseline_pc.csv` — mean & max p_c per (dataset,variable,score) for context.
- `_params.json` — the parameters used.
- `trajectory_regions.py` — the implementation (all params at the top).

## Headline result (this is the important bit)
Co-movement rises **monotonically with literature confidence**, both globally and in the
discovered regions:

| confidence | # regions | mean support |
|---|---|---|
| −2 / −1 (flatness rule: pair inert) | 7 / 8 | ~0.73 |
| +1 … +4 (co-move) | 13 / 2 / 13 / 15 | ~0.69 |
| **+5** (co-move) | **22** | **0.77** |

80 regions total; 57 in Diverse (biggest cohort, most amyloid variables), 17 ROSMAP, 6 Banner.
Positive side: literature-high edges co-move far more, and harder, than low ones (global
co-movement baseline ~0.50 for +1 → 0.58 for +5). Negative side: the "unrelated" pairs stay
flat precisely at **mid-to-late Braak (3–5)** — the same disease-progression window where the
+5 pairs co-move. Both rules agree on where the action is, which is the coherent picture you want.

## Caveats — read before trusting individual regions
- **PMI dominates the strongest +5 regions.** PMI is post-mortem interval — a *technical*
  variable. Proteins co-degrading at certain PMIs is not disease biology. Treat PMI regions as
  a likely technical artifact; the biologically meaningful +5 regions are **Braak, amyThal,
  amyCerad, MMSE**.
- **Coarse ordinals are blocky.** Braak (0–6), CERAD (0–3), amyThal (1–5), amyA (1–3), reag
  (0–3) have few levels, so the rolling slope is coarse (windows spanning one level give no
  slope). Continuous axes (MMSE, PMI, Plaque/Tangle counts) give the cleanest trajectories.
- **3 variables were skipped** — `sex` (binary), `apoe_genotype` (unordered category),
  `amyAny` (0/1 flag): a rolling slope along them is meaningless. They are NOT trajectory axes.
- **Rule is confidence-sign aware, but NOT Sign_80x80 aware yet.** Positive edges are still
  scored by any same-direction co-movement, so a genuine *inhibitory* edge (Sign=−1 but high
  positive confidence — the rubric keeps confidence and sign independent) shows opposite
  directions and is currently undercounted in the + sets. Using `Sign_80x80` so inhibitory
  edges expect *opposite* directions is the next refinement (a per-edge expected-sign test).
- **Not yet the 14-D intersection.** This is steps 1–7 of the plan (per-axis 1-D discovery).
  The multi-variable region construction (step 8+, e.g. {MMSE∈[..], Braak∈[..]}) is the next
  step, to run once you've eyeballed which of these 1-D regions are worth combining.

## Reproduce
`python trajectory_regions.py` from this folder (or with its full path from anywhere). It reads `../DATA/`
(the 14 CSVs) and `../CHRONOS_80x80_confidence.xlsx`, and writes its outputs next to itself; a rerun
reproduces the committed `discovered_regions.csv`, `curves/` and `baseline_pc.csv` exactly. Paths are relative
to the script, so renaming this folder does not break it.

Note: these regions are exploratory. They are defined with the literature matrix C (via its confidence classes),
so they are not independent of C and are not a ground truth for the DAGs.
