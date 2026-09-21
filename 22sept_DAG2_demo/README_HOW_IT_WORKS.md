# CHRONOS causal graphs — how the pipeline works

This document explains, step by step, exactly what the code in `claude_work/` does to turn
raw proteomics tables into one 80 × 80 directed graph per dataset, how it decides which
"trajectory variables" and which "regions" (subgroups of people) matter, and how the three
graphs are compared. It describes the code as it is now — including its known weak spots
(section 9) — not the history of how it got here.

Every rule below is taken directly from the code. Constants are given with their names in
`chronos_regions.py`, so you can find them.

---

## 0. The whole thing in one paragraph

You have a **confidence graph C** (Debbie's prior): for 80 proteins, a number for every ordered
pair saying how confident we are that protein *i* acts on protein *j*. For each proteomics
dataset separately, the pipeline asks: *do the proteins that C says are linked actually rise
and fall together across people's brains, more than proteins C says are unlinked?* It measures
that with one number (the AUROC "signal"). It then tries to **sharpen** the signal by removing
the influence of other variables (cell-type mix, pathology, …) — keeping a variable only if it
helps more than a scrambled copy of itself would. It then looks for **subgroups of people** in
whom C is unusually visible, comparing each subgroup against random groups of the same size.
Finally it writes the **graph**: edges exist only where C allows them; each edge gets a weight
(C × strength of co-movement), a direction, an effect size, a sample size, and a significance
test. `compare_three.py` then checks how much the per-dataset graphs agree beyond chance.

```
                 ┌──────────────────────┐
 raw Banner ───▶ │ prep_banner.py       │──▶ clean Banner matrices (TMT, LFQ)
 files           └──────────────────────┘
                                                     ┌────────────────────────────┐
 ROSMAP matrix ──┐                                   │ chronos_regions.py         │
 Diverse matrix ─┼─▶ + metadata + confidence graph ─▶│  1 load + ID chain         │
 Banner matrix ──┘                                   │  2 covariates (Z)          │
                                                     │  3 signal (AUROC)          │
                                                     │  4 variable selection (Z*) │
                                                     │  5 region search           │
                                                     │  6 graph W, beta, N, p, q  │
                                                     └─────────────┬──────────────┘
                                                                   │ one set of files per dataset
                                                     ┌─────────────▼──────────────┐
                                                     │ compare_three.py           │
                                                     │  agreement vs chance,      │
                                                     │  consensus pairs, TMT/LFQ  │
                                                     └────────────────────────────┘
```

---

## 1. Files

| file | role |
|---|---|
| `prep_banner.py` | Converts the two raw Banner search outputs into clean protein × sample matrices. Only needed for Banner. |
| `chronos_regions.py` | The main engine. One run = one dataset. |
| `compare_three.py` | Reads the per-dataset outputs and compares them. |
| `banner_exclude_possible_diverse_overlap.txt` | 8 Banner donor IDs dropped from Banner so it shares no people with Diverse (section 2.4). |
| `results/` | All outputs. |
| `check_regions_usable_null.py` | Read-only re-check of the region search with a corrected null (section 9.1). |
| `make_results_readme.py` | Regenerates `README_RESULTS.md` from `results/`. |
| `README.md` | Short run instructions. |

Run order:

```bash
DATA=/home/wordson22/projects/CHRONOS/DATA
PRIOR=../CHRONOS_80x80_confidence.xlsx
python prep_banner.py $DATA data
python chronos_regions.py --dataset ROSMAP  --data_dir $DATA --prior $PRIOR --out results
python chronos_regions.py --dataset Diverse --data_dir $DATA --prior $PRIOR --out results \
       --exclude_individuals results/ROSMAP_individuals.txt
python chronos_regions.py --dataset Banner  --data_dir $DATA --prior $PRIOR --out results \
       --prot $PWD/data/Banner_TMT_log2ratio_GIS.csv \
       --exclude_individuals banner_exclude_possible_diverse_overlap.txt
python chronos_regions.py --dataset BannerLFQ --data_dir $DATA --prior $PRIOR --out results \
       --prot $PWD/data/Banner_LFQ_log2.csv \
       --exclude_individuals banner_exclude_possible_diverse_overlap.txt
python compare_three.py --results results
python check_regions_usable_null.py      # optional: corrected region null (section 9.1)
python make_results_readme.py            # optional: rebuild README_RESULTS.md
```

`--fast` cuts the random repeats (200 → 50 region nulls, 10 → 4 variable nulls) for a quick test.
Only Python with numpy, pandas, scipy, openpyxl is needed.

---

## 2. Inputs

### 2.1 The confidence graph C

`CHRONOS_80x80_confidence.xlsx`: an 80 × 80 table, rows = source protein, columns = target
protein. Raw values are on a 0–5 scale.

`load_prior()`:
1. divides everything by `PRIOR_SCALE = 5`, so C is on 0–1;
2. sets the diagonal to 0 (no self-edges).

An ordered pair *i → j* **is a prior edge** if `C[i,j] >= MIN_PRIOR = 0.05`, i.e. a raw value
of at least 0.25. The raw sheet uses 0.1 as a background "nothing" value
(→ C = 0.02), so background entries are not edges.

With the current sheet that gives:
- **688 directed prior edges** (C from 0.2 to 0.9),
- which cover **507 unordered protein pairs** out of 80·79/2 = **3,160**,
- 326 pairs are **one-way** (only *i → j* is allowed),
- 181 pairs are **two-way** (both *i → j* and *j → i* are allowed; 78 of them with equal C both ways).

The one-way / two-way distinction matters a lot later: for a one-way pair the arrow is fixed
by C; for a two-way pair the data picks the arrow.

### 2.2 The proteomics matrices

All four are "rows = proteins, columns = samples" tables whose row names start with a gene
symbol, optionally followed by `|UniProtID` (`fmt="gene_pipe"`). The part before `|` is
matched to the 80 node names.

| dataset | file | what a value means |
|---|---|---|
| ROSMAP | `C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv` | log2 abundance ratio, already batch-corrected upstream (median polish); 400 samples. |
| Diverse | `n1086_residual_log2_batch.csv` | log2 abundance with batch already regressed out upstream ("residual"); 1,086 samples. |
| Banner | `data/Banner_TMT_log2ratio_GIS.csv` (made by `prep_banner.py`) | log2(sample / pooled-standard channel in the same TMT run), per-sample median-centred; 220 samples. |
| BannerLFQ | `data/Banner_LFQ_log2.csv` (made by `prep_banner.py`) | log2 label-free (LFQ) intensity, per-sample median-centred; 225 samples. |

`load_proteomics()` splits every matrix into two parts:
- **X** — the 80 confidence-graph proteins (or as many as the file contains);
- **rest** — every other protein. The rest is never part of the graph, but it is used to
  build trajectory variables (section 4.2).

If a gene appears on several rows, the first row is kept (for Banner, `prep_banner.py` has
already kept the best-measured row per gene).

### 2.3 Metadata and the ID chain

A proteomics column is a *sample*, but clinical data is per *person*. The chain is:

```
proteomics column  ──(assay file: batchChannel)──▶  specimenID  ──(biospecimen file)──▶  individualID  ──▶  clinical row
```

`sample_to_individual()`:
- If a column name is already a specimenID in the biospecimen file, it is used directly (Banner LFQ).
- Otherwise the column name is looked up as a `batchChannel` in the assay file (e.g. `b14.127C`)
  to find the specimenID (ROSMAP, Diverse, Banner TMT). If the same batchChannel appears in more
  than one data release, the entry whose specimen exists in the biospecimen file wins.
- Samples that do not reach an individualID are dropped.

`collapse_to_individuals()` then makes **one row per person**: if a person has several samples
(replicates), their values are averaged protein by protein (missing values ignored).

Per-dataset metadata files:

| dataset | assay | biospecimen | clinical |
|---|---|---|---|
| ROSMAP | `ROSMAP_assay_proteomics_TMTquantitation_metadata.csv` | `ROSMAP_biospecimen_metadata.csv` | `ROSMAP_clinical.csv` |
| Diverse | `AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata_260622.csv` | `AMP-AD_DiverseCohorts_biospecimen_metadata.csv` | `AMP-AD_DiverseCohorts_individual_metadata.csv` |
| Banner | `Banner_TMTquantitation_assay_metadata.csv` | `Banner_biospecimen_metadata.csv` | `Banner_individual_metadata.csv` |
| BannerLFQ | `Banner_proteomics_assay_metadata.csv` | `Banner_biospecimen_metadata.csv` | `Banner_individual_metadata.csv` |

**Batch.** For ROSMAP/Diverse/BannerLFQ the `batch` column of the assay file is attached to each
person. For Banner TMT (`batch_from_channel=True`) the batch is the TMT run taken from the
batchChannel prefix (`b14.127C` → `b14`), because the assay file's `batch` column (b1–b4) is a
sample-preparation batch, not the TMT run.

### 2.4 Keeping the three datasets independent

`--exclude_individuals FILE` drops the listed individualIDs after collapsing to people.
- **Diverse** drops everyone in `ROSMAP_individuals.txt` (ROSMAP donors also appear in Diverse
  under the same IDs): 16 people.
- **Banner** drops the 8 IDs in `banner_exclude_possible_diverse_overlap.txt`. Diverse contains
  donors from the Banner brain bank under a different ID system, with no table linking the two.
  These 8 Banner donors match a Diverse Banner-cohort donor on sex + APOE genotype + age at death +
  Braak stage, so they are treated as possible duplicates.

---

## 3. `prep_banner.py` — making the Banner matrices

### 3.1 TMT (`bannertmt_22batchmulticoncensus_Proteins.txt`)

This is a Proteome Discoverer export: 11,518 protein rows, and for each of 22 TMT runs
(called F1…F22 in the file) 11 channels: 126 (the pooled "GIS" standard) and 127N…131C (10 samples).

For each run *k* and each sample channel *ch*:

```
value = log2( Abundance(Normalized)[F_k, ch] / Abundance(Normalized)[F_k, 126] )
```

- Why divide by channel 126: every run contains the same pooled standard, so dividing by it
  cancels run-to-run differences (the same idea used upstream for ROSMAP).
- Zero or missing abundances become missing.
- Then every sample (column) has its median subtracted, so all samples are centred on 0
  ("equal loading").
- Columns are renamed `F{k}: {ch}` → `b{k:02d}.{ch}` (e.g. F14/127C → `b14.127C`), which is the
  `batchChannel` format of the assay metadata.
- Gene symbol = first entry of the `Gene Symbol` column; row name = `GENE|Accession`.
- If a gene has several rows (isoforms), the row with the most observed samples is kept.

Result: 9,728 proteins × 220 samples, 20.9% missing.

**Checking F_k = b_k.** Nothing in the file states that run F14 is metadata batch b14. The check:
average the Y-chromosome-only proteins (RPS4Y1, DDX3Y, EIF1AY, KDM5D, USP9Y, NLGN4Y) per sample,
and see whether that separates people recorded as male from female. With the F_k = b_k mapping,
a single cut-off classifies 99.0% of 198 matched samples correctly. If the batch labels are
randomly shuffled 200 times, the best achievable accuracy averages 57.8% (maximum 64.1%). So the
mapping is correct.

### 3.2 LFQ (`banner_proteomics_pfc_proteinoutput.txt`)

A MaxQuant-style export: 5,711 protein rows, 225 `LFQ.intensity.<specimen>` columns.

```
value = log2(LFQ intensity),   0 → missing,   then per-sample median centring
```

- Column `LFQ.intensity.b1_002_07` → `b1_002_07_lfq`, which is the specimenID used in the
  biospecimen file.
- 27 of the 225 columns are pooled-standard runs (`…bgis…`, `…egis…`); they have no individualID
  and are dropped by the ID chain.
- Gene = first entry of `GN`; isoform rule as above.

Result: 5,041 proteins × 225 samples, 27.1% missing. 58 of the 80 nodes are present.

---

## 4. `chronos_regions.py` — stage by stage

### 4.1 Cleaning clinical columns into candidate variables (`clean_covariates`)

Every column of the clinical file (plus `batch`) is a candidate trajectory variable **unless**:

- its name matches an identifier/bookkeeping pattern (`NOT_Z`: individualID, specimenID, projid,
  platform, batchChannel, controlType, assay, organ, tissue, species, …), or it is one of the 80 protein names;
- it has fewer than `MIN_REGION_N = 40` non-missing values, or only one distinct value;
- one value covers more than 95% of people.

Then each surviving column is converted to numbers:
- Leading `>=`, `>`, `<=`, `<` and trailing `+` are stripped (so `90+` and `>=90` become 90).
- "nan", "NA", "missing or unknown", "unknown", "not applicable", "" become missing.
- If at least 80% of the non-missing entries are numbers (and there are ≥ 2 distinct numbers),
  the column is **numeric**.
- Otherwise it is **categorical**: its text levels are sorted in "natural" order (numbers by
  value, Roman numerals I–VI by value, so "Stage II" < "Stage IV") and coded 0, 1, 2, …
  A categorical column with more than 12 levels is dropped as free text.
  **Consequence:** batch columns with more than 12 batches are dropped. ROSMAP (64 batch labels in
  the assay file), Diverse (91) and Banner TMT (22 TMT runs) all exceed 12, and the BannerLFQ assay
  file has no batch column, so batch never ends up as a candidate in any of the four runs.

Each variable is also tagged (only for reporting) as `technical` (name contains batch, pmi, ph,
study, cohort, datacontribution), `proteome-derived` (cell-type scores and PCs), or `clinical/pathology`.

### 4.2 Variables built from the rest of the proteome (`proteome_covariates`)

The ~5–7 thousand proteins that are **not** graph nodes carry information about what each brain
sample is made of. Two kinds of variable are built from them, per person:

1. **Cell-type scores.** Each "rest" protein is z-scored across people
   ((value − mean) / sd). For each cell type, the score is the average z-score of its marker
   proteins that are present (at least 2 needed):

   | score | markers |
   |---|---|
   | `celltype_neuron` | SNAP25, SYT1, STMN2, GAP43, SYN1, RBFOX3, NEFL, CAMK2A |
   | `celltype_astrocyte` | GFAP, ALDH1L1, AQP4, SLC1A2, SLC1A3, GJA1 |
   | `celltype_microglia` | AIF1, CD68, ITGAM, CSF1R, P2RY12, CX3CR1, HLA-DRA |
   | `celltype_oligodendrocyte` | MBP, MOG, PLP1, MAG, CNP, MOBP |
   | `celltype_endothelial` | CLDN5, VWF, PECAM1, FLT1, ESAM |

   Intuition: a brain sample with more white matter has more oligodendrocyte proteins; a sample
   with more neuron loss has fewer neuron proteins. These mixture differences move hundreds of
   proteins together, so two proteins can co-vary just because both are made by the same cell
   type — not because one acts on the other.

2. **Proteome principal components (PCs).** Take only "rest" proteins observed in ≥ 95% of
   people, z-score them, fill the few gaps with 0, centre, and take a singular value
   decomposition. The top `N_PROTEOME_PCS = 10` components (score = U·S) become
   `proteomePC1` … `proteomePC10`. Each PC is labelled with the cell-type score it correlates
   with most (Spearman), purely to help you read it; the label is not used in any calculation.
   Intuition: PCs are the biggest "whole-proteome" patterns of variation between people, whatever
   their cause (cell mix, tissue quality, disease, leftover technical effects).

All of these are added to the candidate list. The log prints how much variance each PC explains
and what it tracks.

### 4.3 Co-movement between two proteins: pairwise Spearman (`spearman_pairwise`)

For every pair of proteins (*i*, *j*):
- use only people in whom **both** are measured ("pairwise-complete"); their count is
  **n_ij**, stored in the N matrix;
- Spearman ρ = the correlation of the **ranks** of the two proteins across those people.
  Ranks make it robust to outliers and to any monotone rescaling of the data.

Implementation detail: each protein is ranked once over all its observed people, then a Pearson
correlation of those ranks is computed on the shared people using matrix algebra with a 0/1
"observed" mask (N = MᵀM, sums and sums of squares restricted by the mask). This is exact when a
pair has no extra missingness, and very close to exact otherwise (the ranks are not re-computed
within each pair's subset). It is fast enough to be run thousands of times, which the searches need.

No person is ever dropped because some *other* protein is missing. (If one demanded all 80
proteins, ROSMAP would keep only 40 of 400 people, Banner none.)

### 4.4 The signal score: AUROC against C (`auroc_vs_prior`)

This single number is what everything else tries to improve.

Take all 3,160 unordered protein pairs with n_ij ≥ `MIN_PAIR_N = 20`. Call a pair an
**edge** if C allows it in at least one direction (507 pairs), a **non-edge** otherwise.
Sort all pairs by their correlation ρ_ij. Then

```
AUROC = P( a randomly chosen edge has a higher ρ than a randomly chosen non-edge )
      = ( Σ ranks of edges − n1(n1+1)/2 ) / ( n1 · n0 )        (Mann–Whitney form)
```

where n1 = number of edge pairs, n0 = number of non-edge pairs.

- 0.5 → the data's co-movement tells you nothing about C.
- 1.0 → every C edge co-moves more than every non-edge.
- Below 0.5 → C edges co-move *less* than non-edges.

**Signed, not absolute.** `SIGNED = True`: the ranking uses ρ itself, not |ρ|. So an edge scores
well only if its two proteins rise *together* (positive ρ). This was chosen because C's edges are
overwhelmingly "same direction" relationships, and because signed ρ recovers C much better than
|ρ| in these data (e.g. ROSMAP raw: signed 0.562 vs unsigned 0.515). (The file header still says
"|Spearman rho|"; the code uses signed ρ.)

Small worked example: 3 edge pairs with ρ = 0.30, 0.10, 0.05 and 3 non-edge pairs with
ρ = 0.20, 0.00, −0.10. Of the 9 edge-vs-non-edge comparisons, the edge wins 7 (0.30 beats all
three; 0.10 beats 0.00 and −0.10; 0.05 beats 0.00 and −0.10), so AUROC = 7/9 = 0.78.

### 4.5 Adjusting for variables (`residualise`, `design`)

"Adjusting for Z" = removing from every protein the part that can be predicted linearly from Z,
and then correlating what is left.

For a set of variables Z:
1. `design()` turns them into a numeric matrix: numeric variables as they are, categorical
   variables as 0/1 indicator columns (one per level, first level dropped), plus an intercept.
2. Only people with **every** Z variable observed are used. Everyone else becomes entirely
   missing in the adjusted data. (This is deliberate — raw and adjusted values are never mixed —
   but it means a variable that is missing for many people shrinks the adjusted dataset. See 9.3.)
3. For each protein separately, using those people where that protein is observed
   (at least 20), fit ordinary least squares `protein ≈ a + b·Z` and keep the residual
   `protein − (a + b·Z)`.

The adjusted matrix is then fed back into the pairwise Spearman and the AUROC.

Intuition: if two proteins both track "how much neuron is in this sample", removing neuron
content removes the co-movement they share for that reason, and what remains is closer to
protein-to-protein relationships.

### 4.6 Choosing trajectory variables Z* (`select_variables`)

Greedy forward selection, at most `MAX_Z = 5` variables.

At each step, for every remaining candidate c:
1. **Real score:** adjust for (already selected + c), compute AUROC → s.
2. **Shuffled score:** repeat `N_NULL_VAR = 10` times: take c, randomly permute its values
   among the people who have it (same values, same missingness, but the link to who-is-who is
   broken), adjust for (already selected + shuffled c), compute AUROC. This gives a null mean
   μ and sd σ.
3. **z = (s − μ) / σ.**

A candidate is **eligible** if z > 2 **and** it raises the AUROC by more than 0.0001 over the
current value. Among eligible candidates, the one with the **highest AUROC** (not the highest z)
is added. Selection stops when no candidate is eligible, or after 5 variables.

Why the shuffle comparison is needed: adding *any* variable to a regression removes a little
variation from every protein, and a variable with missing values also changes *which* people are
analysed. Both can move the AUROC even if the variable has nothing to do with biology. The
shuffled copy has exactly the same mechanical side effects, so only the part of the gain that
depends on *which person has which value* is credited.

Reading z: with only 10 shuffles, σ is estimated roughly; when the shuffled AUROCs are nearly
identical (σ tiny), z can be huge (e.g. 775). Treat z as "clearly beats noise" vs "doesn't",
not as a precise size.

Every candidate tried at every step is written to `{ds}_variables.csv` with: step, variable,
adjusted_for, auroc, null_mean, gain_vs_current, z_vs_shuffled, n_subjects (people with the full
Z set), kind.

After selection, the data are adjusted for Z* once more to give **Xfin**, the matrix used for
regions and for the final graph.

### 4.7 Regions: who shows C most strongly (`candidate_regions`, `search_regions`)

A **region** is a subgroup of people defined by one variable (depth 1) or two (depth 2),
evaluated on the Z*-adjusted data.

**Depth-1 candidates**, for every candidate variable except batch and proteome PCs (PCs are not
interpretable as groups of people):
- categorical, or numeric with ≤ 6 distinct values → one region per level (`Braak = Stage IV`);
- other numeric → five regions: low third, middle third, high third (cut at the 1/3 and 2/3
  quantiles), lower half, upper half (cut at the median).

A region is kept only if it has at least `MIN_REGION_N = 40` people and at most 90% of everyone.

**Scoring a region:**
1. AUROC of the adjusted data restricted to the region's people → s.
2. A **same-size random null**: draw `N_NULL_REGION = 200` random subsets of people of the same
   size and compute their AUROC. To save time, sizes are rounded to the nearest 10 (minimum 40)
   and each size's null is reused.
3. z = (s − null mean) / null sd, and p = (1 + #null ≥ s) / 201.

Why same-size: small groups give noisier AUROCs, so a small group can look extreme by luck.
Comparing to random groups of the same size removes that effect. Note p can never be below
1/201 ≈ 0.005.

**Depth 2:** the 5 depth-1 regions with the highest z are each intersected with every other
depth-1 region from a *different* variable (duplicates removed, ≥ 40 people), and scored the same way.

**Multiple testing:** all regions (depth 1 and 2) get Benjamini–Hochberg q-values from their p's.
A region is flagged **STRONG** if q < 0.10 and z > 2. The log also lists the 3 regions with the
lowest z ("C least visible").

If any region is STRONG, the graph is refitted inside the top one and saved as
`{ds}_best_region_W.csv`.

Everything is saved in `{ds}_regions.csv`: region, depth, n, auroc, null_mean, null_sd, z,
p_vs_random_same_size, q_fdr, strong.

### 4.8 Building the graph (`build_graph`)

Computed on Xfin (the adjusted data). For every unordered pair (*i*, *j*):

1. **Correlation test.** ρ_ij and n_ij from the pairwise Spearman. p-value from the t
   approximation: t = ρ·√((n−2)/(1−ρ²)), two-sided with n−2 degrees of freedom. This p is stored
   for both directions.
2. If C allows neither direction → no edge, W = 0 both ways.
3. If n_ij < 20 → a **placeholder**: W = 0.1 × C in the direction C prefers (no data used).
4. Otherwise run the **direction test** (below) on the people who have both proteins, and set W
   by this table (r = |ρ_ij|):

   | C allows | direction test says | result |
   |---|---|---|
   | only *i → j* | *i → j* | W[i,j] = C[i,j] · r |
   | only *i → j* | *j → i* | W[i,j] = 0.3 · C[i,j] · r (kept in C's direction, down-weighted) |
   | both | *i → j* | W[i,j] = C[i,j] · r, W[j,i] = 0 |
   | both | *j → i* | W[j,i] = C[j,i] · r, W[i,j] = 0 |

   So every C pair gets **exactly one arrow**. For one-way pairs, the arrow is always C's
   direction; the data can only shrink the weight. For two-way pairs, the data picks the arrow.

**The direction test (`lingam_direction`)** — a bivariate LiNGAM-style score.
Idea: if *x* causes *y* (y = b·x + noise) and the noise is independent of *x*, then the residual
of regressing y on x is independent of x; the residual of the *wrong* regression (x on y) is not.
With non-Gaussian data this difference is detectable through higher moments.

In code, with x and y standardised and c = mean(x·y):
```
e_xy = y − c·x          (residual of y on x)
e_yx = x − c·y          (residual of x on y)
score = |cov(e_yx, y²)| − |cov(e_xy, x²)|
score > 0  →  x → y   (the backward residual is more dependent on its regressor)
```
It uses only two proteins at a time, ignores all other proteins, and relies on non-Gaussian
distributions. Treat its output as weak evidence.

### 4.9 Effect sizes β (`compute_beta`)

For each target protein *j*, its **parents** are all *i* with W[i,j] > 0. β is the coefficient
of each parent in an ordinary least squares regression

```
protein_j ≈ a + Σ_parents β_ij · protein_i
```

on the people who have the target and **all** its parents. If fewer than (number of parents + 10)
people qualify, β falls back to a one-parent-at-a-time slope on the people who have that pair
(≥ 20). β is in the data's own units (log2 ratio per log2 ratio), not standardised, so it is not
bounded by 1 and not comparable across datasets with different scaling.

### 4.10 The edge table and "data-supported" (`edge_table`)

One row per **directed prior edge** (688 rows): source, target, prior_C, W, beta, n_pairwise,
p_corr, q_fdr, data_supported.
- q_fdr = Benjamini–Hochberg over the 688 p-values (the two directions of a two-way pair share
  one p-value, so it appears twice).
- **data_supported = 1** if q_fdr < 0.05 **and** W > 0 **and** n_pairwise ≥ 20.
  Because only one direction of each pair has W > 0, a pair is counted once.

Proteins missing from a dataset are written as all-zero rows/columns, so every output matrix is 80 × 80.

### 4.11 Outputs per dataset (`results/`)

| file | content |
|---|---|
| `{ds}_causal_W.csv` | 80 × 80 weights, row = source, column = target |
| `{ds}_causal_beta.csv` | 80 × 80 effect sizes |
| `{ds}_causal_N.csv` | 80 × 80 number of people behind each pair |
| `{ds}_causal_edges.csv` | the edge table (4.10) |
| `{ds}_variables.csv` | every variable tried at every step (4.6) |
| `{ds}_regions.csv` | every region tested (4.7) |
| `{ds}_best_region_W.csv` | graph refitted inside the top STRONG region (only if one exists) |
| `{ds}_individuals.txt` | the individualIDs used |
| `{ds}_report.txt`, `{ds}.log` | the printed log |
| `{ds}_summary.json` | headline numbers |
| `{ds}_regions_usable_null_check.csv` | re-scored regions with the corrected null (section 9.1) — ROSMAP, Diverse, Banner |

---

## 5. `compare_three.py` — comparing the graphs

It reads `{ds}_causal_edges.csv` and `{ds}_summary.json` for ROSMAP, Diverse and Banner, plus
BannerLFQ for one extra check. Everything is based on **data-supported edges** (4.10); comparing
all W > 0 edges would be meaningless because C alone decides which ~505 edges get W > 0.

1. **Signal table:** people, raw AUROC, adjusted AUROC, Z* per dataset.
2. **Regions:** number of STRONG regions and the top 3 per dataset.
3. **Supported edges:** count / 688 per dataset.
4. **Pairwise agreement.** Supported edges are turned into unordered pairs. For datasets A and B:
   - Jaccard = |A ∩ B| / |A ∪ B|.
   - **Chance baseline:** U = pairs testable (n ≥ 20) in both; a = |A ∩ U|, b = |B ∩ U|.
     If each dataset's supported pairs were a random pick from U, the expected overlap is
     a·b / |U| and the probability of seeing at least the observed overlap is the
     hypergeometric tail `P(X ≥ shared)`, X ~ Hypergeometric(|U|, a, b). Reported: shared,
     expected, fold = shared / expected, p.
   - **Same direction:** for shared pairs, whether the arrow points the same way in both graphs.
     For one-way C pairs this is automatically yes, so the disagreement comes only from two-way pairs.
5. **Consensus pairs:** unordered pairs supported in ≥ 2 datasets, written to
   `consensus_edges.csv` with each dataset's arrow, W and β, whether the arrows agree
   (`direction_consistent`), and sorted by number of datasets then mean W.
6. **Technical check (Banner TMT vs Banner LFQ):** same people, two measurement methods. Among
   pairs testable in both: overlap of supported pairs vs chance (hypergeometric), and the
   Spearman correlation of |W| across all those pairs.

Outputs: `comparison_report.txt`, `comparison_pairs.csv`, `consensus_edges.csv`.

---

## 6. Parameters

| name | value | used for |
|---|---|---|
| `PRIOR_SCALE` | 5 | raw C / 5 |
| `MIN_PRIOR` | 0.05 | C ≥ this is an edge |
| `MIN_PAIR_N` | 20 | pairs with fewer shared people are not scored / tested |
| `MIN_REGION_N` | 40 | smallest region; also minimum people for adjustment |
| `N_NULL_VAR` | 10 | shuffles per candidate variable |
| `N_NULL_REGION` | 200 | random subsets per region size |
| `MAX_Z` | 5 | maximum number of selected variables |
| `SIGNED` | True | AUROC uses signed ρ |
| `N_PROTEOME_PCS` | 10 | proteome PCs offered as candidates |
| variable z threshold | 2 | selection rule |
| region STRONG | q < 0.10 and z > 2 | region rule |
| edge supported | q < 0.05, W > 0, n ≥ 20 | edge rule |
| down-weight | 0.3 | one-way edge where the direction test disagrees with C |
| placeholder | 0.1 × C | pairs with n < 20 |
| random seed | 22 | `RNG` |

---

## 7. What the graph does and does not learn from data

- **Which edges can exist:** only C decides.
- **Edge weight:** C × |ρ|, so data and C contribute equally (multiplicatively).
- **Arrow direction:** C for one-way pairs (326 of 507), the bivariate direction test for two-way
  pairs (181 of 507).
- **Whether an edge counts as real:** data only (FDR-corrected correlation test).
- **Trajectory variables and regions:** chosen by how much they improve recovery of C.

In other words the output is "which of C's links are visible in this dataset, how strongly, and
in whom", not an independent discovery of causal structure.

---

## 8. Why each dataset gets its own run

The three datasets differ in technology details, preprocessing, cohorts and clinical variables,
and their people do not overlap (section 2.4). Each is fitted from scratch; nothing from one
dataset's fit is used in another. That is what makes their agreement (section 5) informative.

---

## 9. Known weak spots in the current code

These are properties of the code as it stands; none has been changed.

**9.1 Region null includes people who have no adjusted data.** When Z* contains a variable that
is missing for some people, those people are entirely missing in Xfin (4.5). The region
AUROC is computed on whoever in the region has data, but the random same-size subsets are drawn
from *all* people, so a random "group of 100" may contain only ~70 usable people and is
therefore noisier and lower-scoring than it should be. This makes regions look stronger than they
are. It matters only in Diverse, where `amyAny` is in Z* and 288 of 964 people lack it (676
usable). Re-scoring every region with a null drawn only from usable people (and matched to the
region's usable size) is saved as `{ds}_regions_usable_null_check.csv`. Result: Diverse STRONG
regions go from 28 to 0 (best corrected q = 0.12); ROSMAP and Banner stay at 0.

**9.2 Regions with no usable people get p = 0.005.** If a region contains no one with adjusted
data (in Diverse: the Mayo cohort, `mayoDx`, `derivedOutcomeBasedOnMayoDx = TRUE`), its AUROC is
NaN, the comparison `null ≥ NaN` is always false, and p becomes the minimum 1/201. These rows are
not flagged STRONG (z is NaN), but their tiny p-values enter the Benjamini–Hochberg step and lower
everyone else's q. Part of the reason Diverse's original q-values reached 0.04.

**9.3 The adjusted Diverse analysis has 676 people, not 964.** Same cause as 9.1. The summary's
`n_people` reports 964 (people loaded), but everything after `amyAny` was selected — final
adjusted AUROC, regions, and the graph — uses the 676 people with amyloid data (n_pairwise for
most Diverse edges is 676).

**9.4 `Diverse_best_region_W.csv` comes from a region that is not STRONG under 9.1's correction.**

**9.5 Batch is never a candidate.** Categorical columns with more than 12 levels are dropped (4.1);
ROSMAP, Diverse and Banner TMT all have more than 12 batches, and BannerLFQ has no batch column. The data files are already
batch-corrected (ROSMAP, Diverse) or pool-normalised (Banner), so this is mostly harmless, but it
is not tested.

**9.6 Direction is weak.** The direction test is bivariate, ignores other proteins, and assumes
non-Gaussian independent noise. For one-way C pairs it disagrees with C about half the time
(the edge is then kept at 0.3 weight), and for two-way pairs different datasets often pick
different arrows.

**9.7 Small approximations.** Spearman on pre-computed ranks (4.3); only 10 shuffles for the
variable null (z imprecise); region sizes rounded to the nearest 10 for the null; p-values in the
region search can't go below 0.005; the edge FDR counts each two-way pair's p twice; β is
unstandardised.

**9.8 Greedy selection picks the highest AUROC among eligible variables**, not the highest z.
A variable with a large, clearly real gain can be passed over for one with a slightly higher
AUROC, and variables are never removed once added.
