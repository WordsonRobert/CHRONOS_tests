# Per-protein × clinical CSVs (ROSMAP + Diverse Cohorts + Banner-TMT)

One CSV per requested clinical variable. Every row is **one protein in one person**.
Built by `build_protein_clinical.py`. No value is invented; every mapping below is a
fixed, deterministic crosswalk applied to the raw metadata as delivered.

## Columns (identical in every CSV)

| column | meaning |
|---|---|
| `protein` | one of the 80 CHRONOS graph proteins (rows are grouped protein-by-protein) |
| `dataset` | `ROSMAP`, `Diverse`, or `Banner` |
| `individualID` | the person that this protein value came from |
| `expr_raw_log2` | that person's log2 abundance for this protein (one row per person; replicate specimens averaged). ROSMAP/Diverse as delivered by the consortium; Banner recomputed by `prep_banner.py` (log2 of channel ÷ plex reference 126). |
| `expr_z` | `expr_raw_log2` **z-scored within its own dataset, per protein** (mean 0, sd 1). This is the normalization that makes the three datasets comparable, since each has its own log2 centering. Use this column for cross-dataset work. |
| `<variable>` | the clinical variable this file is named after, harmonized to identical labels across datasets (see below) |
| `impaired` | 1 = AD **or any other cognitive decline**, 0 = cognitively normal |

Rows where either the expression **or** the variable value was missing are dropped.

## `impaired` flag (1 = AD/any cognitive decline, 0 = normal)
- **ROSMAP** `cogdx`: 1 → 0 (no impairment); 2,3,4,5,6 → 1 (MCI / AD / other dementia).
- **Diverse** `ADoutcome`: `Control` → 0; `AD` → 1; `Other` → 1; `missing` → dropped.
- **Banner** `diagnosis`: `control` → 0; `Alzheimer Disease` → 1. (Banner recorded **no MCI/intermediate** category, so Banner's `impaired`=1 means AD specifically, vs ROSMAP/Diverse where it also includes MCI/other. Keep this asymmetry in mind.)

## Variable harmonization + which datasets have each

| CSV | datasets | harmonization |
|---|---|---|
| `sex` | R, D, B | ROSMAP `msex` 1→male, 0→female; others already male/female |
| `apoe_genotype` | R, D, B | canonical 2-digit (`33`,`34`,`44`,…). Banner `e3-4`→`34`. **Unordered category, not a number.** |
| `cts_mmse30_lv` | R, B | 0–30. ROSMAP `cts_mmse30_lv`; Banner `lastMMSE`. Same scale. |
| `cts_mmse30_first_ad_dx` | R only | 0–30 (only 126 ROSMAP people were ever AD-diagnosed, so only they have this) |
| `pmi` | R, D, B | hours, numeric. **Ranges differ a lot** (Banner 1–6 h, ROSMAP up to ~98 h). |
| `braaksc` | R, D, B | integer 0–6. Diverse `Stage I..VI`→1..6. **Only ROSMAP has stage 0** (others didn't record it). |
| `ceradsc` | R, B | **severity 0–3** (0=none/no-AD … 3=frequent/definite). ⚠️ **ROSMAP `ceradsc` is REVERSE-coded** (1=definite…4=no-AD) so it was crosswalked 1→3, 2→2, 3→1, 4→0. Banner `CERAD` was already 0–3. Verify this crosswalk before trusting it. |
| `amyThal` | D only | `Phase N`→N (1–5) |
| `amyA` | D only | `Thal Phase 1 or 2`→1, `3`→2, `4 or 5`→3 |
| `amyCerad` | D only | `C0`→0 … `C3`→3 (same 0–3 CERAD concept as `ceradsc`, just Diverse's version) |
| `amyAny` | D only | 0/1 (any amyloid present) |
| `reag` | D only | NIA-Reagan: `No AD`→0, `Low`→1, `Intermediate`→2, `High Likelihood`→3 |
| `PlaqueTotal` | B only | numeric plaque density |
| `TangleTotal` | B only | numeric tangle density |

R = ROSMAP (400 people), D = Diverse (980), B = Banner-TMT (190).
BannerLFQ (same donors, different instrument) was left out — it's a robustness variant, not a 4th cohort.
