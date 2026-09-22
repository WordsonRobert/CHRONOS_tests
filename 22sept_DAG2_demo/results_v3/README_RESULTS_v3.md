# CHRONOS causal graphs — all results (v3 prior)

Everything the pipeline produced, dataset by dataset, with as little interpretation as possible.
How each number is computed is in `README_HOW_IT_WORKS.md`; section numbers like (HIW 4.7) point there.
All tables below were generated from the files in `results_v3/` by `make_results_readme.py`,
using the confidence matrix `CHRONOS_80x80_confidence_v3.xlsx`.

Region results are given twice: the pipeline's original scoring, and a corrected re-scoring that
draws the random comparison groups only from people who have adjusted data (HIW 9.1). Section 5.1
has both counts.

Datasets: **ROSMAP**, **Diverse** (AMP-AD Diverse Cohorts), **Banner** (TMT) are the three
independent graphs. **BannerLFQ** is the same Banner people measured with a second method
(label-free), included as a technical check only.

---

## 1. Data as loaded

| | ROSMAP | Diverse | Banner (TMT) | BannerLFQ |
|---|---|---|---|---|
| proteomics samples in file | 400 | 1,086 | 220 | 225 |
| proteins in file (unique gene symbols) | 8,252 | 9,152 | 9,728 | 5,041 |
| samples matched to a person | 400 | 1,086 | 198 | 198 |
| people after averaging replicates | 400 | 980 | 198 | 198 |
| people excluded as shared with another dataset | 0 | 16 (ROSMAP donors) | 8 (possible Diverse donors) | 8 |
| **people analysed** | **400** | **964** | **190** | **190** |
| people with usable data after adjustment for Z* | 400 | **676** | 189 | 190 |
| graph nodes found (of 80) | 80 | 80 | 80 | 58 |
| missing fraction, 80 nodes × people (raw) | 5.6% | 1.0% | 4.9% | 10.9% |
| nodes with no missing values | 58 | 69 | 66 | 42 |
| people with all nodes measured | 40 | 525 | 0 | 0 |
| median people per protein pair (raw) | 400 | 964 | 190 | 190 |
| smallest people per protein pair (raw) | 120 | 593 | 8 | 0 |
| "rest-of-proteome" proteins used for PCs (≥95% observed) | 5,683 | 7,286 | 6,184 | 2,649 |
| candidate trajectory variables | 29 | 33 | 25 | 25 |

Most-missing graph nodes (fraction of people missing):

| ROSMAP | Diverse | Banner | BannerLFQ |
|---|---|---|---|
| SIGMAR1 0.46 | SREBF2 0.32 | GORASP1 0.73 | SEC23B 0.95 |
| AP4B1 0.42 | AP4M1 0.18 | SREBF2 0.64 | MAPK14 0.94 |
| AP4M1 0.42 | GORASP1 0.17 | AP4M1 0.58 | SEC16A 0.85 |
| AP3B2 0.40 | ABCA7 0.05 | ABCA7 0.40 | IDE 0.73 |
| SREBF2 0.36 | SIGMAR1 0.03 | AP4E1 0.36 | GGA1 0.65 |
| GORASP1 0.34 | AP4E1 0.01 | AP4B1 0.36 | COPG2 0.52 |
| ABCA7 0.30 | MAP2K3 0.01 | MAPKAPK2 0.22 | PREB 0.35 |
| RTN3 0.28 | PLCG2 0.01 | SIGMAR1 0.18 | ITPR2 0.34 |

Nodes absent from BannerLFQ (22): ABCA7, ADAM17, AP4B1, AP4E1, AP4M1, AP4S1, BACE1, ECE1, GGA2,
GORASP1, MAP2K3, MAPKAPK2, NCSTN, PLCG2, PLG, PSEN1, SEC24A, SEC24D, SIGMAR1, SREBF2, SYK, VLDLR.

Candidate trajectory variables (besides 5 cell-type scores and 10 proteome PCs):
- **ROSMAP (14):** Study, msex, educ, apoe_genotype, age_at_visit_max, age_first_ad_dx, age_death, cts_mmse30_first_ad_dx, cts_mmse30_lv, pmi, braaksc, ceradsc, cogdx, dcfdx_lv
- **Diverse (18):** dataContributionGroup, cohort, sex, race, isHispanic, ageDeath, PMI, apoeGenotype, amyThal, amyA, amyCerad, Braak, mayoDx, amyAny, bScore, reag, ADoutcome, derivedOutcomeBasedOnMayoDx
- **Banner, BannerLFQ (10):** sex, ageDeath, apoeGenotype, pmi, diagnosis, CERAD, Braak, PlaqueTotal, TangleTotal, lastMMSE

Cell-type score markers missing from a dataset (all other listed markers were used): ROSMAP microglia
lacked CD68, CSF1R, CX3CR1 (used AIF1, ITGAM, P2RY12, HLA-DRA); Diverse microglia lacked CD68; Banner
neuron lacked RBFOX3 and Banner microglia lacked CD68. In BannerLFQ the microglia score exists for only 32 people and the
endothelial score for 101.

---
## 2. The confidence graph C

| | value |
|---|---|
| file | `CHRONOS_80x80_confidence_v3.xlsx` (first sheet), raw value / 5, diagonal set to 0 |
| directed prior edges (C ≥ 0.05) | 1308 |
| unordered pairs with a prior edge | 654 of 3,160 |
| one-way pairs | 0 |
| two-way pairs | 654 (583 with equal C both ways) |
| C values among edges | 0.2: 403, 0.4: 230, 0.6: 399, 0.8: 123, 1: 153 |
| largest out-degree | APP 44, HSPA5 34, BACE1 34, TMED10 30, SEC13 30 |
| largest in-degree | APP 44, HSPA5 34, BACE1 34, TMED10 30, SEC13 30 |

Provenance of the directed prior edges (from the `Edge_List` sheet):

| provenance | edges |
|---|---|
| AUTO | 712 |
| LIT | 476 |
| COMPLEX | 70 |
| CHRONOS | 50 |

`Sign_80x80` among directed prior edges: +1 30, 0 1203, −1 75. The pipeline does not read this sheet.

Raw off-diagonal values in the matrix: -2: 574, -1: 778, 0: 3660, 1: 403, 2: 230, 3: 399, 4: 123, 5: 153. Values ≤ 0 are non-edges for the pipeline.

---

## 3. Signal: how well each dataset recovers C (AUROC, HIW 4.4)

0.5 = no relation to C. "Signed" (used by the pipeline) ranks pairs by ρ; "unsigned" by abs(ρ).

| | ROSMAP | Diverse | Banner | BannerLFQ |
|---|---|---|---|---|
| signed AUROC, raw | 0.5926 | 0.5767 | 0.5535 | 0.5283 |
| signed AUROC, adjusted for Z* | **0.6232** | **0.6441** | **0.5825** | **0.5369** |
| unsigned AUROC, raw | 0.5471 | 0.5862 | 0.5345 | 0.5213 |
| unsigned AUROC, adjusted | 0.5722 | 0.6120 | 0.5374 | 0.5238 |
| pairs scored (n ≥ 20) | 3,160 | 3,160 | 3,157 | 1,533 |
| prior pairs among them | 654 | 654 | 654 | 386 |

Correlation between prior-linked pairs vs other pairs:

| | ROSMAP raw | ROSMAP adj | Diverse raw | Diverse adj | Banner raw | Banner adj | BannerLFQ raw | BannerLFQ adj |
|---|---|---|---|---|---|---|---|---|
| prior pairs: mean ρ | +0.072 | +0.069 | +0.175 | +0.156 | +0.052 | +0.050 | +0.122 | +0.145 |
| prior pairs: fraction ρ > 0 | 0.641 | 0.680 | 0.746 | 0.804 | 0.587 | 0.598 | 0.668 | 0.687 |
| prior pairs: mean abs(ρ) | 0.138 | 0.122 | 0.278 | 0.190 | 0.152 | 0.130 | 0.268 | 0.276 |
| other pairs: mean ρ | +0.017 | +0.008 | +0.104 | +0.066 | +0.012 | -0.000 | +0.092 | +0.104 |
| other pairs: fraction ρ > 0 | 0.523 | 0.501 | 0.687 | 0.644 | 0.516 | 0.494 | 0.669 | 0.660 |
| other pairs: mean abs(ρ) | 0.114 | 0.088 | 0.221 | 0.130 | 0.130 | 0.109 | 0.258 | 0.263 |

---

## 4. Trajectory variables (HIW 4.6)

### 4.1 Selected Z*, in order

**ROSMAP** — start AUROC 0.5926, end 0.6232

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | proteomePC6 | proteome-derived | 0.5989 | 0.5927 | 0.0063 | 13.13 | 400 |
| 2 | proteomePC2 | proteome-derived | 0.6041 | 0.5989 | 0.0052 | 6.96 | 400 |
| 3 | celltype_neuron | proteome-derived | 0.6122 | 0.6040 | 0.0081 | 23.83 | 400 |
| 4 | celltype_oligodendrocyte | proteome-derived | 0.6186 | 0.6121 | 0.0064 | 10.32 | 400 |
| 5 | proteomePC8 | proteome-derived | 0.6232 | 0.6187 | 0.0047 | 6.21 | 400 |

**Diverse** — start AUROC 0.5767, end 0.6441

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | celltype_oligodendrocyte | proteome-derived | 0.6001 | 0.5767 | 0.0233 | 174.10 | 964 |
| 2 | reag | clinical/pathology | 0.6192 | 0.6139 | 0.0192 | 14.09 | 676 |
| 3 | proteomePC2 | proteome-derived | 0.6346 | 0.6193 | 0.0153 | 58.94 | 676 |
| 4 | proteomePC8 | proteome-derived | 0.6398 | 0.6346 | 0.0053 | 37.66 | 676 |
| 5 | proteomePC6 | proteome-derived | 0.6441 | 0.6398 | 0.0042 | 28.70 | 676 |

**Banner** — start AUROC 0.5535, end 0.5825

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | celltype_oligodendrocyte | proteome-derived | 0.5578 | 0.5529 | 0.0043 | 6.14 | 190 |
| 2 | proteomePC3 | proteome-derived | 0.5634 | 0.5571 | 0.0056 | 4.10 | 190 |
| 3 | celltype_neuron | proteome-derived | 0.5710 | 0.5632 | 0.0076 | 8.44 | 190 |
| 4 | proteomePC8 | proteome-derived | 0.5770 | 0.5704 | 0.0059 | 8.18 | 190 |
| 5 | apoeGenotype | clinical/pathology | 0.5825 | 0.5748 | 0.0055 | 3.38 | 189 |

**BannerLFQ** — start AUROC 0.5283, end 0.5369

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | celltype_neuron | proteome-derived | 0.5349 | 0.5280 | 0.0065 | 9.84 | 190 |
| 2 | proteomePC9 | proteome-derived | 0.5369 | 0.5346 | 0.0020 | 2.42 | 190 |

Proteome PC descriptions (share of rest-of-proteome variance; cell-type score it tracks most, Spearman ρ):

| PC | ROSMAP | Diverse | Banner |
|---|---|---|---|
| 1 | 9.9%, neuron +0.59 | 23.3%, neuron −0.89 | 12.2%, oligodendrocyte −0.59 |
| 2 | 6.1%, neuron +0.27 | 8.8%, astrocyte +0.67 | 6.9%, oligodendrocyte −0.35 |
| 3 | 4.9%, oligodendrocyte −0.47 | 4.0%, oligodendrocyte +0.25 | 5.6%, neuron +0.23 |
| 4 | 4.6%, oligodendrocyte −0.48 | 3.2%, endothelial −0.17 | 3.6%, endothelial −0.38 |
| 5 | 3.2%, astrocyte −0.47 | 3.0%, microglia −0.24 | 3.1%, neuron −0.31 |
| 6 | 2.9%, oligodendrocyte −0.13 | 2.7%, microglia +0.21 | 2.6%, astrocyte −0.20 |
| 7 | 2.4%, microglia −0.23 | 2.0%, oligodendrocyte −0.22 | 2.2%, microglia −0.30 |
| 8 | 2.3%, oligodendrocyte +0.29 | 1.9%, astrocyte +0.14 | 2.1%, oligodendrocyte −0.08 |
| 9 | 2.1%, astrocyte +0.27 | 1.8%, astrocyte −0.37 | 2.0%, endothelial −0.27 |
| 10 | 1.7%, astrocyte +0.32 | 1.5%, astrocyte +0.24 | 1.8%, astrocyte +0.39 |

BannerLFQ PC1: 24.1% of variance, tracks microglia (ρ = −0.42); PC9: 1.2%, microglia (ρ = −0.20).

### 4.2 Every candidate at step 1 (each variable alone vs no adjustment)

`gain` = AUROC change vs unadjusted; `z` = vs 10 shuffled copies; eligible = z > 2 and gain > 0.0001.

**ROSMAP** (unadjusted AUROC 0.5926)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| proteomePC6 | proteome-derived | 0.5989 | 0.5927 | 0.0063 | 13.13 | 400 |
| proteomePC2 | proteome-derived | 0.5987 | 0.5925 | 0.0061 | 10.32 | 400 |
| celltype_oligodendrocyte | proteome-derived | 0.5969 | 0.5928 | 0.0042 | 8.53 | 400 |
| proteomePC8 | proteome-derived | 0.5962 | 0.5926 | 0.0036 | 8.35 | 400 |
| celltype_neuron | proteome-derived | 0.5961 | 0.5928 | 0.0034 | 10.18 | 400 |
| proteomePC7 | proteome-derived | 0.5950 | 0.5923 | 0.0024 | 2.51 | 400 |
| apoe_genotype | clinical/pathology | 0.5941 | 0.5928 | 0.0014 | 1.73 | 400 |
| pmi | technical | 0.5933 | 0.5924 | 0.0007 | 2.13 | 398 |
| celltype_microglia | proteome-derived | 0.5933 | 0.5926 | 0.0007 | 0.82 | 400 |
| msex | clinical/pathology | 0.5932 | 0.5929 | 0.0006 | 0.44 | 400 |
| cts_mmse30_lv | clinical/pathology | 0.5928 | 0.5924 | 0.0002 | 0.60 | 400 |
| educ | clinical/pathology | 0.5926 | 0.5926 | -0.0001 | -0.07 | 400 |
| celltype_astrocyte | proteome-derived | 0.5925 | 0.5926 | -0.0002 | -0.47 | 400 |
| proteomePC10 | proteome-derived | 0.5924 | 0.5928 | -0.0003 | -0.93 | 400 |
| Study | technical | 0.5922 | 0.5926 | -0.0004 | -0.83 | 400 |
| braaksc | clinical/pathology | 0.5922 | 0.5926 | -0.0005 | -0.80 | 400 |
| proteomePC5 | proteome-derived | 0.5920 | 0.5926 | -0.0006 | -0.98 | 400 |
| ceradsc | clinical/pathology | 0.5916 | 0.5927 | -0.0010 | -2.04 | 400 |
| cogdx | clinical/pathology | 0.5915 | 0.5923 | -0.0011 | -1.33 | 400 |
| proteomePC9 | proteome-derived | 0.5915 | 0.5927 | -0.0011 | -2.69 | 400 |
| age_at_visit_max | clinical/pathology | 0.5914 | 0.5926 | -0.0013 | -2.24 | 400 |
| dcfdx_lv | clinical/pathology | 0.5911 | 0.5928 | -0.0016 | -3.39 | 400 |
| age_death | clinical/pathology | 0.5909 | 0.5924 | -0.0018 | -3.30 | 400 |
| celltype_endothelial | proteome-derived | 0.5906 | 0.5924 | -0.0021 | -3.55 | 400 |
| proteomePC3 | proteome-derived | 0.5863 | 0.5932 | -0.0064 | -9.88 | 400 |
| proteomePC1 | proteome-derived | 0.5856 | 0.5927 | -0.0071 | -14.06 | 400 |
| proteomePC4 | proteome-derived | 0.5796 | 0.5927 | -0.0130 | -17.78 | 400 |
| cts_mmse30_first_ad_dx | clinical/pathology | 0.5744 | 0.5728 | -0.0182 | 0.98 | 126 |
| age_first_ad_dx | clinical/pathology | 0.5708 | 0.5732 | -0.0218 | -2.41 | 132 |

**Diverse** (unadjusted AUROC 0.5767)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| celltype_oligodendrocyte | proteome-derived | 0.6001 | 0.5767 | 0.0233 | 174.10 | 964 |
| proteomePC1 | proteome-derived | 0.5954 | 0.5767 | 0.0187 | 162.11 | 964 |
| reag | clinical/pathology | 0.5944 | 0.5906 | 0.0177 | 12.36 | 676 |
| amyCerad | clinical/pathology | 0.5935 | 0.5905 | 0.0168 | 14.85 | 676 |
| amyAny | clinical/pathology | 0.5934 | 0.5906 | 0.0167 | 17.49 | 676 |
| proteomePC2 | proteome-derived | 0.5811 | 0.5768 | 0.0044 | 95.95 | 964 |
| PMI | technical | 0.5811 | 0.5811 | 0.0043 | -0.15 | 818 |
| Braak | clinical/pathology | 0.5808 | 0.5776 | 0.0041 | 12.41 | 908 |
| bScore | clinical/pathology | 0.5804 | 0.5777 | 0.0037 | 18.26 | 908 |
| cohort | technical | 0.5794 | 0.5767 | 0.0027 | 7.37 | 964 |
| ADoutcome | clinical/pathology | 0.5794 | 0.5767 | 0.0027 | 17.54 | 961 |
| proteomePC6 | proteome-derived | 0.5792 | 0.5767 | 0.0025 | 16.08 | 964 |
| celltype_astrocyte | proteome-derived | 0.5786 | 0.5767 | 0.0019 | 35.58 | 964 |
| proteomePC10 | proteome-derived | 0.5777 | 0.5768 | 0.0010 | 11.40 | 964 |
| celltype_endothelial | proteome-derived | 0.5776 | 0.5767 | 0.0009 | 5.41 | 964 |
| dataContributionGroup | technical | 0.5772 | 0.5767 | 0.0005 | 2.01 | 964 |
| isHispanic | clinical/pathology | 0.5770 | 0.5766 | 0.0003 | 5.08 | 963 |
| celltype_neuron | proteome-derived | 0.5768 | 0.5767 | 0.0000 | 0.33 | 964 |
| race | clinical/pathology | 0.5767 | 0.5766 | -0.0000 | 0.43 | 961 |
| sex | clinical/pathology | 0.5766 | 0.5767 | -0.0002 | -0.98 | 964 |
| derivedOutcomeBasedOnMayoDx | clinical/pathology | 0.5765 | 0.5767 | -0.0002 | -1.24 | 964 |
| ageDeath | clinical/pathology | 0.5764 | 0.5767 | -0.0003 | -1.60 | 963 |
| proteomePC9 | proteome-derived | 0.5763 | 0.5768 | -0.0004 | -5.20 | 964 |
| proteomePC8 | proteome-derived | 0.5760 | 0.5767 | -0.0008 | -5.72 | 964 |
| proteomePC4 | proteome-derived | 0.5760 | 0.5767 | -0.0008 | -7.88 | 964 |
| apoeGenotype | clinical/pathology | 0.5757 | 0.5756 | -0.0010 | 0.84 | 895 |
| proteomePC5 | proteome-derived | 0.5756 | 0.5767 | -0.0011 | -16.52 | 964 |
| proteomePC7 | proteome-derived | 0.5754 | 0.5767 | -0.0013 | -9.67 | 964 |
| proteomePC3 | proteome-derived | 0.5726 | 0.5768 | -0.0041 | -33.92 | 964 |
| celltype_microglia | proteome-derived | 0.5723 | 0.5768 | -0.0044 | -54.23 | 964 |
| amyThal | clinical/pathology | 0.5686 | 0.5679 | -0.0081 | 2.38 | 391 |
| amyA | clinical/pathology | 0.5685 | 0.5679 | -0.0083 | 2.41 | 398 |
| mayoDx | clinical/pathology | 0.5592 | 0.5573 | -0.0175 | 5.58 | 285 |

**Banner** (unadjusted AUROC 0.5535)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| celltype_oligodendrocyte | proteome-derived | 0.5578 | 0.5529 | 0.0043 | 6.14 | 190 |
| proteomePC8 | proteome-derived | 0.5572 | 0.5535 | 0.0037 | 5.50 | 190 |
| proteomePC3 | proteome-derived | 0.5571 | 0.5533 | 0.0036 | 9.79 | 190 |
| TangleTotal | clinical/pathology | 0.5556 | 0.5533 | 0.0021 | 2.33 | 190 |
| Braak | clinical/pathology | 0.5553 | 0.5530 | 0.0018 | 1.87 | 190 |
| celltype_neuron | proteome-derived | 0.5553 | 0.5531 | 0.0018 | 2.99 | 190 |
| proteomePC6 | proteome-derived | 0.5545 | 0.5531 | 0.0010 | 1.84 | 190 |
| sex | clinical/pathology | 0.5539 | 0.5527 | 0.0004 | 1.06 | 190 |
| lastMMSE | clinical/pathology | 0.5538 | 0.5531 | 0.0003 | 0.58 | 190 |
| celltype_astrocyte | proteome-derived | 0.5538 | 0.5533 | 0.0003 | 0.56 | 190 |
| PlaqueTotal | clinical/pathology | 0.5532 | 0.5531 | -0.0003 | 0.10 | 190 |
| apoeGenotype | clinical/pathology | 0.5531 | 0.5522 | -0.0004 | 1.04 | 189 |
| proteomePC7 | proteome-derived | 0.5530 | 0.5533 | -0.0005 | -0.44 | 190 |
| proteomePC2 | proteome-derived | 0.5525 | 0.5529 | -0.0010 | -0.51 | 190 |
| pmi | technical | 0.5524 | 0.5534 | -0.0011 | -1.21 | 190 |
| proteomePC9 | proteome-derived | 0.5522 | 0.5532 | -0.0013 | -1.21 | 190 |
| diagnosis | clinical/pathology | 0.5520 | 0.5532 | -0.0015 | -1.11 | 190 |
| ageDeath | clinical/pathology | 0.5517 | 0.5529 | -0.0018 | -1.30 | 190 |
| proteomePC5 | proteome-derived | 0.5515 | 0.5531 | -0.0020 | -2.37 | 190 |
| celltype_endothelial | proteome-derived | 0.5511 | 0.5528 | -0.0024 | -1.77 | 190 |
| proteomePC4 | proteome-derived | 0.5510 | 0.5532 | -0.0025 | -3.93 | 190 |
| celltype_microglia | proteome-derived | 0.5508 | 0.5530 | -0.0027 | -2.95 | 190 |
| CERAD | clinical/pathology | 0.5503 | 0.5532 | -0.0033 | -3.82 | 190 |
| proteomePC10 | proteome-derived | 0.5495 | 0.5532 | -0.0040 | -6.41 | 190 |
| proteomePC1 | proteome-derived | 0.5494 | 0.5531 | -0.0042 | -3.84 | 190 |

**BannerLFQ** (unadjusted AUROC 0.5283)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| celltype_neuron | proteome-derived | 0.5349 | 0.5280 | 0.0065 | 9.84 | 190 |
| celltype_oligodendrocyte | proteome-derived | 0.5303 | 0.5282 | 0.0020 | 1.76 | 190 |
| proteomePC9 | proteome-derived | 0.5294 | 0.5285 | 0.0010 | 0.81 | 190 |
| proteomePC8 | proteome-derived | 0.5290 | 0.5286 | 0.0007 | 0.28 | 190 |
| proteomePC5 | proteome-derived | 0.5287 | 0.5285 | 0.0004 | 0.23 | 190 |
| apoeGenotype | clinical/pathology | 0.5284 | 0.5284 | 0.0001 | -0.01 | 189 |
| celltype_astrocyte | proteome-derived | 0.5280 | 0.5282 | -0.0004 | -0.21 | 190 |
| pmi | technical | 0.5278 | 0.5279 | -0.0006 | -0.15 | 190 |
| Braak | clinical/pathology | 0.5277 | 0.5283 | -0.0006 | -0.59 | 190 |
| proteomePC6 | proteome-derived | 0.5274 | 0.5280 | -0.0009 | -0.88 | 190 |
| proteomePC7 | proteome-derived | 0.5273 | 0.5281 | -0.0011 | -1.52 | 190 |
| ageDeath | clinical/pathology | 0.5269 | 0.5277 | -0.0014 | -1.06 | 190 |
| proteomePC10 | proteome-derived | 0.5266 | 0.5275 | -0.0017 | -0.91 | 190 |
| proteomePC4 | proteome-derived | 0.5266 | 0.5282 | -0.0018 | -1.27 | 190 |
| CERAD | clinical/pathology | 0.5264 | 0.5279 | -0.0019 | -1.66 | 190 |
| PlaqueTotal | clinical/pathology | 0.5263 | 0.5284 | -0.0020 | -1.77 | 190 |
| TangleTotal | clinical/pathology | 0.5260 | 0.5284 | -0.0024 | -2.42 | 190 |
| sex | clinical/pathology | 0.5259 | 0.5280 | -0.0024 | -2.62 | 190 |
| diagnosis | clinical/pathology | 0.5258 | 0.5288 | -0.0025 | -3.06 | 190 |
| proteomePC2 | proteome-derived | 0.5256 | 0.5282 | -0.0028 | -1.84 | 190 |
| lastMMSE | clinical/pathology | 0.5255 | 0.5282 | -0.0029 | -3.38 | 190 |
| proteomePC3 | proteome-derived | 0.5209 | 0.5280 | -0.0074 | -8.37 | 190 |
| celltype_endothelial | proteome-derived | 0.5200 | 0.5187 | -0.0083 | 1.30 | 101 |
| proteomePC1 | proteome-derived | 0.5133 | 0.5282 | -0.0151 | -37.67 | 190 |
| celltype_microglia | proteome-derived | — | — | — | — | 32 |

### 4.3 Clinical/pathology variables at every step

z vs shuffled for each non-proteome variable at each selection step (the full per-step lists are in `{ds}_variables.csv`). Bold = selected at that step.

**ROSMAP**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| Study | -0.83 | -0.95 | -0.75 | -2.54 | -0.42 |
| age_at_visit_max | -2.24 | -4.33 | -0.38 | 0.15 | 0.47 |
| age_death | -3.30 | -3.15 | -3.10 | 0.40 | -0.89 |
| age_first_ad_dx | -2.41 | -1.92 | 0.18 | -0.14 | 0.20 |
| apoe_genotype | 1.73 | 0.85 | 3.43 | 2.56 | 2.39 |
| braaksc | -0.80 | -2.62 | -1.79 | -0.46 | 0.46 |
| ceradsc | -2.04 | -5.98 | -3.89 | -3.33 | -2.84 |
| cogdx | -1.33 | -0.69 | -0.31 | 0.49 | 1.10 |
| cts_mmse30_first_ad_dx | 0.98 | 0.99 | 0.09 | 0.20 | -0.04 |
| cts_mmse30_lv | 0.60 | 0.63 | 0.19 | 1.61 | 2.41 |
| dcfdx_lv | -3.39 | -2.86 | -2.08 | -1.43 | 0.54 |
| educ | -0.07 | -0.93 | 0.21 | -0.11 | -0.12 |
| msex | 0.44 | 1.38 | -0.53 | 1.37 | -0.90 |
| pmi | 2.13 | 0.65 | 3.53 | 0.71 | 1.47 |

**Diverse**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| ADoutcome | 17.54 | 24.76 | 1.48 | 1.64 | 0.40 |
| Braak | 12.41 | 19.94 | 1.72 | 1.88 | 1.90 |
| PMI | -0.15 | 0.35 | 0.16 | -0.86 | 1.39 |
| ageDeath | -1.60 | -1.92 | -1.23 | -1.61 | -2.30 |
| amyA | 2.41 | 0.19 | -0.95 | 0.53 | 0.69 |
| amyAny | 17.49 | 21.95 | 4.42 | 7.46 | 10.28 |
| amyCerad | 14.85 | 13.09 | 1.07 | 2.34 | 2.71 |
| amyThal | 2.38 | 0.25 | -2.13 | -0.46 | -0.76 |
| apoeGenotype | 0.84 | 4.04 | 0.68 | 0.32 | -0.57 |
| bScore | 18.26 | 22.08 | 3.80 | 2.59 | 3.95 |
| cohort | 7.37 | 5.47 | 2.89 | 3.93 | 3.65 |
| dataContributionGroup | 2.01 | 2.85 | 4.59 | 4.22 | 7.09 |
| derivedOutcomeBasedOnMayoDx | -1.24 | 1.38 | 0.06 | -0.22 | 0.34 |
| isHispanic | 5.08 | 2.96 | -0.18 | 0.74 | 0.18 |
| mayoDx | 5.58 | 9.61 | — | — | — |
| race | 0.43 | 5.38 | -0.78 | 0.21 | -1.06 |
| reag | 12.36 | **14.09** | — | — | — |
| sex | -0.98 | -0.61 | -0.51 | 0.36 | 1.02 |

**Banner**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| Braak | 1.87 | 0.93 | 4.65 | 2.08 | 3.00 |
| CERAD | -3.82 | -3.90 | -0.36 | -3.29 | -2.01 |
| PlaqueTotal | 0.10 | -0.09 | 5.53 | -0.54 | 0.49 |
| TangleTotal | 2.33 | 3.22 | 9.14 | 2.81 | 3.27 |
| ageDeath | -1.30 | -0.25 | 0.27 | -0.09 | -0.12 |
| apoeGenotype | 1.04 | 1.80 | 2.16 | 2.52 | **3.38** |
| diagnosis | -1.11 | -0.78 | 6.61 | 0.62 | 1.26 |
| lastMMSE | 0.58 | 1.61 | 11.26 | 1.48 | 5.93 |
| pmi | -1.21 | -0.03 | 0.68 | 0.57 | -1.43 |
| sex | 1.06 | 1.42 | 0.83 | 1.69 | 0.57 |

**BannerLFQ**

| variable | step 1 | step 2 | step 3 |
|---|---|---|---|
| Braak | -0.59 | -1.65 | 1.07 |
| CERAD | -1.66 | -2.79 | -2.86 |
| PlaqueTotal | -1.77 | -3.00 | -3.39 |
| TangleTotal | -2.42 | -1.79 | -1.06 |
| ageDeath | -1.06 | -2.26 | -2.36 |
| apoeGenotype | -0.01 | 0.02 | 0.32 |
| diagnosis | -3.06 | -3.73 | -3.13 |
| lastMMSE | -3.38 | -5.52 | -4.19 |
| pmi | -0.15 | -0.48 | -0.77 |
| sex | -2.62 | -3.06 | -1.17 |

---

## 5. Regions (HIW 4.7)

### 5.1 Overview

| dataset | regions_tested | depth1 | depth2 | z_above_2 | z_below_minus2 | STRONG_original | min_q_original | whole_cohort_AUROC | STRONG_corrected | min_q_corrected | z_above_2_corrected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ROSMAP | 301 | 77 | 224 | 28 | 2 | 0 | 0.264 | 0.623 | 0 | 0.125 | 31 |
| Diverse | 419 | 87 | 332 | 43 | 4 | 34 | 0.051 | 0.644 | 0 | 0.187 | 33 |
| Banner | 126 | 59 | 67 | 3 | 1 | 0 | 0.482 | 0.583 | 0 | 0.524 | 3 |
| BannerLFQ | 139 | 51 | 88 | 14 | 8 | 0 | 0.173 | 0.537 | 0 | 0.148 | 15 |

`original` = the pipeline's output (`{ds}_regions.csv`). `corrected` = same regions re-scored
with random subsets drawn only from people who have adjusted data and matched to the region's
usable size (`{ds}_regions_usable_null_check.csv`, HIW 9.1). People with adjusted data / people analysed:
ROSMAP 400/400, Diverse 676/964, Banner 189/190, BannerLFQ 190/190. Where the two numbers are equal, the two versions differ only by random draws.

### 5.2 ROSMAP

Whole-cohort adjusted AUROC: 0.6232

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_neuron > 0.0321 (upper half)  AND  celltype_endothelial <= -0.043 (lower half) | 2 | 124 | 0.6325 | 0.6041 | 2.59 | 0.2643 |
| celltype_neuron > 0.0321 (upper half)  AND  pmi 5.67-7.61 (middle third) | 2 | 69 | 0.6244 | 0.5908 | 2.58 | 0.2643 |
| celltype_neuron > 0.0321 (upper half)  AND  cts_mmse30_lv <= 26 (lower half) | 2 | 111 | 0.6309 | 0.6023 | 2.57 | 0.2643 |
| celltype_astrocyte > 0.235 (high third)  AND  ceradsc = 2 | 2 | 49 | 0.6234 | 0.5844 | 2.50 | 0.2643 |
| cts_mmse30_lv <= 26 (lower half)  AND  educ <= 14 (low third) | 2 | 85 | 0.6299 | 0.5955 | 2.50 | 0.2643 |
| celltype_neuron > 0.0321 (upper half)  AND  pmi > 6.5 (upper half) | 2 | 99 | 0.6260 | 0.5990 | 2.48 | 0.2643 |
| apoe_genotype = 23  AND  educ <= 16 (lower half) | 2 | 40 | 0.6177 | 0.5794 | 2.47 | 0.2643 |
| celltype_astrocyte > 0.235 (high third)  AND  educ <= 16 (lower half) | 2 | 76 | 0.6291 | 0.5955 | 2.44 | 0.2643 |
| celltype_neuron > 0.0321 (upper half) | 1 | 200 | 0.6313 | 0.6127 | 2.36 | 0.2643 |
| cts_mmse30_lv <= 26 (lower half)  AND  ceradsc = 2 | 2 | 84 | 0.6279 | 0.5955 | 2.35 | 0.2643 |
| celltype_neuron > 0.0321 (upper half)  AND  celltype_microglia > -0.077 (upper half) | 2 | 80 | 0.6274 | 0.5955 | 2.32 | 0.2643 |
| celltype_neuron > 0.0321 (upper half)  AND  educ <= 16 (lower half) | 2 | 123 | 0.6291 | 0.6041 | 2.28 | 0.2643 |
| cts_mmse30_lv <= 26 (lower half)  AND  cogdx = 2 | 2 | 53 | 0.6199 | 0.5844 | 2.27 | 0.2643 |
| apoe_genotype = 23 | 1 | 57 | 0.6209 | 0.5891 | 2.22 | 0.2643 |
| celltype_neuron > 0.0321 (upper half)  AND  braaksc <= 4 (lower half) | 2 | 171 | 0.6311 | 0.6109 | 2.21 | 0.2643 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| cts_mmse30_lv <= 26 (lower half)  AND  educ <= 14 (low third) | 85 | 85 | 0.6299 | 0.5940 | 2.88 | 0.1248 |
| celltype_astrocyte > 0.235 (high third)  AND  educ <= 16 (lower half) | 76 | 76 | 0.6291 | 0.5940 | 2.81 | 0.1248 |
| celltype_astrocyte > 0.235 (high third)  AND  ceradsc = 2 | 49 | 49 | 0.6234 | 0.5819 | 2.76 | 0.1248 |
| celltype_neuron > 0.0321 (upper half)  AND  celltype_endothelial <= -0.043 (lower half) | 124 | 124 | 0.6325 | 0.6032 | 2.75 | 0.1248 |
| cts_mmse30_lv <= 26 (lower half)  AND  ceradsc = 2 | 84 | 84 | 0.6279 | 0.5940 | 2.71 | 0.1248 |
| celltype_neuron > 0.0321 (upper half)  AND  celltype_microglia > -0.077 (upper half) | 80 | 80 | 0.6274 | 0.5940 | 2.68 | 0.1248 |
| celltype_neuron > 0.0321 (upper half)  AND  braaksc <= 4 (lower half) | 171 | 171 | 0.6311 | 0.6094 | 2.64 | 0.1576 |
| celltype_neuron > 0.0321 (upper half)  AND  cts_mmse30_lv <= 26 (lower half) | 111 | 111 | 0.6309 | 0.6019 | 2.63 | 0.1248 |
| apoe_genotype = 23  AND  educ <= 16 (lower half) | 40 | 40 | 0.6177 | 0.5781 | 2.60 | 0.1576 |
| celltype_neuron > 0.0321 (upper half)  AND  pmi 5.67-7.61 (middle third) | 69 | 69 | 0.6244 | 0.5912 | 2.59 | 0.1576 |
| cts_mmse30_lv <= 26 (lower half)  AND  cogdx = 2 | 53 | 53 | 0.6199 | 0.5819 | 2.53 | 0.1248 |
| celltype_neuron > 0.194 (high third)  AND  age_at_visit_max <= 85.7 (low third) | 51 | 51 | 0.6189 | 0.5819 | 2.46 | 0.1248 |
| celltype_neuron > 0.0321 (upper half)  AND  educ <= 16 (lower half) | 123 | 123 | 0.6291 | 0.6032 | 2.43 | 0.1248 |
| celltype_neuron > 0.0321 (upper half) | 200 | 200 | 0.6313 | 0.6135 | 2.42 | 0.1872 |
| cts_mmse30_lv <= 26 (lower half)  AND  celltype_astrocyte > 0.235 (high third) | 80 | 80 | 0.6240 | 0.5940 | 2.41 | 0.1248 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| ceradsc = 1 | 1 | 113 | 0.5762 | 0.6023 | -2.34 | 0.9851 |
| celltype_neuron <= 0.0321 (lower half) | 1 | 200 | 0.5963 | 0.6127 | -2.06 | 0.9734 |
| celltype_astrocyte -0.186-0.235 (middle third) | 1 | 133 | 0.5865 | 0.6066 | -1.97 | 0.9734 |
| cts_mmse30_lv <= 26 (lower half)  AND  dcfdx_lv = 4 | 2 | 102 | 0.5793 | 0.5990 | -1.81 | 0.9734 |
| cts_mmse30_lv <= 26 (lower half)  AND  cogdx = 4 | 2 | 106 | 0.5830 | 0.6023 | -1.73 | 0.9734 |
| celltype_neuron <= -0.149 (low third) | 1 | 134 | 0.5889 | 0.6066 | -1.73 | 0.9734 |
| dcfdx_lv = 4 | 1 | 104 | 0.5802 | 0.5990 | -1.72 | 0.9734 |
| cogdx = 4 | 1 | 109 | 0.5834 | 0.6023 | -1.70 | 0.9734 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| Study = MAP | 252 | 0.6160 | 0.6168 | -0.13 | 0.7299 |
| Study = ROS | 148 | 0.6039 | 0.6072 | -0.32 | 0.7395 |
| age_at_visit_max 85.7-90 (middle third) | 266 | 0.6212 | 0.6184 | 0.50 | 0.5592 |
| age_at_visit_max <= 85.7 (low third) | 134 | 0.5991 | 0.6066 | -0.73 | 0.8812 |
| age_at_visit_max <= 88.9 (lower half) | 200 | 0.6034 | 0.6127 | -1.17 | 0.9386 |
| age_at_visit_max > 88.9 (upper half) | 200 | 0.6232 | 0.6127 | 1.34 | 0.3907 |
| age_death 86.8-90 (middle third) | 266 | 0.6232 | 0.6184 | 0.86 | 0.4260 |
| age_death <= 86.8 (low third) | 134 | 0.5940 | 0.6066 | -1.24 | 0.9386 |
| age_death <= 89.4 (lower half) | 200 | 0.6046 | 0.6127 | -1.01 | 0.9304 |
| age_death > 89.4 (upper half) | 200 | 0.6209 | 0.6127 | 1.04 | 0.4260 |
| age_first_ad_dx 86.5-90 (middle third) | 88 | 0.5958 | 0.5959 | -0.01 | 0.7083 |
| age_first_ad_dx <= 86.5 (low third) | 44 | 0.5646 | 0.5794 | -0.95 | 0.9386 |
| age_first_ad_dx <= 89 (lower half) | 66 | 0.5732 | 0.5908 | -1.36 | 0.9463 |
| age_first_ad_dx > 89 (upper half) | 66 | 0.5969 | 0.5908 | 0.46 | 0.5788 |
| apoe_genotype = 23 | 57 | 0.6209 | 0.5891 | 2.22 | 0.2643 |
| apoe_genotype = 33 | 256 | 0.6129 | 0.6179 | -0.81 | 0.9019 |
| apoe_genotype = 34 | 73 | 0.5946 | 0.5908 | 0.29 | 0.6557 |
| braaksc 3-4 (middle third) | 131 | 0.6049 | 0.6066 | -0.16 | 0.7363 |
| braaksc <= 3 (low third) | 195 | 0.6158 | 0.6127 | 0.39 | 0.6232 |
| braaksc <= 4 (lower half) | 326 | 0.6226 | 0.6208 | 0.46 | 0.5745 |
| braaksc > 4 (high third) | 74 | 0.5889 | 0.5908 | -0.15 | 0.7299 |
| braaksc > 4 (upper half) | 74 | 0.5889 | 0.5908 | -0.15 | 0.7299 |
| ceradsc = 1 | 113 | 0.5762 | 0.6023 | -2.34 | 0.9851 |
| ceradsc = 2 | 138 | 0.6210 | 0.6074 | 1.35 | 0.3120 |
| ceradsc = 3 | 43 | 0.5770 | 0.5794 | -0.15 | 0.7255 |
| ceradsc = 4 | 106 | 0.6029 | 0.6023 | 0.06 | 0.6979 |
| cogdx = 1 | 168 | 0.6170 | 0.6109 | 0.68 | 0.5304 |
| cogdx = 2 | 97 | 0.6145 | 0.5990 | 1.42 | 0.2953 |
| cogdx = 4 | 109 | 0.5834 | 0.6023 | -1.70 | 0.9734 |
| cts_mmse30_first_ad_dx <= 18 (low third) | 48 | 0.5661 | 0.5844 | -1.17 | 0.9386 |
| cts_mmse30_first_ad_dx <= 21 (lower half) | 71 | 0.5929 | 0.5908 | 0.15 | 0.6697 |
| cts_mmse30_first_ad_dx > 21 (upper half) | 55 | 0.5724 | 0.5891 | -1.16 | 0.9304 |
| cts_mmse30_lv 23-28 (middle third) | 174 | 0.6196 | 0.6109 | 0.96 | 0.4260 |
| cts_mmse30_lv <= 23 (low third) | 136 | 0.5983 | 0.6074 | -0.91 | 0.9098 |
| cts_mmse30_lv <= 26 (lower half) | 210 | 0.6277 | 0.6132 | 1.96 | 0.2643 |
| cts_mmse30_lv > 26 (upper half) | 190 | 0.6016 | 0.6125 | -1.36 | 0.9425 |
| cts_mmse30_lv > 28 (high third) | 90 | 0.5890 | 0.5959 | -0.58 | 0.8346 |
| dcfdx_lv = 1 | 174 | 0.6134 | 0.6109 | 0.27 | 0.6557 |
| dcfdx_lv = 2 | 100 | 0.6114 | 0.5990 | 1.13 | 0.3964 |
| dcfdx_lv = 4 | 104 | 0.5802 | 0.5990 | -1.72 | 0.9734 |
| educ 14-18 (middle third) | 168 | 0.6076 | 0.6109 | -0.36 | 0.7488 |
| educ <= 14 (low third) | 158 | 0.6212 | 0.6102 | 1.20 | 0.3702 |
| educ <= 16 (lower half) | 247 | 0.6250 | 0.6168 | 1.35 | 0.3153 |
| educ > 16 (upper half) | 153 | 0.5929 | 0.6072 | -1.41 | 0.9734 |
| educ > 18 (high third) | 74 | 0.5698 | 0.5908 | -1.62 | 0.9726 |
| msex = 0 | 281 | 0.6266 | 0.6184 | 1.43 | 0.3387 |
| msex = 1 | 119 | 0.5910 | 0.6041 | -1.20 | 0.9386 |
| pmi 5.67-7.61 (middle third) | 131 | 0.6136 | 0.6066 | 0.69 | 0.5063 |
| pmi <= 5.67 (low third) | 134 | 0.6076 | 0.6066 | 0.10 | 0.6979 |
| pmi <= 6.5 (lower half) | 201 | 0.6124 | 0.6127 | -0.03 | 0.6979 |
| pmi > 6.5 (upper half) | 197 | 0.6161 | 0.6127 | 0.44 | 0.5918 |
| pmi > 7.61 (high third) | 133 | 0.6000 | 0.6066 | -0.65 | 0.8707 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.186-0.235 (middle third) | 133 | 0.5865 | 0.6066 | -1.97 | 0.9734 |
| celltype_astrocyte <= -0.186 (low third) | 134 | 0.6054 | 0.6066 | -0.11 | 0.7299 |
| celltype_astrocyte <= 0.0643 (lower half) | 200 | 0.6143 | 0.6127 | 0.20 | 0.6621 |
| celltype_astrocyte > 0.0643 (upper half) | 200 | 0.6142 | 0.6127 | 0.19 | 0.6647 |
| celltype_astrocyte > 0.235 (high third) | 133 | 0.6212 | 0.6066 | 1.44 | 0.2953 |
| celltype_endothelial -0.263-0.217 (middle third) | 133 | 0.6075 | 0.6066 | 0.09 | 0.6979 |
| celltype_endothelial <= -0.043 (lower half) | 200 | 0.6154 | 0.6127 | 0.34 | 0.6418 |
| celltype_endothelial <= -0.263 (low third) | 134 | 0.6051 | 0.6066 | -0.15 | 0.7299 |
| celltype_endothelial > -0.043 (upper half) | 200 | 0.6148 | 0.6127 | 0.27 | 0.6557 |
| celltype_endothelial > 0.217 (high third) | 133 | 0.6051 | 0.6066 | -0.15 | 0.7299 |
| celltype_microglia -0.284-0.197 (middle third) | 133 | 0.6066 | 0.6066 | 0.00 | 0.6979 |
| celltype_microglia <= -0.077 (lower half) | 200 | 0.6089 | 0.6127 | -0.48 | 0.8174 |
| celltype_microglia <= -0.284 (low third) | 134 | 0.6079 | 0.6066 | 0.13 | 0.6979 |
| celltype_microglia > -0.077 (upper half) | 200 | 0.6211 | 0.6127 | 1.07 | 0.4260 |
| celltype_microglia > 0.197 (high third) | 133 | 0.6068 | 0.6066 | 0.03 | 0.6979 |
| celltype_neuron -0.149-0.194 (middle third) | 133 | 0.6074 | 0.6066 | 0.09 | 0.6979 |
| celltype_neuron <= -0.149 (low third) | 134 | 0.5889 | 0.6066 | -1.73 | 0.9734 |
| celltype_neuron <= 0.0321 (lower half) | 200 | 0.5963 | 0.6127 | -2.06 | 0.9734 |
| celltype_neuron > 0.0321 (upper half) | 200 | 0.6313 | 0.6127 | 2.36 | 0.2643 |
| celltype_neuron > 0.194 (high third) | 133 | 0.6253 | 0.6066 | 1.84 | 0.2643 |
| celltype_oligodendrocyte -0.371-0.164 (middle third) | 133 | 0.6154 | 0.6066 | 0.87 | 0.4433 |
| celltype_oligodendrocyte <= -0.127 (lower half) | 200 | 0.6194 | 0.6127 | 0.86 | 0.4919 |
| celltype_oligodendrocyte <= -0.371 (low third) | 134 | 0.5971 | 0.6066 | -0.93 | 0.9198 |
| celltype_oligodendrocyte > -0.127 (upper half) | 200 | 0.6078 | 0.6127 | -0.61 | 0.8812 |
| celltype_oligodendrocyte > 0.164 (high third) | 133 | 0.6030 | 0.6066 | -0.35 | 0.7818 |

### 5.3 Diverse

Whole-cohort adjusted AUROC: 0.6441

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| PMI <= 6.83 (lower half)  AND  race = White | 2 | 256 | 0.6615 | 0.6335 | 3.35 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  celltype_microglia <= -0.289 (low third) | 2 | 139 | 0.6633 | 0.6259 | 3.33 | 0.0508 |
| PMI <= 5.17 (low third)  AND  race = White | 2 | 172 | 0.6588 | 0.6284 | 3.12 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  celltype_oligodendrocyte <= -0.11 (lower half) | 2 | 162 | 0.6563 | 0.6275 | 2.80 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  celltype_astrocyte <= -0.213 (low third) | 2 | 108 | 0.6550 | 0.6209 | 2.73 | 0.0508 |
| PMI <= 5.17 (low third)  AND  celltype_astrocyte <= 0.0206 (lower half) | 2 | 116 | 0.6537 | 0.6223 | 2.69 | 0.0719 |
| race = White  AND  sex = female | 2 | 279 | 0.6556 | 0.6348 | 2.67 | 0.0719 |
| celltype_neuron <= 0.103 (lower half)  AND  reag = Low Likelihood | 2 | 64 | 0.6533 | 0.6099 | 2.64 | 0.0893 |
| PMI <= 6.83 (lower half)  AND  celltype_astrocyte <= 0.0206 (lower half) | 2 | 177 | 0.6546 | 0.6296 | 2.61 | 0.0893 |
| celltype_neuron <= 0.103 (lower half)  AND  ADoutcome = Control | 2 | 66 | 0.6471 | 0.6116 | 2.59 | 0.0508 |
| PMI <= 5.17 (low third)  AND  celltype_microglia <= -0.289 (low third) | 2 | 95 | 0.6549 | 0.6203 | 2.55 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  celltype_microglia <= -0.0479 (lower half) | 2 | 202 | 0.6550 | 0.6322 | 2.51 | 0.0893 |
| PMI <= 5.17 (low third)  AND  celltype_astrocyte <= -0.213 (low third) | 2 | 71 | 0.6456 | 0.6116 | 2.48 | 0.0508 |
| PMI <= 5.17 (low third)  AND  celltype_endothelial -0.298-0.239 (middle third) | 2 | 102 | 0.6535 | 0.6203 | 2.45 | 0.0508 |
| PMI <= 5.17 (low third)  AND  sex = female | 2 | 160 | 0.6525 | 0.6275 | 2.43 | 0.0719 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| PMI <= 6.83 (lower half)  AND  race = White | 256 | 236 | 0.6615 | 0.6382 | 3.41 | 0.1866 |
| PMI <= 6.83 (lower half)  AND  celltype_microglia <= -0.289 (low third) | 139 | 109 | 0.6633 | 0.6269 | 3.34 | 0.1964 |
| PMI <= 5.17 (low third)  AND  race = White | 172 | 154 | 0.6588 | 0.6322 | 3.30 | 0.1866 |
| PMI <= 6.83 (lower half)  AND  cohort = Emory | 40 | 40 | 0.6558 | 0.6085 | 3.08 | 0.1866 |
| PMI <= 5.17 (low third)  AND  celltype_microglia <= -0.289 (low third) | 95 | 71 | 0.6549 | 0.6176 | 2.79 | 0.1866 |
| PMI <= 6.83 (lower half)  AND  celltype_astrocyte <= 0.0206 (lower half) | 177 | 150 | 0.6546 | 0.6322 | 2.78 | 0.1866 |
| race = White  AND  sex = female | 279 | 251 | 0.6556 | 0.6384 | 2.74 | 0.1964 |
| celltype_neuron <= -0.0632 (low third)  AND  celltype_microglia <= -0.0479 (lower half) | 77 | 45 | 0.6496 | 0.6085 | 2.67 | 0.1964 |
| celltype_neuron <= 0.103 (lower half)  AND  celltype_oligodendrocyte <= -0.11 (lower half) | 162 | 129 | 0.6563 | 0.6304 | 2.66 | 0.1866 |
| PMI <= 6.83 (lower half)  AND  celltype_astrocyte <= -0.213 (low third) | 108 | 87 | 0.6550 | 0.6251 | 2.65 | 0.1866 |
| PMI <= 5.17 (low third)  AND  celltype_astrocyte <= 0.0206 (lower half) | 116 | 94 | 0.6537 | 0.6251 | 2.54 | 0.1866 |
| PMI <= 5.17 (low third)  AND  celltype_endothelial -0.298-0.239 (middle third) | 102 | 84 | 0.6535 | 0.6230 | 2.52 | 0.1964 |
| celltype_neuron <= 0.103 (lower half)  AND  reag = Low Likelihood | 64 | 64 | 0.6533 | 0.6157 | 2.50 | 0.1964 |
| celltype_neuron <= 0.103 (lower half)  AND  race = White | 224 | 198 | 0.6545 | 0.6357 | 2.46 | 0.1964 |
| PMI <= 5.17 (low third)  AND  isHispanic = FALSE | 228 | 201 | 0.6536 | 0.6357 | 2.33 | 0.2239 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| race = other | 1 | 214 | 0.5767 | 0.6328 | -6.26 | 1.0000 |
| amyThal = Phase 4 | 1 | 52 | 0.5599 | 0.6049 | -2.64 | 0.9924 |
| amyA = Thal Phase 4 or 5 | 1 | 273 | 0.6160 | 0.6353 | -2.44 | 0.9924 |
| celltype_neuron > 0.103 (upper half) | 1 | 482 | 0.6297 | 0.6410 | -2.36 | 0.9924 |
| Braak = Stage VI | 1 | 239 | 0.6182 | 0.6335 | -1.97 | 0.9924 |
| isHispanic = TRUE | 1 | 256 | 0.6176 | 0.6335 | -1.91 | 0.9924 |
| ageDeath <= 83.6 (lower half) | 1 | 482 | 0.6323 | 0.6410 | -1.81 | 0.9842 |
| celltype_endothelial <= -0.298 (low third) | 1 | 322 | 0.6247 | 0.6370 | -1.70 | 0.9586 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| ADoutcome = AD | 560 | 0.6368 | 0.6414 | -1.03 | 0.9002 |
| ADoutcome = Control | 209 | 0.6359 | 0.6328 | 0.35 | 0.5295 |
| ADoutcome = Other | 192 | 0.6304 | 0.6319 | -0.17 | 0.6455 |
| Braak = Stage I | 68 | 0.6158 | 0.6116 | 0.31 | 0.5343 |
| Braak = Stage II | 86 | 0.6090 | 0.6175 | -0.58 | 0.7928 |
| Braak = Stage III | 160 | 0.6282 | 0.6275 | 0.07 | 0.5947 |
| Braak = Stage IV | 164 | 0.6370 | 0.6275 | 0.92 | 0.3522 |
| Braak = Stage V | 191 | 0.6263 | 0.6319 | -0.61 | 0.8101 |
| Braak = Stage VI | 239 | 0.6182 | 0.6335 | -1.97 | 0.9924 |
| PMI 5.17-10 (middle third) | 275 | 0.6400 | 0.6348 | 0.67 | 0.4304 |
| PMI <= 5.17 (low third) | 274 | 0.6526 | 0.6353 | 2.19 | 0.1405 |
| PMI <= 6.83 (lower half) | 409 | 0.6513 | 0.6391 | 2.20 | 0.0719 |
| PMI > 10 (high third) | 269 | 0.6242 | 0.6353 | -1.40 | 0.9586 |
| PMI > 6.83 (upper half) | 409 | 0.6314 | 0.6391 | -1.40 | 0.9424 |
| ageDeath 79-87.7 (middle third) | 303 | 0.6344 | 0.6367 | -0.30 | 0.7404 |
| ageDeath <= 79 (low third) | 339 | 0.6346 | 0.6379 | -0.48 | 0.7712 |
| ageDeath <= 83.6 (lower half) | 482 | 0.6323 | 0.6410 | -1.81 | 0.9842 |
| ageDeath > 83.6 (upper half) | 481 | 0.6463 | 0.6410 | 1.10 | 0.3270 |
| ageDeath > 87.7 (high third) | 321 | 0.6368 | 0.6370 | -0.02 | 0.6313 |
| amyA = Thal Phase 1 or 2 | 64 | 0.5996 | 0.6099 | -0.63 | 0.8001 |
| amyA = Thal Phase 3 | 61 | 0.6336 | 0.6099 | 1.44 | 0.2358 |
| amyA = Thal Phase 4 or 5 | 273 | 0.6160 | 0.6353 | -2.44 | 0.9924 |
| amyAny = 0 | 167 | 0.6365 | 0.6284 | 0.83 | 0.4012 |
| amyAny = 1 | 509 | 0.6415 | 0.6413 | 0.03 | 0.6120 |
| amyCerad = Frequent/Definite/C3 | 326 | 0.6326 | 0.6368 | -0.58 | 0.7941 |
| amyCerad = Moderate/Probable/C2 | 138 | 0.6396 | 0.6259 | 1.22 | 0.2520 |
| amyCerad = None/No AD/C0 | 167 | 0.6365 | 0.6284 | 0.83 | 0.4012 |
| amyCerad = Sparse/Possible/C1 | 45 | 0.6133 | — | — | 0.3840 |
| amyThal = Phase 3 | 61 | 0.6336 | 0.6099 | 1.44 | 0.2358 |
| amyThal = Phase 4 | 52 | 0.5599 | 0.6049 | -2.64 | 0.9924 |
| amyThal = Phase 5 | 217 | 0.6252 | 0.6327 | -0.81 | 0.8699 |
| apoeGenotype = 23 | 73 | 0.6138 | 0.6116 | 0.16 | 0.5770 |
| apoeGenotype = 33 | 449 | 0.6430 | 0.6397 | 0.64 | 0.4624 |
| apoeGenotype = 34 | 286 | 0.6360 | 0.6365 | -0.07 | 0.6536 |
| apoeGenotype = 44 | 58 | 0.5922 | 0.6099 | -1.08 | 0.9091 |
| bScore = Braak Stage I-II | 154 | 0.6221 | 0.6282 | -0.53 | 0.7941 |
| bScore = Braak Stage III-IV | 324 | 0.6401 | 0.6370 | 0.43 | 0.5143 |
| bScore = Braak Stage V-VI | 430 | 0.6315 | 0.6394 | -1.36 | 0.9419 |
| cohort = CLINCOR | 65 | 0.6209 | 0.6099 | 0.66 | 0.4169 |
| cohort = Emory | 97 | 0.6411 | 0.6203 | 1.54 | 0.1978 |
| cohort = MAP | 69 | 0.6262 | 0.6116 | 1.07 | 0.3038 |
| cohort = Mayo Clinic | 248 | — | 0.6343 | — | 0.0508 |
| cohort = Mt Sinai Brain Bank | 205 | 0.6177 | 0.6322 | -1.59 | 0.9613 |
| cohort = ROS | 188 | 0.6363 | 0.6319 | 0.49 | 0.4777 |
| dataContributionGroup = Emory | 141 | 0.6296 | 0.6259 | 0.33 | 0.5144 |
| dataContributionGroup = MSSM | 183 | 0.6166 | 0.6296 | -1.36 | 0.9480 |
| dataContributionGroup = Mayo | 285 | — | 0.6348 | — | 0.0508 |
| dataContributionGroup = Rush | 355 | 0.6431 | 0.6385 | 0.77 | 0.4169 |
| derivedOutcomeBasedOnMayoDx = FALSE | 679 | 0.6441 | 0.6426 | 0.44 | 0.5108 |
| derivedOutcomeBasedOnMayoDx = TRUE | 285 | — | 0.6348 | — | 0.0508 |
| isHispanic = FALSE | 707 | 0.6445 | 0.6427 | 0.55 | 0.4624 |
| isHispanic = TRUE | 256 | 0.6176 | 0.6335 | -1.91 | 0.9924 |
| mayoDx = AD | 181 | — | 0.6296 | — | 0.0508 |
| mayoDx = Other | 78 | — | 0.6145 | — | 0.0508 |
| race = Black or African American | 267 | 0.6270 | 0.6353 | -1.05 | 0.9259 |
| race = White | 465 | 0.6497 | 0.6397 | 1.92 | 0.1588 |
| race = other | 214 | 0.5767 | 0.6328 | -6.26 | 1.0000 |
| reag = High Likelihood | 244 | 0.6276 | 0.6335 | -0.76 | 0.8349 |
| reag = Intermediate Likelihood | 209 | 0.6384 | 0.6328 | 0.63 | 0.4624 |
| reag = Low Likelihood | 191 | 0.6353 | 0.6319 | 0.38 | 0.5003 |
| sex = female | 567 | 0.6445 | 0.6416 | 0.70 | 0.4371 |
| sex = male | 397 | 0.6347 | 0.6388 | -0.71 | 0.8089 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.213-0.236 (middle third) | 321 | 0.6339 | 0.6370 | -0.44 | 0.7683 |
| celltype_astrocyte <= -0.213 (low third) | 322 | 0.6449 | 0.6370 | 1.09 | 0.2857 |
| celltype_astrocyte <= 0.0206 (lower half) | 482 | 0.6414 | 0.6410 | 0.09 | 0.5651 |
| celltype_astrocyte > 0.0206 (upper half) | 482 | 0.6399 | 0.6410 | -0.22 | 0.7026 |
| celltype_astrocyte > 0.236 (high third) | 321 | 0.6364 | 0.6370 | -0.08 | 0.6577 |
| celltype_endothelial -0.298-0.239 (middle third) | 321 | 0.6471 | 0.6370 | 1.40 | 0.1978 |
| celltype_endothelial <= -0.0425 (lower half) | 482 | 0.6348 | 0.6410 | -1.30 | 0.9230 |
| celltype_endothelial <= -0.298 (low third) | 322 | 0.6247 | 0.6370 | -1.70 | 0.9586 |
| celltype_endothelial > -0.0425 (upper half) | 482 | 0.6424 | 0.6410 | 0.29 | 0.5144 |
| celltype_endothelial > 0.239 (high third) | 321 | 0.6345 | 0.6370 | -0.35 | 0.7404 |
| celltype_microglia -0.289-0.179 (middle third) | 321 | 0.6357 | 0.6370 | -0.18 | 0.6929 |
| celltype_microglia <= -0.0479 (lower half) | 482 | 0.6424 | 0.6410 | 0.30 | 0.5144 |
| celltype_microglia <= -0.289 (low third) | 322 | 0.6469 | 0.6370 | 1.37 | 0.2070 |
| celltype_microglia > -0.0479 (upper half) | 482 | 0.6369 | 0.6410 | -0.85 | 0.8349 |
| celltype_microglia > 0.179 (high third) | 321 | 0.6273 | 0.6370 | -1.34 | 0.9424 |
| celltype_neuron -0.0632-0.277 (middle third) | 321 | 0.6324 | 0.6370 | -0.63 | 0.8199 |
| celltype_neuron <= -0.0632 (low third) | 322 | 0.6496 | 0.6370 | 1.75 | 0.1489 |
| celltype_neuron <= 0.103 (lower half) | 482 | 0.6517 | 0.6410 | 2.22 | 0.1097 |
| celltype_neuron > 0.103 (upper half) | 482 | 0.6297 | 0.6410 | -2.36 | 0.9924 |
| celltype_neuron > 0.277 (high third) | 321 | 0.6280 | 0.6370 | -1.25 | 0.9391 |
| celltype_oligodendrocyte -0.482-0.323 (middle third) | 321 | 0.6464 | 0.6370 | 1.30 | 0.2331 |
| celltype_oligodendrocyte <= -0.11 (lower half) | 482 | 0.6395 | 0.6410 | -0.31 | 0.7330 |
| celltype_oligodendrocyte <= -0.482 (low third) | 322 | 0.6331 | 0.6370 | -0.55 | 0.8001 |
| celltype_oligodendrocyte > -0.11 (upper half) | 482 | 0.6402 | 0.6410 | -0.17 | 0.6890 |
| celltype_oligodendrocyte > 0.323 (high third) | 321 | 0.6329 | 0.6370 | -0.57 | 0.8101 |

**The 34 regions flagged STRONG in the original run** (0 remain STRONG with the corrected null)

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| PMI <= 6.83 (lower half)  AND  race = White | 2 | 256 | 0.6615 | 0.6335 | 3.35 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  celltype_microglia <= -0.289 (low third) | 2 | 139 | 0.6633 | 0.6259 | 3.33 | 0.0508 |
| PMI <= 5.17 (low third)  AND  race = White | 2 | 172 | 0.6588 | 0.6284 | 3.12 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  celltype_oligodendrocyte <= -0.11 (lower half) | 2 | 162 | 0.6563 | 0.6275 | 2.80 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  celltype_astrocyte <= -0.213 (low third) | 2 | 108 | 0.6550 | 0.6209 | 2.73 | 0.0508 |
| PMI <= 5.17 (low third)  AND  celltype_astrocyte <= 0.0206 (lower half) | 2 | 116 | 0.6537 | 0.6223 | 2.69 | 0.0719 |
| race = White  AND  sex = female | 2 | 279 | 0.6556 | 0.6348 | 2.67 | 0.0719 |
| celltype_neuron <= 0.103 (lower half)  AND  reag = Low Likelihood | 2 | 64 | 0.6533 | 0.6099 | 2.64 | 0.0893 |
| PMI <= 6.83 (lower half)  AND  celltype_astrocyte <= 0.0206 (lower half) | 2 | 177 | 0.6546 | 0.6296 | 2.61 | 0.0893 |
| celltype_neuron <= 0.103 (lower half)  AND  ADoutcome = Control | 2 | 66 | 0.6471 | 0.6116 | 2.59 | 0.0508 |
| PMI <= 5.17 (low third)  AND  celltype_microglia <= -0.289 (low third) | 2 | 95 | 0.6549 | 0.6203 | 2.55 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  celltype_microglia <= -0.0479 (lower half) | 2 | 202 | 0.6550 | 0.6322 | 2.51 | 0.0893 |
| PMI <= 5.17 (low third)  AND  celltype_astrocyte <= -0.213 (low third) | 2 | 71 | 0.6456 | 0.6116 | 2.48 | 0.0508 |
| PMI <= 5.17 (low third)  AND  celltype_endothelial -0.298-0.239 (middle third) | 2 | 102 | 0.6535 | 0.6203 | 2.45 | 0.0508 |
| PMI <= 5.17 (low third)  AND  sex = female | 2 | 160 | 0.6525 | 0.6275 | 2.43 | 0.0719 |
| celltype_neuron <= -0.0632 (low third)  AND  celltype_microglia <= -0.0479 (lower half) | 2 | 77 | 0.6496 | 0.6145 | 2.41 | 0.0508 |
| celltype_neuron <= -0.0632 (low third)  AND  dataContributionGroup = Rush | 2 | 110 | 0.6506 | 0.6209 | 2.38 | 0.0719 |
| celltype_neuron <= 0.103 (lower half)  AND  race = White | 2 | 224 | 0.6545 | 0.6327 | 2.37 | 0.0719 |
| celltype_neuron <= 0.103 (lower half)  AND  amyAny = 0 | 2 | 51 | 0.6449 | 0.6049 | 2.35 | 0.0719 |
| celltype_neuron <= 0.103 (lower half)  AND  amyCerad = None/No AD/C0 | 2 | 51 | 0.6449 | 0.6049 | 2.35 | 0.0719 |
| PMI <= 5.17 (low third)  AND  isHispanic = FALSE | 2 | 228 | 0.6536 | 0.6333 | 2.30 | 0.0893 |
| PMI <= 6.83 (lower half) | 1 | 409 | 0.6513 | 0.6391 | 2.20 | 0.0719 |
| PMI <= 5.17 (low third)  AND  bScore = Braak Stage III-IV | 2 | 99 | 0.6500 | 0.6203 | 2.19 | 0.0719 |
| PMI <= 5.17 (low third)  AND  derivedOutcomeBasedOnMayoDx = FALSE | 2 | 220 | 0.6526 | 0.6327 | 2.16 | 0.0893 |
| race = White  AND  celltype_microglia <= -0.289 (low third) | 2 | 156 | 0.6496 | 0.6275 | 2.15 | 0.0719 |
| PMI <= 6.83 (lower half)  AND  amyCerad = None/No AD/C0 | 2 | 68 | 0.6410 | 0.6116 | 2.15 | 0.0719 |
| PMI <= 6.83 (lower half)  AND  amyAny = 0 | 2 | 68 | 0.6410 | 0.6116 | 2.15 | 0.0719 |
| PMI <= 5.17 (low third)  AND  reag = Intermediate Likelihood | 2 | 82 | 0.6457 | 0.6145 | 2.14 | 0.0719 |
| PMI <= 6.83 (lower half)  AND  reag = Low Likelihood | 2 | 90 | 0.6489 | 0.6175 | 2.13 | 0.0719 |
| PMI <= 5.17 (low third)  AND  apoeGenotype = 34 | 2 | 80 | 0.6449 | 0.6145 | 2.08 | 0.0719 |
| PMI <= 5.17 (low third)  AND  celltype_microglia <= -0.0479 (lower half) | 2 | 142 | 0.6491 | 0.6259 | 2.07 | 0.0719 |
| celltype_neuron <= 0.103 (lower half)  AND  sex = female | 2 | 295 | 0.6524 | 0.6367 | 2.04 | 0.0893 |
| PMI <= 6.83 (lower half)  AND  ageDeath <= 79 (low third) | 2 | 109 | 0.6461 | 0.6209 | 2.02 | 0.0893 |
| race = White  AND  celltype_neuron <= -0.0632 (low third) | 2 | 129 | 0.6478 | 0.6234 | 2.01 | 0.0893 |

**Regions with no usable people after adjustment** (AUROC NaN, p set to the 1/201 floor; HIW 9.2)

| region | n | auroc | null_mean | p_vs_random_same_size | q_fdr |
|---|---|---|---|---|---|
| dataContributionGroup = Mayo | 285 | — | 0.6348 | 0.0050 | 0.0508 |
| cohort = Mayo Clinic | 248 | — | 0.6343 | 0.0050 | 0.0508 |
| mayoDx = AD | 181 | — | 0.6296 | 0.0050 | 0.0508 |
| mayoDx = Other | 78 | — | 0.6145 | 0.0050 | 0.0508 |
| derivedOutcomeBasedOnMayoDx = TRUE | 285 | — | 0.6348 | 0.0050 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  dataContributionGroup = Mayo | 147 | — | 0.6282 | 0.0050 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  cohort = Mayo Clinic | 129 | — | 0.6234 | 0.0050 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  race = other | 111 | — | 0.6209 | 0.0050 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  mayoDx = AD | 108 | — | 0.6209 | 0.0050 | 0.0508 |
| celltype_neuron <= 0.103 (lower half)  AND  derivedOutcomeBasedOnMayoDx = TRUE | 147 | — | 0.6282 | 0.0050 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  dataContributionGroup = Mayo | 64 | — | 0.6099 | 0.0050 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  cohort = Mayo Clinic | 49 | — | 0.6049 | 0.0050 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  race = other | 49 | — | 0.6049 | 0.0050 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  mayoDx = AD | 43 | — | — | 0.0050 | 0.0508 |
| PMI <= 6.83 (lower half)  AND  derivedOutcomeBasedOnMayoDx = TRUE | 64 | — | 0.6099 | 0.0050 | 0.0508 |
| PMI <= 5.17 (low third)  AND  dataContributionGroup = Mayo | 54 | — | 0.6049 | 0.0050 | 0.0508 |
| PMI <= 5.17 (low third)  AND  cohort = Mayo Clinic | 40 | — | — | 0.0050 | 0.0508 |
| PMI <= 5.17 (low third)  AND  isHispanic = TRUE | 46 | — | 0.6049 | 0.0050 | 0.0508 |
| PMI <= 5.17 (low third)  AND  derivedOutcomeBasedOnMayoDx = TRUE | 54 | — | 0.6049 | 0.0050 | 0.0508 |
| race = White  AND  dataContributionGroup = Mayo | 62 | — | 0.6099 | 0.0050 | 0.0508 |
| race = White  AND  cohort = Mayo Clinic | 62 | — | 0.6099 | 0.0050 | 0.0508 |
| race = White  AND  mayoDx = AD | 43 | — | — | 0.0050 | 0.0508 |
| race = White  AND  derivedOutcomeBasedOnMayoDx = TRUE | 62 | — | 0.6099 | 0.0050 | 0.0508 |
| celltype_neuron <= -0.0632 (low third)  AND  dataContributionGroup = Mayo | 123 | — | 0.6223 | 0.0050 | 0.0508 |
| celltype_neuron <= -0.0632 (low third)  AND  cohort = Mayo Clinic | 111 | — | 0.6209 | 0.0050 | 0.0508 |
| celltype_neuron <= -0.0632 (low third)  AND  race = other | 93 | — | 0.6175 | 0.0050 | 0.0508 |
| celltype_neuron <= -0.0632 (low third)  AND  mayoDx = AD | 92 | — | 0.6175 | 0.0050 | 0.0508 |
| celltype_neuron <= -0.0632 (low third)  AND  derivedOutcomeBasedOnMayoDx = TRUE | 123 | — | 0.6223 | 0.0050 | 0.0508 |

`Diverse_best_region_W.csv` was refitted in `PMI <= 6.83 (lower half)  AND  race = White`, the original top STRONG region.

### 5.4 Banner

Whole-cohort adjusted AUROC: 0.5825

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| PlaqueTotal > 13.5 (high third)  AND  CERAD = 3 | 2 | 50 | 0.5980 | 0.5676 | 2.39 | 0.4822 |
| celltype_microglia > 0.00693 (upper half)  AND  apoeGenotype = e3-3 | 2 | 44 | 0.5973 | 0.5648 | 2.30 | 0.4822 |
| apoeGenotype = e3-3  AND  celltype_neuron <= -1.75e-05 (lower half) | 2 | 41 | 0.5966 | 0.5648 | 2.24 | 0.4822 |
| celltype_microglia > 0.00693 (upper half)  AND  diagnosis = control | 2 | 45 | 0.5912 | 0.5648 | 1.87 | 0.4822 |
| celltype_endothelial > 0.125 (high third) | 1 | 63 | 0.5876 | 0.5687 | 1.81 | 0.4822 |
| celltype_microglia > 0.00693 (upper half)  AND  ageDeath > 86 (upper half) | 2 | 42 | 0.5899 | 0.5648 | 1.77 | 0.4822 |
| celltype_microglia > 0.00693 (upper half)  AND  TangleTotal <= 8 (lower half) | 2 | 40 | 0.5899 | 0.5648 | 1.77 | 0.4822 |
| apoeGenotype = e3-3  AND  CERAD = 3 | 2 | 44 | 0.5894 | 0.5648 | 1.73 | 0.4822 |
| celltype_microglia > 0.00693 (upper half)  AND  celltype_astrocyte <= 0.00686 (lower half) | 2 | 47 | 0.5893 | 0.5676 | 1.71 | 0.4822 |
| celltype_endothelial > 0.125 (high third)  AND  PlaqueTotal > 12.1 (upper half) | 2 | 41 | 0.5878 | 0.5648 | 1.63 | 0.4822 |
| PlaqueTotal > 13.5 (high third)  AND  diagnosis = Alzheimer Disease | 2 | 46 | 0.5882 | 0.5676 | 1.62 | 0.4822 |
| sex = female  AND  TangleTotal <= 8 (lower half) | 2 | 45 | 0.5877 | 0.5648 | 1.62 | 0.4822 |
| celltype_endothelial > 0.125 (high third)  AND  CERAD = 3 | 2 | 41 | 0.5876 | 0.5648 | 1.61 | 0.4822 |
| celltype_microglia > 0.00693 (upper half) | 1 | 95 | 0.5882 | 0.5750 | 1.59 | 0.4822 |
| sex = female  AND  diagnosis = control | 2 | 44 | 0.5866 | 0.5648 | 1.54 | 0.4822 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| PlaqueTotal > 13.5 (high third)  AND  CERAD = 3 | 50 | 50 | 0.5980 | 0.5694 | 2.39 | 0.5237 |
| celltype_microglia > 0.00693 (upper half)  AND  apoeGenotype = e3-3 | 44 | 44 | 0.5973 | 0.5663 | 2.26 | 0.5237 |
| apoeGenotype = e3-3  AND  celltype_neuron <= -1.75e-05 (lower half) | 41 | 41 | 0.5966 | 0.5663 | 2.20 | 0.5237 |
| celltype_microglia > 0.00693 (upper half)  AND  diagnosis = control | 45 | 45 | 0.5912 | 0.5663 | 1.81 | 0.5237 |
| celltype_microglia > 0.00693 (upper half)  AND  ageDeath > 86 (upper half) | 42 | 42 | 0.5899 | 0.5663 | 1.72 | 0.5237 |
| celltype_microglia > 0.00693 (upper half)  AND  TangleTotal <= 8 (lower half) | 40 | 40 | 0.5899 | 0.5663 | 1.71 | 0.5237 |
| apoeGenotype = e3-3  AND  CERAD = 3 | 44 | 44 | 0.5894 | 0.5663 | 1.68 | 0.5237 |
| celltype_microglia > 0.00693 (upper half)  AND  celltype_astrocyte <= 0.00686 (lower half) | 47 | 47 | 0.5893 | 0.5694 | 1.67 | 0.5237 |
| celltype_endothelial > 0.125 (high third) | 63 | 63 | 0.5876 | 0.5679 | 1.66 | 0.5237 |
| PlaqueTotal > 13.5 (high third)  AND  diagnosis = Alzheimer Disease | 46 | 46 | 0.5882 | 0.5694 | 1.57 | 0.5237 |
| celltype_endothelial > 0.125 (high third)  AND  PlaqueTotal > 12.1 (upper half) | 41 | 41 | 0.5878 | 0.5663 | 1.57 | 0.5237 |
| sex = female  AND  TangleTotal <= 8 (lower half) | 45 | 45 | 0.5877 | 0.5663 | 1.55 | 0.5237 |
| celltype_endothelial > 0.125 (high third)  AND  CERAD = 3 | 41 | 41 | 0.5876 | 0.5663 | 1.55 | 0.5237 |
| PlaqueTotal > 13.5 (high third)  AND  TangleTotal > 8 (upper half) | 46 | 46 | 0.5871 | 0.5694 | 1.48 | 0.5237 |
| sex = female  AND  diagnosis = control | 44 | 44 | 0.5866 | 0.5663 | 1.48 | 0.5237 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_endothelial -0.261-0.125 (middle third) | 1 | 63 | 0.5470 | 0.5687 | -2.09 | 0.9851 |
| PlaqueTotal <= 8 (low third) | 1 | 66 | 0.5576 | 0.5710 | -1.34 | 0.9127 |
| celltype_microglia <= 0.00693 (lower half) | 1 | 95 | 0.5642 | 0.5750 | -1.30 | 0.9127 |
| lastMMSE > 25 (upper half) | 1 | 91 | 0.5648 | 0.5750 | -1.24 | 0.9072 |
| TangleTotal > 11.8 (high third) | 1 | 63 | 0.5578 | 0.5687 | -1.05 | 0.8889 |
| celltype_oligodendrocyte <= -0.462 (low third) | 1 | 64 | 0.5585 | 0.5687 | -0.98 | 0.8755 |
| lastMMSE 18-28 (middle third) | 1 | 75 | 0.5643 | 0.5731 | -0.93 | 0.8481 |
| celltype_astrocyte > 0.00686 (upper half) | 1 | 95 | 0.5677 | 0.5750 | -0.88 | 0.8672 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| Braak = 4 | 57 | 0.5797 | 0.5687 | 1.05 | 0.5699 |
| Braak = 5 | 48 | 0.5589 | 0.5676 | -0.69 | 0.8019 |
| CERAD = 3 | 106 | 0.5853 | 0.5773 | 1.00 | 0.5786 |
| PlaqueTotal 8-13.5 (middle third) | 69 | 0.5737 | 0.5710 | 0.27 | 0.6582 |
| PlaqueTotal <= 12.1 (lower half) | 95 | 0.5716 | 0.5750 | -0.40 | 0.7604 |
| PlaqueTotal <= 8 (low third) | 66 | 0.5576 | 0.5710 | -1.34 | 0.9127 |
| PlaqueTotal > 12.1 (upper half) | 95 | 0.5801 | 0.5750 | 0.62 | 0.6510 |
| PlaqueTotal > 13.5 (high third) | 55 | 0.5814 | 0.5687 | 1.21 | 0.4822 |
| TangleTotal 6-11.8 (middle third) | 63 | 0.5765 | 0.5687 | 0.74 | 0.6510 |
| TangleTotal <= 6 (low third) | 64 | 0.5650 | 0.5687 | -0.36 | 0.7414 |
| TangleTotal <= 8 (lower half) | 98 | 0.5756 | 0.5750 | 0.07 | 0.7174 |
| TangleTotal > 11.8 (high third) | 63 | 0.5578 | 0.5687 | -1.05 | 0.8889 |
| TangleTotal > 8 (upper half) | 92 | 0.5742 | 0.5750 | -0.09 | 0.7291 |
| ageDeath 82-88 (middle third) | 54 | 0.5596 | 0.5676 | -0.63 | 0.8019 |
| ageDeath <= 82 (low third) | 73 | 0.5752 | 0.5710 | 0.42 | 0.6582 |
| ageDeath <= 86 (lower half) | 104 | 0.5741 | 0.5750 | -0.11 | 0.7414 |
| ageDeath > 86 (upper half) | 86 | 0.5767 | 0.5750 | 0.21 | 0.6582 |
| ageDeath > 88 (high third) | 63 | 0.5715 | 0.5687 | 0.26 | 0.6582 |
| apoeGenotype = e3-3 | 92 | 0.5874 | 0.5750 | 1.52 | 0.4822 |
| apoeGenotype = e3-4 | 62 | 0.5616 | 0.5687 | -0.68 | 0.8231 |
| diagnosis = Alzheimer Disease | 90 | 0.5734 | 0.5750 | -0.19 | 0.7414 |
| diagnosis = control | 100 | 0.5731 | 0.5750 | -0.23 | 0.7414 |
| lastMMSE 18-28 (middle third) | 75 | 0.5643 | 0.5731 | -0.93 | 0.8481 |
| lastMMSE <= 18 (low third) | 67 | 0.5726 | 0.5710 | 0.16 | 0.6711 |
| lastMMSE <= 25 (lower half) | 99 | 0.5832 | 0.5750 | 1.00 | 0.5731 |
| lastMMSE > 25 (upper half) | 91 | 0.5648 | 0.5750 | -1.24 | 0.9072 |
| lastMMSE > 28 (high third) | 48 | 0.5732 | 0.5676 | 0.44 | 0.6582 |
| pmi 2.5-3.25 (middle third) | 66 | 0.5733 | 0.5710 | 0.24 | 0.6582 |
| pmi <= 2.5 (low third) | 65 | 0.5729 | 0.5687 | 0.40 | 0.6582 |
| pmi <= 3 (lower half) | 117 | 0.5804 | 0.5772 | 0.49 | 0.6582 |
| pmi > 3 (upper half) | 73 | 0.5624 | 0.5710 | -0.87 | 0.8412 |
| pmi > 3.25 (high third) | 59 | 0.5619 | 0.5687 | -0.66 | 0.8231 |
| sex = female | 81 | 0.5850 | 0.5731 | 1.25 | 0.4925 |
| sex = male | 109 | 0.5713 | 0.5773 | -0.74 | 0.8116 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.218-0.213 (middle third) | 63 | 0.5611 | 0.5687 | -0.73 | 0.8322 |
| celltype_astrocyte <= -0.218 (low third) | 64 | 0.5784 | 0.5687 | 0.92 | 0.5786 |
| celltype_astrocyte <= 0.00686 (lower half) | 95 | 0.5809 | 0.5750 | 0.71 | 0.6510 |
| celltype_astrocyte > 0.00686 (upper half) | 95 | 0.5677 | 0.5750 | -0.88 | 0.8672 |
| celltype_astrocyte > 0.213 (high third) | 63 | 0.5711 | 0.5687 | 0.22 | 0.6582 |
| celltype_endothelial -0.261-0.125 (middle third) | 63 | 0.5470 | 0.5687 | -2.09 | 0.9851 |
| celltype_endothelial <= -0.0477 (lower half) | 95 | 0.5719 | 0.5750 | -0.38 | 0.7570 |
| celltype_endothelial <= -0.261 (low third) | 64 | 0.5661 | 0.5687 | -0.25 | 0.7414 |
| celltype_endothelial > -0.0477 (upper half) | 95 | 0.5738 | 0.5750 | -0.14 | 0.7414 |
| celltype_endothelial > 0.125 (high third) | 63 | 0.5876 | 0.5687 | 1.81 | 0.4822 |
| celltype_microglia -0.272-0.252 (middle third) | 63 | 0.5659 | 0.5687 | -0.27 | 0.7414 |
| celltype_microglia <= -0.272 (low third) | 64 | 0.5733 | 0.5687 | 0.44 | 0.6582 |
| celltype_microglia <= 0.00693 (lower half) | 95 | 0.5642 | 0.5750 | -1.30 | 0.9127 |
| celltype_microglia > 0.00693 (upper half) | 95 | 0.5882 | 0.5750 | 1.59 | 0.4822 |
| celltype_microglia > 0.252 (high third) | 63 | 0.5721 | 0.5687 | 0.32 | 0.6582 |
| celltype_neuron -0.147-0.187 (middle third) | 63 | 0.5743 | 0.5687 | 0.54 | 0.6582 |
| celltype_neuron <= -0.147 (low third) | 64 | 0.5663 | 0.5687 | -0.23 | 0.7414 |
| celltype_neuron <= -1.75e-05 (lower half) | 95 | 0.5738 | 0.5750 | -0.15 | 0.7414 |
| celltype_neuron > -1.75e-05 (upper half) | 95 | 0.5777 | 0.5750 | 0.33 | 0.6582 |
| celltype_neuron > 0.187 (high third) | 63 | 0.5759 | 0.5687 | 0.69 | 0.6510 |
| celltype_oligodendrocyte -0.462-0.276 (middle third) | 63 | 0.5683 | 0.5687 | -0.04 | 0.7174 |
| celltype_oligodendrocyte <= -0.129 (lower half) | 95 | 0.5769 | 0.5750 | 0.23 | 0.6582 |
| celltype_oligodendrocyte <= -0.462 (low third) | 64 | 0.5585 | 0.5687 | -0.98 | 0.8755 |
| celltype_oligodendrocyte > -0.129 (upper half) | 95 | 0.5756 | 0.5750 | 0.07 | 0.7174 |
| celltype_oligodendrocyte > 0.276 (high third) | 63 | 0.5804 | 0.5687 | 1.11 | 0.5642 |

### 5.5 BannerLFQ

Whole-cohort adjusted AUROC: 0.5369

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_neuron <= -0.182 (low third)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 2 | 43 | 0.5778 | 0.5358 | 3.79 | 0.1729 |
| celltype_neuron <= 0.0182 (lower half)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 2 | 60 | 0.5671 | 0.5379 | 3.34 | 0.1729 |
| celltype_neuron <= 0.0182 (lower half)  AND  ageDeath <= 86 (lower half) | 2 | 55 | 0.5645 | 0.5379 | 3.05 | 0.1729 |
| celltype_neuron <= 0.0182 (lower half)  AND  pmi <= 3 (lower half) | 2 | 53 | 0.5657 | 0.5363 | 2.99 | 0.1729 |
| celltype_neuron <= 0.0182 (lower half) | 1 | 95 | 0.5573 | 0.5404 | 2.90 | 0.2305 |
| ageDeath <= 86 (lower half) | 1 | 104 | 0.5556 | 0.5404 | 2.62 | 0.2305 |
| ageDeath <= 86 (lower half)  AND  sex = male | 2 | 68 | 0.5561 | 0.5377 | 2.41 | 0.2305 |
| CERAD = 3 | 1 | 106 | 0.5530 | 0.5409 | 2.40 | 0.2470 |
| celltype_neuron <= 0.0182 (lower half)  AND  celltype_oligodendrocyte > 0.391 (high third) | 2 | 45 | 0.5615 | 0.5358 | 2.32 | 0.2305 |
| celltype_neuron <= 0.0182 (lower half)  AND  celltype_astrocyte <= -0.0155 (lower half) | 2 | 45 | 0.5614 | 0.5358 | 2.31 | 0.2305 |
| CERAD = 3  AND  celltype_neuron <= -0.182 (low third) | 2 | 46 | 0.5588 | 0.5363 | 2.29 | 0.2470 |
| ageDeath <= 86 (lower half)  AND  PlaqueTotal > 12.1 (upper half) | 2 | 57 | 0.5575 | 0.5379 | 2.26 | 0.2470 |
| celltype_neuron <= -0.182 (low third) | 1 | 64 | 0.5575 | 0.5379 | 2.25 | 0.2470 |
| ageDeath <= 86 (lower half)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 2 | 57 | 0.5561 | 0.5379 | 2.09 | 0.2470 |
| ageDeath <= 86 (lower half)  AND  TangleTotal > 8 (upper half) | 2 | 62 | 0.5546 | 0.5379 | 1.92 | 0.2634 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| celltype_neuron <= -0.182 (low third)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 43 | 43 | 0.5778 | 0.5353 | 3.93 | 0.1482 |
| celltype_neuron <= 0.0182 (lower half)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 60 | 60 | 0.5671 | 0.5377 | 3.63 | 0.1482 |
| celltype_neuron <= 0.0182 (lower half)  AND  ageDeath <= 86 (lower half) | 55 | 55 | 0.5645 | 0.5377 | 3.31 | 0.1482 |
| celltype_neuron <= 0.0182 (lower half)  AND  pmi <= 3 (lower half) | 53 | 53 | 0.5657 | 0.5354 | 3.22 | 0.1482 |
| celltype_neuron <= 0.0182 (lower half) | 95 | 95 | 0.5573 | 0.5401 | 2.97 | 0.1482 |
| ageDeath <= 86 (lower half) | 104 | 104 | 0.5556 | 0.5401 | 2.68 | 0.1482 |
| CERAD = 3  AND  celltype_neuron <= -0.182 (low third) | 46 | 46 | 0.5588 | 0.5354 | 2.49 | 0.1482 |
| ageDeath <= 86 (lower half)  AND  PlaqueTotal > 12.1 (upper half) | 57 | 57 | 0.5575 | 0.5377 | 2.45 | 0.1482 |
| celltype_neuron <= -0.182 (low third) | 64 | 64 | 0.5575 | 0.5377 | 2.45 | 0.1482 |
| celltype_neuron <= 0.0182 (lower half)  AND  celltype_oligodendrocyte > 0.391 (high third) | 45 | 45 | 0.5615 | 0.5353 | 2.42 | 0.1537 |
| celltype_neuron <= 0.0182 (lower half)  AND  celltype_astrocyte <= -0.0155 (lower half) | 45 | 45 | 0.5614 | 0.5353 | 2.42 | 0.1537 |
| ageDeath <= 86 (lower half)  AND  sex = male | 68 | 68 | 0.5561 | 0.5397 | 2.34 | 0.1482 |
| ageDeath <= 86 (lower half)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 57 | 57 | 0.5561 | 0.5377 | 2.27 | 0.1482 |
| CERAD = 3 | 106 | 106 | 0.5530 | 0.5408 | 2.27 | 0.1482 |
| ageDeath <= 86 (lower half)  AND  TangleTotal > 8 (upper half) | 62 | 62 | 0.5546 | 0.5377 | 2.09 | 0.1482 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| lastMMSE > 25 (upper half) | 1 | 91 | 0.5209 | 0.5410 | -3.06 | 1.0000 |
| celltype_neuron > 0.0182 (upper half) | 1 | 95 | 0.5235 | 0.5404 | -2.91 | 1.0000 |
| lastMMSE > 28 (high third) | 1 | 48 | 0.5094 | 0.5363 | -2.73 | 1.0000 |
| ageDeath > 86 (upper half) | 1 | 86 | 0.5233 | 0.5410 | -2.70 | 1.0000 |
| PlaqueTotal <= 12.1 (lower half) | 1 | 95 | 0.5255 | 0.5404 | -2.56 | 1.0000 |
| PlaqueTotal <= 8 (low third) | 1 | 66 | 0.5206 | 0.5377 | -2.26 | 1.0000 |
| TangleTotal <= 8 (lower half) | 1 | 98 | 0.5282 | 0.5404 | -2.10 | 1.0000 |
| celltype_endothelial > -0.0629 (upper half) | 1 | 50 | 0.5162 | 0.5363 | -2.04 | 1.0000 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| Braak = 4 | 57 | 0.5371 | 0.5379 | -0.09 | 0.7343 |
| Braak = 5 | 48 | 0.5248 | 0.5363 | -1.17 | 0.9760 |
| CERAD = 3 | 106 | 0.5530 | 0.5409 | 2.40 | 0.2470 |
| PlaqueTotal 8-13.5 (middle third) | 69 | 0.5474 | 0.5377 | 1.28 | 0.3804 |
| PlaqueTotal <= 12.1 (lower half) | 95 | 0.5255 | 0.5404 | -2.56 | 1.0000 |
| PlaqueTotal <= 8 (low third) | 66 | 0.5206 | 0.5377 | -2.26 | 1.0000 |
| PlaqueTotal > 12.1 (upper half) | 95 | 0.5500 | 0.5404 | 1.65 | 0.3557 |
| PlaqueTotal > 13.5 (high third) | 55 | 0.5381 | 0.5379 | 0.03 | 0.7343 |
| TangleTotal 6-11.8 (middle third) | 63 | 0.5429 | 0.5379 | 0.58 | 0.5369 |
| TangleTotal <= 6 (low third) | 64 | 0.5217 | 0.5379 | -1.86 | 1.0000 |
| TangleTotal <= 8 (lower half) | 98 | 0.5282 | 0.5404 | -2.10 | 1.0000 |
| TangleTotal > 11.8 (high third) | 63 | 0.5397 | 0.5379 | 0.21 | 0.7327 |
| TangleTotal > 8 (upper half) | 92 | 0.5468 | 0.5410 | 0.89 | 0.4648 |
| ageDeath 82-88 (middle third) | 54 | 0.5326 | 0.5363 | -0.37 | 0.8139 |
| ageDeath <= 82 (low third) | 73 | 0.5493 | 0.5377 | 1.52 | 0.3557 |
| ageDeath <= 86 (lower half) | 104 | 0.5556 | 0.5404 | 2.62 | 0.2305 |
| ageDeath > 86 (upper half) | 86 | 0.5233 | 0.5410 | -2.70 | 1.0000 |
| ageDeath > 88 (high third) | 63 | 0.5233 | 0.5379 | -1.67 | 1.0000 |
| apoeGenotype = e3-3 | 92 | 0.5320 | 0.5410 | -1.36 | 1.0000 |
| apoeGenotype = e3-4 | 62 | 0.5373 | 0.5379 | -0.07 | 0.7343 |
| diagnosis = Alzheimer Disease | 90 | 0.5455 | 0.5410 | 0.69 | 0.5233 |
| diagnosis = control | 100 | 0.5315 | 0.5404 | -1.53 | 1.0000 |
| lastMMSE 18-28 (middle third) | 75 | 0.5458 | 0.5401 | 0.92 | 0.4648 |
| lastMMSE <= 18 (low third) | 67 | 0.5379 | 0.5377 | 0.02 | 0.7343 |
| lastMMSE <= 25 (lower half) | 99 | 0.5507 | 0.5404 | 1.77 | 0.2634 |
| lastMMSE > 25 (upper half) | 91 | 0.5209 | 0.5410 | -3.06 | 1.0000 |
| lastMMSE > 28 (high third) | 48 | 0.5094 | 0.5363 | -2.73 | 1.0000 |
| pmi 2.5-3.25 (middle third) | 66 | 0.5354 | 0.5377 | -0.31 | 0.7740 |
| pmi <= 2.5 (low third) | 65 | 0.5428 | 0.5379 | 0.56 | 0.5389 |
| pmi <= 3 (lower half) | 117 | 0.5396 | 0.5405 | -0.19 | 0.7568 |
| pmi > 3 (upper half) | 73 | 0.5367 | 0.5377 | -0.13 | 0.7343 |
| pmi > 3.25 (high third) | 59 | 0.5345 | 0.5379 | -0.39 | 0.8139 |
| sex = female | 81 | 0.5377 | 0.5401 | -0.37 | 0.8139 |
| sex = male | 109 | 0.5410 | 0.5409 | 0.01 | 0.7343 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.237-0.251 (middle third) | 63 | 0.5475 | 0.5379 | 1.11 | 0.4034 |
| celltype_astrocyte <= -0.0155 (lower half) | 95 | 0.5379 | 0.5404 | -0.43 | 0.8250 |
| celltype_astrocyte <= -0.237 (low third) | 64 | 0.5318 | 0.5379 | -0.70 | 0.8760 |
| celltype_astrocyte > -0.0155 (upper half) | 95 | 0.5394 | 0.5404 | -0.17 | 0.7568 |
| celltype_astrocyte > 0.251 (high third) | 63 | 0.5259 | 0.5379 | -1.38 | 1.0000 |
| celltype_endothelial <= -0.0629 (lower half) | 51 | 0.5369 | 0.5363 | 0.06 | 0.7327 |
| celltype_endothelial > -0.0629 (upper half) | 50 | 0.5162 | 0.5363 | -2.04 | 1.0000 |
| celltype_neuron -0.182-0.263 (middle third) | 63 | 0.5336 | 0.5379 | -0.49 | 0.8406 |
| celltype_neuron <= -0.182 (low third) | 64 | 0.5575 | 0.5379 | 2.25 | 0.2470 |
| celltype_neuron <= 0.0182 (lower half) | 95 | 0.5573 | 0.5404 | 2.90 | 0.2305 |
| celltype_neuron > 0.0182 (upper half) | 95 | 0.5235 | 0.5404 | -2.91 | 1.0000 |
| celltype_neuron > 0.263 (high third) | 63 | 0.5250 | 0.5379 | -1.48 | 1.0000 |
| celltype_oligodendrocyte -0.464-0.391 (middle third) | 63 | 0.5326 | 0.5379 | -0.60 | 0.8674 |
| celltype_oligodendrocyte <= -0.073 (lower half) | 95 | 0.5318 | 0.5404 | -1.48 | 1.0000 |
| celltype_oligodendrocyte <= -0.464 (low third) | 64 | 0.5288 | 0.5379 | -1.04 | 0.9760 |
| celltype_oligodendrocyte > -0.073 (upper half) | 95 | 0.5447 | 0.5404 | 0.74 | 0.5000 |
| celltype_oligodendrocyte > 0.391 (high third) | 63 | 0.5434 | 0.5379 | 0.64 | 0.5233 |

---

## 6. The graphs (HIW 4.8–4.10)

### 6.1 Edge counts and sizes

| dataset | prior_edges | testable_n20 | W_positive | q_below_05_any_direction | data_supported | pct_of_prior_edges | median_n_W_positive | min_n_W_positive |
|---|---|---|---|---|---|---|---|---|
| ROSMAP | 1308 | 1308 | 654 | 502 | 251 | 19.2000 | 400 | 120 |
| Diverse | 1308 | 1308 | 654 | 928 | 464 | 35.5000 | 676 | 455 |
| Banner | 1308 | 1308 | 654 | 286 | 143 | 10.9000 | 189 | 33 |
| BannerLFQ | 1308 | 772 | 409 | 488 | 244 | 18.7000 | 190 | 0 |

| dataset | W_median_all | W_max | W_median_supported | W_min_supported | beta_median_supported | beta_min | beta_max | frac_beta_positive | n_abs_beta_above_2 |
|---|---|---|---|---|---|---|---|---|---|
| ROSMAP | 0.0373 | 0.7317 | 0.0880 | 0.0235 | 0.0795 | -0.7993 | 1.0704 | 0.7170 | 0 |
| Diverse | 0.0592 | 0.8012 | 0.0880 | 0.0174 | 0.0685 | -0.6355 | 2.1548 | 0.7090 | 1 |
| Banner | 0.0384 | 0.8434 | 0.1481 | 0.0379 | 0.1612 | -1.6830 | 2.0510 | 0.7620 | 1 |
| BannerLFQ | 0.0966 | 0.8006 | 0.1500 | 0.0341 | 0.1089 | -1.2499 | 3.0142 | 0.7300 | 1 |

### 6.2 How each prior pair's arrow and weight were set (HIW 4.8)

| | one-way pair, direction test agrees with C (full weight) | one-way pair, direction test disagrees (kept in C's direction at 0.3 weight) | two-way pair (direction test picks the arrow) | n < 20 (placeholder 0.1·C) |
|---|---|---|---|---|
| ROSMAP | 0 | 0 | 654 | 0 |
| Diverse | 0 | 0 | 654 | 0 |
| Banner | 0 | 0 | 654 | 0 |
| BannerLFQ | 0 | 0 | 386 | 23 |

Share of one-way pairs where the direction test agrees with C: ROSMAP — (no one-way pairs), Diverse — (no one-way pairs), Banner — (no one-way pairs), BannerLFQ — (no one-way pairs).

### 6.3 Supported edges by prior confidence

| prior C | prior edges | ROSMAP | Diverse | Banner | BannerLFQ |
|---|---|---|---|---|---|
| 0.2 | 403 | 83 (21%) | 143 (35%) | 32 (8%) | 78 (19%) |
| 0.4 | 230 | 48 (21%) | 76 (33%) | 17 (7%) | 44 (19%) |
| 0.6 | 399 | 55 (14%) | 135 (34%) | 45 (11%) | 70 (18%) |
| 0.8 | 123 | 18 (15%) | 40 (33%) | 15 (12%) | 15 (12%) |
| 1.0 | 153 | 47 (31%) | 70 (46%) | 34 (22%) | 37 (24%) |

### 6.4 Hubs of the supported graph (number of supported edges touching each protein)

- **ROSMAP** top 15: HSPA5 21, COPB1 20, TMED10 20, SEC13 19, COPB2 18, SEC31A 16, COPG1 16, COPE 14, COPA 14, SEC24C 13, TMED2 13, BACE1 12, AP2B1 11, RAB1B 11, RAB1A 10
  - proteins with no supported edge (4): ADAM17, AP4B1, SIGMAR1, VLDLR
- **Diverse** top 15: COPB1 27, SEC13 26, TMED10 24, COPA 23, SEC23A 23, TMED2 23, COPB2 22, HSPA5 22, SEC24C 21, SEC31A 20, SAR1A 20, SEC16A 19, PREB 19, SEC24B 19, RAB1A 18
  - proteins with no supported edge (1): ECE1
- **Banner** top 15: APP 14, COPB1 12, AP2B1 11, COPB2 9, TMED2 9, COPA 9, TMED9 9, HSPA5 9, RAB1A 9, SEC13 9, SEC24C 9, COPG1 8, COPE 8, SEC23A 8, TMED10 8
  - proteins with no supported edge (17): ABCA7, ADAM17, AP4S1, GGA1, IDE, MAP2K3, MAPK14, MAPKAPK2, PLCG2, PLG, SEC16A, SEC23B, SEC24A, SIGMAR1, SREBF2, SYK, VLDLR
- **BannerLFQ** top 15: HSPA5 20, APP 20, COPA 20, TMED10 19, COPB1 18, SEC13 16, RAB1A 16, COPB2 15, SEC24B 14, SORT1 14, COPZ1 13, COPG1 13, RAB1B 13, SAR1B 13, VDAC1 13
  - proteins with no supported edge (26): ABCA7, ADAM17, AP4B1, AP4E1, AP4M1, AP4S1, BACE1, ECE1, EIF2S1, GGA2, GORASP1, MAP2K3, MAPK14, MAPKAPK2, NCSTN, PLCG2, PLG, PSEN1, SEC16A, SEC23B, SEC24A, SEC24D, SIGMAR1, SREBF2, SYK, VLDLR

### 6.5 Supported edges, strongest first (top 40 per dataset; full lists in `{ds}_causal_edges.csv`)

**ROSMAP** (251 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| AP2A1 | AP2B1 | 1.0000 | 0.7317 | 0.5775 | 400 | 1.96e-65 |
| COPB1 | COPG1 | 1.0000 | 0.6281 | 0.4257 | 400 | 6.07e-43 |
| COPB1 | COPE | 1.0000 | 0.5092 | 0.5909 | 400 | 9.76e-26 |
| SEC24C | SEC23A | 1.0000 | 0.5014 | 0.5325 | 400 | 6.19e-25 |
| COPB2 | COPB1 | 1.0000 | 0.4935 | 0.4956 | 400 | 4.47e-24 |
| COPZ1 | COPB1 | 1.0000 | 0.4731 | 0.2813 | 400 | 5.79e-22 |
| COPB2 | COPE | 1.0000 | 0.4676 | 0.5579 | 400 | 2.03e-21 |
| COPE | COPG1 | 1.0000 | 0.4545 | 0.0755 | 400 | 4.12e-20 |
| SEC31A | SEC13 | 1.0000 | 0.4368 | 0.3277 | 400 | 1.88e-18 |
| COPB2 | COPG1 | 1.0000 | 0.4134 | 0.3334 | 400 | 2.00e-16 |
| TMED10 | TMED9 | 0.8000 | 0.4037 | 0.5392 | 400 | 2.97e-25 |
| COPB1 | TMED10 | 1.0000 | 0.3983 | 0.6316 | 400 | 3.46e-15 |
| TMED10 | TMED2 | 1.0000 | 0.3896 | 0.5420 | 400 | 1.65e-14 |
| COPB2 | TMED10 | 1.0000 | 0.3834 | 0.2440 | 400 | 4.91e-14 |
| COPG1 | COPA | 1.0000 | 0.3715 | 0.2016 | 400 | 3.72e-13 |
| COPZ1 | TMED10 | 1.0000 | 0.3565 | 0.0699 | 400 | 4.44e-12 |
| COPZ1 | COPE | 1.0000 | 0.3549 | 0.1813 | 400 | 5.58e-12 |
| COPZ1 | COPB2 | 1.0000 | 0.3410 | 0.1863 | 400 | 4.34e-11 |
| COPB1 | COPA | 1.0000 | 0.3401 | 0.1062 | 400 | 4.85e-11 |
| AP2A1 | AP1G1 | 0.6000 | 0.3275 | 0.5574 | 400 | 3.23e-30 |
| HSPA5 | HSPA9 | 0.8000 | 0.3263 | 0.4301 | 400 | 5.69e-16 |
| COPZ1 | COPG1 | 1.0000 | 0.3192 | 0.0112 | 400 | 1.04e-09 |
| VPS35 | VPS29 | 1.0000 | 0.3186 | 0.1995 | 400 | 1.08e-09 |
| COPB2 | TMED2 | 1.0000 | 0.3098 | 0.4738 | 400 | 3.57e-09 |
| HSPA5 | TMED10 | 0.8000 | 0.3047 | 0.1749 | 400 | 7.35e-14 |
| PICALM | AP2B1 | 1.0000 | 0.2976 | 0.0070 | 400 | 1.77e-08 |
| AP2A1 | PICALM | 1.0000 | 0.2951 | 0.0595 | 400 | 2.38e-08 |
| SEC24C | SEC31A | 0.6000 | 0.2880 | 0.3402 | 400 | 1.26e-22 |
| COPE | TMED10 | 1.0000 | 0.2823 | 0.0350 | 400 | 1.11e-07 |
| SEC31A | SEC23A | 0.8000 | 0.2822 | 0.0759 | 400 | 7.42e-12 |
| COPB2 | TMED9 | 1.0000 | 0.2730 | 0.1459 | 400 | 3.11e-07 |
| COPB1 | COPG2 | 1.0000 | 0.2728 | -0.0673 | 400 | 3.11e-07 |
| COPB2 | COPG2 | 1.0000 | 0.2636 | 0.1295 | 400 | 8.60e-07 |
| EIF2S1 | HSPA5 | 0.6000 | 0.2598 | 0.6119 | 400 | 3.97e-18 |
| HSPA9 | MFN2 | 0.6000 | 0.2559 | 0.3123 | 400 | 1.44e-17 |
| COPE | COPG2 | 1.0000 | 0.2510 | 0.2421 | 400 | 3.05e-06 |
| PLCG2 | SYK | 0.6000 | 0.2461 | 0.3042 | 288 | 8.13e-12 |
| COPA | COPG2 | 1.0000 | 0.2446 | 0.5225 | 400 | 5.76e-06 |
| EIF2S1 | BACE1 | 0.8000 | 0.2360 | -0.1563 | 400 | 2.38e-08 |
| COPB2 | COPA | 1.0000 | 0.2206 | 0.0585 | 400 | 5.55e-05 |

**Diverse** (464 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| AP2B1 | AP2A1 | 1.0000 | 0.8012 | 0.7412 | 676 | 1.16e-149 |
| COPA | COPB1 | 1.0000 | 0.7119 | 0.5629 | 676 | 2.40e-103 |
| SEC24C | SEC23A | 1.0000 | 0.6678 | 0.5419 | 676 | 1.58e-86 |
| COPA | COPG1 | 1.0000 | 0.6643 | 0.4303 | 676 | 2.35e-85 |
| COPB1 | COPG1 | 1.0000 | 0.6619 | 0.5039 | 676 | 1.48e-84 |
| TMED10 | TMED2 | 1.0000 | 0.6516 | 0.5780 | 676 | 3.95e-81 |
| SEC31A | SEC13 | 1.0000 | 0.6451 | 0.3643 | 676 | 4.84e-79 |
| COPB1 | COPB2 | 1.0000 | 0.6386 | 0.4275 | 676 | 5.56e-77 |
| COPB2 | COPA | 1.0000 | 0.6358 | 0.4201 | 676 | 3.76e-76 |
| VPS26A | VPS29 | 1.0000 | 0.5949 | 0.3440 | 676 | 2.40e-64 |
| VPS35 | VPS29 | 1.0000 | 0.5842 | 0.3721 | 676 | 1.40e-61 |
| TMED9 | TMED10 | 0.8000 | 0.5431 | 0.6757 | 676 | 1.82e-90 |
| COPB2 | COPZ1 | 1.0000 | 0.5314 | 0.7592 | 676 | 4.36e-49 |
| COPE | COPA | 1.0000 | 0.5180 | 0.1526 | 676 | 2.97e-46 |
| COPB2 | COPE | 1.0000 | 0.5073 | 0.3412 | 676 | 4.22e-44 |
| COPB1 | COPE | 1.0000 | 0.4989 | 0.3037 | 676 | 1.72e-42 |
| COPG1 | COPB2 | 1.0000 | 0.4965 | 0.1245 | 676 | 4.95e-42 |
| COPZ1 | COPB1 | 1.0000 | 0.4921 | 0.1182 | 676 | 3.06e-41 |
| SEC23A | SEC24B | 1.0000 | 0.4874 | 0.1906 | 676 | 2.40e-40 |
| VPS26A | VPS35 | 1.0000 | 0.4711 | 0.5634 | 676 | 2.06e-37 |
| AP2A1 | AP3B2 | 0.6000 | 0.4611 | 0.8404 | 676 | 2.86e-130 |
| TMED9 | TMED2 | 0.8000 | 0.4544 | 0.2119 | 676 | 1.75e-57 |
| COPB2 | TMED9 | 1.0000 | 0.4199 | 0.3155 | 676 | 3.91e-29 |
| COPZ1 | COPE | 1.0000 | 0.4156 | 0.0941 | 676 | 1.69e-28 |
| AP2A1 | AP1G1 | 0.6000 | 0.4136 | 0.3116 | 676 | 2.64e-94 |
| SEC31A | SEC16A | 0.8000 | 0.4121 | 0.2933 | 676 | 1.12e-45 |
| COPB2 | TMED10 | 1.0000 | 0.4059 | -0.0054 | 676 | 3.89e-27 |
| SEC31A | SEC23A | 0.8000 | 0.4001 | 0.2482 | 676 | 1.05e-42 |
| COPA | COPZ1 | 1.0000 | 0.3999 | 0.1166 | 676 | 2.63e-26 |
| COPA | COPG2 | 1.0000 | 0.3963 | 0.6100 | 676 | 7.92e-26 |
| COPB1 | COPG2 | 1.0000 | 0.3935 | 0.4399 | 676 | 1.88e-25 |
| SEC23A | SEC24A | 1.0000 | 0.3922 | 0.3301 | 676 | 2.79e-25 |
| COPE | COPG1 | 1.0000 | 0.3879 | 0.0455 | 676 | 9.97e-25 |
| SYK | PLCG2 | 1.0000 | 0.3795 | 0.4479 | 664 | 2.84e-23 |
| AP2B1 | AP3B2 | 0.6000 | 0.3736 | -0.0090 | 676 | 3.59e-72 |
| COPE | TMED9 | 1.0000 | 0.3705 | 0.2015 | 676 | 1.53e-22 |
| COPE | TMED10 | 1.0000 | 0.3633 | 0.0259 | 676 | 1.18e-21 |
| COPZ1 | TMED10 | 1.0000 | 0.3601 | 0.1097 | 676 | 2.88e-21 |
| AP3B2 | AP1G1 | 0.6000 | 0.3550 | 0.1455 | 676 | 1.78e-63 |
| MFN2 | VDAC1 | 0.6000 | 0.3541 | 0.9207 | 676 | 4.04e-63 |

**Banner** (143 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| AP2B1 | AP2A1 | 1.0000 | 0.8434 | 0.7454 | 189 | 1.60e-49 |
| COPB2 | COPB1 | 1.0000 | 0.5310 | 0.1839 | 189 | 5.00e-13 |
| TMED2 | TMED10 | 1.0000 | 0.5143 | 0.6470 | 189 | 4.06e-12 |
| COPG1 | COPB2 | 1.0000 | 0.4558 | 0.2496 | 189 | 2.38e-09 |
| COPE | COPB1 | 1.0000 | 0.4479 | 0.1612 | 189 | 4.80e-09 |
| VPS35 | VPS29 | 1.0000 | 0.4418 | 0.3956 | 189 | 7.14e-09 |
| COPG1 | COPB1 | 1.0000 | 0.4058 | 0.0684 | 189 | 2.26e-07 |
| COPE | COPG1 | 1.0000 | 0.4025 | 0.2262 | 189 | 2.91e-07 |
| COPE | COPB2 | 1.0000 | 0.4012 | 0.3086 | 189 | 3.15e-07 |
| COPA | COPB2 | 1.0000 | 0.3909 | 0.0312 | 189 | 7.63e-07 |
| COPZ1 | COPB2 | 1.0000 | 0.3837 | 0.0769 | 189 | 1.32e-06 |
| COPG2 | COPB1 | 1.0000 | 0.3826 | 0.1473 | 189 | 1.39e-06 |
| AP3B2 | AP2A1 | 0.6000 | 0.3818 | 0.1173 | 189 | 2.52e-20 |
| AP3B2 | AP2B1 | 0.6000 | 0.3775 | 0.8163 | 189 | 6.87e-20 |
| LRP1 | APP | 0.8000 | 0.3737 | 1.5741 | 189 | 7.31e-10 |
| COPA | COPB1 | 1.0000 | 0.3717 | 0.0735 | 189 | 3.25e-06 |
| VPS26A | VPS29 | 1.0000 | 0.3701 | 0.2434 | 189 | 3.59e-06 |
| SEC23A | SEC24C | 1.0000 | 0.3666 | 0.2899 | 189 | 4.64e-06 |
| TMED9 | TMED2 | 0.8000 | 0.3558 | 0.5243 | 189 | 6.25e-09 |
| TMED9 | TMED10 | 0.8000 | 0.3535 | 0.1846 | 189 | 7.14e-09 |
| COPB2 | TMED9 | 1.0000 | 0.3441 | 0.3730 | 189 | 2.20e-05 |
| COPA | COPG1 | 1.0000 | 0.3341 | 0.2553 | 189 | 4.30e-05 |
| VPS26A | VPS35 | 1.0000 | 0.3317 | 0.3999 | 189 | 4.98e-05 |
| COPG1 | TMED9 | 1.0000 | 0.3259 | 0.0820 | 189 | 7.21e-05 |
| COPB1 | TMED9 | 1.0000 | 0.3211 | 0.2392 | 189 | 9.62e-05 |
| AP4E1 | AP4B1 | 1.0000 | 0.3023 | 0.6918 | 78 | 0.0357 |
| COPA | TMED9 | 1.0000 | 0.2997 | 0.1347 | 189 | 3.38e-04 |
| TMED10 | COPB2 | 1.0000 | 0.2997 | 0.0221 | 189 | 3.38e-04 |
| SEC31A | SEC13 | 1.0000 | 0.2948 | 0.0081 | 189 | 4.38e-04 |
| AP2B1 | PICALM | 1.0000 | 0.2935 | 0.3876 | 189 | 4.68e-04 |
| CLU | APP | 0.6000 | 0.2913 | 0.5254 | 189 | 1.35e-10 |
| AP2B1 | AP4M1 | 0.8000 | 0.2888 | -1.6830 | 80 | 0.0070 |
| APOE | APP | 0.6000 | 0.2857 | 0.0213 | 189 | 3.44e-10 |
| COPE | TMED9 | 1.0000 | 0.2724 | 0.1679 | 189 | 0.0014 |
| AP1G1 | AP2A1 | 0.6000 | 0.2712 | -0.0976 | 189 | 3.35e-09 |
| CLU | APOE | 0.6000 | 0.2655 | 0.6083 | 189 | 7.14e-09 |
| RTN3 | RTN4 | 0.8000 | 0.2590 | 0.4475 | 189 | 8.18e-05 |
| TMED10 | COPE | 1.0000 | 0.2584 | 0.2756 | 189 | 0.0030 |
| ADAM10 | APP | 1.0000 | 0.2515 | 0.9090 | 189 | 0.0041 |
| HSPA5 | HSPA9 | 0.8000 | 0.2432 | 0.3130 | 189 | 2.71e-04 |

**BannerLFQ** (244 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| TMED10 | COPA | 1.0000 | 0.8006 | 0.0872 | 190 | 7.19e-41 |
| COPZ1 | COPA | 1.0000 | 0.7674 | 0.1205 | 150 | 2.53e-28 |
| TMED10 | COPG1 | 1.0000 | 0.7161 | 0.3567 | 190 | 5.77e-29 |
| COPB1 | COPA | 1.0000 | 0.7131 | 0.1694 | 190 | 1.07e-28 |
| TMED10 | COPZ1 | 1.0000 | 0.6813 | 0.3668 | 150 | 3.27e-20 |
| COPG1 | COPA | 1.0000 | 0.6552 | 0.0013 | 190 | 6.44e-23 |
| COPB1 | COPZ1 | 1.0000 | 0.5628 | 0.1255 | 150 | 9.63e-13 |
| COPE | COPB1 | 1.0000 | 0.5504 | 0.2966 | 190 | 4.05e-15 |
| TMED10 | COPB1 | 1.0000 | 0.5488 | 0.1993 | 190 | 4.98e-15 |
| COPE | COPA | 1.0000 | 0.5141 | 0.2176 | 190 | 5.13e-13 |
| COPZ1 | COPG1 | 1.0000 | 0.5138 | 0.0086 | 150 | 1.92e-10 |
| COPB2 | COPG1 | 1.0000 | 0.5038 | 0.2320 | 190 | 1.77e-12 |
| AP2A1 | AP2B1 | 1.0000 | 0.5002 | 0.5344 | 190 | 2.76e-12 |
| COPB2 | TMED10 | 1.0000 | 0.4903 | 1.0634 | 190 | 8.90e-12 |
| TMED9 | COPB2 | 1.0000 | 0.4784 | 0.2464 | 190 | 3.50e-11 |
| COPE | COPZ1 | 1.0000 | 0.4611 | 0.5880 | 150 | 2.45e-08 |
| COPG1 | COPB1 | 1.0000 | 0.4608 | 0.1383 | 190 | 2.31e-10 |
| TMED9 | COPA | 1.0000 | 0.4519 | 0.0051 | 190 | 5.83e-10 |
| TMED9 | COPG1 | 1.0000 | 0.4490 | -0.0027 | 190 | 7.82e-10 |
| SNX6 | SNX4 | 0.6000 | 0.4313 | 0.3571 | 190 | 3.62e-29 |
| TMED9 | TMED10 | 0.8000 | 0.4279 | 0.5959 | 190 | 3.62e-14 |
| VPS29 | SNX6 | 0.8000 | 0.4132 | 1.2287 | 190 | 3.85e-13 |
| AP1G1 | AP2A1 | 0.6000 | 0.4037 | 0.3003 | 190 | 1.71e-24 |
| TMED10 | SEC24B | 0.6000 | 0.4011 | 0.9936 | 164 | 6.31e-21 |
| RAB1B | COPA | 0.6000 | 0.3937 | -0.1331 | 190 | 6.44e-23 |
| COPB1 | RAB1B | 0.6000 | 0.3931 | 0.1130 | 190 | 6.44e-23 |
| HSPA5 | HSPA9 | 0.8000 | 0.3885 | 0.4639 | 190 | 1.51e-11 |
| SEC24B | RAB1B | 0.6000 | 0.3816 | 0.0893 | 164 | 2.07e-18 |
| PICALM | SORT1 | 0.6000 | 0.3781 | 0.3550 | 190 | 8.75e-21 |
| SEC13 | SEC23A | 0.8000 | 0.3640 | 0.3053 | 190 | 4.26e-10 |
| TMED10 | APP | 0.8000 | 0.3626 | 0.2330 | 190 | 5.12e-10 |
| AP2B1 | AP3B2 | 0.6000 | 0.3596 | 0.5898 | 190 | 2.26e-18 |
| SAR1B | TMED10 | 0.6000 | 0.3595 | 0.5582 | 174 | 6.79e-17 |
| COPB1 | TMED2 | 1.0000 | 0.3562 | -0.3183 | 190 | 2.92e-06 |
| COPB2 | COPA | 1.0000 | 0.3551 | 0.0586 | 190 | 3.09e-06 |
| PREB | SAR1B | 1.0000 | 0.3529 | 0.2337 | 114 | 5.28e-04 |
| APP | LRP1 | 0.8000 | 0.3458 | 0.2145 | 190 | 4.20e-09 |
| RAB1B | COPZ1 | 0.6000 | 0.3425 | 0.5179 | 150 | 3.85e-13 |
| VDAC1 | ADAM10 | 0.6000 | 0.3410 | -0.8810 | 189 | 3.53e-16 |
| VPS29 | SNX4 | 0.6000 | 0.3351 | 0.4043 | 190 | 1.23e-15 |

### 6.6 Selected edges across datasets

| pair | dataset | edge | prior_C | W | beta | n | q | supported |
|---|---|---|---|---|---|---|---|---|
| TMED2–TMED10 | ROSMAP | TMED10→TMED2 | 1.0000 | 0.3896 | 0.5420 | 400 | 1.65e-14 | 1 |
| TMED2–TMED10 | ROSMAP | TMED2→TMED10 | 1.0000 | 0.0000 | 0.0000 | 400 | 1.65e-14 | 0 |
| TMED2–TMED10 | Diverse | TMED10→TMED2 | 1.0000 | 0.6516 | 0.5780 | 676 | 3.95e-81 | 1 |
| TMED2–TMED10 | Diverse | TMED2→TMED10 | 1.0000 | 0.0000 | 0.0000 | 676 | 3.95e-81 | 0 |
| TMED2–TMED10 | Banner | TMED2→TMED10 | 1.0000 | 0.5143 | 0.6470 | 189 | 4.06e-12 | 1 |
| TMED2–TMED10 | Banner | TMED10→TMED2 | 1.0000 | 0.0000 | 0.0000 | 189 | 4.06e-12 | 0 |
| TMED2–TMED10 | BannerLFQ | TMED10→TMED2 | 1.0000 | 0.1182 | -0.0780 | 190 | 0.2355 | 0 |
| TMED2–TMED10 | BannerLFQ | TMED2→TMED10 | 1.0000 | 0.0000 | 0.0000 | 190 | 0.2355 | 0 |
| TMED9–TMED2 | ROSMAP | TMED2→TMED9 | 0.8000 | 0.1338 | -0.0100 | 400 | 0.0033 | 1 |
| TMED9–TMED2 | ROSMAP | TMED9→TMED2 | 0.8000 | 0.0000 | 0.0000 | 400 | 0.0033 | 0 |
| TMED9–TMED2 | Diverse | TMED9→TMED2 | 0.8000 | 0.4544 | 0.2119 | 676 | 1.75e-57 | 1 |
| TMED9–TMED2 | Diverse | TMED2→TMED9 | 0.8000 | 0.0000 | 0.0000 | 676 | 1.75e-57 | 0 |
| TMED9–TMED2 | Banner | TMED9→TMED2 | 0.8000 | 0.3558 | 0.5243 | 189 | 6.25e-09 | 1 |
| TMED9–TMED2 | Banner | TMED2→TMED9 | 0.8000 | 0.0000 | 0.0000 | 189 | 6.25e-09 | 0 |
| TMED9–TMED2 | BannerLFQ | TMED9→TMED2 | 0.8000 | 0.0475 | 0.0289 | 190 | 0.7974 | 0 |
| TMED9–TMED2 | BannerLFQ | TMED2→TMED9 | 0.8000 | 0.0000 | 0.0000 | 190 | 0.7974 | 0 |
| TMED9–TMED10 | ROSMAP | TMED10→TMED9 | 0.8000 | 0.4037 | 0.5392 | 400 | 2.97e-25 | 1 |
| TMED9–TMED10 | ROSMAP | TMED9→TMED10 | 0.8000 | 0.0000 | 0.0000 | 400 | 2.97e-25 | 0 |
| TMED9–TMED10 | Diverse | TMED9→TMED10 | 0.8000 | 0.5431 | 0.6757 | 676 | 1.82e-90 | 1 |
| TMED9–TMED10 | Diverse | TMED10→TMED9 | 0.8000 | 0.0000 | 0.0000 | 676 | 1.82e-90 | 0 |
| TMED9–TMED10 | Banner | TMED9→TMED10 | 0.8000 | 0.3535 | 0.1846 | 189 | 7.14e-09 | 1 |
| TMED9–TMED10 | Banner | TMED10→TMED9 | 0.8000 | 0.0000 | 0.0000 | 189 | 7.14e-09 | 0 |
| TMED9–TMED10 | BannerLFQ | TMED9→TMED10 | 0.8000 | 0.4279 | 0.5959 | 190 | 3.62e-14 | 1 |
| TMED9–TMED10 | BannerLFQ | TMED10→TMED9 | 0.8000 | 0.0000 | 0.0000 | 190 | 3.62e-14 | 0 |
| APP–BACE1 | ROSMAP | BACE1→APP | 1.0000 | 0.0603 | 0.0381 | 312 | 0.4427 | 0 |
| APP–BACE1 | ROSMAP | APP→BACE1 | 0.6000 | 0.0000 | 0.0000 | 312 | 0.4427 | 0 |
| APP–BACE1 | Diverse | APP→BACE1 | 0.6000 | 0.0376 | -0.0040 | 676 | 0.1326 | 0 |
| APP–BACE1 | Diverse | BACE1→APP | 1.0000 | 0.0000 | 0.0000 | 676 | 0.1326 | 0 |
| APP–BACE1 | Banner | BACE1→APP | 1.0000 | 0.2213 | 0.0605 | 181 | 0.0162 | 1 |
| APP–BACE1 | Banner | APP→BACE1 | 0.6000 | 0.0000 | 0.0000 | 181 | 0.0162 | 0 |
| APP–BACE1 | BannerLFQ | APP→BACE1 | 0.6000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| APP–BACE1 | BannerLFQ | BACE1→APP | 1.0000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| LRP1–APP | ROSMAP | LRP1→APP | 0.8000 | 0.0321 | 0.5125 | 312 | 0.6339 | 0 |
| LRP1–APP | ROSMAP | APP→LRP1 | 0.8000 | 0.0000 | 0.0000 | 312 | 0.6339 | 0 |
| LRP1–APP | Diverse | LRP1→APP | 0.8000 | 0.1548 | 2.1548 | 676 | 1.03e-06 | 1 |
| LRP1–APP | Diverse | APP→LRP1 | 0.8000 | 0.0000 | 0.0000 | 676 | 1.03e-06 | 0 |
| LRP1–APP | Banner | LRP1→APP | 0.8000 | 0.3737 | 1.5741 | 189 | 7.31e-10 | 1 |
| LRP1–APP | Banner | APP→LRP1 | 0.8000 | 0.0000 | 0.0000 | 189 | 7.31e-10 | 0 |
| LRP1–APP | BannerLFQ | APP→LRP1 | 0.8000 | 0.3458 | 0.2145 | 190 | 4.20e-09 | 1 |
| LRP1–APP | BannerLFQ | LRP1→APP | 0.8000 | 0.0000 | 0.0000 | 190 | 4.20e-09 | 0 |
| APOE–APP | ROSMAP | APOE→APP | 0.6000 | 0.0317 | -0.2157 | 312 | 0.5113 | 0 |
| APOE–APP | ROSMAP | APP→APOE | 0.2000 | 0.0000 | 0.0000 | 312 | 0.5113 | 0 |
| APOE–APP | Diverse | APOE→APP | 0.6000 | 0.2498 | 0.4941 | 676 | 1.32e-28 | 1 |
| APOE–APP | Diverse | APP→APOE | 0.2000 | 0.0000 | 0.0000 | 676 | 1.32e-28 | 0 |
| APOE–APP | Banner | APOE→APP | 0.6000 | 0.2857 | 0.0213 | 189 | 3.44e-10 | 1 |
| APOE–APP | Banner | APP→APOE | 0.2000 | 0.0000 | 0.0000 | 189 | 3.44e-10 | 0 |
| APOE–APP | BannerLFQ | APP→APOE | 0.2000 | 0.0807 | 0.1599 | 190 | 6.22e-08 | 1 |
| APOE–APP | BannerLFQ | APOE→APP | 0.6000 | 0.0000 | 0.0000 | 190 | 6.22e-08 | 0 |
| CLU–APOE | ROSMAP | CLU→APOE | 0.6000 | 0.1912 | 0.6887 | 400 | 1.08e-09 | 1 |
| CLU–APOE | ROSMAP | APOE→CLU | 0.6000 | 0.0000 | 0.0000 | 400 | 1.08e-09 | 0 |
| CLU–APOE | Diverse | CLU→APOE | 0.6000 | 0.2783 | 1.0294 | 676 | 3.80e-36 | 1 |
| CLU–APOE | Diverse | APOE→CLU | 0.6000 | 0.0000 | 0.0000 | 676 | 3.80e-36 | 0 |
| CLU–APOE | Banner | CLU→APOE | 0.6000 | 0.2655 | 0.6083 | 189 | 7.14e-09 | 1 |
| CLU–APOE | Banner | APOE→CLU | 0.6000 | 0.0000 | 0.0000 | 189 | 7.14e-09 | 0 |
| CLU–APOE | BannerLFQ | CLU→APOE | 0.6000 | 0.3117 | 0.8357 | 190 | 2.65e-13 | 1 |
| CLU–APOE | BannerLFQ | APOE→CLU | 0.6000 | 0.0000 | 0.0000 | 190 | 2.65e-13 | 0 |
| SORL1–APP | ROSMAP | APP→SORL1 | 0.6000 | 0.1523 | 0.0606 | 312 | 3.88e-05 | 1 |
| SORL1–APP | ROSMAP | SORL1→APP | 1.0000 | 0.0000 | 0.0000 | 312 | 3.88e-05 | 0 |
| SORL1–APP | Diverse | APP→SORL1 | 0.6000 | 0.0042 | 0.0397 | 676 | 0.8688 | 0 |
| SORL1–APP | Diverse | SORL1→APP | 1.0000 | 0.0000 | 0.0000 | 676 | 0.8688 | 0 |
| SORL1–APP | Banner | APP→SORL1 | 0.6000 | 0.0167 | -0.0481 | 189 | 0.8284 | 0 |
| SORL1–APP | Banner | SORL1→APP | 1.0000 | 0.0000 | 0.0000 | 189 | 0.8284 | 0 |
| SORL1–APP | BannerLFQ | SORL1→APP | 1.0000 | 0.2771 | 0.1948 | 170 | 0.0011 | 1 |
| SORL1–APP | BannerLFQ | APP→SORL1 | 0.6000 | 0.0000 | 0.0000 | 170 | 0.0011 | 0 |
| PSEN1–NCSTN | ROSMAP | NCSTN→PSEN1 | 1.0000 | 0.1352 | 0.0457 | 384 | 0.0241 | 1 |
| PSEN1–NCSTN | ROSMAP | PSEN1→NCSTN | 1.0000 | 0.0000 | 0.0000 | 384 | 0.0241 | 0 |
| PSEN1–NCSTN | Diverse | PSEN1→NCSTN | 1.0000 | 0.3282 | 0.3095 | 676 | 1.07e-17 | 1 |
| PSEN1–NCSTN | Diverse | NCSTN→PSEN1 | 1.0000 | 0.0000 | 0.0000 | 676 | 1.07e-17 | 0 |
| PSEN1–NCSTN | Banner | NCSTN→PSEN1 | 1.0000 | 0.1378 | 0.0052 | 189 | 0.1783 | 0 |
| PSEN1–NCSTN | Banner | PSEN1→NCSTN | 1.0000 | 0.0000 | 0.0000 | 189 | 0.1783 | 0 |
| PSEN1–NCSTN | BannerLFQ | PSEN1→NCSTN | 1.0000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| PSEN1–NCSTN | BannerLFQ | NCSTN→PSEN1 | 1.0000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| BIN1–BACE1 | ROSMAP | BIN1→BACE1 | 0.8000 | 0.0122 | 0.1002 | 400 | 0.8442 | 0 |
| BIN1–BACE1 | ROSMAP | BACE1→BIN1 | 0.4000 | 0.0000 | 0.0000 | 400 | 0.8442 | 0 |
| BIN1–BACE1 | Diverse | BIN1→BACE1 | 0.8000 | 0.0777 | 0.0185 | 676 | 0.0175 | 1 |
| BIN1–BACE1 | Diverse | BACE1→BIN1 | 0.4000 | 0.0000 | 0.0000 | 676 | 0.0175 | 0 |
| BIN1–BACE1 | Banner | BIN1→BACE1 | 0.8000 | 0.1777 | 0.2948 | 181 | 0.0158 | 1 |
| BIN1–BACE1 | Banner | BACE1→BIN1 | 0.4000 | 0.0000 | 0.0000 | 181 | 0.0158 | 0 |
| BIN1–BACE1 | BannerLFQ | BACE1→BIN1 | 0.4000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| BIN1–BACE1 | BannerLFQ | BIN1→BACE1 | 0.8000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| PICALM–APP | ROSMAP | PICALM→APP | 0.8000 | 0.0431 | 0.4084 | 312 | 0.5001 | 0 |
| PICALM–APP | ROSMAP | APP→PICALM | 0.4000 | 0.0000 | 0.0000 | 312 | 0.5001 | 0 |
| PICALM–APP | Diverse | APP→PICALM | 0.4000 | 0.0331 | 0.0125 | 676 | 0.0447 | 1 |
| PICALM–APP | Diverse | PICALM→APP | 0.8000 | 0.0000 | 0.0000 | 676 | 0.0447 | 0 |
| PICALM–APP | Banner | PICALM→APP | 0.8000 | 0.1266 | 1.7762 | 189 | 0.1064 | 0 |
| PICALM–APP | Banner | APP→PICALM | 0.4000 | 0.0000 | 0.0000 | 189 | 0.1064 | 0 |
| PICALM–APP | BannerLFQ | APP→PICALM | 0.4000 | 0.1217 | 0.2369 | 190 | 9.72e-05 | 1 |
| PICALM–APP | BannerLFQ | PICALM→APP | 0.8000 | 0.0000 | 0.0000 | 190 | 9.72e-05 | 0 |

---

## 7. Agreement between the three graphs (HIW 5)

### 7.1 Pairwise overlap of supported pairs vs chance

| pair | jaccard | shared | expected_by_chance | fold | p_hypergeom | same_direction |
|---|---|---|---|---|---|---|
| ROSMAP-Diverse | 0.4187 | 211 | 178.0795 | 1.1849 | 1.83e-09 | 0.5592 |
| ROSMAP-Banner | 0.2961 | 90 | 54.8823 | 1.6399 | 1.37e-11 | 0.4556 |
| Diverse-Banner | 0.2620 | 126 | 101.4557 | 1.2419 | 5.12e-08 | 0.5556 |

### 7.2 Consensus pairs

- pairs supported in ≥ 2 datasets: **257** (in exactly 2: 172; in all 3: **85**)
- in exactly 2: ROSMAP+Diverse 126, Diverse+Banner 41, ROSMAP+Banner 5
- arrow identical in every dataset that supports the pair: 119 of 257; among the 85 triple pairs: 25 of 85
- split by prior type — one-way C pairs (arrow fixed by C): 0 of 257, 0 consistent; two-way C pairs (arrow picked by the direction test): 257, of which 119 consistent and 138 not
- among the 85 triple pairs: 0 one-way (0 consistent); 85 two-way, of which 25 consistent and 60 not

**All 85 pairs supported in all three datasets**

| protein_a | protein_b | prior_C_a_to_b | prior_C_b_to_a | dir_ROSMAP | W_ROSMAP | dir_Diverse | W_Diverse | dir_Banner | W_Banner | direction_consistent |
|---|---|---|---|---|---|---|---|---|---|---|
| AP2A1 | AP2B1 | 1.000 | 1.000 | AP2A1->AP2B1 | 0.732 | AP2B1->AP2A1 | 0.801 | AP2B1->AP2A1 | 0.843 | 0 |
| COPB1 | COPG1 | 1.000 | 1.000 | COPB1->COPG1 | 0.628 | COPB1->COPG1 | 0.662 | COPG1->COPB1 | 0.406 | 0 |
| COPB1 | COPB2 | 1.000 | 1.000 | COPB2->COPB1 | 0.493 | COPB1->COPB2 | 0.639 | COPB2->COPB1 | 0.531 | 0 |
| TMED10 | TMED2 | 1.000 | 1.000 | TMED10->TMED2 | 0.390 | TMED10->TMED2 | 0.652 | TMED2->TMED10 | 0.514 | 0 |
| SEC23A | SEC24C | 1.000 | 1.000 | SEC24C->SEC23A | 0.501 | SEC24C->SEC23A | 0.668 | SEC23A->SEC24C | 0.367 | 0 |
| COPB1 | COPE | 1.000 | 1.000 | COPB1->COPE | 0.509 | COPB1->COPE | 0.499 | COPE->COPB1 | 0.448 | 0 |
| COPA | COPB1 | 1.000 | 1.000 | COPB1->COPA | 0.340 | COPA->COPB1 | 0.712 | COPA->COPB1 | 0.372 | 0 |
| SEC13 | SEC31A | 1.000 | 1.000 | SEC31A->SEC13 | 0.437 | SEC31A->SEC13 | 0.645 | SEC31A->SEC13 | 0.295 | 1 |
| COPB2 | COPE | 1.000 | 1.000 | COPB2->COPE | 0.468 | COPB2->COPE | 0.507 | COPE->COPB2 | 0.401 | 0 |
| COPA | COPG1 | 1.000 | 1.000 | COPG1->COPA | 0.371 | COPA->COPG1 | 0.664 | COPA->COPG1 | 0.334 | 0 |
| COPB2 | COPG1 | 1.000 | 1.000 | COPB2->COPG1 | 0.413 | COPG1->COPB2 | 0.496 | COPG1->COPB2 | 0.456 | 0 |
| VPS29 | VPS35 | 1.000 | 1.000 | VPS35->VPS29 | 0.319 | VPS35->VPS29 | 0.584 | VPS35->VPS29 | 0.442 | 1 |
| TMED10 | TMED9 | 0.800 | 0.800 | TMED10->TMED9 | 0.404 | TMED9->TMED10 | 0.543 | TMED9->TMED10 | 0.353 | 0 |
| COPB2 | COPZ1 | 1.000 | 1.000 | COPZ1->COPB2 | 0.341 | COPB2->COPZ1 | 0.531 | COPZ1->COPB2 | 0.384 | 0 |
| COPA | COPB2 | 1.000 | 1.000 | COPB2->COPA | 0.221 | COPB2->COPA | 0.636 | COPA->COPB2 | 0.391 | 0 |
| COPE | COPG1 | 1.000 | 1.000 | COPE->COPG1 | 0.455 | COPE->COPG1 | 0.388 | COPE->COPG1 | 0.403 | 1 |
| VPS26A | VPS29 | 1.000 | 1.000 | VPS26A->VPS29 | 0.192 | VPS26A->VPS29 | 0.595 | VPS26A->VPS29 | 0.370 | 1 |
| COPB2 | TMED10 | 1.000 | 1.000 | COPB2->TMED10 | 0.383 | COPB2->TMED10 | 0.406 | TMED10->COPB2 | 0.300 | 0 |
| COPB1 | COPG2 | 1.000 | 1.000 | COPB1->COPG2 | 0.273 | COPB1->COPG2 | 0.394 | COPG2->COPB1 | 0.383 | 0 |
| COPB2 | TMED9 | 1.000 | 1.000 | COPB2->TMED9 | 0.273 | COPB2->TMED9 | 0.420 | COPB2->TMED9 | 0.344 | 1 |
| AP1G1 | AP2A1 | 0.600 | 0.600 | AP2A1->AP1G1 | 0.328 | AP2A1->AP1G1 | 0.414 | AP1G1->AP2A1 | 0.271 | 0 |
| COPE | COPZ1 | 1.000 | 1.000 | COPZ1->COPE | 0.355 | COPZ1->COPE | 0.416 | COPZ1->COPE | 0.218 | 1 |
| COPB1 | TMED10 | 1.000 | 1.000 | COPB1->TMED10 | 0.398 | COPB1->TMED10 | 0.314 | TMED10->COPB1 | 0.237 | 0 |
| TMED2 | TMED9 | 0.800 | 0.800 | TMED2->TMED9 | 0.134 | TMED9->TMED2 | 0.454 | TMED9->TMED2 | 0.356 | 0 |
| COPE | TMED10 | 1.000 | 1.000 | COPE->TMED10 | 0.282 | COPE->TMED10 | 0.363 | TMED10->COPE | 0.258 | 0 |
| SEC23A | SEC24B | 1.000 | 1.000 | SEC24B->SEC23A | 0.202 | SEC23A->SEC24B | 0.487 | SEC23A->SEC24B | 0.200 | 0 |
| AP2B1 | PICALM | 1.000 | 1.000 | PICALM->AP2B1 | 0.298 | AP2B1->PICALM | 0.279 | AP2B1->PICALM | 0.293 | 0 |
| HSPA5 | HSPA9 | 0.800 | 0.800 | HSPA5->HSPA9 | 0.326 | HSPA5->HSPA9 | 0.291 | HSPA5->HSPA9 | 0.243 | 1 |
| COPA | COPG2 | 1.000 | 1.000 | COPA->COPG2 | 0.245 | COPA->COPG2 | 0.396 | COPG2->COPA | 0.206 | 0 |
| HSPA9 | MFN2 | 0.600 | 0.600 | HSPA9->MFN2 | 0.256 | MFN2->HSPA9 | 0.325 | HSPA9->MFN2 | 0.231 | 0 |
| AP2A1 | PICALM | 1.000 | 1.000 | AP2A1->PICALM | 0.295 | AP2A1->PICALM | 0.319 | AP2A1->PICALM | 0.191 | 1 |
| HSPA5 | TMED10 | 0.800 | 0.400 | HSPA5->TMED10 | 0.305 | HSPA5->TMED10 | 0.275 | HSPA5->TMED10 | 0.193 | 1 |
| SEC13 | SEC23A | 0.800 | 0.800 | SEC13->SEC23A | 0.204 | SEC23A->SEC13 | 0.305 | SEC23A->SEC13 | 0.233 | 0 |
| SEC24C | SEC31A | 0.600 | 0.600 | SEC24C->SEC31A | 0.288 | SEC24C->SEC31A | 0.305 | SEC31A->SEC24C | 0.147 | 0 |
| APOE | CLU | 0.600 | 0.600 | CLU->APOE | 0.191 | CLU->APOE | 0.278 | CLU->APOE | 0.266 | 1 |
| AP1G1 | AP2B1 | 0.400 | 0.400 | AP2B1->AP1G1 | 0.205 | AP2B1->AP1G1 | 0.286 | AP2B1->AP1G1 | 0.242 | 1 |
| EIF2S1 | HSPA5 | 0.600 | 0.600 | EIF2S1->HSPA5 | 0.260 | EIF2S1->HSPA5 | 0.272 | EIF2S1->HSPA5 | 0.181 | 1 |
| RTN3 | RTN4 | 0.800 | 0.800 | RTN3->RTN4 | 0.122 | RTN4->RTN3 | 0.323 | RTN3->RTN4 | 0.259 | 0 |
| HSPA9 | VDAC1 | 0.800 | 0.800 | HSPA9->VDAC1 | 0.112 | VDAC1->HSPA9 | 0.350 | VDAC1->HSPA9 | 0.207 | 0 |
| COPG1 | TMED9 | 1.000 | 1.000 | COPG1->TMED9 | 0.129 | COPG1->TMED9 | 0.190 | COPG1->TMED9 | 0.326 | 1 |
| LRP1 | SORL1 | 0.600 | 0.600 | SORL1->LRP1 | 0.152 | SORL1->LRP1 | 0.267 | LRP1->SORL1 | 0.217 | 0 |
| SEC13 | SEC24C | 0.600 | 0.600 | SEC13->SEC24C | 0.213 | SEC13->SEC24C | 0.256 | SEC13->SEC24C | 0.145 | 1 |
| APP | CLU | 0.600 | 0.600 | APP->CLU | 0.094 | CLU->APP | 0.202 | CLU->APP | 0.291 | 0 |
| RAB1A | RAB1B | 0.400 | 0.400 | RAB1A->RAB1B | 0.179 | RAB1A->RAB1B | 0.198 | RAB1B->RAB1A | 0.190 | 0 |
| MIA3 | PREB | 0.600 | 0.600 | MIA3->PREB | 0.162 | PREB->MIA3 | 0.228 | PREB->MIA3 | 0.151 | 0 |
| HSPA5 | MIA3 | 0.400 | 0.400 | MIA3->HSPA5 | 0.192 | MIA3->HSPA5 | 0.142 | HSPA5->MIA3 | 0.175 | 0 |
| GGA2 | SORL1 | 0.800 | 0.800 | GGA2->SORL1 | 0.114 | SORL1->GGA2 | 0.156 | SORL1->GGA2 | 0.221 | 0 |
| BACE1 | SORL1 | 0.400 | 0.800 | SORL1->BACE1 | 0.167 | SORL1->BACE1 | 0.153 | SORL1->BACE1 | 0.164 | 1 |
| RAB1B | TMED2 | 0.600 | 0.600 | RAB1B->TMED2 | 0.182 | RAB1B->TMED2 | 0.162 | RAB1B->TMED2 | 0.114 | 1 |
| SEC24C | TMED10 | 0.600 | 0.600 | SEC24C->TMED10 | 0.133 | TMED10->SEC24C | 0.193 | SEC24C->TMED10 | 0.119 | 0 |
| ADAM10 | VDAC1 | 0.600 | 0.600 | ADAM10->VDAC1 | 0.091 | VDAC1->ADAM10 | 0.140 | ADAM10->VDAC1 | 0.188 | 0 |
| AP1G1 | PICALM | 0.600 | 0.600 | AP1G1->PICALM | 0.129 | AP1G1->PICALM | 0.111 | AP1G1->PICALM | 0.178 | 1 |
| AP2B1 | SEC24C | 0.400 | 0.400 | SEC24C->AP2B1 | 0.140 | AP2B1->SEC24C | 0.138 | AP2B1->SEC24C | 0.113 | 0 |
| AP2A1 | ARF6 | 0.600 | 0.600 | ARF6->AP2A1 | 0.080 | AP2A1->ARF6 | 0.156 | AP2A1->ARF6 | 0.121 | 0 |
| SEC24B | TMED2 | 0.600 | 0.600 | SEC24B->TMED2 | 0.096 | TMED2->SEC24B | 0.134 | SEC24B->TMED2 | 0.125 | 0 |
| AP1G1 | VPS35 | 0.200 | 0.200 | AP1G1->VPS35 | 0.132 | VPS35->AP1G1 | 0.120 | AP1G1->VPS35 | 0.095 | 0 |
| CLU | SORL1 | 0.400 | 0.400 | CLU->SORL1 | 0.066 | CLU->SORL1 | 0.169 | SORL1->CLU | 0.108 | 0 |
| AP2B1 | GGA3 | 0.400 | 0.400 | GGA3->AP2B1 | 0.053 | AP2B1->GGA3 | 0.192 | AP2B1->GGA3 | 0.075 | 0 |
| BIN1 | SNX4 | 0.400 | 0.400 | SNX4->BIN1 | 0.047 | SNX4->BIN1 | 0.138 | BIN1->SNX4 | 0.135 | 0 |
| SAR1B | SEC24B | 0.600 | 0.600 | SAR1B->SEC24B | 0.082 | SAR1B->SEC24B | 0.097 | SEC24B->SAR1B | 0.124 | 0 |
| ITPR2 | MFN2 | 0.400 | 0.400 | MFN2->ITPR2 | 0.071 | MFN2->ITPR2 | 0.101 | MFN2->ITPR2 | 0.127 | 1 |
| AP1G1 | SORT1 | 0.600 | 0.600 | AP1G1->SORT1 | 0.081 | AP1G1->SORT1 | 0.084 | AP1G1->SORT1 | 0.132 | 1 |
| COPA | HSPA5 | 0.400 | 0.400 | COPA->HSPA5 | 0.050 | HSPA5->COPA | 0.123 | HSPA5->COPA | 0.107 | 0 |
| APP | VPS26A | 0.200 | 0.600 | APP->VPS26A | 0.030 | VPS26A->APP | 0.071 | VPS26A->APP | 0.147 | 0 |
| AP2A1 | BIN1 | 0.400 | 0.400 | AP2A1->BIN1 | 0.055 | AP2A1->BIN1 | 0.095 | BIN1->AP2A1 | 0.090 | 0 |
| COPB1 | SEC31A | 0.200 | 0.200 | COPB1->SEC31A | 0.086 | SEC31A->COPB1 | 0.100 | SEC31A->COPB1 | 0.039 | 0 |
| COPA | SEC31A | 0.200 | 0.200 | SEC31A->COPA | 0.042 | SEC31A->COPA | 0.110 | SEC31A->COPA | 0.069 | 1 |
| COPB1 | SEC13 | 0.200 | 0.200 | COPB1->SEC13 | 0.069 | COPB1->SEC13 | 0.080 | SEC13->COPB1 | 0.055 | 0 |
| COPA | SEC13 | 0.200 | 0.200 | COPA->SEC13 | 0.040 | COPA->SEC13 | 0.089 | SEC13->COPA | 0.063 | 0 |
| COPG1 | SEC13 | 0.200 | 0.200 | SEC13->COPG1 | 0.056 | SEC13->COPG1 | 0.063 | COPG1->SEC13 | 0.071 | 0 |
| BACE1 | SORT1 | 0.200 | 0.600 | BACE1->SORT1 | 0.024 | SORT1->BACE1 | 0.126 | BACE1->SORT1 | 0.039 | 0 |
| COPB1 | HSPA5 | 0.200 | 0.200 | HSPA5->COPB1 | 0.039 | HSPA5->COPB1 | 0.067 | HSPA5->COPB1 | 0.076 | 1 |
| BIN1 | SNX6 | 0.200 | 0.200 | SNX6->BIN1 | 0.051 | SNX6->BIN1 | 0.065 | SNX6->BIN1 | 0.057 | 1 |
| SAR1A | TMED10 | 0.200 | 0.200 | SAR1A->TMED10 | 0.050 | SAR1A->TMED10 | 0.077 | SAR1A->TMED10 | 0.045 | 1 |
| COPG2 | SEC13 | 0.200 | 0.200 | SEC13->COPG2 | 0.052 | SEC13->COPG2 | 0.045 | COPG2->SEC13 | 0.064 | 0 |
| COPE | PREB | 0.200 | 0.200 | COPE->PREB | 0.051 | COPE->PREB | 0.047 | PREB->COPE | 0.059 | 0 |
| SAR1A | TMED2 | 0.200 | 0.200 | TMED2->SAR1A | 0.034 | SAR1A->TMED2 | 0.076 | TMED2->SAR1A | 0.047 | 0 |
| COPE | HSPA5 | 0.200 | 0.200 | COPE->HSPA5 | 0.038 | HSPA5->COPE | 0.070 | HSPA5->COPE | 0.048 | 0 |
| COPB1 | MIA3 | 0.200 | 0.200 | MIA3->COPB1 | 0.045 | COPB1->MIA3 | 0.065 | MIA3->COPB1 | 0.043 | 0 |
| COPB1 | SEC23A | 0.200 | 0.200 | SEC23A->COPB1 | 0.047 | SEC23A->COPB1 | 0.064 | SEC23A->COPB1 | 0.042 | 1 |
| AP2B1 | BIN1 | 0.200 | 0.200 | BIN1->AP2B1 | 0.051 | AP2B1->BIN1 | 0.060 | BIN1->AP2B1 | 0.041 | 0 |
| HSPA5 | SAR1B | 0.200 | 0.200 | HSPA5->SAR1B | 0.032 | HSPA5->SAR1B | 0.064 | SAR1B->HSPA5 | 0.050 | 0 |
| COPG1 | SEC23A | 0.200 | 0.200 | COPG1->SEC23A | 0.024 | SEC23A->COPG1 | 0.050 | SEC23A->COPG1 | 0.066 | 0 |
| BACE1 | CLU | 0.200 | 0.200 | BACE1->CLU | 0.035 | BACE1->CLU | 0.038 | BACE1->CLU | 0.057 | 1 |
| APP | RTN4 | 0.200 | 0.200 | APP->RTN4 | 0.038 | APP->RTN4 | 0.026 | RTN4->APP | 0.047 | 0 |

**Pairs supported in exactly two datasets** (172)

| protein_a | protein_b | dir_ROSMAP | W_ROSMAP | dir_Diverse | W_Diverse | dir_Banner | W_Banner |
|---|---|---|---|---|---|---|---|
| COPB1 | COPZ1 | COPZ1->COPB1 | 0.473 | COPZ1->COPB1 | 0.492 | — | — |
| AP2A1 | AP3B2 | — | — | AP2A1->AP3B2 | 0.461 | AP3B2->AP2A1 | 0.382 |
| VPS26A | VPS35 | — | — | VPS26A->VPS35 | 0.471 | VPS26A->VPS35 | 0.332 |
| AP2B1 | AP3B2 | — | — | AP2B1->AP3B2 | 0.374 | AP3B2->AP2B1 | 0.378 |
| COPZ1 | TMED10 | COPZ1->TMED10 | 0.356 | COPZ1->TMED10 | 0.360 | — | — |
| SEC23A | SEC31A | SEC31A->SEC23A | 0.282 | SEC31A->SEC23A | 0.400 | — | — |
| COPA | TMED9 | — | — | COPA->TMED9 | 0.353 | COPA->TMED9 | 0.300 |
| COPE | TMED9 | — | — | COPE->TMED9 | 0.370 | COPE->TMED9 | 0.272 |
| COPB1 | TMED9 | — | — | COPB1->TMED9 | 0.316 | COPB1->TMED9 | 0.321 |
| COPA | COPE | COPA->COPE | 0.119 | COPE->COPA | 0.518 | — | — |
| PLCG2 | SYK | PLCG2->SYK | 0.246 | SYK->PLCG2 | 0.380 | — | — |
| COPB2 | TMED2 | COPB2->TMED2 | 0.310 | COPB2->TMED2 | 0.305 | — | — |
| COPG1 | COPZ1 | COPZ1->COPG1 | 0.319 | COPG1->COPZ1 | 0.266 | — | — |
| COPA | COPZ1 | COPA->COPZ1 | 0.178 | COPA->COPZ1 | 0.400 | — | — |
| AP1G1 | AP3B2 | — | — | AP3B2->AP1G1 | 0.355 | AP3B2->AP1G1 | 0.201 |
| COPB2 | COPG2 | COPB2->COPG2 | 0.264 | COPB2->COPG2 | 0.288 | — | — |
| APOE | APP | — | — | APOE->APP | 0.250 | APOE->APP | 0.286 |
| APP | LRP1 | — | — | LRP1->APP | 0.155 | LRP1->APP | 0.374 |
| COPB1 | TMED2 | COPB1->TMED2 | 0.209 | COPB1->TMED2 | 0.278 | — | — |
| AP4B1 | AP4E1 | — | — | AP4B1->AP4E1 | 0.176 | AP4E1->AP4B1 | 0.302 |
| COPE | TMED2 | TMED2->COPE | 0.185 | COPE->TMED2 | 0.280 | — | — |
| NCSTN | PSEN1 | NCSTN->PSEN1 | 0.135 | PSEN1->NCSTN | 0.328 | — | — |
| AP2B1 | AP4M1 | AP2B1->AP4M1 | 0.170 | — | — | AP2B1->AP4M1 | 0.289 |
| COPZ1 | TMED9 | COPZ1->TMED9 | 0.120 | COPZ1->TMED9 | 0.308 | — | — |
| COPE | COPG2 | COPE->COPG2 | 0.251 | COPE->COPG2 | 0.167 | — | — |
| SEC24B | SEC31A | SEC24B->SEC31A | 0.110 | SEC31A->SEC24B | 0.284 | — | — |
| AP4E1 | AP4M1 | AP4E1->AP4M1 | 0.209 | AP4E1->AP4M1 | 0.178 | — | — |
| COPG1 | TMED10 | COPG1->TMED10 | 0.205 | COPG1->TMED10 | 0.165 | — | — |
| GGA1 | GGA3 | GGA1->GGA3 | 0.116 | GGA3->GGA1 | 0.237 | — | — |
| ADAM10 | APP | — | — | ADAM10->APP | 0.100 | ADAM10->APP | 0.252 |
| SAR1A | SEC24C | SEC24C->SAR1A | 0.118 | SAR1A->SEC24C | 0.228 | — | — |
| MAPK14 | MAPKAPK2 | MAPK14->MAPKAPK2 | 0.134 | MAPKAPK2->MAPK14 | 0.207 | — | — |
| MIA3 | SEC24C | — | — | SEC24C->MIA3 | 0.207 | SEC24C->MIA3 | 0.133 |
| SAR1A | SEC31A | SAR1A->SEC31A | 0.170 | SEC31A->SAR1A | 0.167 | — | — |
| APP | VPS35 | — | — | VPS35->APP | 0.117 | VPS35->APP | 0.219 |
| APBA1 | APBA2 | — | — | APBA1->APBA2 | 0.152 | APBA1->APBA2 | 0.167 |
| SEC13 | TMED10 | SEC13->TMED10 | 0.171 | SEC13->TMED10 | 0.148 | — | — |
| SEC23B | SEC24A | SEC24A->SEC23B | 0.149 | SEC23B->SEC24A | 0.163 | — | — |
| SEC24C | TMED2 | — | — | TMED2->SEC24C | 0.185 | SEC24C->TMED2 | 0.125 |
| AP4E1 | AP4S1 | AP4S1->AP4E1 | 0.162 | AP4S1->AP4E1 | 0.148 | — | — |
| COPG1 | RAB1A | COPG1->RAB1A | 0.090 | — | — | COPG1->RAB1A | 0.219 |
| SEC23A | TMED2 | SEC23A->TMED2 | 0.131 | SEC23A->TMED2 | 0.178 | — | — |
| ITPR1 | MFN2 | — | — | MFN2->ITPR1 | 0.164 | MFN2->ITPR1 | 0.145 |
| AP1G1 | AP4M1 | — | — | AP4M1->AP1G1 | 0.079 | AP1G1->AP4M1 | 0.229 |
| MAP2K3 | MAPKAPK2 | MAPKAPK2->MAP2K3 | 0.184 | MAP2K3->MAPKAPK2 | 0.122 | — | — |
| AP2B1 | BACE1 | — | — | BACE1->AP2B1 | 0.178 | BACE1->AP2B1 | 0.118 |
| RAB1A | SEC23A | — | — | RAB1A->SEC23A | 0.144 | SEC23A->RAB1A | 0.148 |
| SAR1A | SEC13 | SEC13->SAR1A | 0.122 | SAR1A->SEC13 | 0.170 | — | — |
| AP1G1 | GGA1 | AP1G1->GGA1 | 0.103 | AP1G1->GGA1 | 0.186 | — | — |
| SEC31A | TMED10 | SEC31A->TMED10 | 0.164 | SEC31A->TMED10 | 0.122 | — | — |
| PICALM | SEC24C | SEC24C->PICALM | 0.130 | SEC24C->PICALM | 0.151 | — | — |
| SEC24A | TMED10 | SEC24A->TMED10 | 0.100 | SEC24A->TMED10 | 0.177 | — | — |
| AP2B1 | ARF6 | — | — | AP2B1->ARF6 | 0.125 | AP2B1->ARF6 | 0.151 |
| SEC23A | TMED10 | SEC23A->TMED10 | 0.123 | SEC23A->TMED10 | 0.145 | — | — |
| RAB1A | SEC24C | — | — | SEC24C->RAB1A | 0.128 | SEC24C->RAB1A | 0.138 |
| SAR1B | SEC23B | SEC23B->SAR1B | 0.120 | SAR1B->SEC23B | 0.146 | — | — |
| PREB | TMED10 | PREB->TMED10 | 0.140 | PREB->TMED10 | 0.121 | — | — |
| MFN2 | TMED10 | MFN2->TMED10 | 0.124 | MFN2->TMED10 | 0.134 | — | — |
| SEC24A | SEC31A | SEC31A->SEC24A | 0.082 | SEC31A->SEC24A | 0.175 | — | — |
| BACE1 | BIN1 | — | — | BIN1->BACE1 | 0.078 | BIN1->BACE1 | 0.178 |
| EIF2S1 | TMED9 | EIF2S1->TMED9 | 0.097 | EIF2S1->TMED9 | 0.154 | — | — |
| COPA | SEC23A | — | — | SEC23A->COPA | 0.156 | COPA->SEC23A | 0.091 |
| PREB | SEC31A | PREB->SEC31A | 0.047 | PREB->SEC31A | 0.198 | — | — |
| COPZ1 | RAB1B | — | — | COPZ1->RAB1B | 0.071 | RAB1B->COPZ1 | 0.169 |
| RAB1A | TMED9 | — | — | TMED9->RAB1A | 0.105 | RAB1A->TMED9 | 0.133 |
| COPA | SEC24C | COPA->SEC24C | 0.083 | SEC24C->COPA | 0.152 | — | — |
| SNX6 | VPS29 | — | — | VPS29->SNX6 | 0.083 | SNX6->VPS29 | 0.150 |
| SEC23A | TMED9 | — | — | TMED9->SEC23A | 0.122 | SEC23A->TMED9 | 0.107 |
| SEC24A | TMED2 | TMED2->SEC24A | 0.070 | SEC24A->TMED2 | 0.154 | — | — |
| SEC13 | SEC24A | SEC24A->SEC13 | 0.103 | SEC24A->SEC13 | 0.120 | — | — |
| PREB | SEC13 | SEC13->PREB | 0.058 | PREB->SEC13 | 0.164 | — | — |
| MIA3 | SEC24B | MIA3->SEC24B | 0.070 | MIA3->SEC24B | 0.147 | — | — |
| RAB1A | VDAC1 | RAB1A->VDAC1 | 0.100 | RAB1A->VDAC1 | 0.113 | — | — |
| AP2B1 | AP4E1 | AP2B1->AP4E1 | 0.103 | AP2B1->AP4E1 | 0.106 | — | — |
| GORASP1 | TMED2 | GORASP1->TMED2 | 0.134 | GORASP1->TMED2 | 0.073 | — | — |
| GORASP1 | RAB1B | RAB1B->GORASP1 | 0.146 | GORASP1->RAB1B | 0.058 | — | — |
| GGA2 | LRP1 | LRP1->GGA2 | 0.080 | LRP1->GGA2 | 0.124 | — | — |
| SAR1B | TMED2 | TMED2->SAR1B | 0.093 | SAR1B->TMED2 | 0.110 | — | — |
| RAB1B | SEC24C | — | — | RAB1B->SEC24C | 0.059 | SEC24C->RAB1B | 0.141 |
| RAB1A | TMED2 | — | — | TMED2->RAB1A | 0.051 | RAB1A->TMED2 | 0.146 |
| GGA2 | PICALM | GGA2->PICALM | 0.054 | GGA2->PICALM | 0.142 | — | — |
| AP3B2 | AP4E1 | AP4E1->AP3B2 | 0.119 | AP3B2->AP4E1 | 0.075 | — | — |
| RAB1A | SEC13 | — | — | SEC13->RAB1A | 0.071 | SEC13->RAB1A | 0.123 |
| APP | GGA2 | — | — | APP->GGA2 | 0.073 | GGA2->APP | 0.120 |
| ITPR1 | SYK | ITPR1->SYK | 0.094 | SYK->ITPR1 | 0.099 | — | — |
| BIN1 | PICALM | BIN1->PICALM | 0.059 | PICALM->BIN1 | 0.131 | — | — |
| RAB1A | RTN4 | RTN4->RAB1A | 0.114 | — | — | RAB1A->RTN4 | 0.074 |
| COPG1 | RAB1B | RAB1B->COPG1 | 0.089 | COPG1->RAB1B | 0.098 | — | — |
| SEC24D | SEC31A | SEC31A->SEC24D | 0.091 | SEC31A->SEC24D | 0.090 | — | — |
| SAR1A | SAR1B | SAR1B->SAR1A | 0.079 | SAR1A->SAR1B | 0.102 | — | — |
| COPA | RTN4 | COPA->RTN4 | 0.071 | COPA->RTN4 | 0.109 | — | — |
| ARF6 | COPB2 | ARF6->COPB2 | 0.089 | COPB2->ARF6 | 0.090 | — | — |
| HSPA9 | RAB1A | HSPA9->RAB1A | 0.070 | RAB1A->HSPA9 | 0.107 | — | — |
| PICALM | SORL1 | PICALM->SORL1 | 0.088 | PICALM->SORL1 | 0.088 | — | — |
| GORASP1 | TMED10 | GORASP1->TMED10 | 0.113 | GORASP1->TMED10 | 0.062 | — | — |
| SNX6 | VPS35 | SNX6->VPS35 | 0.095 | VPS35->SNX6 | 0.075 | — | — |
| COPG2 | RAB1B | RAB1B->COPG2 | 0.114 | RAB1B->COPG2 | 0.056 | — | — |
| AP2B1 | GGA2 | GGA2->AP2B1 | 0.075 | AP2B1->GGA2 | 0.092 | — | — |
| BACE1 | NCSTN | — | — | BACE1->NCSTN | 0.052 | BACE1->NCSTN | 0.115 |
| AP2A1 | SEC16A | SEC16A->AP2A1 | 0.076 | AP2A1->SEC16A | 0.088 | — | — |
| COPG1 | SEC31A | COPG1->SEC31A | 0.079 | SEC31A->COPG1 | 0.086 | — | — |
| ARF6 | COPB1 | COPB1->ARF6 | 0.077 | ARF6->COPB1 | 0.084 | — | — |
| APP | SORT1 | APP->SORT1 | 0.097 | SORT1->APP | 0.054 | — | — |
| MFN2 | PSEN1 | PSEN1->MFN2 | 0.074 | PSEN1->MFN2 | 0.077 | — | — |
| ARF6 | RAB1B | ARF6->RAB1B | 0.106 | RAB1B->ARF6 | 0.039 | — | — |
| COPZ1 | RTN3 | — | — | RTN3->COPZ1 | 0.041 | COPZ1->RTN3 | 0.102 |
| APP | SNX6 | — | — | SNX6->APP | 0.043 | SNX6->APP | 0.099 |
| RTN4 | VDAC1 | VDAC1->RTN4 | 0.058 | RTN4->VDAC1 | 0.083 | — | — |
| SORL1 | VPS29 | SORL1->VPS29 | 0.088 | SORL1->VPS29 | 0.052 | — | — |
| COPG2 | TMED10 | TMED10->COPG2 | 0.072 | TMED10->COPG2 | 0.068 | — | — |
| COPB1 | SEC24C | COPB1->SEC24C | 0.068 | SEC24C->COPB1 | 0.069 | — | — |
| APP | TMED2 | APP->TMED2 | 0.056 | — | — | APP->TMED2 | 0.079 |
| COPB2 | SEC31A | COPB2->SEC31A | 0.052 | SEC31A->COPB2 | 0.078 | — | — |
| COPB2 | SEC13 | COPB2->SEC13 | 0.056 | SEC13->COPB2 | 0.072 | — | — |
| RAB1B | SAR1B | SAR1B->RAB1B | 0.077 | SAR1B->RAB1B | 0.050 | — | — |
| ABCA7 | SORL1 | ABCA7->SORL1 | 0.084 | ABCA7->SORL1 | 0.036 | — | — |
| PICALM | SEC13 | SEC13->PICALM | 0.037 | PICALM->SEC13 | 0.082 | — | — |
| BACE1 | VPS26A | — | — | BACE1->VPS26A | 0.073 | BACE1->VPS26A | 0.040 |
| MIA3 | SEC31A | MIA3->SEC31A | 0.032 | SEC31A->MIA3 | 0.081 | — | — |
| RAB1B | RTN4 | RAB1B->RTN4 | 0.066 | RTN4->RAB1B | 0.046 | — | — |
| RAB1B | VDAC1 | VDAC1->RAB1B | 0.068 | RAB1B->VDAC1 | 0.044 | — | — |
| COPG1 | SEC24C | SEC24C->COPG1 | 0.056 | SEC24C->COPG1 | 0.054 | — | — |
| MIA3 | SEC13 | MIA3->SEC13 | 0.047 | MIA3->SEC13 | 0.062 | — | — |
| COPE | RTN3 | COPE->RTN3 | 0.073 | COPE->RTN3 | 0.035 | — | — |
| AP1G1 | COPB1 | COPB1->AP1G1 | 0.057 | AP1G1->COPB1 | 0.051 | — | — |
| SEC24B | SREBF2 | SEC24B->SREBF2 | 0.036 | SREBF2->SEC24B | 0.071 | — | — |
| EIF2S1 | HSPA9 | HSPA9->EIF2S1 | 0.046 | EIF2S1->HSPA9 | 0.060 | — | — |
| COPA | PREB | PREB->COPA | 0.029 | PREB->COPA | 0.076 | — | — |
| COPB2 | SAR1B | — | — | SAR1B->COPB2 | 0.058 | SAR1B->COPB2 | 0.044 |
| GGA3 | SEC24C | GGA3->SEC24C | 0.031 | GGA3->SEC24C | 0.072 | — | — |
| OSBP | SEC24D | — | — | SEC24D->OSBP | 0.057 | SEC24D->OSBP | 0.045 |
| COPB2 | SEC24C | SEC24C->COPB2 | 0.054 | COPB2->SEC24C | 0.048 | — | — |
| COPG2 | SAR1A | SAR1A->COPG2 | 0.049 | SAR1A->COPG2 | 0.052 | — | — |
| OSBP | SEC13 | SEC13->OSBP | 0.060 | OSBP->SEC13 | 0.040 | — | — |
| RAB1A | SAR1A | RAB1A->SAR1A | 0.030 | RAB1A->SAR1A | 0.069 | — | — |
| APP | MFN2 | MFN2->APP | 0.030 | — | — | MFN2->APP | 0.070 |
| HSPA5 | MFN2 | HSPA5->MFN2 | 0.054 | MFN2->HSPA5 | 0.044 | — | — |
| COPE | SEC31A | COPE->SEC31A | 0.032 | SEC31A->COPE | 0.067 | — | — |
| COPB2 | PREB | COPB2->PREB | 0.046 | PREB->COPB2 | 0.051 | — | — |
| PREB | TMED9 | PREB->TMED9 | 0.041 | PREB->TMED9 | 0.055 | — | — |
| COPG1 | COPG2 | COPG1->COPG2 | 0.047 | COPG1->COPG2 | 0.049 | — | — |
| COPZ1 | SEC13 | SEC13->COPZ1 | 0.035 | SEC13->COPZ1 | 0.059 | — | — |
| HSPA5 | ITPR1 | HSPA5->ITPR1 | 0.040 | HSPA5->ITPR1 | 0.053 | — | — |
| COPB2 | SEC23A | SEC23A->COPB2 | 0.047 | SEC23A->COPB2 | 0.046 | — | — |
| ARF6 | RAB1A | ARF6->RAB1A | 0.028 | RAB1A->ARF6 | 0.065 | — | — |
| COPE | SEC13 | SEC13->COPE | 0.032 | SEC13->COPE | 0.061 | — | — |
| HSPA9 | ITPR2 | HSPA9->ITPR2 | 0.056 | HSPA9->ITPR2 | 0.034 | — | — |
| COPB1 | SEC24A | SEC24A->COPB1 | 0.036 | COPB1->SEC24A | 0.049 | — | — |
| HSPA5 | SEC13 | SEC13->HSPA5 | 0.047 | SEC13->HSPA5 | 0.038 | — | — |
| APP | VPS29 | — | — | VPS29->APP | 0.028 | APP->VPS29 | 0.056 |
| AP2B1 | COPG1 | AP2B1->COPG1 | 0.028 | AP2B1->COPG1 | 0.055 | — | — |
| CLU | VDAC1 | VDAC1->CLU | 0.026 | VDAC1->CLU | 0.057 | — | — |
| VDAC1 | VPS35 | VPS35->VDAC1 | 0.052 | VPS35->VDAC1 | 0.030 | — | — |
| HSPA5 | VDAC1 | HSPA5->VDAC1 | 0.051 | HSPA5->VDAC1 | 0.030 | — | — |
| HSPA5 | SEC31A | HSPA5->SEC31A | 0.043 | SEC31A->HSPA5 | 0.035 | — | — |
| ARF6 | SYK | ARF6->SYK | 0.047 | SYK->ARF6 | 0.029 | — | — |
| HSPA5 | SAR1A | HSPA5->SAR1A | 0.032 | HSPA5->SAR1A | 0.039 | — | — |
| AP2A1 | COPB1 | COPB1->AP2A1 | 0.035 | AP2A1->COPB1 | 0.033 | — | — |
| LRP1 | RTN4 | — | — | RTN4->LRP1 | 0.029 | RTN4->LRP1 | 0.039 |
| COPB1 | SAR1A | COPB1->SAR1A | 0.030 | COPB1->SAR1A | 0.036 | — | — |
| BIN1 | SORL1 | — | — | SORL1->BIN1 | 0.018 | SORL1->BIN1 | 0.049 |
| COPB1 | GGA2 | COPB1->GGA2 | 0.034 | GGA2->COPB1 | 0.033 | — | — |
| LRP1 | SORT1 | LRP1->SORT1 | 0.024 | SORT1->LRP1 | 0.041 | — | — |
| HSPA5 | RAB1A | HSPA5->RAB1A | 0.031 | RAB1A->HSPA5 | 0.031 | — | — |
| COPB2 | SEC24A | COPB2->SEC24A | 0.029 | COPB2->SEC24A | 0.034 | — | — |
| AP2B1 | COPG2 | AP2B1->COPG2 | 0.038 | COPG2->AP2B1 | 0.023 | — | — |
| COPB2 | GGA1 | COPB2->GGA1 | 0.032 | COPB2->GGA1 | 0.028 | — | — |
| BIN1 | PSEN1 | BIN1->PSEN1 | 0.033 | BIN1->PSEN1 | 0.026 | — | — |
| COPG1 | HSPA5 | HSPA5->COPG1 | 0.024 | COPG1->HSPA5 | 0.034 | — | — |
| MIA3 | SAR1B | — | — | MIA3->SAR1B | 0.019 | SAR1B->MIA3 | 0.038 |
| COPB2 | SEC23B | COPB2->SEC23B | 0.026 | COPB2->SEC23B | 0.022 | — | — |
| BIN1 | VPS35 | BIN1->VPS35 | 0.025 | VPS35->BIN1 | 0.021 | — | — |

---

## 7b. This run vs the run in `results/` (previous prior)

| dataset | AUROC_raw_prev | AUROC_raw_this | AUROC_adj_prev | AUROC_adj_this | prior_edges_prev | prior_edges_this | supported_prev | supported_this | STRONG_orig_prev | STRONG_orig_this |
|---|---|---|---|---|---|---|---|---|---|---|
| ROSMAP | 0.5618 | 0.5926 | 0.5946 | 0.6232 | 688 | 1308 | 172 | 251 | 0 | 0 |
| Diverse | 0.5491 | 0.5767 | 0.6431 | 0.6441 | 688 | 1308 | 354 | 464 | 28 | 34 |
| Banner | 0.5313 | 0.5535 | 0.5811 | 0.5825 | 688 | 1308 | 101 | 143 | 0 | 0 |
| BannerLFQ | 0.5372 | 0.5283 | 0.5896 | 0.5369 | 688 | 1308 | 135 | 244 | 0 | 0 |

| dataset | Z* previous | Z* this run |
|---|---|---|
| ROSMAP | proteomePC2, proteomePC3, proteomePC5, proteomePC10 | proteomePC6, proteomePC2, celltype_neuron, celltype_oligodendrocyte, proteomePC8 |
| Diverse | proteomePC1, proteomePC8, celltype_astrocyte, amyAny, proteomePC9 | celltype_oligodendrocyte, reag, proteomePC2, proteomePC8, proteomePC6 |
| Banner | proteomePC3, proteomePC1, proteomePC5, proteomePC8, celltype_oligodendrocyte | celltype_oligodendrocyte, proteomePC3, celltype_neuron, proteomePC8, apoeGenotype |
| BannerLFQ | proteomePC1, Braak, TangleTotal, PlaqueTotal, proteomePC9 | celltype_neuron, proteomePC9 |

Supported protein pairs (unordered) in the two runs:

| dataset | supported_pairs_prev | supported_pairs_this | supported_in_both | prev_supported_not_prior_now | this_supported_not_prior_before | jaccard |
|---|---|---|---|---|---|---|
| ROSMAP | 172 | 251 | 86 | 52 | 135 | 0.2550 |
| Diverse | 354 | 464 | 168 | 140 | 255 | 0.2580 |
| Banner | 101 | 143 | 53 | 33 | 68 | 0.2770 |
| BannerLFQ | 135 | 244 | 53 | 57 | 139 | 0.1630 |

`prev_supported_not_prior_now` = pairs supported before whose link is no longer a prior edge in this matrix; `this_supported_not_prior_before` = supported pairs that were not prior edges before.

Provenance of each run's supported edges (the arrow chosen, looked up in `Edge_List`):

| dataset | supported | AUTO | LIT | COMPLEX | CHRONOS |
|---|---|---|---|---|---|
| ROSMAP | 251 | 134 | 74 | 28 | 15 |
| Diverse | 464 | 258 | 153 | 32 | 21 |
| Banner | 143 | 63 | 47 | 19 | 14 |
| BannerLFQ | 244 | 136 | 86 | 14 | 8 |

The matrix has 50 off-diagonal edges tagged `CHRONOS`, i.e. set using the earlier 3-cohort run on these same datasets.

---

## 8. Banner-specific checks

### 8.1 TMT batch mapping (HIW 3.1)

Mean of Y-chromosome proteins (RPS4Y1, DDX3Y, EIF1AY, KDM5D, USP9Y, NLGN4Y) per sample, 198 samples with recorded sex:

| | n | mean | min | max |
|---|---|---|---|---|
| female | 86 | −0.541 | −1.817 | −0.072 |
| male | 112 | +0.243 | −0.220 | +0.483 |

- best single cut-off (0.022): 99.0% correctly classified (3 females above the lowest male, 2 males below the highest female)
- with batch labels randomly shuffled (200 times): best accuracy mean 57.8%, maximum 64.1%
- BannerLFQ detects only one Y protein (USP9Y, in 20 samples), so this check is not possible there

### 8.2 Banner TMT vs Banner LFQ (same 190 people)

| | TMT | LFQ |
|---|---|---|
| graph nodes present | 80 | 58 |
| signed AUROC raw → adjusted | 0.5535 → 0.5825 | 0.5283 → 0.5369 |
| Z* | celltype_oligodendrocyte, proteomePC3, celltype_neuron, proteomePC8, apoeGenotype | celltype_neuron, proteomePC9 |
| supported edges | 143 | 244 |
| STRONG regions (original) | 0 | 0 |

From `comparison_report.txt`:

```
6. TECHNICAL CHECK: Banner TMT vs Banner LFQ (same people, different mass-spec method)
   edges testable in both: 386 pairs | supported TMT 124, LFQ 244, both 80 (chance 78.4, p=4.0e-01)
   Spearman of |W| across shared testable edges: -0.008
```

### 8.3 Possible overlap with Diverse

- Diverse contains 43 donors from the Banner cohort; 17 of them are among the 964 Diverse people analysed.
- No ID links the two systems. Matching on sex + APOE + age at death + Braak gave 8 Banner donors
  (09-17, 11-14, 01-16, 06-09, 05-57, 05-35, 13-46, 13-66), matched to 7 Diverse donors; all 8 were removed from Banner and BannerLFQ.

---

## 9. File index (`results_v3/`)

| file | content |
|---|---|
| `{ROSMAP,Diverse,Banner,BannerLFQ}_causal_W.csv` | 80 × 80 weights |
| `…_causal_beta.csv` | 80 × 80 effect sizes |
| `…_causal_N.csv` | 80 × 80 people per pair |
| `…_causal_edges.csv` | 1308 prior edges with W, β, n, p, q, data_supported |
| `…_variables.csv` | every variable at every selection step |
| `…_regions.csv` | every region (original scoring) |
| `{ROSMAP,Diverse,Banner}_regions_usable_null_check.csv` | regions re-scored with the corrected null |
| `Diverse_best_region_W.csv` | graph refitted in the original top Diverse region |
| `…_summary.json`, `…_report.txt`, `….log`, `…_individuals.txt` | run summaries, logs, people used |
| `comparison_report.txt`, `comparison_pairs.csv`, `consensus_edges.csv` | cross-dataset comparison |
