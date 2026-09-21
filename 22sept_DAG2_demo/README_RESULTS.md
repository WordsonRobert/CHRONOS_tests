# CHRONOS causal graphs — all results

Everything the pipeline produced, dataset by dataset, with as little interpretation as possible.
How each number is computed is in `README_HOW_IT_WORKS.md`; section numbers like (HIW 4.7) point there.
All tables below were generated from the files in `results/` by `make_results_readme.py`.

**One correction to results reported earlier.** The first region search reported 28 STRONG
regions in Diverse. Re-scoring with a null drawn only from people who have adjusted data
(HIW 9.1) gives **0 STRONG regions in every dataset**. Both versions are listed in section 5.
Nothing else changes.

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
| people with usable data after adjustment for Z* | 400 | **676** (amyAny missing for 288) | 190 | 190 |
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
| directed prior edges (C ≥ 0.05) | 688 |
| unordered pairs with a prior edge | 507 of 3,160 |
| one-way pairs | 326 |
| two-way pairs | 181 (78 with equal C both ways) |
| C values among edges | 0.2: 148, 0.3: 235, 0.4: 129, 0.5: 88, 0.6: 56, 0.7: 23, 0.8: 7, 0.9: 2 |
| largest out-degree | TMED10 61, TMED2 57, TMED9 45, APOE 28, SEC24D 15 |
| largest in-degree | APP 72, BACE1 52, TMED10 36, TMED2 29, LRP1 26 |

---
## 3. Signal: how well each dataset recovers C (AUROC, HIW 4.4)

0.5 = no relation to C. "Signed" (used by the pipeline) ranks pairs by ρ; "unsigned" by |ρ|.

| | ROSMAP | Diverse | Banner | BannerLFQ |
|---|---|---|---|---|
| signed AUROC, raw | 0.5618 | 0.5491 | 0.5313 | 0.5372 |
| signed AUROC, adjusted for Z* | **0.5946** | **0.6431** | **0.5811** | **0.5896** |
| unsigned AUROC, raw | 0.5151 | 0.5195 | 0.4879 | 0.5209 |
| unsigned AUROC, adjusted | 0.5246 | 0.5933 | 0.5353 | 0.5591 |
| pairs scored (n ≥ 20) | 3,160 | 3,160 | 3,157 | 1,533 |

Correlation between prior-linked pairs vs other pairs:

| | ROSMAP raw | ROSMAP adj | Diverse raw | Diverse adj | Banner raw | Banner adj | LFQ raw | LFQ adj |
|---|---|---|---|---|---|---|---|---|
| prior pairs: mean ρ | +0.057 | +0.056 | +0.155 | +0.132 | +0.038 | +0.048 | +0.131 | +0.086 |
| prior pairs: fraction ρ > 0 | 0.623 | 0.643 | 0.767 | 0.803 | 0.558 | 0.625 | 0.699 | 0.657 |
| prior pairs: mean abs(ρ) | 0.127 | 0.109 | 0.238 | 0.171 | 0.136 | 0.121 | 0.270 | 0.197 |
| other pairs: mean ρ | +0.023 | +0.013 | +0.112 | +0.048 | +0.017 | +0.006 | +0.092 | +0.020 |
| other pairs: fraction ρ > 0 | 0.533 | 0.521 | 0.686 | 0.639 | 0.526 | 0.503 | 0.662 | 0.550 |
| other pairs: mean abs(ρ) | 0.117 | 0.098 | 0.232 | 0.127 | 0.135 | 0.106 | 0.258 | 0.171 |

---
## 4. Trajectory variables (HIW 4.6)

### 4.1 Selected Z*, in order

**ROSMAP** — start AUROC 0.5618, end 0.5946

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | proteomePC2 | proteome-derived | 0.5824 | 0.5617 | 0.0206 | 24.46 | 400 |
| 2 | proteomePC3 | proteome-derived | 0.5892 | 0.5823 | 0.0067 | 26.05 | 400 |
| 3 | proteomePC5 | proteome-derived | 0.5927 | 0.5890 | 0.0035 | 11.69 | 400 |
| 4 | proteomePC10 | proteome-derived | 0.5946 | 0.5924 | 0.0019 | 2.54 | 400 |

**Diverse** — start AUROC 0.5491, end 0.6431

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | proteomePC1 | proteome-derived | 0.6186 | 0.5491 | 0.0695 | 775.31 | 964 |
| 2 | proteomePC8 | proteome-derived | 0.6244 | 0.6186 | 0.0058 | 30.32 | 964 |
| 3 | celltype_astrocyte | proteome-derived | 0.6294 | 0.6244 | 0.0050 | 23.92 | 964 |
| 4 | amyAny | clinical/pathology | 0.6355 | 0.6340 | 0.0060 | 7.01 | 676 |
| 5 | proteomePC9 | proteome-derived | 0.6431 | 0.6356 | 0.0076 | 25.26 | 676 |

**Banner** — start AUROC 0.5313, end 0.5811

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | proteomePC3 | proteome-derived | 0.5555 | 0.5321 | 0.0242 | 32.56 | 190 |
| 2 | proteomePC1 | proteome-derived | 0.5674 | 0.5553 | 0.0118 | 13.00 | 190 |
| 3 | proteomePC5 | proteome-derived | 0.5724 | 0.5665 | 0.0050 | 7.10 | 190 |
| 4 | proteomePC8 | proteome-derived | 0.5775 | 0.5719 | 0.0051 | 6.24 | 190 |
| 5 | celltype_oligodendrocyte | proteome-derived | 0.5811 | 0.5772 | 0.0036 | 3.88 | 190 |

**BannerLFQ** — start AUROC 0.5372, end 0.5896

| step | variable | kind | AUROC_after | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|---|
| 1 | proteomePC1 | proteome-derived | 0.5625 | 0.5374 | 0.0253 | 35.54 | 190 |
| 2 | Braak | clinical/pathology | 0.5831 | 0.5647 | 0.0207 | 12.90 | 190 |
| 3 | TangleTotal | clinical/pathology | 0.5864 | 0.5834 | 0.0033 | 2.28 | 190 |
| 4 | PlaqueTotal | clinical/pathology | 0.5883 | 0.5849 | 0.0019 | 2.62 | 190 |
| 5 | proteomePC9 | proteome-derived | 0.5896 | 0.5863 | 0.0013 | 3.02 | 190 |

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

**ROSMAP** (unadjusted AUROC 0.5618)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| proteomePC2 | proteome-derived | 0.5824 | 0.5617 | 0.0206 | 24.46 | 400 |
| proteomePC5 | proteome-derived | 0.5667 | 0.5614 | 0.0049 | 11.63 | 400 |
| proteomePC3 | proteome-derived | 0.5652 | 0.5619 | 0.0034 | 7.75 | 400 |
| celltype_oligodendrocyte | proteome-derived | 0.5648 | 0.5618 | 0.0030 | 5.67 | 400 |
| proteomePC10 | proteome-derived | 0.5625 | 0.5615 | 0.0007 | 2.40 | 400 |
| proteomePC7 | proteome-derived | 0.5623 | 0.5615 | 0.0005 | 1.05 | 400 |
| age_at_visit_max | clinical/pathology | 0.5622 | 0.5613 | 0.0004 | 1.44 | 400 |
| apoe_genotype | clinical/pathology | 0.5621 | 0.5617 | 0.0003 | 0.63 | 400 |
| age_death | clinical/pathology | 0.5621 | 0.5615 | 0.0003 | 1.20 | 400 |
| msex | clinical/pathology | 0.5621 | 0.5617 | 0.0003 | 0.87 | 400 |
| educ | clinical/pathology | 0.5617 | 0.5613 | -0.0001 | 0.53 | 400 |
| braaksc | clinical/pathology | 0.5616 | 0.5615 | -0.0002 | 0.29 | 400 |
| Study | technical | 0.5614 | 0.5618 | -0.0004 | -1.02 | 400 |
| celltype_microglia | proteome-derived | 0.5610 | 0.5612 | -0.0008 | -0.36 | 400 |
| celltype_neuron | proteome-derived | 0.5609 | 0.5616 | -0.0009 | -1.14 | 400 |
| proteomePC8 | proteome-derived | 0.5607 | 0.5617 | -0.0011 | -2.01 | 400 |
| cts_mmse30_first_ad_dx | clinical/pathology | 0.5605 | 0.5609 | -0.0013 | -0.37 | 126 |
| pmi | technical | 0.5604 | 0.5606 | -0.0014 | -0.47 | 398 |
| ceradsc | clinical/pathology | 0.5603 | 0.5615 | -0.0015 | -1.44 | 400 |
| dcfdx_lv | clinical/pathology | 0.5602 | 0.5617 | -0.0016 | -2.55 | 400 |
| cogdx | clinical/pathology | 0.5601 | 0.5617 | -0.0017 | -2.21 | 400 |
| celltype_astrocyte | proteome-derived | 0.5599 | 0.5616 | -0.0019 | -3.83 | 400 |
| cts_mmse30_lv | clinical/pathology | 0.5596 | 0.5614 | -0.0022 | -3.28 | 400 |
| proteomePC1 | proteome-derived | 0.5593 | 0.5618 | -0.0025 | -4.05 | 400 |
| celltype_endothelial | proteome-derived | 0.5591 | 0.5615 | -0.0027 | -3.30 | 400 |
| age_first_ad_dx | clinical/pathology | 0.5579 | 0.5608 | -0.0039 | -1.50 | 132 |
| proteomePC6 | proteome-derived | 0.5577 | 0.5615 | -0.0041 | -16.54 | 400 |
| proteomePC9 | proteome-derived | 0.5561 | 0.5615 | -0.0057 | -14.74 | 400 |
| proteomePC4 | proteome-derived | 0.5449 | 0.5616 | -0.0169 | -19.45 | 400 |

**Diverse** (unadjusted AUROC 0.5491)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| proteomePC1 | proteome-derived | 0.6186 | 0.5491 | 0.0695 | 775.31 | 964 |
| celltype_neuron | proteome-derived | 0.5894 | 0.5491 | 0.0402 | 249.43 | 964 |
| amyAny | clinical/pathology | 0.5672 | 0.5592 | 0.0181 | 28.43 | 676 |
| amyCerad | clinical/pathology | 0.5664 | 0.5593 | 0.0173 | 17.01 | 676 |
| reag | clinical/pathology | 0.5657 | 0.5592 | 0.0165 | 16.29 | 676 |
| celltype_oligodendrocyte | proteome-derived | 0.5602 | 0.5491 | 0.0111 | 60.20 | 964 |
| ADoutcome | clinical/pathology | 0.5555 | 0.5485 | 0.0064 | 37.32 | 961 |
| Braak | clinical/pathology | 0.5534 | 0.5494 | 0.0043 | 15.24 | 908 |
| bScore | clinical/pathology | 0.5530 | 0.5492 | 0.0039 | 13.59 | 908 |
| proteomePC6 | proteome-derived | 0.5520 | 0.5491 | 0.0029 | 20.96 | 964 |
| PMI | technical | 0.5514 | 0.5514 | 0.0022 | 0.02 | 818 |
| proteomePC8 | proteome-derived | 0.5509 | 0.5491 | 0.0018 | 19.28 | 964 |
| proteomePC4 | proteome-derived | 0.5506 | 0.5492 | 0.0015 | 9.50 | 964 |
| proteomePC3 | proteome-derived | 0.5502 | 0.5492 | 0.0011 | 4.64 | 964 |
| isHispanic | clinical/pathology | 0.5496 | 0.5489 | 0.0005 | 5.19 | 963 |
| ageDeath | clinical/pathology | 0.5493 | 0.5490 | 0.0002 | 1.55 | 963 |
| amyA | clinical/pathology | 0.5492 | 0.5488 | 0.0001 | 0.81 | 398 |
| sex | clinical/pathology | 0.5492 | 0.5492 | 0.0000 | 0.10 | 964 |
| race | clinical/pathology | 0.5490 | 0.5486 | -0.0001 | 1.00 | 961 |
| amyThal | clinical/pathology | 0.5487 | 0.5483 | -0.0004 | 0.95 | 391 |
| proteomePC10 | proteome-derived | 0.5486 | 0.5491 | -0.0005 | -3.83 | 964 |
| apoeGenotype | clinical/pathology | 0.5483 | 0.5479 | -0.0009 | 2.23 | 895 |
| proteomePC9 | proteome-derived | 0.5480 | 0.5492 | -0.0011 | -8.74 | 964 |
| dataContributionGroup | technical | 0.5473 | 0.5490 | -0.0018 | -5.55 | 964 |
| derivedOutcomeBasedOnMayoDx | clinical/pathology | 0.5471 | 0.5491 | -0.0020 | -16.34 | 964 |
| proteomePC5 | proteome-derived | 0.5470 | 0.5490 | -0.0021 | -13.56 | 964 |
| cohort | technical | 0.5467 | 0.5491 | -0.0025 | -5.68 | 964 |
| celltype_microglia | proteome-derived | 0.5454 | 0.5491 | -0.0038 | -32.34 | 964 |
| proteomePC7 | proteome-derived | 0.5434 | 0.5492 | -0.0058 | -33.44 | 964 |
| celltype_astrocyte | proteome-derived | 0.5362 | 0.5491 | -0.0130 | -133.50 | 964 |
| mayoDx | clinical/pathology | 0.5283 | 0.5260 | -0.0208 | 5.98 | 285 |
| celltype_endothelial | proteome-derived | 0.5276 | 0.5492 | -0.0215 | -91.26 | 964 |
| proteomePC2 | proteome-derived | 0.5071 | 0.5491 | -0.0420 | -200.82 | 964 |

**Banner** (unadjusted AUROC 0.5313)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| proteomePC3 | proteome-derived | 0.5555 | 0.5321 | 0.0242 | 32.56 | 190 |
| proteomePC8 | proteome-derived | 0.5354 | 0.5317 | 0.0041 | 5.83 | 190 |
| proteomePC5 | proteome-derived | 0.5344 | 0.5315 | 0.0030 | 5.72 | 190 |
| TangleTotal | clinical/pathology | 0.5343 | 0.5312 | 0.0030 | 2.16 | 190 |
| celltype_neuron | proteome-derived | 0.5342 | 0.5311 | 0.0029 | 3.41 | 190 |
| Braak | clinical/pathology | 0.5337 | 0.5317 | 0.0024 | 3.52 | 190 |
| celltype_astrocyte | proteome-derived | 0.5331 | 0.5316 | 0.0018 | 1.69 | 190 |
| celltype_endothelial | proteome-derived | 0.5325 | 0.5312 | 0.0012 | 1.17 | 190 |
| apoeGenotype | clinical/pathology | 0.5317 | 0.5304 | 0.0004 | 0.79 | 189 |
| sex | clinical/pathology | 0.5307 | 0.5309 | -0.0006 | -0.18 | 190 |
| diagnosis | clinical/pathology | 0.5307 | 0.5308 | -0.0006 | -0.06 | 190 |
| ageDeath | clinical/pathology | 0.5307 | 0.5310 | -0.0007 | -0.40 | 190 |
| proteomePC7 | proteome-derived | 0.5305 | 0.5307 | -0.0008 | -0.19 | 190 |
| proteomePC6 | proteome-derived | 0.5305 | 0.5314 | -0.0008 | -1.54 | 190 |
| proteomePC10 | proteome-derived | 0.5303 | 0.5317 | -0.0010 | -1.32 | 190 |
| lastMMSE | clinical/pathology | 0.5301 | 0.5319 | -0.0013 | -1.65 | 190 |
| proteomePC4 | proteome-derived | 0.5294 | 0.5314 | -0.0019 | -2.82 | 190 |
| pmi | technical | 0.5294 | 0.5314 | -0.0019 | -2.10 | 190 |
| celltype_microglia | proteome-derived | 0.5289 | 0.5306 | -0.0024 | -1.34 | 190 |
| proteomePC1 | proteome-derived | 0.5289 | 0.5312 | -0.0025 | -2.28 | 190 |
| PlaqueTotal | clinical/pathology | 0.5287 | 0.5311 | -0.0026 | -2.86 | 190 |
| celltype_oligodendrocyte | proteome-derived | 0.5251 | 0.5312 | -0.0062 | -6.91 | 190 |
| CERAD | clinical/pathology | 0.5249 | 0.5306 | -0.0065 | -7.48 | 190 |
| proteomePC9 | proteome-derived | 0.5224 | 0.5316 | -0.0089 | -12.71 | 190 |
| proteomePC2 | proteome-derived | 0.5179 | 0.5313 | -0.0134 | -12.96 | 190 |

**BannerLFQ** (unadjusted AUROC 0.5372)

| variable | kind | auroc | shuffled_mean | gain | z | people |
|---|---|---|---|---|---|---|
| proteomePC1 | proteome-derived | 0.5625 | 0.5374 | 0.0253 | 35.54 | 190 |
| Braak | clinical/pathology | 0.5464 | 0.5376 | 0.0092 | 7.73 | 190 |
| diagnosis | clinical/pathology | 0.5420 | 0.5381 | 0.0048 | 2.91 | 190 |
| TangleTotal | clinical/pathology | 0.5420 | 0.5382 | 0.0048 | 3.30 | 190 |
| proteomePC5 | proteome-derived | 0.5403 | 0.5383 | 0.0032 | 1.74 | 190 |
| proteomePC6 | proteome-derived | 0.5403 | 0.5376 | 0.0031 | 1.60 | 190 |
| proteomePC9 | proteome-derived | 0.5396 | 0.5377 | 0.0024 | 1.67 | 190 |
| CERAD | clinical/pathology | 0.5392 | 0.5373 | 0.0020 | 1.43 | 190 |
| proteomePC3 | proteome-derived | 0.5387 | 0.5378 | 0.0016 | 0.66 | 190 |
| PlaqueTotal | clinical/pathology | 0.5385 | 0.5385 | 0.0014 | -0.00 | 190 |
| lastMMSE | clinical/pathology | 0.5376 | 0.5381 | 0.0005 | -0.32 | 190 |
| proteomePC4 | proteome-derived | 0.5374 | 0.5395 | 0.0003 | -1.21 | 190 |
| proteomePC8 | proteome-derived | 0.5374 | 0.5381 | 0.0002 | -0.42 | 190 |
| apoeGenotype | clinical/pathology | 0.5372 | 0.5373 | 0.0001 | -0.04 | 189 |
| pmi | technical | 0.5363 | 0.5376 | -0.0009 | -1.17 | 190 |
| ageDeath | clinical/pathology | 0.5363 | 0.5380 | -0.0009 | -1.05 | 190 |
| sex | clinical/pathology | 0.5360 | 0.5382 | -0.0011 | -2.40 | 190 |
| celltype_oligodendrocyte | proteome-derived | 0.5359 | 0.5383 | -0.0012 | -1.45 | 190 |
| proteomePC10 | proteome-derived | 0.5359 | 0.5384 | -0.0013 | -1.44 | 190 |
| celltype_astrocyte | proteome-derived | 0.5357 | 0.5377 | -0.0015 | -1.96 | 190 |
| proteomePC7 | proteome-derived | 0.5355 | 0.5377 | -0.0016 | -1.69 | 190 |
| celltype_endothelial | proteome-derived | 0.5303 | 0.5288 | -0.0068 | 1.01 | 101 |
| celltype_neuron | proteome-derived | 0.5254 | 0.5388 | -0.0118 | -7.47 | 190 |
| proteomePC2 | proteome-derived | 0.5043 | 0.5389 | -0.0328 | -18.77 | 190 |
| celltype_microglia | proteome-derived | — | — | — | — | 32 |

### 4.3 Clinical/pathology variables at every step

z vs shuffled for each non-proteome variable at each selection step (the full per-step lists are in `{ds}_variables.csv`). Bold = selected at that step.

**ROSMAP**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| Study | -1.02 | -1.15 | -0.60 | -0.94 | -0.17 |
| age_at_visit_max | 1.44 | -0.82 | 0.30 | 1.17 | -0.42 |
| age_death | 1.20 | -0.02 | 0.01 | 1.33 | -1.27 |
| age_first_ad_dx | -1.50 | -2.55 | -2.89 | -1.88 | -0.64 |
| apoe_genotype | 0.63 | 0.24 | 1.04 | 0.29 | -0.36 |
| braaksc | 0.29 | -0.15 | 2.22 | 2.00 | 1.60 |
| ceradsc | -1.44 | -1.12 | 1.86 | 1.69 | 0.29 |
| cogdx | -2.21 | -3.81 | -2.10 | -4.77 | -3.20 |
| cts_mmse30_first_ad_dx | -0.37 | 0.77 | -0.34 | -1.71 | -0.45 |
| cts_mmse30_lv | -3.28 | -3.35 | -3.05 | -1.63 | -1.30 |
| dcfdx_lv | -2.55 | -5.60 | -3.26 | -2.33 | -4.63 |
| educ | 0.53 | -0.94 | 0.08 | -0.17 | -0.08 |
| msex | 0.87 | -0.45 | -0.79 | -2.04 | -0.51 |
| pmi | -0.47 | -1.90 | -2.28 | -4.11 | -1.25 |

**Diverse**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| ADoutcome | 37.32 | 4.07 | 8.24 | 2.14 | -1.00 |
| Braak | 15.24 | -3.74 | -0.28 | -1.72 | -2.59 |
| PMI | 0.02 | 1.65 | 0.65 | -0.22 | 1.56 |
| ageDeath | 1.55 | 4.03 | 3.54 | 7.24 | -2.06 |
| amyA | 0.81 | -5.64 | -5.36 | -2.83 | -2.71 |
| amyAny | 28.43 | 9.02 | 6.66 | **7.01** | — |
| amyCerad | 17.01 | 0.34 | -0.19 | -0.62 | -2.91 |
| amyThal | 0.95 | -5.57 | -8.34 | -6.27 | -1.95 |
| apoeGenotype | 2.23 | -0.33 | -2.19 | -4.75 | -2.16 |
| bScore | 13.59 | -1.24 | 0.03 | -0.68 | -2.80 |
| cohort | -5.68 | 2.38 | 1.87 | -3.07 | -0.03 |
| dataContributionGroup | -5.55 | -5.64 | -3.92 | -2.16 | -2.28 |
| derivedOutcomeBasedOnMayoDx | -16.34 | -2.67 | -6.57 | -4.52 | 0.02 |
| isHispanic | 5.19 | -1.03 | -1.16 | 0.74 | -0.32 |
| mayoDx | 5.98 | -1.72 | -7.65 | -2.28 | — |
| race | 1.00 | -1.76 | -7.38 | -4.69 | -2.50 |
| reag | 16.29 | -1.87 | -7.16 | -3.25 | -4.05 |
| sex | 0.10 | -2.59 | -4.17 | -7.66 | -3.73 |

**Banner**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| Braak | 3.52 | 11.55 | 4.20 | 1.72 | 2.62 |
| CERAD | -7.48 | -3.22 | -8.59 | -6.12 | -6.15 |
| PlaqueTotal | -2.86 | 3.77 | -2.59 | -4.53 | -2.36 |
| TangleTotal | 2.16 | 11.89 | 3.22 | 0.39 | 0.88 |
| ageDeath | -0.40 | -0.37 | -0.66 | -0.16 | -0.37 |
| apoeGenotype | 0.79 | 4.43 | 3.43 | 1.46 | 1.83 |
| diagnosis | -0.06 | 16.96 | -1.08 | -2.06 | -0.67 |
| lastMMSE | -1.65 | 9.38 | 0.43 | -3.56 | 0.65 |
| pmi | -2.10 | -0.61 | -1.57 | -0.31 | -1.87 |
| sex | -0.18 | -0.64 | -1.48 | -1.93 | -2.23 |

**BannerLFQ**

| variable | step 1 | step 2 | step 3 | step 4 | step 5 |
|---|---|---|---|---|---|
| Braak | 7.73 | **12.90** | — | — | — |
| CERAD | 1.43 | 5.17 | 1.12 | 1.65 | 1.24 |
| PlaqueTotal | -0.00 | 1.31 | 1.25 | **2.62** | — |
| TangleTotal | 3.30 | 5.24 | **2.28** | — | — |
| ageDeath | -1.05 | -3.70 | -3.89 | -4.13 | -3.59 |
| apoeGenotype | -0.04 | 0.51 | 1.49 | 0.37 | 1.29 |
| diagnosis | 2.91 | 3.30 | 0.24 | 0.83 | -0.19 |
| lastMMSE | -0.32 | 1.72 | -2.02 | -4.46 | -2.99 |
| pmi | -1.17 | -0.68 | 0.07 | 0.50 | 0.73 |
| sex | -2.40 | -1.82 | -4.16 | -2.82 | -7.38 |

---

## 5. Regions (HIW 4.7)

### 5.1 Overview

| dataset | regions_tested | depth1 | depth2 | z_above_2 | z_below_minus2 | STRONG_original | min_q_original | whole_cohort_AUROC | STRONG_corrected | min_q_corrected | z_above_2_corrected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ROSMAP | 334 | 77 | 257 | 34 | 4 | 0 | 0.151 | 0.595 | 0.0 | 0.138 | 39.0 |
| Diverse | 279 | 87 | 192 | 31 | 2 | 28 | 0.040 | 0.643 | 0.0 | 0.122 | 38.0 |
| Banner | 111 | 59 | 52 | 10 | 3 | 0 | 0.276 | 0.581 | 0.0 | 0.368 | 5.0 |
| BannerLFQ | 84 | 51 | 33 | 0 | 2 | 0 | 0.697 | 0.590 | — | — | — |

`original` = the pipeline's output (`{ds}_regions.csv`). `corrected` = same regions re-scored
with random subsets drawn only from people who have adjusted data and matched to the region's
usable size (`{ds}_regions_usable_null_check.csv`, HIW 9.1). The corrected check was not run for
BannerLFQ (it already had 0 regions with z > 2). For ROSMAP and Banner everyone has adjusted data,
so the two versions differ only by random draws.

### 5.2 ROSMAP

Whole-cohort adjusted AUROC: 0.5946

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  age_death 86.8-90 (middle third) | 2 | 154 | 0.6111 | 0.5817 | 3.36 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  age_at_visit_max 85.7-90 (middle third) | 2 | 152 | 0.6102 | 0.5817 | 3.25 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  pmi > 6.5 (upper half) | 2 | 104 | 0.6099 | 0.5743 | 3.25 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  cts_mmse30_lv <= 26 (lower half) | 2 | 102 | 0.6082 | 0.5743 | 3.09 | 0.1846 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  braaksc <= 4 (lower half) | 2 | 163 | 0.6089 | 0.5836 | 3.03 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  cts_mmse30_lv <= 23 (low third) | 2 | 65 | 0.6023 | 0.5659 | 2.92 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  apoe_genotype = 33 | 2 | 129 | 0.6085 | 0.5809 | 2.83 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  celltype_microglia <= -0.077 (lower half) | 2 | 112 | 0.6069 | 0.5777 | 2.82 | 0.1511 |
| celltype_oligodendrocyte <= -0.127 (lower half) | 1 | 200 | 0.6099 | 0.5874 | 2.80 | 0.1846 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  celltype_microglia <= -0.284 (low third) | 2 | 76 | 0.6050 | 0.5698 | 2.69 | 0.1994 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  educ <= 14 (low third) | 2 | 85 | 0.6036 | 0.5698 | 2.59 | 0.1994 |
| celltype_oligodendrocyte <= -0.371 (low third)  AND  cts_mmse30_lv <= 26 (lower half) | 2 | 68 | 0.6013 | 0.5688 | 2.55 | 0.1511 |
| ceradsc = 4  AND  Study = MAP | 2 | 67 | 0.6010 | 0.5688 | 2.53 | 0.1511 |
| celltype_oligodendrocyte <= -0.371 (low third)  AND  cts_mmse30_lv <= 23 (low third) | 2 | 41 | 0.5956 | 0.5571 | 2.50 | 0.1846 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  age_death > 89.4 (upper half) | 2 | 122 | 0.6020 | 0.5785 | 2.48 | 0.1994 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  apoe_genotype = 33 | 129 | 129 | 0.6085 | 0.5796 | 3.31 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  age_death 86.8-90 (middle third) | 154 | 154 | 0.6111 | 0.5814 | 3.29 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  age_at_visit_max 85.7-90 (middle third) | 152 | 152 | 0.6102 | 0.5814 | 3.19 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half) | 200 | 200 | 0.6099 | 0.5864 | 3.04 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  celltype_microglia <= -0.077 (lower half) | 112 | 112 | 0.6069 | 0.5768 | 3.00 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  pmi > 6.5 (upper half) | 104 | 104 | 0.6099 | 0.5756 | 2.97 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  celltype_microglia <= -0.284 (low third) | 76 | 76 | 0.6050 | 0.5693 | 2.90 | 0.2077 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  cts_mmse30_lv <= 26 (lower half) | 102 | 102 | 0.6082 | 0.5756 | 2.82 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  educ <= 14 (low third) | 85 | 85 | 0.6036 | 0.5693 | 2.79 | 0.2077 |
| celltype_oligodendrocyte <= -0.371 (low third)  AND  cts_mmse30_lv <= 23 (low third) | 41 | 41 | 0.5956 | 0.5574 | 2.73 | 0.1385 |
| educ <= 14 (low third)  AND  braaksc <= 4 (lower half) | 126 | 126 | 0.6032 | 0.5796 | 2.71 | 0.1662 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  braaksc <= 4 (lower half) | 163 | 163 | 0.6089 | 0.5838 | 2.70 | 0.1385 |
| celltype_oligodendrocyte <= -0.127 (lower half)  AND  cts_mmse30_lv <= 23 (low third) | 65 | 65 | 0.6023 | 0.5649 | 2.69 | 0.1385 |
| celltype_oligodendrocyte <= -0.371 (low third)  AND  cts_mmse30_lv <= 26 (lower half) | 68 | 68 | 0.6013 | 0.5691 | 2.65 | 0.1662 |
| ceradsc = 4  AND  Study = MAP | 67 | 67 | 0.6010 | 0.5691 | 2.62 | 0.2216 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_neuron <= 0.0321 (lower half) | 1 | 200 | 0.5685 | 0.5874 | -2.35 | 0.9930 |
| apoe_genotype = 33  AND  educ > 18 (high third) | 2 | 49 | 0.5286 | 0.5609 | -2.27 | 0.9950 |
| celltype_oligodendrocyte > -0.127 (upper half) | 1 | 200 | 0.5693 | 0.5874 | -2.25 | 0.9910 |
| apoe_genotype = 33  AND  celltype_neuron <= 0.0321 (lower half) | 2 | 122 | 0.5590 | 0.5785 | -2.05 | 0.9910 |
| celltype_oligodendrocyte > 0.164 (high third) | 1 | 133 | 0.5634 | 0.5809 | -1.80 | 0.9819 |
| celltype_astrocyte -0.186-0.235 (middle third) | 1 | 133 | 0.5644 | 0.5809 | -1.69 | 0.9697 |
| celltype_microglia > -0.077 (upper half) | 1 | 200 | 0.5752 | 0.5874 | -1.51 | 0.9676 |
| pmi 5.67-7.61 (middle third) | 1 | 131 | 0.5671 | 0.5809 | -1.42 | 0.9147 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| Study = MAP | 252 | 0.5938 | 0.5900 | 0.62 | 0.4621 |
| Study = ROS | 148 | 0.5747 | 0.5817 | -0.81 | 0.8791 |
| age_at_visit_max 85.7-90 (middle third) | 266 | 0.5928 | 0.5907 | 0.39 | 0.5567 |
| age_at_visit_max <= 85.7 (low third) | 134 | 0.5729 | 0.5809 | -0.83 | 0.8791 |
| age_at_visit_max <= 88.9 (lower half) | 200 | 0.5867 | 0.5874 | -0.08 | 0.7101 |
| age_at_visit_max > 88.9 (upper half) | 200 | 0.5890 | 0.5874 | 0.21 | 0.6316 |
| age_death 86.8-90 (middle third) | 266 | 0.5946 | 0.5907 | 0.71 | 0.4621 |
| age_death <= 86.8 (low third) | 134 | 0.5680 | 0.5809 | -1.32 | 0.9147 |
| age_death <= 89.4 (lower half) | 200 | 0.5869 | 0.5874 | -0.06 | 0.7101 |
| age_death > 89.4 (upper half) | 200 | 0.5879 | 0.5874 | 0.06 | 0.6719 |
| age_first_ad_dx 86.5-90 (middle third) | 88 | 0.5683 | 0.5732 | -0.44 | 0.7916 |
| age_first_ad_dx <= 86.5 (low third) | 44 | 0.5542 | 0.5571 | -0.19 | 0.7490 |
| age_first_ad_dx <= 89 (lower half) | 66 | 0.5539 | 0.5688 | -1.18 | 0.9147 |
| age_first_ad_dx > 89 (upper half) | 66 | 0.5666 | 0.5688 | -0.18 | 0.7414 |
| apoe_genotype = 23 | 57 | 0.5531 | 0.5659 | -1.02 | 0.8976 |
| apoe_genotype = 33 | 256 | 0.5979 | 0.5903 | 1.34 | 0.3210 |
| apoe_genotype = 34 | 73 | 0.5641 | 0.5688 | -0.38 | 0.7947 |
| braaksc 3-4 (middle third) | 131 | 0.5854 | 0.5809 | 0.46 | 0.5725 |
| braaksc <= 3 (low third) | 195 | 0.5864 | 0.5874 | -0.13 | 0.7188 |
| braaksc <= 4 (lower half) | 326 | 0.5971 | 0.5924 | 1.30 | 0.3549 |
| braaksc > 4 (high third) | 74 | 0.5543 | 0.5688 | -1.15 | 0.9147 |
| braaksc > 4 (upper half) | 74 | 0.5543 | 0.5688 | -1.15 | 0.9147 |
| ceradsc = 1 | 113 | 0.5675 | 0.5777 | -0.99 | 0.9097 |
| ceradsc = 2 | 138 | 0.5801 | 0.5817 | -0.16 | 0.7188 |
| ceradsc = 3 | 43 | 0.5602 | 0.5571 | 0.21 | 0.6298 |
| ceradsc = 4 | 106 | 0.5970 | 0.5777 | 1.86 | 0.2432 |
| cogdx = 1 | 168 | 0.5910 | 0.5836 | 0.85 | 0.4367 |
| cogdx = 2 | 97 | 0.5819 | 0.5743 | 0.69 | 0.4653 |
| cogdx = 4 | 109 | 0.5750 | 0.5777 | -0.26 | 0.7726 |
| cts_mmse30_first_ad_dx <= 18 (low third) | 48 | 0.5438 | 0.5609 | -1.20 | 0.9144 |
| cts_mmse30_first_ad_dx <= 21 (lower half) | 71 | 0.5612 | 0.5688 | -0.60 | 0.8449 |
| cts_mmse30_first_ad_dx > 21 (upper half) | 55 | 0.5649 | 0.5659 | -0.08 | 0.7188 |
| cts_mmse30_lv 23-28 (middle third) | 174 | 0.5814 | 0.5836 | -0.25 | 0.7490 |
| cts_mmse30_lv <= 23 (low third) | 136 | 0.5773 | 0.5817 | -0.44 | 0.7726 |
| cts_mmse30_lv <= 26 (lower half) | 210 | 0.5904 | 0.5869 | 0.52 | 0.5299 |
| cts_mmse30_lv > 26 (upper half) | 190 | 0.5805 | 0.5863 | -0.74 | 0.8616 |
| cts_mmse30_lv > 28 (high third) | 90 | 0.5830 | 0.5732 | 0.88 | 0.4367 |
| dcfdx_lv = 1 | 174 | 0.5944 | 0.5836 | 1.24 | 0.3693 |
| dcfdx_lv = 2 | 100 | 0.5798 | 0.5743 | 0.50 | 0.5176 |
| dcfdx_lv = 4 | 104 | 0.5680 | 0.5743 | -0.57 | 0.8337 |
| educ 14-18 (middle third) | 168 | 0.5758 | 0.5836 | -0.89 | 0.8791 |
| educ <= 14 (low third) | 158 | 0.5965 | 0.5836 | 1.55 | 0.3003 |
| educ <= 16 (lower half) | 247 | 0.5959 | 0.5900 | 0.96 | 0.4126 |
| educ > 16 (upper half) | 153 | 0.5751 | 0.5817 | -0.76 | 0.8722 |
| educ > 18 (high third) | 74 | 0.5618 | 0.5688 | -0.55 | 0.8394 |
| msex = 0 | 281 | 0.5931 | 0.5909 | 0.40 | 0.5674 |
| msex = 1 | 119 | 0.5689 | 0.5785 | -1.02 | 0.9078 |
| pmi 5.67-7.61 (middle third) | 131 | 0.5671 | 0.5809 | -1.42 | 0.9147 |
| pmi <= 5.67 (low third) | 134 | 0.5828 | 0.5809 | 0.19 | 0.6515 |
| pmi <= 6.5 (lower half) | 201 | 0.5807 | 0.5874 | -0.83 | 0.8791 |
| pmi > 6.5 (upper half) | 197 | 0.5903 | 0.5874 | 0.36 | 0.5749 |
| pmi > 7.61 (high third) | 133 | 0.5786 | 0.5809 | -0.24 | 0.7559 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.186-0.235 (middle third) | 133 | 0.5644 | 0.5809 | -1.69 | 0.9697 |
| celltype_astrocyte <= -0.186 (low third) | 134 | 0.5797 | 0.5809 | -0.13 | 0.7298 |
| celltype_astrocyte <= 0.0643 (lower half) | 200 | 0.5805 | 0.5874 | -0.85 | 0.8791 |
| celltype_astrocyte > 0.0643 (upper half) | 200 | 0.5777 | 0.5874 | -1.21 | 0.9147 |
| celltype_astrocyte > 0.235 (high third) | 133 | 0.5747 | 0.5809 | -0.64 | 0.8394 |
| celltype_endothelial -0.263-0.217 (middle third) | 133 | 0.5682 | 0.5809 | -1.30 | 0.9147 |
| celltype_endothelial <= -0.043 (lower half) | 200 | 0.5862 | 0.5874 | -0.14 | 0.7225 |
| celltype_endothelial <= -0.263 (low third) | 134 | 0.5802 | 0.5809 | -0.08 | 0.7188 |
| celltype_endothelial > -0.043 (upper half) | 200 | 0.5905 | 0.5874 | 0.39 | 0.5674 |
| celltype_endothelial > 0.217 (high third) | 133 | 0.5918 | 0.5809 | 1.11 | 0.3693 |
| celltype_microglia -0.284-0.197 (middle third) | 133 | 0.5807 | 0.5809 | -0.02 | 0.7152 |
| celltype_microglia <= -0.077 (lower half) | 200 | 0.5965 | 0.5874 | 1.14 | 0.3825 |
| celltype_microglia <= -0.284 (low third) | 134 | 0.5887 | 0.5809 | 0.80 | 0.4367 |
| celltype_microglia > -0.077 (upper half) | 200 | 0.5752 | 0.5874 | -1.51 | 0.9676 |
| celltype_microglia > 0.197 (high third) | 133 | 0.5752 | 0.5809 | -0.59 | 0.8337 |
| celltype_neuron -0.149-0.194 (middle third) | 133 | 0.5688 | 0.5809 | -1.25 | 0.9147 |
| celltype_neuron <= -0.149 (low third) | 134 | 0.5720 | 0.5809 | -0.92 | 0.8791 |
| celltype_neuron <= 0.0321 (lower half) | 200 | 0.5685 | 0.5874 | -2.35 | 0.9930 |
| celltype_neuron > 0.0321 (upper half) | 200 | 0.5950 | 0.5874 | 0.95 | 0.4067 |
| celltype_neuron > 0.194 (high third) | 133 | 0.5853 | 0.5809 | 0.45 | 0.5725 |
| celltype_oligodendrocyte -0.371-0.164 (middle third) | 133 | 0.5847 | 0.5809 | 0.38 | 0.5901 |
| celltype_oligodendrocyte <= -0.127 (lower half) | 200 | 0.6099 | 0.5874 | 2.80 | 0.1846 |
| celltype_oligodendrocyte <= -0.371 (low third) | 134 | 0.5951 | 0.5809 | 1.46 | 0.3003 |
| celltype_oligodendrocyte > -0.127 (upper half) | 200 | 0.5693 | 0.5874 | -2.25 | 0.9910 |
| celltype_oligodendrocyte > 0.164 (high third) | 133 | 0.5634 | 0.5809 | -1.80 | 0.9819 |

### 5.3 Diverse

Whole-cohort adjusted AUROC: 0.6431

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| ADoutcome = Other  AND  celltype_neuron <= -0.0632 (low third) | 2 | 57 | 0.6522 | 0.6011 | 3.00 | 0.0397 |
| PMI 5.17-10 (middle third)  AND  amyAny = 1 | 2 | 192 | 0.6580 | 0.6263 | 2.95 | 0.0397 |
| ADoutcome = Other  AND  celltype_endothelial > -0.0425 (upper half) | 2 | 75 | 0.6560 | 0.6080 | 2.90 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  amyCerad = Moderate/Probable/C2 | 2 | 48 | 0.6458 | 0.5971 | 2.77 | 0.0567 |
| ADoutcome = Other  AND  celltype_neuron <= 0.103 (lower half) | 2 | 86 | 0.6571 | 0.6128 | 2.72 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  isHispanic = FALSE | 2 | 252 | 0.6553 | 0.6308 | 2.71 | 0.0397 |
| ADoutcome = Other  AND  sex = female | 2 | 115 | 0.6553 | 0.6175 | 2.70 | 0.0397 |
| PMI 5.17-10 (middle third)  AND  celltype_oligodendrocyte -0.482-0.323 (middle third) | 2 | 107 | 0.6522 | 0.6158 | 2.69 | 0.0567 |
| ADoutcome = Other  AND  amyCerad = Moderate/Probable/C2 | 2 | 49 | 0.6431 | 0.5971 | 2.61 | 0.0731 |
| PMI 5.17-10 (middle third)  AND  celltype_neuron <= -0.0632 (low third) | 2 | 96 | 0.6525 | 0.6152 | 2.59 | 0.0397 |
| amyThal = Phase 3  AND  amyA = Thal Phase 3 | 2 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| amyA = Thal Phase 3 | 1 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| amyThal = Phase 3 | 1 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  derivedOutcomeBasedOnMayoDx = FALSE | 2 | 253 | 0.6532 | 0.6308 | 2.47 | 0.0397 |
| ADoutcome = Other  AND  derivedOutcomeBasedOnMayoDx = FALSE | 2 | 104 | 0.6504 | 0.6152 | 2.45 | 0.0567 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| ADoutcome = Other  AND  sex = female | 115 | 65 | 0.6553 | 0.6088 | 3.43 | 0.1220 |
| ADoutcome = Other  AND  celltype_neuron <= 0.103 (lower half) | 86 | 51 | 0.6571 | 0.6073 | 3.05 | 0.1220 |
| ADoutcome = Other  AND  celltype_endothelial > -0.0425 (upper half) | 75 | 46 | 0.6560 | 0.6073 | 2.98 | 0.1220 |
| PMI 5.17-10 (middle third)  AND  amyAny = 1 | 192 | 192 | 0.6580 | 0.6317 | 2.93 | 0.1220 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  celltype_astrocyte <= -0.213 (low third) | 86 | 64 | 0.6479 | 0.6088 | 2.88 | 0.1220 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  isHispanic = FALSE | 252 | 222 | 0.6553 | 0.6337 | 2.67 | 0.1220 |
| PMI 5.17-10 (middle third)  AND  celltype_neuron <= -0.0632 (low third) | 96 | 82 | 0.6525 | 0.6173 | 2.67 | 0.1220 |
| ADoutcome = Other  AND  amyAny = 1 | 89 | 89 | 0.6525 | 0.6198 | 2.53 | 0.1220 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  ageDeath 79-87.7 (middle third) | 101 | 83 | 0.6501 | 0.6173 | 2.48 | 0.1220 |
| amyA = Thal Phase 3  AND  amyAny = 1 | 40 | 40 | 0.6436 | 0.6002 | 2.40 | 0.1220 |
| amyThal = Phase 3  AND  amyAny = 1 | 40 | 40 | 0.6436 | 0.6002 | 2.40 | 0.1220 |
| amyThal = Phase 3 | 61 | 41 | 0.6435 | 0.6002 | 2.40 | 0.1220 |
| amyA = Thal Phase 3 | 61 | 41 | 0.6435 | 0.6002 | 2.40 | 0.1220 |
| amyA = Thal Phase 3  AND  derivedOutcomeBasedOnMayoDx = FALSE | 41 | 41 | 0.6435 | 0.6002 | 2.40 | 0.1220 |
| amyThal = Phase 3  AND  derivedOutcomeBasedOnMayoDx = FALSE | 41 | 41 | 0.6435 | 0.6002 | 2.40 | 0.1220 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| race = other | 1 | 214 | 0.5889 | 0.6281 | -4.09 | 1.0000 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  isHispanic = TRUE | 2 | 69 | 0.5706 | 0.6059 | -2.04 | 0.9771 |
| celltype_microglia -0.289-0.179 (middle third) | 1 | 321 | 0.6192 | 0.6341 | -1.85 | 0.9836 |
| celltype_astrocyte <= 0.0206 (lower half) | 1 | 482 | 0.6280 | 0.6387 | -1.82 | 0.9727 |
| celltype_oligodendrocyte > 0.323 (high third) | 1 | 321 | 0.6198 | 0.6341 | -1.79 | 0.9757 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  amyThal = Phase 5 | 2 | 72 | 0.5767 | 0.6059 | -1.69 | 0.9727 |
| ADoutcome = Other  AND  celltype_neuron > 0.277 (high third) | 2 | 78 | 0.5802 | 0.6080 | -1.68 | 0.9727 |
| celltype_endothelial <= -0.298 (low third) | 1 | 322 | 0.6206 | 0.6341 | -1.68 | 0.9727 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| ADoutcome = AD | 560 | 0.6393 | 0.6390 | 0.06 | 0.6253 |
| ADoutcome = Control | 209 | 0.6128 | 0.6281 | -1.60 | 0.9742 |
| ADoutcome = Other | 192 | 0.6504 | 0.6263 | 2.25 | 0.0567 |
| Braak = Stage I | 68 | 0.5997 | 0.6059 | -0.36 | 0.8082 |
| Braak = Stage II | 86 | 0.6055 | 0.6128 | -0.45 | 0.8269 |
| Braak = Stage III | 160 | 0.6281 | 0.6216 | 0.55 | 0.4523 |
| Braak = Stage IV | 164 | 0.6384 | 0.6216 | 1.42 | 0.1702 |
| Braak = Stage V | 191 | 0.6181 | 0.6263 | -0.76 | 0.8718 |
| Braak = Stage VI | 239 | 0.6165 | 0.6309 | -1.51 | 0.9727 |
| PMI 5.17-10 (middle third) | 275 | 0.6527 | 0.6324 | 2.44 | 0.0397 |
| PMI <= 5.17 (low third) | 274 | 0.6314 | 0.6318 | -0.04 | 0.6671 |
| PMI <= 6.83 (lower half) | 409 | 0.6411 | 0.6371 | 0.57 | 0.4576 |
| PMI > 10 (high third) | 269 | 0.6188 | 0.6318 | -1.46 | 0.9727 |
| PMI > 6.83 (upper half) | 409 | 0.6381 | 0.6371 | 0.14 | 0.6464 |
| ageDeath 79-87.7 (middle third) | 303 | 0.6404 | 0.6337 | 0.85 | 0.3577 |
| ageDeath <= 79 (low third) | 339 | 0.6294 | 0.6350 | -0.73 | 0.8802 |
| ageDeath <= 83.6 (lower half) | 482 | 0.6315 | 0.6387 | -1.21 | 0.9635 |
| ageDeath > 83.6 (upper half) | 481 | 0.6402 | 0.6387 | 0.27 | 0.5949 |
| ageDeath > 87.7 (high third) | 321 | 0.6294 | 0.6341 | -0.58 | 0.8596 |
| amyA = Thal Phase 1 or 2 | 64 | 0.5944 | 0.6011 | -0.40 | 0.7914 |
| amyA = Thal Phase 3 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| amyA = Thal Phase 4 or 5 | 273 | 0.6184 | 0.6318 | -1.50 | 0.9727 |
| amyAny = 0 | 167 | 0.6105 | 0.6267 | -1.41 | 0.9727 |
| amyAny = 1 | 509 | 0.6477 | 0.6386 | 1.86 | 0.1185 |
| amyCerad = Frequent/Definite/C3 | 326 | 0.6307 | 0.6338 | -0.40 | 0.8082 |
| amyCerad = Moderate/Probable/C2 | 138 | 0.6471 | 0.6209 | 2.07 | 0.0731 |
| amyCerad = None/No AD/C0 | 167 | 0.6105 | 0.6267 | -1.41 | 0.9727 |
| amyCerad = Sparse/Possible/C1 | 45 | 0.6140 | — | — | 0.2639 |
| amyThal = Phase 3 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| amyThal = Phase 4 | 52 | 0.5692 | 0.5971 | -1.59 | 0.9727 |
| amyThal = Phase 5 | 217 | 0.6144 | 0.6279 | -1.38 | 0.9727 |
| apoeGenotype = 23 | 73 | 0.5964 | 0.6059 | -0.55 | 0.8596 |
| apoeGenotype = 33 | 449 | 0.6460 | 0.6378 | 1.21 | 0.2505 |
| apoeGenotype = 34 | 286 | 0.6221 | 0.6330 | -1.26 | 0.9700 |
| apoeGenotype = 44 | 58 | 0.6118 | 0.6011 | 0.62 | 0.4045 |
| bScore = Braak Stage I-II | 154 | 0.6144 | 0.6243 | -0.83 | 0.9123 |
| bScore = Braak Stage III-IV | 324 | 0.6455 | 0.6341 | 1.43 | 0.1666 |
| bScore = Braak Stage V-VI | 430 | 0.6276 | 0.6368 | -1.44 | 0.9727 |
| cohort = CLINCOR | 65 | 0.6231 | 0.6011 | 1.29 | 0.2313 |
| cohort = Emory | 97 | 0.6198 | 0.6152 | 0.32 | 0.5908 |
| cohort = MAP | 69 | 0.6391 | 0.6059 | 1.92 | 0.0941 |
| cohort = Mayo Clinic | 248 | — | 0.6308 | — | 0.0397 |
| cohort = Mt Sinai Brain Bank | 205 | 0.6126 | 0.6281 | -1.35 | 0.9727 |
| cohort = ROS | 188 | 0.6250 | 0.6263 | -0.12 | 0.7295 |
| dataContributionGroup = Emory | 141 | 0.6178 | 0.6209 | -0.24 | 0.7704 |
| dataContributionGroup = MSSM | 183 | 0.6127 | 0.6270 | -1.25 | 0.9668 |
| dataContributionGroup = Mayo | 285 | — | 0.6324 | — | 0.0397 |
| dataContributionGroup = Rush | 355 | 0.6448 | 0.6356 | 1.32 | 0.2180 |
| derivedOutcomeBasedOnMayoDx = FALSE | 679 | 0.6431 | 0.6411 | 0.53 | 0.4776 |
| derivedOutcomeBasedOnMayoDx = TRUE | 285 | — | 0.6324 | — | 0.0397 |
| isHispanic = FALSE | 707 | 0.6405 | 0.6415 | -0.27 | 0.7941 |
| isHispanic = TRUE | 256 | 0.6263 | 0.6321 | -0.62 | 0.8492 |
| mayoDx = AD | 181 | — | 0.6270 | — | 0.0397 |
| mayoDx = Other | 78 | — | 0.6080 | — | 0.0397 |
| race = Black or African American | 267 | 0.6293 | 0.6318 | -0.28 | 0.7484 |
| race = White | 465 | 0.6414 | 0.6376 | 0.67 | 0.3823 |
| race = other | 214 | 0.5889 | 0.6281 | -4.09 | 1.0000 |
| reag = High Likelihood | 244 | 0.6226 | 0.6309 | -0.86 | 0.9123 |
| reag = Intermediate Likelihood | 209 | 0.6445 | 0.6281 | 1.72 | 0.1388 |
| reag = Low Likelihood | 191 | 0.6275 | 0.6263 | 0.11 | 0.6253 |
| sex = female | 567 | 0.6417 | 0.6401 | 0.33 | 0.5466 |
| sex = male | 397 | 0.6330 | 0.6366 | -0.51 | 0.8268 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.213-0.236 (middle third) | 321 | 0.6243 | 0.6341 | -1.22 | 0.9635 |
| celltype_astrocyte <= -0.213 (low third) | 322 | 0.6307 | 0.6341 | -0.43 | 0.8269 |
| celltype_astrocyte <= 0.0206 (lower half) | 482 | 0.6280 | 0.6387 | -1.82 | 0.9727 |
| celltype_astrocyte > 0.0206 (upper half) | 482 | 0.6468 | 0.6387 | 1.38 | 0.1965 |
| celltype_astrocyte > 0.236 (high third) | 321 | 0.6436 | 0.6341 | 1.18 | 0.2349 |
| celltype_endothelial -0.298-0.239 (middle third) | 321 | 0.6379 | 0.6341 | 0.48 | 0.4550 |
| celltype_endothelial <= -0.0425 (lower half) | 482 | 0.6309 | 0.6387 | -1.32 | 0.9668 |
| celltype_endothelial <= -0.298 (low third) | 322 | 0.6206 | 0.6341 | -1.68 | 0.9727 |
| celltype_endothelial > -0.0425 (upper half) | 482 | 0.6352 | 0.6387 | -0.60 | 0.8802 |
| celltype_endothelial > 0.239 (high third) | 321 | 0.6254 | 0.6341 | -1.09 | 0.9550 |
| celltype_microglia -0.289-0.179 (middle third) | 321 | 0.6192 | 0.6341 | -1.85 | 0.9836 |
| celltype_microglia <= -0.0479 (lower half) | 482 | 0.6432 | 0.6387 | 0.77 | 0.3919 |
| celltype_microglia <= -0.289 (low third) | 322 | 0.6456 | 0.6341 | 1.44 | 0.1666 |
| celltype_microglia > -0.0479 (upper half) | 482 | 0.6292 | 0.6387 | -1.60 | 0.9727 |
| celltype_microglia > 0.179 (high third) | 321 | 0.6330 | 0.6341 | -0.13 | 0.7733 |
| celltype_neuron -0.0632-0.277 (middle third) | 321 | 0.6245 | 0.6341 | -1.19 | 0.9635 |
| celltype_neuron <= -0.0632 (low third) | 322 | 0.6399 | 0.6341 | 0.72 | 0.3786 |
| celltype_neuron <= 0.103 (lower half) | 482 | 0.6360 | 0.6387 | -0.45 | 0.8492 |
| celltype_neuron > 0.103 (upper half) | 482 | 0.6407 | 0.6387 | 0.34 | 0.5278 |
| celltype_neuron > 0.277 (high third) | 321 | 0.6334 | 0.6341 | -0.09 | 0.7295 |
| celltype_oligodendrocyte -0.482-0.323 (middle third) | 321 | 0.6532 | 0.6341 | 2.38 | 0.0731 |
| celltype_oligodendrocyte <= -0.11 (lower half) | 482 | 0.6373 | 0.6387 | -0.24 | 0.7914 |
| celltype_oligodendrocyte <= -0.482 (low third) | 322 | 0.6240 | 0.6341 | -1.26 | 0.9635 |
| celltype_oligodendrocyte > -0.11 (upper half) | 482 | 0.6355 | 0.6387 | -0.53 | 0.8697 |
| celltype_oligodendrocyte > 0.323 (high third) | 321 | 0.6198 | 0.6341 | -1.79 | 0.9757 |

**The 28 regions flagged STRONG in the original run** (none remain STRONG with the corrected null)

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| ADoutcome = Other  AND  celltype_neuron <= -0.0632 (low third) | 2 | 57 | 0.6522 | 0.6011 | 3.00 | 0.0397 |
| PMI 5.17-10 (middle third)  AND  amyAny = 1 | 2 | 192 | 0.6580 | 0.6263 | 2.95 | 0.0397 |
| ADoutcome = Other  AND  celltype_endothelial > -0.0425 (upper half) | 2 | 75 | 0.6560 | 0.6080 | 2.90 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  amyCerad = Moderate/Probable/C2 | 2 | 48 | 0.6458 | 0.5971 | 2.77 | 0.0567 |
| ADoutcome = Other  AND  celltype_neuron <= 0.103 (lower half) | 2 | 86 | 0.6571 | 0.6128 | 2.72 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  isHispanic = FALSE | 2 | 252 | 0.6553 | 0.6308 | 2.71 | 0.0397 |
| ADoutcome = Other  AND  sex = female | 2 | 115 | 0.6553 | 0.6175 | 2.70 | 0.0397 |
| PMI 5.17-10 (middle third)  AND  celltype_oligodendrocyte -0.482-0.323 (middle third) | 2 | 107 | 0.6522 | 0.6158 | 2.69 | 0.0567 |
| ADoutcome = Other  AND  amyCerad = Moderate/Probable/C2 | 2 | 49 | 0.6431 | 0.5971 | 2.61 | 0.0731 |
| PMI 5.17-10 (middle third)  AND  celltype_neuron <= -0.0632 (low third) | 2 | 96 | 0.6525 | 0.6152 | 2.59 | 0.0397 |
| amyThal = Phase 3  AND  amyA = Thal Phase 3 | 2 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| amyA = Thal Phase 3 | 1 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| amyThal = Phase 3 | 1 | 61 | 0.6435 | 0.6011 | 2.49 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  derivedOutcomeBasedOnMayoDx = FALSE | 2 | 253 | 0.6532 | 0.6308 | 2.47 | 0.0397 |
| ADoutcome = Other  AND  derivedOutcomeBasedOnMayoDx = FALSE | 2 | 104 | 0.6504 | 0.6152 | 2.45 | 0.0567 |
| PMI 5.17-10 (middle third) | 1 | 275 | 0.6527 | 0.6324 | 2.44 | 0.0397 |
| ADoutcome = Other  AND  amyAny = 1 | 2 | 89 | 0.6525 | 0.6128 | 2.44 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  ageDeath 79-87.7 (middle third) | 2 | 101 | 0.6501 | 0.6152 | 2.42 | 0.0567 |
| PMI 5.17-10 (middle third)  AND  derivedOutcomeBasedOnMayoDx = FALSE | 2 | 246 | 0.6527 | 0.6308 | 2.42 | 0.0567 |
| celltype_oligodendrocyte -0.482-0.323 (middle third) | 1 | 321 | 0.6532 | 0.6341 | 2.38 | 0.0731 |
| PMI 5.17-10 (middle third)  AND  ageDeath <= 79 (low third) | 2 | 78 | 0.6473 | 0.6080 | 2.38 | 0.0731 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  amyAny = 1 | 2 | 194 | 0.6512 | 0.6263 | 2.32 | 0.0567 |
| PMI 5.17-10 (middle third)  AND  celltype_astrocyte > 0.0206 (upper half) | 2 | 141 | 0.6493 | 0.6209 | 2.25 | 0.0567 |
| ADoutcome = Other | 1 | 192 | 0.6504 | 0.6263 | 2.25 | 0.0567 |
| PMI 5.17-10 (middle third)  AND  reag = Intermediate Likelihood | 2 | 68 | 0.6441 | 0.6059 | 2.21 | 0.0397 |
| PMI 5.17-10 (middle third)  AND  celltype_endothelial -0.298-0.239 (middle third) | 2 | 98 | 0.6470 | 0.6152 | 2.20 | 0.0567 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  celltype_astrocyte <= -0.213 (low third) | 2 | 86 | 0.6479 | 0.6128 | 2.16 | 0.0731 |
| amyCerad = Moderate/Probable/C2 | 1 | 138 | 0.6471 | 0.6209 | 2.07 | 0.0731 |

**Regions with no usable people after adjustment** (AUROC NaN, p set to the 1/201 floor; HIW 9.2)

| region | n | auroc | null_mean | p_vs_random_same_size | q_fdr |
|---|---|---|---|---|---|
| dataContributionGroup = Mayo | 285 | — | 0.6324 | 0.0050 | 0.0397 |
| cohort = Mayo Clinic | 248 | — | 0.6308 | 0.0050 | 0.0397 |
| mayoDx = AD | 181 | — | 0.6270 | 0.0050 | 0.0397 |
| mayoDx = Other | 78 | — | 0.6080 | 0.0050 | 0.0397 |
| derivedOutcomeBasedOnMayoDx = TRUE | 285 | — | 0.6324 | 0.0050 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  dataContributionGroup = Mayo | 68 | — | 0.6059 | 0.0050 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  cohort = Mayo Clinic | 59 | — | 0.6011 | 0.0050 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  race = other | 53 | — | 0.5971 | 0.0050 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  mayoDx = AD | 49 | — | 0.5971 | 0.0050 | 0.0397 |
| celltype_oligodendrocyte -0.482-0.323 (middle third)  AND  derivedOutcomeBasedOnMayoDx = TRUE | 68 | — | 0.6059 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  dataContributionGroup = Mayo | 88 | — | 0.6128 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  cohort = Mayo Clinic | 76 | — | 0.6080 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  race = other | 65 | — | 0.6011 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  isHispanic = TRUE | 74 | — | 0.6059 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  mayoDx = Other | 78 | — | 0.6080 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  bScore = Braak Stage I-II | 48 | — | 0.5971 | 0.0050 | 0.0397 |
| ADoutcome = Other  AND  derivedOutcomeBasedOnMayoDx = TRUE | 88 | — | 0.6128 | 0.0050 | 0.0397 |

`Diverse_best_region_W.csv` was refitted in `ADoutcome = Other AND celltype_neuron <= -0.0632 (low third)` (n = 57 listed), the original top STRONG region.

### 5.4 Banner

Whole-cohort adjusted AUROC: 0.5811

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_endothelial <= -0.0477 (lower half)  AND  celltype_oligodendrocyte > -0.129 (upper half) | 2 | 54 | 0.5897 | 0.5604 | 2.51 | 0.2761 |
| apoeGenotype = e3-3  AND  celltype_microglia > 0.00693 (upper half) | 2 | 44 | 0.5890 | 0.5572 | 2.41 | 0.2761 |
| lastMMSE > 28 (high third)  AND  PlaqueTotal <= 12.1 (lower half) | 2 | 44 | 0.5878 | 0.5572 | 2.32 | 0.2761 |
| apoeGenotype = e3-3  AND  celltype_endothelial <= -0.0477 (lower half) | 2 | 47 | 0.5865 | 0.5604 | 2.23 | 0.2761 |
| apoeGenotype = e3-3  AND  celltype_oligodendrocyte > -0.129 (upper half) | 2 | 45 | 0.5849 | 0.5572 | 2.09 | 0.2761 |
| celltype_endothelial <= -0.261 (low third) | 1 | 64 | 0.5875 | 0.5640 | 2.08 | 0.2761 |
| lastMMSE > 28 (high third) | 1 | 48 | 0.5846 | 0.5604 | 2.07 | 0.2761 |
| lastMMSE > 28 (high third)  AND  diagnosis = control | 2 | 48 | 0.5846 | 0.5604 | 2.07 | 0.2761 |
| celltype_endothelial <= -0.0477 (lower half)  AND  TangleTotal <= 6 (low third) | 2 | 40 | 0.5840 | 0.5572 | 2.02 | 0.2761 |
| celltype_endothelial <= -0.261 (low third)  AND  pmi <= 3 (lower half) | 2 | 41 | 0.5839 | 0.5572 | 2.02 | 0.2761 |
| lastMMSE > 28 (high third)  AND  TangleTotal <= 8 (lower half) | 2 | 44 | 0.5827 | 0.5572 | 1.93 | 0.2761 |
| celltype_oligodendrocyte -0.462-0.276 (middle third)  AND  diagnosis = control | 2 | 41 | 0.5824 | 0.5572 | 1.91 | 0.2761 |
| apoeGenotype = e3-3  AND  TangleTotal <= 6 (low third) | 2 | 40 | 0.5822 | 0.5572 | 1.89 | 0.2761 |
| apoeGenotype = e3-3  AND  ageDeath <= 86 (lower half) | 2 | 44 | 0.5820 | 0.5572 | 1.88 | 0.2761 |
| celltype_endothelial <= -0.0477 (lower half)  AND  pmi <= 3 (lower half) | 2 | 57 | 0.5844 | 0.5640 | 1.81 | 0.3797 |

**Top 15 by z (corrected null)** — `n_listed` = people in the region, `n_usable` = those with adjusted data

| region | n_listed | n_usable | auroc | null_usable | z_fixed | q_fixed |
|---|---|---|---|---|---|---|
| apoeGenotype = e3-3  AND  celltype_microglia > 0.00693 (upper half) | 44 | 44 | 0.5890 | 0.5572 | 2.32 | 0.3682 |
| celltype_endothelial <= -0.0477 (lower half)  AND  celltype_oligodendrocyte > -0.129 (upper half) | 54 | 54 | 0.5897 | 0.5615 | 2.31 | 0.3682 |
| lastMMSE > 28 (high third)  AND  PlaqueTotal <= 12.1 (lower half) | 44 | 44 | 0.5878 | 0.5572 | 2.23 | 0.3682 |
| apoeGenotype = e3-3  AND  celltype_endothelial <= -0.0477 (lower half) | 47 | 47 | 0.5865 | 0.5615 | 2.04 | 0.3682 |
| apoeGenotype = e3-3  AND  celltype_oligodendrocyte > -0.129 (upper half) | 45 | 45 | 0.5849 | 0.5572 | 2.02 | 0.3682 |
| celltype_endothelial <= -0.261 (low third) | 64 | 64 | 0.5875 | 0.5639 | 1.98 | 0.3682 |
| celltype_endothelial <= -0.0477 (lower half)  AND  TangleTotal <= 6 (low third) | 40 | 40 | 0.5840 | 0.5572 | 1.95 | 0.3682 |
| celltype_endothelial <= -0.261 (low third)  AND  pmi <= 3 (lower half) | 41 | 41 | 0.5839 | 0.5572 | 1.95 | 0.3682 |
| lastMMSE > 28 (high third) | 48 | 48 | 0.5846 | 0.5615 | 1.89 | 0.3682 |
| lastMMSE > 28 (high third)  AND  diagnosis = control | 48 | 48 | 0.5846 | 0.5615 | 1.89 | 0.3682 |
| lastMMSE > 28 (high third)  AND  TangleTotal <= 8 (lower half) | 44 | 44 | 0.5827 | 0.5572 | 1.86 | 0.3682 |
| celltype_oligodendrocyte -0.462-0.276 (middle third)  AND  diagnosis = control | 41 | 41 | 0.5824 | 0.5572 | 1.84 | 0.3682 |
| apoeGenotype = e3-3  AND  TangleTotal <= 6 (low third) | 40 | 40 | 0.5822 | 0.5572 | 1.82 | 0.3682 |
| apoeGenotype = e3-3  AND  ageDeath <= 86 (lower half) | 44 | 44 | 0.5820 | 0.5572 | 1.81 | 0.3682 |
| celltype_endothelial <= -0.0477 (lower half)  AND  pmi <= 3 (lower half) | 57 | 57 | 0.5844 | 0.5639 | 1.72 | 0.3682 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| lastMMSE <= 18 (low third) | 1 | 67 | 0.5405 | 0.5654 | -2.59 | 0.9950 |
| celltype_oligodendrocyte <= -0.462 (low third) | 1 | 64 | 0.5377 | 0.5640 | -2.34 | 0.9950 |
| lastMMSE <= 25 (lower half) | 1 | 99 | 0.5561 | 0.5732 | -2.10 | 0.9950 |
| celltype_endothelial > -0.0477 (upper half) | 1 | 95 | 0.5570 | 0.5732 | -2.00 | 0.9950 |
| celltype_endothelial -0.261-0.125 (middle third) | 1 | 63 | 0.5436 | 0.5640 | -1.82 | 0.9950 |
| TangleTotal > 8 (upper half) | 1 | 92 | 0.5571 | 0.5709 | -1.57 | 0.9742 |
| diagnosis = Alzheimer Disease | 1 | 90 | 0.5575 | 0.5709 | -1.53 | 0.9742 |
| TangleTotal > 11.8 (high third) | 1 | 63 | 0.5486 | 0.5640 | -1.37 | 0.9742 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| Braak = 4 | 57 | 0.5642 | 0.5640 | 0.02 | 0.6681 |
| Braak = 5 | 48 | 0.5619 | 0.5604 | 0.13 | 0.6420 |
| CERAD = 3 | 106 | 0.5698 | 0.5739 | -0.55 | 0.8713 |
| PlaqueTotal 8-13.5 (middle third) | 69 | 0.5577 | 0.5654 | -0.81 | 0.9166 |
| PlaqueTotal <= 12.1 (lower half) | 95 | 0.5764 | 0.5732 | 0.40 | 0.5522 |
| PlaqueTotal <= 8 (low third) | 66 | 0.5681 | 0.5654 | 0.28 | 0.6311 |
| PlaqueTotal > 12.1 (upper half) | 95 | 0.5626 | 0.5732 | -1.31 | 0.9742 |
| PlaqueTotal > 13.5 (high third) | 55 | 0.5648 | 0.5640 | 0.07 | 0.6420 |
| TangleTotal 6-11.8 (middle third) | 63 | 0.5621 | 0.5640 | -0.17 | 0.7035 |
| TangleTotal <= 6 (low third) | 64 | 0.5765 | 0.5640 | 1.11 | 0.4682 |
| TangleTotal <= 8 (lower half) | 98 | 0.5817 | 0.5732 | 1.05 | 0.4682 |
| TangleTotal > 11.8 (high third) | 63 | 0.5486 | 0.5640 | -1.37 | 0.9742 |
| TangleTotal > 8 (upper half) | 92 | 0.5571 | 0.5709 | -1.57 | 0.9742 |
| ageDeath 82-88 (middle third) | 54 | 0.5566 | 0.5604 | -0.32 | 0.8125 |
| ageDeath <= 82 (low third) | 73 | 0.5721 | 0.5654 | 0.70 | 0.5148 |
| ageDeath <= 86 (lower half) | 104 | 0.5782 | 0.5732 | 0.63 | 0.4844 |
| ageDeath > 86 (upper half) | 86 | 0.5644 | 0.5709 | -0.74 | 0.8930 |
| ageDeath > 88 (high third) | 63 | 0.5571 | 0.5640 | -0.62 | 0.8713 |
| apoeGenotype = e3-3 | 92 | 0.5816 | 0.5709 | 1.21 | 0.4682 |
| apoeGenotype = e3-4 | 62 | 0.5496 | 0.5640 | -1.28 | 0.9742 |
| diagnosis = Alzheimer Disease | 90 | 0.5575 | 0.5709 | -1.53 | 0.9742 |
| diagnosis = control | 100 | 0.5807 | 0.5732 | 0.93 | 0.4682 |
| lastMMSE 18-28 (middle third) | 75 | 0.5689 | 0.5681 | 0.09 | 0.6420 |
| lastMMSE <= 18 (low third) | 67 | 0.5405 | 0.5654 | -2.59 | 0.9950 |
| lastMMSE <= 25 (lower half) | 99 | 0.5561 | 0.5732 | -2.10 | 0.9950 |
| lastMMSE > 25 (upper half) | 91 | 0.5799 | 0.5709 | 1.03 | 0.4682 |
| lastMMSE > 28 (high third) | 48 | 0.5846 | 0.5604 | 2.07 | 0.2761 |
| pmi 2.5-3.25 (middle third) | 66 | 0.5663 | 0.5654 | 0.09 | 0.6420 |
| pmi <= 2.5 (low third) | 65 | 0.5630 | 0.5640 | -0.09 | 0.7035 |
| pmi <= 3 (lower half) | 117 | 0.5796 | 0.5747 | 0.82 | 0.4832 |
| pmi > 3 (upper half) | 73 | 0.5540 | 0.5654 | -1.19 | 0.9742 |
| pmi > 3.25 (high third) | 59 | 0.5562 | 0.5640 | -0.69 | 0.8860 |
| sex = female | 81 | 0.5606 | 0.5681 | -0.79 | 0.9089 |
| sex = male | 109 | 0.5808 | 0.5739 | 0.93 | 0.4682 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.218-0.213 (middle third) | 63 | 0.5702 | 0.5640 | 0.54 | 0.5435 |
| celltype_astrocyte <= -0.218 (low third) | 64 | 0.5594 | 0.5640 | -0.41 | 0.8125 |
| celltype_astrocyte <= 0.00686 (lower half) | 95 | 0.5753 | 0.5732 | 0.26 | 0.6420 |
| celltype_astrocyte > 0.00686 (upper half) | 95 | 0.5667 | 0.5732 | -0.80 | 0.9089 |
| celltype_astrocyte > 0.213 (high third) | 63 | 0.5509 | 0.5640 | -1.17 | 0.9742 |
| celltype_endothelial -0.261-0.125 (middle third) | 63 | 0.5436 | 0.5640 | -1.82 | 0.9950 |
| celltype_endothelial <= -0.0477 (lower half) | 95 | 0.5825 | 0.5732 | 1.16 | 0.4682 |
| celltype_endothelial <= -0.261 (low third) | 64 | 0.5875 | 0.5640 | 2.08 | 0.2761 |
| celltype_endothelial > -0.0477 (upper half) | 95 | 0.5570 | 0.5732 | -2.00 | 0.9950 |
| celltype_endothelial > 0.125 (high third) | 63 | 0.5536 | 0.5640 | -0.92 | 0.9298 |
| celltype_microglia -0.272-0.252 (middle third) | 63 | 0.5655 | 0.5640 | 0.14 | 0.6420 |
| celltype_microglia <= -0.272 (low third) | 64 | 0.5619 | 0.5640 | -0.19 | 0.7082 |
| celltype_microglia <= 0.00693 (lower half) | 95 | 0.5691 | 0.5732 | -0.50 | 0.8713 |
| celltype_microglia > 0.00693 (upper half) | 95 | 0.5764 | 0.5732 | 0.40 | 0.5522 |
| celltype_microglia > 0.252 (high third) | 63 | 0.5754 | 0.5640 | 1.01 | 0.4682 |
| celltype_neuron -0.147-0.187 (middle third) | 63 | 0.5730 | 0.5640 | 0.79 | 0.4682 |
| celltype_neuron <= -0.147 (low third) | 64 | 0.5556 | 0.5640 | -0.75 | 0.8930 |
| celltype_neuron <= -1.75e-05 (lower half) | 95 | 0.5678 | 0.5732 | -0.67 | 0.8930 |
| celltype_neuron > -1.75e-05 (upper half) | 95 | 0.5804 | 0.5732 | 0.90 | 0.4682 |
| celltype_neuron > 0.187 (high third) | 63 | 0.5723 | 0.5640 | 0.73 | 0.4844 |
| celltype_oligodendrocyte -0.462-0.276 (middle third) | 63 | 0.5796 | 0.5640 | 1.38 | 0.4267 |
| celltype_oligodendrocyte <= -0.129 (lower half) | 95 | 0.5625 | 0.5732 | -1.32 | 0.9742 |
| celltype_oligodendrocyte <= -0.462 (low third) | 64 | 0.5377 | 0.5640 | -2.34 | 0.9950 |
| celltype_oligodendrocyte > -0.129 (upper half) | 95 | 0.5806 | 0.5732 | 0.91 | 0.4682 |
| celltype_oligodendrocyte > 0.276 (high third) | 63 | 0.5717 | 0.5640 | 0.68 | 0.4844 |

### 5.5 BannerLFQ

Whole-cohort adjusted AUROC: 0.5896

**Top 15 by z (original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_neuron <= -0.182 (low third)  AND  CERAD = 3 | 2 | 46 | 0.5960 | 0.5722 | 1.98 | 0.6965 |
| celltype_astrocyte <= -0.0155 (lower half)  AND  CERAD = 3 | 2 | 41 | 0.5911 | 0.5652 | 1.90 | 0.6965 |
| sex = female | 1 | 81 | 0.5933 | 0.5770 | 1.69 | 0.6965 |
| celltype_neuron <= -0.182 (low third) | 1 | 64 | 0.5921 | 0.5744 | 1.62 | 0.6965 |
| sex = female  AND  diagnosis = control | 2 | 44 | 0.5871 | 0.5652 | 1.61 | 0.6965 |
| pmi <= 2.5 (low third) | 1 | 65 | 0.5915 | 0.5744 | 1.56 | 0.7164 |
| sex = female  AND  celltype_oligodendrocyte <= -0.073 (lower half) | 2 | 44 | 0.5848 | 0.5652 | 1.44 | 0.6965 |
| sex = female  AND  TangleTotal <= 8 (lower half) | 2 | 45 | 0.5834 | 0.5652 | 1.33 | 0.7684 |
| celltype_astrocyte <= -0.0155 (lower half)  AND  apoeGenotype = e3-3 | 2 | 56 | 0.5888 | 0.5744 | 1.32 | 0.7684 |
| celltype_astrocyte <= -0.0155 (lower half) | 1 | 95 | 0.5880 | 0.5795 | 1.21 | 0.7684 |
| sex = female  AND  pmi <= 3 (lower half) | 2 | 50 | 0.5864 | 0.5722 | 1.18 | 0.7684 |
| sex = female  AND  CERAD = 3 | 2 | 41 | 0.5804 | 0.5652 | 1.12 | 0.7684 |
| sex = female  AND  lastMMSE > 25 (upper half) | 2 | 41 | 0.5800 | 0.5652 | 1.09 | 0.7684 |
| sex = female  AND  celltype_astrocyte <= -0.0155 (lower half) | 2 | 43 | 0.5796 | 0.5652 | 1.05 | 0.7684 |
| celltype_astrocyte <= -0.0155 (lower half)  AND  celltype_oligodendrocyte > -0.073 (upper half) | 2 | 51 | 0.5847 | 0.5722 | 1.04 | 0.7684 |

**Bottom 8 by z (C least visible, original)**

| region | depth | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|---|
| celltype_neuron -0.182-0.263 (middle third) | 1 | 63 | 0.5461 | 0.5744 | -2.59 | 0.9851 |
| sex = male | 1 | 109 | 0.5649 | 0.5814 | -2.51 | 0.9851 |
| pmi > 3.25 (high third) | 1 | 59 | 0.5598 | 0.5744 | -1.33 | 0.9377 |
| celltype_astrocyte <= -0.0155 (lower half)  AND  PlaqueTotal <= 8 (low third) | 2 | 47 | 0.5566 | 0.5722 | -1.30 | 0.9377 |
| celltype_oligodendrocyte -0.464-0.391 (middle third) | 1 | 63 | 0.5608 | 0.5744 | -1.25 | 0.9377 |
| celltype_neuron > 0.0182 (upper half) | 1 | 95 | 0.5712 | 0.5795 | -1.18 | 0.9377 |
| ageDeath 82-88 (middle third) | 1 | 54 | 0.5586 | 0.5722 | -1.14 | 0.9377 |
| celltype_astrocyte <= -0.0155 (lower half)  AND  sex = male | 2 | 52 | 0.5596 | 0.5722 | -1.05 | 0.9377 |

**All single-variable clinical/pathology regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| Braak = 4 | 57 | 0.5778 | 0.5744 | 0.31 | 0.7684 |
| Braak = 5 | 48 | 0.5637 | 0.5722 | -0.71 | 0.9217 |
| CERAD = 3 | 106 | 0.5845 | 0.5814 | 0.48 | 0.7684 |
| PlaqueTotal 8-13.5 (middle third) | 69 | 0.5757 | 0.5750 | 0.07 | 0.7684 |
| PlaqueTotal <= 12.1 (lower half) | 95 | 0.5758 | 0.5795 | -0.52 | 0.8857 |
| PlaqueTotal <= 8 (low third) | 66 | 0.5654 | 0.5750 | -1.01 | 0.9293 |
| PlaqueTotal > 12.1 (upper half) | 95 | 0.5815 | 0.5795 | 0.28 | 0.7684 |
| PlaqueTotal > 13.5 (high third) | 55 | 0.5818 | 0.5744 | 0.68 | 0.7684 |
| TangleTotal 6-11.8 (middle third) | 63 | 0.5796 | 0.5744 | 0.48 | 0.7684 |
| TangleTotal <= 6 (low third) | 64 | 0.5735 | 0.5744 | -0.08 | 0.7684 |
| TangleTotal <= 8 (lower half) | 98 | 0.5782 | 0.5795 | -0.19 | 0.7684 |
| TangleTotal > 11.8 (high third) | 63 | 0.5824 | 0.5744 | 0.74 | 0.7684 |
| TangleTotal > 8 (upper half) | 92 | 0.5811 | 0.5775 | 0.41 | 0.7684 |
| ageDeath 82-88 (middle third) | 54 | 0.5586 | 0.5722 | -1.14 | 0.9377 |
| ageDeath <= 82 (low third) | 73 | 0.5820 | 0.5750 | 0.74 | 0.7684 |
| ageDeath <= 86 (lower half) | 104 | 0.5816 | 0.5795 | 0.29 | 0.7684 |
| ageDeath > 86 (upper half) | 86 | 0.5782 | 0.5775 | 0.08 | 0.7684 |
| ageDeath > 88 (high third) | 63 | 0.5757 | 0.5744 | 0.12 | 0.7684 |
| apoeGenotype = e3-3 | 92 | 0.5716 | 0.5775 | -0.68 | 0.9217 |
| apoeGenotype = e3-4 | 62 | 0.5770 | 0.5744 | 0.24 | 0.7684 |
| diagnosis = Alzheimer Disease | 90 | 0.5810 | 0.5775 | 0.40 | 0.7684 |
| diagnosis = control | 100 | 0.5784 | 0.5795 | -0.15 | 0.7684 |
| lastMMSE 18-28 (middle third) | 75 | 0.5799 | 0.5770 | 0.30 | 0.7684 |
| lastMMSE <= 18 (low third) | 67 | 0.5745 | 0.5750 | -0.05 | 0.7684 |
| lastMMSE <= 25 (lower half) | 99 | 0.5821 | 0.5795 | 0.37 | 0.7684 |
| lastMMSE > 25 (upper half) | 91 | 0.5704 | 0.5775 | -0.81 | 0.9217 |
| lastMMSE > 28 (high third) | 48 | 0.5762 | 0.5722 | 0.33 | 0.7684 |
| pmi 2.5-3.25 (middle third) | 66 | 0.5726 | 0.5750 | -0.25 | 0.7960 |
| pmi <= 2.5 (low third) | 65 | 0.5915 | 0.5744 | 1.56 | 0.7164 |
| pmi <= 3 (lower half) | 117 | 0.5850 | 0.5821 | 0.52 | 0.7684 |
| pmi > 3 (upper half) | 73 | 0.5666 | 0.5750 | -0.89 | 0.9217 |
| pmi > 3.25 (high third) | 59 | 0.5598 | 0.5744 | -1.33 | 0.9377 |
| sex = female | 81 | 0.5933 | 0.5770 | 1.69 | 0.6965 |
| sex = male | 109 | 0.5649 | 0.5814 | -2.51 | 0.9851 |

**All single-variable cell-type regions (original)**

| region | n | auroc | null_mean | z | q_fdr |
|---|---|---|---|---|---|
| celltype_astrocyte -0.237-0.251 (middle third) | 63 | 0.5847 | 0.5744 | 0.94 | 0.7684 |
| celltype_astrocyte <= -0.0155 (lower half) | 95 | 0.5880 | 0.5795 | 1.21 | 0.7684 |
| celltype_astrocyte <= -0.237 (low third) | 64 | 0.5729 | 0.5744 | -0.14 | 0.7684 |
| celltype_astrocyte > -0.0155 (upper half) | 95 | 0.5742 | 0.5795 | -0.75 | 0.9217 |
| celltype_astrocyte > 0.251 (high third) | 63 | 0.5637 | 0.5744 | -0.97 | 0.9293 |
| celltype_endothelial <= -0.0629 (lower half) | 51 | 0.5721 | 0.5722 | -0.01 | 0.7684 |
| celltype_endothelial > -0.0629 (upper half) | 50 | 0.5631 | 0.5722 | -0.76 | 0.9217 |
| celltype_neuron -0.182-0.263 (middle third) | 63 | 0.5461 | 0.5744 | -2.59 | 0.9851 |
| celltype_neuron <= -0.182 (low third) | 64 | 0.5921 | 0.5744 | 1.62 | 0.6965 |
| celltype_neuron <= 0.0182 (lower half) | 95 | 0.5776 | 0.5795 | -0.27 | 0.7972 |
| celltype_neuron > 0.0182 (upper half) | 95 | 0.5712 | 0.5795 | -1.18 | 0.9377 |
| celltype_neuron > 0.263 (high third) | 63 | 0.5696 | 0.5744 | -0.44 | 0.8857 |
| celltype_oligodendrocyte -0.464-0.391 (middle third) | 63 | 0.5608 | 0.5744 | -1.25 | 0.9377 |
| celltype_oligodendrocyte <= -0.073 (lower half) | 95 | 0.5816 | 0.5795 | 0.30 | 0.7684 |
| celltype_oligodendrocyte <= -0.464 (low third) | 64 | 0.5831 | 0.5744 | 0.80 | 0.7684 |
| celltype_oligodendrocyte > -0.073 (upper half) | 95 | 0.5789 | 0.5795 | -0.09 | 0.7684 |
| celltype_oligodendrocyte > 0.391 (high third) | 63 | 0.5752 | 0.5744 | 0.07 | 0.7684 |

---

## 6. The graphs (HIW 4.8–4.10)

### 6.1 Edge counts and sizes

| dataset | prior_edges | testable_n20 | W_positive | q_below_05_any_direction | data_supported | pct_of_688 | median_n_W_positive | min_n_W_positive |
|---|---|---|---|---|---|---|---|---|
| ROSMAP | 688 | 688 | 505 | 264 | 172 | 25.0000 | 400 | 184 |
| Diverse | 688 | 688 | 507 | 483 | 354 | 51.5000 | 676 | 455 |
| Banner | 688 | 688 | 505 | 152 | 101 | 14.7000 | 190 | 34 |
| BannerLFQ | 688 | 404 | 310 | 184 | 135 | 19.6000 | 190 | 0 |

| dataset | W_median_all | W_max | W_median_supported | W_min_supported | beta_median_supported | beta_min | beta_max | frac_beta_positive | n_abs_beta_above_2 |
|---|---|---|---|---|---|---|---|---|---|
| ROSMAP | 0.0160 | 0.4047 | 0.0524 | 0.0077 | 0.0881 | -0.4119 | 1.1784 | 0.7560 | 0 |
| Diverse | 0.0283 | 0.4520 | 0.0461 | 0.0052 | 0.0523 | -1.2715 | 2.2884 | 0.7260 | 1 |
| Banner | 0.0204 | 0.5336 | 0.0750 | 0.0116 | 0.1408 | -1.6528 | 4.1249 | 0.7430 | 3 |
| BannerLFQ | 0.0368 | 0.3466 | 0.0799 | 0.0133 | 0.1259 | -0.7529 | 1.2414 | 0.6960 | 0 |

### 6.2 How each prior pair's arrow and weight were set (HIW 4.8)

| | one-way pair, direction test agrees with C (full weight) | one-way pair, direction test disagrees (kept in C's direction at 0.3 weight) | two-way pair (direction test picks the arrow) | n < 20 (placeholder 0.1·C) |
|---|---|---|---|---|
| ROSMAP | 141 | 185 | 181 | 0 |
| Diverse | 149 | 177 | 181 | 0 |
| Banner | 154 | 172 | 181 | 0 |
| BannerLFQ | 83 | 91 | 115 | 21 |

Share of one-way pairs where the direction test agrees with C: ROSMAP 43%, Diverse 46%, Banner 47%, BannerLFQ 48%.

### 6.3 Supported edges by prior confidence

| prior C | prior edges | ROSMAP | Diverse | Banner | BannerLFQ |
|---|---|---|---|---|---|
| 0.2 | 148 | 39 (26%) | 81 (55%) | 22 (15%) | 32 (22%) |
| 0.3 | 235 | 57 (24%) | 116 (49%) | 27 (11%) | 40 (17%) |
| 0.4 | 129 | 31 (24%) | 65 (50%) | 21 (16%) | 26 (20%) |
| 0.5 | 88 | 26 (30%) | 44 (50%) | 19 (22%) | 19 (22%) |
| 0.6 | 56 | 15 (27%) | 38 (68%) | 8 (14%) | 13 (23%) |
| 0.7 | 23 | 2 (9%) | 8 (35%) | 3 (13%) | 5 (22%) |
| 0.8 | 7 | 2 (29%) | 2 (29%) | 0 (0%) | 0 (0%) |
| 0.9 | 2 | 0 (0%) | 0 (0%) | 1 (50%) | 0 (0%) |

### 6.4 Hubs of the supported graph (number of supported edges touching each protein)

- **ROSMAP** top 15: TMED10 32, BACE1 21, TMED2 20, TMED9 19, MFN2 11, APP 11, EIF2S1 10, LRP1 10, SORL1 9, COPA 8, SAR1A 8, SEC23A 8, SAR1B 7, SEC13 7, SEC31A 6
  - proteins with no supported edge (11): ABCA7, ADAM17, AP4B1, AP4E1, AP4M1, AP4S1, APBA2, COPZ1, PREB, SIGMAR1, SREBF2
- **Diverse** top 15: TMED10 53, TMED2 49, TMED9 40, APP 33, BACE1 33, EIF2S1 20, LRP1 18, SORL1 18, VPS35 17, APOE 16, SAR1A 16, COPA 15, MAPK14 15, SEC23A 12, SAR1B 12
  - proteins with no supported edge (1): PLG
- **Banner** top 15: APP 16, TMED2 15, TMED9 14, TMED10 14, VPS35 8, SORL1 8, MFN2 7, HSPA5 6, LRP1 5, COPA 5, SAR1B 5, BIN1 5, EIF2S1 4, VPS29 4, APOE 4
  - proteins with no supported edge (25): ABCA7, ADAM17, AP1G1, AP4B1, AP4E1, AP4M1, AP4S1, APBA2, COPZ1, GGA3, IDE, MAP2K3, MAPK14, NCSTN, PLCG2, PLG, PREB, PSEN1, RAB1B, SEC16A, SEC23B, SEC24A, SIGMAR1, SREBF2, VLDLR
- **BannerLFQ** top 15: TMED10 26, TMED2 23, APP 23, TMED9 21, APOE 13, PICALM 9, VPS35 9, SEC23A 9, MFN2 8, LRP1 8, SEC13 7, SAR1B 7, EIF2S1 7, OSBP 6, SORT1 6
  - proteins with no supported edge (36): ABCA7, ADAM10, ADAM17, AP1G1, AP4B1, AP4E1, AP4M1, AP4S1, APBA2, ARF6, BACE1, COPG2, COPZ1, ECE1, GGA2, GORASP1, IDE, MAP2K3, MAPK14, MAPKAPK2, NCSTN, PLCG2, PLG, PSEN1, RAB1A, SAR1A, SEC16A, SEC23B, SEC24A, SEC24D, SIGMAR1, SNX4, SORL1, SREBF2, SYK, VLDLR

### 6.5 Supported edges, strongest first (top 40 per dataset; full lists in `{ds}_causal_edges.csv`)

**ROSMAP** (172 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| AP2B1 | AP2A1 | 0.5000 | 0.4047 | 0.8536 | 400 | 1.53e-91 |
| TMED10 | TMED2 | 0.6000 | 0.2699 | 0.5458 | 400 | 2.20e-19 |
| TMED9 | TMED10 | 0.5000 | 0.2694 | 0.1828 | 400 | 2.91e-29 |
| PLCG2 | SYK | 0.5000 | 0.2301 | 0.3743 | 288 | 6.33e-15 |
| SEC24C | SEC23A | 0.5000 | 0.2258 | 0.5004 | 400 | 1.95e-19 |
| VPS29 | VPS35 | 0.5000 | 0.2151 | 0.7672 | 400 | 9.37e-18 |
| COPG1 | COPA | 0.5000 | 0.2075 | 0.2597 | 400 | 1.90e-16 |
| TMED2 | COPB2 | 0.6000 | 0.1895 | 0.0585 | 400 | 2.06e-09 |
| COPB1 | COPA | 0.5000 | 0.1773 | 0.1367 | 400 | 8.50e-12 |
| SEC31A | SEC13 | 0.4000 | 0.1763 | 0.4076 | 400 | 1.35e-18 |
| COPB2 | COPA | 0.5000 | 0.1761 | 0.2827 | 400 | 1.10e-11 |
| MFN2 | TMED10 | 0.5000 | 0.1609 | 0.2712 | 400 | 1.04e-09 |
| EIF2S1 | BACE1 | 0.7000 | 0.1562 | 0.3943 | 400 | 4.44e-05 |
| COPA | COPG2 | 0.6000 | 0.1529 | 0.4817 | 400 | 2.23e-06 |
| MFN2 | HSPA9 | 0.4000 | 0.1474 | 0.5589 | 400 | 9.04e-13 |
| RAB1B | COPA | 0.5000 | 0.1453 | -0.0776 | 400 | 4.49e-08 |
| TMED9 | TMED2 | 0.5000 | 0.1451 | 0.2331 | 400 | 4.53e-08 |
| VPS29 | VPS26A | 0.5000 | 0.1347 | 0.1589 | 400 | 4.97e-07 |
| PICALM | BIN1 | 0.3000 | 0.1316 | 0.7314 | 400 | 1.78e-18 |
| MFN2 | VDAC1 | 0.6000 | 0.1311 | 0.4286 | 400 | 6.51e-05 |
| GORASP1 | TMED2 | 0.4000 | 0.1304 | 0.1138 | 264 | 6.28e-07 |
| MAPK14 | MAPKAPK2 | 0.8000 | 0.1283 | 0.2256 | 352 | 0.0090 |
| GORASP1 | TMED10 | 0.4000 | 0.1233 | 0.0963 | 264 | 2.96e-06 |
| HSPA5 | EIF2S1 | 0.4000 | 0.1202 | 0.0681 | 400 | 1.30e-08 |
| TMED9 | EIF2S1 | 0.5000 | 0.1198 | 0.1374 | 400 | 1.03e-05 |
| BIN1 | BACE1 | 0.5000 | 0.1187 | 0.1163 | 400 | 1.21e-05 |
| HSPA5 | TMED10 | 0.4000 | 0.1183 | 0.1990 | 400 | 2.37e-08 |
| SORL1 | BACE1 | 0.6000 | 0.1180 | -0.3227 | 400 | 3.95e-04 |
| PLCG2 | MAPK14 | 0.5000 | 0.1178 | 0.0563 | 336 | 7.61e-05 |
| PICALM | BACE1 | 0.5000 | 0.1132 | 0.4908 | 400 | 3.36e-05 |
| SEC24C | TMED10 | 0.4000 | 0.1104 | 0.0524 | 400 | 2.38e-07 |
| SORL1 | LRP1 | 0.5000 | 0.1096 | 0.3031 | 400 | 6.18e-05 |
| TMED2 | RAB1B | 0.4000 | 0.1068 | 0.1776 | 400 | 6.28e-07 |
| RTN4 | RTN3 | 0.3000 | 0.1060 | 0.3346 | 288 | 1.11e-08 |
| RAB1A | COPA | 0.5000 | 0.1053 | 0.0019 | 400 | 1.27e-04 |
| SYK | MAPK14 | 0.4000 | 0.1012 | 0.1884 | 312 | 4.20e-05 |
| TMED9 | COPB2 | 0.4000 | 0.0992 | 0.1358 | 400 | 4.66e-06 |
| MFN2 | TMED2 | 0.3000 | 0.0981 | 0.1077 | 400 | 5.16e-10 |
| GGA3 | GGA1 | 0.3000 | 0.0963 | 0.1823 | 400 | 1.11e-09 |
| SEC24B | SEC23A | 0.4000 | 0.0943 | 0.2180 | 400 | 1.38e-05 |

**Diverse** (354 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| COPA | COPB1 | 0.6000 | 0.4520 | 0.6957 | 676 | 1.62e-122 |
| TMED10 | TMED9 | 0.6000 | 0.4273 | 0.7507 | 676 | 1.29e-103 |
| COPA | COPB2 | 0.6000 | 0.4238 | 0.6569 | 676 | 2.72e-101 |
| COPA | COPG1 | 0.6000 | 0.4200 | 0.7921 | 676 | 8.07e-99 |
| TMED10 | TMED2 | 0.6000 | 0.4023 | 0.5801 | 676 | 1.03e-87 |
| AP2B1 | AP2A1 | 0.5000 | 0.3774 | 0.5437 | 676 | 4.78e-123 |
| COPA | COPE | 0.6000 | 0.3557 | 0.7055 | 676 | 7.06e-64 |
| COPA | COPZ1 | 0.6000 | 0.3261 | 0.6051 | 676 | 9.22e-52 |
| VPS29 | VPS26A | 0.5000 | 0.3250 | 0.7610 | 676 | 9.78e-81 |
| VPS29 | VPS35 | 0.5000 | 0.3107 | 0.6277 | 676 | 7.54e-72 |
| TMED9 | TMED2 | 0.5000 | 0.2982 | 0.2264 | 676 | 8.00e-65 |
| EIF2S1 | HSPA5 | 0.6000 | 0.2871 | 0.9095 | 676 | 1.10e-38 |
| PSEN1 | NCSTN | 0.8000 | 0.2700 | 0.3594 | 676 | 1.47e-18 |
| COPA | COPG2 | 0.6000 | 0.2688 | 0.5234 | 676 | 1.76e-33 |
| SEC24C | SEC23A | 0.5000 | 0.2586 | 0.4131 | 676 | 3.44e-46 |
| SEC31A | SEC13 | 0.4000 | 0.2458 | 0.5637 | 676 | 7.76e-70 |
| LRP1 | APP | 0.7000 | 0.2353 | 2.2884 | 676 | 2.00e-18 |
| VPS26A | VPS35 | 0.5000 | 0.2281 | 0.1864 | 676 | 8.60e-35 |
| SYK | PLCG2 | 0.5000 | 0.2132 | 0.4888 | 664 | 1.55e-29 |
| RTN4 | BACE1 | 0.5000 | 0.1854 | -0.3355 | 676 | 1.72e-22 |
| TMED10 | SEC24C | 0.5000 | 0.1824 | 0.1034 | 676 | 9.98e-22 |
| SEC23A | SEC24B | 0.6000 | 0.1821 | 0.2529 | 676 | 4.30e-15 |
| COPA | TMED10 | 0.4000 | 0.1801 | 0.4398 | 676 | 8.11e-34 |
| SAR1B | SEC23B | 0.7000 | 0.1777 | 0.2101 | 676 | 8.77e-11 |
| APOE | APP | 0.4000 | 0.1717 | 0.4024 | 676 | 1.70e-30 |
| SAR1B | SEC24B | 0.6000 | 0.1702 | 0.1184 | 676 | 2.89e-13 |
| TMED2 | COPG2 | 0.6000 | 0.1683 | 0.1044 | 676 | 5.52e-13 |
| CLU | APP | 0.4000 | 0.1632 | 0.2123 | 676 | 2.25e-27 |
| TMED2 | SEC24C | 0.5000 | 0.1623 | 0.0777 | 676 | 3.34e-17 |
| CLU | APOE | 0.3000 | 0.1597 | 1.0311 | 676 | 2.66e-49 |
| SAR1A | SEC24A | 0.6000 | 0.1578 | 0.2995 | 676 | 1.64e-11 |
| BIN1 | PICALM | 0.3000 | 0.1552 | 0.2221 | 676 | 3.44e-46 |
| COPA | TMED2 | 0.4000 | 0.1551 | 0.1396 | 676 | 1.16e-24 |
| MAPKAPK2 | MAPK14 | 0.3000 | 0.1480 | 0.3861 | 667 | 6.37e-41 |
| HSPA5 | TMED10 | 0.4000 | 0.1392 | 0.1990 | 676 | 9.38e-20 |
| EIF2S1 | BACE1 | 0.7000 | 0.1365 | 0.0568 | 676 | 8.85e-07 |
| MFN2 | HSPA9 | 0.4000 | 0.1336 | 0.5823 | 676 | 3.41e-18 |
| TMED10 | COPG2 | 0.5000 | 0.1315 | -0.0423 | 676 | 1.67e-11 |
| MFN2 | VDAC1 | 0.6000 | 0.1313 | 0.3557 | 676 | 2.89e-08 |
| LRP1 | SORL1 | 0.5000 | 0.1306 | 0.3324 | 676 | 2.27e-11 |

**Banner** (101 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| TMED2 | TMED10 | 0.9000 | 0.5336 | 0.6327 | 190 | 2.31e-17 |
| LRP1 | APP | 0.7000 | 0.4212 | 4.1249 | 190 | 7.26e-18 |
| AP2B1 | AP2A1 | 0.5000 | 0.3900 | 0.6026 | 190 | 1.39e-37 |
| VPS35 | VPS29 | 0.5000 | 0.2348 | 0.4788 | 190 | 7.12e-10 |
| COPA | COPB2 | 0.6000 | 0.2152 | 0.1978 | 190 | 8.88e-06 |
| TMED9 | TMED10 | 0.5000 | 0.2138 | 0.4644 | 190 | 4.34e-08 |
| TMED9 | TMED2 | 0.5000 | 0.2118 | 0.1380 | 190 | 5.58e-08 |
| COPA | COPB1 | 0.6000 | 0.2083 | 0.2193 | 190 | 1.77e-05 |
| VPS26A | VPS29 | 0.5000 | 0.1945 | 0.2991 | 190 | 8.38e-07 |
| SEC23A | SEC24C | 0.6000 | 0.1934 | 0.3675 | 190 | 9.04e-05 |
| LRP1 | SORL1 | 0.5000 | 0.1800 | 0.6343 | 190 | 8.63e-06 |
| EIF2S1 | HSPA5 | 0.6000 | 0.1726 | 0.6247 | 190 | 6.90e-04 |
| VPS29 | APP | 0.6000 | 0.1700 | 1.5049 | 190 | 8.16e-04 |
| HSPA9 | MFN2 | 0.4000 | 0.1681 | 0.2626 | 190 | 6.71e-08 |
| APOE | APP | 0.4000 | 0.1658 | 0.7111 | 190 | 1.05e-07 |
| COPG1 | COPA | 0.5000 | 0.1579 | 0.3619 | 190 | 1.35e-04 |
| SEC23A | SEC24B | 0.6000 | 0.1433 | 0.1773 | 190 | 0.0075 |
| SEC31A | SEC13 | 0.4000 | 0.1426 | 0.4028 | 190 | 9.92e-06 |
| TMED10 | COPE | 0.5000 | 0.1423 | 0.1832 | 190 | 7.91e-04 |
| CLU | APP | 0.4000 | 0.1420 | 0.9786 | 190 | 1.05e-05 |
| TMED2 | COPB1 | 0.7000 | 0.1411 | 0.0213 | 190 | 0.0273 |
| OSBP | MFN2 | 0.3000 | 0.1389 | -0.6520 | 190 | 1.18e-09 |
| TMED10 | COPB1 | 0.6000 | 0.1377 | 0.0072 | 190 | 0.0110 |
| VPS26A | VPS35 | 0.5000 | 0.1372 | 0.2330 | 190 | 0.0014 |
| AP2B1 | APP | 0.6000 | 0.1235 | 1.8310 | 190 | 0.0240 |
| TMED2 | MIA3 | 0.3000 | 0.1232 | 0.2865 | 190 | 1.35e-07 |
| CLU | APOE | 0.3000 | 0.1180 | 0.7834 | 190 | 6.20e-07 |
| RAB1A | TMED2 | 0.4000 | 0.1179 | -0.3181 | 190 | 4.69e-04 |
| HSPA5 | TMED10 | 0.4000 | 0.1145 | 0.1778 | 190 | 7.24e-04 |
| SEC24D | GORASP1 | 0.3000 | 0.1138 | 0.5015 | 52 | 0.0279 |
| COPG2 | COPA | 0.5000 | 0.1106 | 0.1268 | 190 | 0.0146 |
| TMED10 | VDAC1 | 0.4000 | 0.1087 | 0.0771 | 190 | 0.0015 |
| TMED10 | MIA3 | 0.3000 | 0.1034 | 0.2457 | 190 | 2.10e-05 |
| SAR1A | TMED10 | 0.4000 | 0.1006 | 0.3106 | 190 | 0.0043 |
| SAR1B | TMED2 | 0.4000 | 0.1003 | 0.1077 | 190 | 0.0043 |
| AP2B1 | PICALM | 0.3000 | 0.0983 | 0.5452 | 190 | 6.60e-05 |
| BIN1 | BACE1 | 0.5000 | 0.0980 | 0.5502 | 182 | 0.0379 |
| VDAC1 | HSPA9 | 0.3000 | 0.0909 | 0.2340 | 190 | 3.00e-04 |
| COPB2 | TMED2 | 0.3000 | 0.0901 | 0.0387 | 190 | 3.40e-04 |
| OSBP | SEC24D | 0.4000 | 0.0901 | 0.6019 | 190 | 0.0131 |

**BannerLFQ** (135 supported)

| source | target | prior_C | W | beta | n_pairwise | q_fdr |
|---|---|---|---|---|---|---|
| AP2A1 | AP2B1 | 0.5000 | 0.3466 | 0.8386 | 190 | 5.33e-26 |
| TMED10 | OSBP | 0.6000 | 0.3155 | 0.3703 | 190 | 3.83e-13 |
| TMED10 | COPB1 | 0.6000 | 0.2886 | -0.2511 | 190 | 7.80e-11 |
| TMED9 | OSBP | 0.4000 | 0.2477 | 0.2918 | 190 | 3.91e-19 |
| TMED9 | TMED10 | 0.5000 | 0.2408 | 0.1259 | 190 | 7.71e-11 |
| PICALM | APP | 0.7000 | 0.2376 | 0.2318 | 190 | 2.12e-05 |
| TMED10 | COPE | 0.5000 | 0.2352 | -0.2652 | 190 | 2.59e-10 |
| VPS26A | VPS29 | 0.5000 | 0.2306 | -0.3506 | 190 | 6.12e-10 |
| TMED9 | COPB2 | 0.4000 | 0.2227 | 0.1790 | 190 | 5.56e-15 |
| VPS35 | APP | 0.7000 | 0.2077 | 0.4644 | 190 | 2.68e-04 |
| VPS35 | VPS26A | 0.5000 | 0.2039 | 0.9371 | 190 | 1.09e-07 |
| AP2A1 | APP | 0.6000 | 0.1936 | 0.7665 | 190 | 6.13e-05 |
| TMED9 | COPE | 0.4000 | 0.1935 | -0.1475 | 190 | 6.73e-11 |
| SEC13 | SEC31A | 0.7000 | 0.1911 | 0.1548 | 190 | 8.87e-04 |
| TMED10 | TMED2 | 0.6000 | 0.1835 | 0.3015 | 190 | 1.59e-04 |
| MFN2 | HSPA9 | 0.4000 | 0.1821 | 0.4257 | 190 | 1.06e-09 |
| TMED10 | SEC23A | 0.4000 | 0.1643 | 0.4984 | 190 | 8.67e-08 |
| SAR1B | SEC23A | 0.7000 | 0.1575 | 0.1679 | 174 | 0.0128 |
| TMED10 | EIF2S1 | 0.5000 | 0.1558 | 0.0803 | 190 | 1.22e-04 |
| TMED9 | TMED2 | 0.5000 | 0.1554 | 0.1874 | 190 | 1.25e-04 |
| SEC24B | SEC23A | 0.4000 | 0.1550 | -0.2849 | 164 | 4.99e-06 |
| MFN2 | TMED10 | 0.5000 | 0.1546 | 0.3477 | 190 | 1.35e-04 |
| TMED10 | COPG1 | 0.5000 | 0.1513 | 0.2866 | 190 | 1.91e-04 |
| TMED10 | SAR1B | 0.4000 | 0.1394 | 0.3489 | 174 | 2.86e-05 |
| OSBP | TMED2 | 0.4000 | 0.1388 | 0.2963 | 190 | 1.32e-05 |
| TMED10 | GGA3 | 0.3000 | 0.1369 | -0.5492 | 190 | 1.03e-09 |
| TMED10 | ITPR2 | 0.4000 | 0.1361 | -0.7529 | 125 | 7.02e-04 |
| RTN3 | APP | 0.4000 | 0.1347 | 0.3841 | 190 | 2.50e-05 |
| SEC23A | SEC24C | 0.6000 | 0.1267 | -0.1047 | 190 | 0.0151 |
| TMED10 | SEC24C | 0.5000 | 0.1252 | -0.1594 | 190 | 0.0028 |
| MIA3 | SAR1B | 0.4000 | 0.1228 | 0.1738 | 158 | 6.40e-04 |
| SEC13 | TMED10 | 0.2000 | 0.1221 | 0.2566 | 190 | 1.24e-18 |
| TMED9 | SEC13 | 0.2000 | 0.1214 | 0.3890 | 190 | 1.93e-18 |
| MFN2 | VDAC1 | 0.6000 | 0.1208 | 0.2391 | 190 | 0.0221 |
| SORT1 | GGA1 | 0.3000 | 0.1199 | 0.7416 | 67 | 0.0044 |
| CLU | APOE | 0.3000 | 0.1190 | 0.7826 | 190 | 2.85e-07 |
| AP2B1 | APP | 0.6000 | 0.1161 | 0.4243 | 190 | 0.0300 |
| TMED10 | SORT1 | 0.2000 | 0.1153 | 0.5834 | 190 | 3.16e-16 |
| SNX6 | APP | 0.5000 | 0.1151 | -0.1929 | 190 | 0.0068 |
| SNX6 | VPS35 | 0.4000 | 0.1118 | -0.0808 | 190 | 6.72e-04 |

### 6.6 Selected edges across datasets

| pair | dataset | edge | prior_C | W | beta | n | q | supported |
|---|---|---|---|---|---|---|---|---|
| TMED2–TMED10 | ROSMAP | TMED10→TMED2 | 0.6000 | 0.2699 | 0.5458 | 400 | 2.20e-19 | 1 |
| TMED2–TMED10 | ROSMAP | TMED2→TMED10 | 0.9000 | 0.0000 | 0.0000 | 400 | 2.20e-19 | 0 |
| TMED2–TMED10 | Diverse | TMED10→TMED2 | 0.6000 | 0.4023 | 0.5801 | 676 | 1.03e-87 | 1 |
| TMED2–TMED10 | Diverse | TMED2→TMED10 | 0.9000 | 0.0000 | 0.0000 | 676 | 1.03e-87 | 0 |
| TMED2–TMED10 | Banner | TMED2→TMED10 | 0.9000 | 0.5336 | 0.6327 | 190 | 2.31e-17 | 1 |
| TMED2–TMED10 | Banner | TMED10→TMED2 | 0.6000 | 0.0000 | 0.0000 | 190 | 2.31e-17 | 0 |
| TMED2–TMED10 | BannerLFQ | TMED10→TMED2 | 0.6000 | 0.1835 | 0.3015 | 190 | 1.59e-04 | 1 |
| TMED2–TMED10 | BannerLFQ | TMED2→TMED10 | 0.9000 | 0.0000 | 0.0000 | 190 | 1.59e-04 | 0 |
| TMED9–TMED2 | ROSMAP | TMED9→TMED2 | 0.5000 | 0.1451 | 0.2331 | 400 | 4.53e-08 | 1 |
| TMED9–TMED2 | ROSMAP | TMED2→TMED9 | 0.7000 | 0.0000 | 0.0000 | 400 | 4.53e-08 | 0 |
| TMED9–TMED2 | Diverse | TMED9→TMED2 | 0.5000 | 0.2982 | 0.2264 | 676 | 8.00e-65 | 1 |
| TMED9–TMED2 | Diverse | TMED2→TMED9 | 0.7000 | 0.0000 | 0.0000 | 676 | 8.00e-65 | 0 |
| TMED9–TMED2 | Banner | TMED9→TMED2 | 0.5000 | 0.2118 | 0.1380 | 190 | 5.58e-08 | 1 |
| TMED9–TMED2 | Banner | TMED2→TMED9 | 0.7000 | 0.0000 | 0.0000 | 190 | 5.58e-08 | 0 |
| TMED9–TMED2 | BannerLFQ | TMED9→TMED2 | 0.5000 | 0.1554 | 0.1874 | 190 | 1.25e-04 | 1 |
| TMED9–TMED2 | BannerLFQ | TMED2→TMED9 | 0.7000 | 0.0000 | 0.0000 | 190 | 1.25e-04 | 0 |
| TMED9–TMED10 | ROSMAP | TMED9→TMED10 | 0.5000 | 0.2694 | 0.1828 | 400 | 2.91e-29 | 1 |
| TMED9–TMED10 | ROSMAP | TMED10→TMED9 | 0.6000 | 0.0000 | 0.0000 | 400 | 2.91e-29 | 0 |
| TMED9–TMED10 | Diverse | TMED10→TMED9 | 0.6000 | 0.4273 | 0.7507 | 676 | 1.29e-103 | 1 |
| TMED9–TMED10 | Diverse | TMED9→TMED10 | 0.5000 | 0.0000 | 0.0000 | 676 | 1.29e-103 | 0 |
| TMED9–TMED10 | Banner | TMED9→TMED10 | 0.5000 | 0.2138 | 0.4644 | 190 | 4.34e-08 | 1 |
| TMED9–TMED10 | Banner | TMED10→TMED9 | 0.6000 | 0.0000 | 0.0000 | 190 | 4.34e-08 | 0 |
| TMED9–TMED10 | BannerLFQ | TMED9→TMED10 | 0.5000 | 0.2408 | 0.1259 | 190 | 7.71e-11 | 1 |
| TMED9–TMED10 | BannerLFQ | TMED10→TMED9 | 0.6000 | 0.0000 | 0.0000 | 190 | 7.71e-11 | 0 |
| APP–BACE1 | ROSMAP | BACE1→APP | 0.9000 | 0.0413 | -0.1399 | 312 | 0.5744 | 0 |
| APP–BACE1 | ROSMAP | APP→BACE1 | 0.4000 | 0.0000 | 0.0000 | 312 | 0.5744 | 0 |
| APP–BACE1 | Diverse | APP→BACE1 | 0.4000 | 0.0717 | 0.0035 | 676 | 6.78e-06 | 1 |
| APP–BACE1 | Diverse | BACE1→APP | 0.9000 | 0.0000 | 0.0000 | 676 | 6.78e-06 | 0 |
| APP–BACE1 | Banner | APP→BACE1 | 0.4000 | 0.0531 | -0.0840 | 182 | 0.1954 | 0 |
| APP–BACE1 | Banner | BACE1→APP | 0.9000 | 0.0000 | 0.0000 | 182 | 0.1954 | 0 |
| APP–BACE1 | BannerLFQ | APP→BACE1 | 0.4000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| APP–BACE1 | BannerLFQ | BACE1→APP | 0.9000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| BACE1–APP | ROSMAP | BACE1→APP | 0.9000 | 0.0413 | -0.1399 | 312 | 0.5744 | 0 |
| BACE1–APP | ROSMAP | APP→BACE1 | 0.4000 | 0.0000 | 0.0000 | 312 | 0.5744 | 0 |
| BACE1–APP | Diverse | APP→BACE1 | 0.4000 | 0.0717 | 0.0035 | 676 | 6.78e-06 | 1 |
| BACE1–APP | Diverse | BACE1→APP | 0.9000 | 0.0000 | 0.0000 | 676 | 6.78e-06 | 0 |
| BACE1–APP | Banner | APP→BACE1 | 0.4000 | 0.0531 | -0.0840 | 182 | 0.1954 | 0 |
| BACE1–APP | Banner | BACE1→APP | 0.9000 | 0.0000 | 0.0000 | 182 | 0.1954 | 0 |
| BACE1–APP | BannerLFQ | APP→BACE1 | 0.4000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| BACE1–APP | BannerLFQ | BACE1→APP | 0.9000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| LRP1–APP | ROSMAP | APP→LRP1 | 0.4000 | 0.0333 | 0.0024 | 312 | 0.2470 | 0 |
| LRP1–APP | ROSMAP | LRP1→APP | 0.7000 | 0.0000 | 0.0000 | 312 | 0.2470 | 0 |
| LRP1–APP | Diverse | LRP1→APP | 0.7000 | 0.2353 | 2.2884 | 676 | 2.00e-18 | 1 |
| LRP1–APP | Diverse | APP→LRP1 | 0.4000 | 0.0000 | 0.0000 | 676 | 2.00e-18 | 0 |
| LRP1–APP | Banner | LRP1→APP | 0.7000 | 0.4212 | 4.1249 | 190 | 7.26e-18 | 1 |
| LRP1–APP | Banner | APP→LRP1 | 0.4000 | 0.0000 | 0.0000 | 190 | 7.26e-18 | 0 |
| LRP1–APP | BannerLFQ | APP→LRP1 | 0.4000 | 0.0610 | -0.0983 | 190 | 0.1230 | 0 |
| LRP1–APP | BannerLFQ | LRP1→APP | 0.7000 | 0.0000 | 0.0000 | 190 | 0.1230 | 0 |
| APOE–APP | ROSMAP | APOE→APP | 0.4000 | 0.0203 | 0.0318 | 312 | 0.5382 | 0 |
| APOE–APP | Diverse | APOE→APP | 0.4000 | 0.1717 | 0.4024 | 676 | 1.70e-30 | 1 |
| APOE–APP | Banner | APOE→APP | 0.4000 | 0.1658 | 0.7111 | 190 | 1.05e-07 | 1 |
| APOE–APP | BannerLFQ | APOE→APP | 0.4000 | 0.0942 | 0.1337 | 190 | 0.0054 | 1 |
| CLU–APOE | ROSMAP | CLU→APOE | 0.3000 | 0.0855 | 0.7426 | 400 | 8.29e-08 | 1 |
| CLU–APOE | ROSMAP | APOE→CLU | 0.4000 | 0.0000 | 0.0000 | 400 | 8.29e-08 | 0 |
| CLU–APOE | Diverse | CLU→APOE | 0.3000 | 0.1597 | 1.0311 | 676 | 2.66e-49 | 1 |
| CLU–APOE | Diverse | APOE→CLU | 0.4000 | 0.0000 | 0.0000 | 676 | 2.66e-49 | 0 |
| CLU–APOE | Banner | CLU→APOE | 0.3000 | 0.1180 | 0.7834 | 190 | 6.20e-07 | 1 |
| CLU–APOE | Banner | APOE→CLU | 0.4000 | 0.0000 | 0.0000 | 190 | 6.20e-07 | 0 |
| CLU–APOE | BannerLFQ | CLU→APOE | 0.3000 | 0.1190 | 0.7826 | 190 | 2.85e-07 | 1 |
| CLU–APOE | BannerLFQ | APOE→CLU | 0.4000 | 0.0000 | 0.0000 | 190 | 2.85e-07 | 0 |
| SORL1–APP | ROSMAP | APP→SORL1 | 0.3000 | 0.0822 | 0.0472 | 312 | 7.45e-06 | 1 |
| SORL1–APP | ROSMAP | SORL1→APP | 0.8000 | 0.0000 | 0.0000 | 312 | 7.45e-06 | 0 |
| SORL1–APP | Diverse | SORL1→APP | 0.8000 | 0.0010 | -0.4243 | 676 | 0.9795 | 0 |
| SORL1–APP | Diverse | APP→SORL1 | 0.3000 | 0.0000 | 0.0000 | 676 | 0.9795 | 0 |
| SORL1–APP | Banner | SORL1→APP | 0.8000 | 0.0485 | 0.2337 | 190 | 0.5950 | 0 |
| SORL1–APP | Banner | APP→SORL1 | 0.3000 | 0.0000 | 0.0000 | 190 | 0.5950 | 0 |
| SORL1–APP | BannerLFQ | APP→SORL1 | 0.3000 | 0.0401 | 0.1974 | 170 | 0.2513 | 0 |
| SORL1–APP | BannerLFQ | SORL1→APP | 0.8000 | 0.0000 | 0.0000 | 170 | 0.2513 | 0 |
| PSEN1–NCSTN | ROSMAP | NCSTN→PSEN1 | 0.6000 | 0.0936 | -0.0006 | 384 | 0.0077 | 1 |
| PSEN1–NCSTN | ROSMAP | PSEN1→NCSTN | 0.8000 | 0.0000 | 0.0000 | 384 | 0.0077 | 0 |
| PSEN1–NCSTN | Diverse | PSEN1→NCSTN | 0.8000 | 0.2700 | 0.3594 | 676 | 1.47e-18 | 1 |
| PSEN1–NCSTN | Diverse | NCSTN→PSEN1 | 0.6000 | 0.0000 | 0.0000 | 676 | 1.47e-18 | 0 |
| PSEN1–NCSTN | Banner | NCSTN→PSEN1 | 0.6000 | 0.0907 | 0.0409 | 190 | 0.1228 | 0 |
| PSEN1–NCSTN | Banner | PSEN1→NCSTN | 0.8000 | 0.0000 | 0.0000 | 190 | 0.1228 | 0 |
| PSEN1–NCSTN | BannerLFQ | NCSTN→PSEN1 | 0.6000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| PSEN1–NCSTN | BannerLFQ | PSEN1→NCSTN | 0.8000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| BIN1–BACE1 | ROSMAP | BIN1→BACE1 | 0.5000 | 0.1187 | 0.1163 | 400 | 1.21e-05 | 1 |
| BIN1–BACE1 | Diverse | BIN1→BACE1 | 0.5000 | 0.0547 | 0.0412 | 676 | 0.0072 | 1 |
| BIN1–BACE1 | Banner | BIN1→BACE1 | 0.5000 | 0.0980 | 0.5502 | 182 | 0.0379 | 1 |
| BIN1–BACE1 | BannerLFQ | BIN1→BACE1 | 0.5000 | 0.0000 | 0.0000 | 0 | 1.0000 | 0 |
| PICALM–APP | ROSMAP | PICALM→APP | 0.7000 | 0.0200 | 0.3949 | 312 | 0.1759 | 0 |
| PICALM–APP | Diverse | PICALM→APP | 0.7000 | 0.0801 | 0.3064 | 676 | 0.0049 | 1 |
| PICALM–APP | Banner | PICALM→APP | 0.7000 | 0.0468 | 0.3640 | 190 | 0.5530 | 0 |
| PICALM–APP | BannerLFQ | PICALM→APP | 0.7000 | 0.2376 | 0.2318 | 190 | 2.12e-05 | 1 |

---

## 7. Agreement between the three graphs (HIW 5)

### 7.1 Pairwise overlap of supported pairs vs chance

| pair | jaccard | shared | expected_by_chance | fold | p_hypergeom | same_direction |
|---|---|---|---|---|---|---|
| ROSMAP-Diverse | 0.3734 | 143 | 120.0947 | 1.1907 | 1.22e-06 | 0.7413 |
| ROSMAP-Banner | 0.2523 | 55 | 34.2643 | 1.6052 | 1.64e-06 | 0.5091 |
| Diverse-Banner | 0.2264 | 84 | 70.5207 | 1.1911 | 5.52e-04 | 0.7500 |

### 7.2 Consensus pairs

- pairs supported in ≥ 2 datasets: **178** (in exactly 2: 126; in all 3: **52**)
- in exactly 2: ROSMAP+Diverse 91, Diverse+Banner 32, ROSMAP+Banner 3
- arrow identical in every dataset that supports the pair: 124 of 178; among the 52 triple pairs: 21 of 52
- split by prior type — one-way C pairs (arrow fixed by C): 93 of 178 consensus pairs, all 93 consistent; two-way C pairs (arrow picked by the direction test): 85, of which 31 consistent and 54 not
- among the 52 triple pairs: 15 one-way (all consistent); 37 two-way, of which 6 consistent (AP2A1–AP2B1, TMED2–TMED9, SEC13–SEC31A, HSPA5–TMED10, APOE–CLU, HSPA5–TMED2) and 31 not

**All 52 pairs supported in all three datasets**

| protein_a | protein_b | prior_C_a_to_b | prior_C_b_to_a | dir_ROSMAP | W_ROSMAP | dir_Diverse | W_Diverse | dir_Banner | W_Banner | direction_consistent |
|---|---|---|---|---|---|---|---|---|---|---|
| TMED10 | TMED2 | 0.600 | 0.900 | TMED10->TMED2 | 0.270 | TMED10->TMED2 | 0.402 | TMED2->TMED10 | 0.534 | 0 |
| AP2A1 | AP2B1 | 0.500 | 0.500 | AP2B1->AP2A1 | 0.405 | AP2B1->AP2A1 | 0.377 | AP2B1->AP2A1 | 0.390 | 1 |
| TMED10 | TMED9 | 0.600 | 0.500 | TMED9->TMED10 | 0.269 | TMED10->TMED9 | 0.427 | TMED9->TMED10 | 0.214 | 0 |
| COPA | COPB1 | 0.600 | 0.500 | COPB1->COPA | 0.177 | COPA->COPB1 | 0.452 | COPA->COPB1 | 0.208 | 0 |
| COPA | COPB2 | 0.600 | 0.500 | COPB2->COPA | 0.176 | COPA->COPB2 | 0.424 | COPA->COPB2 | 0.215 | 0 |
| COPA | COPG1 | 0.600 | 0.500 | COPG1->COPA | 0.207 | COPA->COPG1 | 0.420 | COPG1->COPA | 0.158 | 0 |
| VPS29 | VPS35 | 0.500 | 0.500 | VPS29->VPS35 | 0.215 | VPS29->VPS35 | 0.311 | VPS35->VPS29 | 0.235 | 0 |
| SEC23A | SEC24C | 0.600 | 0.500 | SEC24C->SEC23A | 0.226 | SEC24C->SEC23A | 0.259 | SEC23A->SEC24C | 0.193 | 0 |
| TMED2 | TMED9 | 0.700 | 0.500 | TMED9->TMED2 | 0.145 | TMED9->TMED2 | 0.298 | TMED9->TMED2 | 0.212 | 1 |
| VPS26A | VPS29 | 0.500 | 0.500 | VPS29->VPS26A | 0.135 | VPS29->VPS26A | 0.325 | VPS26A->VPS29 | 0.195 | 0 |
| EIF2S1 | HSPA5 | 0.600 | 0.400 | HSPA5->EIF2S1 | 0.120 | EIF2S1->HSPA5 | 0.287 | EIF2S1->HSPA5 | 0.173 | 0 |
| SEC13 | SEC31A | 0.700 | 0.400 | SEC31A->SEC13 | 0.176 | SEC31A->SEC13 | 0.246 | SEC31A->SEC13 | 0.143 | 1 |
| COPA | COPG2 | 0.600 | 0.500 | COPA->COPG2 | 0.153 | COPA->COPG2 | 0.269 | COPG2->COPA | 0.111 | 0 |
| HSPA9 | MFN2 | 0.400 | 0.400 | MFN2->HSPA9 | 0.147 | MFN2->HSPA9 | 0.134 | HSPA9->MFN2 | 0.168 | 0 |
| LRP1 | SORL1 | 0.500 | 0.500 | SORL1->LRP1 | 0.110 | LRP1->SORL1 | 0.131 | LRP1->SORL1 | 0.180 | 0 |
| SEC23A | SEC24B | 0.600 | 0.400 | SEC24B->SEC23A | 0.094 | SEC23A->SEC24B | 0.182 | SEC23A->SEC24B | 0.143 | 0 |
| COPB2 | TMED2 | 0.300 | 0.600 | TMED2->COPB2 | 0.190 | COPB2->TMED2 | 0.122 | COPB2->TMED2 | 0.090 | 0 |
| SEC24C | TMED10 | 0.400 | 0.500 | SEC24C->TMED10 | 0.110 | TMED10->SEC24C | 0.182 | SEC24C->TMED10 | 0.084 | 0 |
| HSPA5 | TMED10 | 0.400 | 0.400 | HSPA5->TMED10 | 0.118 | HSPA5->TMED10 | 0.139 | HSPA5->TMED10 | 0.115 | 1 |
| APOE | CLU | 0.400 | 0.300 | CLU->APOE | 0.086 | CLU->APOE | 0.160 | CLU->APOE | 0.118 | 1 |
| APP | CLU | 0.200 | 0.400 | APP->CLU | 0.032 | CLU->APP | 0.163 | CLU->APP | 0.142 | 0 |
| MFN2 | TMED10 | 0.500 | 0.300 | MFN2->TMED10 | 0.161 | MFN2->TMED10 | 0.098 | TMED10->MFN2 | 0.056 | 0 |
| SAR1A | TMED10 | 0.400 | 0.400 | TMED10->SAR1A | 0.082 | SAR1A->TMED10 | 0.129 | SAR1A->TMED10 | 0.101 | 0 |
| RTN3 | RTN4 | 0.300 | 0.300 | RTN4->RTN3 | 0.106 | RTN4->RTN3 | 0.124 | RTN3->RTN4 | 0.077 | 0 |
| COPB1 | TMED2 | 0.300 | 0.700 | COPB1->TMED2 | 0.050 | COPB1->TMED2 | 0.111 | TMED2->COPB1 | 0.141 | 0 |
| BACE1 | BIN1 | 0.000 | 0.500 | BIN1->BACE1 | 0.119 | BIN1->BACE1 | 0.055 | BIN1->BACE1 | 0.098 | 1 |
| COPB1 | TMED10 | 0.000 | 0.600 | TMED10->COPB1 | 0.056 | TMED10->COPB1 | 0.077 | TMED10->COPB1 | 0.138 | 1 |
| SAR1A | TMED2 | 0.400 | 0.400 | TMED2->SAR1A | 0.055 | SAR1A->TMED2 | 0.126 | SAR1A->TMED2 | 0.084 | 0 |
| SAR1B | TMED2 | 0.400 | 0.400 | TMED2->SAR1B | 0.048 | TMED2->SAR1B | 0.085 | SAR1B->TMED2 | 0.100 | 0 |
| COPE | TMED10 | 0.000 | 0.500 | TMED10->COPE | 0.027 | TMED10->COPE | 0.060 | TMED10->COPE | 0.142 | 1 |
| COPB2 | TMED9 | 0.000 | 0.400 | TMED9->COPB2 | 0.099 | TMED9->COPB2 | 0.065 | TMED9->COPB2 | 0.050 | 1 |
| SAR1B | TMED10 | 0.400 | 0.400 | SAR1B->TMED10 | 0.058 | TMED10->SAR1B | 0.064 | TMED10->SAR1B | 0.085 | 0 |
| MFN2 | TMED2 | 0.300 | 0.300 | MFN2->TMED2 | 0.098 | MFN2->TMED2 | 0.044 | TMED2->MFN2 | 0.063 | 0 |
| SNX6 | VPS35 | 0.400 | 0.400 | VPS35->SNX6 | 0.069 | SNX6->VPS35 | 0.040 | SNX6->VPS35 | 0.076 | 0 |
| HSPA9 | VDAC1 | 0.300 | 0.300 | HSPA9->VDAC1 | 0.038 | VDAC1->HSPA9 | 0.056 | VDAC1->HSPA9 | 0.091 | 0 |
| COPB2 | TMED10 | 0.000 | 0.500 | TMED10->COPB2 | 0.053 | TMED10->COPB2 | 0.078 | TMED10->COPB2 | 0.050 | 1 |
| HSPA5 | TMED2 | 0.300 | 0.400 | HSPA5->TMED2 | 0.042 | HSPA5->TMED2 | 0.066 | HSPA5->TMED2 | 0.071 | 1 |
| SNX4 | SNX6 | 0.300 | 0.300 | SNX4->SNX6 | 0.039 | SNX4->SNX6 | 0.073 | SNX6->SNX4 | 0.065 | 0 |
| GGA1 | VPS35 | 0.300 | 0.300 | VPS35->GGA1 | 0.059 | VPS35->GGA1 | 0.031 | GGA1->VPS35 | 0.069 | 0 |
| MIA3 | TMED10 | 0.000 | 0.300 | TMED10->MIA3 | 0.018 | TMED10->MIA3 | 0.029 | TMED10->MIA3 | 0.103 | 1 |
| SORL1 | TMED10 | 0.200 | 0.300 | SORL1->TMED10 | 0.048 | TMED10->SORL1 | 0.049 | SORL1->TMED10 | 0.047 | 0 |
| HSPA5 | TMED9 | 0.000 | 0.400 | TMED9->HSPA5 | 0.038 | TMED9->HSPA5 | 0.054 | TMED9->HSPA5 | 0.042 | 1 |
| AP3B2 | APP | 0.400 | 0.000 | AP3B2->APP | 0.034 | AP3B2->APP | 0.016 | AP3B2->APP | 0.083 | 1 |
| TMED10 | VPS35 | 0.200 | 0.200 | VPS35->TMED10 | 0.048 | VPS35->TMED10 | 0.027 | TMED10->VPS35 | 0.042 | 0 |
| HSPA5 | MFN2 | 0.200 | 0.300 | MFN2->HSPA5 | 0.040 | HSPA5->MFN2 | 0.031 | HSPA5->MFN2 | 0.040 | 0 |
| COPG1 | TMED9 | 0.000 | 0.500 | TMED9->COPG1 | 0.020 | TMED9->COPG1 | 0.041 | TMED9->COPG1 | 0.049 | 1 |
| BACE1 | MFN2 | 0.000 | 0.200 | MFN2->BACE1 | 0.026 | MFN2->BACE1 | 0.016 | MFN2->BACE1 | 0.059 | 1 |
| SEC13 | SEC23A | 0.000 | 0.300 | SEC23A->SEC13 | 0.014 | SEC23A->SEC13 | 0.027 | SEC23A->SEC13 | 0.058 | 1 |
| ECE1 | LRP1 | 0.200 | 0.000 | ECE1->LRP1 | 0.009 | ECE1->LRP1 | 0.024 | ECE1->LRP1 | 0.048 | 1 |
| EIF2S1 | LRP1 | 0.200 | 0.000 | EIF2S1->LRP1 | 0.027 | EIF2S1->LRP1 | 0.039 | EIF2S1->LRP1 | 0.013 | 1 |
| MFN2 | TMED9 | 0.000 | 0.200 | TMED9->MFN2 | 0.009 | TMED9->MFN2 | 0.025 | TMED9->MFN2 | 0.018 | 1 |
| MIA3 | TMED9 | 0.000 | 0.200 | TMED9->MIA3 | 0.014 | TMED9->MIA3 | 0.022 | TMED9->MIA3 | 0.013 | 1 |

**Pairs supported in exactly two datasets** (126)

| protein_a | protein_b | dir_ROSMAP | W_ROSMAP | dir_Diverse | W_Diverse | dir_Banner | W_Banner |
|---|---|---|---|---|---|---|---|
| APP | LRP1 | — | — | LRP1->APP | 0.235 | LRP1->APP | 0.421 |
| COPA | COPE | COPA->COPE | 0.088 | COPA->COPE | 0.356 | — | — |
| PLCG2 | SYK | PLCG2->SYK | 0.230 | SYK->PLCG2 | 0.213 | — | — |
| VPS26A | VPS35 | — | — | VPS26A->VPS35 | 0.228 | VPS26A->VPS35 | 0.137 |
| NCSTN | PSEN1 | NCSTN->PSEN1 | 0.094 | PSEN1->NCSTN | 0.270 | — | — |
| APOE | APP | — | — | APOE->APP | 0.172 | APOE->APP | 0.166 |
| BACE1 | EIF2S1 | EIF2S1->BACE1 | 0.156 | EIF2S1->BACE1 | 0.137 | — | — |
| BIN1 | PICALM | PICALM->BIN1 | 0.132 | BIN1->PICALM | 0.155 | — | — |
| APP | VPS29 | — | — | VPS29->APP | 0.107 | VPS29->APP | 0.170 |
| MAPK14 | MAPKAPK2 | MAPK14->MAPKAPK2 | 0.128 | MAPKAPK2->MAPK14 | 0.148 | — | — |
| MFN2 | VDAC1 | MFN2->VDAC1 | 0.131 | MFN2->VDAC1 | 0.131 | — | — |
| COPA | TMED10 | COPA->TMED10 | 0.069 | COPA->TMED10 | 0.180 | — | — |
| BACE1 | PICALM | PICALM->BACE1 | 0.113 | PICALM->BACE1 | 0.119 | — | — |
| GORASP1 | TMED2 | GORASP1->TMED2 | 0.130 | GORASP1->TMED2 | 0.094 | — | — |
| SAR1B | SEC23B | SAR1B->SEC23B | 0.042 | SAR1B->SEC23B | 0.178 | — | — |
| AP2A1 | BIN1 | — | — | BIN1->AP2A1 | 0.117 | BIN1->AP2A1 | 0.085 |
| MAPK14 | SYK | SYK->MAPK14 | 0.101 | SYK->MAPK14 | 0.100 | — | — |
| GORASP1 | TMED10 | GORASP1->TMED10 | 0.123 | GORASP1->TMED10 | 0.078 | — | — |
| MIA3 | TMED2 | — | — | MIA3->TMED2 | 0.074 | TMED2->MIA3 | 0.123 |
| COPA | RAB1B | RAB1B->COPA | 0.145 | RAB1B->COPA | 0.050 | — | — |
| SEC23A | TMED10 | SEC23A->TMED10 | 0.073 | TMED10->SEC23A | 0.118 | — | — |
| OSBP | SEC24D | — | — | OSBP->SEC24D | 0.101 | OSBP->SEC24D | 0.090 |
| SEC23A | TMED2 | SEC23A->TMED2 | 0.058 | TMED2->SEC23A | 0.126 | — | — |
| COPA | RAB1A | RAB1A->COPA | 0.105 | RAB1A->COPA | 0.076 | — | — |
| MFN2 | OSBP | — | — | MFN2->OSBP | 0.042 | OSBP->MFN2 | 0.139 |
| EIF2S1 | TMED9 | TMED9->EIF2S1 | 0.120 | TMED9->EIF2S1 | 0.061 | — | — |
| TMED10 | VDAC1 | — | — | VDAC1->TMED10 | 0.067 | TMED10->VDAC1 | 0.109 |
| SAR1A | SEC24C | SEC24C->SAR1A | 0.054 | SAR1A->SEC24C | 0.118 | — | — |
| SEC23A | SEC24A | SEC23A->SEC24A | 0.083 | SEC24A->SEC23A | 0.089 | — | — |
| RAB1B | TMED2 | TMED2->RAB1B | 0.107 | RAB1B->TMED2 | 0.064 | — | — |
| SEC16A | TMED10 | SEC16A->TMED10 | 0.067 | SEC16A->TMED10 | 0.100 | — | — |
| AP2A1 | APBA1 | APBA1->AP2A1 | 0.092 | APBA1->AP2A1 | 0.074 | — | — |
| PSEN1 | TMED10 | PSEN1->TMED10 | 0.092 | TMED10->PSEN1 | 0.073 | — | — |
| EIF2S1 | TMED10 | EIF2S1->TMED10 | 0.084 | EIF2S1->TMED10 | 0.082 | — | — |
| OSBP | TMED10 | OSBP->TMED10 | 0.070 | TMED10->OSBP | 0.093 | — | — |
| MAPK14 | PLCG2 | PLCG2->MAPK14 | 0.118 | PLCG2->MAPK14 | 0.040 | — | — |
| GGA1 | GGA3 | GGA3->GGA1 | 0.096 | GGA1->GGA3 | 0.060 | — | — |
| SAR1A | SEC24B | SAR1A->SEC24B | 0.072 | SAR1A->SEC24B | 0.084 | — | — |
| BACE1 | SORL1 | SORL1->BACE1 | 0.118 | BACE1->SORL1 | 0.035 | — | — |
| COPG2 | TMED10 | TMED10->COPG2 | 0.018 | TMED10->COPG2 | 0.132 | — | — |
| AP2B1 | PICALM | — | — | AP2B1->PICALM | 0.047 | AP2B1->PICALM | 0.098 |
| APP | VPS35 | — | — | VPS35->APP | 0.096 | VPS35->APP | 0.049 |
| SAR1A | SEC31A | SAR1A->SEC31A | 0.094 | SAR1A->SEC31A | 0.050 | — | — |
| LRP1 | PICALM | LRP1->PICALM | 0.055 | PICALM->LRP1 | 0.089 | — | — |
| BACE1 | LRP1 | LRP1->BACE1 | 0.091 | BACE1->LRP1 | 0.053 | — | — |
| PICALM | SORL1 | PICALM->SORL1 | 0.077 | SORL1->PICALM | 0.067 | — | — |
| APP | EIF2S1 | EIF2S1->APP | 0.090 | EIF2S1->APP | 0.053 | — | — |
| APP | RTN4 | RTN4->APP | 0.064 | RTN4->APP | 0.079 | — | — |
| SAR1A | SEC13 | SAR1A->SEC13 | 0.052 | SAR1A->SEC13 | 0.083 | — | — |
| BACE1 | MAP2K3 | MAP2K3->BACE1 | 0.070 | MAP2K3->BACE1 | 0.060 | — | — |
| ARF6 | RAB1B | ARF6->RAB1B | 0.092 | ARF6->RAB1B | 0.034 | — | — |
| SEC13 | TMED10 | TMED10->SEC13 | 0.059 | SEC13->TMED10 | 0.066 | — | — |
| ITPR1 | ITPR2 | ITPR2->ITPR1 | 0.060 | ITPR2->ITPR1 | 0.065 | — | — |
| BACE1 | TMED2 | TMED2->BACE1 | 0.075 | TMED2->BACE1 | 0.049 | — | — |
| COPG1 | TMED2 | TMED2->COPG1 | 0.079 | TMED2->COPG1 | 0.042 | — | — |
| SEC24B | TMED2 | SEC24B->TMED2 | 0.038 | TMED2->SEC24B | 0.077 | — | — |
| HSPA9 | TMED10 | HSPA9->TMED10 | 0.068 | HSPA9->TMED10 | 0.046 | — | — |
| COPA | TMED9 | — | — | TMED9->COPA | 0.062 | TMED9->COPA | 0.051 |
| APBA1 | SORL1 | — | — | APBA1->SORL1 | 0.036 | APBA1->SORL1 | 0.075 |
| COPB1 | TMED9 | — | — | TMED9->COPB1 | 0.056 | TMED9->COPB1 | 0.051 |
| SEC24D | SEC31A | SEC24D->SEC31A | 0.060 | SEC24D->SEC31A | 0.046 | — | — |
| BACE1 | CLU | BACE1->CLU | 0.043 | CLU->BACE1 | 0.061 | — | — |
| SAR1B | SEC13 | SAR1B->SEC13 | 0.046 | SAR1B->SEC13 | 0.057 | — | — |
| APP | SAR1B | SAR1B->APP | 0.019 | — | — | SAR1B->APP | 0.084 |
| EIF2S1 | SORL1 | EIF2S1->SORL1 | 0.052 | — | — | EIF2S1->SORL1 | 0.049 |
| MIA3 | SAR1B | — | — | MIA3->SAR1B | 0.071 | MIA3->SAR1B | 0.025 |
| MAP2K3 | MAPK14 | MAPK14->MAP2K3 | 0.059 | MAPK14->MAP2K3 | 0.035 | — | — |
| GGA2 | GGA3 | GGA2->GGA3 | 0.056 | GGA3->GGA2 | 0.035 | — | — |
| COPE | TMED2 | TMED2->COPE | 0.033 | TMED2->COPE | 0.057 | — | — |
| EIF2S1 | SYK | — | — | SYK->EIF2S1 | 0.046 | SYK->EIF2S1 | 0.044 |
| BACE1 | SNX4 | SNX4->BACE1 | 0.041 | SNX4->BACE1 | 0.047 | — | — |
| APP | SNX6 | — | — | SNX6->APP | 0.055 | SNX6->APP | 0.032 |
| BACE1 | TMED10 | TMED10->BACE1 | 0.044 | TMED10->BACE1 | 0.042 | — | — |
| COPG1 | TMED10 | TMED10->COPG1 | 0.038 | TMED10->COPG1 | 0.047 | — | — |
| SEC24B | TMED10 | TMED10->SEC24B | 0.012 | TMED10->SEC24B | 0.074 | — | — |
| BACE1 | ECE1 | ECE1->BACE1 | 0.043 | ECE1->BACE1 | 0.042 | — | — |
| APP | GGA1 | — | — | GGA1->APP | 0.054 | GGA1->APP | 0.029 |
| SEC31A | TMED10 | SEC31A->TMED10 | 0.048 | SEC31A->TMED10 | 0.034 | — | — |
| COPE | TMED9 | — | — | TMED9->COPE | 0.043 | TMED9->COPE | 0.036 |
| AP3B2 | BACE1 | AP3B2->BACE1 | 0.050 | AP3B2->BACE1 | 0.027 | — | — |
| SORL1 | TMED2 | — | — | TMED2->SORL1 | 0.054 | TMED2->SORL1 | 0.020 |
| SAR1B | TMED9 | — | — | TMED9->SAR1B | 0.052 | TMED9->SAR1B | 0.020 |
| MAP2K3 | MAPKAPK2 | MAP2K3->MAPKAPK2 | 0.052 | MAP2K3->MAPKAPK2 | 0.019 | — | — |
| APP | MFN2 | MFN2->APP | 0.043 | MFN2->APP | 0.026 | — | — |
| OSBP | TMED9 | TMED9->OSBP | 0.020 | TMED9->OSBP | 0.047 | — | — |
| MAP2K3 | SYK | SYK->MAP2K3 | 0.058 | SYK->MAP2K3 | 0.009 | — | — |
| EIF2S1 | TMED2 | TMED2->EIF2S1 | 0.024 | TMED2->EIF2S1 | 0.042 | — | — |
| APP | COPB1 | COPB1->APP | 0.044 | COPB1->APP | 0.022 | — | — |
| GGA1 | TMED9 | TMED9->GGA1 | 0.047 | TMED9->GGA1 | 0.018 | — | — |
| SEC16A | SEC23A | SEC16A->SEC23A | 0.020 | SEC16A->SEC23A | 0.044 | — | — |
| LRP1 | TMED9 | TMED9->LRP1 | 0.009 | TMED9->LRP1 | 0.051 | — | — |
| SEC13 | TMED9 | TMED9->SEC13 | 0.044 | TMED9->SEC13 | 0.017 | — | — |
| SEC13 | SEC23B | SEC23B->SEC13 | 0.037 | SEC23B->SEC13 | 0.024 | — | — |
| SAR1A | TMED9 | TMED9->SAR1A | 0.038 | TMED9->SAR1A | 0.022 | — | — |
| BIN1 | LRP1 | BIN1->LRP1 | 0.028 | BIN1->LRP1 | 0.032 | — | — |
| GGA2 | SORL1 | — | — | SORL1->GGA2 | 0.038 | SORL1->GGA2 | 0.018 |
| APP | ITPR2 | — | — | ITPR2->APP | 0.040 | ITPR2->APP | 0.013 |
| ADAM10 | SORL1 | ADAM10->SORL1 | 0.042 | ADAM10->SORL1 | 0.009 | — | — |
| TMED9 | VPS35 | TMED9->VPS35 | 0.033 | TMED9->VPS35 | 0.019 | — | — |
| BACE1 | VPS26A | VPS26A->BACE1 | 0.015 | VPS26A->BACE1 | 0.036 | — | — |
| BACE1 | IDE | IDE->BACE1 | 0.031 | IDE->BACE1 | 0.020 | — | — |
| APBA1 | TMED2 | — | — | TMED2->APBA1 | 0.005 | TMED2->APBA1 | 0.044 |
| SORL1 | VPS29 | VPS29->SORL1 | 0.037 | VPS29->SORL1 | 0.012 | — | — |
| SEC23A | SEC31A | SEC23A->SEC31A | 0.025 | SEC23A->SEC31A | 0.023 | — | — |
| GGA1 | TMED10 | TMED10->GGA1 | 0.016 | TMED10->GGA1 | 0.025 | — | — |
| BACE1 | HSPA5 | HSPA5->BACE1 | 0.013 | HSPA5->BACE1 | 0.027 | — | — |
| TMED2 | VDAC1 | — | — | TMED2->VDAC1 | 0.016 | TMED2->VDAC1 | 0.024 |
| BACE1 | SAR1A | SAR1A->BACE1 | 0.013 | SAR1A->BACE1 | 0.024 | — | — |
| BIN1 | SORL1 | BIN1->SORL1 | 0.014 | — | — | BIN1->SORL1 | 0.022 |
| SEC24A | TMED10 | TMED10->SEC24A | 0.011 | TMED10->SEC24A | 0.026 | — | — |
| SORL1 | TMED9 | — | — | TMED9->SORL1 | 0.025 | TMED9->SORL1 | 0.012 |
| APOE | HSPA5 | — | — | APOE->HSPA5 | 0.017 | APOE->HSPA5 | 0.019 |
| APP | OSBP | — | — | OSBP->APP | 0.022 | OSBP->APP | 0.013 |
| CLU | VPS35 | — | — | CLU->VPS35 | 0.021 | CLU->VPS35 | 0.012 |
| GGA3 | TMED9 | TMED9->GGA3 | 0.025 | TMED9->GGA3 | 0.005 | — | — |
| IDE | TMED10 | TMED10->IDE | 0.015 | TMED10->IDE | 0.014 | — | — |
| SNX6 | VPS29 | — | — | VPS29->SNX6 | 0.011 | VPS29->SNX6 | 0.017 |
| BIN1 | TMED9 | — | — | TMED9->BIN1 | 0.013 | TMED9->BIN1 | 0.013 |
| ARF6 | TMED9 | TMED9->ARF6 | 0.015 | TMED9->ARF6 | 0.012 | — | — |
| LRP1 | VLDLR | VLDLR->LRP1 | 0.012 | VLDLR->LRP1 | 0.013 | — | — |
| APOE | VPS35 | — | — | APOE->VPS35 | 0.007 | APOE->VPS35 | 0.016 |
| PLCG2 | TMED2 | TMED2->PLCG2 | 0.009 | TMED2->PLCG2 | 0.013 | — | — |
| BIN1 | TMED10 | TMED10->BIN1 | 0.010 | TMED10->BIN1 | 0.007 | — | — |
| SEC31A | TMED9 | TMED9->SEC31A | 0.009 | TMED9->SEC31A | 0.008 | — | — |
| SYK | TMED2 | TMED2->SYK | 0.008 | TMED2->SYK | 0.009 | — | — |
| BACE1 | SYK | SYK->BACE1 | 0.009 | SYK->BACE1 | 0.007 | — | — |

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
| signed AUROC raw → adjusted | 0.5313 → 0.5811 | 0.5372 → 0.5896 |
| Z* | proteomePC3, proteomePC1, proteomePC5, proteomePC8, celltype_oligodendrocyte | proteomePC1, Braak, TangleTotal, PlaqueTotal, proteomePC9 |
| supported edges | 101 | 135 |
| STRONG regions | 0 | 0 |

- pairs testable in both: 289
- supported: TMT 91, LFQ 135, both 52; expected by chance 42.5; hypergeometric p = 0.011
- Spearman correlation of |W| across the 289 pairs: 0.198

### 8.3 Possible overlap with Diverse

- Diverse contains 43 donors from the Banner cohort; 17 of them are among the 964 Diverse people analysed.
- No ID links the two systems. Matching on sex + APOE + age at death + Braak gave 8 Banner donors
  (09-17, 11-14, 01-16, 06-09, 05-57, 05-35, 13-46, 13-66), matched to 7 Diverse donors; all 8 were removed from Banner and BannerLFQ.

---

## 9. File index (`results/`)

| file | content |
|---|---|
| `{ROSMAP,Diverse,Banner,BannerLFQ}_causal_W.csv` | 80 × 80 weights |
| `…_causal_beta.csv` | 80 × 80 effect sizes |
| `…_causal_N.csv` | 80 × 80 people per pair |
| `…_causal_edges.csv` | 688 prior edges with W, β, n, p, q, data_supported |
| `…_variables.csv` | every variable at every selection step |
| `…_regions.csv` | every region (original scoring) |
| `{ROSMAP,Diverse,Banner}_regions_usable_null_check.csv` | regions re-scored with the corrected null |
| `Diverse_best_region_W.csv` | graph refitted in the original top Diverse region |
| `…_summary.json`, `…_report.txt`, `….log`, `…_individuals.txt` | run summaries, logs, people used |
| `comparison_report.txt`, `comparison_pairs.csv`, `consensus_edges.csv` | cross-dataset comparison |
