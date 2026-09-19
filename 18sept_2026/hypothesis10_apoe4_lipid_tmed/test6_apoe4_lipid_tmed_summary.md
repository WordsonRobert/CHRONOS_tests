# Test 6: APOE4 → Lipid → TMED Bridge
## CHRONOS / ROSMAP

---

## Biological Question

Does APOE4 genotype modulate sphingomyelin biology (SGMS1/2, SMPD1/3, CERT1) in a way that co-varies with TMED2/10/9 expression? Is there evidence for an APOE4-conditioned lipid → TMED signalling axis in human AD brain?

---

## Data

- Bulk RNA: n = 634 (ROSMAP DLPFC logCPM), APOE4 column: `apoe4_count` (0/1/2 alleles)
- Brain lipidomics: n = 384 (ROSMAP DLPFC, Emory LC-MSMS), 92 lipid species total, 6 SM species
- RNA + lipidomics overlap: n = 196 individuals
- Covariates: age_death, msex, pmi, braaksc, cogdx

---

## Results

### Test A: APOE4 → Lipid Genes (RNA, n = 633)

| Gene | ρ | p | FDR |
|---|---|---|---|
| SGMS1 | +0.034 | 0.393 | 0.797 |
| SGMS2 | +0.101 | 0.011 | 0.057 |
| SMPD1 | −0.001 | 0.971 | 0.971 |
| SMPD3 | −0.028 | 0.478 | 0.797 |
| CERT1 | +0.017 | 0.673 | 0.842 |

No lipid gene shows significant association with APOE4 dosage after FDR correction. SGMS2 is the nearest to nominal significance (p = 0.011, FDR = 0.057) with a small positive correlation — more APOE4 alleles, slightly higher SGMS2 (sphingomyelin synthesis enzyme). All others are null.

---

### Test B: APOE4 → TMED Genes (RNA, n = 633)

| Gene | ρ | p | FDR |
|---|---|---|---|
| TMED2 | +0.005 | 0.901 | 0.901 |
| TMED10 | +0.063 | 0.115 | 0.343 |
| TMED9 | −0.011 | 0.783 | 0.901 |

Consistent with prior analysis: APOE4 dosage has no significant association with any TMED gene in ROSMAP DLPFC bulk RNA. All three are null (FDR >> 0.05).

---

### Test C: Lipid Gene ↔ TMED RNA (15 pairwise Spearman, n = 634)

| Lipid gene | TMED gene | ρ | FDR |
|---|---|---|---|
| CERT1 | TMED2 | **+0.887** | <0.001 |
| SMPD1 | TMED9 | **+0.809** | <0.001 |
| SGMS1 | TMED10 | **+0.718** | <0.001 |
| SMPD3 | TMED9 | **+0.718** | <0.001 |
| SGMS1 | TMED2 | **+0.565** | <0.001 |
| SMPD1 | TMED2 | **+0.560** | <0.001 |
| SMPD3 | TMED2 | **+0.580** | <0.001 |
| CERT1 | TMED10 | **+0.574** | <0.001 |
| CERT1 | TMED9 | **+0.461** | <0.001 |
| SMPD1 | TMED10 | **+0.450** | <0.001 |
| SMPD3 | TMED10 | **+0.225** | <0.001 |
| SGMS2 | TMED9 | **−0.340** | <0.001 |
| SGMS2 | TMED2 | **−0.199** | <0.001 |
| SGMS1 | TMED9 | +0.041 | 0.309 |
| SGMS2 | TMED10 | +0.059 | 0.149 |

13 of 15 pairs are FDR-significant. The strongest association in the entire dataset is **CERT1 ↔ TMED2 (ρ = 0.887)** — ceramide transfer protein and TMED2 are almost perfectly co-expressed. SMPD1/3 (sphingomyelin breakdown enzymes) are strongly correlated with TMED9 (ρ = 0.81, 0.72). SGMS2 (sphingomyelin synthesis) shows the only *negative* associations with TMED2 and TMED9.

**Important caveat:** Correlations this strong in bulk RNA (ρ > 0.8) almost always reflect shared cell-type expression patterns rather than direct co-regulation. CERT1 and TMED2 are both active in the secretory pathway and likely co-expressed in the same cell types (neurons, secretory cells). This is co-expression, not a demonstration of functional coupling in the context of APOE4.

---

### Composite Score Associations

| Lipid score | TMED score | ρ | FDR |
|---|---|---|---|
| CERT1 | TMED_score2 | +0.787 | <0.001 |
| CERT1 | TMED_score3 | +0.763 | <0.001 |
| SM_breakdown (SMPD1+SMPD3) | TMED_score3 | +0.673 | <0.001 |
| SM_breakdown | TMED_score2 | +0.509 | <0.001 |
| SM_synthesis (SGMS1+SGMS2) | TMED_score2 | +0.368 | <0.001 |
| SM_synthesis | TMED_score3 | +0.230 | <0.001 |

CERT1 co-expression with TMED is the dominant signal. SM breakdown (SMPD1+SMPD3) also correlates strongly with the TMED composite. SM synthesis (SGMS1+SGMS2) has a weaker positive correlation — note that SGMS2 is individually *negative* with TMED2/9, so the composite mixes directions.

---

### Test D: APOE4 × Lipid-Gene Interaction → TMED (n = 632)

| Lipid gene | TMED gene | β_interaction | p | FDR |
|---|---|---|---|---|
| SGMS1 | TMED10 | +0.113 | 0.078 | 0.555 |
| CERT1 | TMED9 | +0.092 | 0.133 | 0.555 |
| CERT1 | TMED10 | +0.089 | 0.191 | 0.555 |
| … all others | … | … | >0.20 | >0.55 |

No interaction term survives FDR correction. The best p-value is SGMS1×APOE4→TMED10 at p = 0.078 (FDR = 0.555). **The lipid gene / TMED association is not modulated by APOE4 genotype** — the two things correlate strongly (Test C) but that correlation is not different in APOE4 carriers vs non-carriers.

---

### Test E: SM Lipidomics Species ↔ TMED RNA (n = 196, 18 tests)

Only 6 SM species in the lipidomics panel (compared to 62 expected — this dataset has a relatively sparse SM annotation). Results:

| SM species | TMED | ρ | p | FDR |
|---|---|---|---|---|
| SM d38:1 | TMED9 | +0.175 | 0.014 | 0.204 |
| SM d20:0/18:1 | TMED9 | +0.159 | 0.026 | 0.204 |
| SM d36:2-d9 | TMED9 | −0.151 | 0.035 | 0.204 |
| SM d18:1/24:1 | TMED10 | −0.137 | 0.056 | 0.204 |

FDR-significant: **0 of 18**. Three nominal hits, all for TMED9. The lipidomics SM signal for TMED is weak and does not survive multiple testing correction. Sample size (n = 196 overlap) is a limiting factor.

---

### Test F: SM-C18 Species ↔ TMED (n = 196)

One C18-containing SM species identified in the annotation: **SM d20:0/18:1** (sphingomyelin with a C18 acyl chain on a C20 sphingoid base).

| Comparison | ρ | p |
|---|---|---|
| SM d20:0/18:1 ↔ TMED2 | −0.042 | 0.564 |
| SM d20:0/18:1 ↔ TMED10 | −0.101 | 0.158 |
| SM d20:0/18:1 ↔ TMED9 | **+0.159** | **0.026** |
| APOE4 ↔ SM d20:0/18:1 | −0.078 | 0.275 |

The mechanistically relevant SM-C18 species shows no association with TMED2 or TMED10. It has a weak nominal association with TMED9 (p = 0.026) that does not survive correction. APOE4 is not associated with SM-C18 abundance (p = 0.275). **The specific mechanistic prediction (SM-C18 ↔ TMED2) is not supported.**

---

### Test G: Nested Models — APOE4 / SM-C18 / Both → TMED (n = 196)

Using SM d20:0/18:1 as the SM-C18 predictor:

| TMED | Model | Predictor | β | p | R² |
|---|---|---|---|---|---|
| TMED2 | MA: +APOE4 | APOE4 | +0.062 | 0.524 | 0.034 |
| TMED2 | MB: +SM-C18 | SM-C18 | −0.029 | 0.648 | 0.033 |
| TMED2 | MC: +APOE4+SM-C18 | APOE4 | +0.058 | 0.548 | 0.035 |
| TMED2 | MC: +APOE4+SM-C18 | SM-C18 | −0.026 | 0.686 | 0.035 |
| TMED9 | MB: +SM-C18 | SM-C18 | **+0.095** | **0.031** | 0.062 |
| TMED9 | MC: +APOE4+SM-C18 | SM-C18 | +0.096 | 0.030 | 0.062 |

For TMED2 specifically: neither APOE4 nor SM-C18 predicts TMED2 after adjusting for covariates (both p > 0.5, R² ≈ 0.034 — almost entirely driven by covariates). The APOE4 → lipid → TMED2 chain is **not supported** in the lipidomics-RNA overlap sample.

SM-C18 shows a weak independent association with TMED9 (p = 0.030, R² increase from covariates baseline) that is robust to adding APOE4 — but this is for TMED9, not TMED2, and it is marginal.

---

## Summary Verdict

| Finding | Result | Interpretation |
|---|---|---|
| APOE4 → lipid genes (RNA) | All FDR > 0.05; SGMS2 near-nominal | APOE4 does not significantly dysregulate SM gene expression in DLPFC |
| APOE4 → TMED genes | All null | Confirmed prior result |
| Lipid genes ↔ TMED (RNA) | 13/15 pairs FDR < 0.001, ρ up to 0.887 | **Strong co-expression**, but likely cell-type driven |
| APOE4 × lipid interaction | All FDR > 0.55 | APOE4 does not modulate the lipid-TMED relationship |
| SM lipidomics ↔ TMED | 0/18 FDR-significant | No SM abundance signal for TMED in lipidomics |
| SM-C18 ↔ TMED2 | ρ = −0.042, p = 0.564 | **Mechanistic prediction not supported** |
| APOE4 ↔ SM-C18 | ρ = −0.078, p = 0.275 | No APOE4-SM-C18 association |
| Nested models (G) | TMED2 R² ≈ 0.034, all p > 0.5 | SM-C18 and APOE4 are non-predictors of TMED2 |

### Interpretation

The dominant finding in Test 6 is the **very strong co-expression between CERT1, SMPD1/3, and TMED genes** at the RNA level (ρ = 0.46–0.89, all FDR < 0.001). This is biologically coherent: CERT1 (ceramide transfer at ER-Golgi interface), SMPD1/3 (sphingomyelin hydrolysis), and TMED2/9/10 (ER-Golgi vesicle coat) are all part of the secretory pathway and are co-expressed in secretory-active cell types. The co-expression almost certainly reflects shared neuronal/secretory cell expression rather than direct functional coupling.

However, this co-expression is **not APOE4-conditioned** — the interaction test (D) is uniformly null. APOE4 genotype does not significantly change the relationship between lipid metabolism genes and TMED genes.

At the lipidomics level, the actual SM lipid abundance shows no FDR-corrected association with any TMED gene, and the specific mechanistic prediction (SM-C18 ↔ TMED2) is null. The lipidomics sample size (n = 196) limits power, but effect sizes are small (ρ < 0.18).

**What this means for CHRONOS:** The APOE4 → sphingomyelin → TMED2 axis is not supported by ROSMAP transcriptomic or lipidomic data. The strong RNA co-expression between lipid pathway genes and TMED genes is real but is a co-expression phenomenon (shared cell-type biology), not an APOE4-driven signal. For the paper, the honest framing is: lipid metabolism and TMED secretory genes are tightly co-expressed in the brain, consistent with their shared role in the secretory pathway, but APOE4 does not appear to modulate this relationship in DLPFC bulk data.
