# Hypothesis 8: Three-Protein Composite Biomarker Score
## CHRONOS / ROSMAP Bulk RNA Analysis

---

## Biological Claim
If TMED2, TMED10, and TMED9 jointly reflect p24 trafficking dysfunction in AD, combining them into a composite score should outperform any individual gene for AD classification and Braak prediction.

Predicted pattern: TMED2↓ + TMED10↓ + TMED9↑ → composite captures the directional divergence → better signal.

---

## Cohort
- Full ROSMAP bulk RNA: n=634
- AD classification subset (cogdx 1=NCI, ≥4=AD): n=466 (NCI=201, AD=265)

---

## Results

### Test A: H8 Ratio Composite
`H8_ratio = TMED9 − TMED2 − TMED10` (log-scale, equivalent to log(TMED9 / TMED2×TMED10))

| Stat | Value |
|---|---|
| Mean | −5.907 |
| SD | 0.797 |
| corr(ratio, TMED2) | −0.879 (p<0.0001) |
| corr(ratio, TMED10) | −0.862 (p<0.0001) |
| corr(ratio, TMED9) | −0.207 (p<0.0001) |

**Problem:** The ratio is almost entirely driven by TMED2 and TMED10. TMED9's contribution is weak (|r|=0.21). The denominator terms dominate.

---

### Test B: PCA
| Component | Variance explained | TMED2 loading | TMED10 loading | TMED9 loading |
|---|---|---|---|---|
| PC1 | 74.5% | 0.604 | 0.603 | 0.521 |
| PC2 | 18.0% | — | — | — |
| PC3 | 7.5% | — | — | — |

PC1 loads nearly equally on all three genes — it captures general p24 expression level, not the directional divergence (TMED2/TMED10↓ vs TMED9↑) the hypothesis predicts.

---

### Test C: AUC (AD vs NCI)
| Feature | AUC | n |
|---|---|---|
| TMED2 | 0.577 | 466 |
| TMED9 | 0.578 | 466 |
| PC1 | 0.564 | 466 |
| H8_ratio | 0.514 | 466 |
| TMED10 | 0.501 | 466 |

No feature clears a useful threshold. Best individual genes (TMED2, TMED9) ~0.577. The composite ratio is worse than either alone (0.514).

---

### Test D: DeLong Test (H8_ratio vs individual genes)
| Comparison | AUC_ratio | AUC_gene | z | p |
|---|---|---|---|---|
| ratio vs TMED2 | 0.514 | 0.423 | 1.749 | 0.080 |
| ratio vs TMED10 | 0.514 | 0.499 | 0.292 | 0.770 |
| ratio vs TMED9 | 0.514 | 0.422 | 2.196 | 0.028 |

Note: DeLong AUCs here are signed (not flipped), so raw gene AUCs below 0.5 indicate the gene associates in the opposite direction to prediction. The ratio does not significantly outperform TMED2 (p=0.080) and is marginally better than TMED9 raw (p=0.028) but the overall AUCs are all near chance.

---

### Tests E/F: Braak ΔR²
| Feature | R²_null (covariates only) | R²_full | ΔR² |
|---|---|---|---|
| TMED9 | 0.1688 | 0.1786 | **0.0098** |
| PC1 | 0.1688 | 0.1770 | 0.0082 |
| TMED10 | 0.1688 | 0.1742 | 0.0054 |
| TMED2 | 0.1688 | 0.1694 | 0.0006 |
| H8_ratio | 0.1688 | 0.1688 | **0.0000** |

The ratio adds zero Braak variance above covariates. TMED9 individually adds the most (ΔR²=0.0098), still small.

---

### Test G: Three-Gene Joint Model (Braak)
Individual models:
- TMED10: β=−0.365, p=0.045
- TMED9: β=−0.475, p=0.007
- TMED2: β=−0.093, p=0.520

Three-gene joint model (n=633):
- TMED9: β=−0.398, p=0.035 ✓
- TMED10: β=−0.300, p=0.248
- TMED2: β=+0.051, p=0.799

TMED9 is the only term that survives in the joint model. TMED2 loses significance entirely and flips direction — suggests collinearity between TMED2 and TMED10 in the joint model.

---

### Test H: Cross-Validated AUC (5-fold × 10 repeats)
| Feature | CV-AUC | SD |
|---|---|---|
| TMED9 | 0.578 | 0.062 |
| TMED2 | 0.576 | 0.046 |
| PC1 | 0.563 | 0.053 |
| H8_ratio | 0.492 | 0.047 |
| TMED10 | 0.472 | 0.041 |

CV-AUC confirms: the ratio (0.492) performs at chance level. Individual genes slightly better but all near 0.5–0.58, nowhere near clinically useful.

---

### Test I: Permutation Test (H8_ratio, n=500)
- Observed AUC: 0.514
- Permutation p: 0.624

**The composite is not distinguishable from random permutation of labels.**

---

### Test K: Ratio ↔ PC Correlation
- corr(H8_ratio, PC1) = −0.775 (p<0.0001)
- corr(H8_ratio, PC2) = +0.631 (p<0.0001)

The ratio is split across PC1 and PC2, confirming it doesn't align with any clean axis of biological variation captured by PCA of these three genes.

---

## Summary Verdict

**H8 is a clean null.**

| Question | Result |
|---|---|
| Does composite outperform individual genes? | No — ratio AUC 0.514, worse than TMED2/TMED9 alone |
| Is any feature a useful AD classifier? | No — best CV-AUC ~0.578, near chance |
| Does the ratio capture Braak progression? | No — ΔR²=0.000 |
| Is the composite distinguishable from random? | No — permutation p=0.624 |
| Does PCA reveal the predicted divergence pattern? | No — PC1 loads equally on all three genes |

The predicted directional pattern (TMED2/TMED10↓, TMED9↑) does not manifest as a coherent composite signal in ROSMAP bulk RNA. The ratio is dominated by TMED2/TMED10 variance, TMED9 contributes weakly, and combining them adds no discriminative power above individual genes.

TMED9 individually remains the strongest single predictor of Braak in the joint model, but even that effect is small and previously shown to be partially composition-driven (H1).
