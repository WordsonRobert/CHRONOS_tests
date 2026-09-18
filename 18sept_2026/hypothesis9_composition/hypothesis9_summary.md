# Hypothesis 9: Cell-Type Composition Analysis
## CHRONOS / ROSMAP Bulk RNA Analysis

---

## Biological Question
ROSMAP bulk RNA is a tissue mixture. AD brains have fewer neurons, more reactive astrocytes, and altered microglial abundance. Any neuronal gene will appear downregulated with Braak simply because there are fewer neurons — not because of specific biology. This analysis asks: do TMED2, TMED10, and TMED9 have Braak-associated signal *beyond* what is explained by cell-type composition shifts?

---

## Cohort
- n = 634 (all samples with TMED genes + cell scores available)

---

## Results

### Test A: Composition shift across Braak / AD

| Cell type | β Braak | p | β AD | p |
|---|---|---|---|---|
| neuron_score | −0.062 | 0.066 | −0.260 | **0.002** |
| astro_score | +0.121 | **0.0002** | +0.235 | **0.004** |
| micro_score | +0.012 | 0.699 | −0.076 | 0.335 |
| oligo_score | +0.051 | 0.097 | +0.193 | **0.012** |

The expected AD composition shift is confirmed: astrocytes increase significantly with Braak and AD; neurons decrease significantly in AD. Microglia are null in the bulk marker approach — a known limitation of marker-derived deconvolution for microglia. This establishes that bulk RNA composition confounding is real in this dataset.

---

### Test B: Nested TMED models (M1–M4)

| Gene | M1 raw β | M2 +cov β | M3 +neuron β | M4 +all cells β | % change (M1→M4) |
|---|---|---|---|---|---|
| TMED2 | −0.018 (p=0.299) | −0.015 (p=0.436) | +0.014 (p=0.306) | −0.007 (p=0.520) | −60% |
| TMED10 | +0.018 (p=0.231) | +0.008 (p=0.620) | +0.021 (p=0.165) | −0.018 (p=0.045) | −196% |
| TMED9 | −0.030 (p=0.014) | −0.032 (p=0.012) | −0.014 (p=0.162) | −0.025 (p=0.007) | −15% |

**TMED2:** No significant Braak association even before adjustment. Composition adjustment is irrelevant — there is no signal to attenuate.

**TMED10:** No significant raw Braak association. Paradoxically gains marginal significance in M4 (β=−0.018, p=0.045), suggesting suppression effects from correlated cell fractions. Not biologically interpretable as a clean result.

**TMED9:** Only gene with a raw Braak association (p=0.014). Loses significance when neuron_score is added (M3, p=0.162), but survives M4 with all cell types (p=0.007). Partial composition-driven, partial residual signal.

---

### Test C: TMED ↔ cell-fraction correlation matrix

| Gene | r neuron | r astro | r micro | r oligo |
|---|---|---|---|---|
| TMED2 | **0.679** | 0.319 | 0.566 | −0.102 |
| TMED10 | 0.349 | **0.698** | 0.661 | 0.130 |
| TMED9 | **0.640** | 0.234 | 0.423 | −0.321 |

All correlations p<0.001.

- **TMED2** tracks neuronal content most strongly (r=0.679) but also microglia (r=0.566) — its expression is an abundance proxy, not pathway-specific.
- **TMED10** tracks astrocytes most strongly (r=0.698), not neurons. Contradicts the assumption that TMED10 is a neuronal-biology readout.
- **TMED9** tracks neurons most strongly (r=0.640) and anti-correlates with oligodendrocytes (r=−0.321).

---

### Test D: Residualized TMED vs Braak
*(After regressing out all composition + demographic covariates)*

| Gene | r with Braak | p |
|---|---|---|
| TMED2 | −0.023 | 0.556 |
| TMED10 | −0.073 | 0.066 |
| TMED9 | **−0.099** | **0.013** |

After full composition and demographic adjustment, only TMED9 retains a Braak association, and it is weak (r=−0.099). TMED2 is flat. TMED10 is marginal and non-significant.

---

### Test E: Neuronal marker control panel
*(Unrelated neuronal genes as negative controls)*

| Gene | β no cells | p | β +cells | p | % change |
|---|---|---|---|---|---|
| RBFOX3 | −0.066 | 0.008 | −0.019 | 0.314 | −71% |
| SYT1 | −0.066 | 0.015 | −0.009 | 0.205 | −87% |
| SNAP25 | −0.085 | 0.040 | −0.020 | 0.169 | −76% |
| NEFL | −0.095 | 0.002 | −0.038 | 0.009 | −60% |
| MAP2 | −0.046 | 0.097 | −0.004 | 0.732 | −91% |
| SYN1 | −0.087 | 0.009 | −0.029 | 0.004 | −67% |
| CAMK2A | −0.034 | 0.162 | +0.004 | 0.750 | −111% |
| GRIN1 | −0.088 | 0.001 | −0.034 | 0.030 | −62% |
| SLC17A7 | −0.087 | 0.003 | −0.031 | 0.0002 | −65% |

**This is the critical result.** Every unrelated neuronal gene shows a Braak-negative trajectory before adjustment (β ~ −0.07 to −0.09). Most lose significance after cell adjustment. TMED2's unadjusted Braak signal (β=−0.018, p=0.30) is actually *weaker* than these generic neuronal markers — meaning TMED2 is not even a good bulk neuronal marker, let alone showing p24-specific biology.

---

### Test F: Non-neuronal marker controls

| Gene | Class | β raw | β +cells | % change |
|---|---|---|---|---|
| GFAP | Astrocyte | +0.127 | +0.028 | −78% |
| AQP4 | Astrocyte | +0.060 | −0.032 | −153% |
| VIM | Astrocyte | +0.061 | +0.004 | −94% |
| AIF1 | Microglia | −0.044 | −0.068 | +55% |
| CSF1R | Microglia | +0.013 | +0.020 | +52% |
| TMEM119 | Microglia | +0.054 | +0.036 | −35% |

Cell composition adjustment substantially attenuates astrocytic marker Braak associations, confirming the fractions are capturing real biological variation. Microglial markers are inconsistent — consistent with the weak micro_score in Test A.

---

### Test G: Composition decomposition (mediation-style)

| Gene | Total β | Direct β | Indirect β (via neuron) | % via neuron |
|---|---|---|---|---|
| TMED2 | −0.015 | +0.014 | −0.029 | **198%** |
| TMED10 | +0.008 | +0.021 | −0.013 | −163% |
| TMED9 | −0.032 | −0.014 | −0.018 | **57%** |

For TMED2: ~198% of the (already null) total effect passes through neuronal composition — the direct effect actually flips sign. The total effect is near zero because composition-driven and direct paths partially cancel.

For TMED9: ~57% of the total Braak association is mediated via neuronal composition. The remaining ~43% is direct (β=−0.014), which corresponds to the weak residual signal seen in Test D.

---

### Test I: Pre-symptomatic subset (Braak ≤ 2, n=110)

| Gene | β no cells | β +cells | % change |
|---|---|---|---|
| TMED2 | +0.028 | +0.041 | +47% |
| TMED10 | +0.133 | +0.086 | −35% |
| TMED9 | +0.029 | +0.059 | +100% |

In the Braak ≤ 2 subset, all three genes show *positive* β with Braak before and after adjustment — opposite to the predicted direction. This is a small sample (n=110) and the estimates are noisy, but there is no evidence of early TMED downregulation in low-pathology tissue after composition control.

---

### Test J: Composition-adjusted H8 composite

| Model | β | p |
|---|---|---|
| H8_ratio ~ Braak | −0.026 | 0.355 |
| H8_ratio ~ Braak + cells | −0.000 | **0.995** |
| H8_ratio ~ AD | +0.072 | 0.364 |
| H8_ratio ~ AD + cells | +0.085 | 0.128 |

The composite is entirely null before adjustment and achieves β=−0.000 after composition adjustment. There is no residual composite signal.

---

### Test K: TMED vs NeuronMarkerScore (killer control)

| Gene | β alone | p | β + NeuronScore | p |
|---|---|---|---|---|
| TMED2 | −0.067 | 0.436 | +0.178 | 0.122 |
| TMED10 | +0.048 | 0.620 | +0.172 | 0.097 |
| TMED9 | −0.313 | **0.012** | −0.144 | 0.378 |

When NeuronMarkerScore (built from 9 unrelated neuronal genes) is added, TMED9 loses significance (p goes from 0.012 to 0.378). This means TMED9's Braak association is largely explained by the same neuronal depletion signal captured by any set of neuronal genes — it does not add information beyond general neuronal content.

---

## Summary Verdict

**H9 is the most important methodological result in the entire CHRONOS ROSMAP analysis.**

| Finding | Result |
|---|---|
| Is there an AD composition shift in this data? | Yes — astrocyte ↑ (p=0.0002), neuron ↓ in AD (p=0.002) |
| Does TMED2 have Braak signal beyond composition? | No — null even before adjustment; r=−0.023 residual |
| Does TMED10 have Braak signal beyond composition? | No — null raw; unstable across models |
| Does TMED9 have Braak signal beyond composition? | Weak yes — r=−0.099 residual, but lost vs NeuronMarkerScore |
| Is TMED2's Braak pattern distinguishable from generic neuronal depletion? | No — weaker than RBFOX3, SYT1, SNAP25 |
| Does the composite survive composition adjustment? | No — β=0.000, p=0.995 |

**Conclusion:** In ROSMAP bulk RNA, none of the three TMED genes show Braak-associated biology that is distinguishable from cell-type composition shifts. TMED2 and TMED10 are essentially cell-abundance proxies (neuronal and astrocytic respectively). TMED9 has a weak residual signal that is fully absorbed by a generic neuronal depletion score. This does not refute the p24/TMED hypothesis at the cellular level — it means bulk RNA is the wrong assay to test it. Single-nucleus RNA-seq (within cell-type pseudobulk) would be required to properly test the hypothesis.
