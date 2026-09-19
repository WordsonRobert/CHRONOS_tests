# Test 2: RNA vs Protein Discordance
## CHRONOS / ROSMAP

---

## Biological Question

If TMED2/TMED10 are post-transcriptionally degraded in AD (the p24 hypothesis), protein should drop *without* proportional mRNA drop. The standardized protein effect should be substantially larger than the RNA effect, and protein abundance should predict pathology independently of transcript levels.

---

## Cohort

- RNA: n = 634
- Protein: n = 400 individuals
- RNA ∩ Protein ∩ metadata (paired): n = 209
- Complete-case paired (with cell scores): n = 151

---

## Results

### Test B: RNA ~ AD (composition-adjusted)

| Gene | β_AD (simple) | p | β_AD (+cells) | p |
|---|---|---|---|---|
| TMED2 | −0.126 | 0.017 | −0.063 | **0.047** |
| TMED10 | −0.027 | 0.548 | −0.057 | **0.023** |
| TMED9 | −0.081 | 0.021 | −0.035 | 0.183 |

Both TMED2 and TMED10 RNA decrease significantly in AD after composition adjustment.

---

### Test C: Protein ~ AD (composition-adjusted)

| Gene | β_AD (simple) | p | β_AD (+cells) | p |
|---|---|---|---|---|
| TMED2 | −0.042 | 0.095 | −0.045 | 0.083 |
| TMED10 | −0.033 | 0.112 | −0.031 | 0.141 |
| TMED9 | −0.041 | 0.028 | −0.042 | **0.029** |

Protein is weaker than RNA for TMED2/TMED10. Only TMED9 protein reaches significance, but in the same direction as RNA — not the divergence the hypothesis predicts.

---

### Test D: Braak — RNA vs Protein (composition-adjusted)

| Gene | RNA β_Braak | RNA p | Protein β_Braak | Protein p |
|---|---|---|---|---|
| TMED2 | −0.007 | 0.520 | −0.014 | 0.115 |
| TMED10 | −0.018 | **0.045** | −0.008 | 0.272 |
| TMED9 | −0.025 | **0.007** | +0.003 | 0.637 |

RNA shows Braak signal for TMED10 and TMED9. Protein is null for all three. The Braak protein signal is weaker than RNA — opposite to the prediction.

---

### Test E: Standardized Discordance (paired n=151)

| Gene | β_RNA_z | β_Protein_z | Δβ | D = |β_P|/|β_R| |
|---|---|---|---|---|
| TMED2 | −0.021 (p=0.844) | −0.312 (p=0.083) | −0.291 | 15.2 |
| TMED10 | +0.005 (p=0.964) | −0.265 (p=0.141) | −0.269 | 58.9 |
| TMED9 | +0.038 (p=0.787) | −0.390 (p=0.029) | −0.427 | 10.4 |

D ratios are large (10–59×), but this is driven by the RNA effect being essentially zero in the paired cohort (n=151), not by a dramatically large protein effect. The denominator is near zero, so D is misleading here.

---

### Test F: Within-Subject RNA ↔ Protein Spearman (paired n=151)

| Gene | ρ (all) | p | ρ (NCI) | ρ (AD) |
|---|---|---|---|---|
| TMED2 | −0.027 | 0.744 | −0.055 | −0.022 |
| TMED10 | −0.083 | 0.309 | −0.042 | −0.166 |
| TMED9 | +0.135 | 0.099 | +0.103 | +0.140 |

RNA and protein are uncorrelated within individuals for all three genes. This is consistent with post-transcriptional biology, but also consistent with measurement noise and the small paired n.

---

### Test G: Protein ~ RNA + Braak (conditional, paired n=151)

| Gene | Protein ~ Braak β | p | Protein ~ Braak+RNA β_Braak | p |
|---|---|---|---|---|
| TMED2 | −0.0023 | 0.840 | −0.0025 | 0.822 |
| TMED10 | −0.0045 | 0.621 | −0.0046 | 0.617 |
| TMED9 | +0.0066 | 0.431 | +0.0066 | 0.430 |

Adding RNA to the protein model changes nothing. RNA doesn't explain protein (β_RNA ≈ 0 for all genes), and protein doesn't move with Braak regardless of whether RNA is included.

---

### Test H: Incremental R² (paired n=151)

| Gene | R²(M0: Braak+cov) | R²(M1: +RNA) | ΔR² |
|---|---|---|---|
| TMED2 | 0.051 | 0.052 | 0.002 |
| TMED10 | 0.046 | 0.046 | 0.000 |
| TMED9 | 0.043 | 0.045 | 0.003 |

RNA adds essentially zero explanatory power to the protein~Braak model. These two measurements are independent of each other.

---

### Test I: Protein Residual ~ Braak/AD (paired n=151)

| Gene | β_resid~Braak | p | β_resid~AD | p |
|---|---|---|---|---|
| TMED2 | −0.002 | 0.839 | −0.040 | 0.092 |
| TMED10 | −0.004 | 0.652 | −0.027 | 0.154 |
| TMED9 | +0.005 | 0.476 | −0.037 | **0.035** |

After removing RNA from protein, the residual (= protein unexplained by transcript) does not associate with Braak for any gene. TMED9 protein residual associates weakly with AD (p=0.035), but Braak is null.

---

## Summary Verdict

**This is a Negative result against the post-transcriptional discordance hypothesis.**

| Criterion | Result |
|---|---|
| Protein decreases with Braak/AD | Weak: only TMED9 AD p=0.029; Braak null for all |
| RNA substantially weaker than protein | **NO** — RNA is *stronger* than protein for TMED2/TMED10 |
| Protein effect survives conditioning on RNA | Moot — protein had no Braak signal to begin with |
| RNA–protein correlation modest | Yes (ρ ≈ 0), consistent with post-transcriptional biology but also with noise |
| Protein residual associates with pathology | No (Braak); weak TMED9 AD only |

The predicted pattern — protein down without RNA down — is not observed. The RNA signal is actually stronger than the protein signal in the AD direction. The large D ratios in Test E are an artifact of near-zero RNA effects in the smaller paired cohort, not evidence of protein-dominant biology.

**This does not refute the p24/TMED hypothesis at the cellular level.** It means that in ROSMAP bulk proteomics, the signal is not detectable in the expected direction. Possible reasons: small paired n (151 with complete cells), proteomics batch effects with no usable batch correction (rnaBatch unavailable), or the effect being genuinely absent at the bulk tissue level.
