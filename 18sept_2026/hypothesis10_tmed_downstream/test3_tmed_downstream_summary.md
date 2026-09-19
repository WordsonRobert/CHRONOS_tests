# Test 3: TMED → BACE1 / APP / Amyloid (CERAD)
## CHRONOS / ROSMAP

---

## Biological Question

If TMED2/TMED10 regulate APP processing through BACE1, lower TMED should associate with higher BACE1 and APP expression, and lower protein-level TMED should associate with higher amyloid burden (higher CERAD score = more amyloid).

Note: CERAD scores 1–4 where 1 = definite AD neuropathology (most amyloid), 4 = no AD (least amyloid). So the prediction is TMED ↑ → CERAD ↑ (less amyloid), i.e. **positive β** for TMED → CERAD.

---

## Data

- RNA n = 634 (ROSMAP DLPFC logCPM)
- Protein n = 400 (TMT quantitation, paired with clinical)
- Paired RNA+Protein+cells n = 151
- PSEN2 not found in expression matrix — excluded from all RNA tests

---

## Results

### Test A: Raw TMED ↔ Downstream RNA (Spearman)

| TMED | Target | ρ | FDR |
|---|---|---|---|
| TMED2 | BACE1 | **+0.837** | <0.001 |
| TMED2 | APP | **+0.881** | <0.001 |
| TMED2 | PSEN1 | +0.584 | <0.001 |
| TMED2 | ADAM10 | +0.696 | <0.001 |
| TMED10 | BACE1 | +0.652 | <0.001 |
| TMED10 | APP | +0.554 | <0.001 |
| TMED10 | PSEN1 | +0.726 | <0.001 |
| TMED10 | ADAM10 | +0.786 | <0.001 |
| TMED9 | BACE1 | +0.547 | <0.001 |
| TMED9 | APP | +0.568 | <0.001 |
| TMED9 | PSEN1 | +0.259 | <0.001 |
| TMED9 | ADAM10 | +0.203 | <0.001 |

All correlations strongly **positive**. The hypothesis predicted TMED2 ↓ → BACE1 ↑, i.e. a negative correlation. The raw data show the opposite direction for all pairs.

---

### Test B: Composition-Adjusted TMED → Downstream RNA

After controlling for age, sex, PMI, and cell-type composition:

| TMED | Target | β | p | FDR |
|---|---|---|---|---|
| TMED2 | BACE1 | **+0.596** | <0.001 | <0.001 |
| TMED2 | APP | **+0.834** | <0.001 | <0.001 |
| TMED10 | BACE1 | **+0.863** | <0.001 | <0.001 |
| TMED10 | APP | **+0.983** | <0.001 | <0.001 |
| TMED9 | BACE1 | +0.295 | <0.001 | <0.001 |
| TMED9 | ADAM10 | **−0.199** | <0.001 | <0.001 |

All TMED→BACE1 and TMED→APP effects remain strongly **positive** after composition adjustment. The direction is opposite to the mechanistic prediction. TMED9→ADAM10 is the only negative relationship.

---

### Test C: Within-Diagnosis Stratification (Key Pairs)

| TMED | Target | Group | ρ |
|---|---|---|---|
| TMED2 | BACE1 | ALL | +0.837 |
| TMED2 | BACE1 | NCI | +0.830 |
| TMED2 | BACE1 | AD | +0.811 |
| TMED2 | APP | ALL | +0.881 |
| TMED2 | APP | NCI | +0.873 |
| TMED2 | APP | AD | +0.877 |
| TMED9 | BACE1 | ALL | +0.547 |
| TMED9 | BACE1 | NCI | +0.564 |
| TMED9 | BACE1 | AD | +0.544 |

The positive correlations hold within both NCI and AD groups. This is not a between-group artifact — TMED and BACE1/APP co-move within disease state. This is almost certainly a **co-expression / co-regulation signal**, not a TMED-regulates-BACE1 signal.

---

### Test D: TMED RNA ↔ Pathology (Spearman)

| TMED | Outcome | ρ | p | FDR |
|---|---|---|---|---|
| TMED2 | ceradsc | +0.122 | 0.002 | 0.006 |
| TMED2 | braaksc | −0.053 | 0.186 | ns |
| TMED2 | cogdx | −0.118 | 0.003 | 0.007 |
| TMED10 | ceradsc | +0.016 | ns | ns |
| TMED9 | ceradsc | **+0.161** | <0.001 | <0.001 |
| TMED9 | braaksc | −0.155 | <0.001 | <0.001 |
| TMED9 | cogdx | −0.111 | 0.005 | 0.009 |

TMED2 and TMED9 RNA have small but significant associations with CERAD: higher TMED → higher CERAD score → **less amyloid**. TMED9 also shows a negative Braak association (higher TMED9 → lower Braak = less tau pathology). Direction is consistent with the protective hypothesis.

---

### Test E: Composition-Adjusted CERAD Models (RNA)

| TMED | β | p | FDR | R² |
|---|---|---|---|---|
| TMED2 | **+0.370** | 0.009 | 0.009 | 0.102 |
| TMED10 | **+0.665** | <0.001 | <0.001 | 0.112 |
| TMED9 | **+0.517** | 0.003 | 0.004 | 0.105 |

After composition adjustment, all three TMED genes show a **positive, significant association** with CERAD score (= less amyloid). This survives cell-type correction and is the strongest RNA-level result in favor of the hypothesis.

---

### Test F: TMED Protein ↔ CERAD

**Spearman (raw, n=398):**

| TMED | ρ | p |
|---|---|---|
| TMED2 | +0.118 | 0.018 |
| TMED10 | +0.081 | 0.106 |
| TMED9 | +0.036 | ns |

**OLS adjusted (demo + protein batch):**

| TMED | Model | β | p | R² |
|---|---|---|---|---|
| TMED2 | demo+batch | +0.964 | **0.033** | 0.082 |
| TMED2 | demo+batch+cells | +0.323 | 0.630 | 0.090 |
| TMED10 | demo+batch | +1.104 | 0.081 | 0.079 |
| TMED9 | demo+batch | +0.220 | ns | 0.072 |

TMED2 protein reaches significance in the full sample (n=398) but loses it when cell scores are added (paired n=151). Direction is correct (positive = less amyloid) for TMED2 and TMED10.

---

### Test G: Protein ↔ CERAD Within NCI / AD

| TMED | Group | ρ | p |
|---|---|---|---|
| TMED2 | ALL | +0.118 | 0.018 |
| TMED2 | NCI | +0.130 | 0.094 |
| TMED2 | AD | −0.072 | ns |

The ALL-sample TMED2 signal reverses direction in AD (ρ = −0.07) and is only marginal in NCI (ρ = +0.13, p=0.09). This is the pattern of a **between-group artifact**: the positive ALL correlation is mostly driven by NCI vs AD group differences, not a within-group molecular relationship. TMED10 and TMED9 are null in all strata.

---

### Test H: Conditional Protein Model (protein + RNA → CERAD)

| TMED | Model | β_prot | p_prot | R² |
|---|---|---|---|---|
| TMED2 | protein+cov | +0.323 | 0.630 | 0.090 |
| TMED2 | protein+RNA+cov | +0.400 | 0.547 | 0.113 |
| TMED10 | protein+cov | +0.441 | 0.591 | 0.091 |
| TMED9 | protein+cov | −0.791 | 0.378 | 0.094 |

Protein CERAD signal is null in the paired cohort (n=151) whether or not RNA is included. RNA does not confound or explain protein here — both are independently null against CERAD in the paired sample.

---

### Test I: Incremental R² (protein beyond covariates)

Null model (demo + batch): R² = 0.072, n = 398

| TMED | ΔR² | ΔAIC | β | p |
|---|---|---|---|---|
| TMED2 | +0.011 | −2.60 | +0.964 | **0.033** |
| TMED10 | +0.007 | −1.09 | +1.104 | 0.081 |
| TMED9 | +0.000 | +1.89 | +0.220 | ns |

TMED2 protein adds a small but real increment (ΔAIC = −2.6 = better fit). TMED10 is borderline. TMED9 adds nothing.

---

### Test J: BACE1 Protein ↔ TMED Protein

| TMED | ρ | p_sp | β | p_ols |
|---|---|---|---|---|
| TMED2 | −0.091 | 0.071 | −0.091 | 0.091 |
| TMED10 | **−0.257** | <0.001 | **−0.356** | <0.001 |
| TMED9 | −0.163 | 0.001 | −0.200 | 0.012 |

**TMED10 and TMED9 protein are significantly negatively correlated with BACE1 protein** (opposite direction to RNA). This is the mechanistically predicted relationship at the protein level.

**Conditional (protein + RNA → BACE1 protein, paired n=151):**

| TMED | β_prot | p | β_rna | p |
|---|---|---|---|---|
| TMED2 | −0.076 | 0.253 | +0.008 | ns |
| TMED10 | **−0.254** | **0.001** | +0.028 | ns |
| TMED9 | −0.166 | 0.060 | −0.018 | ns |

TMED10 protein → BACE1 protein relationship survives conditioning on TMED10 RNA (β = −0.254, p = 0.001). RNA carries no independent information. This is the cleanest mechanistic signal in the dataset.

---

### Test K: Directionality Matrix

| Relationship | Expected | RNA raw | RNA adj | Prot adj |
|---|---|---|---|---|
| TMED2 → BACE1 | − | **+(***) WRONG** | **+(***) WRONG** | n/a |
| TMED2 → APP | − | **+(***) WRONG** | **+(***) WRONG** | n/a |
| TMED9 → BACE1 | + | +(***) ✓ | +(***) ✓ | n/a |
| TMED2 → CERAD | − (less amyloid) | +(** ) ✓ | +(** ) ✓ | +(* ) ✓ |
| TMED10 → CERAD | − (less amyloid) | +(ns) | +(***) ✓ | +(ns) |
| TMED9 → CERAD | ? | +(***) | +(***) | +(ns) |

---

## Summary Verdict

**Mixed: RNA co-expression is confounded; protein-level signal is mechanistically meaningful.**

| Finding | Result |
|---|---|
| TMED RNA ↔ BACE1/APP RNA direction | **WRONG** — strongly positive, not negative |
| TMED RNA ↔ CERAD (adjusted) | **CORRECT** direction, significant |
| TMED2 protein ↔ CERAD (n=398) | Weak positive (p=0.033), lost with cells |
| TMED2 protein ↔ CERAD within AD | Null / sign reversal |
| TMED10/TMED9 protein ↔ BACE1 protein | **CORRECT** — negative, significant |
| TMED10 protein → BACE1 protein (conditional) | **Survives RNA conditioning** (p=0.001) |

### Interpretation

The TMED↔BACE1/APP RNA correlations are positive and enormous (ρ ≈ 0.65–0.88). This is almost certainly **co-expression**: these genes are all broadly expressed in neurons and co-move with neuronal abundance. It is not evidence for a regulatory relationship.

The protein-level picture is different. TMED10 protein is negatively associated with BACE1 protein (ρ = −0.26, p<0.001) and this holds after conditioning on TMED10 RNA — meaning it is not explained by transcript levels. This is the closest thing in ROSMAP to the mechanistic prediction: higher TMED10 protein → lower BACE1 protein.

The CERAD associations are in the correct direction at both RNA and protein levels, but the stratified analysis (Test G) shows the TMED2 protein signal is mostly between-group (NCI vs AD) rather than a within-disease molecular relationship.

**This does not prove TMED causes amyloid reduction.** ROSMAP is cross-sectional post-mortem data. The correct phrasing is: *Higher TMED10 protein abundance is associated with lower BACE1 protein abundance in bulk DLPFC tissue.*
