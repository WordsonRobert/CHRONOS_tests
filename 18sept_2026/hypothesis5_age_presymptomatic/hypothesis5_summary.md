# Hypothesis 5: Age-dependent Decline in Pre-symptomatic Window
**CHRONOS | Analysis date: September 2026**

---

## Biological claim
Liu 2014 proposes that TMED9/TMED10 expression declines from a postnatal peak
through normal aging, reducing the p24 brake on γ-secretase before AD pathology
develops. Therefore in people with minimal AD pathology (Braak≤2), TMED9 and
TMED10 should decline with age independently of residual pathology and cell
composition.

**Prediction: β_age < 0 for TMED9 and TMED10 in Braak≤2 cohort.**

---

## Data
- ROSMAP bulk DLPFC RNA, n=633 total
- Primary cohort: Braak≤2, n=110 (age 67–90, mean 82.8±6.0)
- Strict cohort: Braak≤2 + NCI (cogdx=1), n=59 (age 67–90, mean 81.9±6.0)
- Pathology covariates: braaksc, ceradsc
- Cell composition: marker-gene estimated scores (neuron, astrocyte, microglia, oligo)

---

## Cohort breakdown (Primary, Braak≤2)

| Braak | n |
|-------|---|
| 0 | 7 |
| 1 | 50 |
| 2 | 53 |

| cogdx | n |
|-------|---|
| NCI (1) | 59 |
| MCI (2) | 26 |
| Other (3-6) | 25 |

---

## Results

### Tests B/C/D: Gene ~ Age (primary cohort)

| Gene | Raw β | Raw p | Adj β | Adj p | FDR | +Path β | +Path p | +Cells β | +Cells p |
|------|-------|-------|-------|-------|-----|---------|---------|----------|----------|
| TMED9 | +0.015 | 0.006 | +0.013 | 0.014 | 0.041 | +0.013 | 0.014 | +0.002 | 0.546 |
| TMED10 | +0.018 | 0.013 | +0.016 | 0.035 | 0.052 | +0.013 | 0.078 | +0.002 | 0.704 |
| TMED2 | +0.015 | 0.079 | +0.013 | 0.135 | 0.135 | +0.012 | 0.155 | −0.006 | 0.247 |

**Direction is OPPOSITE to prediction — all positive slopes, not negative.**

### Test F: Strict NCI-only

| Gene | β | p | n |
|------|---|---|---|
| TMED9 | +0.004 | 0.573 | 59 |
| TMED10 | +0.008 | 0.375 | 59 |
| TMED2 | −0.003 | 0.724 | 59 |

Completely flat in the cleanest cohort.

### Test G: Nested models (cell composition effect)

| Gene | M1 β_age (age only) | M1 p | M2 β_age (+pathology) | M2 p | M3 β_age (+cells) | M3 p |
|------|--------------------|----- |----------------------|------|-------------------|------|
| TMED9 | +0.015 | 0.006 | +0.013 | 0.014 | +0.002 | 0.546 |
| TMED10 | +0.018 | 0.013 | +0.013 | 0.078 | +0.002 | 0.704 |
| TMED2 | +0.015 | 0.079 | +0.012 | 0.155 | −0.006 | 0.247 |

The apparent age signal completely disappears once cell composition is added (M3).
This means the positive age slope was driven by cellular composition differences
across age groups, not by the genes themselves.

### Test H: Nonlinearity
- TMED9: F-test p=0.942 — purely linear (and flat)
- TMED10: F-test p=0.770 — purely linear (and flat)
- TMED2: F-test p=0.261 — no nonlinear component

### Test I: Joint aging signature
- Primary (Braak≤2): aging_p24 ~ age: β=+0.006, p=0.581
- Strict NCI: aging_p24 ~ age: β=+0.013, p=0.410

No joint signal in either cohort.

---

## Verdict

| Claim | Supported? |
|-------|-----------|
| TMED9 declines with age in Braak≤2 | ❌ No — goes slightly UP before composition adjustment |
| TMED10 declines with age in Braak≤2 | ❌ No — same |
| Age effect survives cell composition adjustment | ❌ No — completely disappears in M3 |
| Effect present in strict NCI cohort | ❌ No — flat |
| Joint TMED9/TMED10 score declines with age | ❌ No |

**Overall: Clean null. The pre-symptomatic aging decline hypothesis is not
supported in bulk ROSMAP. The apparent age association (before composition
adjustment) is positive and driven by cellular composition, not biology.
In the cleanest NCI-only cohort, there is no age signal at all.**

---

## Important caveat
ROSMAP participants range from age 67–90+. This dataset cannot capture the
postnatal peak that Liu 2014 describes — it only observes the elderly end
of the lifespan. A true postnatal-peak-to-decline trajectory would require
a lifespan dataset spanning younger ages. The null result here is therefore
specific to the elderly window and does not rule out the existence of a
postnatal peak followed by an earlier decline.

---

## Scripts
- `hypothesis5_age_presymptomatic.py`

## Output files
- `h5_results.csv`
- `h5_nested_models.csv`
- `h5_plots.png`
