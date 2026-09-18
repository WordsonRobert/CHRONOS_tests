# Hypothesis 6: APOE4 Stratification
**CHRONOS | Analysis date: September 2026**

---

## Biological claim
APOE4 drives upstream lipid dysfunction (Castello-Serrano 2025) which impairs
TMED2/TMED10 trafficking and reduces the presynaptic active-zone module score.
In ROSMAP bulk RNA, APOE4 carriers should show lower TMED2, TMED10, and AZscore
even after adjusting for AD pathology and cell composition.

---

## Data
- ROSMAP bulk DLPFC RNA, n=616 with APOE4 status
- APOE4+: n=145 (genotype 34 or 44)
- APOE4−: n=471 (genotype 22, 23, or 33)
- Covariates: age, sex, PMI, Braak, rnaBatch, cell composition (neuron/astro/micro/oligo)

---

## Cohort summary

| | APOE4+ | APOE4− |
|--|--------|--------|
| n | 145 | 471 |
| Age (mean±SD) | 86.0±4.6 | 86.2±4.8 |
| Sex (M%) | 37% | 35% |
| Mean Braak | 3.99 | 3.35 |

APOE4+ carriers have higher Braak on average — as expected biologically.

---

## Results

### Primary tests (full cohort)

| Outcome | Raw β | Raw p | Adj β | Adj p | FDR | Full β | Full p |
|---------|-------|-------|-------|-------|-----|--------|--------|
| TMED2 | −0.025 | 0.619 | −0.015 | 0.767 | 0.827 | −0.043 | 0.175 |
| TMED10 | +0.040 | 0.374 | +0.037 | 0.417 | 0.827 | −0.009 | 0.720 |
| TMED9 | −0.003 | 0.927 | +0.019 | 0.606 | 0.827 | +0.009 | 0.690 |
| AZscore | −0.075 | 0.368 | −0.018 | 0.827 | 0.827 | −0.017 | 0.201 |

All FDR = 0.827. Nothing survives.

### Secondary: Log-ratios

| Ratio | β | p |
|-------|---|---|
| log(TMED2/TMED10) | −0.034 | 0.146 |
| log(TMED9/TMED2) | +0.052 | 0.162 |
| log(TMED9/TMED10) | +0.019 | 0.508 |

### Test H: APOE4 × Braak interaction

| Outcome | Interaction β | p |
|---------|--------------|---|
| TMED2 | −0.004 | 0.881 |
| TMED10 | +0.012 | 0.585 |
| TMED9 | +0.017 | 0.384 |
| AZscore | +0.002 | 0.884 |

No interaction — APOE4 does not modify Braak trajectory for any outcome.

### Test I: Pre-symptomatic (Braak≤2)

| Cohort | n | APOE4+ | Outcomes |
|--------|---|--------|---------|
| Braak≤2 | 106 | 16 | All flat, p>0.5 |
| Braak≤2+NCI | 58 | 8 | All flat, p>0.5 |

Severely underpowered — only 8–16 APOE4+ in these subsets.

---

## Verdict

| Claim | Supported? |
|-------|-----------|
| TMED2 lower in APOE4+ | ❌ No |
| TMED10 lower in APOE4+ | ❌ No |
| AZscore lower in APOE4+ | ❌ No |
| APOE4 × Braak interaction | ❌ No |
| Pre-symptomatic APOE4 effect | ❌ Underpowered — cannot conclude |

**Overall: Clean null. APOE4 status does not associate with TMED2/TMED10/TMED9
or the AZ module score in bulk DLPFC RNA after adjustment for pathology and
cell composition.**

---

## Interpretation
The APOE4 lipid phenotype (from NLA/iPSC astrocyte data) may not propagate to
detectable transcriptomic changes in the presynaptic machinery in human brain.
Two explanations: (1) the effect is fully mediated through Braak/pathology load
which we adjust for, removing it from the residual; (2) the effect exists at
protein/lipid level but not at the mRNA level. The ROSMAP bulk RNA test is
negative but cannot distinguish these.

Evidence streams remain separate:
- NLA/iPSC: APOE4 → CE accumulation (lipid level)
- ROSMAP bulk RNA: APOE4 → no detectable TMED/AZ change (transcriptomic level)

---

## Scripts
- `hypothesis6_apoe4.py`

## Output files
- `h6_primary_results.csv`
- `h6_interactions.csv`
- `h6_presymptomatic.csv`
- `h6_ratios.csv`
- `h6_plots.png`
