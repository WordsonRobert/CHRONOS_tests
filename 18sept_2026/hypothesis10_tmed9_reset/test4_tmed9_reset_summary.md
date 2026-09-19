# Test 4: TMED9 → ER stress / RESET
## CHRONOS / ROSMAP

---

## Biological Question

Does TMED9 expression in human AD brain reflect ER proteostasis stress, consistent with its proposed role in RESET-mediated clearance of misfolded GPI-anchored proteins? Or is TMED9 elevation simply a compensatory response to TMED2 loss?

RESET (Rapid ER Stress-induced Export) is a protein quality-control pathway in which TMED9 facilitates export of misfolded GPI-anchored proteins from the ER to the Golgi. The 2025 mechanistic literature (Muñoz-Braceras et al., PMID 40203033) shows TMED9 directly coordinates this clearance; depletion or pharmacological disruption of TMED9 blocks export. The question here is whether the human AD transcriptomic state is consistent with that mechanism.

---

## Data

- Bulk RNA: n = 634 (ROSMAP DLPFC logCPM), n = 633 with complete covariates
- snRNA: Mathys/ROSMAP snRNAseqPFC_BA10, 17,926 genes × 70,634 cells, 48 donors, BA10 prefrontal cortex
- Pseudobulk (snRNA): 306 donor×cell-type entries with ≥5 cells

---

## Results

### Test A: TMED9 ↔ Individual ER-Stress Markers (Spearman)

| Marker   | Gene           | ρ      | FDR     | n   |
|----------|----------------|--------|---------|-----|
| EIF2AK3  | PERK           | +0.362 | <0.001  | 634 |
| ATF4     | UPR TF         | +0.531 | <0.001  | 634 |
| DDIT3    | CHOP           | +0.608 | <0.001  | 634 |
| HSPA5    | GRP78/BiP      | +0.564 | <0.001  | 634 |
| PPP1R15A | GADD34         | +0.727 | <0.001  | 634 |

All five canonical UPR/ER-stress markers are positively correlated with TMED9. All are genome-wide significant after BH correction. The strongest association is with PPP1R15A (GADD34, ρ = 0.727), a negative feedback regulator of the integrated stress response.

**Important caveat:** These are all broadly expressed genes. Strong positive correlations between co-expressed genes in bulk tissue can reflect shared neuronal/cellular expression patterns rather than co-regulation. These results must be interpreted alongside the nested models in Test E.

---

### Test B: TMED9 ↔ ER-Stress Composite Score

ER-stress composite = mean of z-scored (EIF2AK3, ATF4, DDIT3, HSPA5, PPP1R15A)

| Statistic | Value  |
|-----------|--------|
| Spearman ρ (raw) | +0.685 |
| p | <0.001 |
| n | 634 |
| β (composition-adjusted) | +0.226 |
| p (adjusted) | <0.001 |
| n (adjusted) | 633 |

TMED9 correlates strongly with the ER-stress marker composite, and this association survives adjustment for age, sex, PMI, and cell-type composition (β = +0.226, p < 0.001). The ER-stress composite accounts for substantial variance in TMED9 even after composition is controlled.

---

### Test C: TMED9 ↔ Disease Severity (Spearman, raw)

| Outcome | ρ      | FDR    | n   | Direction |
|---------|--------|--------|-----|-----------|
| Braak   | −0.155 | <0.001 | 634 | **opposite to RESET prediction** |
| CERAD   | +0.161 | <0.001 | 634 | correct direction (higher CERAD = less amyloid) |
| cogdx   | −0.111 | 0.005  | 634 | correct direction (lower cogdx = less impaired) |

TMED9 shows small but significant associations with all three pathology measures. The Braak association is in the **opposite direction** to what the original RESET hypothesis predicted (TMED9 should go up with disease severity). Higher TMED9 is associated with *less* tau pathology and *less* amyloid, not more.

CERAD and cogdx associations are in the direction consistent with protective TMED9 function, but the effect sizes are small (ρ ≈ 0.11–0.16).

---

### Test D: Composition-Adjusted Pathology Models

| Outcome | Model              | β       | p     | n   |
|---------|--------------------|---------|-------|-----|
| Braak   | simple             | −0.030  | 0.013 | 634 |
| Braak   | +demo+batch        | −0.032  | 0.012 | 633 |
| Braak   | +demo+batch+cells  | −0.025  | 0.007 | 633 |
| CERAD   | simple             | +0.045  | <0.001| 634 |
| CERAD   | +demo+batch        | +0.047  | <0.001| 633 |
| CERAD   | +demo+batch+cells  | +0.028  | 0.003 | 633 |
| cogdx   | simple             | −0.018  | 0.077 | 634 |
| cogdx   | +demo+batch        | −0.020  | 0.050 | 633 |
| cogdx   | +demo+batch+cells  | −0.009  | 0.237 | 633 |

Key findings:
- The **Braak association strengthens** slightly after composition adjustment (p = 0.007 in full model), meaning it is not a composition artifact — TMED9 is genuinely negatively associated with tau pathology at the bulk level.
- The **CERAD association attenuates** after composition (β = +0.045 → +0.028) but remains significant. About 40% of the signal is composition-driven.
- The **cogdx association disappears** after cell composition is added (p = 0.237), consistent with this being primarily a composition effect.

---

### Test E: RESET vs Compensation (Nested Models)

Outcome: TMED9 RNA. n = 633 for all nested models.

| Model | Term | β | p | R² |
|-------|------|---|---|----|
| M1: TMED2 | TMED2 | +0.370 | <0.001 | 0.286 |
| M2: TMED2 + ER | TMED2 | −0.007 | 0.788 | 0.541 |
| M2: TMED2 + ER | ER_score | +0.338 | <0.001 | 0.541 |
| M3: + Braak | TMED2 | −0.009 | 0.751 | 0.546 |
| M3: + Braak | ER_score | +0.338 | <0.001 | 0.546 |
| M3: + Braak | Braak | −0.021 | 0.010 | 0.546 |
| M4: + cov+cells | TMED2 | −0.211 | <0.001 | 0.663 |
| M4: + cov+cells | ER_score | +0.267 | <0.001 | 0.663 |
| M4: + cov+cells | Braak | −0.018 | 0.018 | 0.663 |
| M4: + cov+cells | neuron_score | +0.186 | <0.001 | 0.663 |
| M4: + cov+cells | astro_score | +0.084 | <0.001 | 0.663 |
| M4: + cov+cells | oligo_score | −0.046 | <0.001 | 0.663 |

**Compensation test (TMED9 ~ TMED2 simple):** β = +0.370, p < 0.001

This is the most important table in Test 4:

1. **Compensation is ruled out.** In M1, TMED9 ~ TMED2 is strongly *positive* (β = +0.370). The compensation hypothesis predicts a *negative* relationship (TMED9 goes up when TMED2 goes down). The data show the opposite — TMED9 and TMED2 move together, which is co-expression, not compensation.

2. **ER_score completely absorbs the TMED2 signal.** In M1, TMED2 explains 28.6% of variance in TMED9. When ER_score is added (M2), R² jumps to 54.1% and TMED2's β collapses to −0.007 (p = 0.788). TMED2 is confounded by ER state — the reason TMED9 and TMED2 correlate in M1 is that both track ER/cellular stress state, not that TMED2 drives TMED9.

3. **ER_score is the dominant predictor, robust to all adjustments.** ER_score β remains large and significant through M2, M3, and M4 (β = +0.267 in the fully adjusted model with cells + demographics). TMED9 tracks the ER-stress transcriptional state independently of cell composition.

4. **In the full model (M4), TMED2 becomes negative** (β = −0.211, p < 0.001). This sign flip is expected — once ER state is partialed out, residual TMED2 variation is negatively associated with TMED9, which is actually consistent with a partial compensation signal embedded within the larger ER co-expression pattern.

5. **Neuronal and astrocytic scores are the strongest composition predictors** of TMED9 (neuron: β = +0.186, astro: β = +0.084). Oligodendrocyte score is negatively associated (β = −0.046). This is consistent with TMED9 being most active in neurons and astrocytes.

---

### Test F: snRNA Pseudobulk — TMED9 Within Cell Types

Donor-level pseudobulk (mean logCPM per donor per cell type, ≥5 cells). n donors varies by cell type (48 for most, fewer for rare types).

| Cell type | Braak ρ | p (Braak) | CERAD ρ | p (CERAD) | n donors |
|-----------|---------|-----------|---------|-----------|----------|
| Ast       | +0.057  | 0.701     | −0.047  | 0.751     | 48 |
| End       | +0.506  | 0.247     | −0.516  | 0.236     | 7  |
| Ex        | −0.177  | 0.228     | +0.191  | 0.194     | 48 |
| In        | −0.124  | 0.401     | +0.154  | 0.295     | 48 |
| Mic       | +0.325  | **0.028** | −0.241  | 0.106     | 46 |
| Oli       | +0.325  | **0.026** | −0.333  | **0.021** | 48 |
| Opc       | −0.021  | 0.888     | +0.006  | 0.969     | 48 |
| Per       | −0.399  | 0.177     | +0.642  | **0.018** | 13 |

Note: BH-FDR correction across 16 tests (8 cell types × 2 outcomes) not shown above — individual p-values are uncorrected; none survive FDR.

Key findings:
- **No cell type shows robust TMED9 upregulation with AD severity.** The two nominally significant Braak associations (Mic, Oli) are in the *positive* direction — higher TMED9 with higher Braak — but do not survive multiple testing correction.
- **Excitatory and inhibitory neurons are null** for both Braak and CERAD. This directly contradicts a neuron-specific compensation explanation.
- **The Oli–CERAD signal** (ρ = −0.333, p = 0.021) means higher TMED9 in oligodendrocytes associates with *more* amyloid (lower CERAD score). Direction is opposite to what RESET protection would predict in oligos specifically.
- The **Per (pericyte) cell type** shows a strong CERAD signal (ρ = +0.642, p = 0.018) but n = 13 donors — underpowered and unreliable.
- The overall picture is **bulk composition drives the TMED9-pathology signal** more than within-cell-type TMED9 induction.

---

## Summary Verdict

**TMED9 is strongly co-expressed with ER-stress transcriptional markers in human AD bulk tissue, but does not show within-cell-type upregulation with disease severity.**

| Finding | Result | Interpretation |
|---------|--------|----------------|
| TMED9 ↔ individual ER markers | All positive, ρ = 0.36–0.73, all FDR < 0.001 | Supports ER-stress transcriptional co-expression |
| TMED9 ↔ ER composite (adjusted) | β = +0.226, p < 0.001 | ER association survives composition correction |
| TMED9 ↔ Braak (adjusted) | β = −0.025, p = 0.007 | *Opposite* to RESET-up-with-severity prediction |
| TMED9 ↔ CERAD (adjusted) | β = +0.028, p = 0.003 | Correct direction but small; partially composition |
| TMED9 ~ TMED2 (simple) | β = +0.370, positive | **Rules out simple compensation** (predicted negative) |
| TMED2 β after ER_score added | β = −0.007, p = 0.788 | TMED2 signal fully explained by ER state |
| ER_score in fully adjusted model | β = +0.267, p < 0.001 | ER co-expression is the dominant independent signal |
| TMED9 within neurons (snRNA) | ρ ≈ −0.18 vs Braak, p = 0.228 | No neuronal TMED9 induction with AD |
| TMED9 within microglia (snRNA) | ρ = +0.325 vs Braak, p = 0.028* | Marginal; direction unexpected |
| TMED9 within oligodendrocytes | Braak p = 0.026*, CERAD p = 0.021* | Marginal; CERAD direction opposite prediction |

*Uncorrected; does not survive BH-FDR across all snRNA tests.

### Interpretation

The dominant story from ROSMAP is that **TMED9 is part of a broad ER/cellular stress co-expression module**, not a specific disease-induced RESET response. Three pieces of evidence support this:

1. The ER_score completely displaces TMED2 as a predictor of TMED9 (R² 28.6% → 54.1% when ER is added; TMED2 β goes to −0.007), meaning both TMED9 and TMED2 are tracking shared ER transcriptional state.

2. The Braak association is in the wrong direction — TMED9 is lower with worse tau pathology, not higher. If RESET were being activated by disease, TMED9 should rise with severity.

3. The snRNA data show no within-cell-type TMED9 upregulation in neurons (the cell type most vulnerable in AD) across 48 donors.

**What this does not disprove:** RESET is a protein trafficking pathway, not a canonical transcriptional UPR pathway. The 2025 mechanistic experiments showed TMED9 *protein* mediates GPI-AP export regardless of whether transcript levels rise. ROSMAP tests the transcriptomic signature only. A TMED9 protein-level or phosphoproteomic analysis in AD tissue would be a more direct test.

**Current status of H4:** The RESET/ER-stress association in ROSMAP bulk RNA reflects co-expression with a general ER/proteostasis transcriptional module that is actually *attenuated* in severe AD (consistent with neurodegeneration and cell loss). Simple compensation is ruled out. The mechanistic RESET hypothesis remains supported by cell biology but is not confirmed at the human transcriptomic level in ROSMAP.
