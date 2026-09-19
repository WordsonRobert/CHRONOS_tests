# Test 5: eRWR Audit + Normalization Fix
## CHRONOS / ROSMAP

---

## Biological Question

Does the TMED2 secretory-support network score (eRWR) contain a disease-associated signal in ROSMAP DLPFC bulk RNA-seq, once the original normalization bug is corrected?

---

## Background: The Original Problem

The Stage 2 eRWR (expression-weighted Random Walk with Restart) script used **per-sample min-max normalization**: for every individual, the expression vector was independently rescaled to [0,1] before running the walk. This made every patient's network look identical in relative terms — all between-subject variance was erased. The result was a score with SD ≈ 0.000084 and no association with any clinical outcome (NCI vs AD p = 0.179, Braak ρ = +0.043 p = 0.281). That was not a biological null; it was a technical artifact.

---

## Data

- Bulk RNA: n = 634 (ROSMAP DLPFC logCPM), shared with clinical covariates
- Network: STRING, 16,525 nodes; TMED2 has 105 network neighbours
- Subgraph: 106 nodes (TARGET + 105 neighbours), 98/106 with expression data
- Cell-type composition scores: neuron, astro, micro, oligo (z-score marker averages)
- Covariates: age_death, msex, pmi, cell scores

---

## Results

### Test A: Implementation Audit (Synthetic Perturbation)

Artificially increased / decreased TMED2-network gene expression by +2 / −2 logCPM to check whether the eRWR score actually responds.

| Normalization | Low score | Base score | High score | Δhigh | Responds? |
|---|---|---|---|---|---|
| per-sample min-max | 0.007668 | 0.007668 | 0.007668 | +0.000000 | **No** |
| global z-score | 0.007505 | 0.007321 | 0.007630 | +0.000308 | **Yes** |
| global percentile | 0.007976 | 0.007976 | 0.007976 | +0.000000 | **No** |

Per-sample min-max and global percentile are both insensitive to perturbation. Only global z-score responds correctly. **Global percentile fails because rank within a patient is unaffected by a uniform additive shift across all genes** — the relative ordering doesn't change. Global z-score is the only normalization that preserves absolute expression differences between patients.

---

### Test B: Normalization Comparison (n = 634)

| Normalization | Score SD | Braak ρ | p | NCI vs AD p |
|---|---|---|---|---|
| per-sample min-max | 0.000084 | +0.043 | 0.281 | 0.179 |
| global z-score | **0.000332** | **−0.095** | **0.016** | **0.0004** |
| global percentile | 0.000990 | −0.098 | 0.014 | 0.004 |

Global z-score has 4× more variance than the broken version and reveals a real signal. Global percentile has even more variance but fails the synthetic audit (Test A), so it is not the primary score. **Primary normalization for Tests C–G: global z-score.**

---

### Test C: Biological Association (global z-score, n = 634)

| Outcome | ρ | p | FDR | n |
|---|---|---|---|---|
| Braak | −0.095 | 0.016 | 0.016 | 634 |
| CERAD | +0.108 | 0.007 | 0.010 | 634 |
| cogdx | −0.139 | 0.0004 | 0.001 | 634 |

All three pathology/cognitive outcomes are significant after BH correction. Direction is consistent across all three:
- Higher eRWR → lower Braak (less tau)
- Higher eRWR → higher CERAD (less amyloid; CERAD 4 = no neuritic plaques)
- Higher eRWR → lower cogdx (better cognition)

**Higher secretory-support network activity is associated with less AD pathology and better cognition.**

NCI vs AD (Mann-Whitney): p = 0.0004, Cohen's d = 0.36, rank-biserial r = −0.20, n_NCI = 201, n_AD = 220. AD subjects show −1.52% lower mean eRWR score than NCI.

**Direction note vs GSE5281:** The original GSE5281 analysis found TMED2 expression was +30.4% higher in AD vs control (hippocampus). Here, the eRWR *network support score* is lower in AD. These are measuring different things — GSE5281 measured raw TMED2 RNA level (which can go up as a stress response or artefact of cell-type shift); eRWR measures how well the entire secretory network is functioning collectively. Lower eRWR in AD is consistent with the network being less able to support secretory trafficking despite TMED2 transcript levels being elevated.

---

### Test D: Composition-Adjusted Models

| Outcome | Model | β | p | R² | n |
|---|---|---|---|---|---|
| Braak | simple | −3.8×10⁻⁴ | 0.005 | 0.013 | 634 |
| Braak | +demo | −3.7×10⁻⁴ | 0.067 | 0.029 | 363 |
| Braak | +demo+cells | −3.3×10⁻⁴ | 0.305 | 0.440 | 363 |
| CERAD | simple | +5.0×10⁻⁴ | 0.002 | 0.015 | 634 |
| CERAD | +demo | +4.8×10⁻⁴ | 0.038 | 0.031 | 363 |
| CERAD | +demo+cells | +4.2×10⁻⁴ | 0.043 | 0.445 | 363 |
| cogdx | simple | −5.4×10⁻⁴ | 0.0002 | 0.022 | 634 |
| cogdx | +demo | −5.1×10⁻⁴ | 0.045 | 0.030 | 363 |
| cogdx | +demo+cells | −3.7×10⁻⁴ | 0.101 | 0.442 | 363 |

(β values are small because eRWR score is on a 0.006–0.008 scale while outcomes are 0–6.)

Key findings:
- **CERAD survives full composition adjustment** (p = 0.043) — amyloid association is not explained by changing cell-type composition.
- **Braak and cogdx lose significance** after cell composition is added (p = 0.305, 0.101). These associations are partially driven by the fact that AD brains lose neurons and gain glia, which also changes the eRWR score.
- The large R² jump from ~0.03 to ~0.44 when cell scores are added confirms that bulk eRWR reflects cell composition heavily — expected for a gene network active in neurons.

---

### Test E: Residualized Score Sensitivity

After removing the portion of eRWR predicted by cell composition (neuron, astro, micro, oligo scores), the residual score shows:

| Outcome | ρ (raw) | ρ (residualized) | p (residualized) |
|---|---|---|---|
| Braak | −0.095 | −0.024 | 0.551 |
| CERAD | +0.108 | +0.027 | 0.498 |
| cogdx | −0.139 | −0.063 | 0.116 |

After removing composition, the residual eRWR has no significant association with any pathology measure. This confirms that **most of the raw eRWR–pathology association is mediated by cell-type composition shifts** in AD brain, not by an independent network-level signal. The exception is the CERAD association in Test D (which survives in the full OLS model but not the residual approach) — these two analyses differ in how composition is modeled, so the CERAD finding requires cautious interpretation.

---

### Test F: Network Null (500 random networks)

| | Value |
|---|---|
| Null mean score | 0.049 |
| Null SD | 0.050 |
| TMED2 mean score | 0.008 |
| Z | −0.830 |
| p | 0.406 |

TMED2's mean eRWR score does not significantly exceed (or differ from) randomly selected networks. The TMED2 network is not unusual in its average walk-back score relative to other gene networks of comparable size. This is expected — the null tests overall mean score level, not disease association. The biologically relevant question is whether TMED2's *variance* is disease-associated, which Tests C–E address.

---

### Test G: Incremental Validity — Does eRWR Add Beyond TMED2 Alone?

| Outcome | Model | β_eRWR | p_eRWR | R² | ΔAIC |
|---|---|---|---|---|---|
| Braak | M3: covariates + eRWR | −265 | 0.305 | 0.184 | +0.92 |
| Braak | M4: covariates + TMED2 + eRWR | −294 | 0.286 | 0.184 | +2.82 |
| CERAD | M3: covariates + eRWR | +479 | 0.043 | 0.105 | −2.19 |
| CERAD | M4: covariates + TMED2 + eRWR | +294 | 0.241 | 0.117 | −4.86 |
| CERAD | M2: covariates + TMED2 | — | — | 0.113 | −5.45 |
| cogdx | M3: covariates + eRWR | −484 | 0.101 | 0.083 | −0.77 |

For **Braak**: eRWR does not add beyond covariates alone (p = 0.305) and AIC increases — eRWR is not independently predicting tau pathology after composition is controlled.

For **CERAD**: eRWR alone is marginally significant (p = 0.043, ΔAIC = −2.19), but TMED2 alone is a better predictor (ΔAIC = −5.45). When both are in M4, eRWR drops to p = 0.241 — eRWR is not capturing amyloid signal beyond what TMED2 already captures.

For **cogdx**: neither eRWR nor TMED2 reaches significance after adjusting for demographics.

**Summary of Test G:** eRWR does not add meaningful predictive information beyond TMED2 RNA alone for any of the three primary outcomes in ROSMAP DLPFC. The network score does not outperform the single gene.

---

## Summary Verdict

| Finding | Result | Interpretation |
|---|---|---|
| Per-sample min-max | SD = 0.000084, all p > 0.17 | **Technically non-informative — normalization bug confirmed** |
| Global z-score responds to perturbation | Δhigh = +0.000308 | Implementation is correct with global normalization |
| Global percentile does not respond | Δhigh = 0 | Rank is insensitive to uniform additive shifts |
| eRWR (globalz) vs Braak | ρ = −0.095, p = 0.016 | Nominal; partially composition-driven |
| eRWR (globalz) vs CERAD | ρ = +0.108, p = 0.007 | Survives composition adjustment (p = 0.043) |
| eRWR (globalz) vs cogdx | ρ = −0.139, p = 0.0004 | Partially composition-driven; disappears in residual |
| NCI vs AD | p = 0.0004, d = 0.36 | Real signal, AD shows lower secretory-network activity |
| After cell-type residualization | All ρ < 0.07, all p > 0.10 | Most signal is mediated by composition shifts |
| vs null networks | Z = −0.83, p = 0.41 | TMED2 network not unusual in absolute score level |
| eRWR beyond TMED2 (Test G) | p > 0.24 in M4 | Network score adds nothing over TMED2 RNA alone |

### Interpretation

The normalization audit confirms the Stage 2 result was a technical artifact. With global z-score normalization, a real (if modest) disease association emerges: lower eRWR in AD subjects, consistent with reduced secretory network capacity in advanced disease. However, most of this signal is mediated by cell-type composition — AD brains lose neurons, and TMED2-network genes are most active in neurons, so the bulk eRWR score naturally falls with neuronal loss. After residualizing for composition, no independent pathology signal remains.

The network score also does not outperform TMED2 RNA expression alone (Test G), meaning the eRWR walk is not capturing additional network-level information beyond what the single gene already provides in DLPFC bulk data.

**What this does not disprove:** The eRWR approach is conceptually sound and showed regional specificity in GSE5281 (hippocampus > other regions). ROSMAP DLPFC is a single region, so regional specificity cannot be tested here. A multi-region or cell-type-resolved (pseudobulk eRWR per cell type) analysis would be a stronger test.

**Current status of H5 (eRWR):** The original score was non-informative due to normalization. The corrected score shows a real but composition-confounded disease association. The secretory network's apparent decline in AD bulk RNA is largely explained by neuronal loss, not an independent change in secretory pathway gene expression within surviving cells. For the paper, this finding should be framed as: eRWR detects the expected bulk shift consistent with neurodegeneration, but does not provide evidence for an independent secretory network dysregulation signal beyond cell composition in DLPFC.
