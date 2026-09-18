# Hypothesis 7: p38-MAPK Route (Route 4)
**CHRONOS | Analysis date: September 2026**

---

## Biological claim
RESET/secretory-pathway dysfunction activates p38-MAPK (MAPK14) in neurons,
contributing to AD pathology. MAPK14 expression should rise with Braak stage.
Key confound: DAM (disease-associated microglia) also upregulates p38-MAPK
independently, making bulk RNA results difficult to interpret.

---

## Data
- ROSMAP bulk DLPFC RNA, n=633, Braak 0–6
- Genes: MAPK14, MAP2K3, MAP2K6, MAPKAPK2
- Cell type scores: neuron, astrocyte, microglia, oligo (marker-gene estimated)

---

## Results

### Test A: MAPK14 ~ Braak

| Model | β_Braak | p |
|-------|---------|---|
| Raw | −0.014 | 0.285 |
| + age/sex/PMI/batch | −0.006 | 0.656 |
| + full cell composition | −0.017 | **0.042** |

MAPK14 goes slightly **down** with Braak in the full model — opposite to prediction.

### Test B: Nested models — microglial confounding

| Model | β_Braak | p |
|-------|---------|---|
| M1: Braak only | −0.014 | 0.285 |
| M2: + age/sex/PMI/batch | −0.006 | 0.656 |
| M3: + neuron/astro/oligo | −0.002 | 0.820 |
| M4: + microglia | −0.017 | 0.042 |

**Critical finding:** Adding microglia (M3→M4) amplifies the Braak coefficient
rather than attenuating it — a 756% change in magnitude. This is the opposite
of the expected DAM confound pattern. It suggests MAPK14 and microglia score
are negatively correlated, and conditioning on microglia reveals a negative
MAPK14-Braak relationship.

### Test C: MAPK14 ~ microglia score
- r(MAPK14, micro_score) = **−0.040**, p=0.309
- Surprisingly, MAPK14 is **negatively** (not positively) correlated with microglia
- This contradicts the expected DAM upregulation pattern

### Test D: MAPK14 ~ neuron score
- r(MAPK14, neuron_score) = **+0.698**, p=1.3×10⁻⁹³
- MAPK14 is very strongly positively correlated with neuronal abundance
- MAPK14 expression in bulk tissue largely tracks neuronal content

### Test E: AD vs NCI
- MAPK14 ~ AD (adj Braak+cells): β=−0.058, p=0.003
- MAPK14 is **lower** in AD — consistent with neuronal loss driving the signal

### Test F: p38 pathway panel ~ Braak (full model, FDR-corrected)

| Gene | β_Braak | p | FDR |
|------|---------|---|-----|
| MAPK14 | −0.017 | 0.042 | 0.056 |
| MAP2K3 | +0.037 | 0.004 | **0.007** |
| MAP2K6 | +0.013 | 0.368 | 0.368 |
| MAPKAPK2 | +0.043 | 3.0e-04 | **0.001** |

Interesting split: MAPK14 itself goes down, but its upstream activator
(MAP2K3) and downstream substrate (MAPKAPK2) go up. This dissociation
is unusual and potentially reflects post-translational regulation
rather than simple transcriptional upregulation.

### Test H: MAPK14 × microglia interaction
- Braak × micro: β=−0.016, p=0.341 — not significant

---

## Verdict

| Claim | Supported? |
|-------|-----------|
| MAPK14 rises with Braak | ❌ No — goes down |
| Bulk MAPK14 driven by microglia | ❌ No — driven by neurons (r=0.698) |
| p38 pathway panel rises with Braak | ⚠️ Partial — MAP2K3/MAPKAPK2 up, MAPK14 down |

**Overall: Not supported as originally framed. MAPK14 goes down with Braak,
primarily because it tracks neuronal content. However, the dissociation between
MAPK14 (down) and MAP2K3/MAPKAPK2 (up) is potentially interesting — it could
reflect increased MAPK14 kinase activity (post-translational) without increased
mRNA, similar to the BACE1 RNA/protein dissociation in H4.**

---

## Important caveat
RNA-seq measures transcript abundance, not kinase activity. p38 activation is
primarily a phosphorylation event. A flat or declining MAPK14 transcript does
NOT rule out increased p38 kinase activity in specific cell types. The rising
MAPKAPK2 (a direct p38 substrate whose expression can be regulated by p38
signaling) alongside flat MAPK14 is the most interesting observation here.
snRNA-seq (blocked) would be needed to resolve cell-type specificity.

---

## Scripts
- `hypothesis7_p38mapk.py`

## Output files
- `h7_results.csv`
- `h7_pathway_panel.csv`
- `h7_plots.png`
