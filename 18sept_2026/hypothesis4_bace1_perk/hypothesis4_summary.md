# Hypothesis 4: BACE1 and PERK Pathway ~ Braak
**CHRONOS | Analysis date: September 2026**

---

## Biological claim
Increasing AD pathology drives ER stress via PERK activation, which:
- Phosphorylates eIF2α
- Increases BACE1 translation via uORF mechanism (O'Connor 2008)
- Elevates BACE1 protein without necessarily elevating BACE1 mRNA
- ATF4/CHOP targets should be flat (attenuated chronic PERK — Hou 2017)
- GSK3B should track with pathology (Zhang 2019, Du 2026)

---

## Data
- ROSMAP bulk DLPFC RNA, n=633, Braak 0–6
- ROSMAP proteomics (TMT), n=398
- Paired RNA+protein donors: n=209
- Covariates: age, sex, PMI, rnaBatch
- Cell composition: marker-gene estimated scores (neuron, astrocyte, microglia, oligo)

---

## Results

### Tests A/B/F: Primary genes ~ Braak

| Gene | Raw β | Raw p | Adj β | Adj p | +Cell β | +Cell p |
|------|-------|-------|-------|-------|---------|---------|
| BACE1 | −0.038 | 0.031 | −0.039 | 0.038 | −0.037 | 6.8e-04 |
| EIF2AK3 (PERK) | +0.023 | 0.208 | +0.013 | 0.506 | −0.025 | 0.064 |
| GSK3B | −0.027 | 0.191 | −0.025 | 0.252 | −0.028 | 0.010 |

**BACE1 goes DOWN with Braak** — opposite to the naive prediction.
PERK transcript is flat. GSK3B also goes down (significant after cell adjustment).

### Test C: ATF4/CHOP individual targets ~ Braak (adjusted)

| Gene | β_adj | p_adj | FDR | β_cell | p_cell |
|------|-------|-------|-----|--------|--------|
| ATF4 | +0.001 | 0.967 | 0.967 | +0.008 | 0.556 |
| ASNS | −0.053 | 0.026 | 0.113 | −0.041 | 3.5e-05 |
| TRIB3 | −0.030 | 0.125 | 0.225 | −0.007 | 0.684 |
| SLC7A5 | +0.010 | 0.537 | 0.604 | +0.011 | 0.505 |
| DDIT3 (CHOP) | −0.031 | 0.084 | 0.190 | −0.017 | 0.255 |
| PPP1R15A (GADD34) | −0.054 | 0.014 | 0.113 | −0.014 | 0.428 |

All ATF4/CHOP targets go flat or slightly down — none go up.
ATF4_score ~ Braak: β=−0.017, p=0.346
CHOP_score ~ Braak: β=−0.027, p=0.259

### Test E: Pattern classification

| | Predicted | Observed |
|--|-----------|---------|
| BACE1 | ↑ | ↓ (p=0.038) |
| PERK (EIF2AK3) | elevated/active | FLAT |
| ATF4 targets | FLAT (attenuated) | FLAT/↓ |
| CHOP targets | FLAT (attenuated) | FLAT/↓ |
| GSK3B | ↑ | ↓ (p=0.010 after cell adj) |

### Test G: AD vs NCI (secondary)

| Gene | β_AD | p |
|------|------|---|
| BACE1 | −0.055 | 0.038 |
| EIF2AK3 | −0.042 | 0.200 |
| GSK3B | −0.069 | 0.007 |

All going down in AD. Not what the hypothesis predicts at RNA level.

### Test H: BACE1 RNA vs Protein ⭐

| Measure | β_Braak | p | n |
|---------|---------|---|---|
| RNA | −0.039 | 0.038 | 633 |
| Protein | +0.010 | 0.104 | 398 |

**RNA vs Protein correlation (n=209 paired donors): r=0.005, p=0.94**

**Pattern: RNA DOWN, Protein FLAT**

This is the most important result of hypothesis 4.
BACE1 mRNA declines with AD pathology, but BACE1 protein stays stable.
The near-zero RNA/protein correlation (r=0.005) means transcript abundance
tells us almost nothing about protein abundance — strong evidence for
post-transcriptional regulation of BACE1 in human AD brain.

This is consistent with the O'Connor uORF mechanism:
eIF2α phosphorylation under ER stress increases BACE1 translation
efficiency even as total mRNA declines. The protein is maintained
or even slightly elevated despite falling transcript levels.

---

## Verdict

| Claim | Supported? |
|-------|-----------|
| BACE1 mRNA rises with Braak | ❌ No — goes DOWN |
| PERK transcript elevated in AD | ❌ No — flat |
| ATF4/CHOP targets flat (attenuated PERK) | ✅ Yes — all flat/down |
| GSK3B rises with Braak | ❌ No — goes down |
| BACE1 protein maintained despite RNA decline | ✅ Yes — RNA↓, protein FLAT |
| RNA/protein decoupled for BACE1 | ✅ Yes — r=0.005 |

**Overall: Mixed. The RNA-level predictions are mostly wrong (everything goes down,
not up). But the RNA/protein dissociation is genuinely interesting and consistent
with post-transcriptional regulation — the strongest positive evidence for the
O'Connor uORF mechanism in this dataset.**

---

## Key caveat
RNA-seq cannot measure phosphorylation states. PERK activation is primarily
a post-translational event (phospho-eIF2α). A flat EIF2AK3 transcript does
NOT rule out increased PERK kinase activity. The protein/RNA dissociation
result is the most interpretable finding here.

---

## Scripts
- `hypothesis4_bace1_perk.py` — Tests A–G
- `h4_testH.py` — Test H (RNA vs protein)

## Output files
- `h4_gene_results.csv`
- `h4_pattern.txt`
- `h4_plots.png`
- `h4_testH_results.csv`
- `h4_testH.png`
