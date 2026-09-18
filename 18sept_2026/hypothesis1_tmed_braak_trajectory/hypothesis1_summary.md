# Hypothesis 1: TMED2/TMED10/TMED9 Braak Trajectory
**CHRONOS | Analysis date: September 2026**

---

## Biological claim
Increasing AD pathology (Braak 0→6) drives:
- TMED2 ↓
- TMED10 ↓
- TMED9 ↑ (compensatory stress response)

---

## Data
- ROSMAP bulk DLPFC RNA, n=634, 15,582 genes
- Braak 0–6 as continuous pathology variable
- Cell-type composition estimated from marker genes

---

## Results summary

### Tests A–H (Script 1): Individual trajectories

| Gene | Raw β | Raw p | FDR | Adjusted β | Adjusted p | FDR | After neuron adj β | After neuron adj p |
|------|-------|-------|-----|-----------|-----------|-----|-------------------|-------------------|
| TMED2 | -0.018 | 0.298 | 0.298 | -0.015 | 0.436 | 0.620 | +0.014 | 0.306 |
| TMED10 | +0.018 | 0.232 | 0.298 | +0.008 | 0.620 | 0.620 | +0.021 | 0.165 |
| TMED9 | -0.030 | 0.013 | **0.038** | -0.032 | 0.012 | **0.035** | -0.014 | 0.162 |

**Key finding:** TMED9 shows a weak but FDR-significant negative Braak association.
After neuronal composition adjustment, this signal disappears (p=0.16).
TMED2 and TMED10 show no trajectory whatsoever.

### Nonlinear (Test G)
- TMED2: F-test p=0.19 — no nonlinear effect
- TMED10: F-test p=0.20 — no nonlinear effect
- TMED9: F-test p=0.009 — slight nonlinear acceleration in late Braak

### Pre-symptomatic contrast (Test H)
- TMED2: ANOVA p=0.71 — flat
- TMED10: ANOVA p=0.57 — flat
- TMED9: ANOVA p=0.012 — drop only in late group (Braak 5-6), p=0.005 for mid→late

---

### Tests I–L (Script 2): Inter-gene relationships

| Relationship | β | p | Notes |
|-------------|---|---|-------|
| TMED10 ~ TMED2 | 0.689 | 2.8e-130 | Very strong coupling |
| TMED9 ~ TMED2 | 0.363 | 1.5e-46 | Positive — all move together |
| TMED9 ~ TMED10 | 0.415 | 2.3e-47 | Positive — all move together |
| TMED2:Braak → TMED10 (interaction) | 0.003 | 0.85 | No change in coupling with pathology |
| TMED2:Braak → TMED9 (interaction) | -0.019 | 0.29 | No change |
| log(TMED2/TMED10) ~ Braak | -0.023 | 0.054 | Borderline — TMED2 loses slightly vs TMED10 |
| p24 composite ~ Braak (hypothesis) | +0.026 | 0.155 | Not significant |
| p24 composite ~ Braak (data-driven) | -0.032 | 0.280 | Not significant |

---

## Verdict

| Claim | Supported? |
|-------|-----------|
| TMED2 declines with Braak | ❌ No |
| TMED10 declines with Braak | ❌ No |
| TMED9 rises with Braak | ❌ No (goes slightly down, likely neuronal loss) |
| TMED2/TMED10 coupling changes with pathology | ❌ No |
| p24 composite score tracks Braak | ❌ No |

**Overall: Hypothesis 1 bulk RNA trajectory is not supported.**
The three genes are tightly co-expressed but none show the predicted directional change with AD pathology in bulk DLPFC after composition adjustment.

---

## What's still needed (blocked/future)
- Protein-level replication (ROSMAP proteomics — separate script needed)
- Neuronal pseudobulk from snRNA (Mathys dataset — Seurat slot bug blocks this)
- Cell-type-specific Braak trajectory (requires full snRNA)

---

## Scripts
- `tmed_trajectory_script1.py` — Tests A–H
- `tmed_trajectory_script2.py` — Tests I–L

## Output files
- `tmed_testA_stats.csv`
- `tmed_script1_summary.csv`
- `tmed_testG_nonlinear.csv`
- `tmed_testH_presymptomatic.csv`
- `tmed_script2_summary.csv`
- `tmed_trajectory_script1.png`
- `tmed_testH_presymptomatic.png`
- `tmed_script2.png`
