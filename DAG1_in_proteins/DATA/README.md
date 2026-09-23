# DATA — the 14 per-variable CSVs

One CSV per clinical/trajectory variable, stacking the three cohorts (ROSMAP, Diverse Cohorts, Banner-TMT).
Every row is **one protein in one subject**:

`protein, dataset, individualID, expr_raw_log2, expr_z, <variable>, impaired`

- `expr_z` — log2 abundance z-scored within its own dataset and protein (use this across cohorts)
- `<variable>` — the variable the file is named after, harmonised to identical labels across cohorts
- `impaired` — 1 = AD or any cognitive decline, 0 = cognitively normal

| file | variable | range | cohorts |
|---|---|---|---|
| `sex.csv` | sex | male / female | all 3 |
| `apoe_genotype.csv` | APOE genotype | 22, 23, 24, 33, 34, 44 (unordered) | all 3 |
| `cts_mmse30_lv.csv` | last MMSE | 0–30 | ROSMAP, Banner |
| `cts_mmse30_first_ad_dx.csv` | MMSE at first AD diagnosis | 0–30 | ROSMAP |
| `pmi.csv` | post-mortem interval (h) — technical | ~0.75–98 | all 3 |
| `braaksc.csv` | Braak stage | 0–6 | all 3 |
| `ceradsc.csv` | CERAD plaque severity (ROSMAP crosswalked) | 0–3 | ROSMAP, Banner |
| `amyThal.csv` | Thal amyloid phase | 1–5 | Diverse |
| `amyA.csv` | binned Thal phase | 1–3 | Diverse |
| `amyCerad.csv` | CERAD (Diverse coding) | 0–3 | Diverse |
| `amyAny.csv` | any amyloid | 0 / 1 | Diverse |
| `reag.csv` | NIA-Reagan likelihood | 0–3 | Diverse |
| `PlaqueTotal.csv` | plaque load | 0–15 | Banner |
| `TangleTotal.csv` | tangle load | 0–15 | Banner |

Subjects are not stacked 14×: the protein values of a subject are identical across files; each file only adds
its own variable. The pipeline rebuilds one row per subject from these files (`DAG1/methodology/data.py`).

Full harmonisation details (source columns, crosswalks, missing-data handling): `README_protein_clinical_csvs.md`.

Individual-level data derived from AMP-AD (Synapse); project members only, under the AMP-AD data-use terms.
