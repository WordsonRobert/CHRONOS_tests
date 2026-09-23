# One variable at a time

The mirror image of the DAG states. There, one variable was fixed and everything else varied. Here, **one variable
changes and everything else is held the same**, and we look at how each of the 80 proteins moves.

**How:** for each of the 13 variables (same buckets as the DAG states), in each cohort:
protein ~ bucket of that variable + every other clinical variable measured for that person + sex + APOE e4 + log PMI
(+ 5 cell-type scores). The cohorts are pooled by inverse-variance weighting, and the 80 proteins are FDR-corrected per variable.
- Re-codings are never held constant: amyA = binned amyThal, amyAny = (amyCerad > 0), reag = f(Braak, CERAD).
  So **NIA-Reagan is really "Braak + CERAD together"**, not an independent variable.
- `design_check.csv` shows what was held for each variable and how much the variable can still move once the others
  are fixed. In Banner, Braak and tangles barely move independently (15% of their variation is left), so those estimates are weak.
- Braak 0 has only 6 people, so ignore the Braak 0→1 step.

**Main results (cell types held):**
- Once the other pathology is held fixed, most "disease effects" shrink. Braak goes from 30 proteins changed to 7, and the
  amyloid measures from 17–29 to 2–4. Most of what each variable seemed to do is the shared "how sick is this brain" signal.
- The robust independent effects: **APP, LRP1, CLU and APOE go up with amyloid**, and **IDE and PICALM go down with Braak**.
- **TMED2/9/10 per step: no independent effect** of any single pathology variable (all q > 0.2).
- **Cognition step MMSE 18–24 → 24–30** (mild impairment → normal) at the *same* Braak and CERAD:
  TMED10 +0.48, TMED9 +0.48, TMED2 +0.36 SD. They are the **top 3 of 80 proteins** (q = 0.01, 0.01, 0.06), and the
  result repeats in ROSMAP (TMED10/9 are #1/#2) and Banner (TMED2 is #1). It is post-hoc (one of about 40 bucket steps looked at).

**Files:** `one_at_a_time.py` (runs in ~20 s), `output.txt` (printout), `effects_per_step.csv`, `bucket_contrasts.csv`,
`bucket_contrasts_per_stratum.csv`, `effects_matrix_z.csv`, `design_check.csv`, `marginal_vs_held.csv`, `heatmap_cell.png`.
