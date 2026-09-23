"""Build the aligned per-subject table from the 14 CSVs and define trajectory states.

One row per subject (key = dataset + individualID). Nothing is stacked 14x:
the protein levels (expr_z) are identical across the CSVs for a given subject,
each CSV only contributes its own trajectory variable.
"""
import numpy as np
import pandas as pd
from config import CSV_DIR, STATE_VARS, Z_VARS


def build_aligned(csv_dir=CSV_DIR):
    X_parts, Z_cols = [], {}
    for v in Z_VARS:
        d = pd.read_csv(csv_dir / f"{v}.csv", dtype={"individualID": str})
        d["key"] = d["dataset"] + "|" + d["individualID"]
        X_parts.append(d[["key", "protein", "expr_z"]])
        Z_cols[v] = d.drop_duplicates("key").set_index("key")[v]
    X_long = pd.concat(X_parts).drop_duplicates(["key", "protein"])
    X = X_long.pivot(index="key", columns="protein", values="expr_z")
    Z = pd.DataFrame(Z_cols).reindex(X.index)
    meta = pd.DataFrame({"dataset": X.index.str.split("|").str[0],
                         "individualID": X.index.str.split("|", n=1).str[1]}, index=X.index)
    return X, Z, meta


def z_numeric(Z):
    """Numeric encoding of Z for use as regression context."""
    Zn = Z.copy()
    Zn["sex"] = Z["sex"].map({"male": 1.0, "female": 0.0})
    Zn["apoe_genotype"] = Z["apoe_genotype"].map(
        lambda g: float(str(int(g)).count("4")) if pd.notna(g) else np.nan)   # e4 dosage 0/1/2
    return Zn.apply(pd.to_numeric, errors="coerce")


def state_labels(Z):
    """Return {var: Series(label per subject, NaN if subject has no value)}."""
    out = {}
    for v, spec in STATE_VARS.items():
        s = Z[v]
        if spec is None:
            if v == "apoe_genotype":
                lab = s.map(lambda g: str(int(g)) if pd.notna(g) else np.nan)
            elif v == "sex":
                lab = s
            else:
                lab = s.map(lambda x: f"{x:g}" if pd.notna(x) else np.nan)
        else:
            w, M = spec
            x = pd.to_numeric(s, errors="coerce")
            k = np.floor(x / w).clip(upper=M // w - 1)
            lab = k.map(lambda b: f"{int(b*w)}-{int(b*w+w)}" if pd.notna(b) else np.nan)
        out[v] = lab
    return out


def list_states(Z):
    """All (var, label, subject_keys) triples, in a stable order."""
    labs = state_labels(Z)
    states = []
    for v in STATE_VARS:
        lab = labs[v]
        for L in sorted(lab.dropna().unique(), key=_sortkey):
            states.append((v, L, lab.index[lab == L]))
    return states


def _sortkey(L):
    try:
        return (0, float(str(L).split("-")[0]))
    except ValueError:
        return (1, str(L))
