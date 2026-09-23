"""Recursive sink sampling with a REFIT at every step -> one ordering; repeat for N orderings.

Refits for an identical remaining protein set are identical (same subjects, same proteins),
so they are cached and stored once; each ordering keeps a pointer to its refit at every step.
"""
import numpy as np
from probs import edge_probs, sink_probs


def sample_orderings(system, n_orderings, rng):
    P = len(system.proteins)
    cache = {}                     # frozen mask bytes -> refit id
    refit_masks, refit_B = [], []  # unique remaining-sets and their B matrices
    orderings = np.zeros((n_orderings, P), dtype=np.int16)      # source-first protein indices
    step_refit = np.zeros((n_orderings, P), dtype=np.int32)     # refit id used at removal step s
    for m in range(n_orderings):
        remaining = list(range(P))
        sinks = []
        for s in range(P):
            mask = np.zeros(P, bool); mask[remaining] = True
            key = mask.tobytes()
            if key not in cache:
                cache[key] = len(refit_B)
                refit_masks.append(mask)
                refit_B.append(system.fit(remaining).astype(np.float32))
            rid = cache[key]
            step_refit[m, s] = rid
            p = sink_probs(edge_probs(refit_B[rid].astype(float)))
            pick = rng.choice(len(remaining), p=p)
            sinks.append(remaining.pop(pick))
        orderings[m] = sinks[::-1]          # reverse: first = source-most, last = final sink
    return orderings, step_refit, refit_masks, refit_B
