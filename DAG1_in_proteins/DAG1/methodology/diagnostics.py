"""Per-state diagnostics: is the sampler concentrated-but-stochastic?"""
import numpy as np
from scipy.stats import kendalltau
from probs import sink_probs


def first_sink_entropy(A):
    """Normalised entropy of the first-step sink distribution (1 = uniform, 0 = deterministic)."""
    p = sink_probs(A)
    return float(-(p * np.log(p + 1e-300)).sum() / np.log(len(p))), float(p.max())


def ordering_kendall(orderings):
    """Mean pairwise Kendall tau between orderings (1 = identical, 0 = unrelated)."""
    P = orderings.shape[1]
    pos = np.empty_like(orderings)
    for m, o in enumerate(orderings):
        pos[m, o] = np.arange(P)
    M = len(pos)
    taus = [kendalltau(pos[a], pos[b])[0] for a in range(M) for b in range(a + 1, M)]
    return float(np.mean(taus))


def dag_density(G):
    P = G.shape[1]
    return float(G.sum((1, 2)).mean() / (P * (P - 1) / 2))
