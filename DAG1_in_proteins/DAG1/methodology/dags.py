"""Orderings -> Bernoulli DAGs, score S(G;A) = mean_{(i,j) in E_G} A_ij, pick the best."""
import numpy as np


def dag_from_ordering(order, A, rng):
    """Edge i->j allowed only if i precedes j; drawn with probability A_ij. Always acyclic."""
    P = len(order)
    pos = np.empty(P, int); pos[order] = np.arange(P)
    forward = pos[:, None] < pos[None, :]
    return (rng.random((P, P)) < A) & forward


def score(G, A):
    e = G.sum()
    return float(A[G].mean()) if e else -np.inf


def is_acyclic(G, order):
    Gp = G[np.ix_(order, order)]
    return not np.tril(Gp).any()     # strictly upper-triangular in the ordering


def build_and_select(orderings, A, rng):
    G = np.stack([dag_from_ordering(o, A, rng) for o in orderings])
    scores = np.array([score(g, A) for g in G])
    ok = all(is_acyclic(g, o) for g, o in zip(G, orderings))
    return G, scores, int(np.argmax(scores)), ok
