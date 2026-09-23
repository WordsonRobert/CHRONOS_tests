"""B -> edge probabilities A, sink scores, softmax."""
import numpy as np
import config as cfg


def edge_probs(B):
    """A_ij = P(i -> j), from B[target, source] (so B.T is [source, target]).
    threshold: sigmoid((|B_ji| - LAMBDA)/TAU_A)  -- edge existence from effect magnitude
    signed:    sigmoid(B_ji / TAU_SIGMOID)
    The sign of B is kept separately (B_init.csv, best_dag_edges.csv)."""
    Bt = B.T
    if cfg.A_MAP == "threshold":
        A = 1.0 / (1.0 + np.exp(-(np.abs(Bt) - cfg.LAMBDA) / cfg.TAU_A))
    else:
        A = 1.0 / (1.0 + np.exp(-Bt / cfg.TAU_SIGMOID))
    np.fill_diagonal(A, 0.0)
    return A


def resolve_lambda(B0):
    """In quantile mode, set LAMBDA (and TAU_A) for this state from its initial 80x80 B."""
    if cfg.LAMBDA_MODE == "quantile":
        off = ~np.eye(B0.shape[0], dtype=bool)
        cfg.LAMBDA = float(np.quantile(np.abs(B0[off]), cfg.LAMBDA_Q))
        cfg.TAU_A = cfg.TAU_A_RATIO * cfg.LAMBDA
    return cfg.LAMBDA


def sink_probs(A, t=None):
    """p_sink[i] = 1 - max_{j != i} A_ij ;  p = softmax(p_sink / t)."""
    t = cfg.T_SOFTMAX if t is None else t
    k = A.shape[0]
    if k == 1:
        return np.ones(1)
    Am = A.copy()
    np.fill_diagonal(Am, -np.inf)
    psink = 1.0 - Am.max(1)
    z = psink / t
    z -= z.max()
    p = np.exp(z)
    return p / p.sum()
