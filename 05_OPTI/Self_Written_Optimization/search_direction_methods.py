"""Search-direction generators for unconstrained gradient-based optimization.

Each function returns the descent direction ``d`` to use at the current
iterate. The actual step length is determined by a separate line search.

References
----------
- Nocedal & Wright, *Numerical Optimization*, 2nd ed., 2006, §5.2.
- Arora, *Introduction to Optimum Design*, 4th ed., §11.4.
- Rao, *Engineering Optimization: Theory and Practice*, 5th ed., §6.9.
"""
from __future__ import annotations

from typing import Optional

import numpy as np


def steepest_descent(grad: np.ndarray) -> np.ndarray:
    """Steepest-descent direction: ``d = -grad``."""
    return -np.asarray(grad, dtype=float)


def conjugate_gradient(
    grad: np.ndarray,
    grad_prev: Optional[np.ndarray] = None,
    d_prev: Optional[np.ndarray] = None,
    variant: str = "PR+",
    restart_every: Optional[int] = None,
    iter_k: int = 0,
) -> np.ndarray:
    """Nonlinear conjugate-gradient direction.

    Parameters
    ----------
    grad        : current gradient g_k
    grad_prev   : previous gradient g_{k-1}; if None, falls back to steepest descent
    d_prev      : previous direction d_{k-1};   if None, falls back to steepest descent
    variant     : 'FR' (Fletcher-Reeves), 'PR' (Polak-Ribiere),
                  'PR+' (default; beta = max(0, beta_PR), auto-restart on negativity),
                  'HS' (Hestenes-Stiefel)
    restart_every : if set, force a steepest-descent restart every n iterations
    iter_k      : current iteration index (0-based) - used for periodic restart

    Returns
    -------
    d : descent direction
    """
    g = np.asarray(grad, dtype=float)

    # First iteration or forced restart -> steepest descent
    if grad_prev is None or d_prev is None:
        return -g
    if restart_every is not None and iter_k > 0 and iter_k % restart_every == 0:
        return -g

    g_prev = np.asarray(grad_prev, dtype=float)
    d_p = np.asarray(d_prev, dtype=float)

    g_prev_dot = float(g_prev @ g_prev)
    if g_prev_dot < 1e-30:
        return -g

    v = (variant or "PR+").upper()
    if v == "FR":
        beta = float(g @ g) / g_prev_dot
    elif v == "PR":
        beta = float(g @ (g - g_prev)) / g_prev_dot
    elif v == "PR+":
        beta = float(g @ (g - g_prev)) / g_prev_dot
        beta = max(0.0, beta)
    elif v == "HS":
        denom = float(d_p @ (g - g_prev))
        if abs(denom) < 1e-30:
            return -g
        beta = float(g @ (g - g_prev)) / denom
    else:
        raise ValueError(
            f"Unknown CG variant '{variant}'. Choose 'FR', 'PR', 'PR+', or 'HS'."
        )

    d = -g + beta * d_p

    # Safety: if numerical noise makes d an ascent direction, restart
    if float(g @ d) >= 0.0:
        return -g

    return d
