"""Line-search methods for 1D minimization of phi(alpha) = f(x + alpha * d).

Only the inexact Armijo backtracking line search is exposed: it is the
acceptance rule used by both the SD/CG penalized-objective branch in
:mod:`optimization_methods` and (in spirit) by the SLP merit acceptance
in :mod:`misc.slp_MVP_advanced`. Exact line searches (golden section,
bisection, quadratic interpolation) were intentionally removed so that
the project uses a single, consistent line-search rule.

References
----------
- Nocedal & Wright, *Numerical Optimization*, 2nd ed., 2006, Algorithm 3.1.
"""
from __future__ import annotations

from typing import Callable

PhiFn = Callable[[float], float]


# ----------------------------------------------------------------------
# Inexact line search: Armijo backtracking
# ----------------------------------------------------------------------
def armijo_backtracking(
    phi: PhiFn,
    phi0: float,
    dphi0: float,
    alpha0: float = 1.0,
    c1: float = 1e-4,
    rho: float = 0.5,
    max_iter: int = 50,
) -> float:
    """Armijo backtracking line search (Nocedal & Wright, Algorithm 3.1).

    Finds ``alpha`` satisfying the sufficient-decrease (Armijo) condition

        phi(alpha) <= phi0 + c1 * alpha * dphi0

    Parameters
    ----------
    phi    : callable, phi(alpha) -> f value
    phi0   : float,    phi(0)
    dphi0  : float,    phi'(0) = grad(f) . d  (must be < 0 for descent)
    alpha0 : float,    initial trial step
    c1     : float,    sufficient-decrease parameter (typ. 1e-4)
    rho    : float,    backtracking contraction (typ. 0.5)
    max_iter : int,    maximum number of backtracking steps
    """
    if dphi0 >= 0.0:
        # Not a descent direction; fall back to a tiny step
        return max(alpha0 * (rho ** max_iter), 1e-12)

    alpha = float(alpha0)
    for _ in range(max_iter):
        if phi(alpha) <= phi0 + c1 * alpha * dphi0:
            return alpha
        alpha *= rho
        if alpha < 1e-16:
            break
    return max(alpha, 1e-16)
