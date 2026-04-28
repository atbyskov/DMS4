"""Line-search methods for 1D minimization of phi(alpha) = f(x + alpha * d).

All routines are stateless and operate on a callable ``phi`` that maps a scalar
step length ``alpha`` to the (already projected) merit/objective value at
``x + alpha * d``. They return the chosen step length ``alpha``.

References
----------
- Nocedal & Wright, *Numerical Optimization*, 2nd ed., 2006.
- Arora, *Introduction to Optimum Design*, 4th ed., Chapter 10.
- Rao, *Engineering Optimization: Theory and Practice*, 5th ed., Chapter 5.
"""
from __future__ import annotations

from typing import Callable, Tuple

import numpy as np

PhiFn = Callable[[float], float]

# Golden ratio constants
_GR = (np.sqrt(5.0) - 1.0) / 2.0          # ~0.6180339887
_GR2 = 1.0 - _GR                           # ~0.3819660113


# ----------------------------------------------------------------------
# Bracketing helper
# ----------------------------------------------------------------------
def bracket_minimum(
    phi: PhiFn,
    alpha0: float = 1.0,
    growth: float = 2.0,
    max_iter: int = 50,
) -> Tuple[float, float]:
    """Find ``[a, b]`` such that ``phi`` has a minimum inside.

    Starts at ``alpha = 0`` and steps forward, doubling the step until
    ``phi`` stops decreasing. Returns ``(a, b)`` with ``a < b``.
    """
    a = 0.0
    fa = phi(a)
    b = float(alpha0)
    fb = phi(b)

    if fb >= fa:
        # phi already increases at alpha0; shrink instead
        for _ in range(max_iter):
            b *= 0.5
            fb = phi(b)
            if b < 1e-16:
                return 0.0, float(alpha0)
            if fb < fa:
                break
        else:
            return 0.0, float(alpha0)

    # Now fb < fa, expand forward
    a_prev, fa_prev = a, fa
    a_cur, fa_cur = b, fb
    step = b
    for _ in range(max_iter):
        step *= growth
        a_next = a_cur + step
        f_next = phi(a_next)
        if f_next > fa_cur:
            return a_prev, a_next
        a_prev, fa_prev = a_cur, fa_cur
        a_cur, fa_cur = a_next, f_next

    return a_prev, a_cur


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


# ----------------------------------------------------------------------
# Exact line search: Golden section
# ----------------------------------------------------------------------
def golden_section(
    phi: PhiFn,
    a: float | None = None,
    b: float | None = None,
    tol: float = 1e-5,
    max_iter: int = 100,
    alpha0_bracket: float = 1.0,
) -> float:
    """Golden-section search for the minimum of ``phi`` on ``[a, b]``.

    If ``a`` or ``b`` are not provided, an initial bracket is constructed
    automatically by :func:`bracket_minimum`.
    """
    if a is None or b is None:
        a, b = bracket_minimum(phi, alpha0=alpha0_bracket)
    a = float(a)
    b = float(b)
    if b < a:
        a, b = b, a

    # Two interior points
    x1 = a + _GR2 * (b - a)
    x2 = a + _GR * (b - a)
    f1 = phi(x1)
    f2 = phi(x2)

    for _ in range(max_iter):
        if (b - a) < tol:
            break
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = a + _GR2 * (b - a)
            f1 = phi(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + _GR * (b - a)
            f2 = phi(x2)

    return 0.5 * (a + b)


# ----------------------------------------------------------------------
# Exact line search: Bisection on phi'(alpha)
# ----------------------------------------------------------------------
def bisection(
    phi: PhiFn,
    a: float | None = None,
    b: float | None = None,
    tol: float = 1e-5,
    max_iter: int = 100,
    fd_h: float = 1e-4,
    alpha0_bracket: float = 1.0,
) -> float:
    """Bisection on the derivative of ``phi`` (central FD) over ``[a, b]``.

    Looks for a root of ``phi'(alpha) = 0`` inside ``[a, b]``. The derivative
    is approximated by central finite differences with step ``fd_h``, costing
    two extra evaluations of ``phi`` per bisection step.
    """
    if a is None or b is None:
        a, b = bracket_minimum(phi, alpha0=alpha0_bracket)
    a = float(a)
    b = float(b)
    if b < a:
        a, b = b, a

    def dphi(alpha: float) -> float:
        h = max(fd_h, 1e-12)
        # one-sided near a hard zero so we do not query negative alphas
        if alpha - h < 0.0:
            return (phi(alpha + h) - phi(alpha)) / h
        return (phi(alpha + h) - phi(alpha - h)) / (2.0 * h)

    da = dphi(a)
    db = dphi(b)

    if da > 0.0:
        # Minimum is at (or below) the left bracket
        return a
    if db < 0.0:
        # Minimum is at (or beyond) the right bracket
        return b

    for _ in range(max_iter):
        if (b - a) < tol:
            break
        m = 0.5 * (a + b)
        dm = dphi(m)
        if abs(dm) < tol:
            return m
        if dm > 0.0:
            b = m
            db = dm
        else:
            a = m
            da = dm

    return 0.5 * (a + b)


# ----------------------------------------------------------------------
# Polynomial interpolation (safeguarded quadratic)
# ----------------------------------------------------------------------
def quadratic_interpolation(
    phi: PhiFn,
    phi0: float,
    dphi0: float,
    alpha0: float = 1.0,
    c1: float = 1e-4,
    max_iter: int = 20,
    alpha_min: float = 1e-10,
    shrink: float = 0.5,
) -> float:
    """Safeguarded quadratic-interpolation line search.

    Fits a quadratic through ``phi(0)``, ``phi'(0)``, and ``phi(alpha)`` and
    jumps to its minimizer. Falls back to plain backtracking with the Armijo
    sufficient-decrease test.

    Reference: Nocedal & Wright §3.5 ("Interpolation").
    """
    if dphi0 >= 0.0:
        return max(alpha0 * (shrink ** max_iter), alpha_min)

    alpha = float(alpha0)
    f_alpha = phi(alpha)

    for _ in range(max_iter):
        # Armijo acceptance
        if f_alpha <= phi0 + c1 * alpha * dphi0:
            return alpha

        # Quadratic interpolation: minimize the quadratic m(a) defined by
        #   m(0)=phi0, m'(0)=dphi0, m(alpha)=f_alpha
        denom = 2.0 * (f_alpha - phi0 - dphi0 * alpha)
        if abs(denom) < 1e-20:
            alpha_new = alpha * shrink
        else:
            alpha_new = -dphi0 * alpha * alpha / denom

        # Safeguard: keep alpha_new in [shrink*alpha, (1-shrink)*alpha]
        lo = shrink * alpha * shrink     # ~0.25*alpha for shrink=0.5
        hi = (1.0 - shrink) * alpha       # ~0.5*alpha for shrink=0.5
        if not np.isfinite(alpha_new) or alpha_new < lo or alpha_new > hi:
            alpha_new = alpha * shrink

        alpha = max(alpha_new, alpha_min)
        if alpha <= alpha_min:
            return alpha
        f_alpha = phi(alpha)

    return alpha
