"""Sequential Linear Programming (SLP) solver.

Inspired by ``fminslp.py``. Trimmed for simplicity: keeps the adaptive
move-limit trust region and slack-relaxed inequality LP subproblem, drops
the global convergence filter and the augmented-Lagrangian variant.

Problem format (PySLSQP convention)
-----------------------------------
    min   f(x)
    s.t.  c(x) >= 0      (inequalities, callable returns 1D array)
          xl <= x <= xu

References
----------
- Haftka & Gurdal, *Elements of Structural Optimization*, 3rd ed., Ch. 9.
- Rao, *Engineering Optimization*, 5th ed., Ch. 7.
- Companion ``fminslp.py`` in this repository.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Union

import numpy as np
from scipy.optimize import linprog

ArrayLike = Union[Sequence[float], np.ndarray]
ObjectiveFn = Callable[[np.ndarray], float]
ConstraintFn = Callable[[np.ndarray], np.ndarray]


@dataclass
class SLPResult:
    """Return type for :func:`solve_slp` (also reused by ``optimization_methods.minimize``)."""

    x: np.ndarray
    fun: float
    nit: int = 0
    nfev: int = 0
    success: bool = False
    message: str = ""
    history: List[dict] = field(default_factory=list)


# ----------------------------------------------------------------------
# Finite-difference gradients
# ----------------------------------------------------------------------
def _fd_step_array(fd_step, n: int) -> np.ndarray:
    if fd_step is None:
        return np.full(n, np.sqrt(np.finfo(float).eps), dtype=float)
    arr = np.asarray(fd_step, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr), dtype=float)
    if arr.size != n:
        raise ValueError(f"fd_step must be a scalar or have length {n}.")
    return arr.astype(float).ravel()


def _grad_obj_fd(
    obj: ObjectiveFn,
    x: np.ndarray,
    f0: float,
    h: np.ndarray,
    fd_type: str = "forward",
) -> np.ndarray:
    n = x.size
    g = np.zeros(n, dtype=float)
    if fd_type == "central":
        for i in range(n):
            xp = x.copy(); xp[i] += h[i]
            xm = x.copy(); xm[i] -= h[i]
            g[i] = (float(obj(xp)) - float(obj(xm))) / (2.0 * h[i])
    elif fd_type == "backward":
        for i in range(n):
            xm = x.copy(); xm[i] -= h[i]
            g[i] = (f0 - float(obj(xm))) / h[i]
    else:  # forward
        for i in range(n):
            xp = x.copy(); xp[i] += h[i]
            g[i] = (float(obj(xp)) - f0) / h[i]
    return g


def _jac_con_fd(
    con: ConstraintFn,
    x: np.ndarray,
    c0: np.ndarray,
    h: np.ndarray,
    fd_type: str = "forward",
) -> np.ndarray:
    """Return Jacobian J of shape (m, n): J[i, j] = dc_i / dx_j."""
    n = x.size
    m = c0.size
    J = np.zeros((m, n), dtype=float)
    if m == 0:
        return J
    if fd_type == "central":
        for j in range(n):
            xp = x.copy(); xp[j] += h[j]
            xm = x.copy(); xm[j] -= h[j]
            J[:, j] = (np.asarray(con(xp), dtype=float).ravel()
                       - np.asarray(con(xm), dtype=float).ravel()) / (2.0 * h[j])
    elif fd_type == "backward":
        for j in range(n):
            xm = x.copy(); xm[j] -= h[j]
            J[:, j] = (c0 - np.asarray(con(xm), dtype=float).ravel()) / h[j]
    else:  # forward
        for j in range(n):
            xp = x.copy(); xp[j] += h[j]
            J[:, j] = (np.asarray(con(xp), dtype=float).ravel() - c0) / h[j]
    return J


# ----------------------------------------------------------------------
# Adaptive move-limit (trimmed port of fminslp.adaptive_move_limit)
# ----------------------------------------------------------------------
def _adaptive_move_limit(
    x: np.ndarray,
    xLcur: np.ndarray,
    xUcur: np.ndarray,
    xLorg: np.ndarray,
    xUorg: np.ndarray,
    move_limit: float,
    reduce_fac: float,
    expand_fac: float,
    xold1: np.ndarray,
    xold2: np.ndarray,
    reduce_switch: bool,
) -> tuple[np.ndarray, np.ndarray]:
    xLcur = xLcur.copy()
    xUcur = xUcur.copy()

    if reduce_switch:
        expand = reduce_fac
        reduction = reduce_fac
    else:
        reduction = reduce_fac
        expand = expand_fac

    n = x.size
    for i in range(n):
        delta = 0.5 * (xUcur[i] - xLcur[i])
        if abs(x[i] - xold1[i]) > 1.0e-10:
            s1 = (xold1[i] - xold2[i]) / (x[i] - xold1[i])
            if s1 < 0.0:
                delta *= reduction
            else:
                delta *= expand
        else:
            delta *= move_limit

        dmax = (xUorg[i] - xLorg[i]) * move_limit
        if np.isfinite(dmax):
            delta = min(delta, dmax)

        xLcur[i] = max(x[i] - delta, xLorg[i])
        xUcur[i] = min(x[i] + delta, xUorg[i])

        if xLcur[i] > xUcur[i]:
            xLcur[i] = (1.0 - 1.0e-6) * xUcur[i]

    return xLcur, xUcur


# ----------------------------------------------------------------------
# Public entry point
# ----------------------------------------------------------------------
def solve_slp(
    obj: ObjectiveFn,
    con: Optional[ConstraintFn],
    x0: ArrayLike,
    xl: Optional[ArrayLike] = None,
    xu: Optional[ArrayLike] = None,
    *,
    fd_step=None,
    fd_type: str = "forward",
    maxiter: int = 40,
    move_limit: float = 0.1,
    move_limit_reduce: float = 0.5,
    move_limit_expand: float = 1.1,
    infeasibility_penalty: float = 1e3,
    gtol: float = 1e-6,
    ftol: float = 1e-3,
    xtol: float = 1e-8,
    backtrack_max: int = 10,
    display: bool = False,
) -> SLPResult:
    """Sequential Linear Programming with slack-relaxed inequalities.

    See module docstring for the convention used for ``con``.
    """
    x = np.asarray(x0, dtype=float).ravel().copy()
    n = x.size

    if xl is None:
        xL = np.full(n, -np.inf, dtype=float)
    else:
        xL = np.asarray(xl, dtype=float).ravel().copy()
    if xu is None:
        xU = np.full(n, np.inf, dtype=float)
    else:
        xU = np.asarray(xu, dtype=float).ravel().copy()

    h = _fd_step_array(fd_step, n)

    # Initial evaluation
    f = float(obj(x))
    if con is not None:
        c = np.asarray(con(x), dtype=float).ravel()
    else:
        c = np.empty(0, dtype=float)
    m = c.size
    nfev = 1

    def merit(f_val: float, c_val: np.ndarray) -> float:
        if c_val.size == 0:
            return f_val
        viol = np.maximum(0.0, -c_val).sum()
        return f_val + infeasibility_penalty * viol

    P = merit(f, c)
    P_prev = P
    f_prev = f

    xLcur = xL.copy()
    xUcur = xU.copy()
    xold1 = x.copy()
    xold2 = x.copy()

    history: List[dict] = [{"iter": 0, "f": f, "viol": float(np.maximum(0.0, -c).sum() if m else 0.0), "P": P}]

    if display:
        print("*" * 80)
        print("    Custom SLP (slack-relaxed LP, adaptive move limits)")
        print("*" * 80)
        print(f"  {'iter':>5s}  {'f':>12s}  {'viol':>12s}  {'P':>12s}  {'||dx||':>12s}  {'nfev':>6s}")

    success = False
    message = "Maximum iterations reached"
    grad_norm = 0.0
    iter_no = 0

    for iter_no in range(1, maxiter + 1):
        # Gradients
        g = _grad_obj_fd(obj, x, f, h, fd_type=fd_type)
        nfev += n if fd_type != "central" else 2 * n
        if con is not None and m > 0:
            J = _jac_con_fd(con, x, c, h, fd_type=fd_type)
            nfev += n if fd_type != "central" else 2 * n
        else:
            J = np.zeros((0, n), dtype=float)

        grad_norm = float(np.linalg.norm(g))

        # Adaptive move limits
        reduce_switch = False
        xLcur, xUcur = _adaptive_move_limit(
            x, xLcur, xUcur, xL, xU,
            move_limit, move_limit_reduce, move_limit_expand,
            xold1, xold2, reduce_switch,
        )

        # Backtracking on the trust region until the merit decreases
        accepted = False
        x_new = x.copy()
        f_new = f
        c_new = c.copy()
        deltax = np.zeros(n, dtype=float)

        for _bt in range(backtrack_max):
            # LP subproblem variables: z = [dx (n), s (m)]
            #   min   g . dx + M . sum(s)
            #   s.t.  -J dx - s <= c           (linearized c(x_k) + J dx >= -s -> -J dx - s <= c)
            #         s >= 0
            #         xLcur - x_k <= dx <= xUcur - x_k
            n_var = n + m
            c_lp = np.concatenate([g, np.full(m, infeasibility_penalty, dtype=float)])

            if m > 0:
                A_ub = np.zeros((m, n_var), dtype=float)
                A_ub[:, :n] = -J
                A_ub[:, n:] = -np.eye(m)
                b_ub = c.copy()
            else:
                A_ub = None
                b_ub = None

            bounds = [(float(xLcur[j] - x[j]), float(xUcur[j] - x[j])) for j in range(n)] \
                     + [(0.0, None) for _ in range(m)]

            try:
                lp_res = linprog(
                    c=c_lp,
                    A_ub=A_ub, b_ub=b_ub,
                    bounds=bounds,
                    method="highs",
                )
            except Exception as exc:
                message = f"LP solver raised an exception: {exc}"
                break

            if not lp_res.success:
                # Shrink trust region and try again
                xLcur, xUcur = _adaptive_move_limit(
                    x, xLcur, xUcur, xL, xU,
                    move_limit, move_limit_reduce, move_limit_expand,
                    xold1, xold2, reduce_switch=True,
                )
                continue

            z = np.asarray(lp_res.x, dtype=float)
            deltax = z[:n]

            # Trial point (clipped to global bounds for safety)
            x_trial = np.clip(x + deltax, xL, xU)
            f_trial = float(obj(x_trial))
            if con is not None and m > 0:
                c_trial = np.asarray(con(x_trial), dtype=float).ravel()
            else:
                c_trial = np.empty(0, dtype=float)
            nfev += 1

            P_trial = merit(f_trial, c_trial)

            if P_trial < P_prev - 1e-16 * max(1.0, abs(P_prev)):
                # Accept
                x_new = x_trial
                f_new = f_trial
                c_new = c_trial
                accepted = True
                break

            # Reject -> shrink move limits
            xLcur, xUcur = _adaptive_move_limit(
                x, xLcur, xUcur, xL, xU,
                move_limit, move_limit_reduce, move_limit_expand,
                xold1, xold2, reduce_switch=True,
            )

        if not accepted:
            message = "Could not find an acceptable step inside the backtracking budget"
            break

        # Convergence checks
        deltanorm = float(np.linalg.norm(x_new - x))
        deltaf = abs(f_new - f)

        xold2 = xold1.copy()
        xold1 = x.copy()
        x = x_new
        f = f_new
        c = c_new
        P_prev = merit(f, c)
        history.append({
            "iter": iter_no,
            "f": float(f),
            "viol": float(np.maximum(0.0, -c).sum() if m else 0.0),
            "P": float(P_prev),
            "step": deltanorm,
        })

        if display:
            print(f"  {iter_no:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  "
                  f"{P_prev:12.4e}  {deltanorm:12.4e}  {nfev:>6d}")

        if grad_norm <= gtol:
            success = True
            message = f"Converged: ||grad|| <= gtol ({gtol:.2e})"
            break
        if deltaf <= ftol:
            success = True
            message = f"Converged: |df| <= ftol ({ftol:.2e})"
            break
        if deltanorm <= xtol:
            success = True
            message = f"Converged: ||dx|| <= xtol ({xtol:.2e})"
            break

        f_prev = f

    return SLPResult(
        x=x,
        fun=float(f),
        nit=iter_no,
        nfev=nfev,
        success=success,
        message=message,
        history=history,
    )


# ----------------------------------------------------------------------
# Smoke test (only runs when invoked directly)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # min (x0-1)^2 + (x1-2)^2  s.t.  4 - x0^2 - x1^2 >= 0,  x0 + x1 - 1 >= 0
    def f(x):
        return (x[0] - 1.0) ** 2 + (x[1] - 2.0) ** 2

    def c(x):
        return np.array([
            4.0 - x[0] ** 2 - x[1] ** 2,
            x[0] + x[1] - 1.0,
        ])

    res = solve_slp(
        f, c,
        x0=np.array([0.5, 0.5]),
        xl=np.array([-5.0, -5.0]),
        xu=np.array([5.0, 5.0]),
        maxiter=50,
        move_limit=0.25,
        display=True,
    )
    print("\nx*       =", res.x)
    print("f*       =", res.fun)
    print("success  =", res.success)
    print("message  =", res.message)
