"""Unified driver for the self-written optimization schemes.

Public entry point
------------------
``minimize(obj, con, x0, xl, xu, method=..., ...)`` dispatches to one of:

- ``method='steepest_descent'`` : exterior quadratic penalty + Armijo line search
- ``method='conjugate_gradient'`` : same, with nonlinear CG direction (FR/PR/PR+/HS)
- ``method='slp'``               : Sequential Linear Programming (see ``slp.py``)

The line-search rule is fixed to Armijo backtracking
(:func:`line_search_methods.armijo_backtracking`).

Problem format (PySLSQP convention)
-----------------------------------
    min   f(x) = obj(x)
    s.t.  c(x) >= 0          (callable returning a 1D array)
          xl <= x <= xu

Constraint handling for SD/CG (MVP)
-----------------------------------
Exterior quadratic penalty on inequality violations and clipping for the box:

    P(x; mu) = f(x) + mu * sum_i ( max(0, -c_i(x)) )**2
    x_{k+1} = clip(x_k + alpha * d_k, xl, xu)

References
----------
- Nocedal & Wright, *Numerical Optimization*, 2nd ed., 2006 -- Ch. 3, 5, 17.
- Arora, *Introduction to Optimum Design*, 4th ed. -- Ch. 10, 11.
- Rao, *Engineering Optimization*, 5th ed. -- Ch. 5, 6, 7.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

import numpy as np

from .line_search_methods import armijo_backtracking
from .search_direction_methods import conjugate_gradient, steepest_descent
from .misc.slp_MVP_advanced import solve_slp_mvp, SLPResult

ArrayLike = Union[Sequence[float], np.ndarray]
ObjectiveFn = Callable[[np.ndarray], float]
ConstraintFn = Callable[[np.ndarray], np.ndarray]


# ----------------------------------------------------------------------
# Result type
# ----------------------------------------------------------------------
@dataclass
class Result:
    """Common return type for :func:`minimize`."""

    x: np.ndarray
    fun: float
    nit: int = 0
    nfev: int = 0
    success: bool = False
    message: str = ""
    method: str = ""
    history: List[dict] = field(default_factory=list)


# ----------------------------------------------------------------------
# Helpers
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


def _violation(c: np.ndarray) -> np.ndarray:
    """Inequality violation for c(x) >= 0 convention."""
    if c.size == 0:
        return c
    return np.maximum(0.0, -c)


def _penalty_value(f: float, c: np.ndarray, mu: float) -> float:
    if c.size == 0:
        return float(f)
    v = _violation(c)
    return float(f) + float(mu) * float(np.sum(v * v))


# ----------------------------------------------------------------------
# Public entry point
# ----------------------------------------------------------------------
def minimize(
    obj: ObjectiveFn,
    con: Optional[ConstraintFn],
    x0: ArrayLike,
    xl: Optional[ArrayLike] = None,
    xu: Optional[ArrayLike] = None,
    *,
    method: str = "steepest_descent",
    fd_step=None,
    fd_type: str = "forward",
    maxiter: int = 40,
    gtol: float = 1e-6,
    ftol: float = 1e-3,
    xtol: float = 1e-8,
    penalty_weight: float = 1e3,
    cg_variant: str = "PR+",
    cg_restart_every: Optional[int] = None,
    ls_options: Optional[Dict[str, Any]] = None,
    callback: Optional[Callable[[np.ndarray], None]] = None,
    display: bool = True,
    # Forwarded to solve_slp when method='slp'
    slp_options: Optional[Dict[str, Any]] = None,
) -> Result:
    """Minimize ``obj(x)`` subject to ``con(x) >= 0`` and ``xl <= x <= xu``.

    See module docstring for the full convention.
    """
    method_l = (method or "").lower()
    ls_options = ls_options or {}

    # ------------------------------------------------------------------
    # SLP branch -- delegate and wrap
    # ------------------------------------------------------------------
    if method_l == "slp":
        slp_kwargs = dict(
            fd_step=fd_step,
            fd_type=fd_type,
            maxiter=maxiter,
            gtol=gtol,
            ftol=ftol,
            xtol=xtol,
            infeasibility_penalty=penalty_weight,
            display=display,
        )
        if slp_options:
            slp_kwargs.update(slp_options)
        slp_res: SLPResult = solve_slp_mvp(obj, con, x0, xl, xu, **slp_kwargs)
        return Result(
            x=slp_res.x,
            fun=slp_res.fun,
            nit=slp_res.nit,
            nfev=slp_res.nfev,
            success=slp_res.success,
            message=slp_res.message,
            method="slp",
            history=slp_res.history,
        )
    elif method_l == "slp_mvp":
        slp_kwargs = dict(
            fd_step=fd_step,
            fd_type=fd_type,
            maxiter=maxiter,
            gtol=gtol,
            ftol=ftol,
            xtol=xtol,
            infeasibility_penalty=penalty_weight,
            display=display,
            callback=callback,
        )
        if slp_options:
            slp_kwargs.update(slp_options)
        slp_res: SLPResult = solve_slp_mvp(obj, con, x0, xl, xu, **slp_kwargs)
        return Result(
            x=slp_res.x,
            fun=slp_res.fun,
            nit=slp_res.nit,
            nfev=slp_res.nfev,
            success=slp_res.success,
            message=slp_res.message,
            method="slp_mvp",
            history=slp_res.history,
        )

    if method_l not in {"steepest_descent", "conjugate_gradient"}:
        raise ValueError(
            f"Unknown method '{method}'. Choose 'steepest_descent', "
            f"'conjugate_gradient', 'slp', or 'slp_mvp'."
        )

    # ------------------------------------------------------------------
    # Common setup for SD / CG
    # ------------------------------------------------------------------
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
    mu = float(penalty_weight)

    # Counter for true f/c evaluations made by this driver. The underlying
    # ``evaluate_model`` cache de-duplicates repeated calls at the same x, so
    # this is an upper bound on FEA cost.
    nfev = [0]

    def eval_fc(xv: np.ndarray) -> tuple[float, np.ndarray]:
        nfev[0] += 1
        f_val = float(obj(xv))
        if con is None:
            return f_val, np.empty(0, dtype=float)
        c_val = np.asarray(con(xv), dtype=float).ravel()
        return f_val, c_val

    def P_at(xv: np.ndarray) -> tuple[float, float, np.ndarray]:
        f_val, c_val = eval_fc(xv)
        return _penalty_value(f_val, c_val, mu), f_val, c_val

    def grad_P_fd(xv: np.ndarray, P0: float) -> np.ndarray:
        g = np.zeros(n, dtype=float)
        if fd_type == "central":
            for i in range(n):
                xp = xv.copy(); xp[i] += h[i]
                xm = xv.copy(); xm[i] -= h[i]
                Pp, _, _ = P_at(xp)
                Pm, _, _ = P_at(xm)
                g[i] = (Pp - Pm) / (2.0 * h[i])
        elif fd_type == "backward":
            for i in range(n):
                xm = xv.copy(); xm[i] -= h[i]
                Pm, _, _ = P_at(xm)
                g[i] = (P0 - Pm) / h[i]
        else:  # forward
            for i in range(n):
                xp = xv.copy(); xp[i] += h[i]
                Pp, _, _ = P_at(xp)
                g[i] = (Pp - P0) / h[i]
        return g

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------
    P, f, c = P_at(x)
    history: List[dict] = [{
        "iter": 0,
        "f": float(f),
        "viol": float(np.sum(_violation(c))),
        "P": float(P),
    }]

    if display:
        print("*" * 80)
        print(f"    Custom optimizer: method={method_l}, line_search=armijo, "
              f"penalty={mu:g}")
        print("*" * 80)
        print(f"  {'iter':>5s}  {'f':>12s}  {'viol':>12s}  {'P':>12s}  "
              f"{'||g||':>12s}  {'alpha':>10s}  {'nfev':>6s}")
        print(f"  {0:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  "
              f"{P:12.4e}  {'-':>12s}  {'-':>10s}  {nfev[0]:>6d}")

    grad_prev: Optional[np.ndarray] = None
    d_prev: Optional[np.ndarray] = None
    grad_norm = 0.0
    alpha = 0.0
    success = False
    message = "Maximum iterations reached"
    iter_no = 0

    for iter_no in range(1, maxiter + 1):
        # Gradient of penalized objective by FD
        grad = grad_P_fd(x, P)
        grad_norm = float(np.linalg.norm(grad))

        # Convergence on gradient
        if grad_norm <= gtol:
            success = True
            message = f"Converged: ||grad P|| <= gtol ({gtol:.2e})"
            break

        # Search direction
        if method_l == "steepest_descent":
            d = steepest_descent(grad)
        else:  # conjugate_gradient
            d = conjugate_gradient(
                grad, grad_prev, d_prev,
                variant=cg_variant,
                restart_every=cg_restart_every,
                iter_k=iter_no - 1,
            )

        dphi0 = float(grad @ d)
        if dphi0 >= 0.0:
            # Not a descent direction (numerical) -> restart with -grad
            d = -grad
            dphi0 = float(grad @ d)

        # Cap initial step so the first probe stays inside the box for at
        # least one variable -- avoids wasting an FEA on a fully clipped point.
        d_norm = float(np.linalg.norm(d))
        if d_norm > 0.0:
            # Largest alpha that does not push any *active* component beyond
            # its bound (used only as a sanity cap for alpha0).
            with np.errstate(divide="ignore", invalid="ignore"):
                alpha_to_upper = np.where(d > 0, (xU - x) / d, np.inf)
                alpha_to_lower = np.where(d < 0, (xL - x) / d, np.inf)
            alpha_box = float(np.min(np.concatenate([alpha_to_upper, alpha_to_lower])))
            if not np.isfinite(alpha_box) or alpha_box <= 0.0:
                alpha_box = 1.0
        else:
            alpha_box = 1.0

        alpha0 = float(ls_options.get("alpha0", min(1.0, alpha_box)))

        # phi(alpha) clips before evaluating, so we never step outside the box
        def phi(a: float) -> float:
            x_try = np.clip(x + a * d, xL, xU)
            P_try, _, _ = P_at(x_try)
            return P_try

        alpha = armijo_backtracking(
            phi, P, dphi0,
            alpha0=alpha0,
            c1=ls_options.get("c1", 1e-4),
            rho=ls_options.get("rho", 0.5),
            max_iter=ls_options.get("max_iter", 50),
        )

        x_new = np.clip(x + alpha * d, xL, xU)
        P_new, f_new, c_new = P_at(x_new)

        deltax = x_new - x
        deltanorm = float(np.linalg.norm(deltax))
        deltaf = abs(f_new - f)
        deltaP = abs(P_new - P)

        # Update CG memory regardless of acceptance
        grad_prev = grad
        d_prev = d

        x = x_new
        P = P_new
        f = f_new
        c = c_new

        history.append({
            "iter": iter_no,
            "f": float(f),
            "viol": float(np.sum(_violation(c))),
            "P": float(P),
            "alpha": float(alpha),
            "step": float(deltanorm),
            "grad_norm": grad_norm,
        })

        if display:
            print(f"  {iter_no:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  "
                  f"{P:12.4e}  {grad_norm:12.4e}  {alpha:10.3e}  {nfev[0]:>6d}")

        if callback is not None:
            try:
                callback(x.copy())
            except Exception:
                pass

        if deltanorm <= xtol:
            success = True
            message = f"Converged: ||dx|| <= xtol ({xtol:.2e})"
            break
        if deltaf <= ftol and iter_no > 1:
            success = True
            message = f"Converged: |df| <= ftol ({ftol:.2e})"
            break
        if deltaP <= ftol * max(1.0, abs(P)) and iter_no > 1:
            # Penalized-merit stagnation as a backup criterion
            success = True
            message = f"Converged: |dP| <= ftol*max(1,|P|) ({ftol:.2e})"
            break

    if display:
        viol = float(np.sum(_violation(c))) if c.size else 0.0
        print("-" * 80)
        print(f"  Done: success={success}  iters={iter_no}  nfev={nfev[0]}  "
              f"f={f:.6e}  viol={viol:.3e}")
        print(f"  Message: {message}")

    return Result(
        x=x,
        fun=float(f),
        nit=iter_no,
        nfev=nfev[0],
        success=success,
        message=message,
        method=method_l,
        history=history,
    )


# ----------------------------------------------------------------------
# Smoke test
# ----------------------------------------------------------------------
if __name__ == "__main__":
    def f(x):
        return (x[0] - 1.0) ** 2 + (x[1] - 2.0) ** 2

    def c(x):
        return np.array([
            4.0 - x[0] ** 2 - x[1] ** 2,
            x[0] + x[1] - 1.0,
        ])

    for m in ["steepest_descent", "conjugate_gradient", "slp"]:
        res = minimize(
            f, c,
            x0=np.array([0.5, 0.5]),
            xl=np.array([-5.0, -5.0]),
            xu=np.array([5.0, 5.0]),
            method=m,
            maxiter=40,
            fd_step=1e-5,
            penalty_weight=1e3,
        )
        print(f"\n{m}: x*={res.x}, f*={res.fun:.4e}, nfev={res.nfev}")
