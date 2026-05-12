"""
Sequential Linear Programming (SLP) solver implementing the fminslp methodology.

This module is a faithful Python port of the merit-function variant of the
``fminslp`` MATLAB solver written by Rene S. Nielsen.  All algorithmic
ingredients of the reference solver are preserved:

    * Slack variables ``y >= 0`` that relax the linearized inequality
      constraints so that the LP subproblem is always feasible.
    * A quadratic merit objective in the slacks,
      ``f(x) + a_F * sum( y*R + 0.5*y^2 )``,
      which simultaneously punishes infeasibility and reduces to the original
      objective at ``y = 0``.
    * An LP subproblem in the augmented vector ``z = [x, y]`` whose constraints
      are the user-supplied linear constraints together with the slacked
      linearization of the nonlinear constraints.
    * A global convergence filter of the Fletcher--Leyffer--Toint type that
      accepts or rejects each LP trial point by comparing it against a Pareto
      set of past (h, f) pairs.
    * Backtracking driven by both the filter and an Armijo-like sufficient
      decrease test on the merit predicted by the LP.
    * Adaptive move limits that use the sign of the velocity over three
      consecutive iterates to detect oscillation and shrink (or expand) the
      trust region per design variable.
    * First-order, function-value and step-size convergence checks evaluated
      on the merit function in physical units.

The PySLSQP-style convention ``c(x) >= 0`` is exposed at the API boundary; the
solver works internally in the ``g(x) <= 0`` convention by negating ``con(x)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.optimize import linprog


ArrayLike = Union[Sequence[float], np.ndarray]
ObjectiveFn = Callable[[np.ndarray], float]
ConstraintFn = Callable[[np.ndarray], np.ndarray]


# ----------------------------------------------------------------------
# Public result type (signature is unchanged so the rest of the pipeline
# keeps working).
# ----------------------------------------------------------------------
@dataclass
class SLPResult:
    """Return type for :func:`solve_slp_mvp`."""

    x: np.ndarray
    fun: float
    nit: int = 0
    nfev: int = 0
    success: bool = False
    message: str = ""
    history: List[dict] = field(default_factory=list)


# ----------------------------------------------------------------------
# Global convergence filter (Fletcher--Leyffer--Toint style)
# ----------------------------------------------------------------------
@dataclass
class _Filter:
    beta: float
    gamma: float
    sigma: float
    vals: np.ndarray
    n_vals: int
    init_f: float

    def is_accepted(self, h: float, f: float) -> bool:
        for ii in range(self.n_vals):
            hi = self.vals[ii, 0]
            fi = self.vals[ii, 1]
            if (h <= hi * self.beta) or ((f + self.gamma * h) <= fi):
                continue
            return False
        return True

    def add(self, h: float, f: float) -> None:
        kept = []
        for ii in range(self.n_vals):
            hi = self.vals[ii, 0]
            fi = self.vals[ii, 1]
            dominated = (h <= hi) and (f <= fi)
            if not dominated:
                kept.append((hi, fi))
        kept.append((h, f))
        if len(kept) > self.vals.shape[0]:
            new_vals = np.zeros((len(kept) + 32, 2), dtype=float)
            new_vals[: len(kept), :] = np.asarray(kept, dtype=float)
            self.vals = new_vals
        else:
            self.vals[: len(kept), :] = np.asarray(kept, dtype=float)
        self.n_vals = len(kept)


def _init_filter(max_iter: int, max_infeasibility: float) -> _Filter:
    small = 1.0e-6
    vals = np.zeros((max_iter + 4, 2), dtype=float)
    if np.isfinite(max_infeasibility):
        vals[0, 0] = max_infeasibility
        vals[0, 1] = -np.inf
    else:
        vals[0, 0] = np.inf
        vals[0, 1] = np.inf
    return _Filter(
        beta=1.0 - small,
        gamma=small,
        sigma=2.0 * small,
        vals=vals,
        n_vals=1,
        init_f=0.0,
    )


# ----------------------------------------------------------------------
# Bookkeeping helpers
# ----------------------------------------------------------------------
def _fd_step_array(fd_step, n: int) -> np.ndarray:
    if fd_step is None:
        return np.full(n, float(np.sqrt(np.finfo(float).eps)), dtype=float)
    arr = np.asarray(fd_step, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr), dtype=float)
    if arr.size != n:
        raise ValueError(f"fd_step must be a scalar or have length {n}.")
    return np.maximum(np.abs(arr.astype(float).ravel()), 1e-12)


def _as_bounds(x0: np.ndarray, xl, xu) -> Tuple[np.ndarray, np.ndarray]:
    n = x0.size
    xL = np.full(n, -np.inf, dtype=float) if xl is None else np.asarray(xl, dtype=float).ravel()
    xU = np.full(n, np.inf, dtype=float) if xu is None else np.asarray(xu, dtype=float).ravel()
    if xL.size != n or xU.size != n:
        raise ValueError("xl and xu must have the same length as x0.")
    if np.any(xL > xU):
        raise ValueError("Each lower bound must be <= the corresponding upper bound.")
    return xL.copy(), xU.copy()


def _max_violation(g: np.ndarray) -> float:
    """Maximum violation of ``g(x) <= 0`` (positive value means infeasible)."""
    if g.size == 0:
        return 0.0
    return float(np.max(np.maximum(g, 0.0)))


def _sum_violation(g: np.ndarray) -> float:
    if g.size == 0:
        return 0.0
    return float(np.sum(np.maximum(g, 0.0)))


# ----------------------------------------------------------------------
# Finite differences -- obj and constraint Jacobian in a single bound-aware
# loop.  obj and con are queried via wrappers that count true evaluations.
# ----------------------------------------------------------------------
def _fd_perturbations(
    x: np.ndarray,
    i: int,
    h_abs: float,
    xL: np.ndarray,
    xU: np.ndarray,
    fd_type: str,
) -> Tuple[str, Optional[np.ndarray], Optional[np.ndarray], float]:
    h = max(float(abs(h_abs)), 1e-12)
    can_forward = (x[i] + h) <= xU[i] if np.isfinite(xU[i]) else True
    can_backward = (x[i] - h) >= xL[i] if np.isfinite(xL[i]) else True

    if fd_type == "central" and can_forward and can_backward:
        xp = x.copy(); xp[i] += h
        xm = x.copy(); xm[i] -= h
        return "central", xp, xm, 2.0 * h
    if fd_type == "backward" and can_backward:
        xm = x.copy(); xm[i] -= h
        return "backward", None, xm, h
    if can_forward:
        xp = x.copy(); xp[i] += h
        return "forward", xp, None, h
    if can_backward:
        xm = x.copy(); xm[i] -= h
        return "backward", None, xm, h
    return "fixed", None, None, 1.0


def _grad_and_jac_fd(
    obj_eval: Callable[[np.ndarray], Tuple[float, np.ndarray]],
    x: np.ndarray,
    f0: float,
    g0: np.ndarray,
    h: np.ndarray,
    xL: np.ndarray,
    xU: np.ndarray,
    fd_type: str,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """Finite-difference gradient and constraint Jacobian using a single sweep.

    ``obj_eval(x)`` must return ``(f, g)`` where ``g`` is in the ``g <= 0``
    convention.  Returns ``(df, dg, n_evals)`` with shapes ``(n,)`` and
    ``(m, n)``.
    """
    n = x.size
    m = g0.size
    df = np.zeros(n, dtype=float)
    dg = np.zeros((m, n), dtype=float)
    n_evals = 0
    for i in range(n):
        mode, xp, xm, denom = _fd_perturbations(x, i, h[i], xL, xU, fd_type)
        if mode == "central":
            fp, gp = obj_eval(xp)
            fm, gm = obj_eval(xm)
            df[i] = (fp - fm) / denom
            if m:
                dg[:, i] = (gp - gm) / denom
            n_evals += 2
        elif mode == "backward":
            fm, gm = obj_eval(xm)
            df[i] = (f0 - fm) / denom
            if m:
                dg[:, i] = (g0 - gm) / denom
            n_evals += 1
        elif mode == "forward":
            fp, gp = obj_eval(xp)
            df[i] = (fp - f0) / denom
            if m:
                dg[:, i] = (gp - g0) / denom
            n_evals += 1
        else:
            df[i] = 0.0
            if m:
                dg[:, i] = 0.0
    return df, dg, n_evals


# ----------------------------------------------------------------------
# Adaptive move limits -- direct port of fminslp's adaptive_move_limit
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
    min_dv_box_limit: float,
) -> Tuple[np.ndarray, np.ndarray]:
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

        if np.isfinite(xUorg[i]) and np.isfinite(xLorg[i]):
            dmax = (xUorg[i] - xLorg[i]) * move_limit
            delta = min(delta, dmax)

        xLcur[i] = max(x[i] - delta, xLorg[i])
        xUcur[i] = min(x[i] + delta, xUorg[i])

        if xLcur[i] > xUcur[i]:
            if min_dv_box_limit > 0.0:
                mid = 0.5 * (xLcur[i] + xUcur[i])
                xLcur[i] = mid - 0.5 * min_dv_box_limit
                xUcur[i] = mid + 0.5 * min_dv_box_limit
            else:
                xLcur[i] = (1.0 - 1.0e-6) * xUcur[i]

        width = xUcur[i] - xLcur[i]
        if width < min_dv_box_limit:
            target = min_dv_box_limit - width
            grow_left = min(0.5 * target, max(0.0, xLcur[i] - xLorg[i]))
            grow_right = min(0.5 * target, max(0.0, xUorg[i] - xUcur[i]))
            xLcur[i] -= grow_left
            xUcur[i] += grow_right

            remaining = min_dv_box_limit - (xUcur[i] - xLcur[i])
            if remaining > 0:
                extra_right = min(remaining, max(0.0, xUorg[i] - xUcur[i]))
                xUcur[i] += extra_right
                remaining = min_dv_box_limit - (xUcur[i] - xLcur[i])
            if remaining > 0:
                extra_left = min(remaining, max(0.0, xLcur[i] - xLorg[i]))
                xLcur[i] -= extra_left

            xLcur[i] = max(xLcur[i], xLorg[i])
            xUcur[i] = min(xUcur[i], xUorg[i])

    return xLcur, xUcur


def _enforce_box_floor(
    x: np.ndarray,
    xLcur: np.ndarray,
    xUcur: np.ndarray,
    xLorg: np.ndarray,
    xUorg: np.ndarray,
    min_dv_box_limit: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Re-impose ``min_dv_box_limit`` after manual trust-region shrinks.

    The backtracking branches multiply ``(xUcur - xLcur)`` by ``backtrack_shrink``
    directly, bypassing the floor that ``_adaptive_move_limit`` keeps.  Without
    this re-clip the trust region can collapse to zero and the solver reports
    "no acceptable step and trust region is below move_limit_min" spuriously.
    """
    if min_dv_box_limit <= 0.0:
        return xLcur, xUcur
    xLcur = xLcur.copy()
    xUcur = xUcur.copy()
    for i in range(x.size):
        width = xUcur[i] - xLcur[i]
        if width >= min_dv_box_limit:
            continue
        deficit = min_dv_box_limit - width
        grow_left = min(0.5 * deficit, max(0.0, xLcur[i] - xLorg[i]))
        grow_right = min(0.5 * deficit, max(0.0, xUorg[i] - xUcur[i]))
        xLcur[i] -= grow_left
        xUcur[i] += grow_right
        remaining = min_dv_box_limit - (xUcur[i] - xLcur[i])
        if remaining > 0:
            extra_right = min(remaining, max(0.0, xUorg[i] - xUcur[i]))
            xUcur[i] += extra_right
            remaining = min_dv_box_limit - (xUcur[i] - xLcur[i])
        if remaining > 0:
            extra_left = min(remaining, max(0.0, xLcur[i] - xLorg[i]))
            xLcur[i] -= extra_left
        xLcur[i] = max(xLcur[i], xLorg[i])
        xUcur[i] = min(xUcur[i], xUorg[i])
    return xLcur, xUcur


# ----------------------------------------------------------------------
# LP subproblem -- single call to HiGHS via scipy.linprog.
# ----------------------------------------------------------------------
def _solve_lp(
    c: np.ndarray,
    A_ub: Optional[np.ndarray],
    b_ub: Optional[np.ndarray],
    A_eq: Optional[np.ndarray],
    b_eq: Optional[np.ndarray],
    lb: np.ndarray,
    ub: np.ndarray,
) -> Tuple[np.ndarray, bool, str, int]:
    bounds = [(float(lb[i]), float(ub[i])) for i in range(c.size)]
    A_ub_in = A_ub if (A_ub is not None and A_ub.size > 0) else None
    b_ub_in = b_ub if (b_ub is not None and b_ub.size > 0) else None
    A_eq_in = A_eq if (A_eq is not None and A_eq.size > 0) else None
    b_eq_in = b_eq if (b_eq is not None and b_eq.size > 0) else None

    res = linprog(
        c=c,
        A_ub=A_ub_in,
        b_ub=b_ub_in,
        A_eq=A_eq_in,
        b_eq=b_eq_in,
        bounds=bounds,
        method="highs",
    )
    if res.success:
        return np.asarray(res.x, dtype=float), True, str(res.message), int(res.status)
    return np.zeros_like(c), False, str(res.message), int(res.status)


# ----------------------------------------------------------------------
# Main solver
# ----------------------------------------------------------------------
def solve_slp_mvp(
    obj: ObjectiveFn,
    con: Optional[ConstraintFn],
    x0: ArrayLike,
    xl: Optional[ArrayLike] = None,
    xu: Optional[ArrayLike] = None,
    fd_step=None,
    fd_type: str = "forward",
    maxiter: int = 50,
    move_limit: float = 0.10,
    move_limit_min: float = 1e-5,
    move_limit_expand: float = 1.1,
    move_limit_shrink: float = 0.5,
    infeasibility_penalty: float = 1000.0,
    penalty_increase: float = 10.0,
    penalty_max: float = 1e9,
    feasibility_tol: float = 1e-5,
    slack_tol: float = 1e-8,
    xtol: float = 1e-8,
    gtol: float = 1e-6,
    ftol: float = 1e-6,
    backtrack_max: int = 8,
    backtrack_shrink: float = 0.5,
    sufficient_decrease: float = 2.0e-6,
    feasibility_reduction: float = 1e-3,
    objective_limit: float = -1e20,
    max_infeasibility: float = np.inf,
    algorithm: str = "merit",
    display: bool = False,
    callback: Optional[Callable[[np.ndarray], None]] = None,
    **_unused_kwargs,
) -> SLPResult:
    """Solve ``min f(x) s.t. c(x) >= 0, xl <= x <= xu`` with fminslp methodology.

    The PySLSQP convention is used on the public interface (``c >= 0`` means
    feasible).  Internally the algorithm uses ``g(x) = -c(x) <= 0`` and slack
    variables ``y >= 0`` such that ``g(x) - y <= 0`` is always satisfied for
    sufficiently large ``y``.  The merit function in ``y`` then drives the
    slacks toward zero so that the original problem is recovered at a feasible
    optimum.

    Parameters
    ----------
    obj, con :
        Objective and constraint callables.  ``con(x)`` returns a 1-D array of
        constraint values in the ``c(x) >= 0`` convention or ``None``.
    x0, xl, xu :
        Initial guess and bounds (both bounds optional).
    fd_step :
        Scalar or per-variable absolute step size for finite differences.
    fd_type :
        ``'forward'``, ``'backward'`` or ``'central'``.
    maxiter :
        Maximum number of outer iterations.
    move_limit, move_limit_min, move_limit_expand, move_limit_shrink :
        Adaptive move-limit parameters.
    infeasibility_penalty, penalty_increase, penalty_max :
        Quadratic merit penalty ``R``, its multiplicative growth factor, and
        its upper cap.
    feasibility_tol, slack_tol :
        Maximum acceptable nonlinear violation and LP slack for ``success``.
    xtol, gtol, ftol :
        Step, gradient and merit-change tolerances.
    backtrack_max, backtrack_shrink :
        Inner-iteration safeguards: maximum number of LP re-solves per outer
        iteration, and the move-limit shrink factor that is applied each time
        a trial point is rejected (in addition to the adaptive update).
    sufficient_decrease :
        ``sigma`` parameter of the Armijo-like sufficient-decrease test.
    feasibility_reduction :
        Fraction by which the nonlinear violation must drop to accept an
        infeasible step without merit improvement.
    objective_limit :
        Lower stopping value for the merit (``exitflag = 4``).
    max_infeasibility :
        Initial maximum admissible filter ``h`` value.  ``np.inf`` accepts any
        trial point at iteration 0.
    algorithm :
        ``'merit'`` (quadratic merit, default) or ``'al'`` (augmented
        Lagrangian: replaces the linear slack penalty with multipliers).
    display, callback :
        Iteration printout flag and per-accepted-step callback.
    """
    algorithm = str(algorithm).lower()
    if algorithm not in {"merit", "al"}:
        raise ValueError("algorithm must be 'merit' or 'al'.")

    x = np.asarray(x0, dtype=float).ravel().copy()
    n = x.size
    xL, xU = _as_bounds(x, xl, xu)
    x = np.clip(x, xL, xU)

    h = _fd_step_array(fd_step, n)

    # --- True evaluation counters -----------------------------------------
    nfev = [0]

    def obj_only(xv: np.ndarray) -> float:
        nfev[0] += 1
        return float(obj(np.asarray(xv, dtype=float)))

    def obj_and_g(xv: np.ndarray) -> Tuple[float, np.ndarray]:
        nfev[0] += 1
        f_val = float(obj(np.asarray(xv, dtype=float)))
        if con is None:
            return f_val, np.empty(0, dtype=float)
        c_val = np.asarray(con(np.asarray(xv, dtype=float)), dtype=float).ravel()
        # External convention is c >= 0.  Internal convention is g <= 0.
        return f_val, -c_val

    # --- Initial evaluation and slack initialization ----------------------
    f, g = obj_and_g(x)
    m = g.size

    if m > 0:
        y = np.maximum(g, 0.0)
        ylb = np.zeros(m, dtype=float)
        upper_seed = max(float(np.max(np.abs(g))), 1.0) if g.size else 1.0
        yub = np.full(m, upper_seed * 1.0e6, dtype=float)
    else:
        y = np.empty(0, dtype=float)
        ylb = np.empty(0, dtype=float)
        yub = np.empty(0, dtype=float)

    # Augmented-Lagrangian multipliers (zero for pure-merit variant).
    lambda_vec = np.ones(m, dtype=float) if (algorithm == "al" and m > 0) else np.empty(0, dtype=float)

    a_F = max(abs(f), 1.0)
    R = float(infeasibility_penalty)

    def merit_val(f_v: float, y_v: np.ndarray, lam: np.ndarray) -> float:
        if y_v.size == 0:
            return float(f_v)
        if algorithm == "merit":
            return float(f_v + a_F * float(np.sum(y_v * R + 0.5 * y_v * y_v)))
        # algorithm == 'al'
        return float(f_v + a_F * float(np.sum(y_v * lam + 0.5 * R * y_v * y_v)))

    def merit_grad_y(y_v: np.ndarray, lam: np.ndarray) -> np.ndarray:
        if y_v.size == 0:
            return np.empty(0, dtype=float)
        if algorithm == "merit":
            return a_F * (R + y_v)
        return a_F * (lam + R * y_v)

    Phi = merit_val(f, y, lambda_vec)
    init_F = Phi
    filt = _init_filter(max_iter=maxiter, max_infeasibility=max_infeasibility)
    filt.init_f = init_F

    # --- Trust-region bookkeeping ----------------------------------------
    xLcur = xL.copy()
    xUcur = xU.copy()
    xold1 = x.copy()
    xold2 = x.copy()
    min_dv_box_limit = max(2.0 * xtol, move_limit_min)

    # --- History -----------------------------------------------------------
    history: List[dict] = [{
        "iter": 0,
        "accepted": True,
        "f": float(f),
        "viol": _sum_violation(g),
        "max_viol": _max_violation(g),
        "P": float(Phi),
        "penalty": R,
        "step": 0.0,
        "alpha": 0.0,
        "rho": float("nan"),
        "predicted_decrease": 0.0,
        "actual_decrease": 0.0,
        "slack_sum": float(np.sum(y)) if y.size else 0.0,
        "slack_max": float(np.max(y)) if y.size else 0.0,
        "max_delta_y": float(move_limit),
        "lp_status": 0,
        "lp_message": "initial",
    }]

    if display:
        print("*" * 110)
        print("    SLP optimizer  (fminslp methodology: merit + global convergence filter)")
        print(f"    Algorithm: {algorithm}   Constraints: {m}   Variables: {n}")
        print("*" * 110)
        print(
            f"  {'iter':>5s}  {'f':>12s}  {'viol':>12s}  {'P':>12s}  "
            f"{'||dx||':>10s}  {'||g_y||':>10s}  {'rho':>8s}  "
            f"{'sum(s)':>10s}  {'M':>10s}  {'nfev':>6s}"
        )
        print(
            f"  {0:>5d}  {f:12.4e}  {_sum_violation(g):12.4e}  {Phi:12.4e}  "
            f"{0.0:10.3e}  {0.0:10.3e}  {'-':>8s}  "
            f"{(float(np.sum(y)) if y.size else 0.0):10.3e}  {R:10.3e}  {nfev[0]:>6d}"
        )

    if callback is not None:
        try:
            callback(x.copy())
        except Exception:
            pass

    # ---------------------------------------------------------------------
    # Outer iteration loop
    # ---------------------------------------------------------------------
    iter_no = 0
    success = False
    message = f"Number of iterations exceeded the limit: {maxiter:d}"
    optimality_norm = 0.0
    max_viol = _max_violation(g)
    accepted_steps = 0
    rejected_steps = 0
    converged = False

    for iter_no in range(1, maxiter + 1):
        # 1) Sensitivities at the current iterate (forward / central FD).
        df, dg, n_evals_this = _grad_and_jac_fd(
            obj_and_g, x, f, g, h, xL, xU, fd_type=fd_type,
        )
        # nfev was already incremented inside obj_and_g, no extra bookkeeping
        # needed.  (n_evals_this is just diagnostic.)

        # 2) Build the merit gradient w.r.t. z = [x, y].
        dy = merit_grad_y(y, lambda_vec)
        dfmerit = np.concatenate([df, dy]) if m > 0 else df.copy()
        optimality_norm = float(np.linalg.norm(df))

        # 3) Build the LP inequality system.  Linearized constraint:
        #       g(x_k) + dg @ dx + (-I) @ dy_slack <= 0
        #    expanded in z = [x, y]:
        #       dgmerit @ z <= dgmerit @ z_k - (g(x_k) - y_k)
        if m > 0:
            dgmerit = np.zeros((m, n + m), dtype=float)
            dgmerit[:, :n] = dg
            dgmerit[:, n:] = -np.eye(m)
            z_cur = np.concatenate([x, y])
            gmerit_cur = g - y
            A_ub = dgmerit
            b_ub = dgmerit @ z_cur - gmerit_cur
        else:
            A_ub = None
            b_ub = None

        # 4) Adaptive move limits for a fresh outer iteration.
        xLcur, xUcur = _adaptive_move_limit(
            x, xLcur, xUcur, xL, xU,
            move_limit, move_limit_shrink, move_limit_expand,
            xold1, xold2, reduce_switch=False,
            min_dv_box_limit=min_dv_box_limit,
        )
        xold2 = xold1.copy()
        xold1 = x.copy()

        # 5) Backtracking via filter + sufficient decrease.  Each rejected
        #    trial shrinks the move limits and resolves the LP.
        x_new = x.copy()
        y_new = y.copy()
        g_new = g.copy()
        f_new = f
        Phi_new = Phi
        deltax_norm = 0.0
        actual_decrease = 0.0
        predicted_decrease = 0.0
        rho_val = float("nan")
        add_to_filter = False
        backtrack_iter = 0
        accepted = False
        lp_message = ""
        lp_status = 0
        h_trial = 0.0
        f_trial = 0.0

        viol_now = _sum_violation(g)
        viol_trial = viol_now

        for backtrack_iter in range(1, backtrack_max + 1):
            zlb = np.concatenate([xLcur, ylb])
            zub = np.concatenate([xUcur, yub])

            z_new, lp_ok, lp_message, lp_status = _solve_lp(
                dfmerit, A_ub, b_ub, None, None, zlb, zub,
            )
            if not lp_ok:
                xLcur, xUcur = _adaptive_move_limit(
                    x, xLcur, xUcur, xL, xU,
                    move_limit, move_limit_shrink, move_limit_expand,
                    xold1, xold2, reduce_switch=True,
                    min_dv_box_limit=min_dv_box_limit,
                )
                xLcur = x - (x - xLcur) * backtrack_shrink
                xUcur = x + (xUcur - x) * backtrack_shrink
                xLcur, xUcur = _enforce_box_floor(
                    x, xLcur, xUcur, xL, xU, min_dv_box_limit,
                )
                if float(np.max(xUcur - xLcur)) < min_dv_box_limit:
                    break
                continue

            x_new = z_new[:n]
            y_new = z_new[n:] if m > 0 else np.empty(0, dtype=float)

            deltaz = z_new - np.concatenate([x, y]) if m > 0 else (x_new - x)
            deltax = x_new - x
            deltax_norm = float(np.linalg.norm(deltax))

            # Evaluate true nonlinear functions at the trial point.
            f_new, g_new = obj_and_g(x_new)

            # Merit, predicted decrease, filter pair.  ``Phi`` is the merit
            # at the current accepted iterate, so ``actual_decrease`` uses it
            # directly as the baseline.
            Phi_new = merit_val(f_new, y_new, lambda_vec)
            actual_decrease = float(Phi - Phi_new)
            predicted_decrease = float(-dfmerit @ deltaz)

            if m > 0:
                gmerit_new = g_new - y_new
                h_trial = _max_violation(gmerit_new)
            else:
                h_trial = 0.0
            denom = init_F if abs(init_F) > 1.0e-30 else 1.0
            f_trial = float(Phi_new / denom)

            # Feasibility status used by both convergence checks and the
            # secondary acceptance branch.
            current_feasible = max_viol <= feasibility_tol and (
                (y.size == 0) or (float(np.max(y)) <= max(slack_tol, feasibility_tol))
            )

            # Convergence checks evaluated on the trial (mirror fminslp).
            # Success is only reported when the *current accepted* iterate
            # is feasible, so the structural design contract is honoured.
            if optimality_norm <= gtol and current_feasible:
                converged = True
                success = True
                message = f"Converged: ||grad f|| <= gtol ({gtol:.2e})."
                break
            if abs(actual_decrease) <= ftol and current_feasible and iter_no > 1:
                converged = True
                success = True
                message = f"Converged: |dP| <= ftol ({ftol:.2e}) at feasible point."
                break
            if deltax_norm <= xtol and current_feasible:
                converged = True
                success = True
                message = f"Converged: ||dx|| <= xtol ({xtol:.2e}) at feasible point."
                break
            if Phi_new <= objective_limit and iter_no > 1:
                converged = True
                success = True
                message = f"Stopped: merit <= objective_limit ({objective_limit:.2e})."
                break

            filter_accepts = filt.is_accepted(h_trial, f_trial)
            sufficient = (
                actual_decrease >= sufficient_decrease * predicted_decrease
                or predicted_decrease <= 0.0
            )
            # Auxiliary feasibility-progress branch: accept an infeasible
            # step that reduces the nonlinear violation by more than a fixed
            # fraction (feasibility restoration in the filter-SQP spirit).
            viol_trial = _sum_violation(g_new)
            feasibility_progress = (
                viol_now > feasibility_tol
                and viol_trial <= viol_now * (1.0 - feasibility_reduction)
            )

            if filter_accepts and (sufficient or feasibility_progress):
                rho_val = (
                    actual_decrease / predicted_decrease
                    if predicted_decrease > 1.0e-14
                    else float("nan")
                )
                add_to_filter = h_trial > 0.0
                accepted = True
                break

            # Reject -> shrink trust region and re-solve.
            rejected_steps += 1
            xLcur, xUcur = _adaptive_move_limit(
                x, xLcur, xUcur, xL, xU,
                move_limit, move_limit_shrink, move_limit_expand,
                xold1, xold2, reduce_switch=True,
                min_dv_box_limit=min_dv_box_limit,
            )
            xLcur = x - (x - xLcur) * backtrack_shrink
            xUcur = x + (xUcur - x) * backtrack_shrink
            xLcur, xUcur = _enforce_box_floor(
                x, xLcur, xUcur, xL, xU, min_dv_box_limit,
            )

            if float(np.max(xUcur - xLcur)) < min_dv_box_limit:
                break

        # ------------------------------------------------------------------
        # Process the outcome of the inner search.
        # ------------------------------------------------------------------
        if converged:
            history.append({
                "iter": iter_no,
                "accepted": False,
                "f": float(f),
                "viol": _sum_violation(g),
                "max_viol": _max_violation(g),
                "P": float(Phi),
                "penalty": R,
                "step": deltax_norm,
                "alpha": 0.0,
                "rho": float(rho_val) if np.isfinite(rho_val) else float("nan"),
                "predicted_decrease": float(predicted_decrease),
                "actual_decrease": float(actual_decrease),
                "slack_sum": float(np.sum(y)) if y.size else 0.0,
                "slack_max": float(np.max(y)) if y.size else 0.0,
                "max_delta_y": float(np.max(xUcur - xLcur)),
                "lp_status": int(lp_status),
                "lp_message": lp_message,
                "backtracks": int(backtrack_iter),
                "termination": "converged_on_trial",
            })
            if display:
                print(
                    f"  {iter_no:>5d}  {f:12.4e}  {_sum_violation(g):12.4e}  {Phi:12.4e}  "
                    f"{deltax_norm:10.3e}  {optimality_norm:10.3e}  "
                    f"CONVERGED  {(float(np.sum(y)) if y.size else 0.0):10.3e}  {R:10.3e}  {nfev[0]:>6d}"
                )
            break

        if accepted:
            if add_to_filter:
                filt.add(h_trial, f_trial)

            # Commit the accepted step.
            x = x_new.copy()
            f = f_new
            g = g_new
            y = y_new.copy() if m > 0 else y
            Phi = Phi_new
            accepted_steps += 1
            max_viol = _max_violation(g)

            # AL multiplier update.
            if algorithm == "al" and m > 0:
                lambda_vec = np.maximum(lambda_vec + R * y, 1.0)

            # Adaptive penalty: grow R when the slacks remain active or
            # feasibility stalls.
            slack_sum = float(np.sum(y)) if y.size else 0.0
            slack_max = float(np.max(y)) if y.size else 0.0
            if slack_sum > slack_tol or (max_viol > feasibility_tol
                                         and viol_trial >= 0.95 * viol_now):
                new_R = min(R * penalty_increase, penalty_max)
                if new_R > R:
                    R = new_R

            history.append({
                "iter": iter_no,
                "accepted": True,
                "f": float(f),
                "viol": _sum_violation(g),
                "max_viol": max_viol,
                "P": float(Phi),
                "penalty": R,
                "step": deltax_norm,
                "alpha": 1.0,
                "rho": float(rho_val) if np.isfinite(rho_val) else float("nan"),
                "predicted_decrease": float(predicted_decrease),
                "actual_decrease": float(actual_decrease),
                "slack_sum": slack_sum,
                "slack_max": slack_max,
                "max_delta_y": float(np.max(xUcur - xLcur)),
                "lp_status": int(lp_status),
                "lp_message": lp_message,
                "backtracks": int(backtrack_iter - 1),
            })

            if display:
                print(
                    f"  {iter_no:>5d}  {f:12.4e}  {_sum_violation(g):12.4e}  {Phi:12.4e}  "
                    f"{deltax_norm:10.3e}  {optimality_norm:10.3e}  "
                    f"{(rho_val if np.isfinite(rho_val) else 0.0):8.2e}  "
                    f"{slack_sum:10.3e}  {R:10.3e}  {nfev[0]:>6d}"
                )

            if callback is not None:
                try:
                    callback(x.copy())
                except Exception:
                    pass
        else:
            # No acceptable trial point in this outer iteration.
            history.append({
                "iter": iter_no,
                "accepted": False,
                "f": float(f),
                "viol": _sum_violation(g),
                "max_viol": _max_violation(g),
                "P": float(Phi),
                "penalty": R,
                "step": deltax_norm,
                "alpha": 0.0,
                "rho": float(rho_val) if np.isfinite(rho_val) else float("nan"),
                "predicted_decrease": float(predicted_decrease),
                "actual_decrease": float(actual_decrease),
                "slack_sum": float(np.sum(y)) if y.size else 0.0,
                "slack_max": float(np.max(y)) if y.size else 0.0,
                "max_delta_y": float(np.max(xUcur - xLcur)),
                "lp_status": int(lp_status),
                "lp_message": lp_message,
                "backtracks": int(backtrack_iter),
            })

            if display:
                print(
                    f"  {iter_no:>5d}  {f:12.4e}  {_sum_violation(g):12.4e}  {Phi:12.4e}  "
                    f"{deltax_norm:>10.3e}  {optimality_norm:10.3e}  "
                    f"REJECT  "
                    f"{(float(np.sum(y)) if y.size else 0.0):10.3e}  {R:10.3e}  {nfev[0]:>6d}"
                )

            # Bump penalty so the next outer iteration cares more about
            # feasibility, then continue if the inner search merely
            # exhausted backtracks; stop if the trust region collapsed.
            R = min(R * penalty_increase, penalty_max)

            if float(np.max(xUcur - xLcur)) < min_dv_box_limit:
                message = (
                    "Stopped: no acceptable step and trust region is below "
                    f"move_limit_min ({move_limit_min:.2e})."
                )
                break

        if nfev[0] >= 10 * maxiter * max(1, n + 1):
            message = (
                "Stopped: function-evaluation budget exhausted "
                f"(nfev={nfev[0]:d})."
            )
            break

    if display:
        viol = _sum_violation(g)
        print("-" * 110)
        print(
            f"  Done: success={success}  iters={iter_no}  nfev={nfev[0]}  "
            f"f={f:.6e}  viol={viol:.3e}  accepted={accepted_steps}  rejected={rejected_steps}"
        )
        print(f"  Message: {message}")

    return SLPResult(
        x=x,
        fun=float(f),
        nit=int(iter_no),
        nfev=int(nfev[0]),
        success=success,
        message=message,
        history=history,
    )


# ----------------------------------------------------------------------
# Smoke test
# ----------------------------------------------------------------------
if __name__ == "__main__":
    def obj(x: np.ndarray) -> float:
        return (x[0] - 1.0) ** 2 + (x[1] - 2.0) ** 2

    def con(x: np.ndarray) -> np.ndarray:
        return np.array([
            4.0 - x[0] ** 2 - x[1] ** 2,   # circle: 4 - x0^2 - x1^2 >= 0
            x[0] + x[1] - 1.0,             # halfspace: x0 + x1 - 1 >= 0
        ])

    res = solve_slp_mvp(
        obj=obj,
        con=con,
        x0=[0.5, 0.5],
        xl=[-5.0, -5.0],
        xu=[5.0, 5.0],
        fd_step=[1e-5, 1e-5],
        fd_type="central",
        maxiter=80,
        move_limit=0.25,
        display=True,
    )
    print(f"x* = {res.x}")
    print(f"f* = {res.fun}")
    print(f"success = {res.success}")
    print(f"message = {res.message}")
    print(f"nit = {res.nit}, nfev = {res.nfev}")
