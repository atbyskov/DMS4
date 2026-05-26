"""
Faithful Python port of fminslp.m (Rene S. Nielsen).

This module is a literal translation of the merit-function variant of the
MATLAB ``fminslp`` solver. Every algorithmic ingredient and every default
matches the MATLAB reference:

    * Slack variables ``y >= 0`` that relax the linearized inequality
      constraints so the LP subproblem is always feasible.
    * Quadratic merit objective  ``f(x) + a_F * sum(y*R + 0.5*y^2)`` (merit)
      or the augmented-Lagrangian variant ``f(x) + a_F * sum(y*lambda + 0.5*R*y^2)``.
    * LP subproblem in the augmented vector ``z = [x, y]``.
    * Fletcher--Leyffer--Toint global convergence filter
      (``beta = 1 - 1e-6``, ``gamma = 1e-6``, ``sigma = 2e-6``).
    * Sufficient-decrease check ``deltaf < sigma * deltaL`` (with deltaL > 0).
    * Adaptive move limits driven by iterate-history oscillation.
    * Convergence on first-order, function-value, step-size, and objective-limit
      criteria -- evaluated unconditionally in an elif chain, exactly as MATLAB
      does (no feasibility precondition).

External API uses the PySLSQP convention ``c(x) >= 0`` (feasible). Internally
the algorithm uses ``g(x) = -c(x) <= 0`` exactly as in MATLAB.
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
# Result type (public; the rest of the pipeline depends on these fields)
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
# Global convergence filter
#   Mirrors initializeGlobalConvergenceFilter (lines 851-881),
#           EvaluateCurrentDesignPointToFilter (lines 884-899),
#           UpdateFilter (lines 902-936).
# ----------------------------------------------------------------------
@dataclass
class _Filter:
    beta: float
    gamma: float
    sigma: float
    vals: np.ndarray   # (capacity, 2) array of (h, f) pairs
    n_vals: int
    init_f: float
    h: float
    f: float
    accepted: bool

    def evaluate(self) -> None:
        """Port of EvaluateCurrentDesignPointToFilter (lines 884-899)."""
        h = self.h
        f = self.f
        for ii in range(self.n_vals):
            hi = self.vals[ii, 0]
            fi = self.vals[ii, 1]
            if (h <= hi * self.beta) or ((f + self.gamma * h) <= fi):
                self.accepted = True
            else:
                self.accepted = False
                break

    def update(self, hk: float, fk: float) -> None:
        """Port of UpdateFilter (lines 902-936)."""
        kept_h: List[float] = []
        kept_f: List[float] = []
        for ii in range(self.n_vals):
            hi = self.vals[ii, 0]
            fi = self.vals[ii, 1]
            dominated = (hk <= hi) and (fk <= fi)
            if not dominated:
                kept_h.append(hi)
                kept_f.append(fi)
        kept_h.append(hk)
        kept_f.append(fk)
        n_new = len(kept_h)
        if n_new > self.vals.shape[0]:
            new_vals = np.zeros((n_new + 32, 2), dtype=float)
            new_vals[:n_new, 0] = kept_h
            new_vals[:n_new, 1] = kept_f
            self.vals = new_vals
        else:
            self.vals[:n_new, 0] = kept_h
            self.vals[:n_new, 1] = kept_f
        self.n_vals = n_new


def _init_filter(max_iter: int, max_infeasibility: float) -> _Filter:
    """Port of initializeGlobalConvergenceFilter (lines 851-881)."""
    small = 1.0e-6  # SmallVal
    capacity = max(max_iter, 16)
    vals = np.zeros((capacity, 2), dtype=float)
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
        h=1.0e30,
        f=1.0e30,
        accepted=False,
    )


# ----------------------------------------------------------------------
# Adaptive move limit (port of AdaptiveMoveLimit, lines 939-1007)
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

    # Lines 941-947
    if reduce_switch:
        expand = reduce_fac
        reduction = reduce_fac
    else:
        reduction = reduce_fac
        expand = expand_fac

    if min_dv_box_limit is None:
        min_dv_box_limit = 0.0

    n = x.size
    for i in range(n):
        # Line 955: previous allowable half-width
        delta = 0.5 * (xUcur[i] - xLcur[i])

        # Lines 957-967: oscillation detection
        if abs(x[i] - xold1[i]) > 1.0e-10:
            s1 = (xold1[i] - xold2[i]) / (x[i] - xold1[i])
            if s1 < 0.0:
                delta *= reduction
            else:
                delta *= expand
        else:
            delta *= move_limit

        # Lines 968-971: cap delta at the absolute move-limit fraction
        dmax = (xUorg[i] - xLorg[i]) * move_limit
        if delta > dmax:
            delta = dmax

        # Lines 973-977: initial new bounds, then clip to global feasibility
        xLcur[i] = x[i] - delta
        xUcur[i] = x[i] + delta
        xLcur[i] = max(xLcur[i], xLorg[i])
        xUcur[i] = min(xUcur[i], xUorg[i])

        # Lines 980-993: symmetric safety guards (literal port; the
        # second branch in MATLAB has an apparent sign typo, kept as-is
        # for faithfulness -- in practice it is dead code because the
        # first branch already restores L < U).
        if xLcur[i] > xUcur[i]:
            if min_dv_box_limit > 0.0:
                xLcur[i] = xUcur[i] - min_dv_box_limit / 2.0
            else:
                xLcur[i] = (1.0 - 1.0e-6) * xUcur[i]
        if xUcur[i] < xLcur[i]:
            if min_dv_box_limit > 0.0:
                xUcur[i] = xLcur[i] - min_dv_box_limit / 2.0
            else:
                xUcur[i] = (1.0 + 1.0e-6) * xLcur[i]

        # Lines 996-1005: enforce minimum-box-width floor
        if (xUcur[i] - xLcur[i]) < min_dv_box_limit:
            if (xUcur[i] <= xUorg[i] - min_dv_box_limit / 2.0) and \
               (xLcur[i] <= xLorg[i] - min_dv_box_limit / 2.0):
                xUcur[i] = xUcur[i] + min_dv_box_limit / 2.0
                xLcur[i] = xLcur[i] - min_dv_box_limit / 2.0
            elif xUcur[i] <= xUorg[i] - min_dv_box_limit:
                xUcur[i] = xUcur[i] + min_dv_box_limit
            elif xLcur[i] >= xLorg[i] - min_dv_box_limit:
                xLcur[i] = xLcur[i] - min_dv_box_limit

    return xLcur, xUcur


# ----------------------------------------------------------------------
# Finite-difference gradient + constraint Jacobian
#   Mirrors getFunDSA (lines 611-652) + getNonlconDSA (lines 654-745).
#   MATLAB calls them separately; we batch them via obj_and_g for
#   evaluator-cache efficiency.  Result is identical.
# ----------------------------------------------------------------------
def _fd_grad_and_jac(
    obj_and_g: Callable[[np.ndarray], Tuple[float, np.ndarray]],
    x: np.ndarray,
    f0: float,
    g0: np.ndarray,
    h: np.ndarray,
    fd_type: str,
) -> Tuple[np.ndarray, np.ndarray]:
    n = x.size
    m = g0.size
    df = np.zeros(n, dtype=float)
    dg = np.zeros((m, n), dtype=float) if m > 0 else np.zeros((0, n), dtype=float)

    fd = fd_type.lower()
    x_work = x.copy()

    if fd == "forward":
        for i in range(n):
            xin = x_work[i]
            x_work[i] = xin + h[i]
            f_p, g_p = obj_and_g(x_work)
            df[i] = (f_p - f0) / h[i]
            if m > 0:
                dg[:, i] = (g_p - g0) / h[i]
            x_work[i] = xin
    elif fd == "backward":
        for i in range(n):
            xin = x_work[i]
            x_work[i] = xin - h[i]
            f_m, g_m = obj_and_g(x_work)
            df[i] = (f0 - f_m) / h[i]
            if m > 0:
                dg[:, i] = (g0 - g_m) / h[i]
            x_work[i] = xin
    elif fd == "central":
        for i in range(n):
            xin = x_work[i]
            x_work[i] = xin + h[i]
            f_p, g_p = obj_and_g(x_work)
            x_work[i] = xin - h[i]
            f_m, g_m = obj_and_g(x_work)
            df[i] = (f_p - f_m) / (2.0 * h[i])
            if m > 0:
                dg[:, i] = (g_p - g_m) / (2.0 * h[i])
            x_work[i] = xin
    else:
        raise ValueError(f"Unknown fd_type: {fd_type}")

    return df, dg


# ----------------------------------------------------------------------
# LP solver wrapper (MATLAB lpSolver with 'linprog' branch, lines 564-567)
# ----------------------------------------------------------------------
def _solve_lp(
    c: np.ndarray,
    A_ub: Optional[np.ndarray],
    b_ub: Optional[np.ndarray],
    lb: np.ndarray,
    ub: np.ndarray,
) -> Tuple[np.ndarray, bool, str]:
    bounds = [(float(lb[i]), float(ub[i])) for i in range(c.size)]
    A_ub_in = A_ub if (A_ub is not None and A_ub.size > 0) else None
    b_ub_in = b_ub if (b_ub is not None and b_ub.size > 0) else None

    res = linprog(
        c=c,
        A_ub=A_ub_in,
        b_ub=b_ub_in,
        bounds=bounds,
        method="highs",
    )
    if res.success:
        return np.asarray(res.x, dtype=float), True, str(res.message)
    return np.zeros_like(c), False, str(res.message)


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------
def _fd_step_array(fd_step, n: int) -> np.ndarray:
    if fd_step is None:
        # MATLAB default: FiniteDifferenceStepSize = sqrt(eps)
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
    """max(g, 0) in the internal g <= 0 convention."""
    if g.size == 0:
        return 0.0
    return float(np.max(np.maximum(g, 0.0)))


def _sum_violation(g: np.ndarray) -> float:
    if g.size == 0:
        return 0.0
    return float(np.sum(np.maximum(g, 0.0)))


# ----------------------------------------------------------------------
# Main solver  --  literal port of fminslp.solve (lines 196-496)
# ----------------------------------------------------------------------
def solve_slp_mvp(
    obj: ObjectiveFn,
    con: Optional[ConstraintFn],
    x0: ArrayLike,
    xl: Optional[ArrayLike] = None,
    xu: Optional[ArrayLike] = None,
    *,
    # ------- MATLAB-equivalent options (defaults from slpoptions) --------
    fd_step=None,                                 # FiniteDifferenceStepSize
    fd_type: str = "forward",                     # FiniteDifferenceType
    maxiter: int = 1000,                          # MaxIterations
    max_fun_evals: int = 1000,                    # MaxFunctionEvaluations  (counts LP trials, matches MATLAB nFeval)
    infeasibility_penalty: float = 1000.0,        # InfeasibilityPenalization (R)
    gtol: float = 1e-6,                           # OptimalityTolerance
    ftol: float = 1e-6,                           # FunctionTolerance
    xtol: float = 1e-8,                           # StepTolerance
    objective_limit: float = -1e20,               # ObjectiveLimit
    move_limit: float = 0.10,                     # MoveLimit
    move_limit_expand: float = 1.1,               # MoveLimitExpand
    move_limit_shrink: float = 0.5,               # MoveLimitReduce
    max_infeasibility: float = np.inf,            # MaxInfeasibility
    algorithm: str = "merit",                     # Algorithm  ('merit' or 'al')
    display: bool = False,
    callback: Optional[Callable[[np.ndarray], None]] = None,
    # ------- Tolerated for backward compatibility (silently ignored) ----
    # The previous version of this solver had extra heuristics that are
    # not part of fminslp.m. They are accepted but ignored so existing
    # callers do not break.
    **_unused_kwargs: Any,
) -> SLPResult:
    """Sequential Linear Programming solver, faithful port of fminslp.m.

    Solves ``min f(x) s.t. c(x) >= 0, xl <= x <= xu``.
    """
    algorithm = str(algorithm).lower()
    if algorithm not in {"merit", "al"}:
        raise ValueError("algorithm must be 'merit' or 'al'.")

    # ---------------- Setup (MATLAB lines 196-305) ----------------
    x = np.asarray(x0, dtype=float).ravel().copy()
    n = x.size
    xL, xU = _as_bounds(x, xl, xu)
    x = np.clip(x, xL, xU)

    h_fd = _fd_step_array(fd_step, n)

    # Line 225: minDVBoxLimit = 2 * StepTolerance
    min_dv_box_limit = 2.0 * xtol

    # Line 228: x = x0(:);  bounds initialization
    xLcur = xL.copy()
    xUcur = xU.copy()
    xold1 = x.copy()
    xold2 = x.copy()

    # Evaluation counter -- counts true FEA calls (used for SLPResult.nfev,
    # which downstream CSV row count depends on).
    nfev = [0]

    def obj_and_g(xv: np.ndarray) -> Tuple[float, np.ndarray]:
        nfev[0] += 1
        f_val = float(obj(np.asarray(xv, dtype=float)))
        if con is None:
            return f_val, np.empty(0, dtype=float)
        c_val = np.asarray(con(np.asarray(xv, dtype=float)), dtype=float).ravel()
        # PySLSQP c >= 0  ->  MATLAB g <= 0
        return f_val, -c_val

    # Lines 239-256: initial constraint evaluation, slack-variable init.
    # We need one extra (fun, nonlcon) call here because MATLAB calls
    # nonlcon BEFORE fun (line 240), then fun (line 287).  We collapse
    # both into a single obj_and_g call -- algorithmically identical,
    # only saves one evaluation.
    f_cur, g_cur = obj_and_g(x)
    m = g_cur.size

    if m > 0:
        # Line 244:  y = max([gnl(:); gneq(:); -gneq(:)], 0)
        # We have no separate equalities -- y = max(g, 0).
        y = np.maximum(g_cur, 0.0)
        # Line 247:  yub = max(max(g, 1)) * 1e6
        yub_seed = max(float(np.max(g_cur)), 1.0)
        ylb = np.zeros(m, dtype=float)
        yub = np.full(m, yub_seed * 1.0e6, dtype=float)
    else:
        y = np.empty(0, dtype=float)
        ylb = np.empty(0, dtype=float)
        yub = np.empty(0, dtype=float)

    # Lines 275-284: AL multipliers
    R = float(infeasibility_penalty)  # InfeasibilityPenalization, never grows
    if algorithm == "al" and m > 0:
        # Lines 279-281
        lambda_vec = np.ones(m, dtype=float)
        lambda_vec = np.maximum(lambda_vec + R * y, 1.0)
    else:
        lambda_vec = np.empty(0, dtype=float)

    # Lines 287-289: merit scaling factor a_F is set once
    a_F = max(abs(float(f_cur)), 1.0)

    # ---------- Merit helpers (mirror getMeritObj lines 503-531) ----------
    def merit_val(f_v: float, y_v: np.ndarray, lam: np.ndarray) -> float:
        if y_v.size == 0:
            return float(f_v)
        if algorithm == "merit":
            return float(f_v) + a_F * float(np.sum(y_v * R + 0.5 * y_v * y_v))
        # algorithm == 'al'
        return float(f_v) + a_F * float(np.sum(y_v * lam + 0.5 * R * y_v * y_v))

    def merit_grad_y(y_v: np.ndarray, lam: np.ndarray) -> np.ndarray:
        if y_v.size == 0:
            return np.empty(0, dtype=float)
        if algorithm == "merit":
            return a_F * (R + y_v)
        return a_F * (lam + R * y_v)

    # Lines 292-296: initialize filter
    fmerit = merit_val(f_cur, y, lambda_vec)
    filt = _init_filter(max_iter=maxiter, max_infeasibility=max_infeasibility)
    filt.init_f = fmerit
    f_old = fmerit

    # Lines 299-300
    max_inf = float(np.max(y)) if m > 0 else 0.0
    optimality_norm = 0.0

    # ---------- History (preserves the keys the old solver wrote) ----------
    history: List[dict] = [{
        "iter": 0,
        "accepted": True,
        "f": float(f_cur),
        "viol": _sum_violation(g_cur),
        "max_viol": _max_violation(g_cur),
        "P": float(fmerit),
        "penalty": R,
        "step": 0.0,
        "rho": float("nan"),
        "predicted_decrease": 0.0,
        "actual_decrease": 0.0,
        "slack_sum": float(np.sum(y)) if y.size else 0.0,
        "slack_max": float(np.max(y)) if y.size else 0.0,
        "lp_message": "initial",
    }]

    if display:
        print("*" * 105)
        print("                          fminslp optimizer with global convergence filter")
        print("*" * 105)
        print(f"    Algorithm: {algorithm}   nDV: {n}   nConstraints: {m}")
        print(f"    a_F = {a_F:.4e}   R = {R:.4e}   sigma = {filt.sigma:.2e}   "
              f"sqrt(eps) = {float(np.sqrt(np.finfo(float).eps)):.2e}")
        print(f"    StepTol={xtol:.2e}  FunTol={ftol:.2e}  OptTol={gtol:.2e}  "
              f"MaxIter={maxiter}  MaxFevals={max_fun_evals}")
        print("-" * 105)
        print(
            f"  {'f(x)':>14s}  {'Max inf':>12s}  {'Norm dx':>12s}  "
            f"{'nFeval':>8s}  {'IterNo':>8s}"
        )
        print(
            f"  {f_cur:>14.6e}  {max_inf:>12.4e}  {0.0:>12.4e}  "
            f"{nfev[0]:>8d}  {0:>8d}"
        )

    if callback is not None:
        try:
            callback(x.copy())
        except Exception:
            pass

    # ---------- Main loop (MATLAB lines 307-478) ----------
    iter_no = 0
    n_trials = 0           # mirrors MATLAB nFeval semantics (trial counter)
    exit_flag = -1
    success = False
    message = "Maximum iterations reached"

    # Pre-initialize trial-scope variables so they are defined on any exit
    f_new = f_cur
    g_new = g_cur.copy()
    y_new = y.copy()
    fmerit_new = fmerit
    delta_norm = 0.0
    delta_f = 0.0
    delta_L = 0.0
    lp_message = ""
    last_lp_failed = False
    last_trial_accepted = False

    optimize = True
    while optimize:
        iter_no += 1

        # --- Gradients at current x (MATLAB lines 312-314) ---
        df, dg = _fd_grad_and_jac(obj_and_g, x, f_cur, g_cur, h_fd, fd_type)
        dy = merit_grad_y(y, lambda_vec)
        dfmerit = np.concatenate([df, dy]) if m > 0 else df.copy()

        # --- LP linearization (MATLAB lines 316-326) ---
        if m > 0:
            dgmerit = np.zeros((m, n + m), dtype=float)
            dgmerit[:, :n] = dg
            dgmerit[:, n:] = -np.eye(m)
            gmerit_cur = g_cur - y
            z_cur = np.concatenate([x, y])
            A_ub = dgmerit
            b_ub = dgmerit @ z_cur - gmerit_cur
        else:
            A_ub = None
            b_ub = None

        # --- Outer move-limit update (MATLAB lines 328-330) ---
        xLcur, xUcur = _adaptive_move_limit(
            x, xLcur, xUcur, xL, xU,
            move_limit, move_limit_shrink, move_limit_expand,
            xold1, xold2, reduce_switch=False,
            min_dv_box_limit=min_dv_box_limit,
        )

        # --- Update old values BEFORE inner loop (MATLAB lines 333-334) ---
        xold2 = xold1.copy()
        xold1 = x.copy()

        # --- Inner backtrack loop (MATLAB lines 340-447) ---
        backtrack = True
        add_to_filter = False
        last_trial_accepted = False

        while backtrack:
            # MATLAB lines 343-344: LP bounds over [x, y]
            xL_z = np.concatenate([xLcur, ylb])
            xU_z = np.concatenate([xUcur, yub])

            # MATLAB line 347: solve LP
            z_new, lp_ok, lp_message = _solve_lp(
                dfmerit, A_ub, b_ub, xL_z, xU_z,
            )

            if lp_ok:
                # --- Trial point (MATLAB lines 353-359) ---
                x_new = z_new[:n]
                y_new = z_new[n:] if m > 0 else np.empty(0, dtype=float)

                if m > 0:
                    delta_xy = z_new - np.concatenate([x, y])
                else:
                    delta_xy = x_new - x
                delta_x = x_new - x
                delta_norm = float(np.linalg.norm(delta_x))

                # MATLAB line 362: optimality_norm = ||dfmerit(1:nDV)||
                optimality_norm = float(np.linalg.norm(dfmerit[:n]))

                # MATLAB lines 365-370: evaluate constraints + objective at trial
                f_new, g_new = obj_and_g(x_new)
                n_trials += 1
                fmerit_new = merit_val(f_new, y_new, lambda_vec)
                delta_f = f_old - fmerit_new

                # MATLAB lines 376-378: build filter candidate (h, f)
                if m > 0:
                    gmerit_new = g_new - y_new
                else:
                    gmerit_new = np.empty(0, dtype=float)
                filt.h = _max_violation(gmerit_new)
                denom = filt.init_f if abs(filt.init_f) > 1.0e-30 else 1.0
                filt.f = fmerit_new / denom

                # MATLAB line 381: filter test
                filt.evaluate()

                # MATLAB lines 384-403: accept/reject logic
                add_to_filter = False
                if filt.accepted:
                    delta_L = float(-dfmerit @ delta_xy)
                    if (delta_f < filt.sigma * delta_L) and (delta_L > 0.0):
                        # Insufficient decrease relative to predicted -> retry smaller
                        reduce_switch = True
                    else:
                        # Accept the step
                        reduce_switch = False
                        backtrack = False
                        last_trial_accepted = True
                        if filt.h > 0.0:
                            add_to_filter = True
                else:
                    reduce_switch = True

                # MATLAB lines 405-407: shrink move limits if step rejected
                if reduce_switch:
                    xLcur, xUcur = _adaptive_move_limit(
                        x, xLcur, xUcur, xL, xU,
                        move_limit, move_limit_shrink, move_limit_expand,
                        xold1, xold2, reduce_switch=True,
                        min_dv_box_limit=min_dv_box_limit,
                    )

                # --- Convergence checks (MATLAB lines 410-440) ---
                # NOTE: these run REGARDLESS of accept/reject, matching MATLAB.
                # The elif chain mirrors MATLAB's precedence exactly.
                if optimality_norm <= gtol:
                    optimize = False
                    backtrack = False
                    exit_flag = 1
                    success = True
                    message = f"Successfully solved to Optimality Tolerance: <= {gtol:.5e}"
                elif abs(delta_f) <= ftol:
                    optimize = False
                    backtrack = False
                    exit_flag = 2
                    success = True
                    message = f"Successfully solved to Function Tolerance: <= {ftol:.5e}"
                elif delta_norm <= xtol:
                    optimize = False
                    backtrack = False
                    exit_flag = 3
                    success = True
                    message = f"Successfully solved to Step Tolerance: <= {xtol:.5e}"
                elif (fmerit_new <= objective_limit) and iter_no > 1:
                    optimize = False
                    backtrack = False
                    exit_flag = 4
                    success = True
                    message = f"Successfully solved to Objective Limit: <= {objective_limit:.5e}"
                elif iter_no >= maxiter:
                    optimize = False
                    backtrack = False
                    exit_flag = 0
                    success = False
                    message = f"Number of iterations exceeded the limit: {maxiter}"
                elif n_trials >= max_fun_evals:
                    optimize = False
                    backtrack = False
                    exit_flag = 0
                    success = False
                    message = f"Number of function evaluations exceeded the limit: {max_fun_evals}"
            else:
                # --- LP failure (MATLAB lines 441-446) ---
                optimize = False
                backtrack = False
                last_lp_failed = True
                exit_flag = -1
                success = False
                message = f"LP solver failed: {lp_message}"

        # end inner loop

        # MATLAB lines 449-477: only commit + log if NOT terminating
        if optimize:
            # MATLAB lines 451-453: add point to filter if requested
            if add_to_filter:
                filt.update(filt.h, filt.f)

            # MATLAB lines 456-457: commit step
            x = x_new.copy()
            y = y_new.copy() if m > 0 else y
            f_cur = f_new
            g_cur = g_new.copy()
            fmerit = fmerit_new

            # MATLAB lines 459-462: AL multiplier update
            if algorithm == "al" and m > 0:
                lambda_vec = np.maximum(lambda_vec + R * y, 1.0)

            # MATLAB lines 465-471: history
            f_old = fmerit
            max_inf = _max_violation(g_cur)

            rho_val = (delta_f / delta_L) if delta_L > 1.0e-14 else float("nan")
            history.append({
                "iter": iter_no,
                "accepted": True,
                "f": float(f_cur),
                "viol": _sum_violation(g_cur),
                "max_viol": max_inf,
                "P": float(fmerit),
                "penalty": R,
                "step": delta_norm,
                "rho": float(rho_val),
                "predicted_decrease": float(delta_L),
                "actual_decrease": float(delta_f),
                "slack_sum": float(np.sum(y)) if y.size else 0.0,
                "slack_max": float(np.max(y)) if y.size else 0.0,
                "lp_message": lp_message,
            })

            # MATLAB lines 473-476
            if display:
                print(
                    f"  {f_cur:>14.6e}  {max_inf:>12.4e}  {delta_norm:>12.4e}  "
                    f"{nfev[0]:>8d}  {iter_no:>8d}"
                )

            if callback is not None:
                try:
                    callback(x.copy())
                except Exception:
                    pass
        else:
            # Terminating (success, max iter/fevals, or LP failure).
            # Record a final history row so plotters see the trial outcome.
            history.append({
                "iter": iter_no,
                "accepted": bool(last_trial_accepted and not last_lp_failed),
                "f": float(f_new),
                "viol": _sum_violation(g_new),
                "max_viol": _max_violation(g_new),
                "P": float(fmerit_new),
                "penalty": R,
                "step": delta_norm,
                "rho": float("nan"),
                "predicted_decrease": float(delta_L),
                "actual_decrease": float(delta_f),
                "slack_sum": float(np.sum(y_new)) if y_new.size else 0.0,
                "slack_max": float(np.max(y_new)) if y_new.size else 0.0,
                "lp_message": lp_message,
                "termination": message,
            })
            if display:
                tag = "FAILED" if last_lp_failed else ("CONVERGED" if success else "STOPPED")
                print(
                    f"  {f_new:>14.6e}  {_max_violation(g_new):>12.4e}  "
                    f"{delta_norm:>12.4e}  {nfev[0]:>8d}  {iter_no:>8d}  [{tag}]"
                )

    # ----- Final output preparation (MATLAB lines 480-494) -----
    # MATLAB sets fval = freal (the trial's real f at exit).  When the
    # trial was filter-accepted we also commit x to the trial point so
    # that result.x and result.fun are a consistent (x, f(x)) pair --
    # this is a small, deliberate clarification of MATLAB's behaviour
    # which leaves x stale while reporting freal.
    if last_trial_accepted and not last_lp_failed:
        x_out = x_new.copy()
        f_out = float(f_new)
        g_out = g_new
    else:
        x_out = x.copy()
        f_out = float(f_cur)
        g_out = g_cur

    if display:
        print("-" * 105)
        print(
            f"  Done: success={success}  iters={iter_no}  trials={n_trials}  "
            f"nfev={nfev[0]}  f={f_out:.6e}  max_inf={_max_violation(g_out):.3e}"
        )
        print(f"  Exit flag: {exit_flag}")
        print(f"  Message: {message}")

    return SLPResult(
        x=x_out,
        fun=f_out,
        nit=iter_no,
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
        fd_step=1e-5,
        fd_type="forward",
        maxiter=100,
        move_limit=0.25,
        infeasibility_penalty=1000.0,
        display=True,
    )
    print(f"x* = {res.x}")
    print(f"f* = {res.fun}")
    print(f"success = {res.success}")
    print(f"message = {res.message}")
    print(f"nit = {res.nit}, nfev = {res.nfev}")
