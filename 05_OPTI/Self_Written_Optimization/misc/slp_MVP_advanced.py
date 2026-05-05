"""
Readable Sequential Linear Programming (SLP) solver.

Problem format follows PySLSQP:
    min   f(x)
    s.t.  c(x) >= 0
          xl <= x <= xu

The LP subproblem is solved in scaled variables y so variables with different
physical units do not dominate the linear model purely because of their units.
At each iteration, the LP proposes dy and optional linearized-constraint slacks:

    min   grad_y(f)^T dy + M sum(s)
    s.t.  c(x_k) + J_y(x_k) dy + s >= 0
          s >= 0

The nonlinear step is globalized with backtracking on the LP direction and an
L1 merit function. The history intentionally records predicted-vs-actual merit
quality because that is the fastest way to debug rejected SLP steps.
"""
# Information for scaling
#  1. Evaluate gradient/Jacobian in physical space (g_x, J_x).
#  2. Scale them into y-space by multiplying by scale.
#  3. LP solves for a scaled step dy inside the scaled bound box and scaled trust region.
#  4. Backtrack on α in scaled space: y_trial = y + α · dy.
#  5. Unscale to evaluate: x_trial = offset + scale · y_trial, then call obj(x_trial) / con(x_trial).
#  6. Accept/reject, update both y and x.

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
    """Return type for :func:`solve_slp_mvp`."""

    x: np.ndarray
    fun: float
    nit: int = 0
    nfev: int = 0
    success: bool = False
    message: str = ""
    history: List[dict] = field(default_factory=list)


def _fd_step_array(fd_step, n: int) -> np.ndarray:
    if fd_step is None:
        return np.full(n, np.sqrt(np.finfo(float).eps), dtype=float)
    arr = np.asarray(fd_step, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr), dtype=float)
    if arr.size != n:
        raise ValueError(f"fd_step must be a scalar or have length {n}.")
    return np.maximum(np.abs(arr.astype(float).ravel()), 1e-12)


def _as_bounds(x0: np.ndarray, xl: Optional[ArrayLike], xu: Optional[ArrayLike]) -> tuple[np.ndarray, np.ndarray]:
    n = x0.size
    xL = np.full(n, -np.inf, dtype=float) if xl is None else np.asarray(xl, dtype=float).ravel()
    xU = np.full(n, np.inf, dtype=float) if xu is None else np.asarray(xu, dtype=float).ravel()
    if xL.size != n or xU.size != n:
        raise ValueError("xl and xu must have the same length as x0.")
    if np.any(xL > xU):
        raise ValueError("Each lower bound must be <= the corresponding upper bound.")
    return xL.copy(), xU.copy()


def _scaling_from_bounds(x: np.ndarray, xl: np.ndarray, xu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return x = offset + scale*y."""
    finite_range = np.isfinite(xl) & np.isfinite(xu) & ((xu - xl) > 0.0)
    scale = np.where(finite_range, xu - xl, np.maximum(1.0, np.abs(x)))
    offset = np.where(finite_range, xl, 0.0)
    return offset.astype(float), scale.astype(float)


def _to_y(x: np.ndarray, offset: np.ndarray, scale: np.ndarray) -> np.ndarray:
    return (x - offset) / scale


def _to_x(y: np.ndarray, offset: np.ndarray, scale: np.ndarray, xl: np.ndarray, xu: np.ndarray) -> np.ndarray:
    return np.clip(offset + scale * y, xl, xu)


def _constraint_violation(c_val: np.ndarray) -> float:
    if c_val.size == 0:
        return 0.0
    return float(np.sum(np.maximum(0.0, -c_val)))


def _max_constraint_violation(c_val: np.ndarray) -> float:
    if c_val.size == 0:
        return 0.0
    return float(np.max(np.maximum(0.0, -c_val)))


def _merit(f_val: float, c_val: np.ndarray, penalty: float) -> float:
    return float(f_val) + float(penalty) * _constraint_violation(c_val)


def _fd_points(x: np.ndarray, i: int, h_abs: float, xl: np.ndarray, xu: np.ndarray, fd_type: str):
    """Return finite-difference points and denominator while respecting bounds."""
    h = max(float(abs(h_abs)), 1e-12)
    can_forward = x[i] + h <= xu[i] if np.isfinite(xu[i]) else True
    can_backward = x[i] - h >= xl[i] if np.isfinite(xl[i]) else True

    if fd_type == "central" and can_forward and can_backward:
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        return "central", xp, xm, 2.0 * h

    if fd_type == "backward" and can_backward:
        xm = x.copy()
        xm[i] -= h
        return "backward", None, xm, h

    if can_forward:
        xp = x.copy()
        xp[i] += h
        return "forward", xp, None, h

    if can_backward:
        xm = x.copy()
        xm[i] -= h
        return "backward", None, xm, h

    return "fixed", None, None, 1.0


def _grad_obj_fd(
    obj: ObjectiveFn,
    x: np.ndarray,
    f0: float,
    h: np.ndarray,
    xl: np.ndarray,
    xu: np.ndarray,
    fd_type: str = "forward",
) -> tuple[np.ndarray, int]:
    n = x.size
    g = np.zeros(n, dtype=float)
    evals = 0
    for i in range(n):
        mode, xp, xm, denom = _fd_points(x, i, h[i], xl, xu, fd_type)
        if mode == "central":
            g[i] = (float(obj(xp)) - float(obj(xm))) / denom
            evals += 2
        elif mode == "backward":
            g[i] = (f0 - float(obj(xm))) / denom
            evals += 1
        elif mode == "forward":
            g[i] = (float(obj(xp)) - f0) / denom
            evals += 1
    return g, evals


def _jac_con_fd(
    con: ConstraintFn,
    x: np.ndarray,
    c0: np.ndarray,
    h: np.ndarray,
    xl: np.ndarray,
    xu: np.ndarray,
    fd_type: str = "forward",
) -> tuple[np.ndarray, int]:
    n = x.size
    m = c0.size
    J = np.zeros((m, n), dtype=float)
    if m == 0:
        return J, 0
    evals = 0
    for j in range(n):
        mode, xp, xm, denom = _fd_points(x, j, h[j], xl, xu, fd_type)
        if mode == "central":
            J[:, j] = (
                np.asarray(con(xp), dtype=float).ravel()
                - np.asarray(con(xm), dtype=float).ravel()
            ) / denom
            evals += 2
        elif mode == "backward":
            J[:, j] = (c0 - np.asarray(con(xm), dtype=float).ravel()) / denom
            evals += 1
        elif mode == "forward":
            J[:, j] = (np.asarray(con(xp), dtype=float).ravel() - c0) / denom
            evals += 1
    return J, evals


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
    move_limit_expand: float = 1.4,
    move_limit_shrink: float = 0.5,
    infeasibility_penalty: float = 1000.0,
    penalty_increase: float = 10.0,
    penalty_max: float = 1e9,
    feasibility_tol: float = 1e-5,
    slack_tol: float = 1e-8,
    xtol: float = 1e-4,
    gtol: float = 1e-6,
    ftol: float = 1e-3,
    backtrack_max: int = 8,
    backtrack_shrink: float = 0.5,
    sufficient_decrease: float = 1e-4,
    feasibility_reduction: float = 1e-3,
    display: bool = False,
    callback: Optional[Callable[[np.ndarray], None]] = None,
    **kwargs,
) -> SLPResult:
    """
    Sequential Linear Programming with scaled variables, slack-relaxed
    constraints, backtracking, and diagnostic history.
    """
    x = np.asarray(x0, dtype=float).ravel().copy()
    n = x.size
    xL, xU = _as_bounds(x, xl, xu)
    x = np.clip(x, xL, xU)

    h = _fd_step_array(fd_step, n)
    offset, scale = _scaling_from_bounds(x, xL, xU)
    y = _to_y(x, offset, scale)
    yL = _to_y(xL, offset, scale)
    yU = _to_y(xU, offset, scale)

    f = float(obj(x))
    c = np.asarray(con(x), dtype=float).ravel() if con else np.empty(0, dtype=float)
    m = c.size
    nfev = 1
    penalty = float(infeasibility_penalty)
    P = _merit(f, c, penalty)
    viol = _constraint_violation(c)
    max_viol = _max_constraint_violation(c)

    delta_y = np.full(n, float(move_limit), dtype=float)
    delta_y = np.clip(delta_y, move_limit_min, np.inf)

    history = [{
        "iter": 0,
        "accepted": True,
        "f": float(f),
        "viol": viol,
        "max_viol": max_viol,
        "P": float(P),
        "penalty": penalty,
        "step": 0.0,
        "alpha": 0.0,
        "rho": np.nan,
        "predicted_decrease": 0.0,
        "actual_decrease": 0.0,
        "slack_sum": 0.0,
        "slack_max": 0.0,
        "max_delta_y": float(np.max(delta_y)),
        "lp_status": 0,
        "lp_message": "initial",
    }]

    if display:
        print("*" * 104)
        print("    SLP MVP Optimizer (scaled LP, slacks, backtracking)")
        print(f"    Initial constraints: {m}")
        print("*" * 104)
        print(
            f"  {'iter':>5s}  {'f':>12s}  {'viol':>12s}  {'P':>12s}  "
            f"{'||dx||':>10s}  {'alpha':>8s}  {'rho':>8s}  {'sum(s)':>10s}  {'penalty':>10s}  {'nfev':>6s}"
        )
        print(
            f"  {0:>5d}  {f:12.4e}  {viol:12.4e}  {P:12.4e}  "
            f"{0.0:10.3e}  {'-':>8s}  {'-':>8s}  {0.0:10.3e}  {penalty:10.3e}  {nfev:>6d}"
        )

    if callback:
        callback(x.copy())

    success = False
    message = "Maximum iterations reached"
    iter_no = 0

    for iter_no in range(1, maxiter + 1):
        g_x, fe = _grad_obj_fd(obj, x, f, h, xL, xU, fd_type=fd_type)
        nfev += fe
        if m > 0:
            J_x, ce = _jac_con_fd(con, x, c, h, xL, xU, fd_type=fd_type)
            nfev += ce
        else:
            J_x = np.zeros((0, n), dtype=float)

        g_y = g_x * scale
        J_y = J_x * scale[np.newaxis, :]
        grad_norm = float(np.linalg.norm(g_y, ord=np.inf))

        if m > 0:
            c_lp = np.concatenate((g_y, np.full(m, penalty, dtype=float)))
            A_ub = np.hstack((-J_y, -np.eye(m)))
            b_ub = c.copy()
        else:
            c_lp = g_y
            A_ub = None
            b_ub = None

        dy_bounds = [
            (
                float(max(yL[i] - y[i], -delta_y[i])),
                float(min(yU[i] - y[i], delta_y[i])),
            )
            for i in range(n)
        ]
        bounds = dy_bounds + [(0.0, None) for _ in range(m)]

        res = linprog(c_lp, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success:
            delta_y *= move_limit_shrink
            history.append({
                "iter": iter_no,
                "accepted": False,
                "f": float(f),
                "viol": viol,
                "max_viol": max_viol,
                "P": float(P),
                "penalty": penalty,
                "step": 0.0,
                "alpha": 0.0,
                "rho": np.nan,
                "predicted_decrease": 0.0,
                "actual_decrease": 0.0,
                "slack_sum": np.nan,
                "slack_max": np.nan,
                "max_delta_y": float(np.max(delta_y)),
                "lp_status": int(res.status),
                "lp_message": res.message,
            })
            if display:
                print(f"  {iter_no:>5d}  LP FAILED: {res.message}  nfev={nfev}")
            if np.max(delta_y) < move_limit_min:
                message = "Stopped: LP failed and scaled move limits are below move_limit_min."
                break
            continue

        dy = np.asarray(res.x[:n], dtype=float)
        slack = np.asarray(res.x[n:], dtype=float) if m > 0 else np.empty(0, dtype=float)
        slack_sum = float(np.sum(slack)) if slack.size else 0.0
        slack_max = float(np.max(slack)) if slack.size else 0.0

        linear_c_full = c + J_y @ dy
        model_f_full = f + float(g_y @ dy)
        P_model_full = _merit(model_f_full, linear_c_full, penalty)
        predicted_full = max(0.0, P - P_model_full)

        accepted = False
        best_trial = None
        alpha = 1.0
        rho = np.nan
        actual_decrease = 0.0
        predicted_decrease = predicted_full

        for _ in range(backtrack_max):
            y_trial = np.clip(y + alpha * dy, yL, yU)
            x_trial = _to_x(y_trial, offset, scale, xL, xU)
            f_trial = float(obj(x_trial))
            c_trial = np.asarray(con(x_trial), dtype=float).ravel() if con else np.empty(0, dtype=float)
            nfev += 1

            P_trial = _merit(f_trial, c_trial, penalty)
            viol_trial = _constraint_violation(c_trial)
            actual_decrease = P - P_trial

            linear_c_trial = c + J_y @ (alpha * dy)
            model_f_trial = f + float(g_y @ (alpha * dy))
            P_model_trial = _merit(model_f_trial, linear_c_trial, penalty)
            predicted_decrease = max(0.0, P - P_model_trial)
            rho = actual_decrease / predicted_decrease if predicted_decrease > 1e-14 else np.nan

            sufficient_merit = predicted_decrease <= 1e-14 or (
                P_trial <= P - sufficient_decrease * predicted_decrease
            )
            feasibility_progress = (
                viol > feasibility_tol
                and viol_trial <= viol * (1.0 - feasibility_reduction)
            )

            best_trial = (y_trial, x_trial, f_trial, c_trial, P_trial, viol_trial)
            if sufficient_merit or feasibility_progress:
                accepted = True
                break
            alpha *= backtrack_shrink

        y_step_norm = float(np.linalg.norm(alpha * dy))
        x_step_norm = float(np.linalg.norm(best_trial[1] - x)) if best_trial else 0.0

        if accepted and best_trial is not None:
            f_old = f
            P_old = P
            viol_old = viol
            y, x, f, c, P, viol = best_trial
            max_viol = _max_constraint_violation(c)

            if slack_sum > slack_tol or (viol > feasibility_tol and viol >= 0.95 * viol_old):
                penalty = min(penalty * penalty_increase, penalty_max)

            if np.isfinite(rho) and rho > 0.75:
                delta_y = np.minimum(delta_y * move_limit_expand, move_limit)
            elif np.isfinite(rho) and rho < 0.25:
                delta_y *= move_limit_shrink

            history.append({
                "iter": iter_no,
                "accepted": True,
                "f": float(f),
                "viol": viol,
                "max_viol": max_viol,
                "P": float(P),
                "penalty": penalty,
                "step": x_step_norm,
                "scaled_step": y_step_norm,
                "alpha": float(alpha),
                "rho": float(rho) if np.isfinite(rho) else np.nan,
                "predicted_decrease": float(predicted_decrease),
                "actual_decrease": float(P_old - P),
                "slack_sum": slack_sum,
                "slack_max": slack_max,
                "max_delta_y": float(np.max(delta_y)),
                "lp_status": int(res.status),
                "lp_message": res.message,
            })

            if display:
                print(
                    f"  {iter_no:>5d}  {f:12.4e}  {viol:12.4e}  {P:12.4e}  "
                    f"{x_step_norm:10.3e}  {alpha:8.2e}  "
                    f"{(rho if np.isfinite(rho) else 0.0):8.2e}  {slack_sum:10.3e}  "
                    f"{penalty:10.3e}  {nfev:>6d}"
                )

            if callback:
                callback(x.copy())

            if (
                max_viol <= feasibility_tol
                and slack_max <= max(slack_tol, feasibility_tol)
                and x_step_norm <= xtol
            ):
                success = True
                message = f"Converged: feasible and ||dx|| <= xtol ({xtol:.2e})."
                break
            if (
                max_viol <= feasibility_tol
                and slack_max <= max(slack_tol, feasibility_tol)
                and abs(f - f_old) <= ftol
                and grad_norm <= max(gtol, 1e-8)
            ):
                success = True
                message = "Converged: feasible with small objective and scaled-gradient changes."
                break
        else:
            delta_y *= move_limit_shrink
            if slack_sum > slack_tol:
                penalty = min(penalty * penalty_increase, penalty_max)

            history.append({
                "iter": iter_no,
                "accepted": False,
                "f": float(f),
                "viol": viol,
                "max_viol": max_viol,
                "P": float(P),
                "penalty": penalty,
                "step": x_step_norm,
                "scaled_step": y_step_norm,
                "alpha": float(alpha),
                "rho": float(rho) if np.isfinite(rho) else np.nan,
                "predicted_decrease": float(predicted_decrease),
                "actual_decrease": float(actual_decrease),
                "slack_sum": slack_sum,
                "slack_max": slack_max,
                "max_delta_y": float(np.max(delta_y)),
                "lp_status": int(res.status),
                "lp_message": res.message,
            })

            if display:
                print(
                    f"  {iter_no:>5d}  {f:12.4e}  {viol:12.4e}  {P:12.4e}  "
                    f"{x_step_norm:>10.3e}  REJECT  "
                    f"{(rho if np.isfinite(rho) else 0.0):8.2e}  {slack_sum:10.3e}  "
                    f"{penalty:10.3e}  {nfev:>6d}"
                )

            if np.max(delta_y) < move_limit_min:
                message = "Stopped: no acceptable backtracked step and scaled move limits are below move_limit_min."
                break

    return SLPResult(
        x=x,
        fun=float(f),
        nit=iter_no,
        nfev=nfev,
        success=success,
        message=message,
        history=history,
    )


if __name__ == "__main__":
    def obj(x):
        return (x[0] - 1.0) ** 2 + (x[1] - 2.0) ** 2

    def con(x):
        return np.array([
            4.0 - x[0] ** 2 - x[1] ** 2,
            x[0] + x[1] - 1.0,
        ])

    res = solve_slp_mvp(
        obj=obj,
        con=con,
        x0=[0.5, 0.5],
        xl=[-5.0, -5.0],
        xu=[5.0, 5.0],
        fd_step=[1e-5, 1e-5],
        maxiter=50,
        display=True,
    )
    print(f"x* = {res.x}")
    print(f"f* = {res.fun}")
    print(f"success = {res.success}")
    print(f"msg = {res.message}")
