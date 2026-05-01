"""
Minimum Viable Product (MVP) of a Sequential Linear Programming (SLP) solver.

This implementation focuses on maximum readability and simplicity.
It uses:
- Forward difference gradients
- A basic LP subproblem with slack-relaxed linearized constraints
- The simplest possible move limits (trust region) that shrink on rejection
- Simple stopping criteria (step size / move limit size)

Problem format:
    min   f(x)
    s.t.  c(x) >= 0
          xl <= x <= xu

At each SLP iteration the LP subproblem uses variables z = [dx, s]:
    min   grad(f)^T dx + M sum(s)
    s.t.  c(x_k) + J(x_k) dx + s >= 0
          s >= 0
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
    """Return type for :func:`solve_slp_mvp`."""
    x: np.ndarray
    fun: float
    nit: int = 0
    nfev: int = 0
    success: bool = False
    message: str = ""
    history: List[dict] = field(default_factory=list)

# ----------------------------------------------------------------------
# Finite-difference gradients (Ported from slp.py)
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

def _constraint_violation(c_val: np.ndarray) -> float:
    """Total violation for the c(x) >= 0 convention."""
    if c_val.size == 0:
        return 0.0
    return float(np.sum(np.maximum(0.0, -c_val)))

def solve_slp_mvp(
    obj: ObjectiveFn, 
    con: Optional[ConstraintFn], 
    x0: ArrayLike, 
    xl: Optional[ArrayLike] = None, 
    xu: Optional[ArrayLike] = None, 
    fd_step=None,
    fd_type: str = "forward",
    maxiter: int = 50, 
    move_limit: float = 0.1,
    infeasibility_penalty: float = 1000.0,
    xtol: float = 1e-4,
    gtol: float = 1e-6,
    ftol: float = 1e-3,
    display: bool = False,
    callback: Optional[Callable[[np.ndarray], None]] = None,
    **kwargs
) -> SLPResult:
    """
    Simplest possible Sequential Linear Programming (SLP) solver.
    """
    # 1. Initialization
    x = np.array(x0, dtype=float).ravel().copy()
    n = x.size
    if xl is None:
        xl = np.full(n, -np.inf)
    else:
        xl = np.array(xl, dtype=float).ravel()
    if xu is None:
        xu = np.full(n, np.inf)
    else:
        xu = np.array(xu, dtype=float).ravel()
    
    # Finite difference step sizes
    h = _fd_step_array(fd_step, n)

    # Initial function evaluations
    f = float(obj(x))
    c = np.array(con(x), dtype=float).ravel() if con else np.array([])
    m = c.size
    nfev = 1

    def merit(f_val, c_val):
        """Objective plus an exterior L1 penalty on true nonlinear violations."""
        viol = _constraint_violation(c_val)
        return f_val + infeasibility_penalty * viol

    P = merit(f, c)

    # Simplest move limits: A box around the current point.
    delta = move_limit * (xu - xl)
    delta = np.where(np.isinf(delta), 1.0, delta)
    delta = np.clip(delta, 1e-6, np.inf) # Prevent zero move limits if bounds are tight

    history = [{
        "iter": 0,
        "f": f,
        "viol": _constraint_violation(c),
        "P": P,
        "slack_sum": 0.0,
        "slack_max": 0.0,
    }]
    
    success = False
    message = "Maximum iterations reached"

    if display:
        print("*" * 80)
        print("    SLP MVP Optimizer (With Slacks)")
        print(f"    Initial constraints: {m}")
        print("*" * 80)
        print(f"  {'iter':>5s}  {'f':>12s}  {'viol':>12s}  {'P':>12s}  {'||dx||':>12s}  {'sum(s)':>12s}  {'nfev':>6s}")
        print(f"  {0:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  {P:12.4e}  {'-':>12s}  {0.0:12.4e}  {nfev:>6d}")

    # Initial callback
    if callback:
        callback(x.copy())

    needs_gradients = True
    g = None
    J = None

    # 2. Optimization Loop
    for iter_no in range(1, maxiter + 1):
        
        if needs_gradients:
            # Calculate Gradients using Forward Difference
            g = _grad_obj_fd(obj, x, f, h, fd_type=fd_type)
            nfev += n if fd_type != "central" else 2 * n
            if m > 0:
                J = _jac_con_fd(con, x, c, h, fd_type=fd_type)
                nfev += n if fd_type != "central" else 2 * n
            else:
                J = np.zeros((0, n))
            
            needs_gradients = False

        # Setup the Linear Programming (LP) Subproblem
        if m > 0:
            # LP variables are z = [dx (n), s (m)].
            # Slack variables relax the linearized constraints but are
            # penalized, so the LP remains feasible without making
            # infeasibility free.
            c_lp = np.concatenate((g, np.full(m, infeasibility_penalty)))
            # -J*dx - s <= c(x)  <=>  c(x) + J*dx + s >= 0
            A_ub = np.hstack((-J, -np.eye(m)))
            b_ub = c.copy()
        else:
            c_lp = g
            A_ub = None
            b_ub = None

        # Set bounds for the LP variables (dx is bounded by both global bounds and local move limits)
        bounds = [(max(xl[i] - x[i], -delta[i]), min(xu[i] - x[i], delta[i])) for i in range(n)]
        if m > 0:
            # Slack variables are non-negative
            bounds.extend([(0.0, None) for _ in range(m)])

        # Solve the LP
        res = linprog(c_lp, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        
        if not res.success:
            # This should now be very rare since slacks make the LP feasible
            delta *= 0.5
            if display:
                print(f"  {iter_no:>5d}  {' ' * 12}  {' ' * 12}  {' ' * 12}  LP FAILED   {nfev:>6d}")
            
            if np.max(delta) < xtol:
                success = False
                message = f"Stopped (LP solver failed even with slacks and move limits shrunk < xtol: {xtol:.2e})."
                break
            continue
            
        # Extract the proposed step for design variables (first n components)
        dx = res.x[:n]
        slack = res.x[n:] if m > 0 else np.array([])
        slack_sum = float(np.sum(slack)) if slack.size else 0.0
        slack_max = float(np.max(slack)) if slack.size else 0.0
        deltanorm = np.linalg.norm(dx)
        
        # Evaluate the proposed step
        x_new = np.clip(x + dx, xl, xu)
        f_new = float(obj(x_new))
        c_new = np.array(con(x_new), dtype=float).ravel() if con else np.array([])
        nfev += 1
        P_new = merit(f_new, c_new)
        
        # Accept/Reject Step and Update Move Limits
        if P_new < P:
            # Improvement found -> Accept the step
            x = x_new
            f = f_new
            c = c_new
            P = P_new
            needs_gradients = True
            
            history.append({
                "iter": iter_no,
                "f": float(f),
                "viol": _constraint_violation(c),
                "P": float(P),
                "step": float(deltanorm),
                "slack_sum": slack_sum,
                "slack_max": slack_max,
            })

            if display:
                print(f"  {iter_no:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  {P:12.4e}  {deltanorm:12.4e}  {slack_sum:12.4e}  {nfev:>6d}")
            
            # Callback after accepting step
            if callback:
                callback(x.copy())

            # Simple stopping criteria: Converged if step size is very small
            if deltanorm < xtol:
                success = True
                message = f"Converged successfully (Step size < xtol: {xtol:.2e})."
                break
                
            # Slightly expand move limits on success
            delta = np.minimum(delta * 1.2, move_limit * np.where(np.isinf(xu - xl), 1.0, xu - xl))
        else:
            # No improvement -> Reject the step and shrink move limits
            delta *= 0.5
            
            history.append({
                "iter": iter_no,
                "f": float(f),
                "viol": _constraint_violation(c),
                "P": float(P),
                "step": float(deltanorm),
                "slack_sum": slack_sum,
                "slack_max": slack_max,
            })
            if display:
                print(f"  {iter_no:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  {P:12.4e}  REJECTED      {slack_sum:12.4e}  {nfev:>6d}")
            
            # Simple stopping criteria: Converged if move limits become too small
            if np.max(delta) < xtol:
                success = False
                message = f"Stopped (Move limits shrunk < xtol: {xtol:.2e})."
                break

    return SLPResult(
        x=x,
        fun=f,
        nit=iter_no,
        nfev=nfev,
        success=success,
        message=message,
        history=history
    )

if __name__ == "__main__":
    # Simple smoke test
    def obj(x):
        return (x[0] - 1)**2 + (x[1] - 2)**2

    def con(x):
        return np.array([
            4 - x[0]**2 - x[1]**2,  # x0^2 + x1^2 <= 4
            x[0] + x[1] - 1         # x0 + x1 >= 1
        ])

    res = solve_slp_mvp(
        obj=obj,
        con=con,
        x0=[0.5, 0.5],
        xl=[-5, -5],
        xu=[5, 5],
        maxiter=50,
        display=True
    )
    print(f"x* = {res.x}")
    print(f"f* = {res.fun}")
    print(f"msg = {res.message}")
