"""
Minimum Viable Product (MVP) of a Sequential Linear Programming (SLP) solver.

This implementation focuses on maximum readability and simplicity.
It uses:
- Forward difference gradients (reused from slp.py)
- A basic LP subproblem with strict linearized constraints
- The simplest possible move limits (trust region) that shrink on rejection
- Simple stopping criteria (step size / move limit size)
"""

import numpy as np
from scipy.optimize import linprog

# Reuse the gradient and finite difference helpers from the main slp.py
from .slp import _grad_obj_fd, _jac_con_fd, _fd_step_array, SLPResult

def solve_slp_mvp(
    obj, 
    con, 
    x0, 
    xl=None, 
    xu=None, 
    fd_step=None,
    fd_type="forward",
    maxiter=50, 
    move_limit=0.1,
    infeasibility_penalty=1000.0,
    xtol=1e-4,
    gtol=1e-6,
    ftol=1e-3,
    display=False,
    **kwargs
):
    """
    Simplest possible Sequential Linear Programming (SLP) solver.
    
    Parameters:
        obj: Callable returning the objective value (scalar).
        con: Callable returning inequality constraints (numpy array). c(x) >= 0 is feasible.
        x0: Initial guess for design variables.
        xl, xu: Lower and upper bounds for design variables.
        maxiter: Maximum number of iterations.
        move_limit: Fraction of the bound width used for initial move limits.
        infeasibility_penalty: Penalty factor for constraint violations in the merit function.
        xtol: Stopping tolerance for the step size.
        display: Whether to print progress to console.
        
    Returns:
        SLPResult: A dataclass containing x, fun, nit, nfev, success, message, history.
    """
    # 1. Initialization
    x = np.array(x0, dtype=float).copy()
    n = len(x)
    if xl is None:
        xl = np.full(n, -np.inf)
    else:
        xl = np.array(xl, dtype=float)
    if xu is None:
        xu = np.full(n, np.inf)
    else:
        xu = np.array(xu, dtype=float)
    
    # Finite difference step sizes (using helper from slp.py)
    h = _fd_step_array(fd_step, n)

    # Initial function evaluations
    f = float(obj(x))
    c = np.array(con(x), dtype=float) if con else np.array([])
    m = len(c)
    nfev = 1

    def merit(f_val, c_val):
        """Merit function to evaluate step quality: Objective + Penalty * Violations"""
        viol = np.sum(np.maximum(0.0, -c_val)) if len(c_val) > 0 else 0.0
        return f_val + infeasibility_penalty * viol

    P = merit(f, c)

    # Simplest move limits: A box around the current point.
    # Initialized as a fraction of the total variable bounds.
    delta = move_limit * (xu - xl)
    delta = np.where(np.isinf(delta), 1.0, delta)
    delta = np.clip(delta, 1e-6, np.inf) # Prevent zero move limits if bounds are tight

    history = [{"iter": 0, "f": f, "viol": float(np.sum(np.maximum(0.0, -c)) if m > 0 else 0.0), "P": P}]
    
    success = False
    message = "Maximum iterations reached"

    if display:
        print("*" * 80)
        print("    SLP MVP Optimizer (No Slacks)")
        print("*" * 80)
        print(f"  {'iter':>5s}  {'f':>12s}  {'viol':>12s}  {'P':>12s}  {'||dx||':>12s}  {'nfev':>6s}")
        print(f"  {0:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  {P:12.4e}  {'-':>12s}  {nfev:>6d}")

    # 2. Optimization Loop
    for iter_no in range(1, maxiter + 1):
        
        # Calculate Gradients using Forward Difference
        g = _grad_obj_fd(obj, x, f, h, fd_type=fd_type)
        nfev += n if fd_type != "central" else 2 * n
        if m > 0:
            J = _jac_con_fd(con, x, c, h, fd_type=fd_type)
            nfev += n if fd_type != "central" else 2 * n
        else:
            J = np.zeros((0, n))

        # Setup the Linear Programming (LP) Subproblem
        # Variables for LP: [dx_0, ..., dx_{n-1}]
        c_lp = g
        
        # Subject to linearized constraints: c(x) + J * dx >= 0  =>  -J * dx <= c(x)
        if m > 0:
            A_ub = -J
            b_ub = c.copy()
        else:
            A_ub = None
            b_ub = None

        # Set bounds for the LP variables (dx is bounded by both global bounds and local move limits)
        bounds = [(max(xl[i] - x[i], -delta[i]), min(xu[i] - x[i], delta[i])) for i in range(n)]

        # Solve the LP
        res = linprog(c_lp, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        
        if not res.success:
            # LP solver failed (likely infeasible) -> shrink move limits and try again in next iteration
            delta *= 0.5
            if display:
                print(f"  {iter_no:>5d}  {' ' * 12}  {' ' * 12}  {' ' * 12}  INFEASIBLE  {nfev:>6d}")
            
            if np.max(delta) < xtol:
                success = False
                message = f"Stopped (LP infeasible and move limits shrunk < xtol: {xtol:.2e})."
                break
            continue
            
        # Extract the proposed step for design variables
        dx = res.x
        deltanorm = np.linalg.norm(dx)
        
        # Evaluate the proposed step
        x_new = np.clip(x + dx, xl, xu)
        f_new = float(obj(x_new))
        c_new = np.array(con(x_new), dtype=float) if con else np.array([])
        nfev += 1
        P_new = merit(f_new, c_new)
        
        # Accept/Reject Step and Update Move Limits
        if P_new < P:
            # Improvement found -> Accept the step
            x = x_new
            f = f_new
            c = c_new
            P = P_new
            
            history.append({
                "iter": iter_no,
                "f": float(f),
                "viol": float(np.sum(np.maximum(0.0, -c)) if m > 0 else 0.0),
                "P": float(P),
                "step": float(deltanorm)
            })

            if display:
                print(f"  {iter_no:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  {P:12.4e}  {deltanorm:12.4e}  {nfev:>6d}")
            
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
                "viol": float(np.sum(np.maximum(0.0, -c)) if m > 0 else 0.0),
                "P": float(P),
                "step": float(deltanorm)
            })
            if display:
                print(f"  {iter_no:>5d}  {f:12.4e}  {history[-1]['viol']:12.4e}  {P:12.4e}  REJECTED  {nfev:>6d}")
            
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
