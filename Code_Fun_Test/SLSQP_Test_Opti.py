import numpy as np
from scipy import optimize as spo


def run_slsqp_demo():
    """
    A teaching example for following SLSQP behavior in detail.

    What this script shows:
    - Every objective evaluation
    - Every constraint evaluation
    - Small perturbation steps used for finite-difference gradients
    - Larger candidate moves that look like line-search / trial steps
    - Accepted major iterations via callback

    Important:
    SciPy does NOT expose SLSQP's internal line-search alpha directly.
    But by printing every evaluated x, you can still follow the behavior closely.
    """

    # ============================================================
    # Global counters and state for logging
    # ============================================================
    state = {
        "obj_eval": 0,
        "cb_eval": 0,
        "con_eval": {},
        "last_obj_x": None,
        "accepted_iterates": []
    }

    # ------------------------------------------------------------
    # Helper: classify step type
    # ------------------------------------------------------------
    def classify_step(x, prev_x, eps_small=1e-7, eps_fd=1e-5):
        """
        Roughly classify what kind of step this looks like.

        This is only heuristic, but very useful for understanding behavior.
        """
        if prev_x is None:
            return "initial evaluation", None

        dx = x - prev_x
        step_norm = np.linalg.norm(dx)

        if step_norm == 0:
            return "same point repeated", step_norm
        elif step_norm < eps_small:
            return "tiny repeat / numerical noise", step_norm
        elif step_norm < eps_fd:
            return "very small perturbation (likely finite-difference probe)", step_norm
        elif step_norm < 1e-2:
            return "small local probe", step_norm
        else:
            return "larger trial step / candidate move", step_norm

    # ------------------------------------------------------------
    # Objective function
    # ------------------------------------------------------------
    def objective(x):
        state["obj_eval"] += 1

        prev_x = state["last_obj_x"]
        step_type, step_norm = classify_step(x, prev_x)

        # A smooth nonlinear objective
        # Minimum is somewhere near (2, 1), but constraints modify it
        f = (
            (x[0] - 2.0) ** 2
            + 2.0 * (x[1] - 1.0) ** 2
            + 0.3 * np.sin(3.0 * x[0])
            + 0.2 * np.cos(4.0 * x[1])
        )

        print("\n" + "=" * 90)
        print(f"[OBJECTIVE EVALUATION {state['obj_eval']}]")
        print(f"x = {np.array2string(x, precision=10, suppress_small=False)}")

        if prev_x is None:
            print("previous x = None")
            print("step type  = initial evaluation")
        else:
            dx = x - prev_x
            print(f"previous x = {np.array2string(prev_x, precision=10, suppress_small=False)}")
            print(f"dx         = {np.array2string(dx, precision=10, suppress_small=False)}")
            print(f"step norm  = {step_norm:.12e}")
            print(f"step type  = {step_type}")

        print(f"f(x)       = {f:.12f}")

        state["last_obj_x"] = x.copy()
        return f

    # ------------------------------------------------------------
    # Constraint wrapper
    # SciPy inequality constraints must satisfy c(x) >= 0
    # ------------------------------------------------------------
    def make_constraint(name, func):
        state["con_eval"][name] = 0

        def wrapped(x):
            state["con_eval"][name] += 1
            value = func(x)

            print("-" * 90)
            print(f"[CONSTRAINT '{name}' EVALUATION {state['con_eval'][name]}]")
            print(f"x          = {np.array2string(x, precision=10, suppress_small=False)}")
            print(f"value      = {value:.12f}")
            print(f"feasible   = {value >= 0.0}")

            return value

        return {"type": "ineq", "fun": wrapped}

    # ------------------------------------------------------------
    # Callback: called after a major accepted iterate
    # ------------------------------------------------------------
    def callback(xk):
        state["cb_eval"] += 1
        state["accepted_iterates"].append(xk.copy())

        print("\n" + "#" * 90)
        print(f"[CALLBACK / ACCEPTED ITERATE {state['cb_eval']}]")
        print(f"xk = {np.array2string(xk, precision=10, suppress_small=False)}")

        if len(state["accepted_iterates"]) >= 2:
            dx = state["accepted_iterates"][-1] - state["accepted_iterates"][-2]
            print(f"change from previous accepted iterate = {np.array2string(dx, precision=10, suppress_small=False)}")
            print(f"accepted step norm                   = {np.linalg.norm(dx):.12e}")

        print("#" * 90)

    # ============================================================
    # Problem definition
    # ============================================================

    # Start point
    x0 = np.array([0.5, 2.8], dtype=float)

    # Bounds
    bounds = [
        (-1.0, 4.0),   # x0
        (0.0, 4.0),    # x1
    ]

    # Constraints: c(x) >= 0
    constraints = [
        # x0 + x1 >= 2
        make_constraint(
            "x0_plus_x1_ge_2",
            lambda x: x[0] + x[1] - 2.0
        ),

        # x0 <= 3.2  ->  3.2 - x0 >= 0
        make_constraint(
            "x0_le_3p2",
            lambda x: 3.2 - x[0]
        ),

        # nonlinear curved feasible region:
        # x1 >= 0.2 + 0.25*(x0-1.5)^2
        make_constraint(
            "above_parabola",
            lambda x: x[1] - (0.2 + 0.25 * (x[0] - 1.5) ** 2)
        ),

        # x1 <= 3.5 - 0.3*x0
        make_constraint(
            "below_line",
            lambda x: (3.5 - 0.3 * x[0]) - x[1]
        ),
    ]

    # ============================================================
    # Print problem setup
    # ============================================================
    print("=" * 90)
    print("SLSQP LINE-SEARCH / TRIAL-STEP DEMO")
    print("=" * 90)
    print("Objective:")
    print("    f(x0, x1) = (x0-2)^2 + 2*(x1-1)^2 + 0.3*sin(3*x0) + 0.2*cos(4*x1)")
    print("\nInitial guess:")
    print(f"    x0 = {x0}")
    print("\nBounds:")
    for i, b in enumerate(bounds):
        print(f"    x[{i}] in {b}")

    print("\nInequality constraints c(x) >= 0:")
    print("    1) x0 + x1 - 2 >= 0")
    print("    2) 3.2 - x0 >= 0")
    print("    3) x1 - (0.2 + 0.25*(x0-1.5)^2) >= 0")
    print("    4) (3.5 - 0.3*x0) - x1 >= 0")

    print("\nWhat to watch in the output:")
    print("    - very small changes in x: usually finite-difference gradient probes")
    print("    - larger candidate moves: likely trial steps from the QP subproblem / line search")
    print("    - callback hits: accepted major iterates")
    print("=" * 90)

    # ============================================================
    # Run optimizer
    # ============================================================
    result = spo.minimize(
        fun=objective,
        x0=x0,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        callback=callback,
        options={
            "disp": True,
            "maxiter": 50,
            "ftol": 1e-9,
            "eps": 1e-6,   # finite-difference perturbation size
        }
    )

    # ============================================================
    # Final summary
    # ============================================================
    print("\n" + "=" * 90)
    print("FINAL RESULT")
    print("=" * 90)
    print(f"success      = {result.success}")
    print(f"message      = {result.message}")
    print(f"nit          = {result.nit}")
    print(f"nfev         = {result.nfev}")
    print(f"x*           = {np.array2string(result.x, precision=12, suppress_small=False)}")
    print(f"f(x*)        = {result.fun:.12f}")

    print("\nConstraint values at final solution:")
    for con in constraints:
        value = con["fun"](result.x)
        print(f"    {value:.12f}   feasible = {value >= 0.0}")

    print("\nAccepted iterates:")
    if len(state["accepted_iterates"]) == 0:
        print("    No callback iterates recorded.")
    else:
        for i, xk in enumerate(state["accepted_iterates"], start=1):
            print(f"    Iter {i}: {np.array2string(xk, precision=10, suppress_small=False)}")

    print("\nTotal counts:")
    print(f"    Objective evaluations : {state['obj_eval']}")
    print(f"    Callback calls        : {state['cb_eval']}")
    for name, count in state["con_eval"].items():
        print(f"    Constraint '{name}' evaluations: {count}")

    return result


if __name__ == "__main__":
    run_slsqp_demo()