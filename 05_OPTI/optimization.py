import numpy as np
from scipy import optimize as spo
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger


def run_optimization(var, SWcoor, Misc, eps_geom=1, save_folder="Optimization_Logs"):
    x0 = np.array(var, dtype=float)

    bounds = [
        (30.0, 60.0),  # R0
        (31.0, 70.0),  # R1
        (8.0,  20.0),  # R2
        (9.0,  25.0)   # R3
    ]

    def constraint_values(x):
        """
        Return all inequality constraints in the form c(x) >= 0.
        Add as many as you want here.
        """
        return np.array([
            x[1] - x[0] - eps_geom,
            x[3] - x[2] - eps_geom,
            # add more constraints here later if needed
        ], dtype=float)

    constraint_names = [
        "thickness_col",
        "thickness_brace",
        # add more names here if you add more constraints
    ]

    constraints = [
        {"type": "ineq", "fun": lambda x, i=i: constraint_values(x)[i]}
        for i in range(len(constraint_names))
    ]

    options = {
        "disp": True,
        "eps": 0.1,
        "ftol": 1e-3,
        "maxiter": 20
    }

    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        constraint_func=constraint_values,
        constraint_names=constraint_names,
        method="SLSQP",
        options=options,
        save_folder=save_folder
    )

    def objective(x):
        val = RunAPDL(SWcoor, x, Misc)
        logger.log_evaluation(x, val)
        return val

    def callback(xk):
        logger.log_iteration(xk)

    result = spo.minimize(
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        callback=callback,
        options=options
    )

    logger.finalize(result)
    return result, logger.history, logger.txt_path