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

    def constraint_values(x):             # This is the function that is used to calculate the constraint values
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

    constraints = [                 # This is the constraints list that is used to set the constraints for the optimization process
        {"type": "ineq", "fun": lambda x, i=i: constraint_values(x)[i]}
        for i in range(len(constraint_names))
    ]

    options = {                   # This is the options dictionary that is used to set the optimization options
        "disp": True,
        "eps": 0.1,
        "ftol": 1e-3,
        "maxiter": 20
    }

    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="SLSQP",
        options=options,
        save_folder=save_folder,
    )

    def objective(x):                # This is the objective function that is used to optimize the design variables
        val = RunAPDL(SWcoor, x, Misc)
        logger.log_evaluation(x, val)
        return val

    def callback(xk):               # This is the callback function that is called after each iteration in the optimization process
        logger.log_iteration(xk)

    result = spo.minimize(               # This is the optimization function that is used to optimize the design variables
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        callback=callback,
        options=options
    )

    logger.finalize(result)
    return result, logger.txt_path, logger.csv_path