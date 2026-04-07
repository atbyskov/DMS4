import numpy as np
from scipy import optimize as spo
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor


def run_optimization(var, SWcoor, Misc, eps_geom=1, save_folder="Optimization_Logs"):
    x0 = np.array(var, dtype=float)

    bounds = [
        (10.0, 60.0),  # R0
        (11.0, 70.0),  # R1
        (1.0,  20.0),  # R2
        (2.0,  25.0)   # R3
    ]

    def constraint_values(x):
        """
        Return all inequality constraints in the form c(x) >= 0.

        To add a new constraint later, just add one more line
        to the list below.
        """
        utils = PostProcessor()

        Util_LB_values  = utils.Util_LB(x, Misc)
        Util_NF_values  = utils.Util_NF(x, Misc)
        Util_S_values   = utils.Util_S(x, Misc)
        Util_T_values   = utils.Util_T(x, Misc)
        Util_BNS_values = utils.Util_BNS(x, Misc)
        Util_BR_values  = utils.Util_BR(x, Misc)
        Util_IN_values  = utils.Util_IN(x, Misc)

        constraint_list = [
            lambda: x[1] - x[0] - eps_geom,      # thickness column
            lambda: x[3] - x[2] - eps_geom,      # thickness brace

            lambda: 1.0 - Util_LB_values[0],     # local buckling column
            lambda: 1.0 - Util_LB_values[1],     # local buckling brace

            lambda: 1.0 - Util_NF_values[0],     # normal force column
            lambda: 1.0 - Util_NF_values[1],     # normal force brace

            lambda: 1.0 - Util_S_values[0],      # shear column
            lambda: 1.0 - Util_S_values[1],      # shear brace

            lambda: 1.0 - Util_T_values[0],      # torsion column
            lambda: 1.0 - Util_T_values[1],      # torsion brace

            lambda: 1.0 - Util_BNS_values[0],    # bending, normal and shear column
            lambda: 1.0 - Util_BNS_values[1],    # bending, normal and shear brace

            lambda: 1.0 - Util_BR_values[0],     # flexural / torsional buckling column
            lambda: 1.0 - Util_BR_values[1],     # flexural / torsional buckling brace

            lambda: 1.0 - Util_IN_values[0],     # interaction column
            lambda: 1.0 - Util_IN_values[1],     # interaction brace

            # Add future constraints here, example:
            # lambda: 3.0 - some_displacement_value,
            # lambda: 1.0 - utils.Util_SOMETHING(x, Misc)[0],
            # lambda: 1.0 - utils.Util_SOMETHING(x, Misc)[1],
        ]

        return np.array([c() for c in constraint_list], dtype=float)

    n_constraints = len(constraint_values(x0))

    constraints = [
        {"type": "ineq", "fun": lambda x, i=i: constraint_values(x)[i]}
        for i in range(n_constraints)
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
        method="SLSQP",
        options=options,
        save_folder=save_folder,
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
    return result, logger.txt_path, logger.csv_path