import numpy as np
from scipy import optimize as spo
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import Util_LC, Util_NF, Util_S, Util_T, Util_BNS, Util_BR, Util_IN


def run_optimization(var, SWcoor, Misc, eps_geom=1, save_folder="Optimization_Logs"):
    x0 = np.array(var, dtype=float)

    bounds = [
        (30.0, 60.0),  # R0
        (31.0, 70.0),  # R1
        (8.0,  20.0),  # R2
        (9.0,  25.0)   # R3
    ]

    def constraint_values(x):             # This is the function that is used to calculate the constraint values
        Util_LC_values = Util_LC(x, Misc)
        Util_NF_values = Util_NF(x, Misc)
        Util_S_values = Util_S(x, Misc)
        Util_T_values = Util_T(x, Misc)
        Util_BNS_values = Util_BNS(x, Misc)
        Util_BR_values = Util_BR(x, Misc)
        Util_IN_values = Util_IN(x, Misc)
        """
        Return all inequality constraints in the form c(x) >= 0.
        Add as many as you want here.
        """
        return np.array([
            x[1] - x[0] - eps_geom,
            x[3] - x[2] - eps_geom,
            1.0 - Util_LC_values[0],           # local buckling column
            1.0 - Util_LC_values[1],           # local buckling brace
            1.0 - Util_NF_values[0],           # normal force column
            1.0 - Util_NF_values[1],           # normal force brace
            1.0 - Util_S_values[0],           # shear column
            1.0 - Util_S_values[1],           # shear brace
            1.0 - Util_T_values[0],           # Torsion column
            1.0 - Util_T_values[1],           # Torsion brace
            1.0 - Util_BNS_values[0],           # bending, normal and shear column
            1.0 - Util_BNS_values[1],           # bending, normal and shear brace
            1.0 - Util_BR_values[0],           # Flexural and torsional buckling brace
            1.0 - Util_BR_values[1],           # Flexural and torsional buckling brace
            1.0 - Util_IN_values[0],           # Interaction column
            1.0 - Util_IN_values[1],           # Interaction brace
            # add more constraints here later if needed
        ], dtype=float)

    constraint_names = [
        "thickness_col",
        "thickness_brace",
        "local_buckling_column",
        "local_buckling_brace",
        "normal_force_column",
        "normal_force_brace",
        "shear_column",
        "shear_brace",
        "torsion_column",
        "torsion_brace",
        "bending_normal_shear_column",  # Bending, Normal and Shear column
        "bending_normal_shear_brace",
        "flexural_torsional_buckling_column",  # Flexural and Torsional buckling column
        "flexural_torsional_buckling_brace",
        "interaction_column",  # Interaction column
        "interaction_brace",  # Interaction brace
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