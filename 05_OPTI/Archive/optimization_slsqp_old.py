import os

import numpy as np
from scipy import optimize as spo
from pyslsqp import optimize
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor


def run_optimization(
    mapdl,
    var,
    SWcoor,
    Misc,
    eps_geom=0.1,
    save_folder="Optimization_Logs",
):
    x0 = np.array(var, dtype=float)

    bounds = [
        (40.0, 100.0),  # Column Outer Diameter [mm]
        (1.0, 7.0),  # Column Thickness [mm]
        (10.0,  100.0),  # Brace Outer Diameter [mm]
        (0.1,  7.0)   # Brace Thickness [mm]
    ]
    

    def constraint_values(x):
        val = RunAPDL(mapdl, SWcoor, x, Misc)
        logger.log_evaluation(x, val)
        utils = PostProcessor()

        def as_1d_float_array(v):
            if hasattr(v, "to_numpy"):  # pandas Series / DataFrame
                return v.to_numpy(dtype=float).ravel()
            return np.asarray(v, dtype=float).ravel()      

        Util_LB_values_col, Util_LB_values_brace = utils.Util_LB(x, Misc)
        Util_NF_values_col, Util_NF_values_brace = utils.Util_NF(x, Misc)
        Util_S_values_col, Util_S_values_brace = utils.Util_S(x, Misc)
        Util_T_values_col, Util_T_values_brace = utils.Util_T(x, Misc)
        Util_BNS_values_col, Util_BNS_values_brace = utils.Util_BNS(x, Misc)
        Util_BR_values_col, Util_BR_values_brace = utils.Util_BR(x, Misc)
        Util_IN_values_col, Util_IN_values_brace = utils.Util_IN(x, Misc)
        Util_BS_values_brace = utils.Util_BS(x, Misc)
        Util_Class_2_values_col, Util_Class_2_values_brace = utils.Class_2(x, Misc)
        Eigenvalue_1_values = utils.Eigenvalue_1()

        c = np.concatenate([
            np.array([x[1] - eps_geom], dtype=float),
            np.array([x[3] - eps_geom], dtype=float),

            1.0 - as_1d_float_array(Util_LB_values_col),
            1.0 - as_1d_float_array(Util_LB_values_brace),

            1.0 - as_1d_float_array(Util_NF_values_col),
            1.0 - as_1d_float_array(Util_NF_values_brace),

            1.0 - as_1d_float_array(Util_S_values_col),
            1.0 - as_1d_float_array(Util_S_values_brace),

            1.0 - as_1d_float_array(Util_T_values_col),
            1.0 - as_1d_float_array(Util_T_values_brace),

            1.0 - as_1d_float_array(Util_BNS_values_col),
            1.0 - as_1d_float_array(Util_BNS_values_brace),

            1.0 - as_1d_float_array(Util_BR_values_col),
            1.0 - as_1d_float_array(Util_BR_values_brace),

            1.0 - as_1d_float_array(Util_IN_values_col),
            1.0 - as_1d_float_array(Util_IN_values_brace),

            1.0 - as_1d_float_array(Util_BS_values_brace),

            as_1d_float_array(Util_Class_2_values_col),
            as_1d_float_array(Util_Class_2_values_brace),

            as_1d_float_array(Eigenvalue_1_values),
        ])
        print(len(c))   # Print the length of the constraint vector
        return c


    constraints = [
        {"type": "ineq", "fun": constraint_values}
    ]

    options = {                   # This is the options dictionary that is used to set the optimization options
        "disp": True,
        "eps": 0.01,
        "ftol": 1e-3,
        "maxiter": 40
    }

    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="SLSQP",
        options=options,
        save_folder=save_folder,
    )

    def objective(x):                # This is the objective function that is used to optimize the design variables
        val = RunAPDL(mapdl, SWcoor, x, Misc)
        logger.log_evaluation(x, val)

        #utils = PostProcessor()


        return val

    def callback(xk):               # This is the callback function that is called after each iteration in the optimization process
        logger.log_iteration(xk)


    result = spo.minimize (               # This is the optimization function that is used to optimize the design variables 
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        callback=callback,
        options=options
    )
    #result = optimize(x0, objective, constraints=constraints, bounds=bounds, callback=callback)

    logger.finalize(result)
    return result, logger.txt_path, logger.csv_path

