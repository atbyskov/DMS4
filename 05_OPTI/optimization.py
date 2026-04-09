import os

import numpy as np
from scipy import optimize as spo
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

    def constraint_values(x):             # This is the function that is used to calculate the constraint values
        #RunAPDL(mapdl, SWcoor, x, Misc)
        utils = PostProcessor()
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

    
        """
        Return all inequality constraints in the form c(x) >= 0.
        Add as many as you want here.
        """
        return np.array([
            x[1] - eps_geom, # Column Thickness [mm]
            x[3] - eps_geom, # Brace Thickness [mm]
            1.0 - Util_LB_values_col,           # local buckling column
            1.0 - Util_LB_values_brace,           # local buckling brace
            1.0 - Util_NF_values_col,           # normal force column
            1.0 - Util_NF_values_brace,           # normal force brace
            1.0 - Util_S_values_col,           # shear column
            1.0 - Util_S_values_brace,           # shear brace
            1.0 - Util_T_values_col,           # Torsion column
            1.0 - Util_T_values_brace,           # Torsion brace 
            1.0 - Util_BNS_values_col,           # bending, normal and shear column
            1.0 - Util_BNS_values_brace,           # bending, normal and shear brace
            1.0 - Util_BR_values_col,           # Flexural and torsional buckling column
            1.0 - Util_BR_values_brace,           # Flexural and torsional buckling brace
            1.0 - Util_IN_values_col,           # Interaction column
            1.0 - Util_IN_values_brace,           # Interaction brace
            1.0 - Util_BS_values_brace,           # Brace-Step c(x) = sigma_vm/f_y_brace >= 0 "Inequality Constraint" 
            Util_Class_2_values_col,            # Class 2 column c(x) = 70*235/f_y-d0/t0 >= 0 "Inequality Constraint"
            Util_Class_2_values_brace,            # Class 2 brace c(x) = 70*235/f_y_brace-d1/t1 >= 0 "Inequality Constraint"
            Eigenvalue_1_values,            # Eigenvalue 1 c(x)=4.0-alpha_cr >= 0 "Inequality Constraint"
            # add more constraints here later if needed
        ], dtype=float)
    constraint_names = [
        "thickness_column", 
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
        "Brace_Step",
        "Class_2_column",
        "Class_2_brace",
        "Eigenvalue_1_constraint",
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

        utils = PostProcessor()
        Util_list = utils.Util_list(x, Misc)  

        print("\n--- UTILIZATION REPORT ---")
        for key, util in Util_list.items():

            if util is None:
                print(f"{key:10s}  Column:   N/A    Brace:   N/A")
                continue

            util = np.atleast_1d(util)

            col_val = util[0] if len(util) > 0 and np.isfinite(util[0]) else np.nan
            brc_val = util[1] if len(util) > 1 and np.isfinite(util[1]) else np.nan

            print(
                f"{key:10s}  "
                f"Column: {col_val:8.4f}   "
                f"Brace: {brc_val:8.4f}"
            )
        print("------------------")
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

