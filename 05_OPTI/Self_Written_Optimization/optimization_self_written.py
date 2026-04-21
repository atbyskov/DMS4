import os
import numpy as np
from pyslsqp import optimize as pyslsqp_optimize

from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor

# Optimization Function
def run_optimization(mapdl,var,Misc):
    
    # Data saved in this folder:
    save_folder="Optimization_Logs"

    # Design Variables and Bounds
    des_var = [
        ("d0", (40.0,  100)),
        ("t0", (1.0,   7.0)),
        ("d1", (10.0,  100.0)),
        ("t1", (0.1,   7.0)),
        ("rad",(150.0, 350.0)),
    ]
    
    # Set active variables - handle both simple and complex dict formats
    active = []
    for name, default_bnds in des_var:
        if name in var:
            v_data = var[name]
            if isinstance(v_data, dict):
                if v_data.get("active", True):
                    # Use bounds from var dict if available, otherwise use default
                    active.append((name, v_data.get("bounds", default_bnds)))
            else:
                active.append((name, default_bnds))
    
    names = [name for name, _ in active]
    Misc["active_vars"] = names

    # Print active variables
    print(f"{len(names)} Design Variables included: " +
            ", ".join(f"x{i}={n}" for i,n in enumerate(names)))

    # Set initial guess, bounds and ranges to np array
    x0_list = []
    for name, _ in active:
        v_data = var[name]
        if isinstance(v_data, dict):
            x0_list.append(v_data["value"])
        else:
            x0_list.append(v_data)
    x0 = np.array(x0_list, dtype=float)

    bounds = [bnds for _, bnds in active]
    xl, xu = np.array(bounds).T

    # Minimum thickness specification
    eps_geom=0.1

    # Step Options for each variable
    fd_step_options = {
        "d0": 0.01,
        "t0": 0.01,
        "d1": 0.01,
        "t1": 0.01,
        "rad": 2
        }
    fd_step = [fd_step_options[name] for name, _ in active]

    # Maximum objective function tolerance and iterations
    acc = 1e-3
    maxiter = 100

    # Logger Options
    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="PySLSQP",
        options={
            "finite_diff_abs_step": fd_step,
            "acc": acc,
            "maxiter": maxiter,
        },
        save_folder=save_folder,
    )

    # Internal cache so RunAPDL is only executed once per unique x
    cache = {"x": None, "f": None, "c": None,}

    # Prepare initial full var dict (complex format) for use in model evaluation
    initial_full_var = {}
    for name, default_bnds in des_var:
        v_data = var.get(name)
        if v_data is not None:
            if isinstance(v_data, dict):
                initial_full_var[name] = v_data.copy()
            else:
                initial_full_var[name] = {"value": v_data, "active": True}
    
    # Helper
    def arr(v):
        return np.asarray(v, dtype=float).ravel()

    # Model evaluation
    def evaluate_model(x):
        # Read variables and check
        x = arr(x)
        if cache["x"] is not None and np.array_equal(x, cache["x"]):
            return cache["f"], cache["c"]
        
        # Construct full complex var_dict for RunAPDL and simple var_dict for PostProcessor
        full_var_dict = {k: v.copy() for k, v in initial_full_var.items()}
        for name, val in zip(names, x):
            full_var_dict[name]["value"] = float(val)
        
        simple_var_dict = {k: v["value"] for k, v in full_var_dict.items()}

        # Run the Model - passing full_var_dict ensures InputFun has everything it needs
        f = RunAPDL(mapdl, full_var_dict, Misc)
        logger.log_evaluation(x, f)

        # Initiate Post Processing
        pp = PostProcessor(simple_var_dict,Misc)

        # Call Utlization Ratios
        utils = [
            pp.Util_LB(),
            pp.Util_NF(),
            pp.Util_S(),
            pp.Util_T(),
            pp.Util_BNS(),
            pp.Util_BR(),
            pp.Util_IN(),
        ]
        col_brace = [(1-arr(c), 1-arr(b)) for c,b in utils]

        # Set up the constraints
        constraints = []
        # Minimum thickness constraints (only if they are design variables)
        for t_name in ["t0", "t1"]:
            if t_name in names:
                constraints.append(arr(x[names.index(t_name)] - eps_geom))
        
        # Add remaining constraints
        constraints.extend([
            *[v for pair in col_brace for v in pair],
            1 - arr(pp.Util_BS()),
            *map(arr, pp.Class_2()),
            arr(pp.Eigenvalue_1())
        ])
        # Collect them together
        c = np.concatenate(constraints)
        # Print length of constraints
        print(f"Constraint vector length: {len(c)}")

        # Update design variables and constraints
        cache.update(x=x.copy(), f=float(f), c=c)

        return f, c
    
    # Objective function call that only returns mass
    def objective(x):
        """Objective function for PySLSQP."""
        f_val, _ = evaluate_model(x) 
        return f_val 

    # Constrain function call that only returns constraints
    def constraints(x):
        """Constraint vector for PySLSQP."""
        _, c_val = evaluate_model(x) 
        return c_val 

    # Create Folder
    os.makedirs(save_folder, exist_ok=True) 
    save_filename = os.path.join(save_folder, "pyslsqp_history.hdf5") # Save File
    summary_filename = os.path.join(save_folder, "slsqp_summary.out") # Summary File

    # Run custom optimizer (see optimization_methods.minimize)
    #   method      : "steepest_descent" | "conjugate_gradient" | "slp"
    #   line_search : "armijo" | "golden" | "bisection" | "quadratic"
    from .optimization_methods import minimize
    result = minimize(
        obj=objective,
        con=constraints,
        x0=x0,
        xl=xl,
        xu=xu,
        method="slp",
        line_search="golden",
        fd_step=fd_step,
        fd_type="forward",
        maxiter=maxiter,
        ftol=acc,
        penalty_weight=1e3,
        display=True,
    )
    # Touch save_filename/summary_filename so unused-variable linters stay quiet
    _ = (save_filename, summary_filename)
    # return the result, the text path, and the csv path
    return result, logger.txt_path, logger.csv_path 