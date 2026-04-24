# Import packages
import os
import numpy as np
from pyslsqp import optimize as pyslsqp_optimize

# Import Functions
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor

# Optimization Function
def run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings):

    # Set active variables
    active = [(name, data["bounds"]) for name, data in var.items() if data.get("active", True)]
    names = [name for name,_ in active]
    Misc["active_vars"] = names

    # Print active variables
    print(f"{len(names)} Design Variables included: " +
            ", ".join(f"x{i}={n}" for i,n in enumerate(names)), flush=True)

    # Set initial guess, bounds and ranges to np array
    x0 = np.array([var[name]["value"] for name, _ in active], dtype=float)
    bounds = [bnds for _, bnds in active]
    xl, xu = np.array(bounds).T

    # Step Options for each variable (uniform step size for all, except rad)
    fd_step_options = {
        name: 0.01 for name, _ in active
    }
    # Set separate step size for 'rad' if it's an active variable
    if "rad" in fd_step_options:
        fd_step_options["rad"] = 0.5
    fd_step = [fd_step_options[name] for name in names]

    # Logger Options
    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="PySLSQP",
        options={
            "acc": Solver_Settings["acc"],
            "maxiter": Solver_Settings["maxiter"],
            "finite_diff_abs_step": fd_step,
        },
        save_folder=Misc["save_folder"],
    )

    # Internal cache so RunAPDL is only executed once per unique x
    cache = {"x": None, "f": None, "c": None,}

    # Helper
    def arr(v):
        return np.asarray(v, dtype=float).ravel()

    # Model evaluation
    def evaluate_model(x):
        # Read variables and check
        x = arr(x)
        if cache["x"] is not None and np.array_equal(x, cache["x"]):
            return cache["f"], cache["c"]
        
        # Convert variables to dict again
        var_dict = dict(zip(names, x))

        # Run the Model
        f = RunAPDL(mapdl, x, Misc, opti_settings)
        logger.log_evaluation(x, f)

        # Initiate Post Processing
        pp = PostProcessor(var_dict,Misc,opti_settings)

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
        # Ensure all thickness variables remain above the geometric minimum
        thickness_constraints = [
            v - Misc["eps_geom"]
            for name, v in var_dict.items()
            if name.startswith("t0") or name.startswith("t1")
        ]
        
        # Geometric constraints: ensure inner radius remains positive for valid circular tubes
        # NEWLY ADDED
        # For columns: d0/2 - t0 > 0 (inner radius must be positive)
        # For braces: d1/2 - t1 > 0 (inner radius must be positive)
        #geometric_constraints = []
        #for name, value in var_dict.items():
        #    if name.startswith("d0_"):
        #        seg_idx = int(name.split("_")[1])
        #        t_name = f"t0_{seg_idx}"
        #        if t_name in var_dict:
        #            # d0/2 - t0 > 0
        #            geometric_constraints.append(value / 2 - var_dict[t_name])
        #    elif name.startswith("d1_"):
        #        seg_idx = int(name.split("_")[1])
        #        t_name = f"t1_{seg_idx}"
        #        if t_name in var_dict:
        #            # d1/2 - t1 > 0
        #            geometric_constraints.append(value / 2 - var_dict[t_name])
        
        # Combine all inequality constraints for optimization
        constraints = [
            *map(arr, thickness_constraints),
            #*map(arr, geometric_constraints), # NEWLY ADDED
            *[v for pair in col_brace for v in pair],
            1 - arr(pp.Util_BS()),
            *map(arr, pp.Class_2()),
            arr(pp.Eigenvalue_1())
        ]
        # Collect them together
        c = np.concatenate(constraints)

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
    os.makedirs(Misc["save_folder"], exist_ok=True) 
    save_filename = os.path.join(Misc["save_folder"], "pyslsqp_history.hdf5") # Save File
    summary_filename = os.path.join(Misc["save_folder"], "slsqp_summary.out") # Summary File

    # Run PySLSQP
    result = pyslsqp_optimize( 
        x0=x0, 
        obj=objective, 
        con=constraints, 
        meq=0,                          # all constraints are inequalities
        xl=xl,
        xu=xu, 
        finite_diff_abs_step=fd_step,                # finite difference step size for each variable
        maxiter=Solver_Settings["maxiter"], 
        acc=Solver_Settings["acc"],     # Objective Function Tolerance
        iprint=2,                       # print iteration info
        save_itr="major",               # save major iterations
        save_vars=[
            "majiter",                  # major iteration
            "x",                        # design variables
            "objective",                # objective function
            "constraints",              # constraint vector
            "optimality",               # optimality
            "feasibility",              # feasibility
            "step",                     # step
        ],
        save_filename=save_filename, 
        summary_filename=summary_filename, 
        #warm_start=True, # For restarting optimization from prior optimization runs
        #load_filename=save_filename, # Filename for restart file
        visualize=True, 
        visualize_vars=['objective', 'optimality', 'feasibility', 'x[0]', 'gradient[0]', 'constraints[0]', 'multipliers[0]', 'jacobian[0,0]'], 
    ) 
    # return the result, the text path, and the csv path
    return result, logger.txt_path, logger.csv_path 