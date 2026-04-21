# Import packages
import os
import numpy as np
from pyslsqp import optimize as pyslsqp_optimize

# Import Functions
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor
from aggregate import ConstraintAggregate

# Optimization Function
def run_optimization(mapdl,var,Misc,Solver_Settings):

    # Set active variables
    active = [(name, data["bounds"]) for name, data in var.items() if data.get("active", True)]
    names = [name for name,_ in active]
    Misc["active_vars"] = names

    # Print active variables
    print(f"{len(names)} Design Variables included: " +
            ", ".join(f"x{i}={n}" for i,n in enumerate(names)))

    # Set initial guess, bounds and ranges to np array
    x0 = np.array([var[name]["value"] for name, _ in active], dtype=float)
    bounds = [bnds for _, bnds in active]
    xl, xu = np.array(bounds).T

    # Step Options for each variable
    fd_step_options = {
        "d0": 0.01,
        "t0": 0.01,
        "d1": 0.01,
        "t1": 0.01,
        "rad": 2
        }
    fd_step = [fd_step_options[name] for name, _ in active]

    # Logger Options
    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="PySLSQP",
        options={
            "finite_diff_abs_step": fd_step,
            "acc": Solver_Settings["acc"],
            "maxiter": Solver_Settings["maxiter"],
        },
        save_folder=Misc["save_folder"],
    )

    # Internal cache so RunAPDL is only executed once per unique x
    cache = {"x": None, "f": None, "c": None, "util": None}

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
        f = RunAPDL(mapdl, x, Misc)
        logger.log_evaluation(x, f)

        # Initiate Post Processing
        pp = PostProcessor(var_dict,Misc)

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
        col_brace = [(arr(c), arr(b)) for c,b in utils]

        # Set up the constraints
        # Can i only aggregate the utils or do we need to have them all
        c_geom = np.concatenate([
                 arr(x[1] - Misc["eps_geom"]),
                 arr(x[3] - Misc["eps_geom"])])
        
        c_util = np.concatenate([
                *[v for pair in col_brace for v in pair],
                arr(pp.Util_BS())])
        
        c_misc = np.concatenate([
                *map(arr, pp.Class_2()),
                arr(pp.Eigenvalue_1())])
        
        # Handle aggregation of constraints
        agg = ConstraintAggregate(
            method = Solver_Settings["Aggregate"],
            p_value = Solver_Settings["p_value"],
            rho_value = Solver_Settings["rho_value"]
        )

        c_util_agg = agg.agg_output(c_util)
        
        # Collect them together
        c = np.concatenate([c_geom,np.atleast_1d(c_util_agg),c_misc])
        # Print length of constraints
        print(f"Constraint vector length: {len(c)}")

        
        
        util_report = {
                "Util_LB": pp.Util_LB(),
                "Util_NF": pp.Util_NF(),
                "Util_S":  pp.Util_S(),
                "Util_T":  pp.Util_T(),
                "Util_BNS": pp.Util_BNS(),
                "Util_BR": pp.Util_BR(),
                "Util_IN": pp.Util_IN(),
                "Util_BS": pp.Util_BS(),
            }

        # Update design variables and constraints
        cache.update(x=x.copy(), f=float(f), c=c, util = util_report)

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

    
    def iteration_callback(x):
        logger.log_iteration(x)

        util = cache.get("util", None)
        if util is not None:
            logger.log_utilization(util)


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
        callback=iteration_callback,
        finite_diff_abs_step=fd_step, 
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
        visualize=False, 
        visualize_vars=['objective', 'optimality', 'feasibility', 'x[0]', 'gradient[0]', 'constraints[0]', 'multipliers[0]', 'jacobian[0,0]'], 
    ) 
    # return the result, the text path, and the csv path
    return result, logger.txt_path, logger.csv_path 