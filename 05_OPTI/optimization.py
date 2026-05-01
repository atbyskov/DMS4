# Import packages
import os
import numpy as np
from pyslsqp import optimize as pyslsqp_optimize

# Import Functions
from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor
from aggregate import ConstraintAggregate
from acs import ACSclass

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
    fd_step = []

    for name in names:
        xval = var[name]["value"]
        step = 1E-4 * xval
        fd_step.append(step)

    acs = ACSclass(
        alpha=Solver_Settings.get("acs_alpha",1.0),
        c0 = 1.0
    )

    # Logger Options
    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="PySLSQP",
        var_names=names,
        options={
            "acc": Solver_Settings["acc"],
            "finite_diff_abs_step": fd_step, # fd_step
            "maxiter": Solver_Settings["maxiter"],
            "Aggregate": Solver_Settings["Aggregate"],
            "p_value": Solver_Settings["p_value"],
            "rho_value": Solver_Settings["rho_value"],
            "relaxation": Solver_Settings["relaxation"],
        },
        save_folder=Misc["save_folder"],
    )

    # Internal cache so RunAPDL is only executed once per unique x
    cache = {"x": None, "f": None, "c": None, "util": None, "v_agg": None, "g_max": None}

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
        f = RunAPDL(mapdl, x, Misc,opti_settings)

        # Initiate Post Processing
        pp = PostProcessor(var_dict,Misc,opti_settings)

        # Call Utlization Ratios
        #### FIX DOUBLE utils CALL !
        utils = [
            pp.Util_LB(),
            pp.Util_NF(),
            pp.Util_S(),
            pp.Util_T(),
            #pp.Util_BNS(),
            pp.Util_BR(),
            pp.Util_IN(),
        ]
        col_brace = [(1-arr(c), 1-arr(b)) for c,b in utils]

        # Set up the constraints      
        c_util = np.concatenate([
                *[v for pair in col_brace for v in pair],
                1- arr(pp.Util_BS()[0]), # stress constraint here
                1- arr(pp.Util_BS()[1])]) # deflection constraint here
        
        # Handle aggregation of constraints
        agg = ConstraintAggregate(
            method = Solver_Settings["Aggregate"],
            p_value = Solver_Settings["p_value"],
            rho_value = Solver_Settings["rho_value"],
            relaxation = Solver_Settings["relaxation"],
        )

        c_util_agg, v_agg, g_max = agg.agg_output(c_util, c_scale=acs.c)
        
        # Collect them together
        c = np.concatenate([
            np.atleast_1d(c_util_agg),
            *map(arr, pp.Class_2()),
            arr(pp.Eigenvalue_1())
        ])

        print(f"Constraint vector length: {len(c)}")

        logger.log_evaluation(x, f,v_agg=v_agg,g_max=g_max)

        
        
        util_report = {
                "Util_LB": pp.Util_LB(),
                "Util_NF": pp.Util_NF(),
                "Util_S":  pp.Util_S(),
                "Util_T":  pp.Util_T(),
                #"Util_BNS": pp.Util_BNS(),
                "Util_BR": pp.Util_BR(),
                "Util_IN": pp.Util_IN(),
                "Util_BS_sig": pp.Util_BS()[0], # stress
                "Util_BS_defl": pp.Util_BS()[1], # deflection
            }

        # Update design variables and constraints
        cache.update(x=x.copy(), f=float(f), c=c, util = util_report, v_agg = v_agg, g_max=g_max)

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

        if cache["v_agg"] is not None and cache["g_max"] is not None:
            acs.update(cache["v_agg"], cache["g_max"])

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
        finite_diff_rel_step=1E-4, 
        maxiter=Solver_Settings["maxiter"], 
        acc=Solver_Settings["acc"],     # Objective Function Tolerance
        iprint=2,                       # print iteration info
        #obj_scaler=0.8,                # Scale objective function for stability
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
        visualize=False, 
        visualize_vars=['objective', 'optimality', 'feasibility', 'x[0]'], # Extra: 'gradient[0]', 'constraints[0]', 'multipliers[0]', 'jacobian[0,0]'
    ) 
    # return the result, the text path, and the csv path
    return result, logger.txt_path, logger.csv_path 