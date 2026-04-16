import os
import numpy as np
from pyslsqp import optimize as pyslsqp_optimize

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
    """
    Run constrained optimization with PySLSQP.

    Design variables:
        x[0] = Column outer diameter [mm]
        x[1] = Column thickness      [mm]
        x[2] = Brace outer diameter  [mm]
        x[3] = Brace thickness       [mm]

    Constraint convention for PySLSQP:
        - First meq constraints are equalities: c_i(x) = 0
        - Remaining constraints are inequalities: c_i(x) >= 0

    This implementation assumes:
        - all constraints are inequalities
        - utilizations must satisfy utilization <= 1
        - Class_2 and Eigenvalue_1 already return values in >= 0 form
          if they are to be used directly as inequality constraints
    """

    x0 = np.asarray(var, dtype=float)

    bounds = [
        (40.0, 100.0),  # Column outer diameter [mm]
        (1.0,   7.0),   # Column thickness [mm]
        (10.0, 100.0),  # Brace outer diameter [mm]
        (0.1,   7.0),   # Brace thickness [mm]
    ]

    xl = np.array([b[0] for b in bounds], dtype=float)
    xu = np.array([b[1] for b in bounds], dtype=float)

    # Settings that roughly correspond to your old SciPy options
    finite_diff_abs_step = 0.01
    acc = 1e-3
    maxiter = 40

    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="PySLSQP",
        options={
            "finite_diff_abs_step": finite_diff_abs_step,
            "acc": acc,
            "maxiter": maxiter,
        },
        save_folder=save_folder,
    )

    # ------------------------------------------------------------------
    # Internal cache so RunAPDL is only executed once per unique x
    # ------------------------------------------------------------------
    cache = {
        "x": None,
        "f": None,
        "c": None,
    }

    def as_1d_float_array(v):
        #This function is doing defensive programming:
        #"No matter what weird format I get, I want a clean np.ndarray of shape (n,) with floats."
        """Convert input to a flat 1D float NumPy array."""
        if hasattr(v, "to_numpy"): #CHeck if v has a .to_numpy method, and this is typical for pandas.series and pandas.dataframe
            return v.to_numpy(dtype=float).ravel() # convert to numpy array of floats and flatten it to 1D
            # example:
            #v = pd.Series([1, 2, 3])
            #v.to_numpy(dtype=float).ravel() -> array([1., 2., 3.])
        return np.asarray(v, dtype=float).ravel()

    def evaluate_model(x):
        """
        Shared evaluation of objective + constraints at one x.
        Uses caching so repeated calls at the same x do not rerun APDL.
        """
        x = np.asarray(x, dtype=float).ravel()

        if cache["x"] is not None and np.array_equal(x, cache["x"]):
            return cache["f"], cache["c"]

        # Run the expensive model once
        f_val = RunAPDL(mapdl, x, Misc)
        logger.log_evaluation(x, f_val)

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

        c_val = np.concatenate([ #concatenate the arrays into a single array, and return a numpy array of shape (n,) with floats
            # Minimum thickness constraints
            np.array([x[1] - eps_geom], dtype=float), # create a numpy array of shape (1,) with floats
            np.array([x[3] - eps_geom], dtype=float), # create a numpy array of shape (1,) with floats

            # Utilization constraints: require utilization <= 1  ->  1 - util >= 0
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

            # These are assumed already written in c(x) >= 0 form
            as_1d_float_array(Util_Class_2_values_col), 
            as_1d_float_array(Util_Class_2_values_brace),
            as_1d_float_array(Eigenvalue_1_values),
        ]).astype(float, copy=False)

        print(f"Constraint vector length: {len(c_val)}")

        cache["x"] = x.copy() # copy the design variables to the cache
        cache["f"] = float(f_val) # convert the objective function to a float
        cache["c"] = c_val # assign the constraint vector to the cache

        return cache["f"], cache["c"] # return the objective function and the constraint vector

    def objective(x): #This is the objective function for PySLSQP.
        """Objective function for PySLSQP."""
        f_val, _ = evaluate_model(x) # evaluate the model and return the objective function and the constraint vector
        return f_val # return the objective function

    def constraints(x):
        """Constraint vector for PySLSQP."""
        _, c_val = evaluate_model(x) # evaluate the model and return the objective function and the constraint vector
        return c_val # return the constraint vector

    # PySLSQP writes its own files, so keep them in your log folder
    os.makedirs(save_folder, exist_ok=True) # create the save folder if it doesn't exist
    save_filename = os.path.join(save_folder, "pyslsqp_history.hdf5") # create the save filename
    summary_filename = os.path.join(save_folder, "slsqp_summary.out") # create the summary filename

    # ------------------------------------------------------------------
    # Run PySLSQP
    # ------------------------------------------------------------------
    result = pyslsqp_optimize( #This is the optimization function for PySLSQP.
        x0=x0, # initial guess for the design variables
        obj=objective, # objective function
        con=constraints, # constraint function
        meq=0,                         # all constraints are inequalities
        xl=xl, # lower bound for the design variables
        xu=xu, # upper bound for the design variables
        finite_diff_abs_step=finite_diff_abs_step, # finite difference absolute step
        maxiter=maxiter, # maximum number of iterations
        acc=acc, # accuracy Equal to the tolerance for the scipy optimization, tolerance in change in objective function value
        iprint=2, # print iteration info
        save_itr="major", # save major iterations
        save_vars=[
            "majiter", # major iteration
            "x", # design variables
            "objective", # objective function
            "constraints", # constraint vector
            "optimality", # optimality
            "feasibility", # feasibility
            "step", # step
        ],
        save_filename=save_filename, # save the history of the optimization
        summary_filename=summary_filename, # save the summary of the optimization
        visualize=True, # visualize the optimization
        visualize_vars=['objective', 'optimality', 'feasibility', 'x[0]', 'gradient[0]', 'constraints[0]', 'multipliers[0]', 'jacobian[0,0]'], # visualize the optimization variables
    ) 
    #logger.finalize(result) # finalize the optimization
    return result, logger.txt_path, logger.csv_path # return the result, the text path, and the csv path