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
        """Convert input to a flat 1D float NumPy array."""
        if hasattr(v, "to_numpy"):
            return v.to_numpy(dtype=float).ravel()
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
        f_val = RunAPDL(mapdl, SWcoor, x, Misc)
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

        c_val = np.concatenate([
            # Minimum thickness constraints
            np.array([x[1] - eps_geom], dtype=float),
            np.array([x[3] - eps_geom], dtype=float),

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

        cache["x"] = x.copy()
        cache["f"] = float(f_val)
        cache["c"] = c_val

        return cache["f"], cache["c"]

    def objective(x):
        """Objective function for PySLSQP."""
        f_val, _ = evaluate_model(x)
        return f_val

    def constraints(x):
        """Constraint vector for PySLSQP."""
        _, c_val = evaluate_model(x)
        return c_val

    # PySLSQP writes its own files, so keep them in your log folder
    os.makedirs(save_folder, exist_ok=True)
    save_filename = os.path.join(save_folder, "pyslsqp_history.hdf5")
    summary_filename = os.path.join(save_folder, "slsqp_summary.out")

    # ------------------------------------------------------------------
    # Run PySLSQP
    # ------------------------------------------------------------------
    result = pyslsqp_optimize(
        x0=x0,
        obj=objective,
        con=constraints,
        meq=0,                         # all constraints are inequalities
        xl=xl,
        xu=xu,
        finite_diff_abs_step=finite_diff_abs_step,
        maxiter=maxiter,
        acc=acc,
        iprint=2,                     # print iteration info
        save_itr="major",             # save major iterations
        save_vars=[
            "majiter",
            "x",
            "objective",
            "constraints",
            "optimality",
            "feasibility",
            "step",
        ],
        save_filename=save_filename,
        summary_filename=summary_filename,
        visualize=True,
    )

    #logger.finalize(result)
    return result, logger.txt_path, logger.csv_path