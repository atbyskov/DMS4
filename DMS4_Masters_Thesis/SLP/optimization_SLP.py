"""
Optimization_SLP.py
-------------------

This module implements the core optimization workflow for the structural mast
model using the self-written Sequential Linear Programming (SLP) solver. It
serves as the interface between the SLP engine, the finite element solver
(ANSYS MAPDL), and the post-processing routines.

Dependencies:
-------------
- SLP (self-written optimization engine, see SLP/SLP.py)
- ansys.mapdl.core (external, used indirectly via RunAPDL)
- Custom modules:
    * MyAPDLCall   -> FE model execution
    * Post_Process -> Structural response evaluation
    * aggregate    -> Constraint aggregation methods
    * acs          -> Adaptive constraint scaling
    * opt_logger   -> Logging and data storage

Inputs:
-------
- mapdl            : Active MAPDL session
- opti_settings    : Optimization configuration settings
- var              : Dictionary of design variables
- Misc             : Static parameters (loads, material, etc.)
- Solver_Settings  : Optimization and solver parameters

Outputs:
--------
- result           : SLP optimization result object
- txt_path         : Path to the detailed log file
- csv_path         : Path to the evaluation history file
"""

import os
import numpy as np

from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor
from aggregate import ConstraintAggregate
from acs import ACSclass

# Optimization Function
def run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings, method="slp_mvp"):
    
    # Data saved in this folder:
    save_folder = Misc.get("save_folder", "Optimization_Logs")

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

    # Use the same absolute finite-difference convention as the PySLSQP path.
    # This keeps both optimizers comparing the same local model quality.
    fd_step = [
        Solver_Settings.get("finite_diff_rel_step", 1e-3) * max(1.0, abs(var[name]["value"]))
        for name in names
    ]

    # ACS adaptive constraint scaling is only meaningful when an aggregation
    # method is active. If Aggregate is None/"None", agg_output returns the
    # raw constraint vector and ACS has nothing to scale.
    aggregate_method = Solver_Settings.get("Aggregate", "None")
    acs_active = aggregate_method not in (None, "None")
    acs = ACSclass(
        alpha=Solver_Settings.get("acs_alpha", 1.0),
        c0=1.0,
    ) if acs_active else None

    # fminslp.m-equivalent options forwarded to solve_slp_mvp via slp_options.
    # Names match fminslp.slpoptions (lines 794-849).  Extra keys not in
    # fminslp.m are intentionally absent -- the solver accepts and ignores
    # legacy keys for backward compatibility.
    slp_options = {
        "move_limit":        Solver_Settings.get("move_limit", 0.10),
        "move_limit_expand": Solver_Settings.get("move_limit_expand", 1.1),
        "move_limit_shrink": Solver_Settings.get("move_limit_shrink", 0.5),
        "max_fun_evals":     Solver_Settings.get("max_fun_evals", 1000),
        "objective_limit":   Solver_Settings.get("objective_limit", -1e20),
        "max_infeasibility": Solver_Settings.get("max_infeasibility", float("inf")),
    }

    # Logger Options
    logger = OptimizationLogger(
        x0=x0,
        bounds=bounds,
        method="Self-Written-" + method,
        var_names=names,
        options={
            "acc": Solver_Settings["acc"],
            "xtol": Solver_Settings.get("xtol", 1e-8),
            "gtol": Solver_Settings.get("gtol", 1e-6),
            "finite_diff_abs_step": fd_step,
            "maxiter": Solver_Settings["maxiter"],
            "max_fun_evals": slp_options["max_fun_evals"],
            "Aggregate": Solver_Settings.get("Aggregate", "None"),
            "p_value": Solver_Settings.get("p_value", 1),
            "rho_value": Solver_Settings.get("rho_value", 100),
            "relaxation": Solver_Settings.get("relaxation", 0),
            "move_limit": slp_options["move_limit"],
            "move_limit_expand": slp_options["move_limit_expand"],
            "move_limit_shrink": slp_options["move_limit_shrink"],
            "penalty_weight": Solver_Settings.get("penalty_weight", 1000.0),
            "max_infeasibility": slp_options["max_infeasibility"],
        },
        save_folder=save_folder,
    )

    # Internal cache so RunAPDL is only executed once per unique x
    cache = {
        "x": None,
        "f": None,
        "c": None,
        "util": None,
        "v_agg": None,
        "g_max": None,
        "segment_masses": None,
    }

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
        f, segment_masses = RunAPDL(mapdl, x, Misc, opti_settings)

        # Initiate Post Processing
        pp = PostProcessor(var_dict, Misc, opti_settings)

        # Call Utlization Ratios
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
                1 - arr(pp.Util_BS()[0]), # stress constraint here
                1 - arr(pp.Util_BS()[1])]) # deflection constraint here
        
        # Handle aggregation of constraints
        agg = ConstraintAggregate(
            method = Solver_Settings.get("Aggregate", "None"),
            p_value = Solver_Settings.get("p_value", 1),
            rho_value = Solver_Settings.get("rho_value", 100),
            relaxation = Solver_Settings.get("relaxation", 0)
        )

        if acs_active:
            c_util_agg, v_agg, g_max = agg.agg_output(c_util, c_scale=acs.c)
        else:
            c_util_agg, v_agg, g_max = agg.agg_output(c_util)
        
        # Segment mass constraints use the same c(x) >= 0 convention.
        mass_limit = float(opti_settings["segment_mass_limit"])
        mass_constraint = mass_limit - arr(segment_masses)

        # Collect them together
        c = np.concatenate([
            np.atleast_1d(c_util_agg),
            *map(arr, pp.Class_2()),
            arr(pp.Eigenvalue_1()),
            mass_constraint,
        ])

        print(f"Constraint vector length: {len(c)}")

        logger.log_evaluation(
            x, f, segment_masses,
            v_agg=v_agg,
            g_max=g_max,
            c_value=(acs.c if acs_active else None),
        )
        
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
        cache.update(
            x=x.copy(),
            f=float(f),
            c=c,
            util=util_report,
            v_agg=v_agg,
            g_max=g_max,
            segment_masses=segment_masses,
        )

        return f, c
    
    # Objective function call that only returns mass
    def objective(x):
        """Objective function."""
        f_val, _ = evaluate_model(x) 
        return f_val 

    # Constrain function call that only returns constraints
    def constraints(x):
        """Constraint vector."""
        _, c_val = evaluate_model(x) 
        return c_val 

    def iteration_callback(x):
        logger.log_iteration(x)
        if acs_active and cache["v_agg"] is not None and cache["g_max"] is not None:
            acs.update(cache["v_agg"], cache["g_max"])
        util = cache.get("util", None)
        if util is not None:
            logger.log_utilization(util)

    # Create Folder
    os.makedirs(save_folder, exist_ok=True) 

    # Run the fminslp port (Python translation of fminslp.m).  `method` is
    # accepted as a positional argument for backward compatibility but is
    # ignored -- the only available solver is the fminslp port.
    _ = method
    from .SLP import solve_slp_mvp
    result = solve_slp_mvp(
        obj=objective,
        con=constraints,
        x0=x0,
        xl=xl,
        xu=xu,
        fd_step=fd_step,
        fd_type="forward",
        maxiter=Solver_Settings["maxiter"],
        ftol=Solver_Settings["acc"],
        xtol=Solver_Settings.get("xtol", 1e-8),
        gtol=Solver_Settings.get("gtol", 1e-6),
        infeasibility_penalty=Solver_Settings.get("penalty_weight", 1000.0),
        callback=iteration_callback,
        display=True,
        **slp_options,
    )

    # Write a clean FINAL RESULT block (converged x, objective, message) to the TXT log.
    logger.finalize(result)

    # return the result, the text path, and the csv path
    return result, logger.txt_path, logger.csv_path
