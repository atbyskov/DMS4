import os
import numpy as np

from MyAPDLCall import RunAPDL
from opt_logger import OptimizationLogger
from Post_Process import PostProcessor
from aggregate import ConstraintAggregate
from acs import ACSclass


# Order of utilization groups in the aggregated constraint vector.
# Each group is aggregated independently (one P-norm + one ACS per group)
# so failure modes that live on different physical scales (e.g. lateral
# buckling vs. deflection) get their own scaling factor, as recommended by
# Le et al. (2010) section 6 "regional stress measures".
UTIL_GROUPS = ["LB", "NF", "S", "T", "BR", "IN", "BS_sig", "BS_defl"]


def _arr(v):
    return np.asarray(v, dtype=float).ravel()


def _collect_util_groups(pp):
    """Return ordered dict group_name -> 1-D array of raw utilization values."""
    lb_c, lb_b = pp.Util_LB()
    nf_c, nf_b = pp.Util_NF()
    s_c, s_b = pp.Util_S()
    t_c, t_b = pp.Util_T()
    br_c, br_b = pp.Util_BR()
    in_c, in_b = pp.Util_IN()
    bs_sig, bs_defl = pp.Util_BS()
    return {
        "LB":      np.concatenate([_arr(lb_c), _arr(lb_b)]),
        "NF":      np.concatenate([_arr(nf_c), _arr(nf_b)]),
        "S":       np.concatenate([_arr(s_c),  _arr(s_b)]),
        "T":       np.concatenate([_arr(t_c),  _arr(t_b)]),
        "BR":      np.concatenate([_arr(br_c), _arr(br_b)]),
        "IN":      np.concatenate([_arr(in_c), _arr(in_b)]),
        "BS_sig":  _arr(bs_sig),
        "BS_defl": _arr(bs_defl),
    }


# Optimization Function
def run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings, method="slp_mvp"):

    # Data saved in this folder:
    save_folder = Misc.get("save_folder", "Optimization_Logs")

    # Set active variables
    active = [(name, data["bounds"]) for name, data in var.items() if data.get("active", True)]
    names = [name for name, _ in active]
    Misc["active_vars"] = names

    # Print active variables
    print(f"{len(names)} Design Variables included: " +
            ", ".join(f"x{i}={n}" for i, n in enumerate(names)), flush=True)

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

    # ACS is only meaningful when an aggregation method is active.  When
    # Aggregate is None/"None", aggregate_output returns the per-element 1-u
    # vector and ACS is disabled.
    aggregate_method = Solver_Settings.get("Aggregate", "None")
    acs_active = (
        Solver_Settings.get("use_acs", True)
        and aggregate_method not in (None, "None")
    )

    # One aggregator + one ACS instance per failure-mode group, so each
    # failure mode is normalised by its own (sigma_max / sigma_PN) ratio.
    aggregators = {
        g: ConstraintAggregate(
            method=aggregate_method,
            p_value=Solver_Settings.get("p_value", 1),
            rho_value=Solver_Settings.get("rho_value", 100),
            relaxation=Solver_Settings.get("relaxation", 0),
        )
        for g in UTIL_GROUPS
    }
    acs_instances = {
        g: ACSclass(enabled=acs_active, name=g)
        for g in UTIL_GROUPS
    }

    # fminslp.m-equivalent options forwarded to solve_slp_mvp via slp_options.
    slp_options = {
        "algorithm":         Solver_Settings.get("algorithm", "merit"),
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
            "algorithm": slp_options["algorithm"],
            "acc": Solver_Settings["acc"],
            "xtol": Solver_Settings.get("xtol", 1e-8),
            "gtol": Solver_Settings.get("gtol", 1e-6),
            "finite_diff_abs_step": fd_step,
            "maxiter": Solver_Settings["maxiter"],
            "max_fun_evals": slp_options["max_fun_evals"],
            "Aggregate": aggregate_method,
            "p_value": Solver_Settings.get("p_value", 1),
            "rho_value": Solver_Settings.get("rho_value", 100),
            "relaxation": Solver_Settings.get("relaxation", 0),
            "move_limit": slp_options["move_limit"],
            "move_limit_expand": slp_options["move_limit_expand"],
            "move_limit_shrink": slp_options["move_limit_shrink"],
            "penalty_weight": Solver_Settings.get("penalty_weight", 1000.0),
            "max_infeasibility": slp_options["max_infeasibility"],
            "use_acs": acs_active,
            "util_groups": UTIL_GROUPS,
        },
        save_folder=save_folder,
    )

    # Internal cache so RunAPDL is only executed once per unique x
    cache = {
        "x": None,
        "f": None,
        "c": None,
        "util": None,
        "sigma_PN": {},   # group_name -> float (None when method is None)
        "sigma_max": {},  # group_name -> float
        "segment_masses": None,
    }

    # Model evaluation
    def evaluate_model(x):
        # Read variables and check
        x = _arr(x)
        if cache["x"] is not None and np.array_equal(x, cache["x"]):
            return cache["f"], cache["c"]

        # Convert variables to dict again
        var_dict = dict(zip(names, x))

        # Run the Model
        f, segment_masses = RunAPDL(mapdl, x, Misc, opti_settings)

        # Initiate Post Processing
        pp = PostProcessor(var_dict, Misc, opti_settings)

        # Collect raw (non-negative) utilization values per failure mode group.
        # Aggregation operates on the values directly, NOT on the violations --
        # this is what keeps sigma_max / sigma_PN well defined in the feasible
        # region (paper formulation in Le et al. 2010 Eq. 9-10).
        groups_u = _collect_util_groups(pp)

        # Aggregate each group independently
        c_rows = []
        sigma_PN_dict = {}
        sigma_max_dict = {}
        for g in UTIL_GROUPS:
            u_g = groups_u[g]
            c_scale = acs_instances[g].c if acs_active else 1.0
            c_row, sPN, sMax = aggregators[g].agg_output(u_g, c_scale=c_scale)
            c_rows.append(np.atleast_1d(c_row))
            sigma_PN_dict[g] = sPN
            sigma_max_dict[g] = sMax

        c_util_agg = np.concatenate(c_rows) if c_rows else np.empty(0)

        # Worst-group scalars (for the existing CSV columns g_max / v_agg / c_val).
        # We pick the group whose sigma_max is largest so the CSV tracks the
        # binding failure mode.
        if sigma_max_dict:
            worst_g = max(UTIL_GROUPS, key=lambda gg: (sigma_max_dict[gg] or 0.0))
        else:
            worst_g = None
        if worst_g is not None and acs_active:
            worst_sigma_PN = sigma_PN_dict[worst_g]
            worst_sigma_max = sigma_max_dict[worst_g]
            worst_c = acs_instances[worst_g].c
        else:
            worst_sigma_PN = None
            worst_sigma_max = sigma_max_dict[worst_g] if worst_g is not None else None
            worst_c = None

        # Segment mass constraints use the same c(x) >= 0 convention.
        mass_limit = float(opti_settings["segment_mass_limit"])
        mass_constraint = mass_limit - _arr(segment_masses)

        # Collect them together
        c = np.concatenate([
            c_util_agg,
            *map(_arr, pp.Class_2()),
            _arr(pp.Eigenvalue_1()),
            mass_constraint,
        ])

        print(f"Constraint vector length: {len(c)} (worst group: {worst_g})", flush=True)

        logger.log_evaluation(
            x, f, segment_masses,
            v_agg=worst_sigma_PN,
            g_max=worst_sigma_max,
            c_value=worst_c,
        )

        util_report = {
            "Util_LB":      pp.Util_LB(),
            "Util_NF":      pp.Util_NF(),
            "Util_S":       pp.Util_S(),
            "Util_T":       pp.Util_T(),
            "Util_BR":      pp.Util_BR(),
            "Util_IN":      pp.Util_IN(),
            "Util_BS_sig":  groups_u["BS_sig"],
            "Util_BS_defl": groups_u["BS_defl"],
        }

        # Update design variables and constraints
        cache.update(
            x=x.copy(),
            f=float(f),
            c=c,
            util=util_report,
            sigma_PN=sigma_PN_dict,
            sigma_max=sigma_max_dict,
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
        if acs_active:
            sigma_PN_dict = cache.get("sigma_PN", {}) or {}
            sigma_max_dict = cache.get("sigma_max", {}) or {}
            for g in UTIL_GROUPS:
                sPN = sigma_PN_dict.get(g)
                sMax = sigma_max_dict.get(g)
                if sPN is not None and sMax is not None:
                    acs_instances[g].update(sPN, sMax)
        util = cache.get("util", None)
        if util is not None:
            logger.log_utilization(util)

    # Create Folder
    os.makedirs(save_folder, exist_ok=True)

    # Run the fminslp port (Python translation of fminslp.m).  `method` is
    # accepted as a positional argument for backward compatibility but is
    # ignored -- the only available solver is the fminslp port.
    _ = method
    from .slp_MVP_advanced import solve_slp_mvp
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
