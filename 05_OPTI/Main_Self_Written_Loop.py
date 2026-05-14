## MAIN BEAM ELEMENT SLP LOOPED DOCUMENT ##

# Import packages
import sys
import time
from ansys.mapdl.core import launch_mapdl

# For mail
import smtplib
from email.message import EmailMessage

# Import Functions
import Self_Written_Optimization.optimization_self_written as optimization_self_written
from MyAPDLCall import RunAPDL

start_time = time.time()


def read_initial_guesses(filename="search_space.txt"):
    """
    Reads initial guesses from search_space.txt

    Returns
    -------
    names : list[str]
        Variable names (e.g. ["rad", "d0", "t0", "d1", "t1"])
    points : list[list[float]]
        Initial guesses
    """
    names = []
    points = []

    with open(filename, "r") as f:
        lines = f.readlines()

    start = False
    for line in lines:
        line = line.strip()

        if line.startswith("Points:"):
            start = True
            continue

        if not start:
            continue

        # header line with variable names
        if not names and line:
            names = line.split()
            continue

        # data rows
        if line:
            points.append([float(v) for v in line.split()])

    return names, points


def write_search_result(initial, final, names, obj_start, obj_final,
                        filename="search_results_self_written.txt"):
    # Ensure scalar objective values
    if isinstance(obj_start, (list, tuple)):
        obj_start = obj_start[0]
    if hasattr(obj_start, "item"):
        obj_start = obj_start.item()

    if isinstance(obj_final, (list, tuple)):
        obj_final = obj_final[0]
    if hasattr(obj_final, "item"):
        obj_final = obj_final.item()

    with open(filename, "a") as f:
        f.write("\n===============================\n")

        f.write("Initial Guess:\n")
        for n, v in zip(names, initial):
            f.write(f"  {n}: {v:.6f}\n")
        f.write(f"  Objective: {obj_start:.3f}\n")

        f.write("Final Result:\n")
        for n, v in zip(names, final):
            f.write(f"  {n}: {v:.6f}\n")
        f.write(f"  Objective: {obj_final:.3f}\n")


var_names, initial_points = read_initial_guesses("search_space.txt")


# Problem shape -- must be compatible with the variables present in search_space.txt.
# search_space.txt header is: rad d0 t0 d1 t1  -> single-size columns/braces, no brace split.
opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "segment_mass_limit": 23,     # Limits for segment masses [kg]
    "multi_size_columns": False,  # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": False,   # Whether mast segments braces uses different dimensions (True) or not (False)
    "brace_split": False,         # Whether braces are split between horiontal and cross (True) or not (False)
}

# Static variables (matches Main_Self_Written_Test.py)
Misc = {
    "esize": 3,                          # Element Size [mm]
    "Hor_Force": 502.52,                 # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,              # Vertical Force (P_Load_y)   [N]
    "f_y": 700,                          # Column Yield Strength [MPa]
    "f_y_brace": 355,                    # Brace Yield Strength  [MPa]
    "E_mod": 200 * 1E3,                  # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,                # Vertical Force COG (P_COG_y) [N]
    "eps_geom": 0.1,                     # Minimum thickness specification for geometry updates [mm]
    "SW_filename": "LWC_L1_LINES.IGS",   # Filename for IGS File
    "save_folder": "Optimization_Logs"   # Save Folder
}

# Solver Settings -- SLP MVP advanced (matches Main_Self_Written_Test.py)
Solver_Settings = {
    "acc": 1e-3,                  # Maximum objective function tolerance
    "maxiter": 150,               # Maximum iterations                                                         # IMPORTANT TO TUNE FOR FEA
    "Aggregate": None,            # None, "P-norm", "P-norm-mean", "KS", "KS_shift"  (Write exacly)
    "p_value": 8,                 # Value for "P-norm" and "P-norm-mean"
    "rho_value": 100,             # rho value used in KS
    "relaxation": 0,              # Relaxation parameter used in aggregation
    "finite_diff_rel_step": 5e-3, # PySLSQP-like absolute FD step: rel_step * max(1, abs(x))                   # IMPORTANT TO TUNE FOR FEA
    "algorithm": "merit",         # 'merit' (quadratic merit) or 'al' (augmented Lagrangian)
    "penalty_weight": 1,          # Initial L1 merit/slack penalty                                               # 1
    "penalty_increase": 2.0,      # Increase penalty when slacks/infeasibility persist
    "penalty_max": 1e6,           # Upper cap for adaptive penalty
    "move_limit": 0.10,           # Initial/max scaled trust-region size
    "move_limit_min": 1e-3,       # Stop after failed globalization below this scaled size
    "move_limit_expand": 1.1,     # Expand trust region after good predicted-vs-actual agreement
    "move_limit_shrink": 0.5,     # Shrink trust region after rejection or poor agreement
    "backtrack_max": 8,           # Trial alphas along each LP step
    "backtrack_shrink": 0.5,      # alpha <- alpha * shrink during SLP line search
    "feasibility_tol": 1e-5,      # Required max nonlinear constraint violation for convergence
    "slack_tol": 1e-8,            # Required LP slack size for convergence
    "xtol": 1e-4,                 # Physical design step tolerance
    "gtol": 1e-6,                 # Scaled gradient infinity-norm tolerance
    "sufficient_decrease": 1e-4,  # Merit decrease fraction of predicted decrease
    "feasibility_reduction": 1e-3,# Accept infeasible steps that reduce violation by this fraction            # IMPORTANT TO TUNE FOR FEA
    "consecutive_reject_limit": 10, # Abort this run after this many rejected outer iters in a row (0 = off)
}

# Optimization bounds (the bounds passed to the solver, may be wider than search-sampling bounds)
bounds = {
    "rad": (150.0, 350.0),   # Radius Structure [mm]
    "d0":  (48.3, 108.0),    # Column Diameter  [mm]
    "t0":  (2.5, 5.0),       # Column Thickness [mm]
    "d1":  (25.0, 60.0),    # Brace  Diameter  [mm]
    "t1":  (2.0, 6.0),       # Brace  Thickness [mm]
}

# Launch MAPDL once (re-used across every starting guess).
# can vary across calls and introduces ~1e-3..1e-1 noise in mass/stress, which kills
# the SLP's finite-difference gradients and the filter's accept/reject test.
mapdl = launch_mapdl(
    run_location="Ansout",
    log_apdl="apdl_log",
    override=True,
    nproc=8,
    additional_switches="-p ansys -smp",
)

try:
    for i, x0 in enumerate(initial_points, start=1):
        print(f"\n========== Run {i} / {len(initial_points)} ==========", flush=True)

        # Build var dict for this starting guess
        var = {}
        for name, value in zip(var_names, x0):
            if name not in bounds:
                raise KeyError(
                    f"Variable {name!r} from search_space.txt has no bounds entry. "
                    f"Known: {sorted(bounds)}"
                )
            var[name] = {
                "value": value,
                "bounds": bounds[name],
                "active": True,
            }

        # Initial objective for logging (one FEA call at x0)
        f1 = RunAPDL(mapdl, var, Misc, opti_settings)

        # Run SLP optimisation
        result, txt_path, csv_path = optimization_self_written.run_optimization(
            mapdl, opti_settings, var, Misc, Solver_Settings, method="slp_mvp"
        )

        # Extract final design and objective (SLP returns a Result dataclass)
        final_x = result.x
        obj_final = result.fun
        if isinstance(obj_final, (list, tuple)):
            obj_final = obj_final[0]

        obj_start = f1

        # Per-run summary -- highlight the consecutive-rejection abort path
        if "consecutive rejected" in (result.message or ""):
            print(
                f"[Run {i}] FAILED -- aborted by consecutive-rejection guard. "
                f"{result.message} Moving on to next start guess.",
                flush=True,
            )
        else:
            print(
                f"[Run {i}] finished. success={result.success}  nit={result.nit}  "
                f"nfev={result.nfev}  f={obj_final:.3f}  msg: {result.message}",
                flush=True,
            )

        write_search_result(x0, final_x, var_names, obj_start, obj_final)
finally:
    mapdl.exit()

end_time = time.time()


def send_email(subject, body):
    msg = EmailMessage()
    msg.set_content(body)
    msg["Subject"] = subject
    msg["From"] = "atbyskov@gmail.com"
    msg["To"] = "atbyskov@gmail.com"

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login("atbyskov@gmail.com", "zani dtwc hnmw dxpm")
        server.send_message(msg)


tot_time = end_time - start_time

send_email(
    "Python SLP loop job finished",
    "Self-written SLP optimization loop script has completed\n"
    f"Total Time: {tot_time / 3600} Hours"
)
