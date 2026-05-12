## MAIN BEAM ELEMENT DOCUMENT ##

# Import packages
import sys
from ansys.mapdl.core import launch_mapdl
import time

# For mail
import smtplib
from email.message import EmailMessage

# Import Functions
import optimization
import Self_Written_Optimization.optimization_self_written as optimization_self_written
from MyAPDLCall import RunAPDL

start_time = time.time()
def read_initial_guesses(filename="search_space.txt"):
    """
    Reads initial guesses from search_space.txt

    Returns
    -------
    names : list[str]
        Variable names (e.g. ["d0", "t0", "d1", "t1", "rad"])
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

def write_search_result(initial, final, names, obj_start, obj_final, filename="search_results.txt"):
    # Ensure scalar objective values
    if isinstance(obj_start, (list, tuple)):
        obj_start = obj_start[0]
    if hasattr(obj_start, "item"):  # handles numpy scalars
        obj_start = obj_start.item()

    if isinstance(obj_final, (list, tuple)):
        obj_final = obj_final[0]
    if hasattr(obj_final, "item"):  # handles numpy scalars
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



# Solver Settings
Solver_Settings = {
    "acc": 1e-3,                 # Maximum objective function tolerance
    "maxiter": 60,               # Maximum iterations
    "Aggregate": "P-norm",              # None, "P-norm", "P-norm-mean", "KS", "KS_shift"  (Write exacly)
    "p_value": 8,                # Value for "P-norm" and "P-norm-mean"
    "rho_value": 100,            # rho value used in KS
    "relaxation": 0,             # Relaxation parameter used in aggregation
}
# Static variables
Misc = {
    "esize": 3,                     # Element Size [mm]
    "Hor_Force": 502.52,            # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)   [N]
    "f_y": 700 ,                    # Column Yield Strength [MPa]
    "f_y_brace": 355,               # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
    "eps_geom": 0.1,    # Minimum thickness specification for geometry updates [mm]
    "SW_filename": "LWC_L1_LINES.IGS",    # Filename for IGS File
    "save_folder": "Optimization_Logs" # Save Folder
}
# Launch MAPDL
mapdl = launch_mapdl(
    run_location="Ansout",
    log_apdl="apdl_log",
    override=True,
    nproc=6,
    additional_switches="-p ansys -smp",
)

bounds = {
    "d0": (40.0, 100.0),
    "t0": (1.0, 7.0),
    "d1": (10.0, 50.0),
    "t1": (0.1, 4.0),
    "rad": (150.0, 350.0),
}

opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "segment_mass_limit": 23,  # Limits for segment masses [kg]
    "multi_size_columns": False,   # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": False,    # Whether mast segments braces uses different dimensions (True) or not (False)
    "brace_split": False,          # Whether braces are split between horiontal and cross (True) or not (False)
}

for i, x0 in enumerate(initial_points, start=1):
    # Variables. Choose which ones to include by setting "active": True or False.
    # Include Bounds
    var = {}
    for name, value in zip(var_names, x0):
        var[name] = {
            "value": value,
            "bounds": bounds[name],
            "active": True
        }

    # Ensures that MAPDL closes if something chrashes

    f1 = RunAPDL(mapdl,var,Misc,opti_settings)
    result, txt_path, csv_path = optimization.run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings)


    # Append results
    final_x = result["x"]
    print(result)
    obj_final = result["objective"]
    if isinstance(obj_final, (list, tuple)):
        obj_final = obj_final[0]

    obj_start = f1

    write_search_result(x0, final_x, var_names, obj_start, obj_final)
        

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

tot_time = end_time-start_time

send_email(
    "Python job finished ✅",
    "optimization script has completed \n" \
    f"Total Time: {end_time/3600} Hours"
)