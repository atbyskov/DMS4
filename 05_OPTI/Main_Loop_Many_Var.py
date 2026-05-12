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

def write_search_result(initial, final, names, filename="search_results.txt"):
    with open(filename, "a") as f:
        f.write("\n===============================\n")
        f.write("Initial Guess:\n")
        for n, v in zip(names, initial):
            f.write(f"  {n}: {v:.6f}\n")

        f.write("Final Result:\n")
        for n, v in zip(names, final):
            f.write(f"  {n}: {v:.6f}\n")

var_names, initial_points = read_initial_guesses("search_space.txt")

# Variables. Choose which ones to include by setting "active": True or False.
opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "segment_mass_limit": 23,  # Limits for segment masses [kg]
    "multi_size_columns": True,   # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": True,    # Whether mast segments braces uses different dimensions (True) or not (False)
    "brace_split": True,          # Whether braces are split between horiontal and cross (True) or not (False)
}

# Solver Settings
Solver_Settings = {
    "acc": 1e-3,                 # Maximum objective function tolerance
    "maxiter": 100,               # Maximum iterations
    "Aggregate": None,              # None, "P-norm", "P-norm-mean", "KS", "KS_shift"  (Write exacly)
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
    nproc=10,
    additional_switches="-p ansys -smp",
)

# Bounds
column_diameter_bounds = (40.0, 100)   # Column Diameter Bounds [mm]
column_thickness_bounds = (1.0, 7.0)   # Column Thickness Bounds [mm]
brace_diameter_bounds = (10.0, 100.0)   # Brace Diameter Bounds [mm]
brace_thickness_bounds = (1.0, 7.0)    # Brace Thickness Bounds [mm]

# Initial Guess
column_diameter = 76.1 # Column Diameter [mm]
column_thickness = 3.0 # Column Thickness [mm]
brace_diameter = 26.9  # Brace Diameter [mm]
brace_thickness = 2.3  # Brace Thickness [mm]

# Defining variables with bounds and active status
var = {
    "rad": {"value": 202.07, "bounds": (150.0, 350.0), "active": True}, # Radius Structure [mm]
}
if opti_settings["multi_size_columns"]:
    var.update({
        **{f"d0_{i}": {"value": column_diameter, "bounds": column_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},           # Column Diameter  [mm]
        **{f"t0_{i}": {"value": column_thickness,  "bounds": column_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Thickness [mm]
    })
else:
    var.update({
        "d0": {"value": column_diameter, "bounds": column_diameter_bounds, "active": True},           # Column Diameter  [mm]
        "t0": {"value": column_thickness,  "bounds": column_thickness_bounds,  "active": True},       # Column Thickness [mm]
    })
if opti_settings["multi_size_braces"]:
    if opti_settings["brace_split"]:
        var.update({
            **{f"d1_h_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},          # Horizontal Brace Diameter   [mm]
            **{f"t1_h_{i}": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Horizontal Brace Thickness  [mm]
            **{f"d1_c_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},           # Cross Brace Diameter       [mm]
            **{f"t1_c_{i}": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Cross Brace Thickness      [mm]
        })
    else:
        var.update({
            **{f"d1_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},           # Brace Diameter   [mm]
            **{f"t1_{i}": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Brace Thickness  [mm]
        })
else:
    if opti_settings["brace_split"]:
        var.update({
            "d1_h": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},           # Horizontal Brace Diameter   [mm]
            "t1_h": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True},       # Horizontal Brace Thickness  [mm]
            "d1_c": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},           # Cross Brace Diameter       [mm]
            "t1_c": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True},       # Cross Brace Thickness      [mm]
        })
    else:
        var.update({
            "d1": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},           # Brace Diameter   [mm]
            "t1": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True},       # Brace Thickness  [mm]
        })

for i, x0 in enumerate(initial_points, start=1):

    # Map initial guesses to base names
    x0_dict = dict(zip(var_names, x0))

    var = {
        "rad": {
            "value": x0_dict["rad"],
            "bounds": (150.0, 300.0),
            "active": True,
        }
    }

    # Columns -------------------------------------------------
    if opti_settings["multi_size_columns"]:
        for seg in range(1, opti_settings["n_mast_segments"] + 1):
            var[f"d0_{seg}"] = {
                "value": x0_dict["d0"],
                "bounds": column_diameter_bounds,
                "active": True,
            }
            var[f"t0_{seg}"] = {
                "value": x0_dict["t0"],
                "bounds": column_thickness_bounds,
                "active": True,
            }
    else:
        var["d0"] = {
            "value": x0_dict["d0"],
            "bounds": column_diameter_bounds,
            "active": True,
        }
        var["t0"] = {
            "value": x0_dict["t0"],
            "bounds": column_thickness_bounds,
            "active": True,
        }

    # Braces --------------------------------------------------
    if opti_settings["multi_size_braces"]:
        for seg in range(1, opti_settings["n_mast_segments"] + 1):
            if opti_settings["brace_split"]:
                var[f"d1_h_{seg}"] = {
                    "value": x0_dict["d1"],
                    "bounds": brace_diameter_bounds,
                    "active": True,
                }
                var[f"t1_h_{seg}"] = {
                    "value": x0_dict["t1"],
                    "bounds": brace_thickness_bounds,
                    "active": True,
                }
                var[f"d1_c_{seg}"] = {
                    "value": x0_dict["d1"],
                    "bounds": brace_diameter_bounds,
                    "active": True,
                }
                var[f"t1_c_{seg}"] = {
                    "value": x0_dict["t1"],
                    "bounds": brace_thickness_bounds,
                    "active": True,
                }
            else:
                var[f"d1_{seg}"] = {
                    "value": x0_dict["d1"],
                    "bounds": brace_diameter_bounds,
                    "active": True,
                }
                var[f"t1_{seg}"] = {
                    "value": x0_dict["t1"],
                    "bounds": brace_thickness_bounds,
                    "active": True,
                }
    else:
        if opti_settings["brace_split"]:
            var["d1_h"] = {
                "value": x0_dict["d1"],
                "bounds": brace_diameter_bounds,
                "active": True,
            }
            var["t1_h"] = {
                "value": x0_dict["t1"],
                "bounds": brace_thickness_bounds,
                "active": True,
            }
            var["d1_c"] = {
                "value": x0_dict["d1"],
                "bounds": brace_diameter_bounds,
                "active": True,
            }
            var["t1_c"] = {
                "value": x0_dict["t1"],
                "bounds": brace_thickness_bounds,
                "active": True,
            }
        else:
            var["d1"] = {
                "value": x0_dict["d1"],
                "bounds": brace_diameter_bounds,
                "active": True,
            }
            var["t1"] = {
                "value": x0_dict["t1"],
                "bounds": brace_thickness_bounds,
                "active": True,
            }

    # Run optimization
    RunAPDL(mapdl, var, Misc)
    result, txt_path, csv_path = optimization.run_optimization(
        mapdl, var, Misc, Solver_Settings
    )

    write_search_result(x0, result["x"], var_names)
        

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