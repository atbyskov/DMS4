
## MAIN BEAM ELEMENT DOCUMENT ##

# Import packages
import sys
print(sys.version, flush=True)
import time
from ansys.mapdl.core import launch_mapdl

# Import Functions
import optimization
from MyAPDLCall import RunAPDL

# -----------------------------
# SETTINGS
# -----------------------------

opti_settings = {
    "n_mast_segments": 5,
    "mast_segment_height": 810,
    "segment_mass_limit": 23,
    "multi_size_columns": False,
    "multi_size_braces": False,
    "brace_split": False,
}

# Initial Guess (not really used now, but kept)
column_diameter = 76.1
column_thickness = 3
brace_diameter = 26.9
brace_thickness = 2.3

# Bounds
column_diameter_bounds = (48.3, 114.3)   # Column Diameter Bounds [mm]
column_thickness_bounds = (2.5, 5.0)     # Column Thickness Bounds [mm]
brace_diameter_bounds = (10.0, 50.0)     # Brace Diameter Bounds [mm]
brace_thickness_bounds = (1.0, 4.0)      # Brace Thickness Bounds [mm]

# -----------------------------
# VARIABLE DEFINITION
# -----------------------------
var = {
    "rad": {"value": 202.07, "bounds": (150.0, 350.0), "active": True}, # Radius Structure [mm]
}

# Columns
if opti_settings["multi_size_columns"]:
    var.update({
        **{f"d0_{i}": {"value": column_diameter, "bounds": column_diameter_bounds, "active": True}
           for i in range(1, opti_settings["n_mast_segments"]+1)},
        **{f"t0_{i}": {"value": column_thickness, "bounds": column_thickness_bounds, "active": True}
           for i in range(1, opti_settings["n_mast_segments"]+1)},
    })
else:
    var.update({
        "d0": {"value": column_diameter, "bounds": column_diameter_bounds, "active": True},
        "t0": {"value": column_thickness, "bounds": column_thickness_bounds, "active": True},
    })

# Braces
if opti_settings["multi_size_braces"]:
    if opti_settings["brace_split"]:
        var.update({
            **{f"d1_h_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True}
               for i in range(1, opti_settings["n_mast_segments"]+1)},
            **{f"t1_h_{i}": {"value": brace_thickness, "bounds": brace_thickness_bounds, "active": True}
               for i in range(1, opti_settings["n_mast_segments"]+1)},
            **{f"d1_c_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True}
               for i in range(1, opti_settings["n_mast_segments"]+1)},
            **{f"t1_c_{i}": {"value": brace_thickness, "bounds": brace_thickness_bounds, "active": True}
               for i in range(1, opti_settings["n_mast_segments"]+1)},
        })
    else:
        var.update({
            **{f"d1_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True}
               for i in range(1, opti_settings["n_mast_segments"]+1)},
            **{f"t1_{i}": {"value": brace_thickness, "bounds": brace_thickness_bounds, "active": True}
               for i in range(1, opti_settings["n_mast_segments"]+1)},
        })
else:
    if opti_settings["brace_split"]:
        var.update({
            "d1_h": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},
            "t1_h": {"value": brace_thickness, "bounds": brace_thickness_bounds, "active": True},
            "d1_c": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},
            "t1_c": {"value": brace_thickness, "bounds": brace_thickness_bounds, "active": True},
        })
    else:
        var.update({
            "d1": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},
            "t1": {"value": brace_thickness, "bounds": brace_thickness_bounds, "active": True},
        })

# -----------------------------
# STATIC VARIABLES
# -----------------------------
Misc = {
    "esize": 3,
    "Hor_Force": 502.52,
    "Ver_Force": -25.13E+3,
    "f_y": 700,
    "f_y_brace": 235,
    "E_mod": 200*1E3,
    "W_Force": -3.751E+3,
    "SW_filename": "LWC_L1_LINES.IGS",
    "save_folder": "Optimization_Logs"
}
# Solver Settings
Solver_Settings = {
    "acc": 1e-3,                 # Maximum objective function tolerance
    "maxiter": 120,               # Maximum iterations
    "Aggregate": None,       # None, "P-norm", "P-norm-mean", "KS", "KS_shift"  (Write exacly)
    "p_value": 8,               # Value for "P-norm" and "P-norm-mean"
    "rho_value": 100,            # rho value used in KS
    "relaxation": 0,             # Relaxation parameter used in aggregation
    "use_acs": False,              # Toggle ACS to on (True) or off (False)
}


# FUNCTION: READ SEARCH SPACE
def read_search_space(filepath):
    points = []

    with open(filepath, "r") as f:
        lines = f.readlines()

    start = False

    for line in lines:
        # Normalize spaces
        clean_line = " ".join(line.strip().split())

        # Detect header robustly
        if clean_line.startswith("rad d0 t0 d1 t1"):
            start = True
            continue

        if start:
            parts = clean_line.split()

            # Only accept valid numeric rows
            if len(parts) == 5:
                try:
                    values = [float(x) for x in parts]
                    points.append(values)
                except ValueError:
                    # skip malformed lines
                    continue

    return points

search_points = read_search_space("search_space.txt")
print(f"Loaded {len(search_points)} simulation points")

results = []

for i, point in enumerate(search_points):
    print(f"\n--- Running simulation {i+1}/{len(search_points)} ---")
    tic_lic = time.time()
    mapdl = launch_mapdl(
        run_location="Ansout",
        log_apdl="apdl_log",
        override=True,
        nproc=6,
        additional_switches="-p ansys -smp",
)
    toc_lic = time.time()
    print(f"License opened in: {toc_lic - tic_lic:.2f} s")

    rad_val, d0_val, t0_val, d1_val, t1_val = point

    # Update design variables
    var["rad"]["value"] = rad_val
    var["d0"]["value"] = d0_val
    var["t0"]["value"] = t0_val
    var["d1"]["value"] = d1_val
    var["t1"]["value"] = t1_val
    print(f"rad: {var["rad"]["value"]}\nd0: {var["d0"]["value"]}\nt0: {var["t0"]["value"]}\nd1: {var["d1"]["value"]}\nt1: {var["t1"]["value"]}\n  " )

    # Run simulation
    try:
        result = RunAPDL(mapdl, var, Misc, opti_settings)
        result, txt_path, csv_path = optimization.run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings)
    finally:
        mapdl.exit()

print("\n✅ All simulations completed.")

# -----------------------------
# CLOSE MAPDL (OUTSIDE LOOP ✅)
# -----------------------------


