# Main_Test.py
# Script for Standalone Analysis Without Optimization
# Main_Test.py
# Script for Standalone Analysis Without Optimization

# Import packages
import sys
print(sys.version, flush=True)
import time 
import numpy as np
import pandas as pd
from functools import wraps

# Start timing
tic_lic = time.time()

# PyMAPDL Package
from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core.errors import MapdlRuntimeError

# Import Functions
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL



# Collect opti settings
opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "multi_size_columns": False,   # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": False,    # Whether mast segments braces uses different dimensions (True) or not (False)
    "brace_split": False,          # Whether braces are split between horiontal and cross (True) or not (False)
}
# x0 = [194.3214, 101.0977,   3.5106,  54.8785,   2.8011]
# Initial Guess
column_diameter = 101.0977 # Column Diameter [mm]
column_thickness = 3.5106 # Column Thickness [mm]
brace_diameter = 54.8785 # Brace Diameter [mm]
brace_thickness = 2.8011 # Brace Thickness [mm]
# Bounds
column_diameter_bounds = (50.0, 300)   # Column Diameter Bounds [mm]
column_thickness_bounds = (1.0, 7.0)   # Column Thickness Bounds [mm]
brace_diameter_bounds = (10.0, 50.0)   # Brace Diameter Bounds [mm]
brace_thickness_bounds = (1.0, 4.5)    # Brace Thickness Bounds [mm]

# Defining variables with bounds and active status
var = {
    "rad": {"value": 194.3214, "bounds": (150.0, 300.0), "active": True}, # Radius Structure [mm]
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
            **{f"d1_h_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},           # Horizontal Brace Diameter   [mm]
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
    
# Create Misc as dict
Misc = {
    "esize": 3,                          # Element Size                      [mm]
    "Hor_Force": 502.52,                 # Horizontal Force (P_Load_z)       [N]
    "Ver_Force": -25.13E+3,              # Vertical Force (P_Load_y)         [N]
    "f_y": 700 ,                         # Column Yield Strength             [MPa]
    "f_y_brace": 235,                    # Brace Yield Strength              [MPa]
    "E_mod": 200*1E3,                    # Youngs Modulus                    [MPa]
    "W_Force": -3.751E+3,                # Vertical Force COG (P_COG_y)      [N]
    "SW_filename": "LWC_L1_LINES.IGS"    # Filename for IGS File
}

#C:\Program Files\ANSYS Inc\v251\ansys\bin\winx64
# Start License Server
time_license = time.time()
mapdl = launch_mapdl(
    run_location="Ansout", 
    #log_apdl="apdl_logging_test",
    override=True,
    nproc=6,
    additional_switches="-p ansys -smp"
)
time_license2 = time.time()
print(f"Total Time: {time_license2-time_license:.4f} s")


# Run Environment
f = None  # Initialize f to None
try:
    f, segment_masses = RunAPDL(mapdl, var, Misc, opti_settings) # Runs APDL and returns MASS
except MapdlRuntimeError as e:
    # Handle MAPDL runtime errors
    print(f"MAPDL error occurred: {e}")
except Exception as e:
    # Handle other errors
    print(f"Unexpected error: {e}")
finally:
    mapdl.exit()

print(f" -> Mass of Assembly: {f:.2f} kg")

# Pack specifically for compability

active = [(name,data["value"]) for name,data in var.items() if data.get("active", True)]
names = [name for name,_ in active]
x0 = np.array([var[name]["value"] for name, _ in active], dtype=float)
var_dict = dict(zip(names, x0))

utils = PostProcessor(var_dict, Misc, opti_settings)
Util_list = utils.Util_list()  


print("\n--- UTILIZATION REPORT ---")

def extract_max(val):
    if val is None:
        return None

    if isinstance(val, (pd.Series, pd.DataFrame)):
        if len(val) == 0:
            return None
        return np.nanmax(val.to_numpy(dtype=float))

    if np.isscalar(val):
        return float(val)

    return None

for key, util in Util_list.items():

    col_val = None
    brc_val = None

    # Case 1: tuple → (Column, Brace)
    if isinstance(util, tuple):
        col_val = extract_max(util[0])
        brc_val = extract_max(util[1])

    # Case 2: single value (fx Util_BS)
    else:
        brc_val = extract_max(util)

    # Formatting
    col_str = f"{col_val:8.3f}" if col_val is not None else "   N/A  "
    brc_str = f"{brc_val:8.3f}" if brc_val is not None else "   N/A  "

    print(f"{key:15s}  Column: {col_str}  Brace: {brc_str}")

toc_lic = time.time()
print(f"Total Time: {toc_lic-tic_lic:.2f} s")