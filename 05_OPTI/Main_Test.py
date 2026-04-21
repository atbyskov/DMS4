# Main_Test.py
# Script for Standalone Analysis Without Optimization


# Import packages
import time 
import numpy as np
import pandas as pd

# PyMAPDL Package
from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core.errors import MapdlRuntimeError, MapdlException

# Import Functions
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL


# Start timing
tic = time.time()

opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "multi_size_columns": True,   # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": True,    # Whether mast segments braces uses different dimensions (True) or not (False)
}

# Collect variables
# Defining variables with bounds and active status
var = {
    "rad": {"value": 202.07, "bounds": (150.0, 300.0), "active": True}, # Radius Structure [mm]
}
if opti_settings["multi_size_columns"]:
    var.update({
        **{f"d0_{i}": {"value": 76.1, "bounds": (50.0, 100), "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Diameter  [mm]
        **{f"t0_{i}": {"value": 3.0,  "bounds": (1.0, 7.0),  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Thickness [mm]
    })
else:
    var.update({
        "d0": {"value": 76.1, "bounds": (50.0, 100), "active": True},       # Column Diameter  [mm]
        "t0": {"value": 3.0,  "bounds": (1.0, 7.0),  "active": True},       # Column Thickness [mm]
    })
if opti_settings["multi_size_braces"]:
    var.update({
        **{f"d1_{i}": {"value": 26.9, "bounds": (10.0, 40.0), "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},     # Brace Diameter   [mm]
        **{f"t1_{i}": {"value": 2.3,  "bounds": (0.1, 4.5),  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Brace Thickness  [mm]
    })
else:
    var.update({
        "d1": {"value": 26.9, "bounds": (10.0, 40.0), "active": True},     # Brace Diameter   [mm]
        "t1": {"value": 2.3,  "bounds": (0.1, 4.5),  "active": True},       # Brace Thickness  [mm]
    })
    
# Create Misc as dict
Misc = {
    "esize": 3,                     # Element Size                      [mm]
    "Hor_Force": 502.52,            # Horizontal Force (P_Load_z)       [N]
    "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)         [N]
    "f_y": 690 ,                    # Column Yield Strength             [MPa]
    "f_y_brace": 355,               # Brace Yield Strength              [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus                    [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y)      [N]
    "SW_filename": "LWC_L1_LINES.IGS"    # Filename for IGS File
}

#C:\Program Files\ANSYS Inc\v251\ansys\bin\winx64
# Start License Server

tic_lic = time.time()
mapdl = launch_mapdl(
    run_location="Ansout", 
    log_apdl="apdl_logging",
    override=True,
    nproc=2,
    additional_switches="-p ansys -smp"
)
toc_lic = time.time()
print(f"License opened in: {toc_lic-tic_lic:.2f} s")

# Run Environment
try:
    f = RunAPDL(mapdl, var, Misc, opti_settings) # Runs APDL and returns MASS
except MapdlRuntimeError as e:
    # Handle MAPDL runtime errors
    print(f"MAPDL error occured: {e}")
except Exception as e:
    # Handle other errors
    print(f"Unexpected error: {e}")
finally:
    mapdl.exit()

# Print mass of assembly
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

    print(f"{key:10s}  Column: {col_str}  Brace: {brc_str}")