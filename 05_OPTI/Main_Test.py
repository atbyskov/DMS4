# Main_Test.py
# Script for Standalone Analysis Without Optimization


# Import packages
import os
import time 
from ansys.mapdl.core import launch_mapdl
import numpy as np
import pandas as pd

# Import Functions
import SW_Import as SW
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL


# Start timing
tic = time.time()


var = {
    "d0": {"value": 49.4,    "active": True},        # Column Diameter  [mm]
    "t0": {"value": 2.07,     "active": True},        # Column Thickness [mm]
    "d1": {"value": 28.88,    "active": True},        # Brace Diameter   [mm]
    "t1": {"value": 0.62,     "active": True},        # Brace Thickness  [mm]
    "rad": {"value": 314.72, "active": True},       # Radius Structure [mm]
}

# Create Misc as dict
Misc = {
    "esize": 3,                                      # Element Size [mm]
    "Hor_Force": 502.52,                             # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,                          # Vertical Force (P_Load_y)   [N]
    "f_y": 690 ,                                     # Column Yield Strength [MPa]
    "f_y_brace": 355,                                # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,                                # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,                            # Vertical Force COG (P_COG_y) [N]
    "SW_filename": "LWC_LC1.IGS"                     # Filename for IGS File
}

#C:\Program Files\ANSYS Inc\v251\ansys\bin\winx64
# Start License Server

tic_lic = time.time()
mapdl = launch_mapdl(
    run_location="Ansout", 
    log_apdl="apdl_logging",
    override=True,
    nproc=6,
    additional_switches="-p ansys -smp"
)
toc_lic = time.time()
print(f"License opened in: {toc_lic-tic_lic:.2f} s")

# Run Environment
try:
    f = RunAPDL(mapdl, var, Misc) # Runs APDL and returns MASS
finally:
    mapdl.exit()


print(f" -> Mass of Assembly: {f:.2f} kg")


# Pack specifically for compability
active = [(name,data["value"]) for name,data in var.items() if data.get("active", True)]
names = [name for name,_ in active]
x0 = np.array([var[name]["value"] for name, _ in active], dtype=float)
var_dict = dict(zip(names, x0))


utils = PostProcessor(var_dict, Misc)
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