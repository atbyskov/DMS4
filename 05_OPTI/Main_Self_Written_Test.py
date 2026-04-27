## MAIN BEAM ELEMENT DOCUMENT ##

import sys
print(sys.version)

# Import packages
import time 
from ansys.mapdl.core import launch_mapdl
from pyslsqp import optimize

# Import Functions
import optimization
import Self_Written_Optimization.optimization_self_written as optimization_self_written
from MyAPDLCall import RunAPDL


# Variables. Choose which ones to include by setting "active": True or False.
# Include Bounds
var = {
    "d0": {"value": 76.1, "bounds": (40.0, 100), "active": True},       # Column Diameter  [mm]
    "t0": {"value": 3.0,  "bounds": (1.0, 7.0),  "active": True},       # Column Thickness [mm]
    "d1": {"value": 26.9, "bounds": (10.0, 100.0), "active": True},     # Brace Diameter   [mm]
    "t1": {"value": 2.3,  "bounds": (0.1, 7.0),  "active": True},       # Brace Thickness  [mm]
    "rad": {"value": 202.07, "bounds": (150.0, 350.0), "active": True}, # Radius Structure [mm]
}

# Static variables
Misc = {
    "esize": 3,                     # Element Size [mm]
    "Hor_Force": 502.52,            # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)   [N]
    "f_y": 690 ,                    # Column Yield Strength [MPa]
    "f_y_brace": 355,               # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
    "eps_geom": 0.1,    # Minimum thickness specification for geometry updates [mm]
    "SW_filename": "LWC_LC1.IGS",    # Filename for IGS File
    "save_folder": "Optimization_Logs" # Save Folder
}

# Launch MAPDL
mapdl = launch_mapdl(
    run_location="Ansout",
    log_apdl="apdl_log",
    override=True,
    nproc=8,
    additional_switches="-p ansys -smp",
)

# Ensures that MAPDL closes if something chrashes
try:
    f1 = RunAPDL(mapdl,var,Misc)
    result, txt_path, csv_path = optimization_self_written.run_optimization(mapdl, var, Misc, method="slp_mvp")
finally:
    mapdl.exit()

# Print Information
#print("\nOptimal x:", result.x)
#print("Optimal objective:", result.fun)
#print("Message:", result.message)
#print("TXT log file:", txt_path)
#print("Objective CSV:", csv_path)

