## MAIN BEAM ELEMENT DOCUMENT ##

import sys
print(sys.version)

# Import packages
import time 
from ansys.mapdl.core import launch_mapdl
from pyslsqp import optimize

# Import Functions
import optimization
from MyAPDLCall import RunAPDL


# Variables. Uncomment variable to include it.
var = {
    "d0": 76.1,                     # Column Diameter  [mm]
    "t0": 3,                        # Column Thickness [mm]
    "d1": 26.9,                     # Brace Diameter   [mm]
    "t1": 2.3,                      # Brace Thickness  [mm]
    #"rad": 202.07                  # Radius Structure [mm]
}

# Static variables
Misc = {
    "esize": 100,                   # Element Size [mm]
    "Hor_Force": 502.52,            # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)   [N]
    "f_y": 690 ,                    # Column Yield Strength [MPa]
    "f_y_brace": 355,               # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
    "SW_filename": "LWC_LC1.IGS"    # Filename for IGS File
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
    result, txt_path, csv_path = optimization.run_optimization(mapdl, var, Misc)
finally:
    mapdl.exit()

# Print Information
#print("\nOptimal x:", result.x)
#print("Optimal objective:", result.fun)
#print("Message:", result.message)
#print("TXT log file:", txt_path)
#print("Objective CSV:", csv_path)

