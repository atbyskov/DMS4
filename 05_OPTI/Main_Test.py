## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# APDL_Eigen.py         -> Outputs .txt for Eigenbuckling Analysis
# APDL_Nonlin.py        -> Outputs .txt for Nonlinear Analysis
# Post_Process.py       -> Calculates Utilization Ratios


# Import packages
import os
import time 

import builtins
import numpy as np
from optimization import run_optimization

# Import Functions
import SW_Import as SW
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL



tic = time.time()
# Import SW coordinates as list
SW_filename = "LWC.IGS"   # Specify IGES File Name
SW_folder = "IGS"
SWcoor = SW.import_SW(os.path.join(SW_folder,SW_filename))

# Specify tube dimensions
R0 = 70.1/2 # Column Tube inner diameter [mm]
R1 = 76.1/2 # Column Tube outer diameter [mm]
R2 = 22.3/2 # Brace Tube inner diameter  [mm]
R3 = 26.9/2 # Brace Tube outer diameter  [mm]

# Added something here

var = [R0, R1, R2, R3] # Assemble variables

# Other specifications
esize = 100              # Element Size [mm]
Hor_Force = 502.52      # Horizontal Force [N]
Ver_Force = 26400       # Vertical Force   [N]
MomZ = -70364000        # Applied Moment around Z-axis [Nmm]
MomY = 1407140          # Applied Moment around Y-axis [Nmm]
f_y = 690               # Yield Strength of S690 [MPa]
E_mod = 200*1E3         # Youngs Modulus [MPa]

# Create Misc as dict
Misc = {
    "esize": esize,
    "Hor_Force": Hor_Force,
    "Ver_Force": Ver_Force,
    "MomZ": MomZ,
    "MomY": MomY,
    "f_y": f_y,
    "E_mod": E_mod
}

f = RunAPDL(SWcoor,var,Misc) # Runs APDL and returns MASS


print(f"Mass of Assembly: {f:.2f} kg")

utils = PostProcessor()
util_list = utils.Util_list(var,Misc)

# Print the list
print("\n--- UTILIZATION REPORT ---")
for key, val in util_list.items():
    print(f"{key:10s}  Column: {val[0]:8.4f}   Brace: {val[1]:8.4f}")
print("------------------------")
toc = time.time()

print(f"Runtime: {toc-tic:.2f} s \n")