## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# APDL_Eigen.py         -> Outputs .txt for Eigenbuckling Analysis
# APDL_Nonlin.py        -> Outputs .txt for Nonlinear Analysis
# Post_Process.py       -> Calculates Utilization Ratios


# Import packages
import os
import time 

# Import Functions
import SW_Import as SW
from Post_Process import Util_NF
from Post_Process import Util_LC
from Post_Process import print_info
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

Misc = [esize, Hor_Force, Ver_Force, MomZ, MomY, f_y, E_mod]


f = RunAPDL(SWcoor,var,Misc) # Runs APDL and returns MASS


print(f"Mass of Assembly: {f} kg")



# Run Post_Process.py
Util_LC = Util_LC(var,Misc)
Util_NF = Util_NF(var,Misc)

print(f"Util_LC: {Util_LC}")
print(f"Util_NF: {Util_NF}")

toc = time.time()

runtime = toc-tic
print(f"Sim Time: {runtime} s")

print = print_info(var,Misc)

#Hello
