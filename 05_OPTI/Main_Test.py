## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# SW_Import.py          -> Reads IGS File and converts it to Keypoints and Lines for APDL
# MyAPDLCall.py         -> Inputs Swcoor, var and Misc and Runs the APDL Script
#   APDL_Input.py       -> Writes the .txt input file for APDL
# Post_Process.py       -> Reads internal forces from analysis and outputs Utilization Ratios


# Import packages
import os
import time 

# Import Functions
import SW_Import as SW
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL

from ansys.mapdl.core import launch_mapdl

# Start timing
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

# Other specifications``
esize = 100             # Element Size [mm]
Hor_Force = 502.52      # H`orizontal Force [N]
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
#C:\Program Files\ANSYS Inc\v251\ansys\bin\winx64
# Start License Server

tic_lic = time.time()
mapdl = launch_mapdl(
    run_location="Ansout", 
    override=True,
    additional_switches="-p ansys"
)
toc_lic = time.time()
print(f"License opened in: {toc_lic-tic_lic:.2f} s")



# Run Environment
try:
    f = RunAPDL(mapdl, SWcoor, var, Misc) # Runs APDL and returns MASS
finally:
    mapdl.exit()


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