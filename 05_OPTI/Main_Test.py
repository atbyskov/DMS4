## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# SW_Import.py          -> Reads IGS File and converts it to Keypoints and Lines for APDL
# MyAPDLCall.py         -> Inputs Swcoor, var and Misc and Runs the APDL Script
#   APDL_Input.py       -> Writes the .txt input file for APDL
# Post_Process.py       -> Reads internal forces from analysis and outputs Utilization Ratios


# Import packages
import os
import time 
from ansys.mapdl.core import launch_mapdl

# Import Functions
import SW_Import as SW
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL


# Start timing
tic = time.time()

# Import SW coordinates as list
SW_filename = "LWC1.IGS"   # Specify IGES File Name
SW_folder = "IGS"
SWcoor = SW.import_SW(os.path.join(SW_folder,SW_filename))

# Specify tube dimensions
d0 = 76.1       # Column Outer Diameter [mm]
t0 = 3          # Column Thickness      [mm]
d1 = 26.9       # Brace Outer Diameter  [mm]
t1 = 2.3        # Brace Thickness       [mm]

# Collect variables
var = [d0, t0, d1, t1] # Assemble variables

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
    log_apdl="apdl_logging",
    override=True,
    nproc=8,
    additional_switches="-p ansys -smp"
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


"""
# Print the list
print("\n--- UTILIZATION REPORT ---")
for key, val in util_list.items():
    print(f"{key:10s}  Column: {val[0]:8.4f}   Brace: {val[1]:8.4f}")
print("------------------------")
toc = time.time()

print(f"Runtime: {toc-tic:.2f} s \n")
"""