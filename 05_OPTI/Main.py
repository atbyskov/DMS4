## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# APDL_Eigen.py         -> Outputs .txt for Eigenbuckling Analysis
# APDL_Nonlin.py        -> Outputs .txt for Nonlinear Analysis
# Post_Process.py       -> Calculates Utilization Ratios


# Import packages
import os
import time 
from ansys.mapdl.core import launch_mapdl

# Import Functions
from optimization import run_optimization
import SW_Import as SW
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL

# Start global timing
tic = time.time()

# Import SW coordinates as list
SW_filename = "LWC.IGS"   # Specify IGES File Name
SW_folder = "IGS"
SWcoor = SW.import_SW(os.path.join(SW_folder,SW_filename))

# Specify tube dimensions
d0 = 76.1       # Column Outer Diameter [mm]
t0 = 3          # Column Thickness      [mm]
d1 = 26.9       # Brace Outer Diameter  [mm]
t1 = 2.3        # Brace Thickness       [mm]

# Collect variables
var = [d0, t0, d1, t1] # Assemble variables

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

# Launch MAPDL
mapdl = launch_mapdl(
    run_location="Ansout",
    # Uncomment if you want a complete .txt log of everything written to MAPDL
    log_apdl="apdl_log",
    override=True,
    nproc=8,
    additional_switches="-p ansys -smp",
)

# Ensures that MAPDL closes if something chrashes
try:
    f1 = RunAPDL(mapdl,SWcoor,var,Misc)
    result, txt_path, csv_path = run_optimization(mapdl, var, SWcoor, Misc)
finally:
    mapdl.exit()

# Print Information
#print("\nOptimal x:", result.x)
#print("Optimal objective:", result.fun)
#print("Message:", result.message)
#print("TXT log file:", txt_path)
#print("Objective CSV:", csv_path)

