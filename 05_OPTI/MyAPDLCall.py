# MyAPDLCall.py
# When called, this document runs Eigenbuckling Analysis and Nonlinear and Returns mass

# -> INPUT:
#       [SWcoor]    -> Coordinates from Solidworks IGES
#       [var]       -> Radii variables
#       [Misc]      -> Miscellaneous Data (force, mesh etc.) 
# -> OUTPUT:
#       [Mass]      -> Mass of total assembly

# Pseudo code
#   Import [SWcoor], [var] and [Misc]
#   Remove content from "AnsoutEigen" and "AnsoutNonlin" folders
#   Call APDL_Eigen to create .txt input file
#       1. Run Analysis via os.system
#       2. Read first positive eigenvalue
#       3. Calculate Imperfection Force
#   Call APDL_Nonlin to create .txt input file
#       1. Add Imperfection force as input 
#       2. Run Analysis via os.system
#       3. Read Mass of assembly
#   Return Mass
        
# Import Tools
import time
import os

# Import Functions
from APDL_Input import InputFun


def RunAPDL(mapdl,var,Misc):
    import numpy as np
    ans_time_tic = time.time()
   
    if isinstance(var, dict):
        x = np.array([var["d0"], var["t0"], var["d1"], var["t1"]], dtype=float)
    else:
        x = np.asarray(var, dtype=float).ravel()

    var_dict = {
        "d0": x[0],
        "t0": x[1],
        "d1": x[2],
        "t1": x[3],
    }


    mapdl.clear()

    # Create input file for Eigenvalue Analysis
    apdl_cmds = InputFun(var_dict,Misc)

    with mapdl.non_interactive:
        for cmd in apdl_cmds:
            cmd = cmd.strip()
            if cmd:
                mapdl.run(cmd)

    # Clear APDL
    mapdl.finish()
    

    # Read First eigenvalue:
    with open("Ansout/Eigenvalue1.txt") as f:
        eigenvalues = [float(line.strip()) for line in f if line.strip()]
    # Retrieve first positive eigenvalue
    alpha_crit = next(v for v in eigenvalues if v > 0) 

    # Print Information
    print(f"Analysis Complete\n -> Eigenvalue 1: {alpha_crit:.2f}")

    # Open and Read Mass
    with open("Ansout/Mass_Assembly.txt","r") as f:
        Mass = [float(line.strip()) for line in f if line.strip()]

    ans_time_toc = time.time()

    print(f"Sim time: {ans_time_toc-ans_time_tic:.2f} s")

    # Return Mass as float value
    return sum(Mass)