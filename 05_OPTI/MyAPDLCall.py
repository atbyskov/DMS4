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
import os
import math
import shutil

# Import Functions
from APDL_Input import InputFun


def RunAPDL(SWcoor,var,Misc):


    # Clear Ansout folder before running 
    Ansout_dir = "Ansout"
    if os.path.exists(Ansout_dir):
        shutil.rmtree(Ansout_dir)
    os.makedirs(Ansout_dir, exist_ok=True) # Create folder again

    # Create input file for Eigenvalue Analysis
    input_file = InputFun(SWcoor,var,Misc, out_dir = Ansout_dir)

    # Filename for running APDL .bat file
    FileName = os.path.join(Ansout_dir,"APDLRunFile.bat")
    # Open and Edit .bat file
    with open(FileName, 'w') as FileID:
        FileID.write('@echo off\n')
        FileID.write('rem This batch file is placed in your working directory\n')
        FileID.write('SET ANSWAIT=1\n')
        FileID.write('set ANSYS_LOCK=OFF\n')
        FileID.write('rem set ANS_CONSEC=YES\n')
        FileID.write(
            f'"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251"'
            f' -b -p ansys -smp -np 6' 
            f' -i {input_file}'
            f' -dir "{Ansout_dir}"'
            f' -o {os.path.join(Ansout_dir,"AnsysOutputWindow.txt")} \n')

    # Run Eigenvalue Analysis Program
    os.system(f'"{FileName}"')

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

    # Return Mass as float value
    return sum(Mass)