## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py           -> 
# SW_Import.py      -> Used for reading IGES file and converting to APDL Geometry
# APDL_pre          -> Writes all APDL code expect post processing
# APDL_post         -> Write all APDL Post processing code
# Post_Process.py   -> Custom .py for post-processing data
# Plot_tool.py      -> Custom plotting tool

# Import packages
import os
import shutil


# Import Functions
from APDL_EIGEN import Eigen_Fun
from APDL_post import APDL_post_fun
import SW_Import as SW

# Import SW coordinates as list
SW_filename = "SimpleFrame.IGS"   # Specify IGES File Name
SW_folder = "IGS"
SWcoor = SW.import_SW(os.path.join(SW_folder,SW_filename))

# Specify tube dimensions
R0 = 70.1/2 # Corner Tube inner diameter [mm]
R1 = 76.1/2 # Corner Tube outer diameter [mm]
R2 = 22.3/2 # Brace Tube inner diameter  [mm]
R3 = 26.9/2 # Brace Tube outer diameter  [mm]

var = [R0, R1, R2, R3] # Assembly variables

# Other specifications
esize = 10          # Element Size [mm]
Hor_Force = 260     # Horizontal Force at each selection [N]
Ver_Force = 100   # Vertical Force at each selection [N]
f_y = 690           # Yield Strength of S690 [MPa]
E_mod = 200*1E3     # Youngs Modulus [MPa]

Misc = [esize, Hor_Force, Ver_Force, f_y, E_mod]



# Call APDL Pre Function and get CM_dict = [nColumns, nBraces]
CM_dict, pre_file = Eigen_Fun(SWcoor,var,Misc)

# Call APDL Post Function
post_file = APDL_post_fun(CM_dict)

# Create Runfile
RunFile = "Runfile.txt"

# Input pre_file and post_file into the Runfile.txt
with open(RunFile,"w") as f:
    with open(pre_file, "r") as f1:
        f.write(f1.read())
    #with open(post_file, "r") as f2:
    #    f.write(f2.read())

# Filename for running APDL
FileName = "APDLRunFile.bat"

with open(FileName, 'w') as FileID:
    FileID.write('@echo off\n')
    FileID.write('rem This batch file is placed in your working directory\n')
    FileID.write('SET ANSWAIT=1\n')
    FileID.write('set ANSYS_LOCK=OFF\n')
    FileID.write('rem set ANS_CONSEC=YES\n')
    FileID.write('"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251" -b -p ansys -smp -np 8 -i "RunFile.txt" -dir "Ansout" -o "AnsysOutputWindow.txt" \n')

# Clear Ansout folder before running 
#folder = "Ansout"
#if os.path.exists(folder):
#    shutil.rmtree(folder)
#os.makedirs(folder) # Create folder again


# Run Program
os.system(FileName)

# Read eigenvalue:
with open("Ansout/Eigenvalue1.txt","r") as f:
    MS1 = f.readline().strip()

# Print Information
print(f"Eigen Analysis Complete\n -> Eigenvalue 1: {MS1}")

## Post Processing ## 
# Call