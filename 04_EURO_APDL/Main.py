## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# APDL_Eigen.py         -> Outputs .txt for Eigenbuckling Analysis
# APDL_Nonlin.py        -> Outputs .txt for Nonlinear Analysis
# UTIL.py               -> Calculates Utilization Ratios


# Import packages
import os
import shutil
import math

# Import Functions
import SW_Import as SW
from APDL_Eigen import Eigen_Fun
from APDL_Nonlin import Nonlin_Fun

# Import SW coordinates as list
SW_filename = "story"   # Specify IGES File Name
SW_folder = "IGS"
SWcoor = SW.import_SW(os.path.join(SW_folder,SW_filename))

# Specify tube dimensions
R0 = 70.1/2 # Column Tube inner diameter [mm]
R1 = 76.1/2 # Column Tube outer diameter [mm]
R2 = 22.3/2 # Brace Tube inner diameter  [mm]
R3 = 26.9/2 # Brace Tube outer diameter  [mm]

var = [R0, R1, R2, R3] # Assembly variables

# Other specifications
esize = 10              # Element Size [mm]
Hor_Force = 0           # Horizontal Force at each selection [N]
Ver_Force = 1000        # Vertical Force at each selection [N]
f_y = 690               # Yield Strength of S690 [MPa]
E_mod = 200*1E3         # Youngs Modulus [MPa]

Misc = [esize, Hor_Force, Ver_Force, f_y, E_mod]


# Call APDL EIGEN Function and get eigen file to run and CM_dict
CM_dict, eigen_file = Eigen_Fun(SWcoor,var,Misc)


# Filename for running APDL
FileNameEigen = "APDLRunFileEigen.bat"

with open(FileNameEigen, 'w') as FileID:
    FileID.write('@echo off\n')
    FileID.write('rem This batch file is placed in your working directory\n')
    FileID.write('SET ANSWAIT=1\n')
    FileID.write('set ANSYS_LOCK=OFF\n')
    FileID.write('rem set ANS_CONSEC=YES\n')
    FileID.write('"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251" -b -p ansys -smp -np 8 -i "APDL_Eigen.txt" -dir "AnsoutEigen" -o "AnsysOutputWindow.txt" \n')

# Clear Ansout folder before running 
folder = "AnsoutEigen"
if os.path.exists(folder):
    shutil.rmtree(folder)
os.makedirs(folder) # Create folder again


# Run Eigenvalue Analysis Program
os.system(FileNameEigen)

# Read First eigenvalue:
with open("AnsoutEigen/Eigenvalue1.txt") as f:
    eigenvalues = [float(line.strip()) for line in f if line.strip()]

alpha_crit = next(v for v in eigenvalues if v > 0)



# Print Information
print(f"Eigen Analysis Complete\n -> Eigenvalue 1: {alpha_crit}")

# Calculate horizontal equivalent imperfection force
h = max(point[1] for point in SWcoor)
alpha_h = 2/math.sqrt(h)

if alpha_h < 2/3:
    alpha_h = 2/3
elif alpha_h > 1:
    alpha_h = 1

alpha_m = 2 # Assumed for now
imp_ang = 1/200 * alpha_h * alpha_m
imp_force = Ver_Force*imp_ang

print(f"Eq. Imperfection Force: {imp_force:5.5f} N ")


#####  Nonlinear Analysis
Misc = [esize, Hor_Force, Ver_Force, f_y, E_mod, imp_force]   # Add Imperfection Force


CM_dict, Nonlin_file = Nonlin_Fun(SWcoor,var,Misc)            # Create .txt file for APDL




# Filename for running APDL
FileNameNonlin = "APDLRunFileNonlin.bat"

with open(FileNameNonlin, 'w') as FileID:
    FileID.write('@echo off\n')
    FileID.write('rem This batch file is placed in your working directory\n')
    FileID.write('SET ANSWAIT=1\n')
    FileID.write('set ANSYS_LOCK=OFF\n')
    FileID.write('rem set ANS_CONSEC=YES\n')
    FileID.write('"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251" -b -p ansys -smp -np 8 -i "APDL_Nonlin_Input.txt" -dir "AnsoutNonlin" -o "AnsysOutputWindowNonLin.txt" \n')

#os.system(FileNameNonlin)

#print("Nonlinear Analysis Complete")