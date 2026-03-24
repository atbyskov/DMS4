# MyAPDLCall.py
# Runs APDL and outputs weight of assembly

import os
import math
import shutil

from APDL_Eigen import Eigen_Fun
from APDL_Nonlin import Nonlin_Fun

def RunAPDL(SWcoor,var,Misc):

    

    esize, Hor_Force, Ver_Force, f_y, E_mod = Misc


    # Clear Ansout folder before running 
    eigen_dir = "AnsoutEigen"
    nonlin_dir = "AnsoutNonlin"
    for folder in (eigen_dir, nonlin_dir):
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder, exist_ok=True) # Create folder again

    eigen_file = Eigen_Fun(SWcoor,var,Misc, out_dir = eigen_dir)

    # Filename for running APDL
    FileNameEigen = os.path.join(eigen_dir,"APDLRunFileEigen.bat")

    with open(FileNameEigen, 'w') as FileID:
        FileID.write('@echo off\n')
        FileID.write('rem This batch file is placed in your working directory\n')
        FileID.write('SET ANSWAIT=1\n')
        FileID.write('set ANSYS_LOCK=OFF\n')
        FileID.write('rem set ANS_CONSEC=YES\n')
        FileID.write(
            f'"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251"'
            f' -b -p ansys -smp -np 8'
            f' -i {eigen_file}'
            f' -dir "{eigen_dir}"'
            f' -o {os.path.join(eigen_dir,"AnsysOutputWindow.txt")} \n')

    # Run Eigenvalue Analysis Program
    os.system(f'"{FileNameEigen}"')

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

    Nonlin_file = Nonlin_Fun(SWcoor,var,Misc,imp_force,out_dir=nonlin_dir)            # Create .txt file for APDL

    # Filename for running APDL
    FileNameNonlin = os.path.join(nonlin_dir,"APDLRunFileNonlin.bat")

    with open(FileNameNonlin, 'w') as FileID:
        FileID.write('@echo off\n')
        FileID.write('rem This batch file is placed in your working directory\n')
        FileID.write('SET ANSWAIT=1\n')
        FileID.write('set ANSYS_LOCK=OFF\n')
        FileID.write('rem set ANS_CONSEC=YES\n')
        FileID.write(
            f'"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251"'
            f' -b -p ansys -smp -np 8'
            f' -i {Nonlin_file}'
            f' -dir "{nonlin_dir}"'
            f' -o {os.path.join(nonlin_dir,"AnsysOutputWindow.txt")} \n')

    # Run Nonlinear Analysis
    os.system(f"{FileNameNonlin}")
    print("Nonlinear Analysis Complete")


    
    with open("AnsoutEigen/Mass_Assembly.txt","r") as f:
        Mass = [float(line.strip()) for line in f if line.strip()]


    return Mass