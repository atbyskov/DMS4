## MAIN BEAM ELEMENT DOCUMENT ##
# Pseudo code
    # 1. Read Solidworks Geometric Data
    # 2. Convert Coordinates to APDL language
    #   2.1 Create nodes and lines accordingly
    # 3. Set up boundary conditions in APDL
    # 4. Solve the system
    # 5. Post process Results


# Import Functions
from Input import InputFun
from Output import OutputFun
import SW_Import as SW
import os


# Import SW coordinates as list
SW_filename = "SimpleFrame.IGS"
SWcoor = SW.import_SW(SW_filename)

# Specify tube dimensions
R0 = 70.1/2 # Corner Tube inner diameter [mm]
R1 = 76.1/2 # CornSer Tube outer diameter [mm]
R2 = 22.3/2 # Brace Tube inner diameter  [mm]
R3 = 26.9/2 # Brace Tube outer diameter  [mm]

var = [R0, R1, R2, R3] # Assembly variables

esize = 10

# Call Input function and get CM_dict = [nColumns, nBraces]
CM_dict, input_file = InputFun(SWcoor,var,esize)

# Call Output function
output_file = OutputFun(CM_dict)

RunFile = "Runfile.txt"

with open(RunFile,"w") as f:

    with open(input_file, "r") as f1:
        f.write(f1.read())

    with open(output_file, "r") as f2:
        f.write(f2.read())


FileName = "APDLRunFile.bat"

with open(FileName, 'w') as FileID:
    FileID.write('@echo off\n')
    FileID.write('rem This batch file is placed in your working directory\n')
    FileID.write('SET ANSWAIT=1\n')
    FileID.write('set ANSYS_LOCK=OFF\n')
    FileID.write('rem set ANS_CONSEC=YES\n')
    FileID.write('"C:\\Program files\\ANSYS Inc\\v251\\ANSYS\\bin\\winx64\\ansys251" -b -p ansys -smp -np 8 -i "RunFile.txt" -o "AnsysOutputWindow.txt"\n')

os.system(FileName)


# Post Processing 