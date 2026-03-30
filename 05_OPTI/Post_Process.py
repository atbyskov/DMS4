# Post_Process.py
    # Takes [var], [Misc], "APDL_Eigen_Internal.txt", "APDL_Nonlin_Internal" as input
    # 6 Functions for calculating:
        # 1. Local Buckling                 [LC]
        # 2. Normal Force                   [NF]
        # 3. Shear and Torsion              [ST]
        # 4. Bending, Normal and Shear      [BNS]
        # 5. Flexural and Torsional         [FT]
        # 6. Interaction                    [IN]
    # Outputs this as a .txt file for evaluation
    # Output highest utilization factor as list
    #   FORMAT:
    #   [Column NF, Brace NF]

# Import tools
from pathlib import Path
import pandas as pd
import numpy as np
import re
import math
from io import StringIO

# Function to read and parse forces 
# NOTE: I have used Copilot for most of this function, so understanding is low ...
def read_forces(filepath: str) -> pd.DataFrame:
    blocks = []
    current_member = None
    colnames = None
    buffer = []

    with open(filepath, "r") as f:
        lines = f.readlines()

    buffer = []

    for line in lines:
        # Detect new section header
        if line.startswith("NS"):
            if buffer and colnames and current_member:
                df = pd.read_csv(
                    StringIO("".join(buffer)),
                    sep=r"\s+",
                    engine="python",
                    names=colnames
                )
                if "Y_LOC" in df.columns:
                    df = df[df["Y_LOC"] <= 4080]

                df["Member"] = current_member
                blocks.append(df)
            
            buffer = []
            current_member = line.strip().replace("NS ","").replace(".","")
            continue

        if "ElemID" in line and "NF" in line:
            tokens = line.split()
            colnames = [tok for tok in tokens if "[" not in tok]
            continue

        if re.match(r"^\s*\d",line):
            buffer.append(line)

    if buffer and colnames and current_member:
        df = pd.read_csv(
            StringIO("".join(buffer)),
            sep=r"\s+",
            engine="python",
            names=colnames
        )
        df["Member"] = current_member
        blocks.append(df)

    return pd.concat(blocks, ignore_index=True)

# Local Buckling 
def Util_LC(var, Misc):
    # Import Radii
    R0, R1, R2, R3 = var

    # Import Misc
    esize, Hor_Force, Ver_Force, MomZ, MomY, f_y, E_mod = Misc

    # Poissons Ratio
    v = 0.3

    Util_LC = np.zeros(2)


    NonlinData = "AnsoutNonlin/APDL_Nonlin_Internal.txt"

    # Read Data and split into Columns and Braces
    df = read_forces(NonlinData)
    df_col = df[df["Member"].str.startswith("ColMember")].copy()
    df_brace = df[df["Member"].str.startswith("BraceMember")].copy()

    # Critical Force 
    P_cr_col = (2*math.pi*E_mod*(R1-R0)**2)/(math.sqrt(3*(1-v**2)))
    P_cr_brace = (2*math.pi*E_mod*(R3-R2)**2)/(math.sqrt(3*(1-v**2)))

    Util_LC[0] = df_col["NF"].abs().max()/P_cr_col
    Util_LC[1] = df_brace["NF"].abs().max()/P_cr_brace

    return Util_LC

# Function to calcualte Utilization ratios                
def Util_NF(var, Misc):
    # Open and Read Eigenvalue
    with open("AnsoutEigen/Eigenvalue1.txt") as f:
        eigenvalues = [float(line.strip()) for line in f if line.strip()]
    alpha_crit = next(v for v in eigenvalues if v > 0)

    # Initialize
    Util_NF = np.zeros(2)
    # Import Radii
    R0, R1, R2, R3 = var 

    # Areas
    A_Column = np.pi * ((R1**2) - (R0**2)) 
    A_Brace = np.pi * ((R3**2) - (R2**2)) 

    # Import Misc
    esize, Hor_Force, Ver_Force, MomZ, MomY, f_y, E_mod = Misc
    
    # Read Internal Force Data from APDL
    EigenData = "AnsoutEigen/APDL_Eigen_Internal.txt"
    NonlinData = "AnsoutNonlin/APDL_Nonlin_Internal.txt"

    # Read Data and split into Columns and Braces
    df = read_forces(EigenData)
    df_col = df[df["Member"].str.startswith("ColMember")].copy()
    df_brace = df[df["Member"].str.startswith("BraceMember")].copy()

    # Design Resistance
    # N = mm^2 * N/mm^2
    N_Rd_Col = A_Column*f_y 
    N_Rd_Brace = A_Brace*f_y

    # Compression/Tension Utilization [1993-1-1 "6.2.4"]
    df_col["Util_NF"] = df_col["NF"].abs() / N_Rd_Col
    df_brace["Util_NF"] = df_brace["NF"].abs() / N_Rd_Brace

    Util_NF[0] = df_col["Util_NF"].max()
    Util_NF[1] = df_brace["Util_NF"].max()

    # Normal force Utilization
    print(f"Alpha_Crit Value: {alpha_crit}")

    
    return Util_NF

def Util_ST(var, Misc):

    

    return 0

# Update


def print_info(var, Misc):
    """
    Writes the following line to Post_Util.txt:

    E_ID   Y_LOC   UTIL_LC[0]   UTIL_LC[1]   UTIL_NF[0]   UTIL_NF[1]
    """

    # Read nonlinear internal APDL output
    df = read_forces("AnsoutEigen/APDL_Eigen_Internal.txt")

    # Compute utilizations
    util_lc = Util_LC(var, Misc)   # returns [for columns, for braces]
    util_nf = Util_NF(var, Misc)   # same structure

    # Choose representative element: maximum |NF|
    idx = df["NF"].abs().idxmax()
    e_id = int(df.loc[idx, "ElemID"])
    y_loc = float(df.loc[idx, "Y_LOC"])

    # Write results to Post_Util.txt
    with open("Post_Util.txt", "w") as f:
        # Header
        f.write(f"{'E_ID':>8} {'Y_LOC':>10} {'UTIL_LC0':>12} {'UTIL_LC1':>12} "
                f"{'UTIL_NF0':>12} {'UTIL_NF1':>12}\n")

        # Values
        f.write(f"{e_id:8d} {y_loc:10.2f} "
                f"{util_lc[0]:12.5f} {util_lc[1]:12.5f} "
                f"{util_nf[0]:12.5f} {util_nf[1]:12.5f}\n")

