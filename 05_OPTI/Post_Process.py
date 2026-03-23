# Post_Process.py
    # Takes [var], [Misc], "APDL_Eigen_Internal.txt", "APDL_Nonlin_Internal" as input
    # Calculates Utilization Ratios for
        # 1. Local Buckling
        # 2. Normal Force
        # 3. Shear and Torsion
        # 4. Bending, Normal and Shear
        # 5. Flexural and Torsional
        # 6. Interaction
    # Outputs this as a .txt file for evaluation
    # Output highest utilization, place and type

from pathlib import Path
import pandas as pd
import numpy as np
import re
from io import StringIO

# Function to read and parse forces 
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
             
            
def Util_ratio(var, Misc):

    with open("AnsoutEigen/Eigenvalue1.txt") as f:
        eigenvalues = [float(line.strip()) for line in f if line.strip()]

    alpha_crit = next(v for v in eigenvalues if v > 0)

    # Initialize
    Util_list = np.zeros(2)
    # Import Radii
    R0, R1, R2, R3 = var 

    # Areas
    A_Column = np.pi * ((R1**2) - (R0**2)) 
    A_Brace = np.pi * ((R3**2) - (R2**2)) 

    # Import Misc
    esize, Hor_Force, Ver_Force, f_y, E_mod = Misc
    
    # Read Internal Force Data from APDL
    EigenData = "AnsoutEigen/APDL_Eigen_Internal.txt"
    NonlinData = "AnsoutNonlin/APDL_Nonlin_Internal.txt"

    df = read_forces(NonlinData)
    df_col = df[df["Member"].str.startswith("ColMember")].copy()
    df_brace = df[df["Member"].str.startswith("BraceMember")].copy()


    N_Rd_Col = A_Column*f_y 
    N_Rd_Brace = A_Brace*f_y

    df_col["Util_NF"] = df_col["NF"].abs() / N_Rd_Col
    df_brace["Util_NF"] = df_brace["NF"].abs() / N_Rd_Brace

    Util_list[0] = df_col["Util_NF"].max()
    Util_list[1] = df_brace["Util_NF"].max()

    # Normal force Utilization
    print(f"Alpha_Crit Value: {alpha_crit}")

    with open("UTIL_Ratio", "w") as f:
        f.write("Utilization Ratios \n")
        f.write("Util_NF Column \t Util_NF Column \n")

        for i in range(len(df_col)):
            col_val = df_col["Util_NF"].iloc[i]
            brace_val = df_brace["Util_NF"].iloc[i]

            f.write(f"{col_val:20.6f} {brace_val:20.6f}\n")

    
    return Util_list
