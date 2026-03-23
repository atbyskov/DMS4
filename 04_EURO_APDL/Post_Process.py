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

import os

def Util_ratio(var, Misc):
    
    # Read Internal Force Data from APDL
    EigenData = "AnsoutEigen/APDL_Eigen_Internal.txt"
    NonlinData = "AnsoutNonlin/APDL_Nonlin_Internal.txt"

    # Read EigenData
    with open(EigenData,"r") as f:
        f.readlines


    ss = 22
    
    with open("UTIL_Ratio", "w") as f:
        f.write("Utilization Ratio for Model \n")

    
    return ss
