## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# APDL_Eigen.py         -> Outputs .txt for Eigenbuckling Analysis
# APDL_Nonlin.py        -> Outputs .txt for Nonlinear Analysis
# Post_Process.py       -> Calculates Utilization Ratios


# Import packages
import os
import shutil
import math

# Import Functions
import SW_Import as SW
from APDL_Eigen import Eigen_Fun
from APDL_Nonlin import Nonlin_Fun
from Post_Process import Util_ratio
from MyAPDLCall import RunAPDL
import numpy as np
from scipy.optimize import minimize
import scipy.optimize as spo

# Import SW coordinates as list
SW_filename = "story.IGS"   # Specify IGES File Name
SW_folder = "IGS"
SWcoor = SW.import_SW(os.path.join(SW_folder,SW_filename))

# Specify tube dimensions
R0 = 70.1/2 # Column Tube inner diameter [mm]
R1 = 76.1/2 # Column Tube outer diameter [mm]
R2 = 22.3/2 # Brace Tube inner diameter  [mm]
R3 = 26.9/2 # Brace Tube outer diameter  [mm]

var = [R0, R1, R2, R3] # Assemble variables

# Other specifications
esize = 10              # Element Size [mm]
Hor_Force = 0           # Horizontal Force at each selection [N]
Ver_Force = 10000       # Vertical Force at each selection [N]
f_y = 690               # Yield Strength of S690 [MPa]
E_mod = 200*1E3         # Youngs Modulus [MPa]

Misc = [esize, Hor_Force, Ver_Force, f_y, E_mod]


#f = RunAPDL(SWcoor,var,Misc) # Runs APDL and returns MASS

#print(f"Mass of Assembly: {f} kg")

# Run Post_Process.py
#Util_worst = Util_ratio(var,Misc)
#print(Util_worst)


x0 = np.array(var, dtype=float)

def objective(x):
    print(f"\nSciPy tries x = {x}")
    value = RunAPDL(SWcoor, x, Misc)
    print(f"Objective = {value}")
    return float(value)

# Bounds for each variable
bounds = [
    (30.0, 60.0),   # R0
    (31.0, 70.0),   # R1
    (8.0,  20.0),   # R2
    (9.0,  25.0)    # R3
]

# Geometric constraints: outer radius must exceed inner radius
cons = [
    {"type": "ineq", "fun": lambda x: x[1] - x[0] - 1.0},  # column thickness >= 1 mm
    {"type": "ineq", "fun": lambda x: x[3] - x[2] - 1.0},  # brace thickness >= 1 mm
]

result = spo.minimize(
    objective,
    x0,
    method="SLSQP",
    bounds=bounds,
    constraints=cons,
    options={
        "disp": True,
        "eps": 0.1,
        "maxiter": 20
    }
)

print("\nOptimal x:", result.x)
print("Optimal objective:", result.fun)
print("Message:", result.message)


#Optimizing
#result = spo.minimize(lambda x: RunAPDL(SWcoor, var, Misc), var, options={"disp": True})



