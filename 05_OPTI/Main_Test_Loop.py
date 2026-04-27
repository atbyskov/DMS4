## MAIN BEAM ELEMENT DOCUMENT ##
# Main.py               -> Runs the entire script
# SW_Import.py          -> Reads IGS File and converts it to Keypoints and Lines for APDL
# MyAPDLCall.py         -> Inputs Swcoor, var and Misc and Runs the APDL Script
#   APDL_Input.py       -> Writes the .txt input file for APDL
# Post_Process.py       -> Reads internal forces from analysis and outputs Utilization Ratios


# Import packages
import os
import time 
from ansys.mapdl.core import launch_mapdl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import Functions
import SW_Import as SW
from Post_Process import PostProcessor
from MyAPDLCall import RunAPDL


rad_values = []
eigen_values = []
mass_values = []

mapdl = launch_mapdl(
        run_location="Ansout", 
        log_apdl="apdl_logging",
        override=True,
        nproc=6,
        additional_switches="-p ansys -smp"
    )

for rad in range(150,350,1):
    # Collect variables
    var = {
        "d0": {"value": 76.1,    "active": True},        # Column Diameter  [mm]
        "t0": {"value": 3.0,     "active": True},        # Column Thickness [mm]
        "d1": {"value": 26.9,    "active": True},        # Brace Diameter   [mm]
        "t1": {"value": 2.3,     "active": True},        # Brace Thickness  [mm]
        "rad": {"value": rad, "active": True},       # Radius Structure [mm]
    }

    # Create Misc as dict
    Misc = {
        "esize": 3,                     # Element Size [mm]
        "Hor_Force": 502.52,            # Horizontal Force (P_Load_z) [N]
        "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)   [N]
        "f_y": 690 ,                    # Column Yield Strength [MPa]
        "f_y_brace": 355,               # Brace Yield Strength [MPa]
        "E_mod": 200*1E3,               # Youngs Modulus [MPa]
        "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
        "SW_filename": "LWC_LC1.IGS"    # Filename for IGS File
    }

    #C:\Program Files\ANSYS Inc\v251\ansys\bin\winx64
    
    # Run Environment
    mass = RunAPDL(mapdl, var, Misc) # Runs APDL and returns MASS


    print(f" -> Mass of Assembly: {mass:.2f} kg")

    with open("Ansout/Eigenvalue1.txt") as file:
        eigenvalues = [float(line.strip()) for line in file if line.strip()]
    # Retrieve first positive eigenvalue
    alpha_crit = next(v for v in eigenvalues if v > 0) 

    rad_values.append(rad)
    eigen_values.append(alpha_crit)
    mass_values.append(mass)

    print(f"rad = {rad:6.1f} mm | eigenvalue = {alpha_crit:.6f}")

    
    with open("rad_results.txt", "w") as f:
        f.write("rad_mm\teigenvalue\tmass_kg\n")
        for r, e, m in zip(rad_values, eigen_values, mass_values):
            f.write(f"{r:.2f}\t{e:.6f}\t{m:.3f}\n")



    # Pack specifically for compability
    active = [(name,data["value"]) for name,data in var.items() if data.get("active", True)]
    names = [name for name,_ in active]
    x0 = np.array([var[name]["value"] for name, _ in active], dtype=float)
    var_dict = dict(zip(names, x0))


    utils = PostProcessor(var_dict, Misc)
    Util_list = utils.Util_list()  

mapdl.exit()


plt.figure()
plt.plot(rad_values, eigen_values, marker="o")
plt.xlabel("Radius [mm]")
plt.ylabel("Eigenvalue")
plt.title("Eigenvalue vs Radius")
plt.grid(True)
plt.tight_layout()
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.scatter(rad_values, mass_values, eigen_values, c=eigen_values, cmap="viridis")
ax.plot(rad_values, mass_values, eigen_values)

ax.set_xlabel("Radius [mm]")
ax.set_ylabel("Mass [kg]")
ax.set_zlabel("Eigenvalue")
ax.set_title("Eigenvalue vs Radius & Mass")

plt.tight_layout()
plt.show()

