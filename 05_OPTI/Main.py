## MAIN BEAM ELEMENT DOCUMENT ##

# Import packages
import sys
print(sys.version, flush=True)
from ansys.mapdl.core import launch_mapdl

# Import Functions
import optimization
from MyAPDLCall import RunAPDL

# Variables. Choose which ones to include by setting "active": True or False.
# Include Bounds
opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "multi_size_columns": True,   # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": True,    # Whether mast segments braces uses different dimensions (True) or not (False)
    "brace_split": True,          # Whether braces are split between horiontal and cross (True) or not (False)
}

# Defining variables with bounds and active status
var = {
    "rad": {"value": 202.07, "bounds": (150.0, 300.0), "active": True}, # Radius Structure [mm]
}
if opti_settings["multi_size_columns"]:
    var.update({
        **{f"d0_{i}": {"value": 76.1, "bounds": (50.0, 100), "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Diameter  [mm]
        **{f"t0_{i}": {"value": 3.0,  "bounds": (1.0, 7.0),  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Thickness [mm]
    })
else:
    var.update({
        "d0": {"value": 76.1, "bounds": (50.0, 100), "active": True},       # Column Diameter  [mm]
        "t0": {"value": 3.0,  "bounds": (1.0, 7.0),  "active": True},       # Column Thickness [mm]
    })
if opti_settings["multi_size_braces"]:
    var.update({
        **{f"d1_{i}": {"value": 26.9, "bounds": (10.0, 40.0), "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},     # Brace Diameter   [mm]
        **{f"t1_{i}": {"value": 2.3,  "bounds": (0.1, 4.5),  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Brace Thickness  [mm]
    })
else:
    var.update({
        "d1": {"value": 26.9, "bounds": (10.0, 40.0), "active": True},     # Brace Diameter   [mm]
        "t1": {"value": 2.3,  "bounds": (0.1, 4.5),  "active": True},       # Brace Thickness  [mm]
    })

# Static variables
Misc = {
    "esize": 3,                     # Element Size [mm]
    "Hor_Force": 502.52,            # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)   [N]
    "f_y": 690 ,                    # Column Yield Strength [MPa]
    "f_y_brace": 355,               # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
    "eps_geom": 0.1,    # Minimum thickness specification for geometry updates [mm]
    "SW_filename": "LWC_L1_LINES.IGS",    # Filename for IGS File
    "save_folder": "Optimization_Logs" # Save Folder
}

# Solver Settings
Solver_Settings = {
    "acc": 1e-3,       # Maximum objective function tolerance
    "maxiter": 10,     # Maximum iterations
}

# Launch MAPDL
mapdl = launch_mapdl(
    run_location="Ansout",
    log_apdl="apdl_log",
    override=True,
    nproc=2,
    additional_switches="-p ansys -smp",
)

# Ensures that MAPDL closes if something chrashes
try:
    f1 = RunAPDL(mapdl,var,Misc,opti_settings)
    result, txt_path, csv_path = optimization.run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings)
finally:
    mapdl.exit()

# Print Information
print("\n" + "="*60)
print("OPTIMIZATION RESULTS")
print("="*60)
print(f"\nOptimal Mass: {result.fun:.3f} kg")
print(f"\nOptimal Design Variables:")
# Get active variable names in order
active_names = [name for name, data in var.items() if data.get("active", True)]
for name, value in zip(active_names, result.x):
    print(f"  {name:8s}: {value:10.4f}")
print(f"\nOptimization Message: {result.message}")
print(f"\nLog files saved to:")
print(f"  - TXT: {txt_path}")
print(f"  - CSV: {csv_path}")
print("="*60)