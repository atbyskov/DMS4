## MAIN BEAM ELEMENT DOCUMENT ##

# Import packages
import sys
print(sys.version, flush=True)
import time
from ansys.mapdl.core import launch_mapdl

# Import Functions
import Self_Written_Optimization.optimization_self_written as optimization_self_written
from MyAPDLCall import RunAPDL


# Variables. Choose which ones to include by setting "active": True or False.
# Include Bounds
opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "segment_mass_limit": 23,     # Limits for segment masses [kg]
    "multi_size_columns": False,   # Whether mast segments columns uses different dimensions (True) or not (False)
    "multi_size_braces": False,    # Whether mast segments braces uses different dimensions (True) or not (False)
    "brace_split": False,          # Whether braces are split between horiontal and cross (True) or not (False)
}

# Initial Guess
column_diameter = 76.1 # Column Diameter [mm]
column_thickness = 3.0 # Column Thickness [mm]
brace_diameter = 26.9 # Brace Diameter [mm]
brace_thickness = 2.3 # Brace Thickness [mm]

# Bounds
column_diameter_bounds = (40.0, 100) # Column Diameter Bounds [mm]
column_thickness_bounds = (1.0, 7.0)   # Column Thickness Bounds [mm]
brace_diameter_bounds = (10.0, 100.0)   # Brace Diameter Bounds [mm]
brace_thickness_bounds = (1.0, 4.5)    # Brace Thickness Bounds [mm]

# Defining variables with bounds and active status
var = {
    "rad": {"value": 273.7506, "bounds": (150.0, 350.0), "active": True}, # Radius Structure [mm]
}
if opti_settings["multi_size_columns"]:
    var.update({
        **{f"d0_{i}": {"value": column_diameter, "bounds": column_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Diameter  [mm]
        **{f"t0_{i}": {"value": column_thickness,  "bounds": column_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Column Thickness [mm]
    })
else:
    var.update({
        "d0": {"value": column_diameter, "bounds": column_diameter_bounds, "active": True},       # Column Diameter  [mm]
        "t0": {"value": column_thickness,  "bounds": column_thickness_bounds,  "active": True},       # Column Thickness [mm]
    })
if opti_settings["multi_size_braces"]:
    if opti_settings["brace_split"]:
        var.update({
            **{f"d1_h_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},     # Horizontal Brace Diameter   [mm]
            **{f"t1_h_{i}": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Horizontal Brace Thickness  [mm]
            **{f"d1_c_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},     # Cross Brace Diameter       [mm]
            **{f"t1_c_{i}": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Cross Brace Thickness      [mm]
        })
    else:
        var.update({
            **{f"d1_{i}": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},     # Brace Diameter   [mm]
            **{f"t1_{i}": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True} for i in range(1, opti_settings["n_mast_segments"]+1)},       # Brace Thickness  [mm]
        })
else:
    if opti_settings["brace_split"]:
        var.update({
            "d1_h": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},     # Horizontal Brace Diameter   [mm]
            "t1_h": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True},       # Horizontal Brace Thickness  [mm]
            "d1_c": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},     # Cross Brace Diameter       [mm]
            "t1_c": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True},       # Cross Brace Thickness      [mm]
        })
    else:
        var.update({
            "d1": {"value": brace_diameter, "bounds": brace_diameter_bounds, "active": True},     # Brace Diameter   [mm]
            "t1": {"value": brace_thickness,  "bounds": brace_thickness_bounds,  "active": True},       # Brace Thickness  [mm]
        })

# Static variables
Misc = {
    "esize": 3,                     # Element Size [mm]
    "Hor_Force": 502.52,            # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,         # Vertical Force (P_Load_y)   [N]
    "f_y": 700 ,                    # Column Yield Strength [MPa]
    "f_y_brace": 235,               # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
    "eps_geom": 0.1,    # Minimum thickness specification for geometry updates [mm]
    "SW_filename": "LWC_L1_LINES.IGS",    # Filename for IGS File
    "save_folder": "Optimization_Logs" # Save Folder
}

# Solver Settings
Solver_Settings = {
    "acc": 1e-3,                 # Maximum objective function tolerance
    "maxiter": 800,               # Maximum iterations                                                         #IMPORTANT TO TUNE FOR FEA
    "Aggregate": None,           # None, "P-norm", "P-norm-mean", "KS", "KS_shift"  (Write exacly)
    "p_value": 8,                # Value for "P-norm" and "P-norm-mean"
    "rho_value": 100,            # rho value used in KS
    "relaxation": 0,             # Relaxation parameter used in aggregation
    "finite_diff_rel_step": 1e-3, # PySLSQP-like absolute FD step: rel_step * max(1, abs(x))                              #IMPORTANT TO TUNE FOR FEA
    "algorithm": "merit",         # 'merit' (quadratic merit) or 'al' (augmented Lagrangian)
    "penalty_weight": 1e2,        # Initial L1 merit/slack penalty
    "penalty_increase": 2.0,      # Increase penalty when slacks/infeasibility persist
    "penalty_max": 1e6,           # Upper cap for adaptive penalty
    "move_limit": 0.10,           # Initial/max scaled trust-region size
    "move_limit_min": 1e-5,       # Stop after failed globalization below this scaled size
    "move_limit_expand": 1.4,     # Expand trust region after good predicted-vs-actual agreement
    "move_limit_shrink": 0.5,     # Shrink trust region after rejection or poor agreement
    "backtrack_max": 8,           # Trial alphas along each LP step
    "backtrack_shrink": 0.5,      # alpha <- alpha * shrink during SLP line search
    "feasibility_tol": 1e-5,      # Required max nonlinear constraint violation for convergence
    "slack_tol": 1e-8,            # Required LP slack size for convergence
    "xtol": 1e-4,                 # Physical design step tolerance
    "gtol": 1e-6,                 # Scaled gradient infinity-norm tolerance
    "sufficient_decrease": 1e-4,  # Merit decrease fraction of predicted decrease
    "feasibility_reduction": 1e-3,# Accept infeasible steps that reduce violation by this fraction          #IMPORTANT TO TUNE FOR FEA
}

# Launch MAPDL
mapdl = launch_mapdl(
    run_location="Ansout",
    log_apdl="apdl_log",
    override=True,
    nproc=8,
    additional_switches="-p ansys -smp",
)

# Ensures that MAPDL closes if something chrashes
try:
    # f1 = RunAPDL(mapdl,var,Misc,opti_settings) # Test call
    result, txt_path, csv_path = optimization_self_written.run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings, method="slp_mvp")
finally:
    mapdl.exit()

# Print Information
#print("\nOptimal x:", result.x)
#print("Optimal objective:", result.fun)
#print("Message:", result.message)
#print("TXT log file:", txt_path)
#print("Objective CSV:", csv_path)
