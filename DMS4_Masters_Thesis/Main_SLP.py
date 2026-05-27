"""
MAIN BEAM ELEMENT OPTIMIZATION PROGRAM (SLP)
--------------------------------------------

This script serves as the primary driver for structural optimization of a
beam-based mast model using Sequential Linear Programming (SLP) coupled
with ANSYS MAPDL for finite element analysis.

Inputs:
-------
- opti_settings    : Structural configuration parameters
- var              : Design variables (values, bounds, activity flags)
- Misc             : Constants (loads, materials, geometry file, etc.)
- Solver_Settings  : SLP solver configuration

Outputs:
--------
- result           : Optimization result object
- txt_path         : Path to detailed log file
- csv_path         : Path to evaluation history file

External Dependencies:
----------------------
- ansys.mapdl.core        : MAPDL interface (PyMAPDL)
- RunAPDL                 : FE model execution
- optimization_SLP        : SLP optimization engine
- Supporting modules:
    • APDL_Input
    • Post_Process
    • ConstraintAggregate
    • ACS (optional)
    • opt_logger

"""

# Import packages
import sys
print(sys.version, flush=True)
import time
from ansys.mapdl.core import launch_mapdl

# Import Functions
import SLP.optimization_SLP as optimization_SLP
from MyAPDLCall import RunAPDL

start = time.time()
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
column_diameter_bounds = (48.3, 114.3) # Column Diameter Bounds [mm]
column_thickness_bounds = (2.5, 5.0)   # Column Thickness Bounds [mm]
brace_diameter_bounds = (25, 60.0)   # Brace Diameter Bounds [mm]
brace_thickness_bounds = (2.0, 6.0)    # Brace Thickness Bounds [mm]

# Defining variables with bounds and active status
var = {
    "rad": {"value": 202.07, "bounds": (150.0, 350.0), "active": True}, # Radius Structure [mm]
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
    "f_y_brace": 355,               # Brace Yield Strength [MPa]
    "E_mod": 200*1E3,               # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,           # Vertical Force COG (P_COG_y) [N]
    "eps_geom": 0.1,    # Minimum thickness specification for geometry updates [mm]
    "SW_filename": "LWC_L1_LINES.IGS",    # Filename for IGS File
    "save_folder": "Optimization_Logs" # Save Folder
}

# Solver Settings -- match fminslp.m defaults (slpoptions, lines 794-849).
# MATLAB-equivalent option names are noted in the comments.
Solver_Settings = {
    # ---- Aggregation (NOT a solver parameter; consumed by ConstraintAggregate)
    "Aggregate": None,            # None | "P-norm" | "P-norm-mean" |
    "p_value": 8,                 # P-norm exponent
    "relaxation": 0,              # aggregation relaxation
    # ---- Finite differences (MATLAB default: sqrt(eps); too small for FEA)
    "finite_diff_rel_step": 1e-3, # absolute step = rel_step * max(1, |x_i|)                                   #IMPORTANT TO TUNE FOR FEA

    # ---- Penalty
    "penalty_weight": 15,     # InfeasibilityPenalization R (fixed, never grows)

    # ---- Convergence (MATLAB names)
    "acc": 1e-3,                  # FunctionTolerance       -> |delta f_merit|
    "xtol": 1e-6,                 # StepTolerance           -> ||dx||
    "gtol": 1e-6,                 # OptimalityTolerance     -> ||grad f_merit|| over x part
    "maxiter": 5000,              # MaxIterations
    "max_fun_evals": 5000,        # MaxFunctionEvaluations  -- counts LP trial points, mirrors MATLAB nFeval
    "objective_limit": -1e20,     # ObjectiveLimit          -- early stop if f_merit <= this

    # ---- Move limits
    "move_limit": 0.10,           # MoveLimit
    "move_limit_expand": 1.1,     # MoveLimitExpand
    "move_limit_shrink": 0.5,     # MoveLimitReduce
    "use_acs": False,              # Toggle ACS to on (True) or off (False)
    # ---- Filter
    "max_infeasibility": float("inf"),  # MaxInfeasibility  (inf => filter starts unconstrained)
}

# Launch MAPDL
mapdl = launch_mapdl(
    run_location="Ansout",
    log_apdl="apdl_log",
    override=True,
    nproc=10,
    additional_switches="-p ansys -smp",
)

# Ensures that MAPDL closes if something chrashes
try:
    # f1 = RunAPDL(mapdl,var,Misc,opti_settings) # Test call
    result, txt_path, csv_path = optimization_SLP.run_optimization(mapdl, opti_settings, var, Misc, Solver_Settings, method="slp_mvp")
finally:
    mapdl.exit()

end = time.time()
print(f"Optimization completed in {end - start:.2f} seconds.")

# Print Information
print("\nOptimal x:", result.x)
print("Optimal objective:", result.fun)
print("Message:", result.message)
print("TXT log file:", txt_path)
print("Objective CSV:", csv_path)
