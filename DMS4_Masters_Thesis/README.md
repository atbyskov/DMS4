# DMS4 Master's Thesis — Mast Optimization

Python framework that drives ANSYS MAPDL to size a parametric mast structure
(columns, braces, radius) against Eurocode 3 utilization checks.

## Requirements

- Python 3
- ANSYS MAPDL installed and on PATH (the scripts call `launch_mapdl`)
- Python packages: `ansys-mapdl-core`, `numpy`, `pandas`, `scipy`

## How to run

From inside this folder:

```powershell
python Main.py          # optimization with PySLSQP
python Main_SLP.py      # optimization with the custom SLP solver
python BeamModel.py     # single FE evaluation (no optimization)
```

All three scripts launch MAPDL inside the `Ansout/` working folder and write
logs to `Optimization_Logs/`.

> **Note:** If you interrupt or close down an optimization run (e.g. `Main.py`
> or `Main_SLP.py`) before it finishes, the MAPDL process is not always shut
> down cleanly. You will need to manually close the leftover MAPDL window /
> process in Task Manager before starting a new run, otherwise the
> next launch may fail or hang on the `Ansout/` working directory.

Edit the `var`, `Misc`, `opti_settings`, and `Solver_Settings` dictionaries at
the top of the chosen script to change design variables, loads, material, or
solver behavior.

## What's in here

**Entry points**
- `Main.py`         — PySLSQP optimization driver.
- `Main_SLP.py`     — custom Sequential Linear Programming driver.
- `BeamModel.py`    — one-shot analysis of a single design point.

**Core modules**
- `optimization.py` — orchestrates the PySLSQP loop (variables, constraints, aggregation, logging).
- `MyAPDLCall.py`   — bridge between Python and MAPDL; runs the FE analysis and returns masses.
- `APDL_Input.py`   — builds the full APDL command script (mesh, sections, BCs, loads, buckling, nonlinear static).
- `Post_Process.py` — Eurocode 3 utilization checks (local buckling, NF, shear, torsion, BR, interaction, braces).
- `SW_Import.py`    — parses IGES line entities from SolidWorks.
- `aggregate.py`    — P-norm / P-norm-mean constraint aggregation.
- `acs.py`          — Adaptive Constraint Scaling.
- `opt_logger.py`   — writes per-iteration TXT log and per-evaluation CSV history.
- `search_space.py` — sampling utilities (grid, random, LHS, Sobol) for initial points.

**SLP solver (`SLP/`)**
- `SLP/optimization_SLP.py` — SLP-side orchestrator.
- `SLP/SLP.py`              — merit-function SLP solver with LP subproblem, Fletcher–Leyffer–Toint filter, and adaptive move limits.

**Data and output folders**
- `IGS/`                    — mast wireframe geometry (`LWC_L1_LINES.IGS`).
- `Ansout/`                 — MAPDL working directory (created/overwritten at runtime).
- `Optimization_Logs/`      — TXT logs and CSV evaluation histories.

## Output

Each run produces, inside `Optimization_Logs/`:
- `optimization_log_<timestamp>.txt`  — header + per-iteration utilization table.
- `objective_history_<timestamp>.csv` — one row per function evaluation
  (objective, segment masses, aggregated violation, `g_max`).
