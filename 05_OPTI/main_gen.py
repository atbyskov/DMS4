## GENETIC ALGORITHM OPTIMIZATION (main_gen.py) ##
# -----------------------------------------------------------------------------
# Discrete, catalogue-based Genetic Algorithm for the lattice mast.
#
# Why a GA here:
#   The gradient-based SLP/SLSQP paths treat the section sizes as *continuous*
#   variables and then rely on the Class-2 constraint to keep them realistic.
#   In reality only a finite set of production profiles can be bought:
#       * Columns -> Strenx 700 (S700) CHS profiles
#       * Braces  -> S235 CHS profiles
#   Those allowed (d, t) points are exactly the ones drawn in
#   Figures/Class2.py.  A Genetic Algorithm is the natural fit for picking one
#   catalogue entry per member group (a *combinatorial* problem) while still
#   handling the one continuous variable we keep (the structure radius `rad`).
#
# Package:
#   pymoo (https://pymoo.org) -- the de-facto standard evolutionary-optimisation
#   library in Python.  We use its `MixedVariableGA` so the radius stays a real
#   variable while every section is an Integer index into its production
#   catalogue.  Constraints are handled by pymoo's feasibility-first selection.
#
# It reuses the exact same FEA / post-processing pipeline as Main.py:
#       RunAPDL  ->  PostProcessor  ->  utilisation / eigenvalue / mass checks
#
# Run from the 05_OPTI directory (relative imports + Ansout/ paths):
#       cd 05_OPTI
#       python main_gen.py
# -----------------------------------------------------------------------------

# ---- Standard / scientific imports ------------------------------------------
import sys
print(sys.version, flush=True)
import os
import time
import types
import numpy as np

from ansys.mapdl.core import launch_mapdl

# ---- pymoo ------------------------------------------------------------------
try:
    from pymoo.core.problem import ElementwiseProblem
    from pymoo.core.variable import Real, Integer
    from pymoo.core.mixed import MixedVariableGA, MixedVariableSampling
    from pymoo.core.callback import Callback
    from pymoo.optimize import minimize
except ImportError as exc:  # pragma: no cover - guidance only
    raise SystemExit(
        "pymoo is required for main_gen.py. Install it with:\n"
        "    python -m pip install pymoo\n"
        f"(import error: {exc})"
    )

# ---- Project imports (same pipeline as Main.py) -----------------------------
from MyAPDLCall import RunAPDL
from Post_Process import PostProcessor
from opt_logger import OptimizationLogger


# =============================================================================
# 1. PRODUCTION PROFILE CATALOGUES  (copied from Figures/Class2.py)
# =============================================================================
# NOTE: the data is duplicated here on purpose. Importing Figures/Class2.py
# would execute its plotting/`plt.show()` side effects, so we keep the raw
# (d, t) tables local and self-contained.

# High-Strength Steel Columns (S700) -- (outer diameter [mm], wall thickness [mm])
S700_DATA = [
    (48.3, 2.5), (48.3, 3.0), (60.3, 2.5), (60.3, 3.0), (60.3, 4.0),
    (76.1, 3.0), (76.1, 4.0), (76.1, 5.0), (88.9, 3.0), (88.9, 4.0),
    (88.9, 5.0), (101.6, 3.0), (101.6, 4.0), (101.6, 5.0), (108.0, 3.0),
    (108.0, 4.0), (108.0, 5.0), (114.3, 3.0), (114.3, 4.0), (114.3, 5.0),
    (121.0, 4.0), (121.0, 5.0), (127.0, 3.0), (127.0, 4.0), (127.0, 5.0),
    (139.7, 4.0), (139.7, 5.0), (139.7, 6.0), (139.7, 8.0), (139.7, 10.0),
    (168.3, 4.0), (168.3, 5.0), (168.3, 6.0), (168.3, 8.0), (168.3, 10.0),
    (219.1, 5.0), (219.1, 6.0), (219.1, 8.0), (219.1, 10.0), (273.0, 5.0),
    (273.0, 6.0), (273.0, 8.0), (273.0, 10.0), (323.9, 6.0), (323.9, 8.0),
    (323.9, 10.0),
]

# Mild Steel Braces (S235) - defined raw (diameter, [thickness list]) then flattened
S235_RAW = [
    (6, [1, 1.5, 2]),
    (8, [1, 1.5, 2]),
    (10, [1, 1.5, 2, 2.5]),
    (12, [1, 1.5, 2, 2.5]),
    (13, [1, 1.5, 2, 2.5]),
    (14, [1, 1.5, 2, 2.5, 3, 3.5]),
    (15, [1, 1.5, 2, 2.5, 3]),
    (16, [1, 1.5, 2, 2.5, 3, 3.5, 4]),
    (17, [1, 1.5, 2, 2.5]),
    (18, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]),
    (20, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]),
    (21, [1, 1.5, 2, 2.5, 3, 4]),
    (22, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]),
    (23, [1, 1.5, 2]),
    (24, [1, 1.5, 2, 3, 4]),
    (25, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7]),
    (28, [1, 1.5, 2, 2.5, 3, 4, 5, 6]),
    (30, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 7.5]),
    (32, [1.5, 2, 3, 3.5, 4, 5, 6]),
    (34, [1.5, 2, 4]),
    (35, [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 7.5, 8, 9]),
    (36, [1.5, 2, 3, 5, 8]),
    (38, [1.5, 2, 3, 4, 4.5, 5, 6, 7, 8]),
    (40, [1.5, 2, 2.5, 3, 4, 4.5, 5, 5.5, 6, 7, 7.5, 8, 10]),
    (42, [1.5, 2, 3, 4, 5, 7]),
    (45, [1.5, 2, 2.5, 3, 4, 4.5, 5, 6, 7, 10]),
    (50, [1.5, 2, 3, 4, 5, 6, 7, 8, 10]),
    (52, [5]),
    (55, [1.5, 2, 3, 4, 5, 6, 8, 10]),
    (60, [1.5, 2, 3, 4, 5, 6, 8, 9, 10]),
    (62, [6]),
    (65, [4, 6, 10]),
    (70, [2, 3, 4, 5, 6, 8, 10]),
    (80, [2, 3, 4, 5, 6, 8, 10, 12]),
    (90, [5]),
    (100, [3, 4, 5, 6]),
]
S235_DATA = [(float(d), float(t)) for d, t_list in S235_RAW for t in t_list]


def build_catalogue(data, d_bounds, t_bounds, f_y, label):
    """Filter a raw (d, t) production table to the points that are usable.

    A point is kept iff it is
        * inside the project design bounds for diameter and thickness, and
        * Class 2 (or better) per EN 1993-1-1 Table 5.2:  d / t <= 70 * 235 / f_y.

    These are exactly the green/compliant points inside the dashed design-bound
    box in Figures/Class2.py.
    """
    class2_limit = 70.0 * (235.0 / f_y)
    dlo, dhi = d_bounds
    tlo, thi = t_bounds

    catalogue = [
        (float(d), float(t))
        for (d, t) in data
        if dlo - 1e-9 <= d <= dhi + 1e-9
        and tlo - 1e-9 <= t <= thi + 1e-9
        and (d / t) <= class2_limit + 1e-9
    ]
    # Sort by cross-sectional area (light -> heavy) so neighbouring catalogue
    # indices are roughly monotonic in mass; helps the GA's mutation locality.
    catalogue.sort(key=lambda dt: dt[0] ** 2 - (dt[0] - 2 * dt[1]) ** 2)

    if not catalogue:
        raise ValueError(
            f"No {label} catalogue points satisfy the bounds + Class-2 filter "
            f"(d in {d_bounds}, t in {t_bounds}, d/t <= {class2_limit:.3f})."
        )
    print(
        f"{label} catalogue: {len(catalogue)} Class-2 production profiles in bounds "
        f"(d/t <= {class2_limit:.2f}).",
        flush=True,
    )
    return catalogue


# =============================================================================
# 2. PROBLEM SETTINGS  (mirrors Main.py)
# =============================================================================
opti_settings = {
    "n_mast_segments": 5,         # Number of mast segments
    "mast_segment_height": 810,   # Height of each mast segment [mm]
    "segment_mass_limit": 23,     # Limits for segment masses [kg]
    "multi_size_columns": True,   # Columns may differ per segment (True) or not (False)
    "multi_size_braces": True,    # Braces may differ per segment (True) or not (False)
    "brace_split": True,          # Horizontal/Cross braces sized separately (True) or not (False)
}

# Static variables (identical to Main.py: S700 columns, S235 braces)
Misc = {
    "esize": 3,                           # Element Size [mm]
    "Hor_Force": 502.52,                  # Horizontal Force (P_Load_z) [N]
    "Ver_Force": -25.13E+3,               # Vertical Force (P_Load_y)   [N]
    "f_y": 700,                           # Column Yield Strength [MPa]  (S700)
    "f_y_brace": 235,                     # Brace Yield Strength [MPa]   (S235)
    "E_mod": 200 * 1E3,                   # Youngs Modulus [MPa]
    "W_Force": -3.751E+3,                 # Vertical Force COG (P_COG_y) [N]
    "SW_filename": "LWC_L1_LINES.IGS",    # Filename for IGS File
    "save_folder": "Optimization_Logs",   # Save Folder
}

# Bounds -- identical to Main.py
column_diameter_bounds = (48.3, 114.3)   # Column Diameter Bounds [mm]
column_thickness_bounds = (2.5, 5.0)     # Column Thickness Bounds [mm]
brace_diameter_bounds = (10.0, 50.0)     # Brace Diameter Bounds [mm]
brace_thickness_bounds = (1.0, 4.0)      # Brace Thickness Bounds [mm]
rad_bounds = (150.0, 350.0)              # Radius Structure Bounds [mm]

# Initial guess (used only to seed one individual + log x0) -- same as Main.py
init_rad = 202.07
init_col_dt = (76.1, 3.0)
init_brace_dt = (26.9, 2.3)

# GA / solver settings
GA_Settings = {
    "pop_size": 24,      # Number of individuals per generation
    "n_gen": 20,         # Number of generations  (FEA evals ~ pop_size * n_gen)
    "seed": 1,           # RNG seed for reproducibility
    "eigenvalue_min": 4.0,   # Required first positive buckling eigenvalue (a_cr >= 4)
    "fail_penalty": 1e9,     # Objective assigned to FEA failures (kept infeasible)
}


# =============================================================================
# 3. CATALOGUES + GENOME LAYOUT
# =============================================================================
COL_CATALOGUE = build_catalogue(
    S700_DATA, column_diameter_bounds, column_thickness_bounds,
    Misc["f_y"], "Column (S700)",
)
BRACE_CATALOGUE = build_catalogue(
    S235_DATA, brace_diameter_bounds, brace_thickness_bounds,
    Misc["f_y_brace"], "Brace (S235)",
)

N_SEG = int(opti_settings["n_mast_segments"])


def _nearest_index(catalogue, target_dt):
    """Index of the catalogue entry closest (Euclidean in d,t) to target_dt."""
    td, tt = target_dt
    dists = [(d - td) ** 2 + (t - tt) ** 2 for (d, t) in catalogue]
    return int(np.argmin(dists))


# Build the list of genome "genes" (pymoo variable names) for the active layout.
# Columns -------------------------------------------------------------------
if opti_settings["multi_size_columns"]:
    COL_GENES = [f"col_{i}" for i in range(1, N_SEG + 1)]
else:
    COL_GENES = ["col"]

# Braces --------------------------------------------------------------------
if opti_settings["multi_size_braces"]:
    if opti_settings["brace_split"]:
        BRACE_GENES = [f"bh_{i}" for i in range(1, N_SEG + 1)] + \
                      [f"bc_{i}" for i in range(1, N_SEG + 1)]
    else:
        BRACE_GENES = [f"b_{i}" for i in range(1, N_SEG + 1)]
else:
    if opti_settings["brace_split"]:
        BRACE_GENES = ["bh", "bc"]
    else:
        BRACE_GENES = ["b"]

GENE_NAMES = ["rad"] + COL_GENES + BRACE_GENES

# pymoo mixed-variable definition: rad is Real, every section is an Integer
# index into its production catalogue.
PYMOO_VARS = {"rad": Real(bounds=rad_bounds)}
for g in COL_GENES:
    PYMOO_VARS[g] = Integer(bounds=(0, len(COL_CATALOGUE) - 1))
for g in BRACE_GENES:
    PYMOO_VARS[g] = Integer(bounds=(0, len(BRACE_CATALOGUE) - 1))


def decode_to_var_dict(genome):
    """Map a pymoo genome dict -> nested var dict consumed by RunAPDL/PostProcessor.

    Produces exactly the keys APDL_Input.InputFun / PostProcessor expect for the
    current (multi_size_*, brace_split) layout, e.g. d0_1, t0_1, d1_h_1, ...
    Section values come from the catalogues so they are always real, buyable,
    Class-2 profiles.
    """
    var = {"rad": {"value": float(genome["rad"]), "active": True}}

    def put(name, value):
        var[name] = {"value": float(value), "active": True}

    # --- Columns (S700) ----------------------------------------------------
    if opti_settings["multi_size_columns"]:
        for i in range(1, N_SEG + 1):
            d, t = COL_CATALOGUE[int(genome[f"col_{i}"])]
            put(f"d0_{i}", d)
            put(f"t0_{i}", t)
    else:
        d, t = COL_CATALOGUE[int(genome["col"])]
        put("d0", d)
        put("t0", t)

    # --- Braces (S235) -----------------------------------------------------
    if opti_settings["multi_size_braces"]:
        if opti_settings["brace_split"]:
            for i in range(1, N_SEG + 1):
                d, t = BRACE_CATALOGUE[int(genome[f"bh_{i}"])]
                put(f"d1_h_{i}", d)
                put(f"t1_h_{i}", t)
            for i in range(1, N_SEG + 1):
                d, t = BRACE_CATALOGUE[int(genome[f"bc_{i}"])]
                put(f"d1_c_{i}", d)
                put(f"t1_c_{i}", t)
        else:
            for i in range(1, N_SEG + 1):
                d, t = BRACE_CATALOGUE[int(genome[f"b_{i}"])]
                put(f"d1_{i}", d)
                put(f"t1_{i}", t)
    else:
        if opti_settings["brace_split"]:
            d, t = BRACE_CATALOGUE[int(genome["bh"])]
            put("d1_h", d)
            put("t1_h", t)
            d, t = BRACE_CATALOGUE[int(genome["bc"])]
            put("d1_c", d)
            put("t1_c", t)
        else:
            d, t = BRACE_CATALOGUE[int(genome["b"])]
            put("d1", d)
            put("t1", t)

    return var


# Ordered list of *engineering* variable names (for logging x in a stable order)
def engineering_var_names():
    names = ["rad"]
    if opti_settings["multi_size_columns"]:
        names += [f"d0_{i}" for i in range(1, N_SEG + 1)]
        names += [f"t0_{i}" for i in range(1, N_SEG + 1)]
    else:
        names += ["d0", "t0"]
    if opti_settings["multi_size_braces"]:
        if opti_settings["brace_split"]:
            names += [f"d1_h_{i}" for i in range(1, N_SEG + 1)]
            names += [f"t1_h_{i}" for i in range(1, N_SEG + 1)]
            names += [f"d1_c_{i}" for i in range(1, N_SEG + 1)]
            names += [f"t1_c_{i}" for i in range(1, N_SEG + 1)]
        else:
            names += [f"d1_{i}" for i in range(1, N_SEG + 1)]
            names += [f"t1_{i}" for i in range(1, N_SEG + 1)]
    else:
        if opti_settings["brace_split"]:
            names += ["d1_h", "t1_h", "d1_c", "t1_c"]
        else:
            names += ["d1", "t1"]
    return names


ENG_NAMES = engineering_var_names()


def var_dict_to_xarray(var):
    """Flatten a nested var dict into a numeric vector in ENG_NAMES order."""
    return np.array([var[name]["value"] for name in ENG_NAMES], dtype=float)


def initial_genome():
    """Seed genome from Main.py's initial guess (snapped to nearest catalogue)."""
    g = {"rad": float(init_rad)}
    col_idx = _nearest_index(COL_CATALOGUE, init_col_dt)
    brc_idx = _nearest_index(BRACE_CATALOGUE, init_brace_dt)
    for gene in COL_GENES:
        g[gene] = col_idx
    for gene in BRACE_GENES:
        g[gene] = brc_idx
    return g


# =============================================================================
# 4. THE pymoo PROBLEM  (wraps the FEA pipeline)
# =============================================================================
class MastGAProblem(ElementwiseProblem):
    """Single-objective (mass) constrained GA problem.

    Objective : total assembly mass [kg]  (minimise)
    Constraints (pymoo form g(x) <= 0 feasible):
        g_util_struct : max EC3 member utilisation - 1            (LB/NF/S/T/BR/IN)
        g_util_brace  : max brace-step utilisation - 1            (stress + defl)
        g_eig         : eigenvalue_min - a_cr                     (a_cr >= 4)
        g_mass_i      : seg_mass_i - segment_mass_limit           (one per segment)
    Class-2 is guaranteed by the catalogue filter, so it is not re-imposed here.
    """

    def __init__(self, mapdl, logger):
        self.mapdl = mapdl
        self.logger = logger
        self.eig_min = float(GA_Settings["eigenvalue_min"])
        self.mass_limit = float(opti_settings["segment_mass_limit"])
        self.n_g = 3 + N_SEG  # 2 utilisation + 1 eigenvalue + per-segment mass

        # Best feasible individual seen so far (for logging / final report)
        self.best = None  # dict: f, var, segment_masses, util_report, max_util

        super().__init__(vars=PYMOO_VARS, n_obj=1, n_ieq_constr=self.n_g)

    # ---- helpers ----------------------------------------------------------
    @staticmethod
    def _arr(v):
        return np.asarray(v, dtype=float).ravel()

    def _max_util(self, utils_pairs):
        """Max over a list of (column_series, brace_series) utilisation pairs."""
        vals = []
        for pair in utils_pairs:
            for part in pair:
                a = self._arr(part)
                if a.size:
                    vals.append(np.nanmax(a))
        return float(max(vals)) if vals else 0.0

    # ---- main evaluation --------------------------------------------------
    def _evaluate(self, X, out, *args, **kwargs):
        var = decode_to_var_dict(X)
        x_eng = var_dict_to_xarray(var)

        try:
            # 1) FEA: total mass + per-segment masses
            f, segment_masses = RunAPDL(self.mapdl, var, Misc, opti_settings)

            # 2) Post-processing: utilisation ratios (feasible when u <= 1)
            pp = PostProcessor(var, Misc, opti_settings)

            member_utils = [
                pp.Util_LB(), pp.Util_NF(), pp.Util_S(),
                pp.Util_T(), pp.Util_BR(), pp.Util_IN(),
            ]
            g_util_struct = self._max_util(member_utils) - 1.0

            bs_sig, bs_defl = pp.Util_BS()
            g_util_brace = max(
                float(np.nanmax(self._arr(bs_sig))) if self._arr(bs_sig).size else 0.0,
                float(np.nanmax(self._arr(bs_defl))) if self._arr(bs_defl).size else 0.0,
            ) - 1.0

            # 3) Buckling eigenvalue:  a_cr >= eig_min
            a_cr = float(self._arr(pp.Eigenvalue_1())[0]) + self.eig_min  # Eigenvalue_1 returns a_cr - 4
            g_eig = self.eig_min - a_cr

            # 4) Per-segment mass limit
            seg = self._arr(segment_masses)
            seg = (list(seg) + [0.0] * N_SEG)[:N_SEG]
            g_mass = [float(m) - self.mass_limit for m in seg]

            G = [g_util_struct, g_util_brace, g_eig] + g_mass

            out["F"] = float(f)
            out["G"] = G

            # ---- bookkeeping / logging ------------------------------------
            max_util = max(g_util_struct, g_util_brace) + 1.0
            self.logger.log_evaluation(
                x_eng, f, segment_masses,
                v_agg=max_util, g_max=max(G), c_value=None,
            )

            feasible = all(g <= 1e-6 for g in G)
            if feasible and (self.best is None or f < self.best["f"]):
                self.best = {
                    "f": float(f),
                    "var": var,
                    "x_eng": x_eng,
                    "segment_masses": list(seg),
                    "max_util": max_util,
                    "a_cr": a_cr,
                    "util_report": {
                        "Util_LB": pp.Util_LB(),
                        "Util_NF": pp.Util_NF(),
                        "Util_S": pp.Util_S(),
                        "Util_T": pp.Util_T(),
                        "Util_BR": pp.Util_BR(),
                        "Util_IN": pp.Util_IN(),
                        "Util_BS_sig": bs_sig,
                        "Util_BS_defl": bs_defl,
                    },
                }

        except Exception as exc:  # FEA / geometry failure -> infeasible, discarded
            print(f"[GA] Evaluation failed, penalising individual: {exc}", flush=True)
            out["F"] = float(GA_Settings["fail_penalty"])
            out["G"] = [float(GA_Settings["fail_penalty"])] * self.n_g


# =============================================================================
# 5. PER-GENERATION CALLBACK  (iteration logging)
# =============================================================================
class GACallback(Callback):
    def __init__(self, problem):
        super().__init__()
        self.problem = problem

    def notify(self, algorithm):
        logger = self.problem.logger
        gen = algorithm.n_gen
        n_eval = algorithm.evaluator.n_eval

        # opt holds the current best-so-far (feasible if one exists)
        opt = algorithm.opt[0] if algorithm.opt is not None else None
        f_best = float(opt.F[0]) if opt is not None else float("nan")
        cv_best = float(opt.CV[0]) if (opt is not None and opt.CV is not None) else float("nan")

        logger.log_line("-" * 80)
        logger.log_line(f"[GENERATION {gen}]  evals={n_eval}")
        logger.log_line(f"  pop-best objective : {f_best:.3f} kg   (CV={cv_best:.4g})")
        if self.problem.best is not None:
            b = self.problem.best
            logger.log_line(
                f"  best feasible      : {b['f']:.3f} kg   "
                f"max_util={b['max_util']:.3f}  a_cr={b['a_cr']:.3f}"
            )
        logger.log_line("")


# =============================================================================
# 6. RUN
# =============================================================================
def main():
    start = time.time()

    # ---- Logger (reuses the same artefacts as the gradient-based runs) -----
    init_var = decode_to_var_dict(initial_genome())
    x0_eng = var_dict_to_xarray(init_var)
    bounds_for_log = []
    for name in ENG_NAMES:
        if name == "rad":
            bounds_for_log.append(rad_bounds)
        elif name.startswith("d0"):
            bounds_for_log.append(column_diameter_bounds)
        elif name.startswith("t0"):
            bounds_for_log.append(column_thickness_bounds)
        elif name.startswith("d1"):
            bounds_for_log.append(brace_diameter_bounds)
        else:  # t1*
            bounds_for_log.append(brace_thickness_bounds)

    logger = OptimizationLogger(
        x0=x0_eng,
        bounds=bounds_for_log,
        method="GA (pymoo MixedVariableGA)",
        var_names=ENG_NAMES,
        options={
            "pop_size": GA_Settings["pop_size"],
            "n_gen": GA_Settings["n_gen"],
            "seed": GA_Settings["seed"],
            "eigenvalue_min": GA_Settings["eigenvalue_min"],
            "column_catalogue_size": len(COL_CATALOGUE),
            "brace_catalogue_size": len(BRACE_CATALOGUE),
            "multi_size_columns": opti_settings["multi_size_columns"],
            "multi_size_braces": opti_settings["multi_size_braces"],
            "brace_split": opti_settings["brace_split"],
        },
        save_folder=Misc["save_folder"],
    )

    # ---- Launch MAPDL once (same switches as Main.py) ----------------------
    tic_lic = time.time()
    mapdl = launch_mapdl(
        run_location="Ansout",
        log_apdl="apdl_log",
        override=True,
        nproc=7,
        additional_switches="-p ansys -smp",
    )
    print(f"License opened in: {time.time() - tic_lic:.2f} s", flush=True)

    try:
        problem = MastGAProblem(mapdl, logger)

        # Seed the initial population: random mixed-variable individuals with the
        # known Main.py guess injected as member 0 (a sensible, feasible-ish start).
        sampling = MixedVariableSampling()
        init_pop = sampling(problem, GA_Settings["pop_size"])
        init_pop[0].set("X", initial_genome())

        algorithm = MixedVariableGA(pop_size=GA_Settings["pop_size"], sampling=init_pop)

        result = minimize(
            problem,
            algorithm,
            termination=("n_gen", GA_Settings["n_gen"]),
            seed=GA_Settings["seed"],
            callback=GACallback(problem),
            verbose=True,
            save_history=False,
        )
    finally:
        mapdl.exit()

    # ---- Final reporting ---------------------------------------------------
    best = problem.best
    logger.log_line("=" * 80)
    logger.log_line("GENETIC ALGORITHM - FINAL RESULT")
    logger.log_line("=" * 80)

    if best is None:
        logger.log_line("No FEASIBLE design found. Increase pop_size / n_gen or relax bounds.")
        success = False
        best_x = np.full(len(ENG_NAMES), np.nan)
        best_f = float("nan")
        message = "No feasible individual found"
    else:
        success = True
        best_x = best["x_eng"]
        best_f = best["f"]
        message = "Feasible optimum found by GA"
        logger.log_line(f"Best mass            : {best_f:.4f} kg")
        logger.log_line(f"Max utilisation      : {best['max_util']:.4f}")
        logger.log_line(f"First eigenvalue a_cr: {best['a_cr']:.4f}  (limit {GA_Settings['eigenvalue_min']})")
        logger.log_line(f"Segment masses [kg]  : " +
                        ", ".join(f"{m:.3f}" for m in best["segment_masses"]))
        logger.log_line("")
        logger.log_line("Optimal sections (real production profiles):")
        for name in ENG_NAMES:
            logger.log_line(f"  {name:8s} = {best['var'][name]['value']:.3f}")
        logger.log_line("")
        logger.log_line("Per-mode utilisation (max column / max brace):")
        logger.log_utilization(best["util_report"])

    # Build a result-like object so OptimizationLogger.finalize works unchanged
    result_like = types.SimpleNamespace(
        x=best_x,
        fun=best_f,
        nit=int(getattr(result.algorithm, "n_gen", GA_Settings["n_gen"])),
        nfev=int(result.algorithm.evaluator.n_eval),
        success=success,
        message=message,
        status=0 if success else -1,
    )
    logger.finalize(result_like)

    print(f"\nGA completed in {time.time() - start:.2f} seconds.", flush=True)
    print("Best objective:", best_f)
    print("Message:", message)
    print("TXT log file:", logger.txt_path)
    print("Objective CSV:", logger.csv_path)

    return result_like, logger.txt_path, logger.csv_path


if __name__ == "__main__":
    main()
