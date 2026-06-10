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
#   PyGAD (https://pygad.readthedocs.io) -- a widely used, lightweight GA
#   library.  The genome is a flat vector: gene 0 is the continuous radius
#   `rad`; every other gene is an *integer index* into a production catalogue
#   (S700 for columns, S235 for braces).  Per-gene `gene_space` keeps every
#   section choice on a real, buyable, Class-2 profile, and `gene_type` keeps
#   `rad` a float while the index genes stay integers.
#
#   GA configuration (as requested):
#       parent_selection_type = "sss"          (steady-state selection)
#       crossover_type        = "single_point"
#       mutation_type         = "random"
#       mutation_percent_genes = 10
#       stop_criteria         = ["reach_1", "saturate_10"]
#       save_solutions        = True
#
#   Constraint handling: PyGAD maximises a single fitness, so constraints are
#   folded in with an exterior penalty.  Fitness = -(mass + P * violation),
#   where `violation` is the summed positive part of every g(x) <= 0 constraint.
#   Among feasible designs (violation = 0) fitness reduces to -mass, so the GA
#   drives mass down.  NOTE: because fitness is <= 0 for any real design, the
#   "reach_1" stop criterion never fires here -- early stopping is governed by
#   "saturate_10" (no improvement for 10 generations) and num_generations.
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

# ---- PyGAD ------------------------------------------------------------------
try:
    import pygad
except ImportError as exc:  # pragma: no cover - guidance only
    raise SystemExit(
        "PyGAD is required for main_gen.py. Install it with:\n"
        "    python -m pip install pygad\n"
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
    "multi_size_columns": False,  # Columns may differ per segment (True) or not (False)
    "multi_size_braces": False,   # Braces may differ per segment (True) or not (False)
    "brace_split": False,         # Horizontal/Cross braces sized separately (True) or not (False)
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
    "pop_size": 60,           # sol_per_pop  (individuals per generation) -- more diversity
    "n_gen": 120,             # num_generations  (FEA evals <= pop_size * n_gen)
    "seed": 1,                # RNG seed for reproducibility
    "eigenvalue_min": 4.0,    # Required first positive buckling eigenvalue (a_cr >= 4)
    "penalty": 1.0e4,         # Exterior penalty weight P on constraint violation
    "fail_penalty": 1.0e12,   # Fitness penalty for FEA failures (very bad)
    # ---- fixed GA operators (as requested) --------------------------------
    "parent_selection_type": "sss",
    "crossover_type": "single_point",
    "mutation_type": "random",
    "mutation_percent_genes": 20,
    # KEY for coupled moves: force >=2 genes to mutate per offspring so one child
    # can change `rad` AND a section index together. Without this, only 1 of the
    # 3 genes mutates (int(0.20*3)=0 -> clamped to 1), and the optimum at
    # rad=264 / col=(48.3,2.5) is unreachable because each single-gene step is
    # either infeasible or fitness-neutral. Takes precedence over
    # mutation_percent_genes when not None.
    "mutation_num_genes": 2,
    "stop_criteria": ["reach_1", "saturate_40"],   # probe 40 stagnant gens before quitting
    "save_solutions": True,
    # Random mutation step applied to the continuous `rad` gene (index genes
    # mutate by re-sampling their catalogue gene_space instead). Wider = bigger
    # radius jumps to cross the 228 -> 264 mm basin.
    "rad_mutation_step": 40.0,
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


# Build the list of genome "genes" for the active layout. ----------------------
# Columns
if opti_settings["multi_size_columns"]:
    COL_GENES = [f"col_{i}" for i in range(1, N_SEG + 1)]
else:
    COL_GENES = ["col"]

# Braces
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

# Genome order:  rad (float)  +  column index genes (int)  +  brace index genes (int)
GENE_NAMES = ["rad"] + COL_GENES + BRACE_GENES
COL_GENE_SET = set(COL_GENES)
BRACE_GENE_SET = set(BRACE_GENES)


def build_gene_space_and_types():
    """PyGAD per-gene search space + dtype, in GENE_NAMES order.

    rad        -> continuous range {'low', 'high'},  float
    col_*      -> discrete catalogue indices [0..Ncol-1],  int
    b*/bh*/bc* -> discrete catalogue indices [0..Nbrace-1], int
    """
    gene_space = []
    gene_type = []
    for name in GENE_NAMES:
        if name == "rad":
            gene_space.append({"low": rad_bounds[0], "high": rad_bounds[1]})
            gene_type.append(float)
        elif name in COL_GENE_SET:
            gene_space.append(list(range(len(COL_CATALOGUE))))
            gene_type.append(int)
        else:  # brace gene
            gene_space.append(list(range(len(BRACE_CATALOGUE))))
            gene_type.append(int)
    return gene_space, gene_type


def decode_to_var_dict(genome):
    """Map a genome (dict gene_name -> value) -> nested var dict for the pipeline.

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


def solution_to_genome_dict(solution):
    """Zip a flat PyGAD solution vector into a {gene_name: value} dict."""
    return dict(zip(GENE_NAMES, solution))


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
    """Seed genome (as a {gene: value} dict) from Main.py's initial guess,
    snapped to the nearest valid catalogue entry."""
    g = {"rad": float(init_rad)}
    col_idx = _nearest_index(COL_CATALOGUE, init_col_dt)
    brc_idx = _nearest_index(BRACE_CATALOGUE, init_brace_dt)
    for gene in COL_GENES:
        g[gene] = col_idx
    for gene in BRACE_GENES:
        g[gene] = brc_idx
    return g


def build_initial_population(pop_size, rng):
    """Build a (pop_size x n_genes) float array. Row 0 is the Main.py seed;
    the rest are random samples from each gene's space."""
    seed = initial_genome()
    pop = [[float(seed[name]) for name in GENE_NAMES]]
    for _ in range(pop_size - 1):
        row = []
        for name in GENE_NAMES:
            if name == "rad":
                row.append(float(rng.uniform(*rad_bounds)))
            elif name in COL_GENE_SET:
                row.append(float(rng.integers(0, len(COL_CATALOGUE))))
            else:
                row.append(float(rng.integers(0, len(BRACE_CATALOGUE))))
        pop.append(row)
    return np.array(pop, dtype=float)


# =============================================================================
# 4. GA PROBLEM  (wraps the FEA pipeline for PyGAD)
# =============================================================================
class MastGAProblem:
    """Holds run state and provides PyGAD's `fitness_func` and `on_generation`.

    Objective : total assembly mass [kg]  (minimise -> maximise -mass)
    Constraints (g(x) <= 0 feasible):
        g_util_struct : max EC3 member utilisation - 1            (LB/NF/S/T/BR/IN)
        g_util_brace  : max brace-step utilisation - 1            (stress + defl)
        g_eig         : eigenvalue_min - a_cr                     (a_cr >= 4)
        g_mass_i      : seg_mass_i - segment_mass_limit           (one per segment)
    Class-2 is guaranteed by the catalogue filter, so it is not re-imposed.
    Fitness = -(mass + P * sum(max(g_i, 0))).
    """

    def __init__(self, mapdl, logger):
        self.mapdl = mapdl
        self.logger = logger
        self.eig_min = float(GA_Settings["eigenvalue_min"])
        self.mass_limit = float(opti_settings["segment_mass_limit"])
        self.penalty = float(GA_Settings["penalty"])
        self.fail_penalty = float(GA_Settings["fail_penalty"])
        self.n_g = 3 + N_SEG  # 2 utilisation + 1 eigenvalue + per-segment mass

        # Human-readable name for every constraint row (same order as G)
        self.constraint_names = ["util_struct", "util_brace", "eig"] + \
            [f"seg_mass_{i}" for i in range(1, N_SEG + 1)]

        self.cache = {}   # genome-key -> record dict (avoids repeated FEA)
        self.n_eval = 0   # actual FEA evaluations performed
        self.best = None  # best feasible design seen so far

    # ---- helpers ----------------------------------------------------------
    @staticmethod
    def _arr(v):
        return np.asarray(v, dtype=float).ravel()

    @staticmethod
    def _key(solution):
        # rad rounded to a fine grid; index genes are exact integers
        return (round(float(solution[0]), 4),) + tuple(int(round(v)) for v in solution[1:])

    def _max_util(self, utils_pairs):
        vals = []
        for pair in utils_pairs:
            for part in pair:
                a = self._arr(part)
                if a.size:
                    vals.append(np.nanmax(a))
        return float(max(vals)) if vals else 0.0

    # ---- full evaluation of one design ------------------------------------
    def _evaluate(self, solution):
        gd = solution_to_genome_dict(solution)
        var = decode_to_var_dict(gd)
        x_eng = var_dict_to_xarray(var)

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

        # 3) Buckling eigenvalue:  a_cr >= eig_min  (Eigenvalue_1 returns a_cr - 4)
        a_cr = float(self._arr(pp.Eigenvalue_1())[0]) + 4.0
        g_eig = self.eig_min - a_cr

        # 4) Per-segment mass limit
        seg = self._arr(segment_masses)
        seg = (list(seg) + [0.0] * N_SEG)[:N_SEG]
        g_mass = [float(m) - self.mass_limit for m in seg]

        G = [g_util_struct, g_util_brace, g_eig] + g_mass
        max_util = max(g_util_struct, g_util_brace) + 1.0
        return f, G, seg, max_util, a_cr, x_eng, var, pp, bs_sig, bs_defl

    # ---- PyGAD fitness function -------------------------------------------
    def fitness_func(self, ga_instance, solution, solution_idx):
        key = self._key(solution)
        cached = self.cache.get(key)
        if cached is not None:
            return cached["fitness"]

        try:
            (f, G, seg, max_util, a_cr, x_eng, var,
             pp, bs_sig, bs_defl) = self._evaluate(solution)

            violation = sum(max(g, 0.0) for g in G)
            fitness = -(float(f) + self.penalty * violation)  # penalty is the exterior penalty weight P on constraint violation https://pdf.sciencedirectassets.com/271868/1-s2.0-S0045782500X01028/1-s2.0-S0045782599003898/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEAMaCXVzLWVhc3QtMSJHMEUCIDfzwf3l2saCVtOui7FJ7zxSGdGQ31ud2XVtSNkhnV9NAiEAt%2Ffix1WgzzkSOwa9y7KvnEojabmqilZjZqW2hPpfJI0qvAUIzP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDNaLZN5kno%2B3zTwRjCqQBYpdf5rK4E12MdZCktXbRRtzdXt4ODEDLvIouMPVhaC78H2zY6qsX2XZAOjLCKxWGzFtqWw4IEumad8w0c5pHySYrpdFSkzWUT%2BR1torxev%2Fj6cJ9fZuQEY%2BRtJw8cNjJf%2FUlSy%2Bt6%2B0P2sP1kGUAHv9IvB8ZiB0Km%2B5wBXaahVj5dk42%2B83XxYZf2n9qFCsH3pcwhCyjISomAc4rEENsaItsEnfJ%2BkcPl5795764voHZpD8j77SCUGegFbBmGefRmZdL%2BgfAt4%2BJyDBE7%2B6szeeD4M2lOuOG9uOTDfmcgsSsYX0i2mTglBFFQY0Co7xMU%2Fd5U4qlEMYQ8kJQPW4hKBo8DJlmtb5G1q4zLZ9gXMET9VYu%2BBEcKpph%2Fww%2FDnJNEUgkqeSoJROjAtsTx8nGOh320t6G4skWG1cGpIC5UFeSGGo43GOV9G%2F2XNAldj0YeKioFlTSnTRYaDJQaPfZv31zZ2sD%2BPH3WS9xHDXoqwbFUusYC3zWq2WPbtiwKStIuUmRRgTBT9EOJEeH3nvw6xNZQ2lLqHyHRoyOQDBb3ABrJnZsz0q261yn%2BIXOCIhETl63d8XFpb0skoMXJbOzRRfcsok5e0yLEuRtEqiRBEwj2pOmrYcRV41rCLUoAwbZQHUi78mkShc8DGoNoSt%2FnQqFQgXUoipLUXKvBZTVFpL0z55iCHh%2Fbo49hfqejGbvMYbEuE%2B4DVKfhpnBKTc8nGoPbB4v%2BLbdg1WMXZJ%2BAJ8kAnQSRQdTZwsf1S5OaAlyguUk8XHgUyXb5Pujcq7R8aa5x2aqfCAjYCwpifx8Z3nnuOYjWLy24vSNTC5TMQlse3uQwyi6CkvnEqYDdMUhFUluvCGP7FK2K3IejXQBKMhMK3Un9EGOrEBacBhQqMEi7RBJW1tvubCe4xMj19%2Fo4B3pABQpDrazEiiZldYlc48qZuP9tJGgNue%2BuT4ZMJct6M%2Ff6IaFhJkvFHuy3PfApDw%2BUk%2Fr8YIV%2FeLYdWzFCj0UfqCDpYVAC6fP7s8xnTec6zzxy5XvU5%2BchNlLidPY6p16ncQbiUU8jUGMjGDxNhirbvlTmo4vU%2BSmJ2TkSi%2FNBovaxoPgzzF8QqSZpuIBxn5iLccBcJ6c2GL&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20260609T111903Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY4PE3Y34E%2F20260609%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=e6812d8d940d50a93347a9f6acb5dc07e6e2e6a4180c47c9e2c61c9ae01b1692&hash=012ab8b9e8ed2023c4fbcca33c3fadd983f94945bfe14f21b6363535cb8bfd5a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0045782599003898&tid=spdf-b64b4c0e-9839-4b7c-a335-5f0bfe46536f&sid=b2f722af6bd6b64a182a39e81d63b69e799dgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=020c05525b000204515a&rr=a08fbd3508cf8a28&cc=dk
            feasible = violation <= 1e-6

            self.n_eval += 1
            self.logger.log_evaluation(
                x_eng, f, seg, v_agg=max_util, g_max=max(G), c_value=None,
            )

            # Per-solution record (used by on_generation to print g(x) and x)
            record = {
                "fitness": fitness,
                "f": float(f),
                "G": [float(g) for g in G],
                "violation": float(violation),
                "x_eng": x_eng,
                "max_util": float(max_util),
                "a_cr": float(a_cr),
                "feasible": bool(feasible),
            }

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

        except Exception as exc:  # FEA / geometry failure -> very bad fitness
            print(f"[GA] Evaluation failed, penalising individual: {exc}", flush=True)
            gd = solution_to_genome_dict(solution)
            try:
                x_eng = var_dict_to_xarray(decode_to_var_dict(gd))
            except Exception:
                x_eng = None
            record = {
                "fitness": -self.fail_penalty,
                "f": float("nan"),
                "G": None,
                "violation": float("inf"),
                "x_eng": x_eng,
                "max_util": float("nan"),
                "a_cr": float("nan"),
                "feasible": False,
            }

        self.cache[key] = record
        return record["fitness"]

    # ---- PyGAD per-generation callback ------------------------------------
    def on_generation(self, ga_instance):
        gen = ga_instance.generations_completed
        try:
            best_sol, fit, _ = ga_instance.best_solution(ga_instance.last_generation_fitness)
        except Exception:
            best_sol, fit = None, float("nan")

        self.logger.log_line("-" * 80)
        self.logger.log_line(f"[GENERATION {gen}]  FEA evals so far: {self.n_eval}")
        self.logger.log_line(f"  pop-best fitness   : {fit:.4g}")

        # Details (design variables + constraints) of this generation's best
        rec = self.cache.get(self._key(best_sol)) if best_sol is not None else None
        if rec is not None:
            self.logger.log_line(
                f"  objective (mass)   : {rec['f']:.3f} kg   "
                f"max_util={rec['max_util']:.3f}  a_cr={rec['a_cr']:.3f}"
            )

            self.logger.log_line("  design variables:")
            if rec["x_eng"] is not None:
                self.logger.log_line(
                    "    " + ",  ".join(
                        f"{name}={val:.3f}"
                        for name, val in zip(ENG_NAMES, rec["x_eng"])
                    )
                )

            self.logger.log_line("  constraints  g(x) <= 0 feasible:")
            if rec["G"] is not None:
                for name, g in zip(self.constraint_names, rec["G"]):
                    flag = "OK" if g <= 1e-6 else "VIOLATED"
                    self.logger.log_line(f"    {name:12s} = {g:+.4f}   [{flag}]")
                self.logger.log_line(
                    f"    total violation = {rec['violation']:.4g}   "
                    f"({'FEASIBLE' if rec['feasible'] else 'INFEASIBLE'})"
                )
            else:
                self.logger.log_line("    (FEA failed for this individual)")

        if self.best is not None:
            b = self.best
            self.logger.log_line(
                f"  best feasible so far: {b['f']:.3f} kg   "
                f"max_util={b['max_util']:.3f}  a_cr={b['a_cr']:.3f}"
            )
        self.logger.log_line("")


# =============================================================================
# 5. RUN
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
        method="GA (PyGAD)",
        var_names=ENG_NAMES,
        options={
            "pop_size": GA_Settings["pop_size"],
            "n_gen": GA_Settings["n_gen"],
            "seed": GA_Settings["seed"],
            "parent_selection_type": GA_Settings["parent_selection_type"],
            "crossover_type": GA_Settings["crossover_type"],
            "mutation_type": GA_Settings["mutation_type"],
            "mutation_percent_genes": GA_Settings["mutation_percent_genes"],
            "mutation_num_genes": GA_Settings.get("mutation_num_genes"),
            "rad_mutation_step": GA_Settings["rad_mutation_step"],
            "stop_criteria": GA_Settings["stop_criteria"],
            "eigenvalue_min": GA_Settings["eigenvalue_min"],
            "penalty": GA_Settings["penalty"],
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

    ga_instance = None
    try:
        problem = MastGAProblem(mapdl, logger)

        rng = np.random.default_rng(GA_Settings["seed"])
        initial_population = build_initial_population(GA_Settings["pop_size"], rng)
        gene_space, gene_type = build_gene_space_and_types()

        # Mutation count: prefer an explicit gene count (needed for coupled moves
        # on this tiny genome); fall back to the percentage otherwise.
        if GA_Settings.get("mutation_num_genes") is not None:
            mutation_kwargs = {"mutation_num_genes": GA_Settings["mutation_num_genes"]}
        else:
            mutation_kwargs = {"mutation_percent_genes": GA_Settings["mutation_percent_genes"]}

        ga_instance = pygad.GA(
            num_generations=GA_Settings["n_gen"],
            num_parents_mating=max(2, GA_Settings["pop_size"] // 2),
            fitness_func=problem.fitness_func,
            on_generation=problem.on_generation,
            initial_population=initial_population,
            gene_space=gene_space,
            gene_type=gene_type,
            parent_selection_type=GA_Settings["parent_selection_type"],
            crossover_type=GA_Settings["crossover_type"],
            mutation_type=GA_Settings["mutation_type"],
            stop_criteria=GA_Settings["stop_criteria"],
            **mutation_kwargs,
            save_solutions=GA_Settings["save_solutions"],
            random_seed=GA_Settings["seed"],
            # Continuous-gene mutation step (index genes re-sample gene_space)
            random_mutation_min_val=-GA_Settings["rad_mutation_step"],
            random_mutation_max_val=GA_Settings["rad_mutation_step"],
            suppress_warnings=True,
        )

        ga_instance.run()
    finally:
        mapdl.exit()

    # ---- Final reporting ---------------------------------------------------
    best = problem.best
    logger.log_line("=" * 80)
    logger.log_line("GENETIC ALGORITHM (PyGAD) - FINAL RESULT")
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
        logger.log_line("Segment masses [kg]  : " +
                        ", ".join(f"{m:.3f}" for m in best["segment_masses"]))
        logger.log_line("")
        logger.log_line("Optimal sections (real production profiles):")
        for name in ENG_NAMES:
            logger.log_line(f"  {name:8s} = {best['var'][name]['value']:.3f}")
        logger.log_line("")
        logger.log_line("Per-mode utilisation (max column / max brace):")
        logger.log_utilization(best["util_report"])

    n_gen_done = int(getattr(ga_instance, "generations_completed", GA_Settings["n_gen"])) \
        if ga_instance is not None else 0

    # Build a result-like object so OptimizationLogger.finalize works unchanged
    result_like = types.SimpleNamespace(
        x=best_x,
        fun=best_f,
        nit=n_gen_done,
        nfev=int(problem.n_eval),
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
