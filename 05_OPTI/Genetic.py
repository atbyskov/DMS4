## GENETIC ALGORITHM OPTIMIZATION (Genetic.py) ##
# Standard PyGAD implementation for the lattice mast optimisation.
#
# Uses PyGAD with its stock default settings:
#   parent_selection_type = "sss"   (steady-state selection)
#   keep_parents          = 1
#   crossover_type        = "single_point"
#   mutation_type         = "random"
#   mutation_percent_genes= 10
#
# Constraints are folded into the fitness function via an exterior penalty:
#   fitness = -(mass + P * sum(max(g_i, 0)))
# where g_i <= 0 means feasible.  Among feasible designs fitness reduces to
# -mass, so PyGAD (which maximises) drives mass down.
#
# Run from the 05_OPTI directory:
#   cd 05_OPTI
#   python Genetic.py

import sys
print(sys.version, flush=True)
import os
import time
import types
import numpy as np

from ansys.mapdl.core import launch_mapdl

try:
    import pygad
except ImportError as exc:
    raise SystemExit(
        "PyGAD is required.  Install with:\n"
        "    python -m pip install pygad\n"
        f"(import error: {exc})"
    )

from MyAPDLCall import RunAPDL
from Post_Process import PostProcessor
from opt_logger import OptimizationLogger


# =============================================================================
# 1. PRODUCTION PROFILE CATALOGUES
# =============================================================================
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

S235_RAW = [
    (6,   [1, 1.5, 2]),
    (8,   [1, 1.5, 2]),
    (10,  [1, 1.5, 2, 2.5]),
    (12,  [1, 1.5, 2, 2.5]),
    (13,  [1, 1.5, 2, 2.5]),
    (14,  [1, 1.5, 2, 2.5, 3, 3.5]),
    (15,  [1, 1.5, 2, 2.5, 3]),
    (16,  [1, 1.5, 2, 2.5, 3, 3.5, 4]),
    (17,  [1, 1.5, 2, 2.5]),
    (18,  [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]),
    (20,  [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]),
    (21,  [1, 1.5, 2, 2.5, 3, 4]),
    (22,  [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]),
    (23,  [1, 1.5, 2]),
    (24,  [1, 1.5, 2, 3, 4]),
    (25,  [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7]),
    (28,  [1, 1.5, 2, 2.5, 3, 4, 5, 6]),
    (30,  [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 7.5]),
    (32,  [1.5, 2, 3, 3.5, 4, 5, 6]),
    (34,  [1.5, 2, 4]),
    (35,  [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 7.5, 8, 9]),
    (36,  [1.5, 2, 3, 5, 8]),
    (38,  [1.5, 2, 3, 4, 4.5, 5, 6, 7, 8]),
    (40,  [1.5, 2, 2.5, 3, 4, 4.5, 5, 5.5, 6, 7, 7.5, 8, 10]),
    (42,  [1.5, 2, 3, 4, 5, 7]),
    (45,  [1.5, 2, 2.5, 3, 4, 4.5, 5, 6, 7, 10]),
    (50,  [1.5, 2, 3, 4, 5, 6, 7, 8, 10]),
    (52,  [5]),
    (55,  [1.5, 2, 3, 4, 5, 6, 8, 10]),
    (60,  [1.5, 2, 3, 4, 5, 6, 8, 9, 10]),
    (62,  [6]),
    (65,  [4, 6, 10]),
    (70,  [2, 3, 4, 5, 6, 8, 10]),
    (80,  [2, 3, 4, 5, 6, 8, 10, 12]),
    (90,  [5]),
    (100, [3, 4, 5, 6]),
]
S235_DATA = [(float(d), float(t)) for d, t_list in S235_RAW for t in t_list]


def build_catalogue(data, d_bounds, t_bounds, f_y, label):
    """Filter (d, t) table to Class-2 entries inside the design bounds."""
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
    catalogue.sort(key=lambda dt: dt[0] ** 2 - (dt[0] - 2 * dt[1]) ** 2)
    if not catalogue:
        raise ValueError(
            f"No {label} catalogue points satisfy bounds + Class-2 filter "
            f"(d in {d_bounds}, t in {t_bounds}, d/t <= {class2_limit:.3f})."
        )
    print(
        f"{label} catalogue: {len(catalogue)} Class-2 profiles in bounds "
        f"(d/t <= {class2_limit:.2f}).",
        flush=True,
    )
    return catalogue


# =============================================================================
# 2. PROBLEM SETTINGS
# =============================================================================
opti_settings = {
    "n_mast_segments":    5,
    "mast_segment_height": 810,
    "segment_mass_limit": 23,
    "multi_size_columns": False,
    "multi_size_braces":  False,
    "brace_split":        False,
}

Misc = {
    "esize":       3,
    "Hor_Force":   502.52,
    "Ver_Force":   -25.13e3,
    "f_y":         700,
    "f_y_brace":   235,
    "E_mod":       200e3,
    "W_Force":     -3.751e3,
    "SW_filename": "LWC_L1_LINES.IGS",
    "save_folder": "Optimization_Logs",
}

column_diameter_bounds  = (48.3,  114.3)
column_thickness_bounds = (2.5,   5.0)
brace_diameter_bounds   = (10.0,  50.0)
brace_thickness_bounds  = (1.0,   4.0)
rad_bounds              = (150.0, 350.0)

init_rad       = 202.07
init_col_dt    = (76.1, 3.0)
init_brace_dt  = (26.9, 2.3)

GA_Settings = {
    "pop_size":    9,    # sol_per_pop
    "n_gen":       450,    # num_generations
    "seed":        5,
    "eigenvalue_min": 4.0,
    "penalty":     1000,
    "fail_penalty": 1.0e12,
    # --- stock PyGAD defaults ------------------------------------------------
    "parent_selection_type": "sss",
    "keep_parents":           1,
    "crossover_type":         "single_point",
    "mutation_type":          "random",
    "mutation_percent_genes": 34,
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

# Genome: [rad, col_index, brace_index]
# (one column group and one brace group since multi_size_* = False)
GENE_NAMES = ["rad", "col", "b"]

GENE_SPACE = [
    {"low": rad_bounds[0], "high": rad_bounds[1]},  # rad: continuous
    list(range(len(COL_CATALOGUE))),                 # col: integer index
    list(range(len(BRACE_CATALOGUE))),               # b:   integer index
]
GENE_TYPE = [float, int, int]

# Engineering variable names for logging
ENG_NAMES = ["rad", "d0", "t0", "d1", "t1"]


def _nearest_index(catalogue, target_dt):
    td, tt = target_dt
    dists = [(d - td) ** 2 + (t - tt) ** 2 for (d, t) in catalogue]
    return int(np.argmin(dists))


def decode_solution(solution):
    """Map a flat [rad, col_idx, brace_idx] vector to the var dict the pipeline expects."""
    rad   = float(solution[0])
    col_i = int(round(solution[1]))
    brc_i = int(round(solution[2]))

    col_i = max(0, min(col_i, len(COL_CATALOGUE) - 1))
    brc_i = max(0, min(brc_i, len(BRACE_CATALOGUE) - 1))

    d0, t0 = COL_CATALOGUE[col_i]
    d1, t1 = BRACE_CATALOGUE[brc_i]

    return {
        "rad": {"value": rad,        "active": True},
        "d0":  {"value": float(d0),  "active": True},
        "t0":  {"value": float(t0),  "active": True},
        "d1":  {"value": float(d1),  "active": True},
        "t1":  {"value": float(t1),  "active": True},
    }


def var_to_x(var):
    return np.array([var[n]["value"] for n in ENG_NAMES], dtype=float)


def initial_solution():
    """Return a flat genome vector seeded from the Main.py initial guess."""
    col_i = _nearest_index(COL_CATALOGUE,   init_col_dt)
    brc_i = _nearest_index(BRACE_CATALOGUE, init_brace_dt)
    return [float(init_rad), float(col_i), float(brc_i)]


# =============================================================================
# 4. FITNESS FUNCTION  (stateful via a mutable container)
# =============================================================================
# PyGAD requires a module-level or instance-level callable; we use a simple
# class to carry the shared MAPDL instance, logger, cache, and best result.

class _State:
    def __init__(self, mapdl, logger):
        self.mapdl    = mapdl
        self.logger   = logger
        self.cache    = {}
        self.n_eval   = 0
        self.best     = None
        self.eig_min  = float(GA_Settings["eigenvalue_min"])
        self.mass_lim = float(opti_settings["segment_mass_limit"])
        self.penalty  = float(GA_Settings["penalty"])
        self.fail_pen = float(GA_Settings["fail_penalty"])

    @staticmethod
    def _key(sol):
        return (round(float(sol[0]), 4), int(round(sol[1])), int(round(sol[2])))

    def _arr(self, v):
        return np.asarray(v, dtype=float).ravel()

    def _max_util(self, utils_pairs):
        vals = []
        for pair in utils_pairs:
            for part in pair:
                a = self._arr(part)
                if a.size:
                    vals.append(float(np.nanmax(a)))
        return max(vals) if vals else 0.0

    def evaluate(self, solution):
        key    = self._key(solution)
        cached = self.cache.get(key)
        if cached is not None:
            return cached

        try:
            var    = decode_solution(solution)
            x_eng  = var_to_x(var)

            f, segment_masses = RunAPDL(self.mapdl, var, Misc, opti_settings)

            pp = PostProcessor(var, Misc, opti_settings)

            member_utils = [
                pp.Util_LB(), pp.Util_NF(), pp.Util_S(),
                pp.Util_T(),  pp.Util_BR(), pp.Util_IN(),
            ]
            g_util_struct = self._max_util(member_utils) - 1.0

            bs_sig, bs_defl = pp.Util_BS()
            bs_sig_arr  = self._arr(bs_sig)
            bs_defl_arr = self._arr(bs_defl)
            g_util_brace = max(
                float(np.nanmax(bs_sig_arr))  if bs_sig_arr.size  else 0.0,
                float(np.nanmax(bs_defl_arr)) if bs_defl_arr.size else 0.0,
            ) - 1.0

            a_cr  = float(self._arr(pp.Eigenvalue_1())[0]) + 4.0
            g_eig = self.eig_min - a_cr

            seg = list((list(self._arr(segment_masses)) + [0.0] * N_SEG)[:N_SEG])
            g_mass = [float(m) - self.mass_lim for m in seg]

            G         = [g_util_struct, g_util_brace, g_eig] + g_mass
            violation = sum(max(g, 0.0) for g in G)
            fitness   = -(float(f) + self.penalty * violation)
            feasible  = violation <= 1e-6
            max_util  = max(g_util_struct, g_util_brace) + 1.0

            self.n_eval += 1
            self.logger.log_evaluation(
                x_eng, f, seg, v_agg=max_util, g_max=max(G), c_value=None,
            )

            record = {
                "fitness":  fitness,
                "f":        float(f),
                "G":        [float(g) for g in G],
                "violation": float(violation),
                "x_eng":    x_eng,
                "max_util": float(max_util),
                "a_cr":     float(a_cr),
                "feasible": bool(feasible),
            }

            if feasible and (self.best is None or f < self.best["f"]):
                self.best = {
                    "f":              float(f),
                    "var":            var,
                    "x_eng":          x_eng,
                    "segment_masses": seg,
                    "max_util":       max_util,
                    "a_cr":           a_cr,
                    "util_report": {
                        "Util_LB":      pp.Util_LB(),
                        "Util_NF":      pp.Util_NF(),
                        "Util_S":       pp.Util_S(),
                        "Util_T":       pp.Util_T(),
                        "Util_BR":      pp.Util_BR(),
                        "Util_IN":      pp.Util_IN(),
                        "Util_BS_sig":  bs_sig,
                        "Util_BS_defl": bs_defl,
                    },
                }

        except Exception as exc:
            print(f"[GA] Evaluation failed: {exc}", flush=True)
            record = {
                "fitness":   -self.fail_pen,
                "f":         float("nan"),
                "G":         None,
                "violation": float("inf"),
                "x_eng":     None,
                "max_util":  float("nan"),
                "a_cr":      float("nan"),
                "feasible":  False,
            }

        self.cache[key] = record
        return record


# Module-level state object; assigned in main() before GA runs.
_state: _State = None


def fitness_func(ga_instance, solution, solution_idx):
    return _state.evaluate(solution)["fitness"]


def on_generation(ga_instance):
    gen = ga_instance.generations_completed
    try:
        best_sol, fit, _ = ga_instance.best_solution(ga_instance.last_generation_fitness)
    except Exception:
        best_sol, fit = None, float("nan")

    _state.logger.log_line("-" * 70)
    _state.logger.log_line(f"[GEN {gen}]  evals: {_state.n_eval}  pop-best fitness: {fit:.4g}")

    rec = _state.cache.get(_state._key(best_sol)) if best_sol is not None else None
    if rec is not None and rec["x_eng"] is not None:
        _state.logger.log_line(
            f"  mass={rec['f']:.3f} kg  max_util={rec['max_util']:.3f}  a_cr={rec['a_cr']:.3f}"
        )
        _state.logger.log_line(
            "  x: " + "  ".join(f"{n}={v:.3f}" for n, v in zip(ENG_NAMES, rec["x_eng"]))
        )
        if rec["G"] is not None:
            constraint_names = (
                ["util_struct", "util_brace", "eig"]
                + [f"seg_mass_{i}" for i in range(1, N_SEG + 1)]
            )
            for cname, g in zip(constraint_names, rec["G"]):
                flag = "OK" if g <= 1e-6 else "VIOLATED"
                _state.logger.log_line(f"    {cname:12s} = {g:+.4f}  [{flag}]")
            _state.logger.log_line(
                f"    violation = {rec['violation']:.4g}  "
                f"({'FEASIBLE' if rec['feasible'] else 'INFEASIBLE'})"
            )

    if _state.best is not None:
        b = _state.best
        _state.logger.log_line(
            f"  best feasible: {b['f']:.3f} kg  "
            f"max_util={b['max_util']:.3f}  a_cr={b['a_cr']:.3f}"
        )
    _state.logger.log_line("")


# =============================================================================
# 5. MAIN
# =============================================================================
def main():
    global _state
    start = time.time()

    # ---- Logger ---------------------------------------------------------------
    init_var  = decode_solution(initial_solution())
    x0_eng    = var_to_x(init_var)
    bounds_for_log = [
        rad_bounds,
        column_diameter_bounds,
        column_thickness_bounds,
        brace_diameter_bounds,
        brace_thickness_bounds,
    ]
    logger = OptimizationLogger(
        x0=x0_eng,
        bounds=bounds_for_log,
        method="GA (PyGAD stock)",
        var_names=ENG_NAMES,
        options={
            "pop_size":               GA_Settings["pop_size"],
            "n_gen":                  GA_Settings["n_gen"],
            "seed":                   GA_Settings["seed"],
            "parent_selection_type":  GA_Settings["parent_selection_type"],
            "keep_parents":           GA_Settings["keep_parents"],
            "crossover_type":         GA_Settings["crossover_type"],
            "mutation_type":          GA_Settings["mutation_type"],
            "mutation_percent_genes": GA_Settings["mutation_percent_genes"],
            "eigenvalue_min":         GA_Settings["eigenvalue_min"],
            "penalty":                GA_Settings["penalty"],
            "column_catalogue_size":  len(COL_CATALOGUE),
            "brace_catalogue_size":   len(BRACE_CATALOGUE),
        },
        save_folder=Misc["save_folder"],
    )

    # ---- Launch MAPDL ---------------------------------------------------------
    tic = time.time()
    mapdl = launch_mapdl(
        run_location="Ansout",
        log_apdl="apdl_log",
        override=True,
        nproc=7,
        additional_switches="-p ansys -smp",
    )
    print(f"License opened in {time.time() - tic:.2f} s", flush=True)

    ga_instance = None
    try:
        _state = _State(mapdl, logger)

        # Seed the first individual from the Main.py initial guess; the rest
        # are sampled randomly by PyGAD from gene_space.
        seed_row   = initial_solution()
        rng        = np.random.default_rng(GA_Settings["seed"])
        init_pop   = [seed_row]
        for _ in range(GA_Settings["pop_size"] - 1):
            row = [
                float(rng.uniform(*rad_bounds)),
                float(rng.integers(0, len(COL_CATALOGUE))),
                float(rng.integers(0, len(BRACE_CATALOGUE))),
            ]
            init_pop.append(row)
        initial_population = np.array(init_pop, dtype=float)

        # Print catalogues
        print("\nColumn catalogue (S700):", flush=True)
        print(f"  {'idx':>4}  {'d0':>6}  {'t0':>4}  {'d/t':>5}  {'area_proxy':>10}", flush=True)
        for i, (d, t) in enumerate(COL_CATALOGUE):
            print(f"  {i:>4}  {d:>6.1f}  {t:>4.1f}  {d/t:>5.1f}  {d**2-(d-2*t)**2:>10.1f}", flush=True)

        print("\nBrace catalogue (S235):", flush=True)
        print(f"  {'idx':>4}  {'d1':>6}  {'t1':>4}  {'d/t':>5}  {'area_proxy':>10}", flush=True)
        for i, (d, t) in enumerate(BRACE_CATALOGUE):
            print(f"  {i:>4}  {d:>6.1f}  {t:>4.1f}  {d/t:>5.1f}  {d**2-(d-2*t)**2:>10.1f}", flush=True)

        # Print initial population
        print("\nInitial population:", flush=True)
        print(f"  {'#':>3}  {'rad':>8}  {'col_idx':>7}  {'d0':>6}  {'t0':>4}  {'brace_idx':>9}  {'d1':>6}  {'t1':>4}", flush=True)
        for i, row in enumerate(initial_population):
            col_i = int(round(row[1]))
            brc_i = int(round(row[2]))
            d0, t0 = COL_CATALOGUE[col_i]
            d1, t1 = BRACE_CATALOGUE[brc_i]
            tag = "  <-- seed" if i == 0 else ""
            print(f"  {i:>3}  {row[0]:>8.2f}  {col_i:>7}  {d0:>6.1f}  {t0:>4.1f}  {brc_i:>9}  {d1:>6.1f}  {t1:>4.1f}{tag}", flush=True)
        print(flush=True)

        # ------------------------------------------------------------------
        # Stock PyGAD GA — only mandatory parameters + gene_space/type;
        # everything else is the library default.
        # ------------------------------------------------------------------
        ga_instance = pygad.GA(
            num_generations=GA_Settings["n_gen"],
            num_parents_mating=max(2, GA_Settings["pop_size"] // 2),
            sol_per_pop=GA_Settings["pop_size"],
            fitness_func=fitness_func,
            on_generation=on_generation,
            initial_population=initial_population,
            gene_space=GENE_SPACE,
            gene_type=GENE_TYPE,
            # --- stock defaults (explicit for clarity) --------------------
            parent_selection_type=GA_Settings["parent_selection_type"],
            keep_parents=GA_Settings["keep_parents"],
            crossover_type=GA_Settings["crossover_type"],
            mutation_type=GA_Settings["mutation_type"],
            mutation_percent_genes=GA_Settings["mutation_percent_genes"],
            random_seed=GA_Settings["seed"],
            suppress_warnings=True,
        )

        ga_instance.run()

    finally:
        mapdl.exit()

    # ---- Final report ---------------------------------------------------------
    best    = _state.best
    n_eval  = _state.n_eval

    logger.log_line("=" * 70)
    logger.log_line("GENETIC ALGORITHM (PyGAD stock) - FINAL RESULT")
    logger.log_line("=" * 70)

    if best is None:
        logger.log_line("No feasible design found.  Increase pop_size / n_gen or relax bounds.")
        success  = False
        best_x   = np.full(len(ENG_NAMES), np.nan)
        best_f   = float("nan")
        message  = "No feasible individual found"
    else:
        success  = True
        best_x   = best["x_eng"]
        best_f   = best["f"]
        message  = "Feasible optimum found by GA"
        logger.log_line(f"Best mass            : {best_f:.4f} kg")
        logger.log_line(f"Max utilisation      : {best['max_util']:.4f}")
        logger.log_line(f"First eigenvalue a_cr: {best['a_cr']:.4f}  (limit {GA_Settings['eigenvalue_min']})")
        logger.log_line(
            "Segment masses [kg]  : "
            + ", ".join(f"{m:.3f}" for m in best["segment_masses"])
        )
        logger.log_line("")
        logger.log_line("Optimal sections:")
        for name in ENG_NAMES:
            logger.log_line(f"  {name:4s} = {best['var'][name]['value']:.3f}")
        logger.log_line("")
        logger.log_line("Per-mode utilisation:")
        logger.log_utilization(best["util_report"])

    n_gen_done = int(getattr(ga_instance, "generations_completed", GA_Settings["n_gen"])) \
        if ga_instance is not None else 0

    result_like = types.SimpleNamespace(
        x=best_x,
        fun=best_f,
        nit=n_gen_done,
        nfev=n_eval,
        success=success,
        message=message,
        status=0 if success else -1,
    )
    logger.finalize(result_like)

    elapsed = time.time() - start
    print(f"\nGA completed in {elapsed:.2f} s", flush=True)
    print(f"Best objective : {best_f}")
    print(f"Message        : {message}")
    print(f"TXT log        : {logger.txt_path}")
    print(f"Objective CSV  : {logger.csv_path}")

    return result_like, logger.txt_path, logger.csv_path


if __name__ == "__main__":
    main()
