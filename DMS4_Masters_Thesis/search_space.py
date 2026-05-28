
"""
Search_space.py
----------------

This module generates quasi-random and structured sampling points for initializing
design variables in an optimization problem. It enables exploration of the design
space prior to or alongside optimization, improving robustness and reducing the
risk of local minima.

Sampling Methods:
   - Grid Sampling (Equal Spacing):
        Structured, deterministic grid across all dimensions.
   - Random Sampling:
        Uniform random distribution within bounds.
   - Latin Hypercube Sampling (LHS):
        Stratified sampling ensuring better space coverage.
   - Sobol Sampling:
        Quasi-random low-discrepancy sequence for optimal uniformity.


Outputs:
--------
- search_space.txt : Text file containing sampled design points
- Console output   : Discrepancy scores for all sampling methods
- Plots            : Visual comparison of sampling distributions

Notes:
------
- Bounds must match those used in the main optimization script.
- This module is intended for pre-processing and is not part of the main optimization loop.

"""

# Import Tools 
import numpy as np
import time
from functools import wraps
from scipy.stats import qmc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

# Write Bounds, must match exacly Main.py
var_bounds = {
    "rad": (150.0,  350.0),
    "d0":  (48.3,  108),
    "t0":  (2.5,   5),
    "d1":  (25,    60),
    "t1":  (2,     6)
}


# Equal Spacing
def equal_spacing(var_bounds,n_per_dim):
    
    names = list(var_bounds.keys())

    axes = []
    for name in names:
        # Read lower and upper bounds
        lower,upper = var_bounds[name]

        if isinstance(n_per_dim, dict):
            num = int(n_per_dim[name])
        else:
            num = int(n_per_dim)

        # Create values and append
        values = np.linspace(lower, upper, num=num)
        axes.append(values)

    # Create mesh
    mesh = np.meshgrid(*axes, indexing="ij")

    grid = np.stack(mesh, axis=-1).reshape(-1, len(names))

    return grid, names

# Random Spacing
def random_spacing(var_bounds, n_points, seed=None):
    """
    Random points uniformly inside bounds.
    Output format matches equal_spacing: (points, names)
    """
    rng = np.random.default_rng(seed)

    names = list(var_bounds.keys())
    lowers = np.array([var_bounds[name][0] for name in names], dtype=float)
    uppers = np.array([var_bounds[name][1] for name in names], dtype=float)

    u = rng.random((n_points, len(names)))   # in [0,1)
    points = lowers + u * (uppers - lowers)  # scale to [lower, upper]

    return points, names

# LHS Spacing
def latin_spacing(var_bounds, n_points, seed=None, scramble=True):
    names = list(var_bounds.keys())
    d = len(names)

    # LHS in the unit hypercube [0,1]^d
    sampler = qmc.LatinHypercube(d=d, scramble=scramble, seed=seed,optimization="random-cd")
    sample_unit = sampler.random(n=n_points)

    # Scale to physical bounds
    lower = np.array([var_bounds[name][0] for name in names], dtype=float)
    upper = np.array([var_bounds[name][1] for name in names], dtype=float)
    points = qmc.scale(sample_unit, lower, upper)

    return points, names

# Sobol Spacing
def sobol_spacing(var_bounds, n_points, seed=None, scramble=True, power2=True):
    
    names = list(var_bounds.keys())
    d = len(names)

    sampler = qmc.Sobol(d=d, scramble=scramble, seed=seed)

    if power2:
        # Find smallest 2^m >= n_points
        m = int(np.ceil(np.log2(n_points)))
        sample_unit = sampler.random_base2(m=m)[:n_points]
    else:
        sample_unit = sampler.random(n=n_points)

    lower = np.array([var_bounds[name][0] for name in names], dtype=float)
    upper = np.array([var_bounds[name][1] for name in names], dtype=float)
    points = qmc.scale(sample_unit, lower, upper)

    return points, names


n_per_dim_equal = 2

# Call grid function with options
grid, names = equal_spacing(var_bounds, n_per_dim=n_per_dim_equal)
N = grid.shape[0]

# Call random samplint (Seed is set for reproducability)
rand_points,_ = random_spacing(var_bounds, n_points=N, seed=42)

# Use LHS 
lhs_points, _ = latin_spacing(var_bounds, n_points=N, seed=42, scramble=True)

# Use Sobol
sobol_points, _ = sobol_spacing(var_bounds, n_points=N, seed=42, scramble=True, power2=True)


# Evaluate Discrepancy Function
def evaluate_spacing(points, var_bounds, names, method="CD"):
    """
    Scale points to [0,1]^d using var_bounds and return discrepancy.
    """
    lowers = np.array([var_bounds[n][0] for n in names], dtype=float)
    uppers = np.array([var_bounds[n][1] for n in names], dtype=float)

    points_unit = (points - lowers) / (uppers - lowers)
    return qmc.discrepancy(points_unit, method=method)


# Call evaluations
grid_m = evaluate_spacing(grid,          var_bounds, names)
rand_m = evaluate_spacing(rand_points,   var_bounds, names)
LHS_m = evaluate_spacing(lhs_points,     var_bounds, names)
Sobol_m = evaluate_spacing(sobol_points, var_bounds, names)

data = {"GRID": {"points": grid, "names": names, "discrepancy": grid_m},
        "RAND": {"points": rand_points, "names": names, "discrepancy": rand_m},
        "LHS" : {"points": lhs_points, "names": names, "discrepancy": LHS_m},
        "SOBOL":{"points": sobol_points, "names": names, "discrepancy": Sobol_m}

}

# Print results
print("!---- Results ---- !")
print("Number of simulations :", grid.shape[0])
print(f"GRID: {grid_m:.6f}")
print(f"RAND: {rand_m:.6f}")
print(f"LHS:  {LHS_m:.6f}")
print(f"SOBOL: {Sobol_m:.6f}")


method_to_print = "LHS"
# Write to file
def read_to_file(var_bounds, data, method_key, n_per_dim_equal):
    points = data[method_key]["points"]
    names = data[method_key]["names"]
    DS = data[method_key]["discrepancy"]
    n_points = points.shape[0]

    with open("search_space.txt","w") as f:
        f.write("! Search Space Points \n")
        f.write(f"Method: {method_key} \n")
        f.write(f"Points per Dimension: {n_per_dim_equal} \n")
        f.write(f"Number of Simulations: {n_points} \n")
        f.write(f"discrepancy score: {DS:.4f} \n \n")
        f.write(f"Active Design Variables and Bounds: \n")
        for name, (lower, upper) in var_bounds.items():
            f.write(f" {name}: [{lower}, {upper}] \n")

        f.write("\nPoints: \n")
        f.write("  "+ "  ".join(names) + "\n")
        for row in points:
            f.write("  "+"  ".join(f"{val:.4f}" for val in row) + "\n")


read_to_file(var_bounds,data,method_to_print,n_per_dim_equal)


# Plotting 2D
def plotting2D():

    fig, axs = plt.subplots(1, 4, figsize=(14, 4))

    # Data + labels
    methods = ["Grid", "Random", "LHS", "Sobol"]
    scores  = [grid_m, rand_m, LHS_m, Sobol_m]
    data    = [grid, rand_points, lhs_points, sobol_points]

    # Bounds
    xmin, xmax = var_bounds["d0"]
    ymin, ymax = var_bounds["t0"]

    for i, ax in enumerate(axs):

        # Scatter
        ax.scatter(data[i][:,0], data[i][:,1])

        # Bounding box
        ax.plot([xmin, xmax, xmax, xmin, xmin],
                [ymin, ymin, ymax, ymax, ymin], color="black")

        # Titles
        ax.set_title(f"{methods[i]}\nCD = {scores[i]:.4f}")

        # Labels
        ax.set_xlabel(r"$d_0$ [mm]")
        ax.set_ylabel(r"$t_0$ [mm]")

        # Limits (important for comparison)
        ax.set_xlim(xmin-5, xmax+5)
        ax.set_ylim(ymin-1, ymax+1)

        # Grid
        ax.grid(True)

    # Main title
    fig.suptitle(f"Search Methods ({n_per_dim_equal} values per dimension)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()



plotting2D()

# Aggregate Output
