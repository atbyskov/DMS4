# Search_space.py
# This script creates quasi-random points for inital values to the design variables
# These points can then be transfered and run in main_loop.py or similar.
# Pseudocode for the script:
#   1. Write design variable names and bounds
#   2. Create Grid sampling
#   3. Create Random sampling
#   4. Create Latin Hypercube sampling
#   5. Create Sobol samplint
#   6. Calculate Disceptancy for each method and print
#   7. write Latin Hypercube sampling points to a .txt file: "seach_space.txt"

# Import Tools 
import numpy as np
import time
from functools import wraps
from scipy.stats import qmc
import matplotlib.pyplot as plt

# Write Bounds, must match exacly Main.py
var_bounds = {
    "d0":  (40.0,  100.0),
    "t0":  (1.0,     7.0),
    #"d1":  (10.0,  100.0),
    #"t1":  (0.1,     7.0),
    #"rad": (150.0,  350.0)
}

# Timer Decorator
def timer(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        tic = time.perf_counter()
        result = fn(*args, **kwargs)
        toc = time.perf_counter()
        print(f"Time for {fn.__name__}: {toc - tic:.6f} s")
        return result
    return wrapper

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
    sampler = qmc.LatinHypercube(d=d, scramble=scramble, seed=seed)
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

# Use non-equal spacing (remove =2 from n_per_dim=2 in)
n_per_dim = {
    "d0": 3,
    "t0": 3, 
    #"d1": 3,
    #"t1": 3,
    #"rad": 3,
}
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

# Print results
print("!---- Results ---- !")
print("Number of simulations :", grid.shape[0])
print(f"GRID: {grid_m:.6f}")
print(f"RAND: {rand_m:.6f}")
print(f"LHS:  {LHS_m:.6f}")
print(f"SOBOL: {Sobol_m:.6f}")

# Plotting
def plotting():

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
                [ymin, ymin, ymax, ymax, ymin])

        # Titles
        ax.set_title(f"{methods[i]}\nCD = {scores[i]:.4f}")

        # Labels
        ax.set_xlabel("d0 [mm]")
        ax.set_ylabel("t0 [mm]")

        # Limits (important for comparison)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        # Grid
        ax.grid(True)

    # Main title
    fig.suptitle(f"Search Methods ({n_per_dim_equal} values per dimension)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


plotting()