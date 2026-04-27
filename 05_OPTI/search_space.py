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

import numpy as np
import time
from functools import wraps
from scipy.stats import qmc
import matplotlib.pyplot as plt

var_bounds = {
    "d0":  (40.0,  100.0),
    "t0":  (1.0,     7.0),
    "d1":  (10.0,  100.0),
    "t1":  (0.1,     7.0),
    "rad": (150.0,  350.0)
}

n = 30      # Sampling size
d = len(var_bounds)


# Timer Function
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

        # Step size
        num = int(n_per_dim[name])

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


# Use equal spacing
n_per_dim = {
    "d0": 3,
    "t0": 3, 
    "d1": 3,
    "t1": 3,
    "rad": 3,
}

grid, names = equal_spacing(var_bounds, n_per_dim)

# Use random
N = grid.shape[0]
rand_points,_ = random_spacing(var_bounds, n_points=N, seed=42)

# Use LHS
lhs_points, _ = latin_spacing(var_bounds, n_points=N, seed=42, scramble=True)


## Plotting
# def plot equal spacing

def grid_3dplot():
    x = grid[:, 0]
    y = grid[:, 1]
    z = grid[:, 2]

    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(x,y,z,s=20)

    ax.set_xlabel(names[0])
    ax.set_ylabel(names[1])
    ax.set_zlabel(names[2])
    ax.set_title("3D Grid Points")

    plt.grid(True)
    plt.show()

def grid_2dplot():

    plt.figure(figsize=(8,5))
    plt.scatter(grid[:, 0], grid[:, 1], s=20)
    plt.xlabel(names[0])
    plt.ylabel(names[1])
    plt.title("2D Grid points")
    plt.grid(True)
    plt.show()

# def plot random spacing
def rand_2dplot():
    plt.figure(figsize=(8,5))
    plt.scatter(rand_points[:,0],rand_points[:, 1], s=20)
    plt.xlabel(names[0])
    plt.ylabel(names[1])
    plt.title("2D Random Points")
    plt.grid(True)
    plt.show()

# Evaluate Discrepancy
def evaluate_spacing(points, var_bounds, names, method="CD"):
    """
    Scale points to [0,1]^d using var_bounds and return discrepancy.
    """
    lowers = np.array([var_bounds[n][0] for n in names], dtype=float)
    uppers = np.array([var_bounds[n][1] for n in names], dtype=float)

    points_unit = (points - lowers) / (uppers - lowers)
    return qmc.discrepancy(points_unit, method=method)


grid_m = evaluate_spacing(grid, var_bounds, names)
rand_m = evaluate_spacing(rand_points, var_bounds, names)
LHS_m = evaluate_spacing(lhs_points, var_bounds, names)


print("!---- Results ---- !")
print(f"GRID: {grid_m:.4f}")
print(f"RAND: {rand_m:.4f}")
print(f"LHS:  {LHS_m:.4f}")
print("Number of simulations :", grid.shape[0])

#grid_2dplot()
#rand_2dplot()

