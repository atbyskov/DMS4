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
from scipy.stats import qmc
import matplotlib.pyplot as plt

var_bounds = {
    "d0":  (40.0,  100.0),
    "t0":  (1.0,     7.0),
    "d1":  (10.0,  100.0),
    #"t1":  (0.1,     7.0),
    #"rad": (150.0,  350.0)
}

n = 30      # Sampling size
d = len(var_bounds)


# Variables and Bounds
def equal_spacing(var_bounds,step):
    
    names = list(var_bounds.keys())

    axes = []
    for name in names:
        # Read lower and upper bounds
        lower,upper = var_bounds[name]

        # Step size
        step_size = float(step[name])

        # Create values and append
        values = np.arange(lower,upper + 0.5*step_size, step_size)
        axes.append(values)

    # Create mesh
    mesh = np.meshgrid(*axes, indexing="ij")

    grid = np.stack(mesh, axis=-1).reshape(-1, len(names))

    return grid, names



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


# Use equal spacing
step = {
    "d0": 10,
    "t0": 1, 
    "d1": 10,
    "t1": 1,
    "rad": 10,
}
# Use grid
grid, names = equal_spacing(var_bounds, step)

# Use random
N = grid.shape[0]
rand_points,_ = random_spacing(var_bounds, n_points=N, seed=42)


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

#grid_2dplot()
rand_2dplot()
grid_2dplot()

# We scale
def scale_to_unit(points, var_bounds, names):
    lowers = np.array([var_bounds[n][0] for n in names])
    uppers = np.array([var_bounds[n][1] for n in names])

