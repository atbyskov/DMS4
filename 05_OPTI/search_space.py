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

var_bounds = {
    "d0":  (40.0,  100.0),
    "t0":  (1.0,     7.0),
    #"d1":  (10.0,  100.0),
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

# Use equal spacing
step = {
    "d0": 1,
    "t0": 1, 
    "d1": 10,
    "t1": 1,
    "rad": 10,
}
grid, names = equal_spacing(var_bounds, step)


print("Variable order:", names)
print("Grid shape:", grid.shape)
print("First 5 points:\n", grid[:5])
print("Last 5 points:\n", grid[-5:])



