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

# Random
def random(n,d,seed):
    rng = np.random.default_rng(seed)
    return rng.random((n,d))

ran=random(n,d,47)

print(random)

