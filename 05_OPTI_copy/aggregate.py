
"""
Aggregate.py
------------

This module provides functionality for aggregating multiple constraint values
into a single representative constraint for use in gradient-based optimization.

Implemented Methods:
-------------------
1. None (No Aggregation):
   - Returns all constraints directly (no reduction).

2. P-norm:
   - Smooth approximation of the maximum violated constraint.
   - Defined as:
         v = (Σ g_k^p)^(1/p)
   - As p → ∞, approaches max(g_k).

3. P-norm-mean:
   - Normalized version of P-norm:
         v = (1/n * Σ g_k^p)^(1/p)
   - Reduces sensitivity to number of constraints.

Inputs:
-------
- g         : Array-like constraint vector (g_i ≥ 0 formulation)
- method    : Aggregation method ("P-norm", "P-norm-mean", or None)
- p_value   : Order of P-norm (controls smoothness vs. max approximation)
- rho_value : (Reserved for KS-type methods, not currently implemented)
- c_scale   : External scaling factor (typically from ACS)

Outputs:
--------
- c_agg     : Aggregated constraint value (scalar)
- v         : Raw aggregated violation measure
- g_max     : Maximum individual violation

"""

import numpy as np

# Define the class
class ConstraintAggregate:
    def __init__(self, method=None, p_value = None, rho_value = None, relaxation = None):
        # Read Method
        self.method = method
        self.rho = rho_value

        # Read P value
        self.p = p_value

    # Function for handling aggregation
    def agg_output(self,g, c_scale = 1.0):
        # Relaxation parameter for P-norm methods (0<eps<1)
        relaxation = 0
        g = np.asarray(g, dtype = float)                            # Read constraints
        g_k = np.maximum(-(g - relaxation), 0.0)                    # Only violated
        g_max = np.max(g_k)                                         # Maximum

        # No Aggregate Method 
        if self.method is None:
            return g, None, None
        
        # P-norm method
        if self.method == "P-norm":
            p = float(self.p)                                       # Read p-value
            v = (np.sum(g_k**p))**(1/p)                             # Formula
            print(f"Max Const: {g_max:.3f} | P-norm: {v:.3}")       # Print

        # P-norm mean method
        if self.method == "P-norm-mean":
            p = float(self.p)
            n = np.size(g)
            g_k = np.maximum(-(g - relaxation), 0.0)
            v = (1/n*np.sum(g_k**p))**(1/p)


            print(f"Max Const: {g_max:.3f} | P-norm-mean: {v:.3}")

    
        
        v_scaled = c_scale * v
        
        return -v_scaled, v, g_max
  

