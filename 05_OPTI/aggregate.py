# Aggregate.py
# This script handles the aggregation of utilization constraints for optimization
# ----------------------------------------------------
# INPUT  <- Utilization Constraints as list (c_util_agg)
#        <- Method: is chosen in Main
#        <- P-value: is chosen in Main
# OUTPUT -> Aggregated Utilization Constraints 
# ----------------------------------------------------
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
    def agg_output(self,g):
        # Relaxation parameter for P-norm methods (0<eps<1)
        relaxation = 0
        # Read constraints
        g = np.asarray(g, dtype = float)

        # No Aggregate Method 
        if self.method is None:
            return g
        
        # P-norm method
        if self.method == "P-norm":
            p = float(self.p)
            
            g_k = np.maximum(-(g - relaxation), 0.0)
            v_agg = (np.sum(g_k**p))**(1/p)

            print(f"P-norm: {v_agg:.2f}")

            return -v_agg
        
        # P-norm mean method
        if self.method == "P-norm-mean":
            p = float(self.p)
            con_constraints = g - relaxation
            n0 = con_constraints.size

            return (np.sum((con_constraints)**p) / n0)**(1 / p)
        
        if self.method == "KS":
            # rho parameter (typically between 2 and 200)
            rho = float(self.rho)
            g = g - relaxation
            g_max = np.max(g)

            return ((1/rho)*np.log(np.sum(np.exp(rho*g))))
        
        if self.method == "KS_shift":
            # rho parameter (typically between 2 and 200)
            rho = float(self.rho)
            g = 1.0 - g - relaxation
            g_max = np.max(g)

            return (g_max + (1/rho)*np.log(np.sum(np.exp(rho * (g-g_max)))))

