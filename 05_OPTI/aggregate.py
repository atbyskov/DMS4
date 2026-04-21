# Aggregate.py
# This script handles the aggregation of utilization constraints for optimization
# ----------------------------------------------------
# INPUT  <- Utilization Constraints as list (c_util_agg)
# OUTPUT -> Aggregated utilization Constraints as list
# ----------------------------------------------------
import numpy as np

class ConstraintAggregate:
    def __init__(self, method=None, p_value = None):
        # Read Method
        self.method = method
        print(method)
        # Read P value
        self.p = p_value

        

    
    def agg_output(self,g):
        eps = 0
        
        g = np.asarray(g, dtype = float)

        # No Aggregate Method 
        if self.method is None:
            return 1.0 - g
        
        # P-norm method
        if self.method == "P-norm":
            p = float(self.p)
            con_constraints = 1 - g - eps
            
            return (np.sum((con_constraints)**p))**(1 / p)
        
        # P-norm mean method
        if self.method == "P-norm-mean":
            p = float(self.p)
            con_constraints = 1 - g - eps
            n0 = con_constraints.size

            
            return (np.sum((con_constraints)**p) / n0)**(1 / p)
        
        if self.method == "KS":

            return "hey ho lets go"

