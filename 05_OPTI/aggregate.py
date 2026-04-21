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
        
        g = np.asarray(g, dtype = float)

        # No Aggregate Method 
        if self.method is None:
            print("Aggregate first value: ")
            print(g)
            return 1.0 - g
        
        # P-norm method
        if self.method == "P-norm":
            pass

