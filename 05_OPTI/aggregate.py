# Aggregate.py
# P-norm / P-norm-mean aggregation of constraint *values* (e.g. utilization
# ratios), paired with the ACS scaling factor c from acs.py.
#
# Follows Le et al. (2010) Eq. 8-9 and Oest & Lund (2017) Eq. 15:
#     sigma_PN(x)  = ( sum  u_i(x)^p )^(1/p)         (or (1/n * sum)^(1/p) for mean)
#     sigma_max(x) = max  u_i(x)
#     aggregated constraint:  c * sigma_PN(x) <= 1
#                  in c(x) >= 0 form:   1 - c * sigma_PN(x) >= 0
#
# This is *NOT* a P-norm of violations -- it is a P-norm of the actual
# utilization values themselves.  That is what makes ACS well defined in the
# feasible region (the ratio sigma_max/sigma_PN stays bounded away from 0/0).
# ----------------------------------------------------
# INPUT  <- raw utilization values u (non-negative array)
#        <- method, p_value chosen in Main_Self_Written_Test.py
#        <- c_scale: the ACS scaling factor (default 1.0)
# OUTPUT -> c_row     : scalar 1 - c_scale*sigma_PN (or per-element 1-u when method=None)
#        -> sigma_PN  : the unscaled P-norm value
#        -> sigma_max : the true max of u
# ----------------------------------------------------
import numpy as np


class ConstraintAggregate:
    def __init__(self, method=None, p_value=None, rho_value=None, relaxation=None):
        self.method = method
        self.p = p_value
        self.rho = rho_value
        # relaxation kept for backward compatibility but not used here.

    def agg_output(self, u, c_scale=1.0):
        u = np.asarray(u, dtype=float).ravel()
        # Negative utilization values are not physical; clamp at zero so the
        # P-norm is well defined and faithful to the Le et al. formulation.
        u_pos = np.maximum(u, 0.0)
        sigma_max = float(np.max(u_pos)) if u_pos.size else 0.0

        # ---- No aggregation: pass through per-element 1 - u in c(x)>=0 form ----
        if self.method is None or self.method == "None":
            return 1.0 - u_pos, None, sigma_max

        # ---- P-norm of values ----
        if self.method == "P-norm":
            p = float(self.p)
            sigma_PN = float(np.sum(u_pos ** p) ** (1.0 / p))

        elif self.method == "P-norm-mean":
            p = float(self.p)
            n = max(int(u_pos.size), 1)
            sigma_PN = float((np.sum(u_pos ** p) / n) ** (1.0 / p))

        else:
            raise ValueError(f"Unknown aggregation method: {self.method!r}")

        # Aggregated constraint in c(x) >= 0 form:  1 - c_scale * sigma_PN >= 0
        c_row = 1.0 - float(c_scale) * sigma_PN
        return c_row, sigma_PN, sigma_max
