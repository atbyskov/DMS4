# ACS.py
# Adaptive Constraint Scaling (ACS) following:
#   Le et al. (2010) - Stress-based topology optimization for continua, Eqs. 9-10
#   Oest & Lund (2017) - Topology optimization with finite-life fatigue constraints, Fig. 2
#
# The scaling factor c corrects the over-conservative P-norm aggregate so that
# c * sigma_PN -> sigma_max as iterations stabilise.  c is bounded in (0, 1]
# because sigma_PN >= sigma_max for nonnegative values with p >= 1.


class ACSclass:
    def __init__(self, enabled=True, name=""):
        self.enabled = enabled
        self.name = name

        # State (matches T.c_Le / T.alpha_Le / T.c_Le_old / T.c_Le_oldold in ACS.m)
        self.alpha = 1.0
        self.c = 1.0          # paper starts with no scaling; learns from first iter
        self.c_old = None
        self.c_oldold = None

        # Iteration counter (matches IterNo in ACS.m, 0-based)
        self.iter = 0

    def update(self, sigma_PN, sigma_max):
        """Update c using the most recent iteration's (sigma_PN, sigma_max).

        Both inputs must come from the *previous* iterate's aggregate so the
        next forward call (which uses self.c) is consistent with Eq. 10:
            c^I = alpha^I * (sigma_max^{I-1} / sigma_PN^{I-1})
                  + (1 - alpha^I) * c^{I-1}
        """
        # ---- ACS OFF ----
        if not self.enabled:
            self.c = 1.0
            return self.c

        # Safety floor: sigma_PN can never be exactly 0 if there is any nonzero
        # utilization, but clamp anyway for robustness against the all-zero corner.
        sigma_PN = max(float(sigma_PN), 1e-12)
        sigma_max = max(float(sigma_max), 0.0)
        ratio = sigma_max / sigma_PN  # in (0, 1] for valid P-norm aggregates

        if self.iter == 0:
            # ACS.m line 22-27
            self.c = ratio
            self.alpha = 1.0

        elif self.iter == 1:
            # ACS.m line 30-37
            self.c_old = self.c
            self.c = ratio

        else:
            # ACS.m line 40-88
            self.c_oldold = self.c_old
            self.c_old = self.c

            scale_up = 1.20
            scale_down = 0.80
            decimal = 2

            # Tentative c with the *current* alpha (used for oscillation check only)
            c_trial = self.alpha * ratio + (1.0 - self.alpha) * self.c_old

            c_oldold_r = round(self.c_oldold, decimal)
            c_old_r = round(self.c_old, decimal)
            c_trial_r = round(c_trial, decimal)

            oscillating = (
                (c_oldold_r < c_old_r and c_old_r > c_trial_r) or
                (c_oldold_r > c_old_r and c_old_r < c_trial_r)
            )

            if oscillating:
                self.alpha = max(0.5, self.alpha * scale_down)
            else:
                self.alpha = min(1.0, self.alpha * scale_up)

            # Recompute c with possibly updated alpha (ACS.m line 73)
            self.c = self.alpha * ratio + (1.0 - self.alpha) * self.c_old

        self.iter += 1

        print(
            f"[ACS:{self.name:>8s}] iter={self.iter:3d} | "
            f"sigma_max={sigma_max:.4f} | sigma_PN={sigma_PN:.4f} | "
            f"ratio={ratio:.4f} | alpha={self.alpha:.3f} | c={self.c:.4f}",
            flush=True,
        )

        return self.c
