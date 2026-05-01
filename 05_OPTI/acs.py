# ACS.py
# Script for handling Adaptive Constraint Scaling Method
# Is only called per iteration

class ACSclass:
    def __init__(self, alpha=1.0, c0=1.0):
        self.alpha = alpha
        self.c = c0
        self.prev_v = None
        self.prev_gmax = None
        self.iter = 0


    def update(self,v,gmax):
        # Called once per major iteration
        if self.prev_v is not None and self.prev_v > 0.0:
            ratio = self.prev_gmax / self.prev_v
            self.c = self.alpha * ratio + (1.0 - self.alpha) * self.c

        self.prev_v = v
        self.prev_gmax = gmax
        self.iter += 1

        print(
            f"[ACS] iter= {self.iter:3d} |"
            f"g_max = {gmax:.4f} |"
            f"v = {v:.4f} |"
            f"c={self.c:.4f}"
        )