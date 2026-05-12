# ACS.py

class ACSclass:
    def __init__(self, alpha=1.0, c0=1.0):
        self.alpha = float(alpha)
        self.c = float(c0)

        self.c_old = None
        self.c_oldold = None

        self.iter = 0

    def update(self, g, gmax):
        """
        g     : aggregated constraint (scalar)
        gmax  : max constraint before aggregation (scalar)
        """

        # Safety
        g = max(float(g), 1e-12)
        gmax = float(gmax)

        # ---- Iteration 0 ----
        if self.iter == 0:
            self.c = gmax / g
            self.alpha = 1.0

        # ---- Iteration 1 ----
        elif self.iter == 1:
            self.c_old = self.c
            self.c = gmax / g

        # ---- Iteration > 1 ----
        else:
            self.c_oldold = self.c_old
            self.c_old = self.c

            scale_up = 1.20
            scale_down = 0.80
            decimal = 2

            ratio = gmax / g

            # Temporary c
            c_new = self.alpha * ratio + (1 - self.alpha) * self.c_old

            # Rounded for oscillation detection
            c_oldold_r = round(self.c_oldold, decimal)
            c_old_r = round(self.c_old, decimal)
            c_new_r = round(c_new, decimal)

            oscillating = (
                (c_oldold_r < c_old_r and c_old_r > c_new_r) or
                (c_oldold_r > c_old_r and c_old_r < c_new_r)
            )

            if oscillating:
                self.alpha = max(0.5, self.alpha * scale_down)
            else:
                self.alpha = min(1.0, self.alpha * scale_up)

            # Final update
            self.c = self.alpha * ratio + (1 - self.alpha) * self.c_old

        self.iter += 1

        # Debug print
        print(
            f"[ACS] iter={self.iter:3d} | "
            f"g={g:.4f} | gmax={gmax:.4f} | "
            f"alpha={self.alpha:.3f} | c={self.c:.4f}"
        )

        return self.c