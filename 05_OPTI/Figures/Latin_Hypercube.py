import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
n = 3
np.random.seed(2)

# --- Generate Latin Hypercube ---
perm_x = np.random.permutation(n)
perm_y = np.random.permutation(n)

x = (perm_x + np.random.rand(n)) / n
y = (perm_y + np.random.rand(n)) / n

# --- Plot ---
fig, ax = plt.subplots(figsize=(6, 6))

# Draw grid
for i in range(n + 1):
    ax.axhline(i / n, color="gray", lw=0.8, alpha=0.6)
    ax.axvline(i / n, color="gray", lw=0.8, alpha=0.6)

# Plot points (blue dots)
ax.scatter(x, y, color="blue", s=40, zorder=3)

# Add grid coordinate labels
for i in range(n):
    gx = perm_x[i] + 1
    gy = perm_y[i] + 1

    cx = (perm_x[i] + 0.5) / n
    cy = (perm_y[i] + 0.5) / n

    ax.text(
        cx,
        cy,
        f"({gx},{gy})",
        ha="center",
        va="center",
        fontsize=11
    )

# Formatting
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")

# Bigger axis labels
ax.set_xlabel(r"$x_1$", fontsize=16)
ax.set_ylabel(r"$x_2$", fontsize=16)

ax.set_title("Latin Hypercube Sampling", fontsize=14)

plt.tight_layout()
plt.show()