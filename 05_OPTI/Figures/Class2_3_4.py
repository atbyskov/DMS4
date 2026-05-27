import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 1200
matplotlib.rcParams['text.antialiased'] = True
matplotlib.rcParams['lines.antialiased'] = True
matplotlib.rcParams['patch.antialiased'] = True
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['CMU Serif', 'Computer Modern Roman', 'Times New Roman', 'DejaVu Serif']
matplotlib.rcParams['axes.linewidth'] = 1.2
matplotlib.rcParams['xtick.major.width'] = 1.0
matplotlib.rcParams['ytick.major.width'] = 1.0
matplotlib.rcParams['xtick.minor.width'] = 0.6
matplotlib.rcParams['ytick.minor.width'] = 0.6
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# ---------------------------------------------------------
# 1. Production Profile Dataset (D, t) — S700 Columns
# ---------------------------------------------------------

s700_data = [
    (48.3, 2.5), (48.3, 3.0), (60.3, 2.5), (60.3, 3.0), (60.3, 4.0),
    (76.1, 3.0), (76.1, 4.0), (76.1, 5.0), (88.9, 3.0), (88.9, 4.0),
    (88.9, 5.0), (101.6, 3.0), (101.6, 4.0), (101.6, 5.0), (108.0, 3.0),
    (108.0, 4.0), (108.0, 5.0), (114.3, 3.0), (114.3, 4.0), (114.3, 5.0),
    (121.0, 4.0), (121.0, 5.0), (127.0, 3.0), (127.0, 4.0), (127.0, 5.0),
    (139.7, 4.0), (139.7, 5.0), (139.7, 6.0), (139.7, 8.0), (139.7, 10.0),
    (168.3, 4.0), (168.3, 5.0), (168.3, 6.0), (168.3, 8.0), (168.3, 10.0),
    (219.1, 5.0), (219.1, 6.0), (219.1, 8.0), (219.1, 10.0), (273.0, 5.0),
    (273.0, 6.0), (273.0, 8.0), (273.0, 10.0), (323.9, 6.0), (323.9, 8.0),
    (323.9, 10.0)
]

# ---------------------------------------------------------
# 2. EN 1993-1-1 Table 5.2 — CHS Class Limits
# ---------------------------------------------------------

fy = 700
eps_sq = 235 / fy  # ε² = 235 / f_y

limit_class1 = 50 * eps_sq   # d/t ≤ 50ε²  → 16.786
limit_class2 = 70 * eps_sq   # d/t ≤ 70ε²  → 23.500
limit_class3 = 90 * eps_sq   # d/t ≤ 90ε²  → 30.214

# Specified project design bounds
bound_col_t = (2.5, 5.0)
bound_col_d = (48.3, 114.3)

# ---------------------------------------------------------
# 3. Classify Each Production Profile
# ---------------------------------------------------------

def classify(d, t):
    ratio = d / t
    if ratio <= limit_class1:
        return 1
    elif ratio <= limit_class2:
        return 2
    elif ratio <= limit_class3:
        return 3
    else:
        return 4

class_points = {1: ([], []), 2: ([], []), 3: ([], []), 4: ([], [])}
for d_val, t_val in s700_data:
    cls = classify(d_val, t_val)
    class_points[cls][0].append(t_val)
    class_points[cls][1].append(d_val)

# ---------------------------------------------------------
# 4. Colour Palette
# ---------------------------------------------------------

region_colors = {
    1: '#1b7a3d',   # forest green
    2: '#2e86c1',   # steel blue
    3: '#d4a017',   # dark gold
    4: '#c0392b',   # crimson
}

region_alphas = {1: 0.14, 2: 0.14, 3: 0.14, 4: 0.14}

scatter_face = {
    1: '#27ae60',
    2: '#3498db',
    3: '#f1c40f',
    4: '#e74c3c',
}

scatter_edge = {
    1: '#145a2a',
    2: '#1a5276',
    3: '#7d6608',
    4: '#78281f',
}

boundary_colors = {
    1: '#145a2a',
    2: '#1a5276',
    3: '#7d6608',
}

# ---------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------

max_t = 11
max_d = 340

fig, ax = plt.subplots(figsize=(11, 8.5))

t_space = np.linspace(0.01, max_t + 5, 5000)
d_class1 = limit_class1 * t_space
d_class2 = limit_class2 * t_space
d_class3 = limit_class3 * t_space

# Region fills (bottom to top)
ax.fill_between(t_space, 0, d_class1,
                color=region_colors[1], alpha=region_alphas[1], label='Class 1')
ax.fill_between(t_space, d_class1, d_class2,
                color=region_colors[2], alpha=region_alphas[2], label='Class 2')
ax.fill_between(t_space, d_class2, d_class3,
                color=region_colors[3], alpha=region_alphas[3], label='Class 3')
ax.fill_between(t_space, d_class3, max_d + 150,
                color=region_colors[4], alpha=region_alphas[4], label='Class 4')

# Boundary lines
ax.plot(t_space, d_class1, color=boundary_colors[1], linestyle='-', linewidth=2.5,
        label=rf'$d/t = 50\varepsilon^2 = {limit_class1:.1f}$', antialiased=True)
ax.plot(t_space, d_class2, color=boundary_colors[2], linestyle='--', linewidth=2.5,
        label=rf'$d/t = 70\varepsilon^2 = {limit_class2:.1f}$', antialiased=True)
ax.plot(t_space, d_class3, color=boundary_colors[3], linestyle='-.', linewidth=2.5,
        label=rf'$d/t = 90\varepsilon^2 = {limit_class3:.1f}$', antialiased=True)

# Design bounds box
rect = patches.Rectangle(
    (bound_col_t[0], bound_col_d[0]),
    bound_col_t[1] - bound_col_t[0],
    bound_col_d[1] - bound_col_d[0],
    linewidth=2.5, edgecolor='#2c3e50', facecolor='none',
    linestyle='--', label='Design Bounds', zorder=6
)
ax.add_patch(rect)

# Scatter production profiles per class
for cls in [1, 2, 3, 4]:
    t_pts, d_pts = class_points[cls]
    if t_pts:
        ax.scatter(t_pts, d_pts,
                   color=scatter_face[cls], edgecolors=scatter_edge[cls],
                   s=65, linewidths=1.2, zorder=7, label=f'Class {cls} Profile')

# Formatting
ax.set_xlim(0, max_t)
ax.set_ylim(0, max_d)
ax.set_xlabel(r'Wall Thickness, $t_0$ [mm]', fontsize=17, labelpad=10)
ax.set_ylabel(r'Outer Diameter, $d_0$ [mm]', fontsize=17, labelpad=10)
ax.set_title(
    r'Strenx 700 Columns — Cross-Section Classification (EN 1993-1-1)'
    '\n'
    rf'$\varepsilon^2 = 235\,/\,f_y = {eps_sq:.4f}$',
    fontsize=19, pad=18
)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=11)
ax.grid(True, linestyle=':', alpha=0.5, linewidth=0.7)
ax.legend(loc='upper left', fontsize=13, framealpha=0.95, fancybox=True,
          edgecolor='#555555', shadow=True)

fig.tight_layout()

save_dir = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(save_dir, 'Class2_3_4.png'),
            dpi=1200, bbox_inches='tight', pad_inches=0.15,
            facecolor='white', edgecolor='none')
fig.savefig(os.path.join(save_dir, 'Class2_3_4.pdf'),
            bbox_inches='tight', pad_inches=0.15,
            facecolor='white', edgecolor='none')
fig.savefig(os.path.join(save_dir, 'Class2_3_4.svg'),
            bbox_inches='tight', pad_inches=0.15,
            facecolor='white', edgecolor='none')

plt.show()
