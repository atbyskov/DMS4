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
# 1. Production Profile Datasets (D, t)
# ---------------------------------------------------------

# High-Strength Steel Columns (S700)
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

# Structural Steel Braces (S355)
s355_data = [
    (25, 2), (25, 3), (25, 4), (25, 5), (25, 6),
    (30, 2), (30, 3), (30, 4), (30, 5), (30, 6),
    (32, 5), (38, 4),
    (40, 2), (40, 3), (40, 4), (40, 5), (40, 6), (40, 7), (40, 8),
    (50, 2), (50, 3), (50, 4), (50, 5), (50, 6), (50, 9),
    (55, 6),
    (60, 2), (60, 3), (60, 4), (60, 5), (60, 6)
]

# Mild Steel Braces (S235) - Defined raw and flattened
s235_raw = [
    (6, [1, 1.5, 2]),
    (8, [1, 1.5, 2]),
    (10, [1, 1.5, 2, 2.5]),
    (12, [1, 1.5, 2, 2.5]),
    (13, [1, 1.5, 2, 2.5]),
    (14, [1, 1.5, 2, 2.5, 3, 3.5]),
    (15, [1, 1.5, 2, 2.5, 3]),
    (16, [1, 1.5, 2, 2.5, 3, 3.5, 4]),
    (17, [1, 1.5, 2, 2.5]),
    (18, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]),
    (20, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]),
    (21, [1, 1.5, 2, 2.5, 3, 4]),
    (22, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]),
    (23, [1, 1.5, 2]),
    (24, [1, 1.5, 2, 3, 4]),
    (25, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7]),
    (28, [1, 1.5, 2, 2.5, 3, 4, 5, 6]),
    (30, [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 7.5]),
    (32, [1.5, 2, 3, 3.5, 4, 5, 6]),
    (34, [1.5, 2, 4]),
    (35, [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 7.5, 8, 9]),
    (36, [1.5, 2, 3, 5, 8]),
    (38, [1.5, 2, 3, 4, 4.5, 5, 6, 7, 8]),
    (40, [1.5, 2, 2.5, 3, 4, 4.5, 5, 5.5, 6, 7, 7.5, 8, 10]),
    (42, [1.5, 2, 3, 4, 5, 7]),
    (45, [1.5, 2, 2.5, 3, 4, 4.5, 5, 6, 7, 10]),
    (50, [1.5, 2, 3, 4, 5, 6, 7, 8, 10]),
    (52, [5]),
    (55, [1.5, 2, 3, 4, 5, 6, 8, 10]),
    (60, [1.5, 2, 3, 4, 5, 6, 8, 9, 10]),
    (62, [6]),
    (65, [4, 6, 10]),
    (70, [2, 3, 4, 5, 6, 8, 10]),
    (80, [2, 3, 4, 5, 6, 8, 10, 12]),
    (90, [5]),
    (100, [3, 4, 5, 6])
]

s235_data = [(d, t) for d, t_list in s235_raw for t in t_list]

# ---------------------------------------------------------
# 2. Material Limits & Design Bounds
# ---------------------------------------------------------

# EN 1993-1-1 Class 2 limits (d/t <= 70 * eps^2, where eps^2 = 235 / fy)
limit_700 = 70 * (235 / 700)   # 23.500
limit_355 = 70 * (235 / 355)   # 46.338
limit_235 = 70 * (235 / 235)   # 70.000

# Specified project bounds
bound_col_t = (2.5, 5.0)
bound_col_d = (48.3, 114.3)

bound_brc_t = (1.0, 4.0)
bound_brc_d = (10.0, 50.0)

# ---------------------------------------------------------
# 3. Plotting Execution Engine
# ---------------------------------------------------------

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 8.5))

def plot_panel_hybrid(ax, data, limit, title, max_t, max_d, t_label, d_label, design_bounds=None, legend_loc='upper left', legend_fontsize=11.0):
    t_space = np.linspace(0.01, max_t + 5, 5000)
    d_limit = limit * t_space

    ax.fill_between(t_space, 0, d_limit, color='#2ecc71', alpha=0.12, label='Class 1 & 2')
    ax.fill_between(t_space, d_limit, max_d + 150, color='#e74c3c', alpha=0.12, label='Class 3 & 4')

    ax.plot(t_space, d_limit, color='#d35400', linestyle='-', linewidth=2.8, antialiased=True)

    if design_bounds:
        bound_t, bound_d = design_bounds
        rect = patches.Rectangle(
            (bound_t[0], bound_d[0]),
            bound_t[1] - bound_t[0],
            bound_d[1] - bound_d[0],
            linewidth=2.8,
            edgecolor='#2980b9',
            facecolor='none',
            linestyle='--',
            label='Specified Design Bounds',
            zorder=6
        )
        ax.add_patch(rect)

    comp_t, comp_d = [], []
    slend_t, slend_d = [], []
    for d_val, t_val in data:
        if d_val / t_val <= limit:
            comp_t.append(t_val)
            comp_d.append(d_val)
        else:
            slend_t.append(t_val)
            slend_d.append(d_val)

    if comp_t:
        ax.scatter(comp_t, comp_d, color='#1abc9c', edgecolors='#116b57',
                   s=60, linewidths=1.2, zorder=7, label='Compliant Production Profile')
    if slend_t:
        ax.scatter(slend_t, slend_d, color='#962d22', edgecolors='#4a110b',
                   s=60, linewidths=1.2, zorder=7, label='Non-Compliant Production Profile')

    ax.set_xlim(0, max_t)
    ax.set_ylim(0, max_d)
    ax.set_xlabel(t_label, fontsize=17, fontweight='normal', labelpad=10)
    ax.set_ylabel(d_label, fontsize=17, fontweight='normal', labelpad=10)
    ax.set_title(title, fontsize=19, fontweight='normal', pad=18)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=11)
    ax.grid(True, linestyle=':', alpha=0.5, linewidth=0.7)
    ax.legend(loc=legend_loc, fontsize=legend_fontsize, framealpha=0.95, fancybox=True,
              edgecolor='#555555', shadow=True)

# Generate Plot 1: Columns S700
plot_panel_hybrid(
    ax1, s700_data, limit_700, 
    'Strenx 700 Columns\n(Production Profiles & Initial Bounds)', 
    11, 340, r'Wall Thickness, $t_0$ [mm]', r'Outer Diameter, $d_0$ [mm]', 
    design_bounds=(bound_col_t, bound_col_d),
    legend_fontsize=10.5
)

# Generate Plot 2: Braces S355
plot_panel_hybrid(
    ax2, s355_data, limit_355, 
    'S355 Braces\n(Production Profiles & Initial Bounds)', 
    10, 65, r'Wall Thickness, $t_1$ [mm]', r'Outer Diameter, $d_1$ [mm]', 
    design_bounds=(bound_brc_t, bound_brc_d),
    legend_loc='lower right'
)

# Generate Plot 3: Braces S235
plot_panel_hybrid(
    ax3, s235_data, limit_235, 
    'S235 Braces\n(Production Profiles & Initial Bounds)', 
    13, 110, r'Wall Thickness, $t_1$ [mm]', r'Outer Diameter, $d_1$ [mm]', 
    design_bounds=(bound_brc_t, bound_brc_d),
    legend_loc='lower right'
)

# Additional figure: only S700 and S355
fig2, (ax4, ax5) = plt.subplots(1, 2, figsize=(16, 8.5))
plot_panel_hybrid(
    ax4, s700_data, limit_700,
    'Strenx 700 Columns\n(Production Profiles & Initial Bounds)', 
    11, 340, r'Wall Thickness, $t_0$ [mm]', r'Outer Diameter, $d_0$ [mm]', 
    design_bounds=(bound_col_t, bound_col_d),
    legend_fontsize=11.0
)
plot_panel_hybrid(
    ax5, s355_data, limit_355, 
    'S355 Braces\n(Production Profiles & Initial Bounds)', 
    10, 65, r'Wall Thickness, $t_1$ [mm]', r'Outer Diameter, $d_1$ [mm]', 
    design_bounds=(bound_brc_t, bound_brc_d),
    legend_loc='lower right'
)

# Additional figure: only S355 and S235
fig3, (ax6, ax7) = plt.subplots(1, 2, figsize=(16, 8.5))
plot_panel_hybrid(
    ax6, s355_data, limit_355,
    'S355 Braces\n(Production Profiles & Initial Bounds)', 
    13, 110, r'Wall Thickness, $t_1$ [mm]', r'Outer Diameter, $d_1$ [mm]', 
    design_bounds=(bound_brc_t, bound_brc_d),
    legend_loc='lower right'
)
plot_panel_hybrid(
    ax7, s235_data, limit_235,
    'S235 Braces\n(Production Profiles & Initial Bounds)', 
    13, 110, r'Wall Thickness, $t_1$ [mm]', r'Outer Diameter, $d_1$ [mm]', 
    design_bounds=(bound_brc_t, bound_brc_d),
    legend_loc='lower right'
)

fig.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

save_dir = os.path.dirname(os.path.abspath(__file__))
for name, f in [('Class2_all3', fig), ('Class2_S700_S355', fig2), ('Class2_S355_S235', fig3)]:
    f.savefig(os.path.join(save_dir, f'{name}.png'),
              dpi=1200, bbox_inches='tight', pad_inches=0.15,
              facecolor='white', edgecolor='none')
    f.savefig(os.path.join(save_dir, f'{name}.pdf'),
              bbox_inches='tight', pad_inches=0.15,
              facecolor='white', edgecolor='none')
    f.savefig(os.path.join(save_dir, f'{name}.svg'),
              bbox_inches='tight', pad_inches=0.15,
              facecolor='white', edgecolor='none')

plt.show()