# Plotting.py
# This script is for plotting various items in the optimization scheme
# Make functions that can be called, this is easier to uncomment 

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.ticker import MaxNLocator

# Plotting ACS
#   Plots ACS development over time
def plot_ACS():
    # ---- Load CSV ----
    filename = "Optimization_Logs/objective_history_2026-05-12_08-48-10.csv"
    filename = "Optimization_Logs/objective_history_2026-05-12_14-46-23.csv"
    #filename = "Optimization_Logs/objective_history_2026-05-12_14-59-18.csv"
    
    df = pd.read_csv(filename, comment='#', skiprows=26)   
    df.columns = df.columns.str.strip()                         # Clean 
    # ---- Extract relevant columns ----
    iter_col = df["eval_index"]
    g_max = df["g_max"]
    v_agg = df["v_agg"]
    # ---- Get LAST value per iteration ----
    df_last = df.groupby("iteration").last().reset_index()
    iterations = df_last["eval_index"]
    g_max_last = df_last["g_max"]
    v_agg_last = df_last["v_agg"]

    # ---- Plot ----
    plt.figure(figsize=(8, 5))
    plt.plot(iterations, v_agg_last, 'k-', linewidth=2, label=r"$\sigma_{PN}$")
    plt.plot(iterations, g_max_last, 'k--', linewidth=2, label=r"$\sigma_{max}$")
    # Axis labels
    plt.xlabel("Iteration", fontsize=12)
    plt.ylabel("Stress measure", fontsize=12)
    # Style
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(frameon=True)
    # Optional: match paper style
    plt.xlim(left=0)
    plt.tight_layout()
    # Save figure
    # plt.savefig("acs_plot.png", dpi=300)

    # Show
    plt.show()

#plot_ACS()

# Objective functions over iterations

def plot_obj(inp):
    csv_files = [f for f in os.listdir("LHS_results") if f.endswith(".csv")]

    # ✅ Mapping from user input → column in CSV + label
    

    col_map = {
        "objective": ("objective", "Mass [kg]", "objective"),
        "d0": ("x1", r"$d_0$", r"$d_0$"),
        "t0": ("x2", r"$t_0$", r"$t_0$"),
        "d1": ("x3", r"$d_1$", r"$d_1$"),
        "t1": ("x4", r"$t_1$", r"$t_1$"),
        "rad": ("x0", "rad", "rad"),
    }



    if inp not in col_map:
        raise ValueError("Choose from: obj, d0, t0, d1, t1, rad")

    y_col, ylabel, title_label = col_map[inp]
    x_col = "iteration"

    fig = plt.figure(figsize=(8,4))

    for file in csv_files:
        path = os.path.join("LHS_results_SLP_Self_written", file)
        df = pd.read_csv(path, skiprows=lambda x: x < 26)
        df.columns = df.columns.str.strip()

        # Ensure numeric
        df["iteration"] = df["iteration"].astype(int)
        df["eval_index"] = df["eval_index"].astype(int)
        df[y_col] = pd.to_numeric(df[y_col], errors="coerce")

        # 🔥 Move last evaluation to next iteration
        last_idx = df.groupby("iteration")["eval_index"].idxmax()
        df.loc[last_idx, "iteration"] += 1

        # Sort
        df = df.sort_values(by=["iteration", "eval_index"])

        # Plot
        plt.plot(df[x_col], df[y_col], label=file)

        # Mark final point
        x_last = df[x_col].iloc[-1]
        y_last = df[y_col].iloc[-1]


        plt.scatter(x_last, y_last, s=50, marker="x")

    plt.grid(True)
    plt.xlabel("Iterations")
    plt.ylabel(f"{ylabel} [mm]")
    #title_label = title_label.replace("_", r"\_")
    plt.title(f"{title_label} for 32 sampling points")
    # plt.legend()
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()


# Plot all simulation design variables and objective function over evaluations
def plot_D_O_All(x_axis):

    if x_axis == "eval":
        x_col = "eval_index"
        x_col_lab = "evaluations"
        title = "Design variables for n=32 sampling points"
    elif x_axis == "iter":
        x_col = "iteration"
        x_col_lab = "iterations"
        title = "Design variables for n=32 sampling points"
    else:
        raise ValueError("x_axis must be either 'eval' or 'iter'")

    folder = "LHS_results_SLP_Self_written"
    csv_files = [f for f in os.listdir(folder) if f.endswith(".csv")]

    fig = plt.figure(figsize=(14, 8))

    ax1 = plt.subplot2grid((2, 6), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
    ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
    ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)
    ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)

    axes = [ax1, ax2, ax3, ax4, ax5]

    cols = ["x0", "x1", "x2", "x3", "x4"]
    plot_labels = ["rad", "d1", "t1", "d0", "t0"]

    for file in csv_files:
        path = os.path.join(folder, file)

        df = pd.read_csv(path, skiprows=lambda x: x < 26)
        df.columns = df.columns.str.strip()

        # ✅ Ensure numeric types
        df["iteration"] = df["iteration"].astype(int)
        df["eval_index"] = df["eval_index"].astype(int)

        # 🔥 Move last evaluation to next iteration
        last_idx = df.groupby("iteration")["eval_index"].idxmax()
        df.loc[last_idx, "iteration"] += 1

        # Optional (recommended for clean plots)
        df = df.sort_values(by=["iteration", "eval_index"])

        for i, col in enumerate(cols):
            axes[i].plot(df[x_col].astype(float),
                         df[col].astype(float),
                         label=file)

    # Formatting
    for i, ax in enumerate(axes):
        ax.set_title(plot_labels[i])
        ax.set_xlabel(x_col_lab)
        ax.grid(True)
        ax.set_ylabel("[mm]")

    fig.suptitle(title)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.show()

# Select "eval" or "iter" for x-axis and mode = "all" or "start_end" to plot all or only start/end values
def plot_D_O_start_end_simple():

    folder = "LHS_results_SLP_Self_written"
    csv_files = [f for f in os.listdir(folder) if f.endswith(".csv")]

    title = "Start vs End Design Variables and Objective"

    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()

    cols = ["x0", "x1", "x2", "x3", "x4", "objective"]
    plot_labels = ["rad", "d0", "t0", "d1", "t1", "objective"]

    for file in csv_files:
        path = os.path.join(folder, file)

        df = pd.read_csv(path, skiprows=lambda x: x < 26)
        df.columns = df.columns.str.strip()

        # ✅ Extract start and end values
        start = df.iloc[0]
        end = df.iloc[-1]

        # ✅ Fake x-axis: 0=start, 1=end
        x_vals = [0, 1]

        for i, col in enumerate(cols):
            y_vals = [float(start[col]), float(end[col])]

            axes[i].plot(x_vals, y_vals, marker='o', label=file)

    # Formatting
    for i, ax in enumerate(axes):
        ax.set_title(plot_labels[i])

        # ✅ Force x-axis to be exactly [0,1]
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["start", "end"])

        ax.grid(True)

    # ✅ FIX: put BEFORE show()
    fig.suptitle(title)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.show()

def plot_max_util(folder="LHS_results_SLP_Self_written"):

    txt_files = [f for f in os.listdir(folder) if f.endswith(".txt")]

    plt.figure(figsize=(10, 6))

    for file in txt_files:
        path = os.path.join(folder, file)

        with open(path, "r") as f:
            lines = f.readlines()

        iterations = []
        max_utils = []

        current_iter = None
        collecting = False
        util_values = []

        for line in lines:

            # Detect iteration
            match = re.search(r"\[ITERATION (\d+)\]", line)
            if match:
                if current_iter is not None and util_values:
                    iterations.append(current_iter)
                    max_utils.append(max(util_values))

                current_iter = int(match.group(1))
                collecting = True
                util_values = []
                continue

            # Extract utilizations
            if collecting and "Util_" in line:
                nums = re.findall(r"\d+\.\d+", line)
                nums = [float(n) for n in nums]
                util_values.extend(nums)

        # Save last iteration
        if current_iter is not None and util_values:
            iterations.append(current_iter)
            max_utils.append(max(util_values))

        # 🔥 Shift iteration index (final value belongs to next iteration)
        iterations = [i + 1 for i in iterations]

        # Plot
        plt.plot(iterations, max_utils, marker='o', label=file)

    plt.xlabel("Iterations")
    plt.ylabel("Maximum Utilization")
    plt.title("Maximum Utilization per Iteration")
    plt.grid(True)

    plt.tight_layout()
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()

# Select what plots to use
#plot_obj("rad")                          # objective, d0, t0, d1, t1, rad
#plot_D_O_All(x_axis="iter")
#plot_D_O_start_end_simple()
plot_max_util()



def objective_summary_table(folder="LHS_results_SLP_Self_written", start_mass=107.89):

    csv_files = [f for f in os.listdir(folder) if f.endswith(".csv")]

    results = []

    for file in csv_files:
        path = os.path.join(folder, file)

        # Load data
        df = pd.read_csv(path, skiprows=lambda x: x < 26)
        df.columns = df.columns.str.strip()

        # Ensure numeric types
        df["iteration"] = df["iteration"].astype(int)
        df["eval_index"] = df["eval_index"].astype(int)
        df["objective"] = df["objective"].astype(float)

        # Design variables
        design_cols = ["x0", "x1", "x2", "x3", "x4"]
        for col in design_cols:
            df[col] = df[col].astype(float)



def objective_summary_table(folder="LHS_results_SLP_Self_written",
                            start_mass=107.89,
                            output_file="objective_summary.txt"):

    csv_files = [f for f in os.listdir(folder) if f.endswith(".csv")]

    results = []

    for file in csv_files:
        path = os.path.join(folder, file)

        # Load data
        df = pd.read_csv(path, skiprows=lambda x: x < 26)
        df.columns = df.columns.str.strip()

        # Ensure numeric types
        df["iteration"] = df["iteration"].astype(int)
        df["eval_index"] = df["eval_index"].astype(int)
        df["objective"] = df["objective"].astype(float)

        design_cols = ["x0", "x1", "x2", "x3", "x4"]
        for col in design_cols:
            df[col] = df[col].astype(float)

        # 🔥 Move last evaluation to next iteration
        last_idx = df.groupby("iteration")["eval_index"].idxmax()
        df.loc[last_idx, "iteration"] += 1

        # Sort
        df = df.sort_values(by=["iteration", "eval_index"]).reset_index(drop=True)

        # Extract designs
        initial_design = df.iloc[0][design_cols].values
        final_design = df.iloc[-1][design_cols].values

        # Objectives
        start_val = start_mass
        end_val = df["objective"].iloc[-1]
        improvement = 100 * (start_val - end_val) / start_val

        results.append(
            [start_val, end_val, improvement] +
            list(initial_design) +
            list(final_design)
        )

    # Column names
    columns = (
        ["Start [kg]", "End [kg]", "Improvement [%]"] +
        [f"{c}_start" for c in design_cols] +
        [f"{c}_end" for c in design_cols]
    )

    summary = pd.DataFrame(results, columns=columns)

    # Sort by best result
    summary = summary.sort_values("End [kg]").reset_index(drop=True)

    # Round values
    summary = summary.round(3)

    # ✅ Convert to string
    table_str = summary.to_string(index=False)

    # ✅ Write to txt file
    with open(output_file, "w") as f:
        f.write("=== Objective + Design Summary ===\n\n")
        f.write(table_str)

    # Also print to console (optional)
    print("\n=== Objective + Design Summary (saved to file) ===")
    print(table_str)

    return summary



objective_summary_table()