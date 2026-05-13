# Plotting.py
# This script is for plotting various items in the optimization scheme
# Make functions that can be called, this is easier to uncomment 

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# Plotting ACS
#   Plots ACS development over time
def plot_ACS():
    # ---- Load CSV ----
    filename = "Optimization_Logs/objective_history_2026-05-12_08-48-10.csv"
    filename = "Optimization_Logs/objective_history_2026-05-12_14-46-23.csv"
    #filename = "Optimization_Logs/objective_history_2026-05-12_14-59-18.csv"
    
    df = pd.read_csv(filename, comment='#', skiprows=20)   
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
def plot_obj():
    csv_files = [f for f in os.listdir("Optimization_Logs") if f.endswith(".csv")]

    cols = ["objective"]
    x_col = "iteration"

    fig = plt.figure(figsize=(8,4))


    for file in csv_files:
        path = os.path.join("Optimization_Logs", file)
        df = pd.read_csv(path, skiprows=lambda x: x < 20)
        df.columns = df.columns.str.strip()

        for i, col in enumerate(cols):
            plt.plot(df[x_col].astype(float),
                        df[col].astype(float),
                        label=file)
            
    plt.grid(True)
    plt.xlabel("Iterations")
    plt.ylabel("Mass [kg]")
    plt.title("Objective functions for n=32 sampling points")
    plt.show()
#plot_obj()

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

    folder = "Optimization_Logs"
    csv_files = [f for f in os.listdir(folder) if f.endswith(".csv")]

    fig = plt.figure(figsize=(14, 8))

    # ✅ Create centered layout (your idea, cleaned up)
    ax1 = plt.subplot2grid((2, 6), (0, 0), colspan=2)  # rad
    ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)  # d1
    ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)  # t1

    ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)  # d0 (centered left)
    ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)  # t0 (centered right)

    axes = [ax1, ax2, ax3, ax4, ax5]

    cols = ["x0", "x1", "x2", "x3", "x4"]
    plot_labels = ["rad", "d1", "t1", "d0", "t0"]

    for file in csv_files:
        path = os.path.join(folder, file)

        df = pd.read_csv(path, skiprows=lambda x: x < 20)
        df.columns = df.columns.str.strip()

        for i, col in enumerate(cols):
            axes[i].plot(df[x_col].astype(float),
                         df[col].astype(float),
                         label=file)

    # Formatting
    for i, ax in enumerate(axes):
        ax.set_title(plot_labels[i])
        ax.set_xlabel(x_col_lab)
        ax.grid(True)

    # ✅ Title with spacing fix
    fig.suptitle(title)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.show()
#plot_D_O_All(x_axis="iter")

# Select "eval" or "iter" for x-axis and mode = "all" or "start_end" to plot all or only start/end values

def plot_D_O_start_end_simple():

    folder = "Optimization_Logs"
    csv_files = [f for f in os.listdir(folder) if f.endswith(".csv")]

    title = "Start vs End Design Variables and Objective"

    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()

    cols = ["x0", "x1", "x2", "x3", "x4", "objective"]
    plot_labels = ["rad", "d0", "t0", "d1", "t1", "objective"]

    for file in csv_files:
        path = os.path.join(folder, file)

        df = pd.read_csv(path, skiprows=lambda x: x < 20)
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

    plt.tight_layout()
    plt.show()
#plot_D_O_start_end_simple()


def plot_max_util(folder="Optimization_Logs"):

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
                # Save previous iteration
                if current_iter is not None and util_values:
                    iterations.append(current_iter)
                    max_utils.append(max(util_values))

                # Start new iteration
                current_iter = int(match.group(1))
                collecting = True
                util_values = []
                continue

            # If inside iteration block, extract utilizations
            if collecting and "Util_" in line:

                # Extract all numbers in line (ignore N/A automatically)
                nums = re.findall(r"\d+\.\d+", line)

                # Convert to float
                nums = [float(n) for n in nums]

                util_values.extend(nums)

        # Save last iteration
        if current_iter is not None and util_values:
            iterations.append(current_iter)
            max_utils.append(max(util_values))

        # Plot for this file
        plt.plot(iterations, max_utils, marker='o', label=file)

    # Formatting
    plt.xlabel("Iterations")
    plt.ylabel("Maximum Utilization")
    plt.title("Maximum Utilization per Iteration")
    plt.grid(True)
    #plt.legend(fontsize=7)

    plt.tight_layout()
    plt.show()

plot_max_util()