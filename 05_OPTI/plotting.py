# Plotting.py
# This script is for plotting various items in the optimization scheme
# Make functions that can be called, this is easier to uncomment 

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import pandas as pd
import matplotlib.pyplot as plt

# Function for plotting Aggregate
def plot_aggregate(csv_path):
    aggregate = "Unknown"
    rho_value = "N/A"

    header_lines = []
    data_start_line = None

    # ---------
    # Read file line-by-line
    # ---------
    with open(csv_path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            header_lines.append(line.strip())

            # Extract metadata
            if "Aggregate" in line:
                aggregate = line.split(":")[-1].strip()
            elif "rho_value" in line:
                rho_value = line.split(":")[-1].strip()

            # Detect start of CSV table
            if line.startswith("iteration,"):
                data_start_line = i
                break

    if data_start_line is None:
        raise ValueError("Could not find CSV header line")

    # ---------
    # Read numerical data ONLY
    # ---------
    df = pd.read_csv(
        csv_path,
        skiprows=data_start_line
    )

    g_max = df["g_max"].to_numpy()
    v_agg = df["v_agg"].to_numpy()

    avg_g = np.average(g_max)
    avg_v = np.average(v_agg)

    ratio = avg_v/avg_g

    # ---------
    # Plot
    # ---------
    plt.figure(figsize=(8, 6))
    plt.scatter(g_max, v_agg, s=30, alpha=0.7)


    plt.xlabel("g_max (true max constraint violation)")
    plt.ylabel("v_agg (aggregated constraint)")
    plt.title(f"{aggregate} aggregation (rho = {rho_value}) | average (avg_v/avg_g) = {ratio:.3f}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Function to handle various optimization results
def time_step():
    # Reads "ResultsOptimization.xlsx" 
    # Selects Time step page and gatheres information for each step size
    exc_file = r"C:\Users\atbys\OneDrive - Aalborg Universitet\DMS4 - General\03_Subjects\05_Optimization\ResultsOptimization.xlsx"

    df = pd.read_excel(exc_file,sheet_name="Time Step", header=None)

    start_row1 = 37
    end_row1 = 1250
    start_col1 = 5

    cols = ["it", "eval", "ob", "rad", "d0", "t0", "d1", "t1", "gmax", "vagg"]
    data1 = df.iloc[start_row1:end_row1,start_col1:start_col1+10]

    data_out1 = {
    f"{name}_1": data1.iloc[:, i]
    for i, name in enumerate(cols)
}
    return data_out1
def plot_ts3D(x_key="d0_1", y_key="rad_1", z_key="ob_1"):
    data = time_step()

    df_plot = pd.DataFrame({
        "it": pd.to_numeric(data["it_1"], errors="coerce"),
        "x": pd.to_numeric(data[x_key], errors="coerce"),
        "y": pd.to_numeric(data[y_key], errors="coerce"),
        "z": pd.to_numeric(data[z_key], errors="coerce"),
    })

    df_plot = df_plot.dropna()
    df_plot = df_plot.sort_values("it")

    # 👇 THIS is the key line
    df_first = df_plot.groupby("it", as_index=False).first()

    x = df_first["x"]
    y = df_first["y"]
    z = df_first["z"]

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    sc = ax.scatter(x, y, z, c=z)
    for xi, yi, zi, iti in zip(x, y, z, df_first["it"]):
        ax.text(xi, yi, zi, f"{int(iti)}", fontsize=8)

    fig.colorbar(sc, label=z_key)

    ax.set_xlabel(x_key)
    ax.set_ylabel(y_key)
    ax.set_zlabel(z_key)

    ax.set_title("3D Optimization Plot (First per Iteration)")

    plt.show()

plot_ts3D("d0_1", "rad_1", "ob_1")
