# Plotting.py
# This script is for plotting various items in the optimization scheme
# Make functions that can be called, this is easier to uncomment 

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import pandas as pd
import matplotlib.pyplot as plt


def plot_aggregate(csv_path):
    aggregate = "Unknown"
    p_value = "N/A"

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
            elif "p_value" in line:
                p_value = line.split(":")[-1].strip()

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
    plt.title(f"{aggregate} aggregation (p = {p_value}) | average (avg_v/avg_g) = {ratio:.3f}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# With p-norm
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_18-47-50.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_19-09-29.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_19-53-15.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_20-09-59.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_21-15-31.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_22-13-22.csv")

# With p-norm-mean
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_22-42-06.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_22-59-39.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_23-27-24.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-23_23-56-21.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-24_00-48-02.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-24_01-03-34.csv")
