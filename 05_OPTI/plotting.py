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

# With KS
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_12-54-09.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_14-08-22.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_15-22-39.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_16-49-51.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_17-59-45.csv")


# With KS Stable
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_19-17-51.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_20-32-20.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_21-53-59.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-28_23-11-02.csv")
plot_aggregate("Optimization_Logs/objective_history_2026-04-29_00-22-51.csv")


