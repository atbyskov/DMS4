# Video Tracking Main.py

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

path = "VideoTrackData.txt"

# Load Data
df = pd.read_csv(path,sep="\t",decimal=",")
df.columns = df.columns.str.strip() # Remove Header
df = df.loc[:,~df.columns.str.contains("^Unnamed")]

# Extract Variables
t = df["t"].values           # [s]
F = df["F_Cor"].values   # [N]
x = df["x"].values           # [m]

# Remove NaN values
mask = ~np.isnan(x) & ~np.isnan(F)
x_clean = x[mask]
F_clean = F[mask]

# Convert units
x_mm = x_clean*1000                # [mm]
F_kN = F_clean/1000        # [kN]

# Add Failure Point
x_last = x_mm[-1]
F_last = F_kN[-1]
fail_point = x_last, F_last

# Fit Polynomial to data
degree = 3
coeffs = np.polyfit(x_mm,F_kN,degree)
p = np.poly1d(coeffs)
r2 = 1 - np.sum((F_kN - p(x_mm))**2) / np.sum((F_kN - np.mean(F_kN))**2) # Calculate R^2
# generate values for plotting
x_fit = np.linspace(x_mm.min(),x_mm.max(),100)
y_fit = p(x_fit)


# Print Results
print("\nResults from Analysis")
print(f"R² Value: {r2} ")
print(f"P_cr: {F_last}")
print(f"x_max: {x_last}\n")


# Response curve
plt.figure(figsize=(10,6))

plt.plot(x_mm,F_kN,".",color="black",label="Data Points")
plt.plot(x_last, F_last,"x",color="red",markersize=9,markeredgewidth=2,label="Failure Point")
plt.plot(x_fit,y_fit,"-",color="blue",label=f"Polynomial Fit: $R^2$ = {r2:.3f}")

plt.grid()
plt.xlabel("x-displacement [mm]")
plt.ylabel("Force [kN]")
plt.legend()

plt.title("Load Response Curve for Destructive Test (5-section)")

plt.show()




