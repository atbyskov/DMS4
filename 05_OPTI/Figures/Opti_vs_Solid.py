import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = "LC1_load_displacement.csv"
opti_path = "LC1_opti.csv"
data = pd.read_csv(file_path)
data_opti = pd.read_csv(opti_path)

# Display column names (optional, for debugging)
print(data.columns)
print(data_opti.columns)

# Extract columns
x1 = data['X-Directional Deformation']
y1 = data['Load'] / 1000  # Convert N → kN
x2 = data_opti["x-directional deformation"]
y2 = data_opti["load"] / 1000 # Convert N -> kN


# Plot
plt.figure()
plt.plot(x1, y1, label="Solid Model")
plt.plot(x2, y2, label ="Optimized Beam Model")
plt.plot([0,200],[25.3,25.3],color="black",linestyle="--")
plt.text(110,26,"Vertical Load used in optimization")
plt.legend()

plt.xlim([0,200])

# Add title and labels
plt.title(f"Load Response Curve for LC1 \n Solid and optimized beam model")
plt.xlabel("x-displacement [mm]")
plt.ylabel("Force [kN]")

# Show grid (optional but useful)
plt.grid()

# Show plot
plt.show()
