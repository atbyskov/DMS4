import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = "LC1_load_displacement.csv"
data = pd.read_csv(file_path)

# Display column names (optional, for debugging)
print(data.columns)

# Extract columns
x = data['X-Directional Deformation']
y = data['Load'] / 1000  # Convert N → kN

# Plot
plt.figure()
plt.plot(x, y)

# Add title and labels
plt.title("Load Response Curve for LC1 in FE")
plt.xlabel("x-displacement [mm]")
plt.ylabel("Force [kN]")

# Show grid (optional but useful)
plt.grid()

# Show plot
plt.show()
