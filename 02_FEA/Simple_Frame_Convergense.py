# Simple Frame Buckling Analysis
# Convergense Study


import matplotlib.pyplot as plt
import numpy as np

Esize = [60, 30, 15]                # Global Element Size
gamma = [3.7823, 3.5401, 3.4886]    # Load Multiplier from Buckling Analysis

per = np.zeros(len(gamma))

# Relative percentage change
per[0] = 0
per[1] = (gamma[1]-gamma[0])/gamma[0]*100
per[2] = (gamma[2]-gamma[1])/gamma[1]*100

# Print Statistics
print(f"Percentage change: \n {per[0]} \n {per[1]} \n {per[2]}")

# Plot Statistics

plt.figure(figsize=(8,6))

plt.plot(Esize,gamma,"--*",color="Black",label="Data Points")
plt.text(Esize[1] - 2, gamma[1] + 0.005, f"{per[1]:.2f}%", 
             color="red", fontsize=10, fontweight='bold', ha='center')
plt.text(Esize[2] - 0, gamma[2] + 0.01, f"{per[2]:.2f}%", 
             color="red", fontsize=10, fontweight='bold', ha='center')


plt.grid()
plt.gca().invert_xaxis()
plt.xlabel("Global Element Size [mm]")
plt.ylabel("Load Multiplier")
plt.title("Simple Frame Convergense Study")

# Save Figure
plt.savefig(r"02_FEA\Mesh_Convergence_Study.png", dpi=300, bbox_inches='tight')

# Show Figure
plt.show()



