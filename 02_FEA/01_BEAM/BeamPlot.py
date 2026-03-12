# BEAM Plot py
# Reads "ResultsBEAM.txt" and plots the response curve

import pandas as pd
import matplotlib.pyplot as plt

# Filepath
file_path = "ResultsBEAM.txt"

df = pd.read_csv(file_path, sep='\s+',decimal=',',engine='python',encoding='utf-16')

disp = df['Disp'].values
force = df['Force'].values
force = abs(force)

term_disp = disp[-1]
term_force = force[-1]

plt.figure(figsize=(10,6))
plt.plot(disp, force,linewidth=1.5)

plt.plot(term_disp,term_force, 'rx',markersize=6,label="Termination Point")

plt.xlabel("x-displacement [mm]")
plt.ylabel("Vertical Force (ABS) [N]")
plt.grid()
plt.legend(loc='center right')
plt.title("BEAM189 Simple Frame Response Curve")

plt.show()