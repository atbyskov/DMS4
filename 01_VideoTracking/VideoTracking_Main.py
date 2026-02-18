# Video Tracking Main.py

import pandas as pd

path = "VideoTrackData.txt"

origin = pd.read_csv(path,sep="\t",decimal=",")
origin = origin.columns.str.strip()

# Assigning values
t = origin.iloc[:, 0] # [s]
F_Cor = origin.iloc[:, 1] # [N]
x = origin.iloc[:, 2] 


