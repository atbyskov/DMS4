# Video Tracking Main.py

import pandas as pd

path = "VideoTrackData.txt"

origin = pd.read_csv(path,sep="\t",decimal=",")
origin = origin.columns.str.strip()

print(origin)

