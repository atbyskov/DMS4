# Function for importing Solidworks coordinates

import os
import pandas as pd

def read_iges(SW_filename):

    lines110 = []
    data = []

    with open(SW_filename) as f:
        raw = f.readlines()

        for line in raw:
            filter_line = line.strip()

            if filter_line.startswith("110") and "P" in filter_line:
                lines110.append(filter_line)
        return lines110
    
def create_data(lines110):
    coordinates = []

    for line in lines110:
        line = line.replace(";", ",")
        parts = line.split(",")
        coor_only = parts[1:7]

        nums = [float(x.strip()) for x in coor_only]

        coordinates.append(nums)

    return coordinates
    


def import_SW(SW_filename):
    lines110 = read_iges(SW_filename)
    coords = create_data(lines110)
    return coords