# LWC Coordinate File
# INPUT -> .txt file "SW_Coordinates.txt"
# OUTPUT -> Ansys Coordinate File for Geometry Creation (Multiply for each section)

import numpy as np

nSegments = 8
y_offset = 0.810

# Read SW_Coordinates.txt
def load_SW(file_path):

    try:
        with open(file_path, "r") as f:

            next(f) # Skip header
            data = [line.replace(",",".").split() for line in f if line.strip]

            SWxyz = np.array(data, dtype=np.float64)
            return SWxyz

    except FileNotFoundError:
        print(f"ERROR: {file_path} not found")

# Format as
# GroupNumber IDNumber X Y Z
def write_Ans(out_file, base_data, Segments, offset):
    global_id = 1
    
    with open(out_file, "w") as f:
        for s in range(Segments):
            if s == 0:
                current_data = base_data
            else:
                current_data = base_data[3:]

            for coords in current_data:
            
                x = coords[0]
                y = coords[1] + (s*offset)
                z = coords[2]

                x_str = f"{x:.5f}".replace('.', ',')
                y_str = f"{y:.5f}".replace('.', ',')
                z_str = f"{z:.5f}".replace('.', ',')

                f.write(f"{1} \t {global_id} \t {x_str} \t {y_str} \t {z_str} \n")

                global_id += 1


out_file = "Ans_Coordinates.txt"
file_name = "SW_Coordinates.txt"
SWxyz = load_SW(file_name)
write_Ans(out_file, SWxyz, nSegments, y_offset)



