# Function for importing Solidworks coordinates

def read_iges(filename):
    param_lines = []

    with open(filename, "r") as f:
        for raw_line in f:
            data = raw_line[:72].strip()

            section = raw_line[72:73]

            if section == "P" and data.startswith("110"):
                param_lines.append(data)
    
    return param_lines

def Get_coor(params):
    j



Filename = "SimpleFrameBrace.IGS"
params = read_iges(Filename)

for p in params:
    print(p)