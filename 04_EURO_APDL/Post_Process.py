# PostProcess.py
# Reads "APDL_out.txt" and parameterises values
# Uses conventional text book Beam equations to evaulue failure and utilization ratios
# Stores and prints results for evaluation and later use 

# Import packages
import re
import math


# Import from main
from Main import Misc 
from Main import var

R0, R1, R2, R3 = var
esize, Hor_Force, Ver_Force, f_y, E_mod = Misc

f_y_shear = f_y/math.sqrt(3)

# Read "APDL_out.txt"
# FORMAT: [ElemID, NF, My, Mz, Vy, Vz, T]
def parse_apdl_output(file_path):
    member_data = {}
    current_member = None

    header_pattern = re.compile(r"NS\s+([\w_]+)\.")

    with open(file_path,"r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Check for header
            header_match = header_pattern.search(line)
            if header_match:
                current_member = header_match.group(1)
                member_data[current_member] = []
                continue
            
            if "ElemID" in line:
                continue

            if current_member:
                parts = line.split()
                row_values = [int(float(parts[0]))] + [float(p) for p in parts[1:]]
                member_data[current_member].append(row_values)
                continue
    return member_data

data = parse_apdl_output("Ansout/APDL_out.txt")

# Initialize post_eval
post_eval = {}

# Local Buckling (Timoshenko)
#   P_cr = (2*pi*E*t^2)/(sqrt(3*(1-v^2)))
P_cr_Column = (2*math.pi*E_mod*(R1-R0)**2)/(math.sqrt(3*(1-0.3**2)))
P_cr_Brace  = (2*math.pi*E_mod*(R3-R2)**2)/(math.sqrt(3*(1-0.3**2)))

A_Col = math.pi * (R1**2 - R0**2)               # mm^2
A_Brace = math.pi * (R3**2 - R2**2)             # mm^2

I_Col = (math.pi/4) * (R1**4 - R0**4)           # mm^4
I_Brace = (math.pi/4) * (R3**4 - R2**4)         # mm^4

J_Col = math.pi/2 * (R1**4 - R0**4)             # mm^4
J_Brace = math.pi/2 * (R3**4 - R2**4)           # mm^4

for member_name, element_rows in data.items():
    post_eval[member_name] = []

    if "Col" in member_name:
        P_cr, A, I, y_max, J = P_cr_Column, A_Col, I_Col, R1, J_Col
    elif "Brace" in member_name:
        P_cr, A, I, y_max, J = P_cr_Brace, A_Brace, I_Brace, R3, J_Brace
    else:
        continue

    for row in element_rows:
        elem_id = row[0]
        N = row[1]
        My = row[2] 
        Mz = row[3] 
        Vy = row[4]
        Vz = row[5]
        T = row[6]

        # Local Buckling
        util_Pcr = abs(N) / P_cr if N < 0 else 0

        # Normal Stress
        sig_My = (N / A) - (My * y_max / I)
        sig_Mz = (N / A) - (Mz * y_max / I)

        # Utilization for Normal Force
        util_Ny = abs(sig_My) / f_y 
        util_Nz = abs(sig_Mz) / f_y 

        # Shear Stress from Bending 4V/3A
        tau_max_y = 4*Vy / 3*A
        tau_max_z = 4*Vz / 3*A

        # Util Shear Stress
        util_tau_y =tau_max_y / f_y_shear
        util_tau_z =tau_max_z / f_y_shear

        # Shear Stress Torsion Tr/J
        tau_tor = T*y_max/J

        # Shear Stress Torsion Utilization 
        util_tor = tau_tor / f_y_shear
        
        post_eval[member_name].append([elem_id, util_Pcr, util_Ny, util_Nz, util_tau_y, util_tau_z, util_tor])

# Normal Stress sig_x = N/A - My/I

# Shear and Torsional Stress 

# Write to txt file
output_file = "Post_Eval_Results.txt"

with open(output_file, "w") as f:
    f.write("Evaluation of Beam elements utilization \n")
    f.write("======================================= \n")

    for member_name, results in post_eval.items():
        # Write the Member Header (e.g., NS ColMember_1)
        f.write(f"NS{member_name}\n")
        
        # Write the Column Headers
        # <15 ensures alignment with a 15-character width
        f.write(f"{'ElemID':<10} {'Local Buckling Util':<25} {'sig_y Util':<15} {'sig_z Util':<15} {'tau_y Util'} {'tau_z Util'} {'tau_tor Util'}\n")
        
        # Iterate through the [elem_id, util] pairs
        for row in results:
            elem_id = row[0]
            util_Pcr = row[1]
            util_Ny = row[2]
            util_Nz = row[3]
            util_tau_y = row[4]
            util_tau_z = row[5]
            util_tor = row[6]

            
            # Write the data with specific formatting
            # .0f for ElemID (no decimals), .4f for Utilization
            f.write(f"{elem_id:<10.0f} {util_Pcr:<25.3f} {util_Ny:<15.3f} {util_Nz:<15.3f} {util_tau_y:<15.3f} {util_tau_z:<15.3f} {util_tor:<15.3f}\n")
        
        # Add a newline for visual separation between members
        f.write("\n")