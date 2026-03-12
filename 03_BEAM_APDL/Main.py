## MAIN BEAM ELEMENT DOCUMENT ##
# Pseudo code
    # 1. Read Solidworks Geometric Data
    # 2. Convert Coordinates to APDL language
    #   2.1 Create nodes and lines accordingly
    # 3. Set up boundary conditions in APDL
    # 4. Solve the system
    # 5. Post process Results


# Import tools
from GeoCreate import GeoCreateFun
import SW_Import as SW


# Import SW coordinates as list
SW_filename = "NotSoSimpleFrameBrace.IGS"
SWcoor = SW.import_SW(SW_filename)

# Specify tube dimensions
R0 = 70.1/2 # Corner Tube inner diameter [mm]
R1 = 76.1/2 # CornSer Tube outer diameter [mm]
R2 = 22.3/2 # Brace Tube inner diameter  [mm]
R3 = 26.9/2 # Brace Tube outer diameter  [mm]

var = [R0, R1, R2, R3] # Assembly variables

esize = 100
# Call GeoCreateFun
GeoCreateFun(SWcoor,var,esize)