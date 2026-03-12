### GeoCreate.py ###
# -> INPUT:
#       - SW_coor (Coordinates from Soldiworks)
# -> OUTPUT:
#       - .txt file with node and line information to APDL

# Pseudo Code
    # Input is coordinate list with format [(x1,y1,z1),(x2,y2,z3)]
    # Create Two nodes with these coordinates and create a line with them
    # Sort lines into two groups:
    #   1. Vertical / Corner beams, that only varies in y direction
    #   2. Brace beams, so everyone else
    # Define element type with ET,1,BEAM189
    # Define cross section with 
    #   SECTYPE,1,BEAM,CIRC and SECDATA,R1,R0,N
    #   SECTYPE,2,BEAM,CIRC and SECDATA,R3,R2,N
    # Define material properties (S690 Nonlinear)
    # Apply SECTYPE,1 to vertical beams
    # Apply SECTYPE,2 to brace beams
    # Apply same material and element type to all lines

# Format to be written
# N, Node, X,Y,Z
# L, P1, P2

# Must connect 

