### GeoCreate.py ###
# -> INPUT:
#       - SW_coor (Coordinates from Soldiworks)
# -> OUTPUT:
#       - .txt Input file for APDL

# Pseudo Code
    # Input is coordinate list with format [(x1,y1,z1),(x2,y2,z3)]
    # Create Two Keypoints with these coordinates and create a line with them
    #   K,ID,X,Y,Z
    #   L,P1,P2
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



def PreFun(SWcoor, var, Misc):

    # Import Radii
    R0, R1, R2, R3 = var 

    # Import Misc
    esize, Hor_Force, Ver_Force, f_y, E_mod = Misc

    # Function to group lines
    def beam_class(p1, p2):
        # Corner beam if x2 == x1 AND z2 == z1, i.e. these are unchanged. 
        (x1, y1, z1) = p1
        (x2, y2, z2) = p2

        if x1 == x2 and z1 == z2:
            return "corner"
        else:
            return "brace"
    
    pre_file = "APDL_Pre.txt"

    with open(pre_file, "w") as f:
    # SETUP
        f.write("! ===== APDL INPUT FILE ====== ! \n \n")
        f.write("/UNITS,MPa \n")
        f.write("/PREP7 \n \n")
        f.write("ET,1,BEAM189 ! Use BEAM189 \n \n")
        # CROSS SECTION
        f.write("! CROSS SECTION !\n")
        f.write("! Corner Type (SECTYPE = 1) \n")
        f.write("SECTYPE,1,BEAM,CTUBE \n")
        f.write(f"SECDATA,{R0},{R1},8 \n")
        f.write("! Brace Section Type (SECTYPE = 2) \n")
        f.write("SECTYPE,2,BEAM,CTUBE \n")
        f.write(f"SECDATA,{R2},{R3},8 \n \n")
        # MATERIAL
        f.write("! MATERIAL DATA\n")
        f.write(f"MP,EX,1,{E_mod} ! [MPa]\n")
        f.write("MP,PRXY,1,0.3 \n")
        f.write("MP,DENS,1.7850E-6 ! [kg/mm^3] \n \n")

        # NODES
        f.write("! KEYPOINT AND LINES ! \n")
        
        key_id = 1
        line_id = 1
        corner_lines = []
        brace_lines = []
        corner_id = 1
        brace_id = 1

        kp_dict = {}
        CM_Brace_dict = 0
        CM_Column_dict = 0

        for x1, y1, z1, x2, y2, z2 in SWcoor:
            
            

            # First point
            p1 = (x1, y1, z1)
            if p1 in kp_dict:
                kp1 = kp_dict[p1]
            else:
                kp1 = key_id
                kp_dict[p1] = kp1
                f.write(f"K,{key_id}, {x1:.3f}, {y1:.3f}, {z1:.3f} \n")
                key_id += 1

            # Second Point
            
            p2 = (x2, y2, z2)
            if p2 in kp_dict:
                kp2 = kp_dict[p2]
            else:
                kp2 = key_id
                kp_dict[p2] = kp2
                f.write(f"K,{kp2}, {x2:.3f}, {y2:.3f}, {z2:.3f}\n")
                key_id += 1

            group = beam_class(p1,p2)

            f.write(f"L,{kp1},{kp2} \n")
            f.write(f"LSEL,S,LINE,,{line_id}\n")

            if group == "corner":
                corner_lines.append(line_id)
                f.write(f"CM,COLUMN_{corner_id},LINE \n")
                corner_id += 1
                CM_Column_dict += 1
            else:
                brace_lines.append(line_id)
                f.write(f"CM,BRACE_{brace_id},LINE \n")
                brace_id += 1
                CM_Brace_dict += 1
            
            
            
            f.write("LSEL,ALL \n")   # Reset

            line_id += 1
        f.write("\n")

        # ELEMENT DEFINITION
        f.write("! ELEMENT SIZE ! \n")
        f.write(f"ESIZE,{esize} \n\n")
        
        def group_mesh(block_name, secnum, line_ids):
            
            f.write(f"! {block_name} ! \n")
            f.write(f"SECNUM,{secnum} \n")
            f.write("LSEL,ALL\n")
            f.write("LSEL,NONE\n")


            first = True
            for lid in line_ids:
                if first:
                    f.write(f"LSEL,S,LINE,,{lid}\n")
                    first = False
                else:
                    f.write(f"LSEL,A,LINE,,{lid}\n")
            
            f.write("LMESH,ALL \n\n")
            
        
        group_mesh("Meshing CORNER Beams (SECNUM=1)",1,corner_lines)

        group_mesh("Meshing BRACE beam (SECNUM=2)",2, brace_lines)

        f.write("/ESHAPE,1 ! Display Cross Section\n")

        # BOUNDARY CONDITIONS
        # FORCE
        f.write("\n! -- BOUNDARY CONDITIONS -- ! \n \n! Force \n")
        f.write("ALLSEL,ALL \n")
        f.write("SELTOL,1.0E-6 \n")
        #f.write("*GET, NodeYMax, NODE, 0, MXLOC, Y \n")
        #f.write("NSEL,S,LOC,Y,NodeYMax \n")
        #f.write("F,ALL,FX,300000 !N \n")
        #f.write("ALLSEL,ALL\n\n")
        
        # ONLY FOR SIMPLEFRAME
        f.write("*GET, NodeXMax, NODE, 0, MXLOC, X \n")
        f.write("*GET, NodeXMin, NODE, 0, MNLOC, X \n")
        f.write("*GET, NodeYMax, NODE, 0, MXLOC, Y \n")
        f.write("NSEL,S,LOC,X,NodeXMin \n")
        f.write("NSEL,R,LOC,Y,NodeYMax \n")
        f.write("F,ALL,FY,-300000 \n")
        f.write("ALLSEL,ALL \n")
        f.write("NSEL,S,LOC,X,NodeXMax \n")
        f.write("NSEL,R,LOC,Y,NodeYMax \n")
        f.write(f"F,ALL,FY,{-Ver_Force} \n")
        f.write("ALLSEL,ALL \n")

        f.write("NSEL,S,LOC,X,NodeXMin \n")
        f.write("NSEL,R,LOC,Y,NodeYMax \n")
        f.write(f"F,ALL,FX,{Hor_Force} \n")
        f.write("ALLSEL,ALL \n")

        # FIXED DISPLACEMENT
        f.write("! Displacement ! \n")
        f.write("ALLSEL,ALL \n")
        f.write("*GET, NodeYMin, NODE, 0, MNLOC, Y \n")
        f.write("NSEL,S,LOC,Y,NodeYMin \n")
        f.write("D,ALL,ALL,0 \n")
        f.write("ALLSEL,ALL\n\n")

        # RUN STATIC ANALYSIS
        f.write("! SOLUTION ! \n")
        f.write("/SOLU \n")
        f.write("ANTYPE, STATIC \n")
        f.write("NLGEOM,ON \n")
        f.write("ARCLEN,ON \n")
        f.write("ARCTRM,L \n ")

        # SETTINGS
        f.write("! Solver Settings !\n")
        f.write("AUTOTS,OFF \n")
        f.write("TIME,1 \n")
        f.write("NSUBST,80,100,40 \n")
        f.write("SOLVE \n \n")
        
        # FINISH
        f.write("FINISH \n \n")

        CM_dict = [CM_Column_dict, CM_Brace_dict]
        
        
        print(f"Input written to:  {pre_file}")
        
        return CM_dict, pre_file, kp_dict



    


