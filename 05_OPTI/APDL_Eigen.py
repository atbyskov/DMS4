### ADPL_Eigen.py ###
# -> INPUT:
#       - [SW_coor]   -> Coordinates from Soldiworks
#       - [var]       -> Radius variables
#       - [Misc]      -> Miscellaneous Data (force, mesh etc.) 
# -> OUTPUT:
#       - .txt Input file for APDL Eigenbuckling Analysis

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

import os

def Eigen_Fun(SWcoor, var, Misc, out_dir = "AnsoutEigen"):
    
    # Import Radii
    R0, R1, R2, R3 = var 

    # Import Misc
    esize, Hor_Force, Ver_Force, Mom, f_y, E_mod = Misc

    # Function to group lines
    def beam_class(p1, p2):
        
        (x1, y1, z1) = p1
        (x2, y2, z2) = p2

        # Corner beam if x2 == x1 AND z2 == z1, i.e. these are unchanged. 
        if x1 == x2 and z1 == z2:
            return "corner"
        else:
            return "brace"
        
    # Create eigen_file    
    eigen_file = os.path.join(out_dir, "APDL_Eigen.txt")
    
    # Open and Edit .txt file
    with open(eigen_file, "w") as f:

                            # SETUP
        f.write("! ===== APDL INPUT FILE ====== ! \n")
        f.write("!   Eigenbuckling Analysis     ! \n")
        f.write("! ============================ ! \n\n")
        f.write("/UNITS,MPa ! Set units [mm,Mg,s,C] \n\n")

        f.write("! Static Structural Analysis ! \n")
        f.write("/PREP7 \n")
        f.write("ET,1,BEAM189 ! Use BEAM189 Elements \n \n")
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
        # Initialize values used later        
        key_id = 1
        line_id = 1
        corner_lines = []
        brace_lines = []
        corner_id = 1
        brace_id = 1

        kp_dict = {}
        CM_Brace_dict = 0
        CM_Column_dict = 0

        # Loop through and create lines
        for x1, y1, z1, x2, y2, z2 in SWcoor:
            
            # First point 
            p1 = (x1, y1, z1)
            if p1 in kp_dict:   # (We check if the point already exists)
                kp1 = kp_dict[p1]
            else:
                kp1 = key_id
                kp_dict[p1] = kp1
                f.write(f"K,{key_id}, {x1:.3f}, {y1:.3f}, {z1:.3f} \n")
                key_id += 1

            # Second Point
            p2 = (x2, y2, z2)
            if p2 in kp_dict:   # (We check if the point already exists)
                kp2 = kp_dict[p2]
            else:
                kp2 = key_id
                kp_dict[p2] = kp2
                f.write(f"K,{kp2}, {x2:.3f}, {y2:.3f}, {z2:.3f}\n")
                key_id += 1

            # Group the lines points
            group = beam_class(p1,p2)
            # Create line 
            f.write(f"L,{kp1},{kp2} \n")
            f.write(f"LSEL,S,LINE,,{line_id}\n")

            # Split lines into corner or brace 
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
            
            # Reset
            f.write("LSEL,ALL \n")   

            line_id += 1
        f.write("\n")

        # ELEMENT DEFINITION
        f.write("! ELEMENT SIZE ! \n")
        f.write(f"ESIZE,{esize} \n\n")
        
        # Function to mesh each group seperatly
        # Makes sure each section has the correct cross section
        def group_mesh(block_name, secnum, line_ids):
            
            f.write(f"! {block_name} ! \n")
            f.write(f"SECNUM,{secnum} \n")
            f.write("LSEL,ALL\n")

            first = True
            for lid in line_ids:
                if first:
                    f.write(f"LSEL,S,LINE,,{lid}\n")
                    first = False
                else:
                    f.write(f"LSEL,A,LINE,,{lid}\n")
            
            f.write("LMESH,ALL \n\n")
            
        # Run the function for Corner and Brace
        group_mesh("Meshing CORNER Beams (SECNUM=1)",1,corner_lines)
        group_mesh("Meshing BRACE beam (SECNUM=2)",2, brace_lines)

        # Display Cross section
        f.write("/ESHAPE,1 ! Display Cross Section\n")

        # Create and save .png of the mesh
        f.write("/SHOW,PNG,,0 \n")
        f.write("/RGB,INDEX,100,100,100,0 \n")
        f.write("/RGB,INDEX,80,80,80,13 \n")
        f.write("/RGB,INDEX,60,60,60,14 \n")
        f.write("/RGB,INDEX,0,0,0,15 \n")
        f.write("/TYPE,,4 \n")
        f.write("/VIEW,,0,0,1 \n")
        f.write("/ANGLE,,30,YM \n")
        f.write("EPLOT \n")
        f.write("/SHOW,close \n")
        f.write("/SHOW,TERM \n")

        # RUN STATIC ANALYSIS
        # We use sparse solver with pre-stress on
        f.write("! SOLUTION ! \n")
        f.write("/SOLU \n")
        f.write("ANTYPE, STATIC \n")
        f.write("EQSLV,SPARSE \n")
        f.write("PSTRES,ON \n")

        # BOUNDARY CONDITIONS
        f.write("\n! -- BOUNDARY CONDITIONS -- ! \n")
        f.write("ALLSEL,ALL \n")
        f.write("SELTOL,1.0E-6 \n") # Important for node selection

        # Get top and bottom nodes
        f.write("*GET, NodeXMax, NODE, 0, MXLOC, X \n")
        f.write("*GET, NodeXMin, NODE, 0, MNLOC, X \n")
        f.write("*GET, NodeYMax, NODE, 0, MXLOC, Y \n")
        f.write("*GET, NodeYMin, NODE, 0, MNLOC, Y \n")
        # Apply force at top node
        f.write("NSEL,S,LOC,Y,NodeYMax \n")
        f.write(f"F,ALL,FY,{-Ver_Force} \n")
        f.write(f"F,ALL,FX,{Hor_Force} \n")
        f.write(f"F,ALL,MZ,{Mom} \n")
        f.write("ALLSEL,ALL \n")

        # Fixed displacement at bottom nodes
        f.write("! Displacement ! \n")
        f.write("ALLSEL,ALL \n")
        f.write("NSEL,S,LOC,Y,NodeYMin \n")
        f.write("NSEL,R,LOC,X,202.07 \n")
        f.write("D,ALL,UY,0 \n")

        f.write("ALLSEL,ALL \n")
        f.write("NSEL,S,LOC,Y,NodeYMin \n")
        f.write("NSEL,R,LOC,Z,-175 \n")
        f.write("NSEL,R,LOC,X,-101.04 \n")

        f.write("ALLSEL,ALL \n")
        f.write("NSEL,S,LOC,Y,NodeYMin \n")
        f.write("NSEL,R,LOC,Z,175 \n")
        f.write("NSEL,R,LOC,X,-101.04 \n")

        

        # SOLVE
        f.write("! Solve the System \n")
        f.write("SOLVE \n")
        f.write("FINISH \n\n")
    
        # Save column and brace list
        CM_dict = [CM_Column_dict, CM_Brace_dict]

        f.write("! ===== APDL OUTPUT FILE ===== ! \n")

        f.write("/POST1 \n")
        f.write("SET,LAST \n")
        f.write("ALLSEL,ALL \n \n")

        # ONLY SELECT BEAM189 ELEMENTS
        f.write("*GET,E_COUNT,ELEM,0,COUNT \n \n") 
        f.write(f"! Number of Columns: {CM_dict[0]} \n")
        f.write(f"! Number of Braces : {CM_dict[1]} \n \n")
        # SET OUTPUT FILE
        f.write("! Open file to write \n")
        f.write("*CFOPEN, APDL_Eigen_Internal,txt \n")
        # LOOP OVER COLUMNS
        f.write("! Loop over Columns \n")
        f.write(f"*DO,ii,1,{CM_dict[0]},1 \n")
        f.write("   CMSEL,S,COLUMN_%ii% \n")
        f.write("   ESLL,S \n") 
        f.write("   ESEL,R,ENAME,,189 \n")
        # FORMAT
        f.write("   *IF,ii,LT,10,THEN \n")
        f.write("       *VWRITE,ii \n")
        f.write('       ("NS ColMember_",F2.0) \n')
        f.write("   *ELSE \n")
        f.write("       *VWRITE,ii \n")
        f.write('       ("NS ColMember_",F3.0) \n')
        f.write("   *ENDIF \n")
        # RESULT 
        f.write("*GET,nElem,ELEM,0,COUNT \n")
        f.write("   *VWRITE,'ElemID','NF [N]','My [Nmm]','Mz [Nmm]','Vy [N]','Vz [N]','T [N/mm]' \n")
        f.write("   (A12,A20,A20,A20,A20,A20,A20) \n")
        f.write("   ELEM = 0 \n")
        f.write("   *DO,jj,1,nElem,1 \n")
        f.write("       ELEM = ELNEXT(ELEM) \n")
        f.write("       *GET,NF,ELEM,ELEM,SMISC,1 \n")
        f.write("       *GET,MY,ELEM,ELEM,SMISC,2 \n")
        f.write("       *GET,MZ,ELEM,ELEM,SMISC,3 \n")
        f.write("       *GET,VY,ELEM,ELEM,SMISC,4 \n")
        f.write("       *GET,VZ,ELEM,ELEM,SMISC,5 \n")
        f.write("       *GET,TQ,ELEM,ELEM,SMISC,6 \n")
        f.write("       *VWRITE,ELEM,NF,MY,MZ,VY,VZ,TQ \n")
        f.write("       (F12.0,6E20.8) \n")
        f.write("   *ENDDO \n \n")
        f.write("*ENDDO")

        # LOOP OVER BRACES
        f.write("! Loop over Braces \n")
        f.write(f"*DO,ii,1,{CM_dict[1]},1 \n")
        f.write("   CMSEL,S,BRACE_%ii%,LINE \n")
        f.write("   ESLL,S \n")
        f.write("   ESEL,R,ENAME,,189 \n")
        # FORMAT
        f.write("   *IF,ii,LT,10,THEN \n")
        f.write("       *VWRITE,ii \n")
        f.write('       ("NS BraceMember_",F2.0) \n')
        f.write("   *ELSE \n")
        f.write("       *VWRITE,ii \n")
        f.write('       ("NS BraceMember_",F3.0) \n')
        f.write("   *ENDIF \n")
        # RESULT 
        f.write("   *GET,nElem,ELEM,0,COUNT \n")
        f.write("   *VWRITE,'ElemID','NF [N]','My [Nmm]','Mz [Nmm]','Vy [N]','Vz [N]','T [N/mm]' \n")
        f.write("   (A12,6A20) \n")
        f.write("   elem = 0 \n")
        f.write("   *DO,jj,1,nElem,1 \n")
        f.write("       ELEM = ELNEXT(ELEM) \n")
        f.write("       *GET,NF,ELEM,ELEM,SMISC,1 \n")
        f.write("       *GET,MY,ELEM,ELEM,SMISC,2 \n")
        f.write("       *GET,MZ,ELEM,ELEM,SMISC,3 \n")
        f.write("       *GET,VY,ELEM,ELEM,SMISC,4 \n")
        f.write("       *GET,VZ,ELEM,ELEM,SMISC,5 \n")
        f.write("       *GET,TQ,ELEM,ELEM,SMISC,6 \n")
        f.write("       *VWRITE,ELEM,NF,MY,MZ,VY,VZ,TQ \n")
        f.write("       (F12.0,6E20.8) \n")
        f.write("   *ENDDO \n")
        f.write("*ENDDO \n")
        f.write("*CFCLOS \n \n")

        # MASS OF ASSEMBLY
        f.write("! Get and Print Mass \n")
        f.write("ALLSEL \n")
        f.write("*GET,ecnt,ELEM,0,COUNT \n")
        f.write("*GET,enum,ELEM,0,NUM,MIN \n")
        f.write("totvol = 0 \n")
        # Loop over each element and get volume
        f.write("*DO,i,1,ecnt \n")
        f.write("   *GET,ev,ELEM,enum,VOLU \n")
        f.write("   totvol = totvol + ev \n")
        f.write("   enum = ELNEXT(enum) \n")
        f.write("*ENDDO \n")
        # Calculate Mass
        f.write("dens = 7.85E-6 ! kg/mm^3 \n") # Desnity
        f.write("Comp_mass = dens*totvol \n \n")
        # Open and write to file
        f.write("*CFOPEN,MASS_assembly,txt \n")
        f.write("   *VWRITE,Comp_mass \n")
        f.write("   (F12.5) \n")
        f.write("*CFCLOS \n")
        f.write("FINISH \n")
        f.write("ALLSEL,ALL \n \n")

        # EIGENBUCKLING 
        f.write("! Eigenbuckling Solution! \n")
        f.write("/SOLU \n")
        f.write("ANTYPE,BUCKLE \n")
        f.write("BUCOPT,LANB,10 \n")
        f.write("MXPAND,ALL \n")
        f.write("OUTRES,ALL,ALL \n")
        f.write("SOLVE \n")
        f.write("FINISH \n \n")
        
        f.write("! Retrieve first 10 eigenvalues \n")
        f.write("/POST1 \n")
        f.write("*CFOPEN,Eigenvalue1,txt \n")
        f.write("*DO,jj,1,10,1 \n")
        f.write("   *GET,MS%jj%,MODE,jj,FREQ \n")
        f.write("   *VWRITE,MS%jj% \n")
        f.write('   (F10.5) \n')
        f.write("*ENDDO \n")
        f.write("*CFCLOS \n")

        
    return eigen_file



    