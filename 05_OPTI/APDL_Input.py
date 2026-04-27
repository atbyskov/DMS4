### ADPL_Input.py ###
# -> INPUT:
#       - [SW_coor]   -> Coordinates from Soldiworks
#       - [var]       -> Radius variables
#       - [Misc]      -> Miscellaneous Data (force, mesh etc.) 
# -> OUTPUT:
#       - .txt Input file for APDL Analysis (Eigenbuckling and Nonlinear)

# Pseudo Code
    # Input is coordinate list with format [(x1,y1,z1),(x2,y2,z3)]
    # Create Two Keypoints with these coordinates and create a line with them
    #   K,ID,X,Y,Z
    #   L,P1,P2
    # Sort lines into two groups:
    #   1. Vertical / Corner beams, that only varies in y direction
    #   2. Brace beams, so everyone else
    # Define element type with ET,1,BEAM189
    # Define Cross Section
    # Define material properties
    # Apply SECTYPE,1 to vertical beams
    # Apply SECTYPE,2 to brace beams
    # Apply SECTYPE,3 to Top (constant)

# Import packages
import os
import SW_Import as SW


def InputFun(var,Misc,opti_settings):
    SW_filename = Misc["SW_filename"]

    # Initialize APDL Command for PyMAPDL
    SWcoor = SW.import_SW(os.path.join("IGS",SW_filename))
    ap = []
    

    # Only use Rad if selected in main
    rad = var["rad"]["value"] if "rad" in var and var["rad"].get("active", True) else None

    n = opti_settings["n_mast_segments"]
    mast_height = opti_settings["mast_segment_height"]
    brace_start = n + 1
    top_sectype = 99

    # Convert to Radii
    if opti_settings.get("multi_size_columns", True):
        R0_list = [round(var[f"d0_{i}"]["value"]/2 - var[f"t0_{i}"]["value"], 4) for i in range(1, n+1)]
        R1_list = [round(var[f"d0_{i}"]["value"]/2, 4) for i in range(1, n+1)]
    else:
        R0 = round(var["d0"]["value"]/2 - var["t0"]["value"], 4)
        R1 = round(var["d0"]["value"]/2, 4)

    if opti_settings.get("multi_size_braces", True):
        if opti_settings.get("brace_split", False):
            # Horizontal braces
            R2_horiz_list = [round(var[f"d1_h_{i}"]["value"]/2 - var[f"t1_h_{i}"]["value"], 4) for i in range(1, n+1)]
            R3_horiz_list = [round(var[f"d1_h_{i}"]["value"]/2, 4) for i in range(1, n+1)]
            # Cross braces
            R2_cross_list = [round(var[f"d1_c_{i}"]["value"]/2 - var[f"t1_c_{i}"]["value"], 4) for i in range(1, n+1)]
            R3_cross_list = [round(var[f"d1_c_{i}"]["value"]/2, 4) for i in range(1, n+1)]
        else:
            # Combined braces (same size for both horizontal and cross)
            R2_list = [round(var[f"d1_{i}"]["value"]/2 - var[f"t1_{i}"]["value"], 4) for i in range(1, n+1)]
            R3_list = [round(var[f"d1_{i}"]["value"]/2, 4) for i in range(1, n+1)]
    else:
        if opti_settings.get("brace_split", False):
            R2_horiz = round(var["d1_h"]["value"]/2 - var["t1_h"]["value"], 4)
            R3_horiz = round(var["d1_h"]["value"]/2, 4)
            R2_cross = round(var["d1_c"]["value"]/2 - var["t1_c"]["value"], 4)
            R3_cross = round(var["d1_c"]["value"]/2, 4)
        else:
            R2 = round(var["d1"]["value"]/2 - var["t1"]["value"], 4)
            R3 = round(var["d1"]["value"]/2, 4)

    # Change Radius Logic
    if rad is not None:

        import math

        # Check if Column Member 
        def is_column(x1, z1, x2, z2): 
            return (x1 == x2) and (z1 == z2)
        
        # Only take Column Members
        col_loc = {(x1, z1) for x1, y1, z1, x2, y2, z2 in SWcoor if is_column(x1, z1, x2, z2)}

        loc_to_new = {}
        for (x, z) in col_loc: 
            # r0 = \sqrt(x^2+z^2) 
            r0 = math.hypot(x, z) 
            # Scale radius 
            s = rad / 202.07         # NEEDS TO CHANGE IF ANOTHER IGS FILE IS USED
            # Apply scaling to x and z values 
            loc_to_new[(x, z)] = (x * s, z * s)
            
        def adjust(x, y, z):
            if (x,z) in loc_to_new:
                nx,nz = loc_to_new[(x,z)]
                return nx,y,nz
            return x, y, z

        SWcoor = [
            (*adjust(x1, y1, z1), *adjust(x2, y2, z2))
            for x1, y1, z1, x2, y2, z2 in SWcoor
        ]


    # Function to group lines
    def beam_class(p1, p2):
        
        (x1, y1, z1) = p1
        (x2, y2, z2) = p2
        
        tol = 1e-3  # Tolerance for floating point comparison

        # Corner beam if x2 == x1 AND z2 == z1, i.e. these are unchanged. 
        if abs(x1 - x2) < tol and abs(z1 - z2) < tol:
            return "corner"
        # Horizontal brace if y coordinates are the same
        elif abs(y1 - y2) < tol:
            return "horizontal"
        # Cross brace if y coordinates differ
        else:
            return "cross"
        

                        # SETUP
    ap.append("! ===== APDL INPUT FILE ====== !  ")
##################################################################
##################### Eigenbuckling Analysis #####################
##################################################################
    ap.append("!   Eigenbuckling Analysis     !  ")
    ap.append("! ============================ !   ")
    ap.append("/UNITS,MPa ! Set units [mm,Mg,s,C]   ")

    ap.append("! Static Structural Analysis !  ")
    ap.append("/PREP7  ")
    ap.append("ET,1,BEAM189 ! Use BEAM189 Elements    ")
    ## CROSS SECTION
    ap.append("! CROSS SECTION ! ")
    if opti_settings.get("multi_size_columns", True):
        for i in range(1, n+1):
            ap.append(f"SECTYPE,{i},BEAM,CTUBE  ")
            ap.append(f"SECDATA,{R0_list[i-1]},{R1_list[i-1]},8  ")
    else:
        ap.append("SECTYPE,1,BEAM,CTUBE  ")
        ap.append(f"SECDATA,{R0},{R1},8  ")
    
    if opti_settings.get("multi_size_braces", True):
        if opti_settings.get("brace_split", False):
            # Horizontal braces
            for i in range(1, n+1):
                secnum = brace_start + 2*(i-1)
                ap.append(f"SECTYPE,{secnum},BEAM,CTUBE  ")
                ap.append(f"SECDATA,{R2_horiz_list[i-1]},{R3_horiz_list[i-1]},8  ")
            # Cross braces
            for i in range(1, n+1):
                secnum = brace_start + 2*(i-1) + 1
                ap.append(f"SECTYPE,{secnum},BEAM,CTUBE  ")
                ap.append(f"SECDATA,{R2_cross_list[i-1]},{R3_cross_list[i-1]},8  ")
        else:
            # Combined braces (same size for both horizontal and cross)
            for i in range(1, n+1):
                secnum = brace_start + i - 1
                ap.append(f"SECTYPE,{secnum},BEAM,CTUBE  ")
                ap.append(f"SECDATA,{R2_list[i-1]},{R3_list[i-1]},8  ")
    else:
        if opti_settings.get("brace_split", False):
            # Horizontal braces
            ap.append(f"SECTYPE,{brace_start},BEAM,CTUBE  ")
            ap.append(f"SECDATA,{R2_horiz},{R3_horiz},8  ")
            # Cross braces
            ap.append(f"SECTYPE,{brace_start+1},BEAM,CTUBE  ")
            ap.append(f"SECDATA,{R2_cross},{R3_cross},8  ")
        else:
            ap.append(f"SECTYPE,{brace_start},BEAM,CTUBE  ")
            ap.append(f"SECDATA,{R2},{R3},8  ")
    
    ap.append(f"SECTYPE,{top_sectype},BEAM,CTUBE  ")
    ap.append("SECDATA,35.05,38.05,8  ")

    # MATERIAL
    ap.append("! MATERIAL DATA ")
    ap.append(f"MP,EX,1,{Misc['E_mod']} ! [MPa]")
    ap.append("MP,PRXY,1,0.3  ")
    ap.append("MP,DENS,1,7.850E-6 ! [kg/mm^3]")

    # STIFF Material
    ap.append("! INF STIFNESS MATERIAL REGION ABOVE Y=4078 ")
    ap.append("MP,EX,2,2E+09 ")
    ap.append("MP,PRXY,2,0.3 ")
    ap.append("MP,DENS,2,7.850E-6 ")

    # NODES
    ap.append("! KEYPOINT AND LINES !  ")
    # Initialize values used later        
    key_id = 2
    line_id = 1
    if opti_settings.get("multi_size_columns", True):
        corner_lines = [[] for _ in range(n)]
    else:
        corner_lines = []
    if opti_settings.get("multi_size_braces", True):
        horiz_lines = [[] for _ in range(n)]
        cross_lines = [[] for _ in range(n)]
    else:
        horiz_lines = []
        cross_lines = []
    Top_lines = []
    corner_id = 1
    horiz_id = 1
    cross_id = 1
    Top_lines = []        

    kp_dict = {}
    CM_Brace_dict = 0
    CM_Column_dict = 0

    # For Remote Loading
    #ap.append("K,1,0,4179.14,0  ")
    

    # Loop through and create lines
    for x1, y1, z1, x2, y2, z2 in SWcoor:
        
        # First point 
        p1 = (x1, y1, z1)
        if p1 in kp_dict:   # (We check if the point already exists)
            kp1 = kp_dict[p1]
        else:
            kp1 = key_id
            kp_dict[p1] = kp1
            ap.append(f"K,{key_id}, {x1:.6f}, {y1:.6f}, {z1:.6f}  ")
            key_id += 1

        # Second Point
        p2 = (x2, y2, z2)
        if p2 in kp_dict:   # (We check if the point already exists)
            kp2 = kp_dict[p2]
        else:
            kp2 = key_id
            kp_dict[p2] = kp2
            ap.append(f"K,{kp2}, {x2:.6f}, {y2:.6f}, {z2:.6f} ")
            key_id += 1

        # Group the lines points
        group = beam_class(p1,p2)
        # Create line 
        ap.append(f"L,{kp1},{kp2}  ")
        ap.append(f"LSEL,S,LINE,,{line_id} ")

        
        # Split lines into corner, brace or top section
        if y1 > 4078 and y2 > 4078:
            # This is a TOP beam!
            Top_lines.append(line_id)
            ap.append(f"CM,TOPMAT,LINE ")

        else:
            # Not a top beam → classify as corner, horizontal or cross brace
            if group == "corner":
                if opti_settings.get("multi_size_columns", True):
                    segment = int(min(y1, y2) // mast_height) + 1
                    segment = min(max(segment, 1), n)
                    corner_lines[segment-1].append(line_id)
                else:
                    corner_lines.append(line_id)
                ap.append(f"CM,COLUMN_{corner_id},LINE ")
                corner_id += 1
                CM_Column_dict += 1
            elif group == "horizontal":
                if opti_settings.get("multi_size_braces", True):
                    segment = int(min(y1, y2) // mast_height) + 1
                    segment = min(max(segment, 1), n)
                    horiz_lines[segment-1].append(line_id)
                else:
                    horiz_lines.append(line_id)
                ap.append(f"CM,HORIZ_{horiz_id},LINE ")
                horiz_id += 1
                CM_Brace_dict += 1
            else:  # cross
                if opti_settings.get("multi_size_braces", True):
                    segment = int(min(y1, y2) // mast_height) + 1
                    segment = min(max(segment, 1), n)
                    cross_lines[segment-1].append(line_id)
                else:
                    cross_lines.append(line_id)
                ap.append(f"CM,CROSS_{cross_id},LINE ")
                cross_id += 1
                CM_Brace_dict += 1
        # Reset
        ap.append("LSEL,ALL  ")   

        line_id += 1
    ap.append(" ")

    # Printing lines for debug (for debug)
    #print("Corner Lines")
    #print(corner_lines)
    #print("Horizontal Lines")
    #print(horiz_lines)
    #print("Cross Lines")
    #print(cross_lines)

    # ELEMENT DEFINITION
    ap.append("! ELEMENT SIZE !  ")
    ap.append(f"ESIZE,,{Misc['esize']}   ")
    
    # Function to mesh each group seperatly
    # Makes sure each section has the correct cross section
    def group_mesh(block_name, secnum, line_ids):
        
        ap.append(f"! {block_name} !  ")
        ap.append(f"SECNUM,{secnum}  ")
        ap.append("LSEL,ALL ")

        first = True
        for lid in line_ids:
            if first:
                ap.append(f"LSEL,S,LINE,,{lid} ")
                first = False
            else:
                ap.append(f"LSEL,A,LINE,,{lid} ")
        
        ap.append("LMESH,ALL   ")
            
        
    # Save brace counts
    CM_dict = [CM_Column_dict, horiz_id - 1, cross_id - 1]  # Brace counts (horizontal, cross)
            
        
    # Run the function for Corner and Brace
    if opti_settings.get("multi_size_columns", True):
        for i in range(n):
            if corner_lines[i]:
                group_mesh(f"Meshing CORNER Beams Segment {i+1} (SECNUM={i+1})", i+1, corner_lines[i])
    else:
        group_mesh("Meshing CORNER Beams (SECNUM=1)", 1, corner_lines)
    
    if opti_settings.get("multi_size_braces", True):
        if opti_settings.get("brace_split", False):
            for i in range(n):
                # Horizontal braces
                if horiz_lines[i]:
                    secnum = brace_start + 2*i
                    group_mesh(f"Meshing HORIZ Beams Segment {i+1} (SECNUM={secnum})", secnum, horiz_lines[i])
                # Cross braces
                if cross_lines[i]:
                    secnum = brace_start + 2*i + 1
                    group_mesh(f"Meshing CROSS Beams Segment {i+1} (SECNUM={secnum})", secnum, cross_lines[i])
        else:
            # Combined braces - same section number for both horizontal and cross
            for i in range(n):
                if horiz_lines[i]:
                    secnum = brace_start + i
                    group_mesh(f"Meshing HORIZ Beams Segment {i+1} (SECNUM={secnum})", secnum, horiz_lines[i])
                if cross_lines[i]:
                    secnum = brace_start + i
                    group_mesh(f"Meshing CROSS Beams Segment {i+1} (SECNUM={secnum})", secnum, cross_lines[i])
    else:
        if opti_settings.get("brace_split", False):
            group_mesh(f"Meshing HORIZ Beams (SECNUM={brace_start})", brace_start, horiz_lines)
            group_mesh(f"Meshing CROSS Beams (SECNUM={brace_start+1})", brace_start+1, cross_lines)
        else:
            # Combine all braces (horizontal + cross) into one section
            combined_lines = horiz_lines + cross_lines
            group_mesh(f"Meshing ALL Braces (SECNUM={brace_start})", brace_start, combined_lines)
    
    group_mesh(f"Meshing TOP Beam (SECNUM={top_sectype})", top_sectype, Top_lines)

    # Select all lines
    ap.append("LSEL,ALL ")

    # Select all 'top material' lines again
    first = True
    for lid in Top_lines:
        if first:
            ap.append(f"LSEL,S,LINE,,{lid} ")
            first = False
        else:
            ap.append(f"LSEL,A,LINE,,{lid} ")

    ap.append("ESLL,S  ! select elements on these lines ")
    ap.append("EMODIF,ALL,MAT,2  ! modify selected elements to material 2 ")
    ap.append("ALLSEL,ALL  ")

    # Select all nodes with x=0
    ap.append("SELTOL,1.0E-6  ") # Important for node selection
    ap.append("NSEL,S,LOC,X,0")
    ap.append("NSEL,R,LOC,Z,0")
    ap.append("NSEL,U,LOC,Y,4283,6000") # Old value: 4282
    ap.append("*GET,SlaveNum,NODE,0,COUNT")

    ap.append("*DIM,SlaveIDs,ARRAY,SlaveNum")
    ap.append("*VGET,SlaveIDs(1),NODE,,NLIST") # Stores all node IDs

    # Create Master / Independent Node 
    #ap.append("N,99999,0,4.179140091E+03,0")
    ap.append("N,99999,0,4180.14,0 ") # Old value: 4182.1384
    ap.append("*SET,tid,11")
    ap.append("*SET,cid,10")
    ap.append("ET,cid,175")
    ap.append("ET,tid,170")
    ap.append("KEYO,tid,2,1")
    ap.append("KEYO,tid,4,0")
    ap.append("KEYO,cid,12,5")
    ap.append("KEYO,cid,4,0")
    ap.append("KEYO,cid,2,2")
    ap.append("MAT,10")
    ap.append("REAL,10")
    ap.append("TYPE,10")

    ap.append("*DO,ii,1,SlaveNum,1")
    ap.append("    *SET,elemID,8999+ii")
    ap.append("    *SET,nodeID,SlaveIDs(ii)")
    ap.append("    EN,elemID,nodeID")
    ap.append("*ENDDO")

    # Pilot Node Options
    ap.append("*SET,_npilot,99999")
    ap.append("_npilot1=_npilot")
    ap.append("TYPE,tid")
    ap.append("MAT,cid")
    ap.append("REAL,cid")
    ap.append("TSHAPE,PILO")
    ap.append("EN,79999,_npilot")
    ap.append("TSHAPE")
    # Display Cross section
    ap.append("ALLSEL")

    # Select all nodes with x=0
    ap.append("SELTOL,1.0E-6  ") # Important for node selection
    ap.append("NSEL,S,LOC,X,0")
    ap.append("NSEL,R,LOC,Z,0")
    ap.append("NSEL,U,LOC,Y,4283,6000") # Old value: 4282
    ap.append("*GET,SlaveNum,NODE,0,COUNT")

    ap.append("*DIM,SlaveIDs,ARRAY,SlaveNum")
    ap.append("*VGET,SlaveIDs(1),NODE,,NLIST") # Stores all node IDs

    # Create Master / Independent Node 
    #ap.append("N,99999,0,4.179140091E+03,0")
    ap.append("N,99999,0,4180.14,0 ") # Old value: 4182.1384
    ap.append("*SET,tid,11")
    ap.append("*SET,cid,10")
    ap.append("ET,cid,175")
    ap.append("ET,tid,170")
    ap.append("KEYO,tid,2,1")
    ap.append("KEYO,tid,4,0")
    ap.append("KEYO,cid,12,5")
    ap.append("KEYO,cid,4,0")
    ap.append("KEYO,cid,2,2")
    ap.append("MAT,10")
    ap.append("REAL,10")
    ap.append("TYPE,10")

    ap.append("*DO,ii,1,SlaveNum,1")
    ap.append("    *SET,elemID,8999+ii")
    ap.append("    *SET,nodeID,SlaveIDs(ii)")
    ap.append("    EN,elemID,nodeID")
    ap.append("*ENDDO")

    # Pilot Node Options
    ap.append("*SET,_npilot,99999")
    ap.append("_npilot1=_npilot")
    ap.append("TYPE,tid")
    ap.append("MAT,cid")
    ap.append("REAL,cid")
    ap.append("TSHAPE,PILO")
    ap.append("EN,79999,_npilot")
    ap.append("TSHAPE")
    # Display Cross section
    ap.append("ALLSEL")

    # REMOTE POINT APPLICATION
    # Select all nodes with x=0
    ap.append("SELTOL,1.0E-6  ") # Important for node selection
    ap.append("NSEL,S,LOC,X,0")
    ap.append("NSEL,R,LOC,Z,0")
    ap.append("NSEL,U,LOC,Y,4284,6000")
    ap.append("*GET,SlaveNum,NODE,0,COUNT")

    ap.append("*DIM,SlaveIDs,ARRAY,SlaveNum")
    ap.append("*VGET,SlaveIDs(1),NODE,,NLIST") # Stores all node IDs

    # Create Master / Independent Node 
    #ap.append("N,99999,0,4.179140091E+03,0")
    ap.append("N,99999,0,4182.1384,0 ")
    ap.append("*SET,tid,11")
    ap.append("*SET,cid,10")
    ap.append("ET,cid,175")
    ap.append("ET,tid,170")
    ap.append("KEYO,tid,2,1")
    ap.append("KEYO,tid,4,0")
    ap.append("KEYO,cid,12,5")
    ap.append("KEYO,cid,4,0")
    ap.append("KEYO,cid,2,2")
    ap.append("MAT,10")
    ap.append("REAL,10")
    ap.append("TYPE,10")

    # Create slave elements
    #ap.append("*CFOPEN,SlaveNodes")
    #ap.append("*SET,firstnode,SlaveIDs(1)")
    #ap.append("*VWRITE,firstnode")
    #ap.append("(F15.0)")
    #ap.append("*CFCLOSE")
    
    ap.append("*DO,ii,1,SlaveNum,1")
    ap.append("    *SET,elemID,8999+ii")
    ap.append("    *SET,nodeID,SlaveIDs(ii)")
    ap.append("    EN,elemID,nodeID")
    ap.append("*ENDDO")

    # Pilot Node Options
    ap.append("*SET,_npilot,99999")
    ap.append("_npilot1=_npilot")
    ap.append("TYPE,tid")
    ap.append("MAT,cid")
    ap.append("REAL,cid")
    ap.append("TSHAPE,PILO")
    ap.append("EN,79999,_npilot")
    ap.append("TSHAPE")
    # Display Cross section
    #ap.append("/ESHAPE,1 ! Display Cross Section ")
        
    #Create and save .png of the mesh
    ap.append("/SHOW,PNG,,0  ")
    ap.append("/RGB,INDEX,100,100,100,0  ")
    ap.append("/RGB,INDEX,80,80,80,13  ")
    ap.append("/RGB,INDEX,60,60,60,14  ")
    ap.append("/RGB,INDEX,0,0,0,15  ")
    ap.append("/TYPE,,4  ")
    ap.append("/VIEW,,0,0,1  ")
    ap.append("/ANGLE,,30,YM  ")
    ap.append("EPLOT  ")
    ap.append("/SHOW,close  ")
    ap.append("/SHOW,TERM  ")
    

    # RUN STATIC ANALYSIS
    # We use sparse solver with pre-stress on
    ap.append("! SOLUTION !  ")
    ap.append("/SOLU  ")
    ap.append("ANTYPE, STATIC  ")
    ap.append("EQSLV,SPARSE  ")
    ap.append("PSTRES,ON  ")

    # BOUNDARY CONDITIONS
    ap.append(" ! -- BOUNDARY CONDITIONS -- !  ")
    ap.append("ALLSEL,ALL  ")
    

    # Get top and bottom nodes
    ap.append("*GET, NodeYMax, NODE, 0, MXLOC, Y  ")
    ap.append("*GET, NodeYMin, NODE, 0, MNLOC, Y  ")
    ap.append("*GET, NodeXMax, NODE, 0, MXLOC, X  ")

    ap.append("NSEL,S,LOC,Y,NodeYMax")
    ap.append("NSEL,R,LOC,X,NodeXMax")
    ap.append(f"F,ALL,FY,{Misc['Ver_Force']}")
    ap.append(f"F,ALL,FZ,{Misc['Hor_Force']}")
    ap.append("NSEL,ALL")

    # Force at x = 0 / Weight
    ap.append("NSEL,S,NODE,,99999")
    ap.append(f"F,ALL,FY,{Misc['W_Force']}")

    # Fixed displacement at bottom nodes
    ap.append("! Displacement !  ")
    ap.append("ALLSEL,ALL  ")
    ap.append("NSEL,S,LOC,Y,NodeYMin  ")
    ap.append("D,ALL,ALL,0  ")
    ap.append("ALLSEL  ")

    # SOLVE
    ap.append("! Solve the System  ")
    ap.append("SOLVE  ")
    ap.append("FINISH   ")

    # MASS OF ASSEMBLY
    ap.append("! Get and Print Mass  ")
    ap.append("ALLSEL  ")
    ap.append("NSEL,S,LOC,Y,,4050  ")
    ap.append("ESLN  ")
    ap.append("*GET,ecnt,ELEM,0,COUNT  ")
    ap.append("*GET,enum,ELEM,0,NUM,MIN  ")
    ap.append("totvol = 0  ")

    # Loop over each element and get volume
    ap.append("*DO,i,1,ecnt  ")
    ap.append("   *GET,ev,ELEM,enum,VOLU  ")
    ap.append("   totvol = totvol + ev  ")
    ap.append("   enum = ELNEXT(enum)  ")
    ap.append("*ENDDO  ")
    
    # Calculate Mass
    ap.append("dens = 7.85E-6 ! kg/mm^3  ") # Density
    ap.append("Comp_mass = dens*totvol    ")

    # Open and write to file
    ap.append("*CFOPEN,MASS_assembly,txt  ")
    ap.append("   *VWRITE,Comp_mass  ")
    ap.append("   (F12.5)  ")
    ap.append("*CFCLOS  ")
    ap.append("FINISH  ")
    ap.append("ALLSEL,ALL    ")

    # EIGENBUCKLING 
    ap.append("! Eigenbuckling Solution!  ")
    ap.append("/SOLU  ")
    ap.append("ANTYPE,BUCKLE  ")
    ap.append("BUCOPT,LANB,10  ")
    ap.append("MXPAND,ALL  ")
    ap.append("OUTRES,ALL,ALL  ")
    ap.append("SOLVE  ")
    ap.append("FINISH    ")
    
    ap.append("! Retrieve first 10 eigenvalues  ")      
    ap.append("/POST1  ")
    ap.append("*CFOPEN,Eigenvalue1,txt  ")
    ap.append("*DO,jj,1,10,1  ")
    ap.append("   *GET,MS%jj%,MODE,jj,FREQ  ")
    ap.append("   *VWRITE,MS%jj%  ")
    ap.append('   (F15.5)  ')
    ap.append("*ENDDO  ")
    ap.append("*CFCLOS  ")

    # Calculate Imperfection Force
    ap.append("ALLSEL ")
    ap.append("*GET,N_HIGHEST_ret,NODE,0,MXLOC,Y  ")
    ap.append("N_HIGHEST = N_HIGHEST_ret*1E-3  ")
    ap.append("alpha_h = 2/sqrt(N_HIGHEST)  ")
    ap.append("*IF,alpha_h,LT,0.66,THEN  ")
    ap.append("   alpha_h = 0.66  ")
    ap.append("*ELSEIF,alpha_h,GT,1,THEN  ")
    ap.append("   alpha_h = 1  ")
    ap.append("*ENDIF  ")
    ap.append("alpha_m = 2 ! Assumed for now  ")
    ap.append("imp_ang = 1/200 * alpha_h * alpha_m  ")
    ap.append(f"FORCE_IMP = {Misc['Ver_Force']}*imp_ang  ")

    #ap.append("ALLSEL,ALL   FDELE,ALL,ALL   DDELE,ALL,ALL  ")

##################################################################
##################### Nonlinear Analysis #########################
##################################################################

    # NONLINEAR ANALYSIS SETTINGS
    ap.append("! SOLUTION !  ")
    ap.append("/SOLU  ")
    ap.append("ANTYPE, STATIC  ")
    ap.append("NROPT,FULL")
    ap.append("EQSLV,SPARSE")
    ap.append("NLGEOM,ON  ")
    #ap.append("cnvtol,f,,0.005,,0.01")
    ap.append("ARCLEN,ON")
    #ap.append("ARCTRM,L")
    #ap.append("AUTOTS,OFF  ")
    #ap.append("NSUBST,100,100,100") 
    #ap.append("time,1.")
    #ap.append("NEQIT,200")

    # Get top and bottom nodes
    ap.append("*GET, NodeYMax, NODE, 0, MXLOC, Y  ")
    ap.append("*GET, NodeYMin, NODE, 0, MNLOC, Y  ")
    ap.append("*GET, NodeXMax, NODE, 0, MXLOC, X  ")


    ap.append("NSEL,S,LOC,Y,NodeYMax")
    ap.append("NSEL,R,LOC,X,NodeXMax")
    ap.append(f"F,ALL,FY,{Misc['Ver_Force']}")
    ap.append(f"F,ALL,FZ,{Misc['Hor_Force']}")
    ap.append("NSEL,ALL")

    # Force at x = 0 / Weight
    ap.append("NSEL,S,NODE,,99999")
    ap.append(f"F,ALL,FY,{Misc['W_Force']}")
    ap.append(f"F,ALL,FX,FORCE_IMP")

    # Fixed Displacement at Bottom Nodes
    ap.append("! Displacement !  ")
    ap.append("ALLSEL,ALL  ")
    ap.append("NSEL,S,LOC,Y,NodeYMin  ")
    ap.append("D,ALL,ALL,0  ")
    ap.append("ALLSEL  ")

    # SOLVE
    ap.append("! Solve the System  ")
    ap.append("SOLVE  ")
    ap.append("FINISH   ")

    ap.append("! ===== APDL OUTPUT FILE ===== !  ")

    ap.append("/POST1  ")
    ap.append("SET,LAST  ")
    ap.append("ALLSEL,ALL    ")

    # ONLY SELECT BEAM189 ELEMENTS
    ap.append("*GET,E_COUNT,ELEM,0,COUNT    ") 
    ap.append(f"! Number of Columns: {CM_dict[0]}  ")
    ap.append(f"! Number of Horiz Braces : {CM_dict[1]}    ")
    ap.append(f"! Number of Cross Braces : {CM_dict[2]}    ")
    # SET OUTPUT FILE
    ap.append("! Open file to write  ")
    ap.append("*CFOPEN, APDL_Nonlin_Internal,txt  ")
    # LOOP OVER COLUMNS
    ap.append("! Loop over Columns  ")
    ap.append(f"*DO,ii,1,{CM_dict[0]},1  ")
    ap.append("   CMSEL,S,COLUMN_%ii%  ")
    ap.append("   ESLL,S  ") 
    ap.append("   ESEL,R,ENAME,,189  ")
    # FORMAT
    ap.append("   *IF,ii,LT,10,THEN  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS ColMember_",F3.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS ColMember_",F4.0)  ')
    ap.append("   *ENDIF  ")
    # RESULT 
    ap.append("   *GET,nElem,ELEM,0,COUNT  ")
    ap.append("       *VWRITE,'ElemID','NF [N]','My [Nmm]','Mz [Nmm]','Vy [N]','Vz [N]','T [N/mm]','Y_LOC'  ")
    ap.append("       (A12,A20,A20,A20,A20,A20,A20,A20)  ")
    ap.append("       ELEM = 0  ")
    ap.append("       *DO,jj,1,nElem,1  ")
    ap.append("           ELEM = ELNEXT(ELEM)  ")
    ap.append("           *GET,NF,ELEM,ELEM,SMISC,1  ")
    ap.append("           *GET,MY,ELEM,ELEM,SMISC,2  ")
    ap.append("           *GET,MZ,ELEM,ELEM,SMISC,3  ")
    ap.append("           *GET,VY,ELEM,ELEM,SMISC,6  ")
    ap.append("           *GET,VZ,ELEM,ELEM,SMISC,5  ")
    ap.append("           *GET,TQ,ELEM,ELEM,SMISC,4  ")
    ap.append("           NSLE  ")
    ap.append("           *GET,Y_LOC,NODE,0,MNLOC,Y  ")
    ap.append("           *VWRITE,ELEM,NF,MY,MZ,VY,VZ,TQ,Y_LOC  ")
    ap.append("           (F12.0,7E20.8)  ")
    ap.append("       *ENDDO    ")
    ap.append("*ENDDO")

    # LOOP OVER HORIZONTAL BRACES
    ap.append("! Loop over Horizontal Braces  ")
    ap.append(f"*DO,ii,1,{CM_dict[1]},1  ")
    ap.append("   CMSEL,S,HORIZ_%ii%,LINE  ")
    ap.append("   ESLL,S  ")
    ap.append("   ESEL,R,ENAME,,189  ")
    # FORMAT
    ap.append("   *IF,ii,LT,10,THEN  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS HorizMember_",F2.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS HorizMember_",F3.0)  ')
    ap.append("   *ENDIF  ")
    # RESULT 
    ap.append("   *GET,nElem,ELEM,0,COUNT  ")
    ap.append("   *VWRITE,'ElemID','NF [N]','My [Nmm]','Mz [Nmm]','Vy [N]','Vz [N]','T [N/mm]','Y_LOC'  ")
    ap.append("   (A12,7A20)  ")
    ap.append("   elem = 0  ")
    ap.append("   *DO,jj,1,nElem,1  ")
    ap.append("       ELEM = ELNEXT(ELEM)  ")
    ap.append("       *GET,NF,ELEM,ELEM,SMISC,1  ")
    ap.append("       *GET,MY,ELEM,ELEM,SMISC,2  ")
    ap.append("       *GET,MZ,ELEM,ELEM,SMISC,3  ")
    ap.append("       *GET,VY,ELEM,ELEM,SMISC,6  ")
    ap.append("       *GET,VZ,ELEM,ELEM,SMISC,5  ")
    ap.append("       *GET,TQ,ELEM,ELEM,SMISC,4  ")
    ap.append("       NSLE  ")
    ap.append("       *GET,Y_LOC,NODE,0,MNLOC,Y  ")
    ap.append("       *VWRITE,ELEM,NF,MY,MZ,VY,VZ,TQ,Y_LOC  ")
    ap.append("       (F12.0,7E20.8)  ")
    ap.append("   *ENDDO  ")
    ap.append("*ENDDO  ")

    # LOOP OVER CROSS BRACES
    ap.append("! Loop over Cross Braces  ")
    ap.append(f"*DO,ii,1,{CM_dict[2]},1  ")
    ap.append("   CMSEL,S,CROSS_%ii%,LINE  ")
    ap.append("   ESLL,S  ")
    ap.append("   ESEL,R,ENAME,,189  ")
    # FORMAT
    ap.append("   *IF,ii,LT,10,THEN  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS CrossMember_",F2.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS CrossMember_",F3.0)  ')
    ap.append("   *ENDIF  ")
    # RESULT 
    ap.append("   *GET,nElem,ELEM,0,COUNT  ")
    ap.append("   *VWRITE,'ElemID','NF [N]','My [Nmm]','Mz [Nmm]','Vy [N]','Vz [N]','T [N/mm]','Y_LOC'  ")
    ap.append("   (A12,7A20)  ")
    ap.append("   elem = 0  ")
    ap.append("   *DO,jj,1,nElem,1  ")
    ap.append("       ELEM = ELNEXT(ELEM)  ")
    ap.append("       *GET,NF,ELEM,ELEM,SMISC,1  ")
    ap.append("       *GET,MY,ELEM,ELEM,SMISC,2  ")
    ap.append("       *GET,MZ,ELEM,ELEM,SMISC,3  ")
    ap.append("       *GET,VY,ELEM,ELEM,SMISC,6  ")
    ap.append("       *GET,VZ,ELEM,ELEM,SMISC,5  ")
    ap.append("       *GET,TQ,ELEM,ELEM,SMISC,4  ")
    ap.append("       NSLE  ")
    ap.append("       *GET,Y_LOC,NODE,0,MNLOC,Y  ")
    ap.append("       *VWRITE,ELEM,NF,MY,MZ,VY,VZ,TQ,Y_LOC  ")
    ap.append("       (F12.0,7E20.8)  ")
    ap.append("   *ENDDO  ")
    ap.append("*ENDDO  ")
    ap.append("*CFCLOS    ")


        
    return ap



