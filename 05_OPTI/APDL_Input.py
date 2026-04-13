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
    # 


def InputFun(SWcoor, var, Misc):

    # Initialize APDL Command for PyMAPDL
    ap = []
    
    # Import variables
    d0, t0, d1, t1 = var

    # Convert to Radii
    R0 = d0/2 - t0       # Column Inner Radius [mm]
    R1 = d0/2            # Column Outer Radius [mm]
    R2 = d1/2 - t1       # Brace Inner Radius  [mm]
    R3 = d1/2            # Brace Outer Radius  [mm]

    R0 = round(R0,4)
    R1 = round(R1,4)
    R2 = round(R2,4)
    R3 = round(R3,4)

    # Import Misc
    esize     = Misc["esize"]
    Hor_Force = Misc["Hor_Force"]
    Ver_Force = Misc["Ver_Force"]
    MomZ      = Misc["MomZ"]
    MomY      = Misc["MomY"]
    f_y       = Misc["f_y"]
    E_mod     = Misc["E_mod"]


    # Function to group lines
    def beam_class(p1, p2):
        
        (x1, y1, z1) = p1
        (x2, y2, z2) = p2

        # Corner beam if x2 == x1 AND z2 == z1, i.e. these are unchanged. 
        if x1 == x2 and z1 == z2:
            return "corner"
        else:
            return "brace"
        

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
    # CROSS SECTION
    ap.append("! CROSS SECTION ! ")
    ap.append("! Corner Type (SECTYPE = 1)  ")
    ap.append("SECTYPE,1,BEAM,CTUBE  ")
    ap.append(f"SECDATA,{R0},{R1},8  ")
    ap.append("! Brace Section Type (SECTYPE = 2)  ")
    ap.append("SECTYPE,2,BEAM,CTUBE  ")
    ap.append(f"SECDATA,{R2},{R3},8    ")
    ap.append("! Top Section Type (SECTYPE= 3 )  ")
    ap.append("SECTYPE,3,BEAM,CTUBE  ")
    ap.append("SECDATA,35.05,38.05,8  ")
    # MATERIAL
    ap.append("! MATERIAL DATA ")
    ap.append(f"MP,EX,1,{E_mod} ! [MPa]")
    ap.append("MP,PRXY,1,0.3  ")
    ap.append("MP,DENS,1,1.7850E-6 ! [kg/mm^3]")
    # STIFF Material
    ap.append("! INF STIFNESS MATERIAL REGION ABOVE Y=4070 ")
    ap.append("MP,EX,2,2E+09 ")
    ap.append("MP,PRXY,2,0.3 ")
    ap.append("MP,DENS,2,1.7850E-6 ")

    # NODES
    ap.append("! KEYPOINT AND LINES !  ")
    # Initialize values used later        
    key_id = 1
    line_id = 1
    corner_lines = []
    brace_lines = []
    corner_id = 1
    brace_id = 1
    Top_lines = []        

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
            ap.append(f"K,{key_id}, {x1:.3f}, {y1:.3f}, {z1:.3f}  ")
            key_id += 1

        # Second Point
        p2 = (x2, y2, z2)
        if p2 in kp_dict:   # (We check if the point already exists)
            kp2 = kp_dict[p2]
        else:
            kp2 = key_id
            kp_dict[p2] = kp2
            ap.append(f"K,{kp2}, {x2:.3f}, {y2:.3f}, {z2:.3f} ")
            key_id += 1

        # Group the lines points
        group = beam_class(p1,p2)
        # Create line 
        ap.append(f"L,{kp1},{kp2}  ")
        ap.append(f"LSEL,S,LINE,,{line_id} ")

        
        # Split lines into corner, brace or top section
        if y1 > 4070 and y2 > 4070:
            # This is a TOP beam!
            Top_lines.append(line_id)
            ap.append(f"CM,TOPMAT,LINE ")

        else:
            # Not a top beam → classify as corner or brace
            if group == "corner":
                corner_lines.append(line_id)
                ap.append(f"CM,COLUMN_{corner_id},LINE ")
                corner_id += 1
                CM_Column_dict += 1
            else:
                brace_lines.append(line_id)
                ap.append(f"CM,BRACE_{brace_id},LINE ")
                brace_id += 1
                CM_Brace_dict += 1
        
        # Reset
        ap.append("LSEL,ALL  ")   

        line_id += 1
    ap.append(" ")

    # ELEMENT DEFINITION
    ap.append("! ELEMENT SIZE !  ")
    ap.append(f"ESIZE,{esize}   ")
    
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
            
        
    # Run the function for Corner and Brace
    group_mesh("Meshing CORNER Beams (SECNUM=1)",1, corner_lines)
    group_mesh("Meshing BRACE Beams  (SECNUM=2)",2, brace_lines)
    group_mesh("Meshing TOP Beam     (SECNUM=3)",3, Top_lines)

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

    # Display Cross section
    #ap.append("/ESHAPE,1 ! Display Cross Section ")

    # Remote Point for Moment Application
    
    # Select all nodes with x=0
    ap.append("NSEL,S,LOC,X,0")
    ap.append("NSEL,U,LOC,Y,4286,5000")
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
        
    #Create and save .png of the mesh
    #ap.append("/SHOW,PNG,,0  ")
    #ap.append("/RGB,INDEX,100,100,100,0  ")
    #ap.append("/RGB,INDEX,80,80,80,13  ")
    #ap.append("/RGB,INDEX,60,60,60,14  ")
    #ap.append("/RGB,INDEX,0,0,0,15  ")
    #ap.append("/TYPE,,4  ")
    #ap.append("/VIEW,,0,0,1  ")
    #ap.append("/ANGLE,,30,YM  ")
    #ap.append("EPLOT  ")
    #ap.append("/SHOW,close  ")
    #ap.append("/SHOW,TERM  ")

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
    ap.append("SELTOL,1.0E-6  ") # Important for node selection

    # Get top and bottom nodes
    ap.append("*GET, NodeYMax, NODE, 0, MXLOC, Y  ")
    ap.append("*GET, NodeYMin, NODE, 0, MNLOC, Y  ")

    ap.append("NSEL,S,LOC,Y,NodeYMax")
    ap.append(f"F_HOR = {Hor_Force}  ")
    ap.append(f"F_VER = {-Ver_Force}  ")
    ap.append("NSEL,ALL")

    # Moment Application
    ap.append("ALLSEL")
    ap.append("NSEL,S,NODE,,99999")
    ap.append(f"F,ALL,MX,0  ")
    ap.append(f"F,ALL,MY,{MomY}  ")
    ap.append(f"F,ALL,MZ,{MomZ}  ")
    ap.append("NSEL,ALL")

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

    # Save column and brace list
    CM_dict = [CM_Column_dict, CM_Brace_dict]

    ap.append("! ===== APDL OUTPUT FILE ===== !  ")

    ap.append("/POST1  ")
    ap.append("SET,LAST  ")
    ap.append("ALLSEL,ALL    ")

    # ONLY SELECT BEAM189 ELEMENTS
    ap.append("*GET,E_COUNT,ELEM,0,COUNT    ") 
    ap.append(f"! Number of Columns: {CM_dict[0]}  ")
    ap.append(f"! Number of Braces : {CM_dict[1]}    ")
    # SET OUTPUT FILE
    ap.append("! Open file to write  ")
    ap.append("*CFOPEN, APDL_Eigen_Internal,txt  ")
    # LOOP OVER COLUMNS
    ap.append("! Loop over Columns  ")
    ap.append(f"*DO,ii,1,{CM_dict[0]},1  ")
    ap.append("   CMSEL,S,COLUMN_%ii%  ")
    ap.append("   ESLL,S  ") 
    ap.append("   ESEL,R,ENAME,,189  ")
    # FORMAT
    ap.append("   *IF,ii,LT,10,THEN  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS ColMember_",F2.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS ColMember_",F3.0)  ')
    ap.append("   *ENDIF  ")
    # RESULT 
    ap.append("   *GET,nElem,ELEM,0,COUNT  ")
    ap.append("   *VWRITE,'ElemID','NF [N]','My [Nmm]','Mz [Nmm]','Vy [N]','Vz [N]','T [N/mm]','Y_LOC'  ")
    ap.append("   (A12,A20,A20,A20,A20,A20,A20,A20)  ")
    ap.append("   ELEM = 0  ")
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
    ap.append("   *ENDDO    ")
    ap.append("*ENDDO")

    # LOOP OVER BRACES
    ap.append("! Loop over Braces  ")
    ap.append(f"*DO,ii,1,{CM_dict[1]},1  ")
    ap.append("   CMSEL,S,BRACE_%ii%,LINE  ")
    ap.append("   ESLL,S  ")
    ap.append("   ESEL,R,ENAME,,189  ")
    # FORMAT
    ap.append("   *IF,ii,LT,10,THEN  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS BraceMember_",F2.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS BraceMember_",F3.0)  ')
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

    # MASS OF ASSEMBLY
    ap.append("! Get and Print Mass  ")
    ap.append("ALLSEL  ")
    ap.append("NSEL,S,LOC,Y,,4080  ")
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
    ap.append(f"FORCE_IMP = {Ver_Force}*imp_ang  ")

    #ap.append("ALLSEL,ALL   FDELE,ALL,ALL   DDELE,ALL,ALL  ")

##################################################################
##################### Nonlinear Analysis #########################
##################################################################

    # NONLINEAR ANALYSIS SETTINGS
    ap.append("! SOLUTION !  ")
    ap.append("/SOLU  ")
    ap.append("ANTYPE, STATIC  ")
    ap.append("NLGEOM,ON  ")
    ap.append("ARCLEN,ON  ")
    ap.append("ARCTRM,L  ")
    ap.append("AUTOTS,OFF  ")
    #ap.append("NSUBST,30,100,10  ") ######### Doesnt work ######## 

    # Apply Force
    ap.append("NSEL,S,LOC,X,0   ")
    ap.append("*GET,N_LOW,NODE,,MNLOC,Y  ")
    ap.append("*GET,n_load_c,NODE,0,COUNT  ")

    ap.append("NSEL,S,LOC,Y,NodeYMax")
    ap.append(f"F_HOR = {Hor_Force}  ")
    ap.append(f"F_VER = {-Ver_Force}  ")
    ap.append("NSEL,ALL")

    # Moment Application
    ap.append("ALLSEL")
    ap.append("NSEL,S,NODE,,99999")
    ap.append(f"F,ALL,MX,0  ")
    ap.append(f"F,ALL,MY,{MomY}  ")
    ap.append(f"F,ALL,MZ,{MomZ}  ")
    ap.append("NSEL,ALL")

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
    ap.append(f"! Number of Braces : {CM_dict[1]}    ")
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
    ap.append('       ("NS ColMember_",F2.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS ColMember_",F3.0)  ')
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

    # LOOP OVER BRACES
    ap.append("! Loop over Braces  ")
    ap.append(f"*DO,ii,1,{CM_dict[1]},1  ")
    ap.append("   CMSEL,S,BRACE_%ii%,LINE  ")
    ap.append("   ESLL,S  ")
    ap.append("   ESEL,R,ENAME,,189  ")
    # FORMAT
    ap.append("   *IF,ii,LT,10,THEN  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS BraceMember_",F2.0)  ')
    ap.append("   *ELSE  ")
    ap.append("       *VWRITE,ii  ")
    ap.append('       ("NS BraceMember_",F3.0)  ')
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



    