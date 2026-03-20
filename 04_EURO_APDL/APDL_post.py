# OUTPUT.py # 
# Controls output parameters, such as axial forces etc.

def APDL_post_fun(CM_dict):
    post_file = "APDL_post.txt"

    with open(post_file, "w") as f:
        f.write("! ===== APDL OUTPUT FILE ===== ! \n \n")

        f.write("/POST1 \n")
        f.write("SET,LAST \n \n")
        f.write("ALLSEL,ALL \n")

        # ONLY IN GUI MODE
        # f.write("/INPUT,APDL_out,txt \n")

        # ONLY SELECT BEAM189 ELEMENTS
        f.write("*GET,E_COUNT,ELEM,0,COUNT \n \n") 
        f.write(f"! Number of Columns: {CM_dict[0]} \n")
        f.write(f"! Number of Braces : {CM_dict[1]} \n \n")
        # SET OUTPUT FILE
        f.write("*CFOPEN, APDL_out,txt \n \n")
        # LOOP OVER COLUMNS
        f.write("! Loop over Columns \n")
        f.write(f"*DO,ii,1,{CM_dict[0]},1 \n")
        f.write("   CMSEL,S,COLUMN_%ii% \n")

        f.write("ESLL,S \n")
        f.write("ESEL,R,ENAME,,189 \n")
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
        f.write("*GET,nElem,ELEM,0,COUNT \n")
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

        f.write("*CFCLOS \n")


        return post_file
        # 


    