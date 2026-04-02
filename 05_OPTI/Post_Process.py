# Post_Process.py
    # Takes [var], [Misc], "APDL_Eigen_Internal.txt", "APDL_Nonlin_Internal" as input
    # 6 Functions for calculating:
        # 1. Local Buckling                 [LB]
        # 2. Normal Force                   [NF]
        # 3. Shear                          [S]
        # 4. Torsion                        [T]
        # 5. Bending, Normal and Shear      [BNS]
        # 6. Buckling Resistance            [BR]
        # 7. Interaction                    [IN]
    # Outputs this as a .txt file for evaluation
    # Output highest utilization factor as list
    #   FORMAT:
    #   [Column NF, Brace NF]

# Import tools
from pathlib import Path
import pandas as pd
import numpy as np
import re
import math
from io import StringIO

class PostProcessor:

    def __init__(self):
        # self.df_eigen  = self.read_forces("AnsoutEigen/APDL_Eigen_Internal.txt")
        self.df_nonlin = self.read_forces("Ansout/APDL_Nonlin_Internal.txt")
        self.df_col = self.df_nonlin[self.df_nonlin["Member"].str.startswith("ColMember")].copy()
        self.df_brace = self.df_nonlin[self.df_nonlin["Member"].str.startswith("BraceMember")].copy()

    # Function to read and parse forces 
    # NOTE: I have used Copilot for most of this function, so understanding is low ...
    def read_forces(self,filepath: str) -> pd.DataFrame:
        blocks = []
        current_member = None
        colnames = None
        buffer = []

        with open(filepath, "r") as f:
            lines = f.readlines()

        buffer = []

        for line in lines:
            # Detect new section header
            if line.startswith("NS"):
                if buffer and colnames and current_member:
                    df = pd.read_csv(
                        StringIO("".join(buffer)),
                        sep=r"\s+",
                        engine="python",
                        names=colnames
                    )
                    if "Y_LOC" in df.columns:
                        df = df[df["Y_LOC"] <= 4080]

                    df["Member"] = current_member
                    blocks.append(df)
                
                buffer = []
                current_member = line.strip().replace("NS ","").replace(".","")
                continue

            if "ElemID" in line and "NF" in line:
                tokens = line.split()
                colnames = [tok for tok in tokens if "[" not in tok]
                continue

            if re.match(r"^\s*\d",line):
                buffer.append(line)

        if buffer and colnames and current_member:
            df = pd.read_csv(
                StringIO("".join(buffer)),
                sep=r"\s+",
                engine="python",
                names=colnames
            )
            df["Member"] = current_member
            blocks.append(df)

        return pd.concat(blocks, ignore_index=True)

    # Function to evaluate all constraints at once
    def Util_list(self, var, Misc):

            util_data = {
                "Util_LB": self.Util_LB(var, Misc),
                "Util_NF": self.Util_NF(var, Misc),
                "Util_S": self.Util_S(var, Misc),
                "Util_T": self.Util_T(var, Misc),
                "Util_BNS": self.Util_BNS(var, Misc),
                "Util_BR": self.Util_BR(var, Misc),
                "Util_IN": self.Util_IN(var, Misc)
            }

            return util_data
    
    # Local Buckling  [Timeshenko p. 458]
    def Util_LB(self, var, Misc ):
        # Import Radii
        R0, R1, R2, R3 = var

        # Import E_mod 
        E_mod = Misc["E_mod"]

        # Poissons Ratio
        v = 0.3

        # Initialize np.zeros
        Util_LB = np.zeros(2)
        
        def LocalBucklingFun(df_member,Ro,Ri):

            # Only compressive members are considered
            df_member = df_member[df_member["NF"] < 0].copy()

            # Switch to positve values (following Eurocode)
            df_member["NF"] = df_member["NF"].abs()

            # Critical Force [Timoshenko p. 460]
            P_cr = (2*math.pi*E_mod*(Ro-Ri)**2)/(math.sqrt(3*(1-v**2)))

            # Utilization Ratio
            Util_local = df_member["NF"].max()/P_cr

            # Max Utilization
            util_max = Util_local.max()

            return util_max
        
        Util_LB[0] = LocalBucklingFun(self.df_col,R1,R0)
        Util_LB[1] = LocalBucklingFun(self.df_brace,R3,R2)

        return Util_LB

    # Normal Force [6.2.4]                
    def Util_NF(self, var, Misc):

        # Initialize
        Util_NF = np.zeros(2)

        # Import Radii
        R0, R1, R2, R3 = var 

        # Import E-modulus
        f_y = Misc["f_y"]

        # Function to calculate Util
        def NormalForceFun(df_member,Ro,Ri):

            # Area
            A = np.pi * ((Ro**2) - (Ri**2))  

            # Design Resistance [6.2.4 (6.10)]
            N_Rd = A*f_y 

            # Compression/Tension Utilization (1993-1-1 [6.2.3] and [6.2.4])
            df_member["Util_NF"] = df_member["NF"].abs() / N_Rd

            # Utilization Ratio [1993-1-1 [6.2.4] -  p. 49 (6.9)]
            # This takes maximum normal force of each member
            util_max = df_member["Util_NF"].max()

            # Return maximum value
            return util_max

        Util_NF[0] = NormalForceFun(self.df_col,R1,R0)
        Util_NF[1] = NormalForceFun(self.df_brace,R3,R2)
        
        return Util_NF

    # Shear Force [6.2.6]
    def Util_S(self, var, Misc):
        
        # Initialize
        Util_S = np.zeros(2)

        # Import Radii
        R0, R1, R2, R3 = var        

        # Import Misc
        f_y = Misc["E_mod"]

        # Function        
        def shearFun(df_member,Ro,Ri):

            A = np.pi * ((Ro**2) - (Ri**2)) 
            
            # Av Area [6.2.6 (3) with option (g)]
            Av = 2/math.pi * A

            # Design Plastic Shear Resistance [6.2.6 (6.18)]
            V_cRd = Av*(f_y/math.sqrt(3))

            # Utilization Ratio [6.2.6 (6.17)]
            # We have two directions (y,z), so we take whatever value in each element, that is highest
            df_member["Util_S"] = df_member[["Vy","Vz"]].abs().max(axis=1) / V_cRd

            # Take max 
            util_s_max = df_member["Util_S"].max()

            return util_s_max
        
        Util_S[0] = shearFun(self.df_col,R1,R0)
        Util_S[1] = shearFun(self.df_brace,R3,R2)


        return Util_S

    # Torsion [6.2.7]
    def Util_T(self, var, Misc):

        # Initialize Util_T
        Util_T = np.zeros(2)

        # Import Radii
        R0, R1, R2, R3 = var

        # Import Misc
        f_y = Misc["f_y"]

        # Function to handle col and brace
        
        def torsionFun(df_member,Ro,Ri):

            D0 = Ro*2 # Outer Diameter
            Di = Ri*2 #Inner Diameter

            # Torsion Design Resistance [6.2.7 (6.23)]
            T_Rd = (math.pi/16 * (D0**4 - Di**4)/D0)*f_y/math.sqrt(3)

            # Torsion Check [6.2.7 (6.23)]
            df_member["Util_T"] = df_member["T"].abs() / T_Rd
  
            # Maximum Value
            Util_T_max = df_member["Util_T"].max()

            # Return Max Value
            return Util_T_max

        Util_T[0] = torsionFun(self.df_col,R1,R0)
        Util_T[1] = torsionFun(self.df_brace,R3,R2)

        return Util_T

    # Bending, Normal and Shear [6.2.9]
    def Util_BNS(self, var, Misc):

        # Initialize Util_BNS
        Util_BNS = np.zeros(2)

        # Import Radii
        R0, R1, R2, R3 = var

        # Import f_y
        f_y = Misc["f_y"]

        def bnsFun(df_member,Ro,Ri):
        
            D0 = Ro*2 # Outer Diameter
            Di = Ri*2 # Inner Diameter

            # Area
            A= np.pi * ((Ro**2) - (Ri**2)) 

            # Design Moment Resistance
            M_Rd = (D0**3 - Di**3)/6 * f_y

            # Design Normal Force Resistanc
            N_Rd = A * f_y
  
            # Forces (Here read as abs values)
            N = df_member["NF"].abs()
            My = df_member["My"].abs()
            Mz = df_member["Mz"].abs()

            # Fraction to be used in Utilization Ratio
            red_col = 1 - N / N_Rd

            # Utilization Ratio
            df_member["Util_BNS"] = (My / (M_Rd * red_col))**2 + (Mz / (M_Rd * red_col))**2

            Util_BNS_max = df_member["Util_BNS"].max()

            return Util_BNS_max
        
        Util_BNS[0] = bnsFun(self.df_col,R1,R0)
        Util_BNS[1] = bnsFun(self.df_brace,R3,R2)


        return Util_BNS

    # Buckling Resistance [6.3.1]
    def Util_BR(self, var, Misc):

        # Open and Read Eigenvalue
        with open("Ansout/Eigenvalue1.txt") as f:
            eigenvalues = [float(line.strip()) for line in f if line.strip()]
        a_cr = next(v for v in eigenvalues if v > 0)

        # Utilize Util_IN
        Util_IN = np.zeros(2)

        # Import Radii
        R0, R1, R2, R3 = var

        # Import Misc
        f_y = Misc["f_y"]

        # Function to handle columns and brace
        def bucklingResFun(df_member,Ro,Ri):
            
            # Only consider compressive members
            df_member = df_member[df_member["NF"] < 0].copy()

            # Convert NF to positive values (following Eurocode)
            df_member["NF"] = df_member["NF"].abs()

            # Area
            A = np.pi * ((Ro**2) - (Ri**2))

            # N_cr (a_imp taken from Table 6.2 with Hollow Section Cold Formed)
            a_imp = 0.49
            N_cr = a_cr * df_member["NF"].max()

            # Slenderness ([6.3.1.2] (6.49))
            slen = np.sqrt((A*f_y)/N_cr)

            # Phi ([6.3.1.2] (6.49))
            Phi = 0.5 * (1+a_imp*(slen-0.2)+slen**2)

            # Chi ([6.3.1.2] (6.49))
            Chi = 1/(Phi+np.sqrt(Phi**2-slen**2))

            # Design Buckling Resistance ([6.3.1.1] (6.47))
            N_bRd = Chi*A*f_y

            # Buckling Resistance Check ([6.3.1.1] (6.46))
            df_member["Util_IN"] = df_member["NF"]/N_bRd

            # Take maximum value
            util_max = df_member["Util_IN"].max()

            # Return Max Value
            return util_max
        
        Util_IN[0] = bucklingResFun(self.df_col,R1,R0)
        Util_IN[1] = bucklingResFun(self.df_brace,R3,R2)

        return Util_IN

    # Interaction Force [6.3.3]
    def Util_IN(self, var, Misc):

        # Initialize Util_IN
        Util_IN = np.zeros(2)

        # Import Radii
        R0, R1, R2, R3 = var

        # Imperfection Factor
        a_imp = 0.49

        # Import Misc
        f_y = Misc["f_y"]

        # Open and Read Eigenvalue
        with open("Ansout/Eigenvalue1.txt") as f:
            eigenvalues = [float(line.strip()) for line in f if line.strip()]
        a_cr = next(v for v in eigenvalues if v > 0)

        def interaction(df_member,Ro,Ri):

            # Consider only compression members
            df_member = df_member[df_member["NF"] < 0].copy()

            # Use positive value (following Eurocode)
            df_member["NF"] = df_member["NF"].abs()

            # Diameters
            Do = 2 * Ro
            Di = 2 * Ri

            # Area
            A = math.pi * (Ro**2-Ri**2)

            # N_cr with eigenvalue 
            N_cr = a_cr * df_member["NF"].max()

            # Slenderness, Phi and Chi ([6.3.1.2] (6.49))
            slen = np.sqrt(A*f_y/N_cr)
            Phi = 0.5 * (1+a_imp*(slen-0.2)+slen**2)
            Chi = 1/(Phi+np.sqrt(Phi**2-slen**2))

            # mu factor
            mu = (1-df_member["NF"]/N_cr)/(1-Chi*(df_member["NF"])/N_cr)


            # Determining Psi (Table A.2)
            M_start_y = df_member["My"].iloc[0]
            M_end_y = df_member["My"].iloc[-1]

            if abs(M_start_y) >= abs(M_end_y):
                M1_y = M_start_y
                M2_y = M_end_y
            else:
                M1_y = M_end_y
                M2_y = M_start_y
            
            Psi_y = M2_y / M1_y
            Psi_y = max(min(Psi_y, 1.0),-1.0)

            # Determining Psi (Table A.2)
            M_start_z = df_member["Mz"].iloc[0]
            M_end_z = df_member["Mz"].iloc[-1]

            if abs(M_start_z) >= abs(M_end_z):
                M1_z = M_start_z
                M2_z = M_end_z
            else:
                M1_z = M_end_z
                M2_z = M_start_z
            
            # Psi (Table A.2)
            Psi_z = M2_z / M1_z
            Psi_z = max(min(Psi_z, 1.0),-1.0)

            # C values (Table A.1)
            Cmy = 0.79 + 0.21 * Psi_y + 0.36 * (Psi_y - 0.33)*df_member["NF"]/N_cr
            CmLT = 1
            Cmz = 0.79 + 0.21 * Psi_z + 0.36 * (Psi_z - 0.33)*df_member["NF"]/N_cr

            # K values (Table A.1)    
            k_yy = Cmy*CmLT * (mu/(1-df_member["NF"]/N_cr))
            k_yz = Cmz*CmLT * (mu/(1-df_member["NF"]/N_cr))

            # Design Resistance Values
            N_Rk = A*f_y
            M_Rk = (math.pi*(Do**4-Di**4))/(32*Do) * f_y

            # Utilization [6.3.3 (6.61)]
            util = df_member["NF"]/(Chi*N_Rk) + k_yy * (df_member["My"].abs().max())/(M_Rk) + k_yz*(df_member["Mz"].abs().max())/(M_Rk)
            
            return util.max()

        Util_IN[0] = interaction(self.df_col,R1,R0)
        Util_IN[1] = interaction(self.df_brace,R3,R2)

        return Util_IN

        """
        Writes the following line to Post_Util.txt:

        E_ID   Y_LOC   UTIL_LC[0]   UTIL_LC[1]   UTIL_NF[0]   UTIL_NF[1]
        """

        # Read nonlinear internal APDL output
        df = read_forces("AnsoutEigen/APDL_Eigen_Internal.txt")

        # Compute utilizations
        util_lc = Util_LC(var, Misc)   # returns [for columns, for braces]
        util_nf = Util_NF(var, Misc)   # same structure

        # Choose representative element: maximum |NF|
        idx = df["NF"].abs().idxmax()
        e_id = int(df.loc[idx, "ElemID"])
        y_loc = float(df.loc[idx, "Y_LOC"])

        # Write results to Post_Util.txt
        with open("Post_Util.txt", "w") as f:
            # Header
            f.write(f"{'E_ID':>8} {'Y_LOC':>10} {'UTIL_LC0':>12} {'UTIL_LC1':>12} "
                    f"{'UTIL_NF0':>12} {'UTIL_NF1':>12}\n")

            # Values
            f.write(f"{e_id:8d} {y_loc:10.2f} "
                    f"{util_lc[0]:12.5f} {util_lc[1]:12.5f} "
                    f"{util_nf[0]:12.5f} {util_nf[1]:12.5f}\n")


  

