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
        # 8. Brace-Step                     [BS]
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
import SW_Import as SW


# Function to calculate the span of horizontal brace - It works, but i don't understand this
def _brace_span_mm(var, misc):
    co = SW.import_SW(str(Path("IGS") / misc["SW_filename"]))
    if (rad := var.get("rad")) is not None:
        col = {(x1, z1) for x1, y1, z1, x2, y2, z2 in co if x1 == x2 and z1 == z2}
        ln = {p: (p[0]*rad/h, p[1]*rad/h) for p in col for h in [math.hypot(p[0], p[1])+1e-5]}
        adj = lambda x,y,z: (ln[(x,z)][0],y,ln[(x,z)][1]) if (x,z) in ln else (x,y,z)
        co = [(*adj(x1,y1,z1),*adj(x2,y2,z2)) for x1,y1,z1,x2,y2,z2 in co]
    return next(math.hypot(x2-x1,y2-y1,z2-z1) for x1,y1,z1,x2,y2,z2 in co if x1!=x2 or z1!=z2)
class PostProcessor:

    def __init__(self,var,Misc):
        # Store inputs
        self.var = var
        self.Misc = Misc

        # Unpack variables
        self.d0 = var["d0"]
        self.t0 = var["t0"]
        self.d1 = var["d1"]
        self.t1 = var["t1"]
        self.E_mod = Misc["E_mod"]
        self.f_y = Misc["f_y"]
        self.f_y_brace = Misc["f_y_brace"]


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
            if "Y_LOC" in df.columns:
                        df = df[df["Y_LOC"] <= 4080]
            df["Member"] = current_member
            blocks.append(df)

        return pd.concat(blocks, ignore_index=True)

    # Function to evaluate all constraints at once
    def Util_list(self):

            util_data = {
                "Util_LB": self.Util_LB(),
                "Util_NF": self.Util_NF(),
                "Util_S": self.Util_S(),
                "Util_T": self.Util_T(),
                "Util_BNS": self.Util_BNS(),
                "Util_BR": self.Util_BR(),
                "Util_IN": self.Util_IN(),
                "Util_BS": self.Util_BS()
            }

            return util_data
    
    # Local Buckling  [Timeshenko p. 458]
    def Util_LB(self):

        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]

        # Import E_mod 
        E_mod = self.E_mod

        # Poissons Ratio
        v = 0.3

        # Initialize np.zeros
        Util_LB_col = np.zeros(len(self.df_col))
        Util_LB_brace = np.zeros(len(self.df_brace))
        

        def LocalBucklingFun(df_member, Ro, Ri):
            df = df_member.copy()

            P_cr = (2 * math.pi * E_mod * (Ro - Ri)**2) / math.sqrt(3 * (1 - v**2))

            N_raw = df["NF"].to_numpy(dtype=float)
            N_comp = np.maximum(-N_raw, 0.0)   # only compression counts

            util_local = N_comp / P_cr

            return pd.Series(util_local, index=df.index, name="Util_LB")
                
        Util_LB_col = LocalBucklingFun(self.df_col,R1,R0)
        Util_LB_brace = LocalBucklingFun(self.df_brace,R3,R2)

        return Util_LB_col, Util_LB_brace

    # Normal Force [6.2.4]                
    def Util_NF(self):

        # Initialize
        Util_NF = np.zeros(2)


        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]

        # Import E-modulus
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        # Function to calculate Util
        def NormalForceFun(df_member,Ro,Ri,f_y):


            # Area
            A = np.pi * ((Ro**2) - (Ri**2)) 

            # Design Resistance [6.2.4 (6.10)] 

            N_Rd = A*f_y 


            # Compression/Tension Utilization (1993-1-1 [6.2.3] and [6.2.4])
            df_member["Util_NF"] = df_member["NF"].abs() / N_Rd

            # Utilization Ratio [1993-1-1 [6.2.4] -  p. 49 (6.9)]
            # This takes Utilization Ratio of each member
            Util_NF = df_member["Util_NF"]

            # Return Utilization Ratio
            return Util_NF

        Util_NF_col = NormalForceFun(self.df_col,R1,R0,f_y) # Column
        Util_NF_brace = NormalForceFun(self.df_brace,R3,R2,f_y_brace) # Brace
        
        return Util_NF_col, Util_NF_brace

    # Shear Force [6.2.6]
    def Util_S(self):
        
        # Initialize
        Util_S = np.zeros(2)

        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]      

        # Import Misc
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        # Function        
        def shearFun(df_member,Ro,Ri,f_y):

            A = np.pi * ((Ro**2) - (Ri**2)) 
            
            # Av Area [6.2.6 (3) with option (g)]
            Av = 2/math.pi * A 

            # Design Plastic Shear Resistance [6.2.6 (6.18)]
            V_cRd = Av*(f_y/math.sqrt(3))

            # Utilization Ratio [6.2.6 (6.17)]
            # We have two directions (y,z), so we take whatever value in each element, that is highest
            df_member["Util_S"] = df_member[["Vy","Vz"]].abs().max(axis=1) / V_cRd

            # Take max Utilization Ratio
            Util_S = df_member["Util_S"]

            return Util_S
        
        Util_S_col = shearFun(self.df_col,R1,R0,f_y) # Column
        Util_S_brace = shearFun(self.df_brace,R3,R2,f_y_brace) # Brace


        return Util_S_col, Util_S_brace

    # Torsion [6.2.7]
    def Util_T(self):

        # Initialize Util_T
        Util_T_col = np.zeros(len(self.df_col))
        Util_T_brace = np.zeros(len(self.df_brace))


        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]

        # Import Misc
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        # Function to handle col and brace
        
        def torsionFun(df_member,Ro,Ri,f_y):

            D0 = Ro*2 # Outer Diameter
            Di = Ri*2 #Inner Diameter

            # Torsion Design Resistance [6.2.7 (6.23)]
            T_Rd = (math.pi/16 * (D0**4 - Di**4)/D0)*f_y/math.sqrt(3)

            # Torsion Check [6.2.7 (6.23)]
            df_member["Util_T"] = df_member["T"].abs() / T_Rd
  
            # Utilization Ratio
            Util_T = df_member["Util_T"]

            # Return Utilization Ratio
            return Util_T

        Util_T_col = torsionFun(self.df_col,R1,R0,f_y) # Column
        Util_T_brace = torsionFun(self.df_brace,R3,R2,f_y_brace) # Brace

        return Util_T_col, Util_T_brace

    # Bending, Normal and Shear [6.2.9]
    def Util_BNS(self):

        # Initialize Util_BNS
        Util_BNS_col = np.zeros(len(self.df_col))
        Util_BNS_brace = np.zeros(len(self.df_brace))

        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]]

        # Import f_y
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        def bnsFun(df_member,Ro,Ri,f_y):
        
            D0 = Ro*2 # Outer Diameter
            Di = Ri*2 # Inner Diameter
            t = (D0-Di)/2   # Thickness

            # Area
            A = np.pi * ((Ro**2) - (Ri**2)) 
            Aw = (A-2*D0*t)/A
            Aw = np.clip(Aw,0,0.5)

            # Design Moment Resistance
            M_Rd = (D0**3 - Di**3)/6 * f_y

            # Design Normal Force Resistanc
            N_Rd = A * f_y
  
            # Forces (Here read as abs values)
            N = df_member["NF"].abs() 
            My = df_member["My"].abs()
            Mz = df_member["Mz"].abs()

            # Fraction to be used in Utilization Ratio
            red_col = (1 - (N / N_Rd)**1.7)

            # Utilization Ratio
            df_member["Util_BNS"] = (My / (M_Rd * red_col))**2 + (Mz / (M_Rd * red_col))**2

            Util_BNS = df_member["Util_BNS"]

            return Util_BNS
        
        Util_BNS_col = bnsFun(self.df_col,R1,R0,f_y) # Column
        Util_BNS_brace = bnsFun(self.df_brace,R3,R2,f_y_brace) # Brace

        return Util_BNS_col, Util_BNS_brace

    # Buckling Resistance [6.3.1]
    def Util_BR(self):

        # Open and Read Eigenvalue
        with open("Ansout/Eigenvalue1.txt") as f:
            eigenvalues = [float(line.strip()) for line in f if line.strip()]
        a_cr = next(v for v in eigenvalues if v > 0)

        # Utilize Util_IN
        Util_BR_col = np.zeros(len(self.df_col))
        Util_BR_brace = np.zeros(len(self.df_brace))

        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]

        # Import Misc
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        # Function to handle columns and brace
        def bucklingResFun(df_member, Ro, Ri, f_y):
            df = df_member.copy()
            util = np.zeros(len(df), dtype=float)

            A = np.pi * (Ro**2 - Ri**2)

            N_raw = df["NF"].to_numpy(dtype=float)
            N_comp = np.maximum(-N_raw, 0.0)

            active = N_comp > 0.0

            if np.any(active):
                a_imp = 0.49
                N_cr = a_cr * N_comp[active]

                slen = np.sqrt((A * f_y) / N_cr)
                Phi = 0.5 * (1 + a_imp * (slen - 0.2) + slen**2)
                Chi = 1.0 / (Phi + np.sqrt(Phi**2 - slen**2))
                N_bRd = Chi * A * f_y

                util[active] = N_comp[active] / N_bRd

            return pd.Series(util, index=df.index, name="Util_BR")
        
        Util_BR_col = bucklingResFun(self.df_col,R1,R0,f_y) # Column
        Util_BR_brace = bucklingResFun(self.df_brace,R3,R2,f_y_brace) # Brace

        return Util_BR_col, Util_BR_brace

    # Interaction Force [6.3.3]
    def Util_IN(self):

        # Initialize Util_IN
        Util_IN_col = np.zeros(len(self.df_col))
        Util_IN_brace = np.zeros(len(self.df_brace))


        # Convert to Radii
        R0 = self.d0/2 - self.t0       # Column Inner Radius [mm]
        R1 = self.d0/2            # Column Outer Radius [mm]
        R2 = self.d1/2 - self.t1       # Brace Inner Radius  [mm]
        R3 = self.d1/2            # Brace Outer Radius  [mm]

        # Imperfection Factor
        a_imp = 0.49

        # Import Misc
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        # Open and Read Eigenvalue
        with open("Ansout/Eigenvalue1.txt") as f:
            eigenvalues = [float(line.strip()) for line in f if line.strip()]
        a_cr = next(v for v in eigenvalues if v > 0)

        def interaction(df_member, Ro, Ri, f_y):
            df = df_member.copy()
            util = np.zeros(len(df), dtype=float)

            Do = 2 * Ro
            Di = 2 * Ri
            A = math.pi * (Ro**2 - Ri**2)

            N_raw = df["NF"].to_numpy(dtype=float)
            N_comp = np.maximum(-N_raw, 0.0)
            active = N_comp > 0.0

            # Psi should come from the full member data, not filtered compression-only data
            M_start_y = df["My"].iloc[0]
            M_end_y = df["My"].iloc[-1]
            M1_y, M2_y = (M_start_y, M_end_y) if abs(M_start_y) >= abs(M_end_y) else (M_end_y, M_start_y)
            Psi_y = 0.0 if abs(M1_y) < 1e-12 else np.clip(M2_y / M1_y, -1.0, 1.0)

            M_start_z = df["Mz"].iloc[0]
            M_end_z = df["Mz"].iloc[-1]
            M1_z, M2_z = (M_start_z, M_end_z) if abs(M_start_z) >= abs(M_end_z) else (M_end_z, M_start_z)
            Psi_z = 0.0 if abs(M1_z) < 1e-12 else np.clip(M2_z / M1_z, -1.0, 1.0)

            if np.any(active):
                N_cr = a_cr * N_comp[active]

                slen = np.sqrt(A * f_y / N_cr)
                Phi = 0.5 * (1 + a_imp * (slen - 0.2) + slen**2)
                Chi = 1.0 / (Phi + np.sqrt(Phi**2 - slen**2))

                mu = (1 - N_comp[active] / N_cr) / (1 - Chi * N_comp[active] / N_cr)

                Cmy = 0.79 + 0.21 * Psi_y + 0.36 * (Psi_y - 0.33) * N_comp[active] / N_cr
                CmLT = 1.0
                Cmz = 0.79 + 0.21 * Psi_z + 0.36 * (Psi_z - 0.33) * N_comp[active] / N_cr

                k_yy = Cmy * CmLT * (mu / (1 - N_comp[active] / N_cr))
                k_yz = Cmz * CmLT * (mu / (1 - N_comp[active] / N_cr))

                N_Rk = A * f_y
                M_Rk = (math.pi * (Do**4 - Di**4)) / (32 * Do) * f_y

                My = df["My"].abs().to_numpy(dtype=float)
                Mz = df["Mz"].abs().to_numpy(dtype=float)

                util[active] = (
                    N_comp[active] / (Chi * N_Rk)
                    + k_yy * My[active] / M_Rk
                    + k_yz * Mz[active] / M_Rk
                )

            return pd.Series(util, index=df.index, name="Util_IN")

        Util_IN_col = interaction(self.df_col,R1,R0,f_y) # Column
        Util_IN_brace = interaction(self.df_brace,R3,R2,f_y_brace) # Brace

        return Util_IN_col, Util_IN_brace


    def Util_BS(self): 

        # Import variables
        d0 = self.d0
        t0 = self.t0
        d1 = self.d1
        t1 = self.t1

        # Convert to Radii
        R0 = d0/2 - t0       # Column Inner Radius [mm]
        R1 = d0/2            # Column Outer Radius [mm]
        R2 = d1/2 - t1       # Brace Inner Radius  [mm]
        R3 = d1/2            # Brace Outer Radius  [mm]

        # Import Misc
        f_y_brace = self.f_y_brace

        # Force and Length
        P = 200 * 9.82      # [N]
        #L = 350             # [mm]
        L = _brace_span_mm(self.var, self.Misc)  # [mm] first horizontal brace (same for all)
        print("Brace Length:", L)
        # Max Moment
        M = 1/8*P*L         # [Nmm]

        # Moment of Intertia and Area
        I = np.pi/64*(d1**4-(d1-2*t1)**4)   # [mm^4]
        A = np.pi/4*(d1**2-(d1-2*t1)**2)    # [mm^2]

        # Sectional Modulus
        W = I/(d1/2)                        # [mm^3]

        # Bending Stress
        sig_b = M/W                         # [N/mm^2]
        
        # Vertical Force
        V = P/2                             # [N]

        # Average and Max Shear Stress
        tau_avg = V/A                       # [N/mm^2]
        tau_max = 2*tau_avg                 # [N/mm^2]

        # Von Misses
        sig_vm = np.sqrt(sig_b**2+3*tau_max**2) # [N/mm^2]

        # Write to Utilization Ratio
        Util_BS_brace = sig_vm/f_y_brace                 # [Na] Brace

        return Util_BS_brace

    def Class_2(self): 



        # For this function we use it for a constraint, such that we always have a Class 2 section. or below according to
        # Eurocode 3 Section Table 5.2
        # And for the optimization scheme we need to implement the form:
        # c(x) >= 0 "Inequality Constraint"
        # Therefore we get:
        # 70*235/f_y-dw/tw >=0
        
        d0 = self.d0
        t0 = self.t0
        d1 = self.d1
        t1 = self.t1

        #Yield Strength of Columns
        f_y = self.f_y
        # Yield Strength of Braces
        f_y_brace = self.f_y_brace

        Util_Class_2 = np.zeros(2)
        Util_Class_2[0] = 70*(235/f_y)-(d0/t0) #column
        Util_Class_2[1] = 70*(235/f_y_brace)-(d1/t1) #brace
        return Util_Class_2


    def Eigenvalue_1(self):

        # We implement a constrain, such that the first positive eigenvalue is greater than or equal to 4.0
        # And for the optimization scheme we need to implement the form:
        # c(x) >= 0 "Inequality Constraint"
        # Therefore we get:
        # a_cr >= 4.0

        # Open and Read Eigenvalue
        with open("Ansout/Eigenvalue1.txt") as f:
            eigenvalues = [float(line.strip()) for line in f if line.strip()]
        a_cr = next(v for v in eigenvalues if v > 0)

        # Return Constraint Value
        # SLSQP inequality form: c(x) >= 0
        return a_cr - 4.0