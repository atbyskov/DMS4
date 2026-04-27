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


# Function to calculate the span of horizontal brace - It works, but i (Christian) don't understand this
def _brace_span_mm(var, misc, tol=1e-3):
    co = SW.import_SW(str(Path("IGS") / misc["SW_filename"]))

    # If a target radius is provided, move the column-connected points
    # radially in the x-z plane while keeping y unchanged
    rad = var.get("rad")
    if isinstance(rad, dict):
        rad = rad.get("value")
    if rad is not None:
        # Find vertical column footprints in the x-z plane
        col = {
            (x1, z1)
            for x1, y1, z1, x2, y2, z2 in co
            if abs(x1 - x2) < tol and abs(z1 - z2) < tol
        }

        # Map each footprint to the same direction but new radius = rad
        ln = {}
        for x, z in col:
            h = math.hypot(x, z)
            if h < tol:
                ln[(x, z)] = (x, z)  # avoid dividing by ~0
            else:
                ln[(x, z)] = (x * rad / h, z * rad / h)

        def adj(x, y, z):
            if (x, z) in ln:
                x_new, z_new = ln[(x, z)]
                return (x_new, y, z_new)
            return (x, y, z)

        co = [
            (*adj(x1, y1, z1), *adj(x2, y2, z2))
            for x1, y1, z1, x2, y2, z2 in co
        ]

    # Find the first member that is horizontal:
    #   - same y-coordinate
    #   - but different x/z coordinates
    for x1, y1, z1, x2, y2, z2 in co:
        is_horizontal = abs(y2 - y1) < tol
        has_horizontal_span = abs(x2 - x1) >= tol or abs(z2 - z1) >= tol

        if is_horizontal and has_horizontal_span:
            return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

class PostProcessor:

    def __init__(self,var,Misc,opti_settings):
        # Store inputs
        self.var = var
        self.Misc = Misc
        self.opti_settings = opti_settings

        self.n = int(opti_settings.get("n_mast_segments"))
        self.mast_height = float(opti_settings.get("mast_segment_height"))
        self.multi_size_columns = bool(opti_settings.get("multi_size_columns"))
        self.multi_size_braces = bool(opti_settings.get("multi_size_braces"))
        self.brace_split = bool(opti_settings.get("brace_split", False))
        self.E_mod = Misc["E_mod"]
        self.f_y = Misc["f_y"]
        self.f_y_brace = Misc["f_y_brace"]

        # Build per-segment section radii for columns and braces
        self.column_R0, self.column_R1 = self._build_section_radii("d0", "t0", self.multi_size_columns)
        self.column_d0_t0 = self._build_section_dimensions("d0", "t0", self.multi_size_columns)
        
        # Handle braces based on split setting
        if self.brace_split:
            # Horizontal braces
            self.brace_horiz_R2, self.brace_horiz_R3 = self._build_section_radii("d1_h", "t1_h", self.multi_size_braces)
            # Cross braces
            self.brace_cross_R2, self.brace_cross_R3 = self._build_section_radii("d1_c", "t1_c", self.multi_size_braces)
            # Build per-segment diameter/thickness data for class checks
            self.brace_horiz_d1_t1 = self._build_section_dimensions("d1_h", "t1_h", self.multi_size_braces)
            self.brace_cross_d1_t1 = self._build_section_dimensions("d1_c", "t1_c", self.multi_size_braces)
            self.brace_d1_t1 = self.brace_horiz_d1_t1 + self.brace_cross_d1_t1
        else:
            # Combined braces (same size for both)
            self.brace_R2, self.brace_R3 = self._build_section_radii("d1", "t1", self.multi_size_braces)
            self.brace_d1_t1 = self._build_section_dimensions("d1", "t1", self.multi_size_braces)

        # Load nonlinear results and separate columns from braces
        self.df_nonlin = self.read_forces("Ansout/APDL_Nonlin_Internal.txt")
        self.df_col = self.df_nonlin[self.df_nonlin["Member"].str.startswith("ColMember")].copy()
        
        if self.brace_split:
            self.df_horiz = self.df_nonlin[self.df_nonlin["Member"].str.startswith("HorizMember")].copy()
            self.df_cross = self.df_nonlin[self.df_nonlin["Member"].str.startswith("CrossMember")].copy()
        else:
            # If not split, combine all braces
            self.df_brace = self.df_nonlin[
                (self.df_nonlin["Member"].str.startswith("HorizMember")) |
                (self.df_nonlin["Member"].str.startswith("CrossMember"))
            ].copy()

        self._annotate_member_sections()
        self.df_brace = self._combined_brace_dataframe()

    # Extract numeric value from the var dictionary, supporting nested dict shape
    def _get_value(self, key):
        value = self.var.get(key)
        if isinstance(value, dict):
            return float(value["value"])
        return float(value)

    # Build inner and outer radii for one or multiple section sizes
    def _build_section_radii(self, d_key, t_key, multi_size):
        if multi_size:
            inner = []
            outer = []
            for i in range(1, self.n + 1):
                d = self._get_value(f"{d_key}_{i}")
                t = self._get_value(f"{t_key}_{i}")
                outer.append(d / 2.0)
                inner.append(d / 2.0 - t)
            return inner, outer

        d = self._get_value(d_key)
        t = self._get_value(t_key)
        return [d / 2.0 - t], [d / 2.0]

    # Build diameter/thickness tuples for class checks and brace sizing
    def _build_section_dimensions(self, d_key, t_key, multi_size):
        if multi_size:
            return [
                (self._get_value(f"{d_key}_{i}"), self._get_value(f"{t_key}_{i}"))
                for i in range(1, self.n + 1)
            ]
        return [(self._get_value(d_key), self._get_value(t_key))]

    # Convert Y location to mast segment indices for section lookup
    def _segment_series(self, y_series):
        segments = np.ceil(y_series / self.mast_height).astype(int)
        segments = np.maximum(segments, 1)
        segments = np.minimum(segments, self.n)
        return segments

    # Assign section radii to each member row based on its segment position
    def _annotate_member_sections(self):
        if not self.df_col.empty:
            self.df_col["segment"] = self._segment_series(self.df_col["Y_LOC"])
            if len(self.column_R0) > 1:
                self.df_col["R0"] = self.df_col["segment"].map(lambda s: self.column_R0[s - 1])
                self.df_col["R1"] = self.df_col["segment"].map(lambda s: self.column_R1[s - 1])
            else:
                self.df_col["R0"] = self.column_R0[0]
                self.df_col["R1"] = self.column_R1[0]

        if self.brace_split:
            # Handle separate horizontal and cross braces
            if not self.df_horiz.empty:
                self.df_horiz["segment"] = self._segment_series(self.df_horiz["Y_LOC"])
                if len(self.brace_horiz_R2) > 1:
                    self.df_horiz["R2"] = self.df_horiz["segment"].map(lambda s: self.brace_horiz_R2[s - 1])
                    self.df_horiz["R3"] = self.df_horiz["segment"].map(lambda s: self.brace_horiz_R3[s - 1])
                else:
                    self.df_horiz["R2"] = self.brace_horiz_R2[0]
                    self.df_horiz["R3"] = self.brace_horiz_R3[0]

            if not self.df_cross.empty:
                self.df_cross["segment"] = self._segment_series(self.df_cross["Y_LOC"])
                if len(self.brace_cross_R2) > 1:
                    self.df_cross["R2"] = self.df_cross["segment"].map(lambda s: self.brace_cross_R2[s - 1])
                    self.df_cross["R3"] = self.df_cross["segment"].map(lambda s: self.brace_cross_R3[s - 1])
                else:
                    self.df_cross["R2"] = self.brace_cross_R2[0]
                    self.df_cross["R3"] = self.brace_cross_R3[0]
        else:
            # Handle combined braces
            if not self.df_brace.empty:
                self.df_brace["segment"] = self._segment_series(self.df_brace["Y_LOC"])
                if len(self.brace_R2) > 1:
                    self.df_brace["R2"] = self.df_brace["segment"].map(lambda s: self.brace_R2[s - 1])
                    self.df_brace["R3"] = self.df_brace["segment"].map(lambda s: self.brace_R3[s - 1])
                else:
                    self.df_brace["R2"] = self.brace_R2[0]
                    self.df_brace["R3"] = self.brace_R3[0]

    def _combined_brace_dataframe(self):
        if self.brace_split:
            frames = [df for df in (self.df_horiz, self.df_cross) if not df.empty]
            if not frames:
                return pd.DataFrame()
            # Combine horizontal and cross braces into one continuous DataFrame.
            # Use a fresh index to avoid separate original index spaces.

            return pd.concat(frames, ignore_index=True)
        return self.df_brace

    def _brace_utilization(self, func):
        if self.brace_split:
            values = [
                func(df, "R3", "R2", self.f_y_brace)
                for df in (self.df_horiz, self.df_cross)
                if not df.empty
            ]
            if values:
                return pd.concat(values, ignore_index=False)
            return pd.Series(dtype=float)
        return func(self.df_brace, "R3", "R2", self.f_y_brace)

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
                #"Util_BNS": self.Util_BNS(),
                "Util_BR": self.Util_BR(),
                "Util_IN": self.Util_IN(),
                "Util_BS_sig": self.Util_BS()[0], # stress
                "Util_BS_defl": self.Util_BS()[1], # deflection
            }

            return util_data
    
    # Local Buckling  [Timeshenko p. 458]
    def Util_LB(self):
        # Import E_mod 
        E_mod = self.E_mod

        # Poissons Ratio
        v = 0.3

        def LocalBucklingFun(df_member, Ro_col, Ri_col):
            df = df_member.copy()

            P_cr = (2 * math.pi * E_mod * (df[Ro_col] - df[Ri_col])**2) / math.sqrt(3 * (1 - v**2))

            N_raw = df["NF"].to_numpy(dtype=float)
            N_comp = np.maximum(-N_raw, 0.0)   # only compression counts

            return pd.Series(N_comp / P_cr, index=df.index, name="Util_LB")

        Util_LB_col = LocalBucklingFun(self.df_col, "R1", "R0")
        Util_LB_brace = LocalBucklingFun(self.df_brace, "R3", "R2")

        return Util_LB_col, Util_LB_brace

    # Normal Force [6.2.4]                
    def Util_NF(self):

        # Import f_y
        f_y = self.f_y

        def NormalForceFun(df_member, Ro_col, Ri_col, f_y):
            # Area
            A = np.pi * ((df_member[Ro_col]**2) - (df_member[Ri_col]**2))

            # Compression/Tension Utilization (1993-1-1 [6.2.3] and [6.2.4])
            df_member["Util_NF"] = df_member["NF"].abs() / (A * f_y)

            # Utilization Ratio [1993-1-1 [6.2.4] - p. 49 (6.9)]
            return df_member["Util_NF"]

        Util_NF_col = NormalForceFun(self.df_col, "R1", "R0", f_y) # Column
        Util_NF_brace = self._brace_utilization(NormalForceFun) # Brace
        
        return Util_NF_col, Util_NF_brace

    # Shear Force [6.2.6]
    def Util_S(self):

        # Import Misc
        f_y = self.f_y

        def shearFun(df_member, Ro_col, Ri_col, f_y):
            A = np.pi * ((df_member[Ro_col]**2) - (df_member[Ri_col]**2))

            # Av Area [6.2.6 (3) with option (g)]
            Av = 2/math.pi * A

            # Design Plastic Shear Resistance [6.2.6 (6.18)]
            V_cRd = Av * (f_y / math.sqrt(3))

            # Utilization Ratio [6.2.6 (6.17)]
            # We have two directions (y,z), so we take whatever value in each element, that is highest
            df_member["Util_S"] = df_member[["Vy","Vz"]].abs().max(axis=1) / V_cRd

            return df_member["Util_S"]

        Util_S_col = shearFun(self.df_col, "R1", "R0", f_y) # Column
        Util_S_brace = self._brace_utilization(shearFun) # Brace

        return Util_S_col, Util_S_brace

    # Torsion [6.2.7]
    def Util_T(self):

        # Import Misc
        f_y = self.f_y

        def torsionFun(df_member, Ro, Ri, f_y):
            D0 = df_member[Ro] * 2 # Outer Diameter
            Di = df_member[Ri] * 2 # Inner Diameter

            # Torsion Design Resitance [6.2.7 (6.23)]
            T_Rd = (math.pi/16 * (D0**4 - Di**4) / D0) * f_y / math.sqrt(3)

            # Torsion Check [6.2.7 (6.23)]
            df_member["Util_T"] = df_member["T"].abs() / T_Rd

            # Return Utilization Ratio
            return df_member["Util_T"]

        Util_T_col = torsionFun(self.df_col, "R1", "R0", f_y) # Column
        Util_T_brace = self._brace_utilization(torsionFun) # Brace

        return Util_T_col, Util_T_brace

    # Bending, Normal and Shear [6.2.9]
    def Util_BNS(self):
        # Import f_y
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        def bnsFun(df_member, Ro, Ri, f_y):
            D0 = df_member[Ro] * 2
            Di = df_member[Ri] * 2
            t = (D0 - Di) / 2

            # Area
            A = np.pi * ((df_member[Ro]**2) - (df_member[Ri]**2))
            Aw = (A - 2 * D0 * t) / A
            Aw = np.clip(Aw, 0, 0.5)

            # Design Moment Resitance
            M_Rd = (D0**3 - Di**3) / 6 * f_y

            # Design Normal Force Resistance
            N_Rd = A * f_y

            # Forces (Here read as values)
            N = df_member["NF"].abs()
            My = df_member["My"].abs()
            Mz = df_member["Mz"].abs()

            # Design normal force
            N_Rd = A*f_y 
            n = N/N_Rd

            # Fraction to be used in Utilization Ratio
            red = np.clip(1 - (N / N_Rd)**1.7,0.1,1.0)

            # Utilization Ratio
            df_member["Util_BNS"] = (My / (M_Rd * red))**2 + (Mz / (M_Rd * red))**2

            return df_member["Util_BNS"]

        Util_BNS_col = bnsFun(self.df_col, "R1", "R0", f_y) # Column
        Util_BNS_brace = bnsFun(self.df_brace, "R3", "R2", f_y_brace) # Brace

        return Util_BNS_col, Util_BNS_brace

    # Buckling Resistance [6.3.1]
    def Util_BR(self):

        # Open and Read Eigenvalue
        with open("Ansout/Eigenvalue1.txt") as f:
            eigenvalues = [float(line.strip()) for line in f if line.strip()]
        a_cr = next(v for v in eigenvalues if v > 0)

        # Import Misc
        f_y = self.f_y
        f_y_brace = self.f_y_brace

        def bucklingResFun(df_member, Ro, Ri, f_y):
            df = df_member.copy()
            util = np.zeros(len(df), dtype=float)

            A = np.pi * (df[Ro]**2 - df[Ri]**2)

            N_raw = df["NF"].to_numpy(dtype=float)
            N_comp = np.maximum(-N_raw, 0.0)

            active = N_comp > 0.0

            if np.any(active):
                a_imp = 0.49
                N_cr = a_cr * N_comp[active]

                slen = np.sqrt((A[active] * f_y) / N_cr)
                Phi = 0.5 * (1 + a_imp * (slen - 0.2) + slen**2)
                Chi = 1.0 / (Phi + np.sqrt(Phi**2 - slen**2))
                N_bRd = Chi * A[active] * f_y

                util[active] = N_comp[active] / N_bRd

            return pd.Series(util, index=df.index, name="Util_BR")
        
        Util_BR_col = bucklingResFun(self.df_col, "R1", "R0", f_y) # Column
        Util_BR_brace = bucklingResFun(self.df_brace, "R3", "R2", f_y_brace) # Brace

        return Util_BR_col, Util_BR_brace

    # Interaction Force [6.3.3]
    def Util_IN(self):

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
            df = df_member.copy().reset_index(drop=True)
            util = np.zeros(len(df), dtype=float)

            # Geometry (row-wise)
            Do = 2 * df[Ro]
            Di = 2 * df[Ri]
            A = math.pi * (df[Ro]**2 - df[Ri]**2)

            N_raw = df["NF"].to_numpy(dtype=float)
            N_comp = np.maximum(-N_raw, 0.0)

            # Loop per structural member
            for member, g in df.groupby("Member"):
                idx = g.index.to_numpy()

                N_c = N_comp[idx]
                active = N_c > 0.0
                if not np.any(active):
                    continue

                # ---- END MOMENTS (per member) ---------------------------------
                My_start, My_end = g["My"].iloc[0], g["My"].iloc[-1]
                M1y, M2y = (My_start, My_end) if abs(My_start) >= abs(My_end) else (My_end, My_start)
                Psi_y = 0.0 if abs(M1y) < 1e-12 else np.clip(M2y / M1y, -1.0, 1.0)

                Mz_start, Mz_end = g["Mz"].iloc[0], g["Mz"].iloc[-1]
                M1z, M2z = (Mz_start, Mz_end) if abs(Mz_start) >= abs(Mz_end) else (Mz_end, Mz_start)
                Psi_z = 0.0 if abs(M1z) < 1e-12 else np.clip(M2z / M1z, -1.0, 1.0)

                # ---- BUCKLING / INTERACTION ----------------------------------- (Calculations per element)
                N_cr = a_cr * N_c[active]

                slen = np.sqrt(A[idx][active] * f_y / N_cr)
                Phi = 0.5 * (1 + a_imp * (slen - 0.2) + slen**2)
                Chi = 1.0 / (Phi + np.sqrt(Phi**2 - slen**2))

                mu = (1 - N_c[active] / N_cr) / (1 - Chi * N_c[active] / N_cr)

                Cmy = 0.79 + 0.21 * Psi_y + 0.36 * (Psi_y - 0.33) * N_c[active] / N_cr
                Cmz = 0.79 + 0.21 * Psi_z + 0.36 * (Psi_z - 0.33) * N_c[active] / N_cr
                CmLT = 1.0

                k_yy = Cmy * CmLT * (mu / (1 - N_c[active] / N_cr))
                k_yz = Cmz * CmLT * (mu / (1 - N_c[active] / N_cr))

                N_Rk = A[idx][active] * f_y
                M_Rk = (math.pi * (Do[idx][active]**4 - Di[idx][active]**4)) / (32 * Do[idx][active]) * f_y

                My = df.loc[idx, "My"].abs().to_numpy()[active]
                Mz = df.loc[idx, "Mz"].abs().to_numpy()[active]

                util_idx = idx[active]
                util[util_idx] = (
                    N_c[active] / (Chi * N_Rk)
                    + k_yy * My / M_Rk
                    + k_yz * Mz / M_Rk
                )

            return pd.Series(util, index=df.index, name="Util_IN")

        Util_IN_col = interaction(self.df_col, "R1", "R0", f_y) # Column
        Util_IN_brace = interaction(self.df_brace, "R3", "R2", f_y_brace) # Brace

        return Util_IN_col, Util_IN_brace

    # Brace Step [Machine Design]
    # Brace-Step
    def Util_BS(self): 

        brace_sections = self.brace_horiz_d1_t1 if self.brace_split else self.brace_d1_t1

        # Import Misc
        f_y_brace = self.f_y_brace
        E_mod = self.E_mod

        # Force and Length
        P = 150 * 9.82      # [N] 150 kg from standard
        L = _brace_span_mm(self.var, self.Misc)  # [mm] first horizontal brace (same for all)
        M = 1/8 * P * L     # Max Moment [Nmm]

        # Calculate utilization for all horizontal brace sections
        util_values_sig = []
        util_values_defl = []
        
        for d1, t1 in brace_sections:
            # Moment of Inertia and Area
            I = np.pi / 64 * (d1**4 - (d1 - 2 * t1)**4) # [mm^4]
            A = np.pi / 4 * (d1**2 - (d1 - 2 * t1)**2)  # [mm^2]

            # Sectional Modulus
            W = I / (d1 / 2)                            # [mm^3]

            # Bending Stress
            sig_b = M / W                               # [N/mm^2]

            # Vertical Force
            V = P / 2                                   # [N]

            # Average and Max Shear Stress
            tau_avg = V / A                             # [N/mm^2]
            tau_max = 2 * tau_avg                       # [N/mm^2]

            # Von Mises
            sig_vm = np.sqrt(sig_b**2 + 3 * tau_max**2) # [N/mm^2]

            util_values_sig.append(sig_vm / f_y_brace)      # [Na] Brace

            # Deflection (u_max = (1/192) * ((P * L^3) / (E * I)))
            u_max = (P * L**3) / (192 * E_mod * I)         # [mm]

            util_values_defl.append(u_max / (L / 200))          # [Na] Brace Max deflection from standard
        # Return as Series to allow extract_max to work with multiple values
        return pd.Series(util_values_sig, name="Util_BS_sig"), pd.Series(util_values_defl, name="Util_BS_defl")

    def Class_2(self): 

        # For this function we use it for a constraint, such that we always have a Class 2 section or below according to
        # Eurocode 3 Section Table 5.2
        # And for the optimization scheme we need to implement the form:
        # c(x) >= 0 "Inequality Constraint"
        # Therefore we get:
        # 70*235/f_y-dw/tw >=0

        f_y = self.f_y
        f_y_brace = self.f_y_brace

        if self.multi_size_columns:
            col_ratio = max(d / t for d, t in self.column_d0_t0)
        else:
            col_ratio = self.column_d0_t0[0][0] / self.column_d0_t0[0][1]

        if self.multi_size_braces:
            brace_ratio = max(d / t for d, t in self.brace_d1_t1)
        else:
            brace_ratio = self.brace_d1_t1[0][0] / self.brace_d1_t1[0][1]

        Util_Class_2 = np.zeros(2)
        Util_Class_2[0] = 70 * (235 / f_y) - col_ratio
        Util_Class_2[1] = 70 * (235 / f_y_brace) - brace_ratio
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
