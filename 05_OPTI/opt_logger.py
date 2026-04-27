import os
from datetime import datetime
from collections.abc import Sequence

import numpy as np

# The reason we also have log_iteration function is because the Evaluation functions from the optimization scheme is
# used to estimate the gradidents, and then from the gradients, the step size is calculated, and then for every new
# Step size, there is an iteration happening, that's why we both track evaluations, and iterations.


class OptimizationLogger:
    """Text log of the run, plus a small CSV of each objective evaluation (FEA call)."""

    def __init__(
        self,
        x0,
        bounds,
        method="SLSQP",
        options=None,
        var_names = None,
        save_folder="Optimization_Logs",
    ):
        
        self.options = options or {} # This line stores the optimization options.
        self.eval_counter = 0 # This line stores the evaluation counter.
        self.iter_counter = 0 # This line stores the iteration counte
        self.start_dt = datetime.now() #Get the time once now, and then we store it in different formats underneath
        self.start_str = self.start_dt.strftime("%Y-%m-%d %H:%M:%S") # THis is a human readable format.
        ts = self.start_dt.strftime("%Y-%m-%d_%H-%M-%S") # And this one is the timestamp for the file name, so for example optimization_log_2026-03-31_31_10-23-25.txt

        os.makedirs(save_folder, exist_ok=True) # This line creates the folder "Optimization_Logs" if it doesn't exist.
        self.txt_path = os.path.join(save_folder, f"optimization_log_{ts}.txt") # This line creates the path to the optimization log file.
        self.csv_path = os.path.join(save_folder, f"objective_history_{ts}.csv") # This line creates the path to the objective history file.
        self._write_header(x0, bounds, method) # This line writes the header to the optimization log file.

        n_vars = len(x0)
        var_names = var_names or [f"x{i}" for i in range(n_vars)]
        x_headers = ",".join(f"x{i}" for i in range(n_vars))

        # Header
        with open(self.csv_path, "w", encoding="utf-8") as f: # This line writes the header to the objective history file.
            f.write("Optimization Data\n \n")
            f.write(f"Start time: {self.start_str} \n \n")
            f.write(f"Solver Options:\n")
            f.write(f"Method: {method} \n")
            for key, value in self.options.items():
                f.write(f"# {key}: {value} \n")

            f.write("Design Variables and Bounds \n")
            f.write("# x0 = [" + ", ".join(f"{v:.2f}" for v in x0) + "]\n")
            for i, (name,bnd) in enumerate(zip(var_names, bounds)):
                f.write(f" x[{i}] ({name}) bounds = [{bnd[0]}, {bnd[1]}] \n")

            f.write("\n\n")
            f.write(f"iteration,eval_index,objective,{x_headers},g_max,v_agg\n")

        
        

    def log_line(self, text="", mode="a"): # Help Function to help write to the Optimization_Log folder files. It is instead of using print statements.
        import sys
        print(text, flush=True)  # Flush stdout immediately
        with open(self.txt_path, mode, encoding="utf-8") as f:
            f.write(text + "\n")
            f.flush()  # Ensure file is written to disk

    def _write_header(self, x0, bounds, method): # Function to write the initial Header in the Optimization_Log folder files., using the log_line function.
        self.log_line("=" * 80, mode="w")
        self.log_line("OPTIMIZATION LOG")
        self.log_line("=" * 80)
        self.log_line(f"Run start date/time : {self.start_str}")
        self.log_line(f"Log file            : {self.txt_path}")
        self.log_line(f"Objective CSV       : {self.csv_path}")
        self.log_line("")
        self.log_line("Initial design variables:")
        self.log_line(
            f"x0 = {np.array2string(np.asarray(x0), precision=6, separator=', ')}"
        )
        self.log_line("")
        self.log_line("Bounds:")
        for i, bnd in enumerate(bounds):
            self.log_line(f"  x[{i}] in [{bnd[0]}, {bnd[1]}]")
        self.log_line("")
        self.log_line("Optimizer settings:")
        self.log_line(f"  Method  : {method}")
        for key, value in self.options.items():
            self.log_line(f"  {key:<8}: {value}")
        self.log_line("")
        self.log_line("=" * 80)
        self.log_line("")

    def log_evaluation(self, x, fun, v_agg=None,g_max=None,p_value=None): # When called, it will log the evaluation of the objective function and design variables, and write it inside the Optimization_Log folder files.
        x = np.asarray(x, dtype=float) # Shows the design variables as a numpy array.
        fun = float(fun) # Shows the objective function as a float.
        self.eval_counter += 1 # Increments the evaluation counter.
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S") # Shows the current date and time in a human readable format.

        #NOTE: Relative time that it takes per evaluation is not implemented yet.
        #Implement, such that it is relative time that it takes per evaluation ############################################

        
# Convert x into CSV format
        x_csv = ",".join(f"{v:.2f}" for v in x)

        v_agg_str = f"{v_agg:.4f}" if v_agg is not None else "NA"
        g_max_str = f"{g_max:.4f}" if g_max is not None else "NA"

        with open(self.csv_path, "a", encoding="utf-8") as f:
            f.write(
                f"{self.iter_counter},"
                f"{self.eval_counter},"
                f"{fun:.2f},"
                f"{x_csv},"
                f"{g_max_str},"
                f"{v_agg_str}\n"
            )




        self.log_line(f"[EVALUATION {self.eval_counter}]") # Writes the evaluation number to the optimization log file.
        self.log_line(f"  Date/time : {now}") # Writes the current date and time to the optimization log file.
        x_str = ", ".join(f"{v:.2f}" for v in x) ######### Maybe doesnt work
        self.log_line(f"  x         : [{x_str}]")
        # Writes the design variables to the optimization log file.
        self.log_line(f"  Objective : {fun:.3f}") # Writes the objective function to the optimization log file.
        self.log_line("") # Writes an empty line to the optimization log file.

    def log_iteration(self, x): # When called, it will log the iteration of the design variables, and write it inside the Optimization_Log folder files.
        self.iter_counter += 1 # Increments the iteration counter.
        x = np.asarray(x, dtype=float) # Shows the design variables as a numpy array.
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S") # Shows the current date and time in a human readable format.
        self.log_line(f"[ITERATION {self.iter_counter}]") # Writes the iteration number to the optimization log file.
        self.log_line(f"  Date/time : {now}") # Writes the current date and time to the optimization log file.
        self.log_line(f"  x         : {np.array2string(x, precision=6, separator=', ')}") # Writes the design variables to the optimization log file.
        self.log_line("") # Writes an empty line to the optimization log file.

    def finalize(self, result): # When called, it will log the final result of the optimization, and write it inside the Optimization_Log folder files. This is only used, when convergence is met.
        end = datetime.now() # Gets the current date and time.
        dt_sec = (end - self.start_dt).total_seconds() # Calculates the total runtime of the optimization.
        xopt = np.asarray(result.x, dtype=float) # Shows the optimal design variables as a numpy array.

        self.log_line("=" * 80) # ======= lines
        self.log_line("FINAL RESULT") 
        self.log_line("=" * 80)
        self.log_line(f"Run end date/time   : {end.strftime('%Y-%m-%d %H:%M:%S')}") # Writes the end date and time to the optimization log file.
        self.log_line(f"Total runtime [s]   : {dt_sec:.3f}") # Writes the total runtime to the optimization log file.
        self.log_line(f"Success             : {result.success}") # Writes the success status to the optimization log file.
        self.log_line(f"Status              : {result.status}") # Writes the status to the optimization log file.
        self.log_line(f"Message             : {result.message}") # Writes the message to the optimization log file.
        self.log_line(f"Iterations          : {result.nit}") # Writes the number of iterations to the optimization log file.
        self.log_line(f"Function evaluations: {result.nfev}") # Writes the number of function evaluations to the optimization log file.
        self.log_line(f"Optimal x           : {np.array2string(xopt, precision=6, separator=', ')}") # Writes the optimal design variables to the optimization log file.
        self.log_line(f"Optimal objective   : {float(result.fun):.6f}") # Writes the optimal objective function to the optimization log file.
        self.log_line("") # Writes an empty line to the optimization log file.
        self.log_line("=" * 80) # ======= lines
        self.log_line("END OF LOG") 
        self.log_line("=" * 80) # ======= lines

    
    def _extract_max(self, val):
        if val is None:
            return None
        arr = np.asarray(val, dtype=float)
        if arr.size == 0:
            return None
        return np.nanmax(arr)


    def log_utilization(self, util_dict):

        for key, util in util_dict.items():

            col_val = None
            brc_val = None

            # Column / Brace case
            if isinstance(util, Sequence) and len(util) == 2:
                col_val = self._extract_max(util[0])
                brc_val = self._extract_max(util[1])

            # Brace-only case (e.g. Util_BS)
            else:
                brc_val = self._extract_max(util)

            col_str = f"{col_val:8.3f}" if col_val is not None else "   N/A  "
            brc_str = f"{brc_val:8.3f}" if brc_val is not None else "   N/A  "

            self.log_line(
                f"{key:10s}  Column: {col_str}  Brace: {brc_str}"
            )

