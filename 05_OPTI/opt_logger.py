import os
from datetime import datetime
import numpy as np


class OptimizationLogger:               # This is the class that is used to log the optimization process
    def __init__(               # This is the constructor for the OptimizationLogger class
        self,
        x0,               # Initial design variables
        bounds,           # Bounds for the design variables
        constraint_func=None,           # Function to calculate the constraint values
        constraint_names=None,           # Names for the constraints
        method="SLSQP",           # Method for the optimization
        options=None,           # Options for the optimization
        save_folder="Optimization_Logs"           # Folder to save the log file
    ):               # This is the constructor for the OptimizationLogger class
        self.start_dt = datetime.now()               # Start date and time
        self.start_str = self.start_dt.strftime("%Y-%m-%d %H:%M:%S")               # Start date and time string
        self.timestamp_file = self.start_dt.strftime("%Y-%m-%d_%H-%M-%S")               # Timestamp for the log file

        os.makedirs(save_folder, exist_ok=True)               # Create the save folder if it doesn't exist
        self.txt_path = os.path.join(save_folder, f"optimization_log_{self.timestamp_file}.txt")               # Create the text path
        self.history_path = os.path.join(
        save_folder,               # Save folder
        f"objective_history_{self.timestamp_file}.csv"
    )               # Create the history path

        with open(self.history_path, "w", encoding="utf-8") as f:               # Open the history file
            f.write("type,index,objective\n")               # Write the header to the history file

        self.constraint_func = constraint_func               # Function to calculate the constraint values
        self.constraint_names = constraint_names               # Names for the constraints
        self.options = options if options is not None else {}               # Options for the optimization

        self.history = {               # History of the optimization process
            "run_start": self.start_str,
            "eval_num": [],               # Number of evaluations
            "eval_x": [],               # Design variables
            "eval_fun": [],               # Objective function
            "eval_time": [],               # Time of the evaluation
            "iter_num": [],               # Number of iterations
            "iter_x": [],               # Design variables
            "iter_fun": [],               # Objective function
            "iter_time": [],               # Time of the iteration
            "iter_step_norm": [],               # Step norm
            "eval_constraints": [],
            "iter_constraints": []               # Constraints
        }

        self.eval_counter = 0               # Counter for the evaluations
        self.prev_iter_x = None               # Previous iteration design variables
        self.latest_eval_x = None               # Latest evaluation design variables
        self.latest_eval_fun = None

        self._write_header(x0, bounds, method)               # Write the header

    def log_line(self, text="", mode="a"):               # This is the function that is used to log a line      
        print(text)               # Print the text
        with open(self.txt_path, mode, encoding="utf-8") as f:               # Open the text file
            f.write(text + "\n")               # Write the text to the file
    def log_objective_history(self, kind, index, objective):
        """
        Writes a compact line to the objective history file.

        Parameters
        ----------
        kind : str
            "evaluation" or "iteration"
        index : int
            Evaluation number or iteration number
        objective : float
            Objective function value
        """
        with open(self.history_path, "a", encoding="utf-8") as f:
            f.write(f"{kind},{index},{objective:.12f}\n")    

    def get_constraint_values(self, x):               # This is the function that is used to get the constraint values
        if self.constraint_func is None:               # If the constraint function is None, return an empty list
            return []               # Return an empty list

        vals = self.constraint_func(np.asarray(x, dtype=float))               # Calculate the constraint values

        if vals is None:
            return []               # Return an empty list  

        vals = np.atleast_1d(vals).astype(float).tolist()               # Convert the constraint values to a list
        return vals               # Return the constraint values

    def get_constraint_names(self, n_constraints):
        if self.constraint_names is not None:               # If the constraint names are not None, check if the length of the constraint names is equal to the number of constraints
            if len(self.constraint_names) != n_constraints:               # If the length of the constraint names is not equal to the number of constraints, raise an error
                raise ValueError(
                    f"constraint_names has length {len(self.constraint_names)}, "               # The constraint names have the wrong length
                    f"but there are {n_constraints} constraints."               # The number of constraints is wrong
                )               # Raise an error
            return self.constraint_names               # Return the constraint names

        return [f"c{i+1}" for i in range(n_constraints)]               # Return the constraint names

    def _write_header(self, x0, bounds, method):               # This is the function that is used to write the header
        self.log_line("=" * 80, mode="w")               # Write the header
        self.log_line("OPTIMIZATION LOG")               # Write the optimization log
        self.log_line("=" * 80)
        self.log_line(f"Run start date/time : {self.start_str}")               # Write the run start date and time
        self.log_line(f"Log file            : {self.txt_path}")               # Write the log file
        self.log_line("")
        self.log_line("Initial design variables:")               # Write the initial design variables   
        self.log_line(f"x0 = {np.array2string(np.asarray(x0), precision=6, separator=', ')}")               # Write the initial design variables
        self.log_line("")
        self.log_line("Bounds:")               # Write the bounds
        for i, bnd in enumerate(bounds):               # Write the bounds
            self.log_line(f"  x[{i}] in [{bnd[0]}, {bnd[1]}]")               # Write the bounds
        self.log_line("")
        self.log_line("Optimizer settings:")
        self.log_line(f"  Method  : {method}")               # Write the method
        for key, value in self.options.items():               # Write the optimization options
            self.log_line(f"  {key:<8}: {value}")               # Write the optimization options
        self.log_line("")

        initial_constraints = self.get_constraint_values(x0)               # Get the initial constraint values  
        if initial_constraints:
            names = self.get_constraint_names(len(initial_constraints))               # Get the initial constraint names
            self.log_line("Constraints detected:")               # Write the constraints detected
            for name in names:               # Write the constraint names
                self.log_line(f"  {name}(x) >= 0")               # Write the constraint names
            self.log_line("")

        self.log_line("=" * 80)               # Write the header
        self.log_line("")               # Write the empty line

    def log_evaluation(self, x, fun):               # This is the function that is used to log the evaluation process
        x = np.asarray(x, dtype=float)               # Convert the design variables to a numpy array
        fun = float(fun)               # Convert the objective function to a float

        self.eval_counter += 1               # Increment the evaluation counter
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")               # Get the current date and time
        constraints = self.get_constraint_values(x)               # Get the constraint values
        names = self.get_constraint_names(len(constraints)) if constraints else []               # Get the constraint names

        self.history["eval_num"].append(self.eval_counter)               # Append the evaluation number to the history
        self.history["eval_x"].append(x.tolist())
        self.history["eval_fun"].append(fun)               # Append the objective function to the history
        self.history["eval_time"].append(now_str)               # Append the current date and time to the history
        self.history["eval_constraints"].append(constraints)               # Append the constraint values to the history

        self.latest_eval_x = np.copy(x)               # Set the latest evaluation design variables to the current design variables
        self.latest_eval_fun = fun               # Set the latest evaluation function to the current objective function
        self.log_objective_history("evaluation", self.eval_counter, fun)

        self.log_line(f"[EVALUATION {self.eval_counter}]")               # Log the evaluation number
        self.log_line(f"  Date/time : {now_str}")               # Log the current date and time
        self.log_line(f"  x         : {np.array2string(x, precision=6, separator=', ')}")               # Log the design variables
        self.log_line(f"  Objective : {fun:.6f}")               # Log the objective function

        for name, value in zip(names, constraints):               # Log the constraint names and values
            self.log_line(f"  {name:<9}: {value:.6f}")               # Log the constraint names and values

        self.log_line("")               # Log the empty line

    def log_iteration(self, x):               # This is the function that is used to log the iteration process
        x = np.asarray(x, dtype=float)          # Convert the design variables to a numpy array
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")          # Get the current date and time
        iter_num = len(self.history["iter_num"]) + 1          # Get the number of iterations
        constraints = self.get_constraint_values(x)          # Get the constraint values
        names = self.get_constraint_names(len(constraints)) if constraints else []          # Get the constraint names

        if self.prev_iter_x is None:          # If the previous iteration is None, set the step norm to nan
            step_norm = np.nan          # Step norm is the norm of the difference between the current and previous design variables
        else:
            step_norm = float(np.linalg.norm(x - self.prev_iter_x))          # Step norm is the norm of the difference between the current and previous design variables

        if self.latest_eval_x is not None and np.allclose(x, self.latest_eval_x, rtol=1e-10, atol=1e-10):          # If the latest evaluation is not None and the current and latest evaluation are close, set the objective function to the latest evaluation function
            fun = self.latest_eval_fun          # Objective function is the value of the objective function
        else:
            fun = np.nan          # Objective function is nan
        if not np.isnan(fun):
            self.log_objective_history("iteration", iter_num, fun)
        self.history["iter_num"].append(iter_num)          # Append the iteration number to the history
        self.history["iter_x"].append(x.tolist())          # Append the design variables to the history
        self.history["iter_fun"].append(fun)          # Append the objective function to the history
        self.history["iter_time"].append(now_str)          # Append the current date and time to the history
        self.history["iter_step_norm"].append(step_norm)          # Append the step norm to the history
        self.history["iter_constraints"].append(constraints)          # Append the constraint values to the history

        self.prev_iter_x = np.copy(x)          # Set the previous iteration design variables to the current design variables

        self.log_line(f"[ITERATION {iter_num}]")          # Log the iteration number
        self.log_line(f"  Date/time : {now_str}")          # Log the current date and time
        self.log_line(f"  x         : {np.array2string(x, precision=6, separator=', ')}")          # Log the design variables
        self.log_line(f"  Objective : {fun:.6f}" if not np.isnan(fun) else "  Objective : nan")          # Log the objective function
        self.log_line(f"  Step norm : {step_norm:.6f}" if not np.isnan(step_norm) else "  Step norm : nan")          # Log the step norm

        for name, value in zip(names, constraints):          # Log the constraint names and values
            self.log_line(f"  {name:<9}: {value:.6f}")

        self.log_line("")

    def finalize(self, result):               # This is the function that is used to finalize the optimization process
        end_dt = datetime.now()
        runtime_sec = (end_dt - self.start_dt).total_seconds()

        final_constraints = self.get_constraint_values(result.x)          # Get the final constraint values
        names = self.get_constraint_names(len(final_constraints)) if final_constraints else []          # Get the final constraint names

        self.log_line("=" * 80)
        self.log_line("FINAL RESULT")          # Log the final result
        self.log_line("=" * 80)
        self.log_line(f"Run end date/time   : {end_dt.strftime('%Y-%m-%d %H:%M:%S')}")          # Log the final date and time   
        self.log_line(f"Total runtime [s]   : {runtime_sec:.3f}")          # Log the total runtime
        self.log_line(f"Success             : {result.success}")          # Log the success
        self.log_line(f"Status              : {result.status}")          # Log the status
        self.log_line(f"Message             : {result.message}")          # Log the message
        self.log_line(f"Iterations          : {getattr(result, 'nit', 'N/A')}")          # Log the number of iterations
        self.log_line(f"Function evaluations: {getattr(result, 'nfev', 'N/A')}")          # Log the number of function evaluations
        self.log_line(f"Optimal x           : {np.array2string(np.asarray(result.x), precision=6, separator=', ')}")          # Log the optimal design variables
        self.log_line(f"Optimal objective   : {float(result.fun):.6f}")          # Log the optimal objective function
        self.log_line("")

        if final_constraints:          # If the final constraints are not None, log the final constraint values
            self.log_line("Final constraint values:")          # Log the final constraint values
            for name, value in zip(names, final_constraints):          # Log the final constraint names and values
                self.log_line(f"  {name:<9}= {value:.6f}")          # Log the final constraint names and values
            self.log_line("")

        self.log_line("=" * 80) 
        self.log_line("END OF LOG")
        self.log_line("=" * 80)