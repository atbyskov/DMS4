import os
from datetime import datetime
import numpy as np


class OptimizationLogger:
    def __init__(
        self,
        x0,
        bounds,
        constraint_func=None,
        constraint_names=None,
        method="SLSQP",
        options=None,
        save_folder="Optimization_Logs"
    ):
        self.start_dt = datetime.now()
        self.start_str = self.start_dt.strftime("%Y-%m-%d %H:%M:%S")
        self.timestamp_file = self.start_dt.strftime("%Y-%m-%d_%H-%M-%S")

        os.makedirs(save_folder, exist_ok=True)
        self.txt_path = os.path.join(save_folder, f"optimization_log_{self.timestamp_file}.txt")

        self.constraint_func = constraint_func
        self.constraint_names = constraint_names
        self.options = options if options is not None else {}

        self.history = {
            "run_start": self.start_str,
            "eval_num": [],
            "eval_x": [],
            "eval_fun": [],
            "eval_time": [],
            "iter_num": [],
            "iter_x": [],
            "iter_fun": [],
            "iter_time": [],
            "iter_step_norm": [],
            "eval_constraints": [],
            "iter_constraints": []
        }

        self.eval_counter = 0
        self.prev_iter_x = None
        self.latest_eval_x = None
        self.latest_eval_fun = None

        self._write_header(x0, bounds, method)

    def log_line(self, text="", mode="a"):
        print(text)
        with open(self.txt_path, mode, encoding="utf-8") as f:
            f.write(text + "\n")

    def get_constraint_values(self, x):
        if self.constraint_func is None:
            return []

        vals = self.constraint_func(np.asarray(x, dtype=float))

        if vals is None:
            return []

        vals = np.atleast_1d(vals).astype(float).tolist()
        return vals

    def get_constraint_names(self, n_constraints):
        if self.constraint_names is not None:
            if len(self.constraint_names) != n_constraints:
                raise ValueError(
                    f"constraint_names has length {len(self.constraint_names)}, "
                    f"but there are {n_constraints} constraints."
                )
            return self.constraint_names

        return [f"c{i+1}" for i in range(n_constraints)]

    def _write_header(self, x0, bounds, method):
        self.log_line("=" * 80, mode="w")
        self.log_line("OPTIMIZATION LOG")
        self.log_line("=" * 80)
        self.log_line(f"Run start date/time : {self.start_str}")
        self.log_line(f"Log file            : {self.txt_path}")
        self.log_line("")
        self.log_line("Initial design variables:")
        self.log_line(f"x0 = {np.array2string(np.asarray(x0), precision=6, separator=', ')}")
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

        initial_constraints = self.get_constraint_values(x0)
        if initial_constraints:
            names = self.get_constraint_names(len(initial_constraints))
            self.log_line("Constraints detected:")
            for name in names:
                self.log_line(f"  {name}(x) >= 0")
            self.log_line("")

        self.log_line("=" * 80)
        self.log_line("")

    def log_evaluation(self, x, fun):
        x = np.asarray(x, dtype=float)
        fun = float(fun)

        self.eval_counter += 1
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        constraints = self.get_constraint_values(x)
        names = self.get_constraint_names(len(constraints)) if constraints else []

        self.history["eval_num"].append(self.eval_counter)
        self.history["eval_x"].append(x.tolist())
        self.history["eval_fun"].append(fun)
        self.history["eval_time"].append(now_str)
        self.history["eval_constraints"].append(constraints)

        self.latest_eval_x = np.copy(x)
        self.latest_eval_fun = fun

        self.log_line(f"[EVALUATION {self.eval_counter}]")
        self.log_line(f"  Date/time : {now_str}")
        self.log_line(f"  x         : {np.array2string(x, precision=6, separator=', ')}")
        self.log_line(f"  Objective : {fun:.6f}")

        for name, value in zip(names, constraints):
            self.log_line(f"  {name:<9}: {value:.6f}")

        self.log_line("")

    def log_iteration(self, x):
        x = np.asarray(x, dtype=float)
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        iter_num = len(self.history["iter_num"]) + 1
        constraints = self.get_constraint_values(x)
        names = self.get_constraint_names(len(constraints)) if constraints else []

        if self.prev_iter_x is None:
            step_norm = np.nan
        else:
            step_norm = float(np.linalg.norm(x - self.prev_iter_x))

        if self.latest_eval_x is not None and np.allclose(x, self.latest_eval_x, rtol=1e-10, atol=1e-10):
            fun = self.latest_eval_fun
        else:
            fun = np.nan

        self.history["iter_num"].append(iter_num)
        self.history["iter_x"].append(x.tolist())
        self.history["iter_fun"].append(fun)
        self.history["iter_time"].append(now_str)
        self.history["iter_step_norm"].append(step_norm)
        self.history["iter_constraints"].append(constraints)

        self.prev_iter_x = np.copy(x)

        self.log_line(f"[ITERATION {iter_num}]")
        self.log_line(f"  Date/time : {now_str}")
        self.log_line(f"  x         : {np.array2string(x, precision=6, separator=', ')}")
        self.log_line(f"  Objective : {fun:.6f}" if not np.isnan(fun) else "  Objective : nan")
        self.log_line(f"  Step norm : {step_norm:.6f}" if not np.isnan(step_norm) else "  Step norm : nan")

        for name, value in zip(names, constraints):
            self.log_line(f"  {name:<9}: {value:.6f}")

        self.log_line("")

    def finalize(self, result):
        end_dt = datetime.now()
        runtime_sec = (end_dt - self.start_dt).total_seconds()

        final_constraints = self.get_constraint_values(result.x)
        names = self.get_constraint_names(len(final_constraints)) if final_constraints else []

        self.log_line("=" * 80)
        self.log_line("FINAL RESULT")
        self.log_line("=" * 80)
        self.log_line(f"Run end date/time   : {end_dt.strftime('%Y-%m-%d %H:%M:%S')}")
        self.log_line(f"Total runtime [s]   : {runtime_sec:.3f}")
        self.log_line(f"Success             : {result.success}")
        self.log_line(f"Status              : {result.status}")
        self.log_line(f"Message             : {result.message}")
        self.log_line(f"Iterations          : {getattr(result, 'nit', 'N/A')}")
        self.log_line(f"Function evaluations: {getattr(result, 'nfev', 'N/A')}")
        self.log_line(f"Optimal x           : {np.array2string(np.asarray(result.x), precision=6, separator=', ')}")
        self.log_line(f"Optimal objective   : {float(result.fun):.6f}")
        self.log_line("")

        if final_constraints:
            self.log_line("Final constraint values:")
            for name, value in zip(names, final_constraints):
                self.log_line(f"  {name:<9}= {value:.6f}")
            self.log_line("")

        self.log_line("=" * 80)
        self.log_line("END OF LOG")
        self.log_line("=" * 80)