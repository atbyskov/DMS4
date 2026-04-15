from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.optimize import linprog

ArrayLike = Union[Sequence[float], np.ndarray]
ObjectiveFn = Callable[[np.ndarray], Union[float, Tuple[float, np.ndarray]]]
ConstraintFn = Callable[
    [np.ndarray],
    Union[
        Tuple[ArrayLike, ArrayLike],
        Tuple[ArrayLike, ArrayLike, np.ndarray, np.ndarray],
    ],
]


@dataclass
class SLPOptions:
    """Options for the Sequential Linear Programming solver.

    The names intentionally mirror the original MATLAB version closely.
    """

    Algorithm: str = "merit"  # 'merit' or 'al'
    Solver: str = "scipy"  # only scipy/highs is implemented here
    Display: str = "off"  # 'off' or 'iter'

    # Convergence parameters
    MaxFunctionEvaluations: int = 1000
    MaxIterations: int = 1000
    InfeasibilityPenalization: float = 1000.0
    OptimalityTolerance: float = 1e-6
    FunctionTolerance: float = 1e-6
    StepTolerance: float = 1e-8
    ObjectiveLimit: float = -1e20

    # Gradient parameters
    FiniteDifferenceType: str = "forward"  # 'forward', 'backward', 'central'
    FiniteDifferenceStepSize: float = np.sqrt(np.finfo(float).eps)
    SpecifyConstraintGradient: bool = False
    SpecifyObjectiveGradient: bool = False
    CheckGradients: bool = False

    # Move-limit parameters
    MoveLimitMethod: str = "adaptive"
    MoveLimit: float = 0.1
    MoveLimitExpand: float = 1.1
    MoveLimitReduce: float = 0.5

    # Global convergence filter parameters
    MaxInfeasibility: float = np.inf

    def __post_init__(self) -> None:
        self.Algorithm = self.Algorithm.lower()
        self.Solver = self.Solver.lower()
        self.Display = self.Display.lower()
        self.FiniteDifferenceType = self.FiniteDifferenceType.lower()
        self.MoveLimitMethod = self.MoveLimitMethod.lower()

        if self.Algorithm not in {"merit", "al"}:
            raise ValueError("Algorithm must be 'merit' or 'al'.")
        if self.Solver not in {"scipy", "linprog"}:
            raise ValueError("Only Solver='scipy'/'linprog' is supported in this Python version.")
        if self.Display not in {"off", "iter"}:
            raise ValueError("Display must be 'off' or 'iter'.")
        if self.FiniteDifferenceType not in {"forward", "backward", "central"}:
            raise ValueError("FiniteDifferenceType must be 'forward', 'backward', or 'central'.")
        if self.MoveLimitMethod not in {"adaptive"}:
            raise ValueError("Only MoveLimitMethod='adaptive' is supported.")

        positive_fields = [
            "MaxFunctionEvaluations",
            "MaxIterations",
            "InfeasibilityPenalization",
            "OptimalityTolerance",
            "FunctionTolerance",
            "StepTolerance",
            "FiniteDifferenceStepSize",
            "MoveLimit",
            "MoveLimitExpand",
            "MoveLimitReduce",
        ]
        for name in positive_fields:
            value = getattr(self, name)
            if not (np.isfinite(value) and value > 0):
                raise ValueError(f"{name} must be positive and finite.")


@dataclass
class GlobalConvergenceFilter:
    SmallVal: float = 1.0e-6
    gamma: float = 1.0e-6
    beta: float = 1.0 - 1.0e-6
    sigma: float = 2.0e-6
    delta: float = 1.0e-6
    vals: np.ndarray = field(default_factory=lambda: np.zeros((1, 2), dtype=float))
    nVals: int = 1
    PointAcceptedByFilter: bool = False
    h: float = 1.0e30
    f: float = 1.0e30
    initF: float = 0.0


@dataclass
class IterationHistory:
    f: List[float] = field(default_factory=list)
    xnorm: List[float] = field(default_factory=list)
    constrviolation: List[float] = field(default_factory=list)
    maxInf: List[float] = field(default_factory=list)


@dataclass
class SLPOutput:
    iterations: int = 0
    funcCount: int = 0
    constrviolation: float = 0.0
    firstorderopt: float = 0.0
    message: str = ""
    iterHistory: IterationHistory = field(default_factory=IterationHistory)
    lambda_: Optional[np.ndarray] = None
    trueObjectiveCalls: int = 0
    trueConstraintCalls: int = 0
    acceptedSteps: int = 0
    rejectedSteps: int = 0


class FminSLP:
    """Sequential Linear Programming optimizer.

    This is a faithful Python rewrite of the supplied MATLAB class, with a few
    defensive fixes and a SciPy LP backend.

    Expected callback signatures
    ----------------------------
    fun(x)
        If SpecifyObjectiveGradient=False:
            return f
        If SpecifyObjectiveGradient=True:
            return f, grad

    nonlcon(x)
        If SpecifyConstraintGradient=False:
            return g, geq
        If SpecifyConstraintGradient=True:
            return g, geq, dg, dgeq

        where:
            g    : inequality constraints, satisfied when g <= 0
            geq  : equality constraints, satisfied when geq == 0
            dg   : shape (n_dv, n_ineq)
            dgeq : shape (n_dv, n_eq)
    """

    name = "fminslp"
    version = "v1.5-python"

    def __init__(
        self,
        fun: ObjectiveFn,
        x0: ArrayLike,
        A: Optional[ArrayLike] = None,
        b: Optional[ArrayLike] = None,
        Aeq: Optional[ArrayLike] = None,
        beq: Optional[ArrayLike] = None,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
        nonlcon: Optional[ConstraintFn] = None,
        **options: Any,
    ) -> None:
        if not callable(fun):
            raise TypeError("fun must be callable.")
        self.fun = fun

        self.x0 = self._col_vec_required(x0, "x0")
        self.nDV = self.x0.size

        self.A = self._matrix_or_empty(A)
        self.b = self._col_vec_or_empty(b)
        self.Aeq = self._matrix_or_empty(Aeq)
        self.beq = self._col_vec_or_empty(beq)
        self.lb = self._col_vec_with_default(lb, -np.inf, self.nDV, "lb")
        self.ub = self._col_vec_with_default(ub, np.inf, self.nDV, "ub")
        self.nonlcon = nonlcon

        if np.any(self.lb > self.ub):
            raise ValueError("Each lower bound must be <= its upper bound.")
        if np.any(self.x0 < self.lb) or np.any(self.x0 > self.ub):
            raise ValueError("x0 must lie within [lb, ub].")

        if self.A.size > 0:
            if self.b.size == 0:
                raise ValueError("b must be provided when A is provided.")
            if self.A.shape[0] != self.b.size:
                raise ValueError("A and b must have the same number of rows.")
            if self.A.shape[1] != self.nDV:
                raise ValueError("A must have nDV columns.")
        elif self.b.size > 0:
            raise ValueError("A must be provided when b is provided.")

        if self.Aeq.size > 0:
            if self.beq.size == 0:
                raise ValueError("beq must be provided when Aeq is provided.")
            if self.Aeq.shape[0] != self.beq.size:
                raise ValueError("Aeq and beq must have the same number of rows.")
            if self.Aeq.shape[1] != self.nDV:
                raise ValueError("Aeq must have nDV columns.")
        elif self.beq.size > 0:
            raise ValueError("Aeq must be provided when beq is provided.")

        if nonlcon is not None and not callable(nonlcon):
            raise TypeError("nonlcon must be callable or None.")

        self.options = SLPOptions(**options)
        self.filter = self.initialize_global_convergence_filter(self.options)

        self.initialized = True
        self.f0: float = 0.0
        self.nGnl: int = 0
        self.aFac: float = 1.0

        self._true_objective_calls = 0
        self._true_constraint_calls = 0

        if self.options.CheckGradients:
            self.check_user_supplied_gradients()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def solve(self) -> Tuple[np.ndarray, float, int, SLPOutput]:
        exitflag = -1
        output = SLPOutput()

        if self.options.Display == "iter":
            print("*" * 121)
            print("\t\t\t\t\t fminslp optimizer with global convergence filter")
            print("*" * 121)
            print(f"\t {'f(x)':>10s} \t\t {'Max inf':>10s} \t\t {'Norm dx':>10s} \t\t {'nFeval':>10s} \t\t {'IterNo':>10s}")

        min_dv_box_limit = 2.0 * self.options.StepTolerance

        x = self.x0.copy()
        xLcur = self.lb.copy()
        xUcur = self.ub.copy()
        xold1 = x.copy()
        xold2 = x.copy()

        y, ylb, yub = self._initialize_slack_variables(x)

        if self.Aeq.size > 0:
            Ameq = np.zeros((self.Aeq.shape[0], self.Aeq.shape[1] + self.nGnl), dtype=float)
            Ameq[:, : self.Aeq.shape[1]] = self.Aeq
        else:
            Ameq = np.empty((0, self.nDV + self.nGnl), dtype=float)

        if self.A.size > 0:
            Am = np.zeros((self.A.shape[0] + self.nGnl, self.A.shape[1] + self.nGnl), dtype=float)
            Am[: self.A.shape[0], : self.A.shape[1]] = self.A
            bm = np.zeros(self.b.size + self.nGnl, dtype=float)
            bm[: self.b.size] = self.b
        else:
            Am = np.empty((0, self.nDV + self.nGnl), dtype=float)
            bm = np.empty(0, dtype=float)

        if self.options.Algorithm == "al":
            lambda_vec = np.ones(self.nGnl, dtype=float)
            lambda_vec = self.update_lambda(y, lambda_vec)
        else:
            lambda_vec = np.empty(0, dtype=float)

        self.f0 = float(self._eval_objective_value(x))
        self.aFac = max(abs(self.f0), 1.0)

        fmerit, freal, _ = self.get_merit_obj(x, y, lambda_vec)
        self.filter.initF = float(fmerit)
        fOld = float(fmerit)

        maxInf = float(max(np.max(y), 0.0)) if y.size else 0.0
        optimalityNorm = 0.0

        nFeval = 1
        iterNo = 0
        optimize = True

        while optimize:
            iterNo += 1

            _, _, dfmerit = self.get_merit_obj(x, y, lambda_vec)
            gmerit_cur, _, dgmerit = self.get_merit_constraints(x, y)

            if self.nonlcon is not None and self.A.size > 0:
                if self.nGnl > 0:
                    Am[self.A.shape[0] :, :] = dgmerit
                    bm[self.b.size :] = dgmerit @ np.concatenate([x, y]) - gmerit_cur
            elif self.nonlcon is not None and self.A.size == 0:
                Am = dgmerit.copy()
                bm = dgmerit @ np.concatenate([x, y]) - gmerit_cur
            else:
                Am = np.empty((0, self.nDV + self.nGnl), dtype=float)
                bm = np.empty(0, dtype=float)

            reduceSwitch = False
            xLcur, xUcur = self.adaptive_move_limit(
                x,
                xLcur,
                xUcur,
                self.lb,
                self.ub,
                self.options.MoveLimit,
                self.options.MoveLimitReduce,
                self.options.MoveLimitExpand,
                xold1,
                xold2,
                reduceSwitch,
                min_dv_box_limit,
            )

            xold2 = xold1.copy()
            xold1 = x.copy()

            backtrack = True
            AddToFilter = False
            xNew = x.copy()
            yNew = y.copy()

            while backtrack:
                zlb = np.concatenate([xLcur, ylb])
                zub = np.concatenate([xUcur, yub])
                z0 = np.concatenate([x, y])

                z_new, lp_exitflag, lp_message = self.lp_solver(
                    dfmerit,
                    Am,
                    bm,
                    Ameq,
                    self.beq,
                    zlb,
                    zub,
                    z0,
                )

                if lp_exitflag == 1:
                    xNew = z_new[: self.nDV]
                    yNew = z_new[self.nDV :]

                    deltaxy = z_new - z0
                    deltax = xNew - x
                    deltanorm = float(np.linalg.norm(deltax))
                    optimalityNorm = float(np.linalg.norm(dfmerit[: self.nDV]))

                    gmerit_new, greal_new, _ = self.get_merit_constraints(xNew, yNew)
                    nFeval += 1
                    fmerit, freal, _ = self.get_merit_obj(xNew, yNew, lambda_vec)
                    deltaf = float(fOld - fmerit)

                    self.filter.h = float(max(np.max(gmerit_new), 0.0)) if gmerit_new.size else 0.0
                    denom = self.filter.initF if abs(self.filter.initF) > 1e-30 else 1.0
                    self.filter.f = float(fmerit / denom)
                    self.filter = self.evaluate_current_design_point_to_filter(self.filter)

                    AddToFilter = False
                    if self.filter.PointAcceptedByFilter:
                        deltaL = float(-dfmerit @ deltaxy)
                        if (deltaf < self.filter.sigma * deltaL) and (deltaL > 0.0):
                            reduceSwitch = True
                            output.rejectedSteps += 1
                        else:
                            reduceSwitch = False
                            backtrack = False
                            if self.filter.h > 0.0:
                                AddToFilter = True
                            output.acceptedSteps += 1
                    else:
                        reduceSwitch = True
                        output.rejectedSteps += 1

                    if reduceSwitch:
                        xLcur, xUcur = self.adaptive_move_limit(
                            x,
                            xLcur,
                            xUcur,
                            self.lb,
                            self.ub,
                            self.options.MoveLimit,
                            self.options.MoveLimitReduce,
                            self.options.MoveLimitExpand,
                            xold1,
                            xold2,
                            reduceSwitch,
                            min_dv_box_limit,
                        )

                    if optimalityNorm <= self.options.OptimalityTolerance:
                        optimize = False
                        backtrack = False
                        exitflag = 1
                        output.message = (
                            f"Successfully solved to Optimality Tolerance: <= {self.options.OptimalityTolerance:0.5e}"
                        )
                    elif abs(deltaf) <= self.options.FunctionTolerance:
                        optimize = False
                        backtrack = False
                        exitflag = 2
                        output.message = (
                            f"Successfully solved to Function Tolerance: <= {self.options.FunctionTolerance:0.5e}"
                        )
                    elif deltanorm <= self.options.StepTolerance:
                        optimize = False
                        backtrack = False
                        exitflag = 3
                        output.message = (
                            f"Successfully solved to Step Tolerance: <= {self.options.StepTolerance:0.5e}"
                        )
                    elif (fmerit <= self.options.ObjectiveLimit) and (iterNo > 1):
                        optimize = False
                        backtrack = False
                        exitflag = 4
                        output.message = (
                            f"Successfully solved to Objective Limit: <= {self.options.ObjectiveLimit:0.5e}"
                        )
                    elif iterNo >= self.options.MaxIterations:
                        optimize = False
                        backtrack = False
                        exitflag = 0
                        output.message = (
                            f"Number of iterations exceeded the limit: {self.options.MaxIterations:d}"
                        )
                    elif nFeval >= self.options.MaxFunctionEvaluations:
                        optimize = False
                        backtrack = False
                        exitflag = 0
                        output.message = (
                            f"Number of function evaluations exceeded the limit: {self.options.MaxFunctionEvaluations:d}"
                        )
                else:
                    optimize = False
                    backtrack = False
                    output.message = f"LP solver failed: {lp_message}"

            if optimize:
                if AddToFilter:
                    self.filter = self.update_filter(self.filter, self.filter.h, self.filter.f)

                x = xNew.copy()
                y = yNew.copy()

                if self.options.Algorithm == "al":
                    lambda_vec = self.update_lambda(y, lambda_vec)

                fOld = float(fmerit)
                maxInf = float(max(np.max(greal_new), 0.0)) if greal_new.size else 0.0

                output.iterHistory.f.append(float(freal))
                output.iterHistory.xnorm.append(float(deltanorm))
                output.iterHistory.constrviolation.append(float(maxInf))
                if self.nGnl > 0:
                    output.iterHistory.maxInf.append(float(maxInf))

                if self.options.Display == "iter":
                    print(
                        f"\t {freal:6.4e} \t\t {maxInf:6.4e} \t\t {deltanorm:6.4e}"
                        f" \t\t {nFeval:10d} \t\t {iterNo:10d}"
                    )

        fval = float(self._eval_objective_value(x))

        output.constrviolation = float(maxInf)
        output.iterations = int(iterNo)
        output.funcCount = int(nFeval)
        output.firstorderopt = float(optimalityNorm)
        output.trueObjectiveCalls = int(self._true_objective_calls)
        output.trueConstraintCalls = int(self._true_constraint_calls)

        if self.options.Algorithm == "al":
            output.lambda_ = lambda_vec / self.aFac

        return x, fval, exitflag, output

    # ------------------------------------------------------------------
    # Merit / constraints
    # ------------------------------------------------------------------
    def get_merit_obj(
        self,
        x: np.ndarray,
        y: np.ndarray,
        lambda_vec: np.ndarray,
    ) -> Tuple[float, float, np.ndarray]:
        if self.options.SpecifyObjectiveGradient:
            fval, df = self._eval_objective_and_gradient(x)
        else:
            fval = self._eval_objective_value(x)
            df = self.get_fun_dsa(x, fin=fval)

        if self.options.Algorithm == "merit":
            dy = self.aFac * (self.options.InfeasibilityPenalization + y)
            fmerit = fval + self.aFac * np.sum(
                y * self.options.InfeasibilityPenalization + 0.5 * y**2
            )
        elif self.options.Algorithm == "al":
            dy = self.aFac * (lambda_vec + self.options.InfeasibilityPenalization * y)
            fmerit = fval + self.aFac * np.sum(
                y * lambda_vec + self.options.InfeasibilityPenalization * 0.5 * y**2
            )
        else:
            raise ValueError(f"Unknown algorithm '{self.options.Algorithm}'.")

        dfmerit = np.concatenate([df, dy])
        return float(fmerit), float(fval), dfmerit

    def get_merit_constraints(
        self,
        x: np.ndarray,
        y: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.nonlcon is None:
            return np.empty(0, dtype=float), np.empty(0, dtype=float), np.empty((0, self.nDV + self.nGnl), dtype=float)

        if self.options.SpecifyConstraintGradient:
            gn, gneq, dgnl, dgneq = self._eval_constraints_and_gradients(x)
        else:
            gn, gneq = self._eval_constraints(x)
            dgnl, dgneq = self.get_nonlcon_dsa(x, gnl0=gn, gnleq0=gneq)

        gn = self._as_1d(gn)
        gneq = self._as_1d(gneq)
        greal = np.concatenate([gn, gneq, -gneq])
        gmerit = greal - y

        dgmerit = np.zeros((self.nGnl, self.nDV + self.nGnl), dtype=float)
        if self.nGnl > 0:
            parts = []
            if gn.size > 0:
                parts.append(self._normalize_constraint_gradient(dgnl, self.nDV, gn.size).T)
            if gneq.size > 0:
                dgeq_norm = self._normalize_constraint_gradient(dgneq, self.nDV, gneq.size).T
                parts.append(dgeq_norm)
                parts.append(-dgeq_norm)
            if parts:
                dgmerit[:, : self.nDV] = np.vstack(parts)
            dgmerit[:, self.nDV :] = -np.eye(self.nGnl)

        return gmerit, greal, dgmerit

    def update_lambda(self, y: np.ndarray, lambda_vec: np.ndarray) -> np.ndarray:
        return np.maximum(
            lambda_vec + self.options.InfeasibilityPenalization * y,
            1.0,
        )

    # ------------------------------------------------------------------
    # LP subproblem
    # ------------------------------------------------------------------
    def lp_solver(
        self,
        c: np.ndarray,
        A_ub: np.ndarray,
        b_ub: np.ndarray,
        A_eq: np.ndarray,
        b_eq: np.ndarray,
        lb: np.ndarray,
        ub: np.ndarray,
        x0: np.ndarray,
    ) -> Tuple[np.ndarray, int, str]:
        del x0  # Reserved for future LP solvers that accept warm starts.

        nvar = c.size
        bounds = [(float(lb[i]), float(ub[i])) for i in range(nvar)]

        A_ub_in = None if A_ub.size == 0 else A_ub
        b_ub_in = None if b_ub.size == 0 else b_ub
        A_eq_in = None if A_eq.size == 0 else A_eq
        b_eq_in = None if b_eq.size == 0 else b_eq

        result = linprog(
            c=c,
            A_ub=A_ub_in,
            b_ub=b_ub_in,
            A_eq=A_eq_in,
            b_eq=b_eq_in,
            bounds=bounds,
            method="highs",
        )

        if result.success:
            return np.asarray(result.x, dtype=float), 1, result.message
        return np.zeros_like(c), 0, result.message

    # ------------------------------------------------------------------
    # Finite differences
    # ------------------------------------------------------------------
    def get_fun_dsa(self, x: np.ndarray, fin: Optional[float] = None) -> np.ndarray:
        df = np.zeros(self.nDV, dtype=float)
        h = self.options.FiniteDifferenceStepSize
        fd_type = self.options.FiniteDifferenceType

        if fd_type == "forward":
            if fin is None:
                fin = self._eval_objective_value(x)
            for i in range(self.nDV):
                xp = x.copy()
                xp[i] += h
                p1 = self._eval_objective_value(xp)
                df[i] = (p1 - fin) / h
        elif fd_type == "backward":
            if fin is None:
                fin = self._eval_objective_value(x)
            for i in range(self.nDV):
                xm = x.copy()
                xm[i] -= h
                m1 = self._eval_objective_value(xm)
                df[i] = (fin - m1) / h
        elif fd_type == "central":
            for i in range(self.nDV):
                xp = x.copy()
                xm = x.copy()
                xp[i] += h
                xm[i] -= h
                p1 = self._eval_objective_value(xp)
                m1 = self._eval_objective_value(xm)
                df[i] = (p1 - m1) / (2.0 * h)
        else:
            raise ValueError("Unknown FiniteDifferenceType.")

        return df

    def get_nonlcon_dsa(
        self,
        x: np.ndarray,
        gnl0: Optional[np.ndarray] = None,
        gnleq0: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if self.nonlcon is None:
            return np.empty((self.nDV, 0), dtype=float), np.empty((self.nDV, 0), dtype=float)

        h = self.options.FiniteDifferenceStepSize
        fd_type = self.options.FiniteDifferenceType

        if gnl0 is None or gnleq0 is None:
            gnl0, gnleq0 = self._eval_constraints(x)
        gnl0 = self._as_1d(gnl0)
        gnleq0 = self._as_1d(gnleq0)

        ngnl = gnl0.size
        ngnleq = gnleq0.size
        dg = np.zeros((self.nDV, ngnl), dtype=float)
        dgeq = np.zeros((self.nDV, ngnleq), dtype=float)

        if fd_type == "forward":
            for i in range(self.nDV):
                xp = x.copy()
                xp[i] += h
                gp1, geqp1 = self._eval_constraints(xp)
                if ngnl:
                    dg[i, :] = (gp1 - gnl0) / h
                if ngnleq:
                    dgeq[i, :] = (geqp1 - gnleq0) / h
        elif fd_type == "backward":
            for i in range(self.nDV):
                xm = x.copy()
                xm[i] -= h
                gm1, geqm1 = self._eval_constraints(xm)
                if ngnl:
                    dg[i, :] = (gnl0 - gm1) / h
                if ngnleq:
                    dgeq[i, :] = (gnleq0 - geqm1) / h
        elif fd_type == "central":
            for i in range(self.nDV):
                xp = x.copy()
                xm = x.copy()
                xp[i] += h
                xm[i] -= h
                gp1, geqp1 = self._eval_constraints(xp)
                gm1, geqm1 = self._eval_constraints(xm)
                if ngnl:
                    dg[i, :] = (gp1 - gm1) / (2.0 * h)
                if ngnleq:
                    dgeq[i, :] = (geqp1 - geqm1) / (2.0 * h)
        else:
            raise ValueError("Unknown FiniteDifferenceType.")

        return dg, dgeq

    def check_user_supplied_gradients(self) -> None:
        original_fd_type = self.options.FiniteDifferenceType
        self.options.FiniteDifferenceType = "central"
        try:
            if self.options.SpecifyObjectiveGradient:
                _, df_user = self._eval_objective_and_gradient(self.x0)
                df_fd = self.get_fun_dsa(self.x0)
                max_diff = float(np.max(np.abs(df_user - df_fd))) if df_user.size else 0.0
                print("\n\tDerivative Check Information")
                print(
                    f" Objective function derivatives:\n"
                    f" Maximum difference between user-supplied and finite-difference derivatives = {max_diff:0.5e}"
                )

            if self.options.SpecifyConstraintGradient and self.nonlcon is not None:
                _, _, dgnl_user, dgneq_user = self._eval_constraints_and_gradients(self.x0)
                gn0, geq0 = self._eval_constraints(self.x0)
                dgnl_fd, dgneq_fd = self.get_nonlcon_dsa(self.x0, gn0, geq0)

                if dgnl_user.size:
                    dgnl_user_n = self._normalize_constraint_gradient(dgnl_user, self.nDV, self._as_1d(gn0).size)
                    dsa_diff = np.abs(dgnl_user_n - dgnl_fd)
                    idx = int(np.argmax(dsa_diff))
                    dv_no, g_no = np.unravel_index(idx, dsa_diff.shape)
                    max_diff = float(dsa_diff[dv_no, g_no])
                    print("\n\tDerivative Check Information")
                    print(
                        f" Nonlinear inequality constraint derivatives:\n"
                        f" Maximum difference between user-supplied and finite-difference derivatives = {max_diff:0.5e}\n"
                        f"\t User-supplied constraint derivative element ({dv_no},{g_no}): {dgnl_user_n[dv_no, g_no]:0.5e}\n"
                        f"\t Finite-difference constraint derivative element ({dv_no},{g_no}): {dgnl_fd[dv_no, g_no]:0.5e}"
                    )

                if dgneq_user.size:
                    dgneq_user_n = self._normalize_constraint_gradient(dgneq_user, self.nDV, self._as_1d(geq0).size)
                    dsa_diff = np.abs(dgneq_user_n - dgneq_fd)
                    idx = int(np.argmax(dsa_diff))
                    dv_no, g_no = np.unravel_index(idx, dsa_diff.shape)
                    max_diff = float(dsa_diff[dv_no, g_no])
                    print("\n\tDerivative Check Information")
                    print(
                        f" Nonlinear equality constraint derivatives:\n"
                        f" Maximum difference between user-supplied and finite-difference derivatives = {max_diff:0.5e}\n"
                        f"\t User-supplied constraint derivative element ({dv_no},{g_no}): {dgneq_user_n[dv_no, g_no]:0.5e}\n"
                        f"\t Finite-difference constraint derivative element ({dv_no},{g_no}): {dgneq_fd[dv_no, g_no]:0.5e}"
                    )
        finally:
            self.options.FiniteDifferenceType = original_fd_type

    # ------------------------------------------------------------------
    # Global convergence filter
    # ------------------------------------------------------------------
    @staticmethod
    def initialize_global_convergence_filter(options: SLPOptions) -> GlobalConvergenceFilter:
        small = 1.0e-6
        filt = GlobalConvergenceFilter(
            SmallVal=small,
            gamma=small,
            beta=1.0 - small,
            sigma=2.0 * small,
            delta=small,
            vals=np.zeros((options.MaxIterations + 1, 2), dtype=float),
            nVals=1,
            PointAcceptedByFilter=False,
            h=1.0e30,
            f=1.0e30,
            initF=0.0,
        )

        if np.isfinite(options.MaxInfeasibility):
            filt.vals[0, 0] = options.MaxInfeasibility
            filt.vals[0, 1] = -np.inf
        else:
            filt.vals[0, 0] = np.inf
            filt.vals[0, 1] = np.inf
        return filt

    @staticmethod
    def evaluate_current_design_point_to_filter(filter_obj: GlobalConvergenceFilter) -> GlobalConvergenceFilter:
        h = filter_obj.h
        f = filter_obj.f
        accepted = True
        for ii in range(filter_obj.nVals):
            hi = filter_obj.vals[ii, 0]
            fi = filter_obj.vals[ii, 1]
            if (h <= hi * filter_obj.beta) or ((f + filter_obj.gamma * h) <= fi):
                continue
            accepted = False
            break
        filter_obj.PointAcceptedByFilter = accepted
        return filter_obj

    @staticmethod
    def update_filter(
        filter_obj: GlobalConvergenceFilter,
        hk: float,
        fk: float,
    ) -> GlobalConvergenceFilter:
        new_filter = GlobalConvergenceFilter(
            SmallVal=filter_obj.SmallVal,
            gamma=filter_obj.gamma,
            beta=filter_obj.beta,
            sigma=filter_obj.sigma,
            delta=filter_obj.delta,
            vals=np.zeros_like(filter_obj.vals),
            nVals=0,
            PointAcceptedByFilter=False,
            h=filter_obj.h,
            f=filter_obj.f,
            initF=filter_obj.initF,
        )

        for ii in range(filter_obj.nVals):
            hi = filter_obj.vals[ii, 0]
            fi = filter_obj.vals[ii, 1]
            dominated = (hk <= hi) and (fk <= fi)
            if not dominated:
                new_filter.vals[new_filter.nVals, 0] = hi
                new_filter.vals[new_filter.nVals, 1] = fi
                new_filter.nVals += 1

        new_filter.vals[new_filter.nVals, 0] = hk
        new_filter.vals[new_filter.nVals, 1] = fk
        new_filter.nVals += 1
        return new_filter

    # ------------------------------------------------------------------
    # Move limits
    # ------------------------------------------------------------------
    @staticmethod
    def adaptive_move_limit(
        x: np.ndarray,
        xLcur: np.ndarray,
        xUcur: np.ndarray,
        xLorg: np.ndarray,
        xUorg: np.ndarray,
        move_limit: float,
        reduce_fac: float,
        expand_fac: float,
        xold1: np.ndarray,
        xold2: np.ndarray,
        reduce_switch: bool,
        min_dv_box_limit: float = 0.0,
    ) -> Tuple[np.ndarray, np.ndarray]:
        xLcur = xLcur.copy()
        xUcur = xUcur.copy()

        if reduce_switch:
            expand = reduce_fac
            reduction = reduce_fac
        else:
            reduction = reduce_fac
            expand = expand_fac

        n_dv = x.size
        for i in range(n_dv):
            delta = 0.5 * (xUcur[i] - xLcur[i])

            if abs(x[i] - xold1[i]) > 1.0e-10:
                s1 = (xold1[i] - xold2[i]) / (x[i] - xold1[i])
                if s1 < 0.0:
                    delta *= reduction
                else:
                    delta *= expand
            else:
                delta *= move_limit

            dmax = (xUorg[i] - xLorg[i]) * move_limit
            if np.isfinite(dmax):
                delta = min(delta, dmax)

            xLcur[i] = max(x[i] - delta, xLorg[i])
            xUcur[i] = min(x[i] + delta, xUorg[i])

            if xLcur[i] > xUcur[i]:
                if min_dv_box_limit > 0:
                    mid = 0.5 * (xLcur[i] + xUcur[i])
                    xLcur[i] = mid - 0.5 * min_dv_box_limit
                    xUcur[i] = mid + 0.5 * min_dv_box_limit
                else:
                    xLcur[i] = (1.0 - 1.0e-6) * xUcur[i]

            width = xUcur[i] - xLcur[i]
            if width < min_dv_box_limit:
                target = min_dv_box_limit - width
                grow_left = min(0.5 * target, max(0.0, xLcur[i] - xLorg[i]))
                grow_right = min(0.5 * target, max(0.0, xUorg[i] - xUcur[i]))

                xLcur[i] -= grow_left
                xUcur[i] += grow_right

                remaining = min_dv_box_limit - (xUcur[i] - xLcur[i])
                if remaining > 0:
                    extra_right = min(remaining, max(0.0, xUorg[i] - xUcur[i]))
                    xUcur[i] += extra_right
                    remaining = min_dv_box_limit - (xUcur[i] - xLcur[i])
                if remaining > 0:
                    extra_left = min(remaining, max(0.0, xLcur[i] - xLorg[i]))
                    xLcur[i] -= extra_left

                xLcur[i] = max(xLcur[i], xLorg[i])
                xUcur[i] = min(xUcur[i], xUorg[i])

        return xLcur, xUcur

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
    def _initialize_slack_variables(
        self, x: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.nonlcon is None:
            self.nGnl = 0
            return np.empty(0, dtype=float), np.empty(0, dtype=float), np.empty(0, dtype=float)

        gnl, gnleq = self._eval_constraints(x)
        gnl = self._as_1d(gnl)
        gnleq = self._as_1d(gnleq)

        self.nGnl = gnl.size + 2 * gnleq.size
        y_seed = np.concatenate([gnl, gnleq, -gnleq]) if self.nGnl > 0 else np.empty(0, dtype=float)
        y = np.maximum(y_seed, 0.0)
        ylb = np.zeros(self.nGnl, dtype=float)
        if y_seed.size > 0:
            upper_seed = max(float(np.max(y_seed)), 1.0)
        else:
            upper_seed = 1.0
        yub = np.full(self.nGnl, upper_seed * 1.0e6, dtype=float)
        return y, ylb, yub

    def _eval_objective_value(self, x: np.ndarray) -> float:
        self._true_objective_calls += 1
        out = self.fun(np.asarray(x, dtype=float))
        if isinstance(out, tuple):
            return float(out[0])
        return float(out)

    def _eval_objective_and_gradient(self, x: np.ndarray) -> Tuple[float, np.ndarray]:
        self._true_objective_calls += 1
        out = self.fun(np.asarray(x, dtype=float))
        if not isinstance(out, tuple) or len(out) != 2:
            raise ValueError(
                "When SpecifyObjectiveGradient=True, fun(x) must return (f, grad)."
            )
        f, grad = out
        grad_arr = self._as_1d(grad)
        if grad_arr.size != self.nDV:
            raise ValueError("Objective gradient must have length nDV.")
        return float(f), grad_arr

    def _eval_constraints(self, x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if self.nonlcon is None:
            return np.empty(0, dtype=float), np.empty(0, dtype=float)
        self._true_constraint_calls += 1
        out = self.nonlcon(np.asarray(x, dtype=float))
        if not isinstance(out, tuple) or len(out) < 2:
            raise ValueError("nonlcon(x) must return at least (g, geq).")
        g, geq = out[:2]
        return self._as_1d(g), self._as_1d(geq)

    def _eval_constraints_and_gradients(
        self, x: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if self.nonlcon is None:
            return (
                np.empty(0, dtype=float),
                np.empty(0, dtype=float),
                np.empty((self.nDV, 0), dtype=float),
                np.empty((self.nDV, 0), dtype=float),
            )
        self._true_constraint_calls += 1
        out = self.nonlcon(np.asarray(x, dtype=float))
        if not isinstance(out, tuple) or len(out) != 4:
            raise ValueError(
                "When SpecifyConstraintGradient=True, nonlcon(x) must return (g, geq, dg, dgeq)."
            )
        g, geq, dg, dgeq = out
        g_arr = self._as_1d(g)
        geq_arr = self._as_1d(geq)
        dg_arr = self._normalize_constraint_gradient(dg, self.nDV, g_arr.size)
        dgeq_arr = self._normalize_constraint_gradient(dgeq, self.nDV, geq_arr.size)
        return g_arr, geq_arr, dg_arr, dgeq_arr

    @staticmethod
    def _normalize_constraint_gradient(
        grad: Any,
        n_dv: int,
        n_con: int,
    ) -> np.ndarray:
        if n_con == 0:
            return np.empty((n_dv, 0), dtype=float)
        arr = np.asarray(grad, dtype=float)
        if arr.shape == (n_dv, n_con):
            return arr
        if arr.shape == (n_con, n_dv):
            return arr.T
        raise ValueError(
            f"Constraint gradient must have shape ({n_dv}, {n_con}) or ({n_con}, {n_dv})."
        )

    @staticmethod
    def _as_1d(x: Any) -> np.ndarray:
        if x is None:
            return np.empty(0, dtype=float)
        arr = np.asarray(x, dtype=float)
        if arr.size == 0:
            return np.empty(0, dtype=float)
        return np.ravel(arr).astype(float)

    @staticmethod
    def _matrix_or_empty(x: Optional[ArrayLike]) -> np.ndarray:
        if x is None:
            return np.empty((0, 0), dtype=float)
        arr = np.asarray(x, dtype=float)
        if arr.size == 0:
            return np.empty((0, 0), dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.ndim != 2:
            raise ValueError("Expected a 2D array for a matrix input.")
        return arr

    @staticmethod
    def _col_vec_or_empty(x: Optional[ArrayLike]) -> np.ndarray:
        if x is None:
            return np.empty(0, dtype=float)
        arr = np.asarray(x, dtype=float)
        if arr.size == 0:
            return np.empty(0, dtype=float)
        return np.ravel(arr).astype(float)

    @staticmethod
    def _col_vec_required(x: ArrayLike, name: str) -> np.ndarray:
        arr = np.asarray(x, dtype=float)
        if arr.size == 0:
            raise ValueError(f"{name} must be non-empty.")
        return np.ravel(arr).astype(float)

    @staticmethod
    def _col_vec_with_default(
        x: Optional[ArrayLike],
        default_value: float,
        size: int,
        name: str,
    ) -> np.ndarray:
        if x is None:
            return np.full(size, default_value, dtype=float)
        arr = np.ravel(np.asarray(x, dtype=float)).astype(float)
        if arr.size != size:
            raise ValueError(f"{name} must have length {size}.")
        return arr


# ----------------------------------------------------------------------
# Minimal example
# ----------------------------------------------------------------------
if __name__ == "__main__":
    def fun(x: np.ndarray) -> float:
        return (x[0] - 1.0) ** 2 + (x[1] - 2.0) ** 2

    def nonlcon(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        g = np.array([
            x[0] ** 2 + x[1] ** 2 - 4.0,   # inside circle
            -x[0] - x[1] + 1.0,            # x0 + x1 >= 1
        ])
        geq = np.array([])
        return g, geq

    solver = FminSLP(
        fun=fun,
        x0=np.array([0.5, 0.5]),
        lb=np.array([-5.0, -5.0]),
        ub=np.array([5.0, 5.0]),
        nonlcon=nonlcon,
        Display="iter",
        MaxIterations=100,
        MaxFunctionEvaluations=500,
        MoveLimit=0.25,
        FiniteDifferenceType="central",
    )

    x_opt, f_opt, exitflag, out = solver.solve()
    print("\nDone")
    print("x*       =", x_opt)
    print("f*       =", f_opt)
    print("exitflag =", exitflag)
    print("message  =", out.message)
