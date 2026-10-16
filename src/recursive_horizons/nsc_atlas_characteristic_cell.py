"""One bounded implicit-midpoint cell of the normalized mixed equations.

Corners use X=(r,L,c), L=log(C). SE is displaced in the first null
coordinate and NW in the second. This constructs no grid or wall data;
an algebraic root is not a continuum or null-constraint error bound.
"""
from collections.abc import Mapping

import numpy as np

from .nsc_atlas_characteristic_initialization import Budget
from .nsc_atlas_nonlinear_bulk import bulk_operator


NAMES = ("r", "L", "c")


class CellSolveError(RuntimeError):
    """A bounded local solve stopped; diagnostics preserve its last state."""

    def __init__(self, reason, diagnostics):
        self.reason = reason
        self.diagnostics = dict(diagnostics, reason=reason)
        super().__init__("characteristic cell: " + reason)


def _number(value, name):
    raw = np.asarray(value)
    if raw.shape != () or raw.dtype.kind not in "iuf":
        raise ValueError(name + " must be a finite real scalar")
    value = float(raw)
    if not np.isfinite(value):
        raise ValueError(name + " must be a finite real scalar")
    return value


def _state(value, name):
    if not isinstance(value, Mapping) or set(value) != set(NAMES):
        raise ValueError(name + " must contain exactly r, L, c")
    return np.array([_number(value[key], name + "." + key) for key in NAMES])


def _triple(value, name, *, positive):
    if isinstance(value, (list, tuple)) and any(isinstance(x, (bool, np.bool_)) for x in value):
        raise ValueError(name + " must contain three finite real values")
    raw = np.asarray(value)
    if raw.shape != (3,) or raw.dtype.kind not in "iuf":
        raise ValueError(name + " must contain three finite real values")
    result = raw.astype(float)
    if not np.all(np.isfinite(result)) or np.any(result <= 0 if positive else result < 0):
        raise ValueError(name + " has invalid values")
    return result


def _mapping(values):
    return dict(zip(NAMES, map(float, values)))


class _Cell:
    def __init__(self, sw, nw, se, du, dv, alpha, region, tube, budget):
        self.sw, self.nw, self.se = (_state(v, n) for v, n in ((sw, "SW"), (nw, "NW"), (se, "SE")))
        self.du, self.dv = _number(du, "du"), _number(dv, "dv")
        self.alpha = _number(alpha, "alpha")
        if self.du <= 0 or self.dv <= 0 or self.alpha < 0:
            raise ValueError("require du,dv>0 and alpha>=0")
        self.area = self.du*self.dv
        if not np.isfinite(self.area) or self.area <= 0:
            raise ValueError("cell area must be positive and representable")
        if not isinstance(region, str) or region not in ("cap", "exterior"):
            raise ValueError("region must be cap or exterior")
        self.region = region
        if not isinstance(tube, Mapping) or set(tube) != set(NAMES):
            raise ValueError("tube must contain finite bounds for r, L, c")
        limits = []
        for key in NAMES:
            pair = tube[key]
            if not isinstance(pair, (list, tuple, np.ndarray)) or np.shape(pair) != (2,):
                raise ValueError("tube bounds must be pairs")
            low, high = (_number(x, "tube." + key) for x in pair)
            if not low < high:
                raise ValueError("tube bounds must be strictly increasing")
            limits.append((low, high))
        self.limits = np.array(limits)
        if self.limits[0, 0] <= 0:
            raise ValueError("tube radius floor must be positive")
        with np.errstate(over="ignore", under="ignore"):
            lapses = np.exp(self.limits[1])
        if not np.all(np.isfinite(lapses)) or np.any(lapses <= 0):
            raise ValueError("tube lapse bounds must be positive and representable")
        for value in (self.sw, self.nw, self.se):
            if not self.inside(value):
                raise ValueError("known corner exits the declared tube")
        self.base = self.nw+self.se-self.sw
        if not np.all(np.isfinite(self.base)):
            raise ValueError("corner arithmetic is not representable")
        self.budget = budget if budget is not None else Budget(cpu_seconds=30.)
        if not callable(getattr(self.budget, "check", None)):
            raise ValueError("budget must expose check(count)")
        self.rhs_evaluations = 0
        self.jacobian_evaluations = 0
        try:
            self.budget.check()
        except RuntimeError as error:
            raise CellSolveError("budget_exhausted", dict(self.counts(), detail=str(error))) from error

    def inside(self, value):
        return (np.all(np.isfinite(value)) and np.all(value >= self.limits[:, 0])
                and np.all(value <= self.limits[:, 1]))

    def counts(self):
        return {"rhs_evaluations": self.rhs_evaluations,
                "jacobian_evaluations": self.jacobian_evaluations}

    def evaluate(self, ne, *, jacobian=False):
        if not self.inside(ne):
            raise ValueError("NE exits the declared tube")
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            mid = (ne+self.nw+self.se+self.sw)/4
            xu = (ne+self.se-self.nw-self.sw)/(2*self.du)
            xv = (ne+self.nw-self.se-self.sw)/(2*self.dv)
        if not self.inside(mid) or not np.all(np.isfinite([*xu, *xv])):
            raise ValueError("midpoint or first jets exit the finite cell domain")
        r, L, c = mid
        C = float(np.exp(L))
        jets = dict(r_U=float(xu[0]), r_V=float(xv[0]), c_U=float(xu[2]), c_V=float(xv[2]))
        self.rhs_evaluations += 1
        self.budget.check(1)
        row = bulk_operator(dict(r=r, C=C, c=c), jets, alpha=self.alpha, region=self.region)
        rhs = np.array([row["r_UV"], row["log_C_UV"], row["c_UV"]], dtype=float)
        residual = ne-self.base-self.area*rhs
        if not np.all(np.isfinite(residual)):
            raise ValueError("nonfinite algebraic residual")
        result = {"ne": _mapping(ne), "midpoint": _mapping(mid), "midpoint_C": C,
                  "first_jets": dict(jets, log_C_U=float(xu[1]), log_C_V=float(xv[1])),
                  "mixed_rhs": _mapping(rhs), "residual": residual}
        if jacobian:
            self.jacobian_evaluations += 1
            self.budget.check(1)
            result["jacobian"] = np.eye(3)-self.area*self._rhs_jacobian(mid, xu, xv, row)
        return result

    def _rhs_jacobian(self, mid, xu, xv, row):
        """Forward directional chain rule with each NE component as a seed."""
        r, L, c = mid
        C = np.exp(L)
        rU, rV, cU, cV = xu[0], xv[0], xu[2], xv[2]
        alpha, h2 = self.alpha, int(self.region == "cap")
        u = c/r
        uU, uV = (cU-u*rU)/r, (cV-u*rV)/r
        q, D = alpha*u*u/6, row["D"]
        increment = row["potential_increment"]
        product = rU*rV
        vac = product+C*(1-h2*r*r)/4
        derivatives = np.empty((3, 3))
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            for j in range(3):
                seed = np.eye(3)[j]
                dr, dL, dc = seed/4
                dC = C*dL
                drU, drV = seed[0]/(2*self.du), seed[0]/(2*self.dv)
                dcU, dcV = seed[2]/(2*self.du), seed[2]/(2*self.dv)
                du = (dc-u*dr)/r
                duU = (dcU-du*rU-u*drU-uU*dr)/r
                duV = (dcV-du*rV-u*drV-uV*dr)/r
                dD = alpha*u*du/3
                d_increment = 6*h2*(1+q)*dD
                dsigma = alpha*((duU*uV+uU*duV)/D-uU*uV*dD/(D*D))
                dproduct = drU*rV+rU*drV
                dsum = dproduct+dC/4
                drhs = -dsum/r+(product+C/4)*dr/(r*r)
                drhs += (dC*r+C*dr)*(3*h2+increment)/4+C*r*d_increment/4
                dLrhs = 2*dsum/(r*r)-4*(product+C/4)*dr/(r**3)-dsigma
                dvac = dproduct+dC*(1-h2*r*r)/4-C*h2*r*dr/2
                dK = dvac/(r*r)-2*vac*dr/(r**3)
                dK -= (dC*increment+C*d_increment)/12+dsigma/6
                derivatives[:, j] = (drhs, dLrhs, -c*dK-row["K"]*dc)
        if not np.all(np.isfinite(derivatives)):
            raise ValueError("nonfinite analytic Jacobian")
        return derivatives


def cell_evaluation(ne, sw, nw, se, *, du, dv, alpha, region, tube, budget=None,
                    jacobian=False):
    """Evaluate the local equation; one RHS and optionally one Jacobian charge.

    No solve or grid is started. The returned mixed first jets do not supply
    the pure second derivatives needed to evaluate null constraints.
    """
    cell = _Cell(sw, nw, se, du, dv, alpha, region, tube, budget)
    result = cell.evaluate(_state(ne, "NE"), jacobian=jacobian)
    return dict(result, **cell.counts())


def solve_characteristic_cell(sw, nw, se, *, du, dv, alpha, region, tube,
                              atol=(1e-12, 1e-12, 1e-12), rtol=1e-12,
                              reference_scales=(1., 1., 1.), max_updates=12,
                              max_backtracks=8, maximum_condition=1e12,
                              budget=None, initial_guess=None):
    """Solve exactly one cell, or raise CellSolveError with bounded diagnostics.

    Fixed scales are atol+rtol*max(reference_scales,abs(known corners)).
    Each accepted update uses a full analytic Newton direction with at most
    max_backtracks trial lengths 1,1/2,... . A tube never changes in a solve.
    The default affine predictor is used only if it is admissible; otherwise
    a recorded convex average of NW and SE supplies a positive initial guess.
    Invalid input raises ValueError; numerical failure is explicit and never
    triggers cell subdivision, a changed amplitude or a different equation.
    """
    cell = _Cell(sw, nw, se, du, dv, alpha, region, tube, budget)
    atol = _triple(atol, "atol", positive=True)
    reference = _triple(reference_scales, "reference_scales", positive=True)
    rtol = _number(rtol, "rtol")
    maximum_condition = _number(maximum_condition, "maximum_condition")
    if rtol < 0 or maximum_condition < 1:
        raise ValueError("require rtol>=0 and maximum_condition>=1")
    for value, label, ceiling in ((max_updates, "max_updates", 12), (max_backtracks, "max_backtracks", 8)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not 1 <= value <= ceiling:
            raise ValueError(label + " exceeds the admitted positive iteration range")
    scales = atol+rtol*np.maximum(reference, np.max(np.abs([cell.sw, cell.nw, cell.se]), axis=0))
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("fixed error scales must be positive and finite")
    if initial_guess is None:
        ne = cell.base.copy()
        seed_kind = "affine_predictor"
        if not cell.inside(ne):
            ne = (cell.nw+cell.se)/2
            seed_kind = "convex_corner_seed"
    else:
        ne = _state(initial_guess, "initial_guess")
        if not cell.inside(ne):
            raise ValueError("initial_guess exits the declared tube")
        seed_kind = "supplied"
    diagnostics = dict(updates=0, backtrack_trials=0, step_factors=[], seed_kind=seed_kind,
                       scales=scales.tolist(), ne=_mapping(ne), residual_norm=None,
                       correction_norm=None, condition_number=None)

    def fail(reason):
        raise CellSolveError(reason, dict(diagnostics, **cell.counts()))

    try:
        for iteration in range(max_updates+1):
            row = cell.evaluate(ne, jacobian=True)
            residual = row["residual"]/scales
            norm = float(np.linalg.norm(residual, ord=np.inf))
            J = row["jacobian"]*scales[None, :]/scales[:, None]
            condition = float(np.linalg.cond(J))
            diagnostics.update(ne=_mapping(ne), residual_norm=norm, condition_number=condition)
            if not np.isfinite(condition) or condition > maximum_condition:
                fail("ill_conditioned_jacobian")
            direction_scaled = np.linalg.solve(J, -residual)
            cell.budget.check()
            correction = float(np.linalg.norm(direction_scaled, ord=np.inf))
            diagnostics["correction_norm"] = correction
            if norm <= 1 and correction <= 1:
                result = {key: value for key, value in row.items() if key != "jacobian"}
                result["residual"] = row["residual"].tolist()
                return dict(result, **diagnostics, **cell.counts())
            if iteration == max_updates:
                fail("newton_update_limit")
            accepted = False
            for trial in range(max_backtracks):
                factor = 2.**(-trial)
                candidate = ne+factor*scales*direction_scaled
                diagnostics["backtrack_trials"] += 1
                if not cell.inside(candidate):
                    cell.budget.check()
                    continue
                try:
                    trial_row = cell.evaluate(candidate)
                    trial_norm = float(np.linalg.norm(trial_row["residual"]/scales, ord=np.inf))
                except (ValueError, FloatingPointError):
                    continue
                if trial_norm <= (1-1e-4*factor)*norm or (trial_norm <= 1 and factor*correction <= 1):
                    ne = candidate
                    diagnostics["updates"] += 1
                    diagnostics["step_factors"].append(factor)
                    accepted = True
                    break
            if not accepted:
                fail("backtracking_limit")
    except CellSolveError:
        raise
    except np.linalg.LinAlgError:
        fail("linear_solve_failure")
    except (ValueError, FloatingPointError, OverflowError) as error:
        diagnostics["detail"] = str(error)
        fail("nonfinite_or_invalid_evaluation")
    except RuntimeError as error:
        diagnostics["detail"] = str(error)
        fail("budget_exhausted")
