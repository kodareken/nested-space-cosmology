"""Finite forced scalar evolution on the original clipped exterior grid.

The state is z=y-g1(V)-B(U), so its moving-wall derivative vanishes. Known
forcing coefficients may be cached in a moving eta coordinate; the unknown
and its Volterra integral remain in the original U,V chart. No physical
preparation, metric query or direct boundary quadrature occurs in this solver.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_field_preparation import exterior_coefficient, radius_from_null


SCHEMA = "nsc-atlas-scalar-evolution-v1"
SCOPE = "formal exterior cubic coefficient; no central waveform or finite-amplitude remainder"
DOMAIN_KEYS = ("U0", "Ua", "Ub", "Va", "Vb")
MAX_GRID_NODES = 257
MAX_QUERY = 513
MAX_RHS_CALLS = 6000


def _array(value, name, maximum=MAX_QUERY):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError(name+" must contain finite real numbers")
    result = np.asarray(raw, dtype=float)
    if result.ndim > 1 or not 1 <= result.size <= maximum or not np.all(np.isfinite(result)):
        raise ValueError(name+" must be a bounded finite scalar or one-dimensional array")
    return result


def _number(value, name):
    result = _array(value, name)
    if result.ndim != 0:
        raise ValueError(name+" must be scalar")
    return float(result)


def _time(value):
    value = _number(value, "normalized time")
    if not 0 <= value <= 1:
        raise ValueError("scalar evolution query exits the finite time interval")
    return value


def reference_kernel(U, V):
    """The retained leading exterior coefficient, evaluated only inside the wall."""
    return exterior_coefficient(radius_from_null(U, V))


class ScalarEvolutionSolution:
    """Dense lifted solution with source and material wall readout.

    Construct with solve_scalar_correction. The forcing object supplies
    eval(U,t), domain, geometry['wall'], training_times, and the shared metric
    budget; eval returns U,F,G,g1,g1_V,Fsource. The kernel can be injected for
    analytic controls. All kernel values share the preparation budget.
    """

    def __init__(self, forcing, lift, U_grid, *, tolerance, atol, max_step,
                 rhs_calls_max, kernel):
        self.forcing, self.lift = forcing, lift
        self.budget, self.geometry = forcing.budget, forcing.geometry
        self.domain = {key: _number(forcing.domain[key], key) for key in DOMAIN_KEYS}
        d = self.domain
        if not (d["U0"] < d["Ua"] < d["Ub"] and d["Va"] < d["Vb"]):
            raise ValueError("invalid finite scalar evolution domain")
        if (lift.budget is not self.budget or any(not math.isclose(
                _number(lift.domain[key], "lift "+key), d[key], rel_tol=0, abs_tol=1e-13)
                for key in DOMAIN_KEYS)):
            raise ValueError("forcing and lift must share their finite domain and cumulative budget")
        self.external = getattr(self.budget, "external", None)
        if (self.external is None or self.external.started != self.budget.started
                or not callable(getattr(self.external, "check", None))):
            raise ValueError("scalar evolution needs the original shared preparation budget")
        self.grid = np.array(_array(U_grid, "correction grid", MAX_GRID_NODES), copy=True)
        if (self.grid.ndim != 1 or self.grid.size < 2 or np.any(np.diff(self.grid) <= 0)
                or self.grid[0] != d["U0"] or self.grid[-1] != d["Ub"]):
            raise ValueError("correction grid must increase strictly from U0 to Ub")
        self.grid.setflags(write=False)
        self.dV = d["Vb"]-d["Va"]
        self.kernel = reference_kernel if kernel is None else kernel
        if not callable(self.kernel):
            raise ValueError("kernel must be callable")
        tolerance, atol, max_step = (_number(value, label) for value, label in
            ((tolerance, "relative tolerance"), (atol, "absolute tolerance"), (max_step, "maximum step")))
        if not (0 < tolerance <= 1e-5 and 0 < atol <= 1e-7 and 0 < max_step <= .125):
            raise ValueError("invalid bounded scalar ODE controls")
        if (isinstance(rhs_calls_max, (bool, np.bool_)) or not isinstance(rhs_calls_max, (int, np.integer))
                or not 1 <= rhs_calls_max <= MAX_RHS_CALLS):
            raise ValueError("scalar RHS ceiling must be between 1 and 6000")
        self.rhs_calls_max = int(rhs_calls_max)
        self.settings = dict(grid_nodes=int(self.grid.size), tolerance=tolerance, atol=atol,
                             max_step=max_step, rhs_calls_max=self.rhs_calls_max)
        self.rhs_calls = self.readout_rhs_calls = self.kernel_values = 0
        self.budget.check()
        endpoints = [self._wall(t)[1] for t in (0., 1.)]
        endpoint_tolerance = 8*np.finfo(float).eps*max(1., abs(d["Ua"]), abs(d["Ub"]))
        if not np.allclose(endpoints, [d["Ua"], d["Ub"]], rtol=0, atol=endpoint_tolerance):
            raise ValueError("wall endpoints do not cover the correction grid")
        self.B_grid = self._lift_values(self.grid)
        if np.any(self.B_grid[self.grid <= d["Ua"]] != 0):
            raise ValueError("lift must preserve the exact quiet bottom anchor")
        times = [0., .25, .75, 1.]
        for label, raw in (("forcing time knots", getattr(forcing, "training_times", (0., 1.))),
                           ("lift time knots", getattr(lift, "knot_times", (0., 1.)))):
            values = _array(raw, label, 257)
            if values.ndim != 1 or np.any(values < 0) or np.any(values > 1):
                raise ValueError(label+" exit the original finite interval")
            times.extend(values.tolist())
        births = []
        for U in self.grid[(self.grid > d["Ua"]) & (self.grid < d["Ub"])]:
            root = brentq(lambda t: self._wall(float(t))[1]-U, 0., 1., xtol=1e-14)
            births.append(float(root))
        # Avoid zero-length DOP853 segments when two independently obtained
        # breakpoints differ only by floating-point representation.
        times.extend(births)
        ordered = sorted(set(times))
        merged = [ordered[0]]
        merge_width = 0.
        for t in ordered[1:]:
            if t-merged[-1] <= 8*np.finfo(float).eps:
                merge_width = max(merge_width, t-merged[-1])
                if t in (0., .25, .75, 1.):
                    merged[-1] = t
            else:
                merged.append(t)
        self.cuts = np.asarray(merged)
        self.wall_birth_times = births
        self.maximum_split_merge_width = merge_width
        state = np.zeros(self.grid.size)
        initial = self._inputs(0.)
        if abs(initial["g1"]) > 1e-24:
            raise ValueError("source correction must have the quiet bottom anchor")
        self.segments = []
        for a, b in zip(self.cuts[:-1], self.cuts[1:]):
            solved = solve_ivp(self._rate, (a, b), state, method="DOP853", rtol=tolerance,
                atol=atol, first_step=min(max_step, b-a), max_step=min(max_step, b-a), dense_output=True)
            if not solved.success or solved.t[-1] != b or not np.all(np.isfinite(solved.y[:, -1])):
                raise RuntimeError("scalar correction failed to reach its finite endpoint")
            state = solved.y[:, -1]
            self.segments.append(solved.sol)
        self.dense_data = {"grid": self.grid, "cuts": self.cuts, "segments": tuple(self.segments)}
        self.budget.check()

    def _wall(self, t):
        V = min(self.domain["Vb"], max(self.domain["Va"], self.domain["Va"]+self.dV*t))
        raw = _number(self.geometry["wall"](V), "wall coordinate")
        tolerance = 8*np.finfo(float).eps*max(1., abs(self.domain["Ua"]), abs(self.domain["Ub"]))
        if not self.domain["Ua"]-tolerance <= raw <= self.domain["Ub"]+tolerance:
            raise ValueError("wall exits the retained correction domain")
        wall = min(self.domain["Ub"], max(self.domain["Ua"], raw))
        derivative = _number(self.geometry["wall"](V, 1), "wall derivative")
        if derivative <= 0:
            raise ValueError("scalar correction needs an increasing wall map")
        return V, wall, derivative

    def _lift_values(self, U):
        result = np.asarray(self.lift.B(U), dtype=float)
        if result.shape != np.asarray(U).shape or not np.all(np.isfinite(result)):
            raise ValueError("invalid boundary lift values")
        return result

    def _inputs(self, t):
        self.budget.check()
        V, wall, wall_V = self._wall(t)
        clipped = np.minimum(self.grid, wall)
        row = self.forcing.eval(clipped, t)
        result = {"V": V, "wall": wall, "wall_V": wall_V, "U": clipped}
        for key in ("F", "G", "U"):
            array = np.asarray(row[key], dtype=float)
            if array.shape != clipped.shape or not np.all(np.isfinite(array)):
                raise ValueError("cached scalar forcing has invalid "+key)
            if key == "U":
                if not np.allclose(array, clipped, rtol=0, atol=8*np.finfo(float).eps):
                    raise ValueError("cached forcing changed the original clipped coordinates")
            else:
                result[key] = array
        for key in ("g1", "g1_V", "Fsource"):
            result[key] = _number(row[key], "cached "+key)
        if result["F"][0] != result["Fsource"]:
            raise ValueError("source subtraction must use the same cached F(U0,V)")
        result["B_wall"] = float(self._lift_values(np.asarray(wall)))
        return result

    def _rate(self, t, state, *, counted=True, return_inputs=False):
        t = _time(t)
        if counted:
            self.rhs_calls += 1
            if self.rhs_calls > self.rhs_calls_max:
                raise RuntimeError("scalar correction RHS budget exhausted")
        else:
            self.readout_rhs_calls += 1
        data = self._inputs(t)
        if state.shape != self.grid.shape or not np.all(np.isfinite(state)):
            raise ValueError("invalid lifted scalar state")
        active = self.grid < data["wall"]
        z = np.where(active, state, state[-1])
        B = np.where(active, self.B_grid, data["B_wall"])
        y = z+data["g1"]+B
        self.external.check(int(self.grid.size))
        self.kernel_values += self.grid.size
        K = np.asarray(self.kernel(data["U"], data["V"]), dtype=float)
        if K.shape != self.grid.shape or not np.all(np.isfinite(K)):
            raise ValueError("invalid retained exterior kernel")
        integrand = K*y+data["G"]
        prefix = np.r_[0., np.cumsum(.5*(integrand[:-1]+integrand[1:])*np.diff(data["U"]))]
        rate_V = data["F"]-data["Fsource"]-prefix
        self.budget.check()
        if return_inputs:
            return self.dV*rate_V, {**data, "z": z, "y": y, "rate_V": rate_V,
                "source_state_defect": float(state[0]),
                "inactive_state_spread": float(np.max(abs(state[~active]-state[-1])))}
        return self.dV*rate_V

    def _state(self, t):
        self.budget.check()
        index = min(int(np.searchsorted(self.cuts, t, side="right")-1), len(self.segments)-1)
        return np.asarray(self.segments[index](t), dtype=float)

    def sample(self, t):
        """One wall row; y_wall_V is the total derivative along the wall."""
        t = _time(t)
        _, row = self._rate(t, self._state(t), counted=False, return_inputs=True)
        b = _number(self.lift.b_at_U(row["wall"]), "boundary derivative")
        return {"t": t, "V": row["V"], "U_wall": row["wall"], "U_wall_V": row["wall_V"],
            "g1": row["g1"], "g1_V": row["g1_V"], "B_wall": row["B_wall"], "b_wall": b,
            "z_wall": float(row["z"][-1]), "y_wall": float(row["y"][-1]),
            "z_wall_V": float(row["rate_V"][-1]),
            "y_wall_V": float(row["g1_V"]+row["rate_V"][-1]+row["wall_V"]*b),
            "source_defect": row["source_state_defect"], "source_rate_defect": float(row["rate_V"][0]),
            "wall_flux": float(row["F"][-1]),
            "inactive_state_spread": row["inactive_state_spread"]}

    def wall(self, times):
        times = _array(times, "wall times")
        rows = [self.sample(_time(t)) for t in times.ravel()]
        return {key: np.asarray([row[key] for row in rows]).reshape(times.shape) for key in rows[0]}

    def bulk(self, U, t):
        """Piecewise-linear retained readout of z and y; no derivative claim."""
        t = _time(t)
        points = _array(U, "bulk U query")
        V, wall, _ = self._wall(t)
        if np.any(points < self.domain["U0"]) or np.any(points > wall):
            raise ValueError("scalar bulk readout exits the retained exterior")
        state = self._state(t)
        clipped = np.minimum(self.grid, wall)
        z_grid = np.where(self.grid < wall, state, state[-1])
        z = np.interp(points, clipped, z_grid)
        g1 = _number(self.forcing.eval(np.array([self.domain["U0"]]), t)["g1"], "cached source")
        return {"U": points.copy(), "V": np.full_like(points, V), "z": z,
                "y": z+g1+self._lift_values(points)}

    def summary(self, times=np.linspace(0., 1., 9)):
        rows = [self.sample(_time(t)) for t in _array(times, "summary times").ravel()]
        return {"schema": SCHEMA, "scope": SCOPE, "domain": dict(self.domain), "settings": dict(self.settings),
            "wall_samples": rows, "rhs_calls": self.rhs_calls, "readout_rhs_calls": self.readout_rhs_calls,
            "kernel_values": int(self.kernel_values), "segments": len(self.segments),
            "split_times": self.cuts.tolist(), "maximum_split_merge_width": self.maximum_split_merge_width,
            "maximum_absolute_source_defect": max(abs(row["source_defect"]) for row in rows),
            "maximum_absolute_source_rate_defect": max(abs(row["source_rate_defect"]) for row in rows),
            "maximum_absolute_wall_flux": max(abs(row["wall_flux"]) for row in rows)}


def solve_scalar_correction(forcing, lift, U_grid=None, *, tolerance=2e-10, atol=1e-12,
                            max_step=1/64, rhs_calls_max=MAX_RHS_CALLS, kernel=None):
    """Evolve the formal cubic exterior coefficient using already checked inputs."""
    if U_grid is None:
        U_grid = np.linspace(forcing.domain["U0"], forcing.domain["Ub"], 129)
    return ScalarEvolutionSolution(forcing, lift, U_grid, tolerance=tolerance, atol=atol,
        max_step=max_step, rhs_calls_max=rhs_calls_max, kernel=kernel)
