"""Bounded source snapshots for subsequent characteristic evolution.

The compact reduced scalar and its derivatives remain analytic. Only the
integrated Jordan radius correction and transverse jets are interpolated.
No source is integrated, no regional equation is evolved, and no terminal
rectangle trace is extrapolated or substituted for a moving-wall solution.
"""
from collections.abc import Mapping
import copy
import math
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from .nsc_atlas_characteristic_initialization import Budget, source_jet
from .nsc_atlas_validation import source_jet_roundoff_tolerances
from .nsc_atlas_nonlinear_bulk import bulk_operator
from .nsc_atlas_nonlinear_source import source_amplitude_bound


MAX_QUERY_POINTS = 513
MAX_EVALUATIONS = 20000
TRANSVERSE = {"r_Ubar": "r_UV", "log_C_Ubar": "log_C_UV", "c_Ubar": "c_UV"}
REQUIRED_ROWS = {"x", "h", "h_x", "rJ", "rJ_x", "r", "C", "log_C", "c",
                 "r_x", "log_C_x", "c_x", *TRANSVERSE}


def _array(value, label, shape=None):
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.size == 0:
            raise ValueError
        if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
            raise ValueError
        result = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(result)) or (shape is not None and result.shape != shape):
            raise ValueError
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(label+" must contain finite real values of the required shape") from error
    return result


def _budget(value):
    if value is None:
        value = Budget(maximum_evaluations=MAX_EVALUATIONS, cpu_seconds=30.)
    if (not isinstance(value, Budget) or type(value.maximum_evaluations) is not int
            or type(value.evaluations) is not int or value.maximum_evaluations <= 0
            or not 0 <= value.evaluations <= value.maximum_evaluations
            or isinstance(value.cpu_seconds, bool) or not math.isfinite(value.cpu_seconds)
            or value.cpu_seconds <= 0 or not math.isfinite(value.started)
            or value.started > time.process_time()):
        raise ValueError("use a live characteristic Budget with finite positive ceilings")
    value.check()
    return value


def _mixed(row, alpha):
    return bulk_operator({key: row[key] for key in ("r", "C", "c")},
        {"r_U": row["r_Ubar"], "r_V": row["r_x"],
         "c_U": row["c_Ubar"], "c_V": row["c_x"]}, alpha=alpha, region="exterior")


class SourceSnapshot:
    """Immutable interpolation data with a separately live evaluation budget.

    ``source_rows`` is the first-jet mapping returned by ExactSource.sample
    on ``x_nodes``. Saved rows must be authenticated by their caller. Inputs
    are copied. Knots cover the full original source interval, with 4..513
    points. They need not be uniform; the compact profile is never fitted.

    h, h_x and h_xx are derivatives of one cubic Hermite polynomial. Its
    actual h_xx supplies r_xx and the reported null constraint. The source
    ODE curvature target is reported separately, not substituted to force
    that constraint to vanish. Transverse jets use their saved values and
    mixed-equation slopes at the knots; their interpolation defects and the
    numerical error of the original source integration remain distinct.

    Each snapshot has a 20,000-evaluation ceiling (three reference checks,
    knot construction points, and queried points), with at most 513 points
    per query. The default Budget permits 30 CPU seconds. A supplied shared
    Budget keeps its original clock and ceilings; snapshot evaluation does
    not consult or reset an expired ExactSource child budget.
    """

    def __init__(self, geometry, *, alpha, x_nodes, source_rows, budget=None):
        self.budget = _budget(budget)
        self.evaluations = 0
        self.query_points = 0
        self.query_calls = 0
        self.origin = "saved_first_jet_rows"
        self._geometry = copy.deepcopy(dict(geometry))
        self._alpha = float(_array(alpha, "alpha", ()))
        self._charge(3)
        self._bounds = source_amplitude_bound(self._geometry)
        if not 0 <= self._alpha < self._bounds["alpha_star"]:
            raise ValueError("snapshot alpha exits the original source-only threshold")
        nodes = _array(x_nodes, "x_nodes")
        if (nodes.ndim != 1 or not 4 <= nodes.size <= MAX_QUERY_POINTS
                or np.any(np.diff(nodes) <= 0)
                or nodes[0] != self._geometry["x_corner"] or nodes[-1] != self._geometry["x_last"]):
            raise ValueError("snapshot knots must increasingly cover the full original source interval")
        if not isinstance(source_rows, Mapping) or not REQUIRED_ROWS <= set(source_rows):
            raise ValueError("source_rows must include the original source values and first jets")
        rows = {key: _array(source_rows[key], "source_rows."+key, nodes.shape) for key in REQUIRED_ROWS}
        if not np.array_equal(rows["x"], nodes):
            raise ValueError("source rows and snapshot knots differ")
        self._charge(int(nodes.size))
        frame = self._frame(nodes, rows["h"], rows["h_x"], np.zeros_like(nodes))
        pulse_atol = dict(zip(("c", "c_x"), source_jet_roundoff_tolerances(nodes, geometry)))
        for key in ("rJ", "rJ_x", "r", "C", "log_C", "c", "r_x", "log_C_x", "c_x"):
            if not np.allclose(rows[key], frame[key], rtol=2e-13, atol=2e-13+pulse_atol.get(key, 0.)):
                raise ValueError("source frame or compact profile mismatch: "+key)
        # The common quiet corner is fixed by the original affine gauge.
        a0 = (.8/nodes[0]-1)/4
        corner = {"h": 0., "h_x": 0., "r_Ubar": a0,
                  "c_Ubar": 0., "log_C_Ubar": -(1/.8+1/nodes[0])*a0}
        if any(abs(rows[key][0]-value) > 5e-13 for key, value in corner.items()):
            raise ValueError("source snapshot is incompatible with the quiet corner")
        self._nodes = nodes
        self._nodes.flags.writeable = False
        self._h = CubicHermiteSpline(nodes, rows["h"], rows["h_x"], extrapolate=False)
        rhs = _mixed(rows, self.alpha)
        self._transverse = {key: CubicHermiteSpline(nodes, rows[key], rhs[output], extrapolate=False)
                            for key, output in TRANSVERSE.items()}
        self.budget.check()

    @classmethod
    def from_live(cls, source, x_nodes, *, budget):
        """Capture one live source sample up front; never reset its child clock.

        The supplied budget must be the source budget or its declared parent,
        so the live sample and snapshot work have cumulative accounting.
        Call this before the original source child expires. Later sample()
        calls use only the retained polynomials and analytic compact source.
        """
        budget = _budget(budget)
        child = getattr(source, "budget", None)
        if child is not budget and getattr(child, "parent", None) is not budget:
            raise ValueError("live capture must share the source budget or its parent")
        if not callable(getattr(source, "sample", None)):
            raise ValueError("live source must provide sample(x)")
        nodes = _array(x_nodes, "x_nodes")
        if nodes.ndim != 1 or not 4 <= nodes.size <= MAX_QUERY_POINTS:
            raise ValueError("live capture requires 4..513 source knots")
        rows = source.sample(nodes)
        result = cls(source.geometry, alpha=source.alpha, x_nodes=nodes,
                     source_rows=rows, budget=budget)
        result.origin = "single_upfront_live_source_sample"
        return result

    @property
    def geometry(self):
        return copy.deepcopy(self._geometry)

    @property
    def alpha(self):
        return self._alpha

    @property
    def x_nodes(self):
        return self._nodes.copy()

    @property
    def domain(self):
        CJ0 = -3.2/self._geometry["U0"]
        return {**{key: self._geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")},
                "x_a": self._nodes[0], "x_b": self._nodes[-1], "Ubar_a": 0.,
                "Ubar_b": CJ0*(self._geometry["Ua"]-self._geometry["U0"]), "CJ0": CJ0}

    def _charge(self, count):
        self.evaluations += count
        self.budget.check(count)
        if self.evaluations > MAX_EVALUATIONS:
            raise RuntimeError("source snapshot evaluation budget exhausted")

    def _frame(self, x, h, hx, hxx):
        p, px, pxx = source_jet(x, self._geometry)
        quiet = (x <= self._geometry["support_x"][0]) | (x >= self._geometry["support_x"][1])
        g, gx, gxx = (np.where(quiet, 0., value) for value in (x*p, p+x*px, 2*px+x*pxx))
        alpha = self.alpha
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            rJ, rJx, rJxx = x+alpha*h, 1+alpha*hx, alpha*hxx
            F = 1-alpha*g*g/(6*rJ*rJ)
            if (np.any(rJ <= self._bounds["r_min"]) or np.any(abs(rJx) >= self._bounds["v_max"])
                    or np.any(F <= self._bounds["F_min"])):
                raise ValueError("source snapshot exits the source continuation guards")
            r = rJ*np.sqrt(F)
            rx = (rJ*rJx-alpha*g*gx/6)/r
            rxx = (rJx*rJx+rJ*rJxx-alpha*(gx*gx+g*gxx)/6-rx*rx)/r
            pJ, pJx = g/rJ, (gx-g*rJx/rJ)/rJ
            pJxx = (gxx-2*pJx*rJx-pJ*rJxx)/rJ
            Lx = -alpha*pJ*pJx/(3*F)
            ux = (gx-g*rx/r)/r
            Q = (g*gxx-2*gx*gx)/(6*rJ)+g*gx*rJx/(3*rJ*rJ)
            result = {"x": x, "h": h, "h_x": hx, "h_xx": hxx,
                "rJ": rJ, "rJ_x": rJx, "rJ_xx": rJxx, "r": r, "r_x": rx, "r_xx": rxx,
                "C": F, "log_C": np.log1p(-alpha*g*g/(6*rJ*rJ)), "log_C_x": Lx,
                "c": g, "c_x": gx, "c_xx": gxx,
                "source_ode_h_xx": Q, "source_h_xx_defect": hxx-Q,
                "source_constraint": rxx-Lx*rx+alpha*r*F*ux*ux/2,
                "source_constraint_jordan": np.sqrt(F)*rJxx+alpha*rJ*(2*pJx*pJx-pJ*pJxx)/(6*np.sqrt(F)),
                "source_constraint_from_curvature_defect": alpha*(hxx-Q)/np.sqrt(F)}
        return result

    def sample(self, x):
        """Return consistent snapshot values/jets; no extrapolation or fallback.

        Knot second derivatives use SciPy's right interval convention, with
        the final knot using the left interval. h is C1, not generally C2.
        source_constraint is the actual interpolant residual; mixed_*_defect
        similarly compares differentiated transverse polynomials with the
        mixed field equations. Neither residual bounds source ODE error.
        """
        x = _array(x, "source snapshot query")
        if (x.ndim > 1 or x.size > MAX_QUERY_POINTS
                or np.any(x < self._nodes[0]) or np.any(x > self._nodes[-1])):
            raise ValueError("source snapshot query exits the closed original interval or batch limit")
        self._charge(int(x.size))
        self.query_points += int(x.size)
        self.query_calls += 1
        try:
            row = self._frame(x, self._h(x), self._h(x, 1), self._h(x, 2))
            row.update({key: polynomial(x) for key, polynomial in self._transverse.items()})
            rhs = _mixed(row, self.alpha)
            for key, output in TRANSVERSE.items():
                row[key+"_x"] = self._transverse[key](x, 1)
                row["mixed_"+key+"_defect"] = row[key+"_x"]-rhs[output]
            row["m"] = rhs["m"]
        except (FloatingPointError, OverflowError) as error:
            raise ValueError("nonfinite source snapshot arithmetic") from error
        if not all(np.all(np.isfinite(value)) for value in row.values()):
            raise ValueError("nonfinite source snapshot output")
        self.budget.check()
        return row

    def summary(self):
        return {"schema": "nsc-atlas-source-snapshot-v1", "origin": self.origin,
                "alpha": self.alpha, "domain": self.domain, "knot_count": len(self._nodes),
                "evaluations": self.evaluations, "maximum_evaluations": MAX_EVALUATIONS,
                "query_points": self.query_points, "query_calls": self.query_calls,
                "profile": "original analytic compact g,g_x,g_xx",
                "interpolation": "cubic Hermite h and transverse source jets",
                "source_numerical_error_is_enclosed": False,
                "interpolation_error_is_enclosed": False,
                "source_integrations": 0, "regional_marches": 0}
