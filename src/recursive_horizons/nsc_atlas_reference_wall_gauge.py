"""Fixed reference wall labels in the original affine source coordinates.

This is a bounded quadrature/inversion of the known vacuum wall, not a source
or spacetime solve. W and eta prescribe a coordinate gauge; returned vacuum
r,L,T are controls, not replacement fields for finite-amplitude evolution.
"""
from collections import OrderedDict
from collections.abc import Mapping
import copy
import math
import time

import numpy as np
from scipy.optimize import brentq

from .nsc_atlas_characteristic_initialization import Budget, validate_geometry


MAX_QUERY_POINTS = 513
MAX_CACHE_POINTS = 4096
INVERSE_MAX_ITERATIONS = 32
R_MIN, R_MAX = .45, .55


def _array(value, label):
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.size == 0:
            raise ValueError
        if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
            raise ValueError
        result = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(result)):
            raise ValueError
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(label+" must contain finite real values") from error
    return result


def _F(radius):
    return (1-radius/.8)*np.exp(radius/.8)


def _state(radius):
    """Analytic signed wall kinematics, regular through beta_plus=0."""
    R = np.asarray(radius)
    bm, bp = 1/(5*R**2)+3*R/4, 1/(5*R**2)-5*R/4
    q = 25*R**2/16-1+3/(10*R)+1/(25*R**4)
    v = np.sqrt(q)
    v_R = (25*R/8-3/(10*R**2)-4/(25*R**5))/(2*v)
    bp_R = -2/(5*R**3)-5/4
    difference = v-bp
    advanced_R = -1/(v*difference)
    log_advanced_R = -v_R/v-(v_R-bp_R)/difference
    return bm, bp, v, advanced_R, log_advanced_R


class ReferenceWallGauge:
    """Evaluate the reference W(x), its analytic jets and relative clocks.

    The supplied geometry is the unchanged source_geometry(64) dictionary.
    I(R)=integral_.45^R[-1/(v*(v-beta_plus))]dR is inverted against
    1.6*log(V(x)/Va), V(x)=F(x)/U0, on the closed .45..55 wall segment.
    Fixed Gauss order 32 or 64 controls this smooth integral; derivatives
    use analytic wall equations, never differences of sampled coordinates.

    s=T=0 at x=.4. T is the reference static cap time, equal to the central
    proper clock there. An old absolute clock requires its inherited offset.
    The original Va,Vb and U endpoints are not fitted or rescaled. Endpoint
    quadrature mismatch is retained as an inversion residual.

    The default cumulative Budget is 30 CPU seconds/20,000 evaluations.
    Explicit shared budgets retain their own clock and ceilings; no extra
    lifetime cap is imposed. Every query and reference quadrature node is
    charged. Exact repeated x queries reuse at most 4,096 cached rows.
    """

    def __init__(self, geometry, *, order=64, budget=None):
        if type(order) is not int or order not in (32, 64):
            raise ValueError("reference quadrature order must be 32 or 64")
        if not isinstance(geometry, Mapping):
            raise ValueError("supply the original source geometry")
        self._geometry = copy.deepcopy(dict(geometry))
        self.order = order
        self.budget = Budget(maximum_evaluations=20000, cpu_seconds=30.) if budget is None else budget
        if (not isinstance(self.budget, Budget) or type(self.budget.maximum_evaluations) is not int
                or type(self.budget.evaluations) is not int or self.budget.maximum_evaluations <= 0
                or not 0 <= self.budget.evaluations <= self.budget.maximum_evaluations
                or isinstance(self.budget.cpu_seconds, bool) or not math.isfinite(self.budget.cpu_seconds)
                or self.budget.cpu_seconds <= 0 or not math.isfinite(self.budget.started)
                or self.budget.started > time.process_time()):
            raise ValueError("reference gauge requires a live cumulative characteristic Budget")
        self.evaluations = 0
        self.quadrature_nodes = 0
        self.query_points = 0
        self.cache_hits = 0
        self.inverse_calls = 0
        self.maximum_inverse_iterations = 0
        self.maximum_inverse_residual = 0.
        self._cache = OrderedDict()
        self._charge(3)
        validate_geometry(self._geometry)
        if type(self._geometry["order"]) is not int or self._geometry["order"] != 64:
            raise ValueError("reference gauge requires the original order-64 source geometry")
        self.CJ0 = -3.2/self._geometry["U0"]
        self._lower, self._upper = .4, float(self._geometry["x_last"])
        self._nodes, self._weights = np.polynomial.legendre.leggauss(order)
        if not math.isclose(self._geometry["U0"]*self._geometry["Va"],
                            float(_F(.4)), rel_tol=0, abs_tol=3e-14):
            raise ValueError("source corner does not have the original affine radius")
        for R, U, V in ((R_MIN, self._geometry["Ua"], self._geometry["Va"]),
                        (R_MAX, self._geometry["Ub"], self._geometry["Vb"])):
            if not math.isclose(U*V, float(_F(R)), rel_tol=0, abs_tol=3e-14):
                raise ValueError("geometry endpoints do not match the original reference wall")
        self._integral_end = self._integrals(R_MAX)
        end_target = self._target(self._upper)
        self.endpoint_integral_defect = self._integral_end[0]-end_target
        if abs(self.endpoint_integral_defect) > 3e-13:
            raise ValueError("reference integral is incompatible with original incoming endpoint")
        self.budget.check()

    @property
    def geometry(self):
        return copy.deepcopy(self._geometry)

    @property
    def domain(self):
        return {"x_a": self._lower, "x_b": self._upper, "r_a": R_MIN, "r_b": R_MAX,
                "W_a": self.CJ0*(self._geometry["Ua"]-self._geometry["U0"]),
                "W_b": self.CJ0*(self._geometry["Ub"]-self._geometry["U0"])}

    def _charge(self, count=0):
        self.evaluations += count
        self.budget.check(count)

    def _integrals(self, R):
        if R == R_MIN:
            return np.zeros(3)
        self._charge(self.order)
        self.quadrature_nodes += self.order
        radius = (R+R_MIN)/2+(R-R_MIN)*self._nodes/2
        bm, _, speed, advanced_R, _ = _state(radius)
        return (R-R_MIN)/2*np.array([
            self._weights @ advanced_R, self._weights @ (1/speed),
            self._weights @ (bm/((1-radius*radius)*speed))])

    def _target(self, x):
        if x == self._lower:
            return 0.
        V = float(_F(x))/self._geometry["U0"]
        return 1.6*math.log1p((V-self._geometry["Va"])/self._geometry["Va"])

    def _invert(self, x):
        target = self._target(x)
        if x == self._lower:
            return R_MIN, np.zeros(3), 0, 0.
        if x == self._upper:
            return R_MAX, self._integral_end.copy(), 0, self.endpoint_integral_defect
        endpoints = {R_MIN: np.zeros(3), R_MAX: self._integral_end}
        values = dict(endpoints)
        def residual(R):
            if R not in values:
                values[R] = self._integrals(R)
            return float(values[R][0]-target)
        if target < self._integral_end[0]:
            # Only inherited endpoint roundoff can put an in-domain label
            # beyond the smooth integral's bracket; never extend the domain.
            if abs(target-self._integral_end[0]) <= 3e-13 and self._upper-x <= 8*np.spacing(self._upper):
                return R_MAX, self._integral_end.copy(), 0, self._integral_end[0]-target
            raise ValueError("incoming label lies outside the reference inversion bracket")
        self.inverse_calls += 1
        R, result = brentq(residual, R_MIN, R_MAX, xtol=5e-15,
                          rtol=4*np.finfo(float).eps, maxiter=INVERSE_MAX_ITERATIONS,
                          full_output=True, disp=False)
        if not result.converged:
            raise RuntimeError("reference wall inversion did not converge within 32 iterations")
        integrals = values.get(R)
        if integrals is None:
            integrals = self._integrals(R)
        error = float(integrals[0]-target)
        if abs(error) > 3e-13:
            raise RuntimeError("reference wall inverse residual exceeds its fixed tolerance")
        return R, integrals, int(result.iterations), error

    def _evaluate(self, x):
        R, integrals, iterations, residual = self._invert(x)
        bm, bp, v, advanced_R, log_advanced_R = map(float, _state(R))
        lam = -x/(.8*(.8-x))
        lam_x = -1/(.8-x)**2
        R_x = 1.6*lam/advanced_R
        R_xx = R_x*(lam_x/lam-log_advanced_R*R_x)
        pR, pR_R = -R/(.8*(.8-R)), -1/(.8-R)**2
        Z = pR*R_x-lam
        Z_x = pR_R*R_x*R_x+pR*R_xx-lam_x
        V = float(_F(x))/self._geometry["U0"]
        U = float(_F(R))/V
        if x == self._lower:
            U = self._geometry["Ua"]
        elif x == self._upper:
            U = self._geometry["Ub"]
        W = self.CJ0*(U-self._geometry["U0"])
        W_x, W_xx = self.CJ0*U*Z, self.CJ0*U*(Z*Z+Z_x)
        C = x/R*math.exp((x-R)/.8)
        L = math.log(x/R)+(x-R)/.8
        L_x_total = 1/x+1/.8-(1/R+1/.8)*R_x
        r_Ubar = C*(.8-x)/(4*x)
        r_x = (1-W/3.2)*C
        L_Ubar, L_x = -(1/R+1/.8)*r_Ubar, 1/x+1/.8-(1/R+1/.8)*r_x
        s_x = R_x/v
        if min(W, W_x, C, s_x) <= 0:
            raise ValueError("reference wall leaves the positive, future-oriented gauge")
        eta_x = -.5*W_xx/W_x
        row = {"x": x, "U": U, "V": V, "W": W, "W_x": W_x, "W_xx": W_xx,
            "eta_plus": -.5*math.log(W_x), "eta_plus_x": eta_x,
            "r": R, "r_x_total": R_x, "r_xx_total": R_xx,
            "C": C, "L": L, "L_x_total": L_x_total,
            "r_Ubar": r_Ubar, "r_x": r_x, "log_C_Ubar": L_Ubar, "log_C_x": L_x,
            "c": 0., "c_Ubar": 0., "c_x": 0.,
            "s": float(integrals[1]), "s_x": s_x, "x_s": 1/s_x, "Ubar_s": W_x/s_x,
            "T": float(integrals[2]), "T_x": bm*s_x/(1-R*R), "T_s": bm/(1-R*R),
            "r_s": v, "beta_minus": bm, "beta_plus": bp,
            "inverse_iterations": iterations, "inverse_residual": residual,
            "normalization_defect": C*W_x/(s_x*s_x)-1}
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError("reference wall derived quantities must be finite")
        self.maximum_inverse_iterations = max(self.maximum_inverse_iterations, iterations)
        self.maximum_inverse_residual = max(self.maximum_inverse_residual, abs(residual))
        return row

    def sample(self, x):
        """Return scalar-shaped or one-dimensional rows on the closed interval."""
        query = _array(x, "reference wall query")
        if (query.ndim > 1 or query.size > MAX_QUERY_POINTS
                or np.any(query < self._lower) or np.any(query > self._upper)):
            raise ValueError("reference query exits its original interval or batch limit")
        self._charge(int(query.size))
        self.query_points += int(query.size)
        rows = []
        for point in query.reshape(-1):
            key = float(point)
            if key in self._cache:
                self.cache_hits += 1
                row = self._cache.pop(key)
            else:
                row = self._evaluate(key)
            self._cache[key] = row
            if len(self._cache) > MAX_CACHE_POINTS:
                self._cache.popitem(last=False)
            rows.append(row)
        self.budget.check()
        return {key: np.array([row[key] for row in rows]).reshape(query.shape) for key in rows[0]}

    def summary(self):
        return {"schema": "nsc-atlas-reference-wall-gauge-v1", "domain": self.domain,
                "quadrature_order": self.order, "evaluations": self.evaluations,
                "quadrature_nodes": self.quadrature_nodes, "query_points": self.query_points,
                "cache_hits": self.cache_hits, "inverse_calls": self.inverse_calls,
                "maximum_inverse_iterations": self.maximum_inverse_iterations,
                "maximum_inverse_residual": self.maximum_inverse_residual,
                "endpoint_integral_defect": self.endpoint_integral_defect,
                "clock_origin": "T=s=0 at the quiet x=.4, r=.45 wall anchor",
                "comparison_is_rigorous_enclosure": False,
                "source_integrations": 0, "regional_marches": 0}
