"""Finite linear Einstein propagation on the genuine exterior rectangle.

The leading scalar is supplied by an independently solved preparation. The
metric uses y=r0*R_E1 and z=delta log(C_E*r_E), whose mixed equations contain
no spatial derivatives. No metric is continued past the original wall.
"""
from dataclasses import dataclass, field as dataclass_field
import math
import time

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_characteristic_initialization import source_jet, validate_geometry
from .nsc_atlas_field_preparation import radius_from_null


SCHEMA = "nsc-atlas-bulk-response-v1"
CPU_SECONDS = 120
MAX_RHS = 12000
MAX_METRIC_VALUES = 30_000_000
M = .4
A_RADIUS = .8
GRID_SIZES = (33, 65, 129, 257)


@dataclass
class BulkBudget:
    cpu_seconds: float = CPU_SECONDS
    started: float = dataclass_field(default_factory=time.process_time)
    metric_values: int = 0
    source_nodes: int = 0
    rhs_calls: int = 0
    external: object = None

    def __post_init__(self):
        self.cpu_seconds = _number(self.cpu_seconds, "CPU budget")
        self.started = _number(self.started, "budget start")
        if not 0 < self.cpu_seconds <= CPU_SECONDS:
            raise ValueError("CPU budget must be positive and at most 120 seconds")
        for name in ("metric_values", "source_nodes", "rhs_calls"):
            self._count(getattr(self, name))

    @staticmethod
    def _count(value):
        if (isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer)) or value < 0):
            raise ValueError("budget counts must be nonnegative integers")
        return int(value)

    def check(self, *, metric=0, source=0, rhs=False):
        if not isinstance(rhs, (bool, np.bool_)):
            raise ValueError("RHS accounting flag must be boolean")
        self.metric_values += self._count(metric)
        self.source_nodes += self._count(source)
        self.rhs_calls += int(rhs)
        if (self.metric_values+self.source_nodes > MAX_METRIC_VALUES
                or self.rhs_calls > MAX_RHS
                or time.process_time()-self.started > self.cpu_seconds):
            raise RuntimeError("bulk response cumulative compute budget exhausted")
        if self.external is not None:
            self.external.check(0)


def _number(value, label):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(label+" must be a finite real number")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(label+" must be a finite real number") from error
    if not math.isfinite(value):
        raise ValueError(label+" must be a finite real number")
    return value


def background(U, V):
    """Independent smooth Schwarzschild-IV radius, lapse and analytic jets."""
    U, V = np.broadcast_arrays(np.asarray(U, dtype=float), np.asarray(V, dtype=float))
    r = radius_from_null(U, V)
    Fp = -r*np.exp(r/A_RADIUS)/A_RADIUS**2
    rU, rV = V/Fp, U/Fp
    C = 4*A_RADIUS**3*np.exp(-r/A_RADIUS)/r
    log_slope = -(1/r+1/A_RADIUS)
    rUV = -A_RADIUS**4*np.exp(-r/A_RADIUS)/r**3
    return {"r": r, "r_U": rU, "r_V": rV, "C": C,
            "logC_U": log_slope*rU, "logC_V": log_slope*rV,
            "r_UV": rUV, "r_UU": log_slope*rU*rU,
            "r_VV": log_slope*rV*rV}


def source_traces(V, geometry, *, order=64, source_scale=1., budget=None):
    """Exact-source jets with bounded moments, in the original (U,V) gauge.

    source_scale=0 is the explicit zero-source software control. The physical
    preparation is source_scale=1; this number is not a chosen finite epsilon.
    """
    V = _number(V, "source coordinate")
    if not isinstance(order, (int, np.integer)) or isinstance(order, (bool, np.bool_)) or order not in (32, 64):
        raise ValueError("source moments use order 32 or 64")
    source_scale = _number(source_scale, "source scale")
    if source_scale not in (0., 1.):
        raise ValueError("source scale is the original 1 or zero-source control 0")
    if not geometry["Va"] <= V <= geometry["Vb"]:
        raise ValueError("source trace exits the fixed characteristic")
    g = background(geometry["U0"], V)
    x, xV = float(g["r"]), float(g["r_V"])
    # Preserve exact source endpoint membership despite Lambert roundoff.
    x = min(geometry["x_last"], max(geometry["x_corner"], x))
    p, px, pxx = (float(value) for value in source_jet(x, geometry, source_scale))
    lower, support_upper = geometry["support_x"]
    upper = min(x, support_upper)
    I0 = I1 = I2 = Im = 0.
    if upper > lower and source_scale:
        nodes, weights = np.polynomial.legendre.leggauss(order)
        points = (lower+upper)/2+(upper-lower)*nodes/2
        weights = weights*(upper-lower)/2
        if budget is not None:
            budget.check(source=order)
        values, derivatives, _ = source_jet(points, geometry, source_scale)
        I0 = float(weights @ (values*values))
        I1 = float(weights @ (points*derivatives*derivatives))
        I2 = float(weights @ ((x-points)*points*derivatives*derivatives))
        Im = float(weights @ (points*(A_RADIUS-points)*derivatives*derivatives))
    h = -I2/2+x*p*p/12-I0/6
    hx = -I1/2+x*p*px/6-p*p/12
    H = h-x*p*p/12
    Hx = hx-p*p/12-x*p*px/6
    y = x*H
    yV = (H+x*Hx)*xV
    z = h/x-p*p/4
    zV = (hx/x-h/x**2-p*px/2)*xV
    muE = -Im/4
    # Z=C/(2*r_U); the original source gauge differs from affine E gauge.
    zetaE = h/x-p*p/4-I0/(6*(A_RADIUS-x))
    kE = -float(g["r_U"])/float(g["C"])*zetaE
    return {"x": x, "x_V": xV, "p": p, "p_V": px*xV,
            "p_x": px, "p_xx": pxx, "h_J": h, "h_J_x": hx,
            "H_E": H, "H_E_V": Hx*xV, "c_E": -p*p/6,
            "y": y, "z": z, "y_V": yV, "z_V": zV,
            "mu_E": muE, "zeta_E": zetaE, "k_E": kE,
            "integral_p_squared": I0}


def mixed_sources(r, C, y, z, p_U, p_V):
    """Derivative-free linear mixed equations, with canonical scalar source."""
    return C*(y/r**2-z)/4, C*(z/r**2-3*y/r**4)/4-p_U*p_V


def _prefix(values, coordinates):
    return np.r_[0., np.cumsum((values[:-1]+values[1:])*np.diff(coordinates)/2)]


def _field_values(prepared_field, U, t, order):
    values = prepared_field.bulk(U, float(t), order=order)
    result = {}
    for key in ("p", "p_U", "p_V"):
        if key not in values:
            raise ValueError("prepared field is missing "+key)
        array = np.asarray(values[key], dtype=float)
        if array.shape != U.shape or not np.all(np.isfinite(array)):
            raise ValueError("prepared field has invalid "+key+" samples")
        result[key] = array
    return result


class BulkSolution:
    """Dense finite rectangle solution; samples never extend into the wall strip."""
    def __init__(self, *, U, geometry, prepared_field, segments, cuts, rate,
                 budget, source_order, field_order, source_scale):
        self.U, self.geometry, self.prepared_field = U, dict(geometry), prepared_field
        self.domain = {key: geometry[key] for key in ("U0", "Ua", "Va", "Vb")}
        self.segments, self.cuts, self._rate = segments, cuts, rate
        self.budget = budget
        self.source_order, self.field_order, self.source_scale = source_order, field_order, source_scale
        self.nodes = len(U)

    def _state(self, t):
        t = _number(t, "normalized time")
        if not 0 <= t <= 1:
            raise ValueError("metric query exits the solved rectangle time interval")
        index = min(np.searchsorted(self.cuts, t, side="right")-1, len(self.segments)-1)
        return t, self.segments[index](t).reshape(6, self.nodes)

    def sample(self, t, *, independent_order=None):
        """Metric, analytic first jets and independent integral-constraint defects.

        An optional Gauss integral checks endpoint mass with a different
        spatial rule. Dense output is transient and is not saved by this class.
        """
        t, state = self._state(t)
        _, diagnostic = self._rate(t, state.ravel(), counted=False, diagnostics=True)
        g, source, scalar = diagnostic["background"], diagnostic["source"], diagnostic["scalar"]
        r, rU, rV, C = (g[key] for key in ("r", "r_U", "r_V", "C"))
        y, z = state[0]+source["y"], state[1]+source["z"]
        yU, zU = state[2], state[3]
        yV, zV = diagnostic["y_V"], diagnostic["z_V"]
        H, cE = y/r, z-y/r**2
        HU, HV = yU/r-y*rU/r**2, yV/r-y*rV/r**2
        cU = zU-yU/r**2+2*y*rU/r**3
        cV = zV-yV/r**2+2*y*rV/r**3
        p, pU, pV = (scalar[key] for key in ("p", "p_U", "p_V"))
        h = H+r*p*p/12
        hU, hV = HU+(rU*p*p+2*r*p*pU)/12, HV+(rV*p*p+2*r*p*pV)/12
        cJ, cJU, cJV = cE+p*p/6, cU+p*pU/3, cV+p*pV/3
        muE = M*H/r+2*r*(rV*HU+rU*HV-rU*rV*cE)/C
        muJ = M*h/r+2*r*(rV*hU+rU*hV-rU*rV*cJ)/C
        frame = M*p*p/12+r*r*p*(rU*pV+rV*pU)/(3*C)
        muU = source["mu_E"]+_prefix(-r*r*rV*pU*pU/C, self.U)
        kE = (HU-rU*cE)/C
        kExpected = source["k_E"]+_prefix(-r*pU*pU/(2*C), self.U)
        lE = (HV-rV*cE)/C
        zeta = cE-HU/rU
        zeta_radial = source["zeta_E"]+_prefix(r*pU*pU/(2*rU), self.U)
        result = {
            "U": self.U.copy(), "V": diagnostic["V"], "r0": r,
            "y": y, "z": z, "y_U": yU, "z_U": zU, "y_V": yV, "z_V": zV,
            "H_E": H, "H_E_U": HU, "H_E_V": HV,
            "c_E": cE, "c_E_U": cU, "c_E_V": cV,
            "h_J": h, "h_J_U": hU, "h_J_V": hV,
            "c_J": cJ, "c_J_U": cJU, "c_J_V": cJV,
            "mu_E": muE, "mu_J": muJ, "frame_mass_difference": frame,
            "mu_E_flux_U": muU, "mu_E_flux_V": state[5],
            "mass_defect_U": muE-muU, "mass_defect_V": muE-state[5],
            "frame_mass_defect": muJ-muE-frame,
            "null_U_integral_defect": kE-kExpected,
            "null_V_integral_defect": lE-state[4],
            "zeta_E": zeta, "zeta_E_radial": zeta_radial,
            "areal_null_lapse_defect": zeta-zeta_radial,
            "p": p, "p_U": pU, "p_V": pV,
            "source": source,
        }
        if independent_order is not None:
            if (not isinstance(independent_order, (int, np.integer))
                    or isinstance(independent_order, (bool, np.bool_))
                    or independent_order not in (32, 64)):
                raise ValueError("independent mass quadrature uses order 32 or 64")
            nodes, weights = np.polynomial.legendre.leggauss(independent_order)
            U0, Ua = self.domain["U0"], self.domain["Ua"]
            points = (U0+Ua)/2+(Ua-U0)*nodes/2
            bg = background(points, diagnostic["V"])
            sf = _field_values(self.prepared_field, points, t, self.field_order)
            self.budget.check(metric=independent_order)
            integrand = -bg["r"]**2*bg["r_V"]*sf["p_U"]**2/bg["C"]
            independent = source["mu_E"]+(Ua-U0)/2*float(weights @ integrand)
            result["mu_E_flux_U_gauss_endpoint"] = independent
            result["mass_defect_U_gauss_endpoint"] = float(muE[-1]-independent)
        return result

    def summary(self, times=(0., .25, .5, .75, 1.), *, independent_order=32):
        """Compact finite coefficients and defects, with no dense trial arrays."""
        keys = ("mass_defect_U", "mass_defect_V", "null_U_integral_defect",
                "null_V_integral_defect", "areal_null_lapse_defect", "frame_mass_defect")
        maxima = dict.fromkeys(keys, 0.)
        rows = []
        for t in times:
            sampled = self.sample(t, independent_order=independent_order)
            for key in keys:
                maxima[key] = max(maxima[key], float(np.max(np.abs(sampled[key]))))
            rows.append({"t": float(t), **{key: float(sampled[key][-1]) for key in
                         ("h_J", "c_J", "h_J_U", "h_J_V", "c_J_U", "c_J_V", "mu_E", "mu_J")},
                         "mass_defect_U_gauss": sampled["mass_defect_U_gauss_endpoint"]})
        return {"schema": SCHEMA, "domain": self.domain, "metric_nodes": self.nodes,
                "source_order": self.source_order, "field_order": self.field_order,
                "rhs_calls": self.budget.rhs_calls, "metric_values": self.budget.metric_values,
                "source_nodes": self.budget.source_nodes,
                "process_cpu_seconds": time.process_time()-self.budget.started,
                "maximum_absolute_defects": maxima, "transverse_leg_samples": rows,
                "scope": "formal exterior rectangle; U=Ua data for subsequent wall matching"}


def solve_bulk_rectangle(prepared_field, geometry, *, nodes=33, source_order=64,
                         field_order=32, source_scale=1., tolerance=2e-9, budget=None):
    """Propagate the two linear metric fields from their genuine initial legs."""
    validate_geometry(geometry)
    if (not isinstance(nodes, (int, np.integer))
            or isinstance(nodes, (bool, np.bool_)) or nodes not in GRID_SIZES):
        raise ValueError("metric grid must have 33, 65, 129 or 257 nodes")
    if any(not isinstance(value, (int, np.integer)) or isinstance(value, (bool, np.bool_))
           or value not in (32, 64) for value in (source_order, field_order)):
        raise ValueError("quadrature orders must be 32 or 64")
    source_scale = _number(source_scale, "source scale")
    if source_scale not in (0., 1.):
        raise ValueError("source scale is the original 1 or zero-source control 0")
    tolerance = _number(tolerance, "ODE tolerance")
    if not 1e-12 <= tolerance <= 1e-6:
        raise ValueError("ODE tolerance is outside the admitted useful range")
    domain = getattr(prepared_field, "domain", None)
    if not isinstance(domain, dict) or any(
            key not in domain or not math.isclose(float(domain[key]), geometry[key], rel_tol=0, abs_tol=2e-12)
            for key in ("U0", "Ua", "Va", "Vb")):
        raise ValueError("field and metric characteristic domains do not agree")
    if budget is None:
        external = getattr(prepared_field, "budget", None)
        budget = BulkBudget(external=external)
    budget.check()
    U = np.linspace(geometry["U0"], geometry["Ua"], nodes)
    Va, Vb = geometry["Va"], geometry["Vb"]
    span = Vb-Va

    # The data comparison is a physical compatibility check, not a source fit.
    for t in (0., .5, 1.):
        V = Va+span*t
        src = source_traces(V, geometry, order=source_order, source_scale=source_scale, budget=budget)
        given = _field_values(prepared_field, U[:1], t, field_order)
        if (not np.isclose(given["p"][0], src["p"], rtol=2e-8, atol=2e-10)
                or not np.isclose(given["p_V"][0], src["p_V"], rtol=2e-8, atol=2e-10)):
            raise ValueError("prepared scalar and gravitational source data are incompatible")
    quiet = _field_values(prepared_field, U, 0., field_order)
    if any(np.any(np.abs(quiet[key]) > 2e-12) for key in ("p", "p_U", "p_V")):
        raise ValueError("bottom field data are not quiet")

    def rate(t, vector, *, counted=True, diagnostics=False):
        budget.check(metric=nodes, rhs=counted)
        V = min(Vb, max(Va, Va+span*float(t)))
        state = np.asarray(vector).reshape(6, nodes)
        src = source_traces(V, geometry, order=source_order, source_scale=source_scale, budget=budget)
        bg = background(U, V)
        sf = _field_values(prepared_field, U, min(1., max(0., t)), field_order)
        if (not np.isclose(sf["p"][0], src["p"], rtol=2e-8, atol=2e-10)
                or not np.isclose(sf["p_V"][0], src["p_V"], rtol=2e-8, atol=2e-10)):
            raise ValueError("prepared scalar and gravitational source data are incompatible")
        r, C = bg["r"], bg["C"]
        y, z = state[0]+src["y"], state[1]+src["z"]
        Fy, Fz = mixed_sources(r, C, y, z, sf["p_U"], sf["p_V"])
        y_rate, z_rate = _prefix(Fy, U), _prefix(Fz, U)
        lv_rate = -r*sf["p_V"]**2/(2*C)
        mass_rate = -r*r*bg["r_U"]*sf["p_V"]**2/C
        derivative = span*np.array([y_rate, z_rate, Fy, Fz, lv_rate, mass_rate]).ravel()
        if diagnostics:
            return derivative, {"V": V, "background": bg, "source": src, "scalar": sf,
                                "y_V": y_rate+src["y_V"], "z_V": z_rate+src["z_V"]}
        return derivative

    cuts = np.array([0., .25, .75, 1.])
    state = np.zeros(6*nodes)
    segments = []
    for lower, upper in zip(cuts[:-1], cuts[1:]):
        solved = solve_ivp(rate, (lower, upper), state, method="DOP853", dense_output=True,
                           rtol=tolerance, atol=1e-11, first_step=min(1/32, upper-lower), max_step=1/32)
        if not solved.success or solved.t[-1] != upper:
            raise RuntimeError("bulk rectangle propagation did not reach its finite endpoint")
        state = solved.y[:, -1]
        if not np.all(np.isfinite(state)):
            raise RuntimeError("bulk rectangle propagation produced nonfinite coefficients")
        segments.append(solved.sol)
    budget.check()
    return BulkSolution(U=U, geometry=geometry, prepared_field=prepared_field,
                        segments=segments, cuts=cuts, rate=rate, budget=budget,
                        source_order=source_order, field_order=field_order, source_scale=source_scale)
