"""Two-level spherical nested initial slice with supplied vacuum curvatures.

The quiet Einstein and conformal-scalar action has xi=1/6 and phi=0, and both
walls use the same tension k. Regional Hubble values are supplied. This is one
initial slice, not a derived Lambda, a formation history, or a common sibling
parent. Static horizon clocks and a joint worldtube spacetime are not built.
"""
from __future__ import annotations

import math
import time

import numpy as np
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_connection_slice import schwarzschild_bridge, smooth_step
from recursive_horizons.nsc_atlas_geometric_average import KAPPA, scalar_average_sources


IDENTITY = 1e-8
QUADRATURE_FACTOR = 1e-7
SAMPLE_LIMIT = 100_000
CPU_SECONDS = 30
INNER_WALL = 0.8
HOST_STATIC_LEFT = 0.1
HOST_STATIC_RIGHT = 0.5
HOST_EXPANDING_START = 0.9
OUTER_WALL = 1.5
OUTER_TAPER_STOP = 1.0
OUTER_CUT = 4.0
HOST_MASS = 0.01
OUTER_MASS = 0.4
TENSION = 2.0
INTERIOR_HUBBLE = 1.0
EXTERIOR_HUBBLE = 0.0


class EvaluationBudget:
    def __init__(self, limit=SAMPLE_LIMIT):
        self.limit = int(limit)
        self.count = 0

    def add(self, number):
        self.count += int(number)
        if self.count > self.limit:
            raise RuntimeError(
                f"geometric evaluation budget exhausted at {self.count} points")


def _fail(message):
    raise ValueError(message)


def _require(name, residual, scale=1.0):
    value = float(residual)
    if not math.isfinite(value) or abs(value) > IDENTITY * scale:
        _fail(f"{name} residual {value} exceeds the {IDENTITY} identity gate")
    return value


def metric_lapse(radius, mass, hubble):
    radius = np.asarray(radius, dtype=float)
    return 1 - 2*mass/radius - hubble**2 * radius**2


def metric_lapse_radius(radius, mass, hubble):
    radius = np.asarray(radius, dtype=float)
    return 2*mass/radius**2 - 2*hubble**2 * radius


def metric_lapse_second(radius, mass, hubble):
    radius = np.asarray(radius, dtype=float)
    return -4*mass/radius**3 - 2*hubble**2


def kottler_horizons(mass, hubble, counter=None):
    """Black-hole and cosmological roots of f. Mass is not adjusted to the root."""
    if not math.isfinite(mass) or mass <= 0 or not math.isfinite(hubble) or hubble < 0:
        _fail("Kottler black-hole roots require finite positive mass and nonnegative H")
    if hubble == 0:
        if mass == 0:
            _fail("flat exterior has no black-hole root")
        return {"rb": 2*mass, "rc": None, "nariai_parameter": 0.0,
                "root_f_residual": 0.0, "cosmological_root_f_residual": None,
                "mass_identity_residual": 0.0}

    mu = mass*hubble
    nariai = 27*mu**2
    if nariai >= 1:
        _fail("Kottler data are not sub-Nariai")

    # Solve in x=H*R, preserving units and root precision under ell rescaling.
    def lapse(x):
        if counter is not None:
            counter.add(1)
        return float(1-2*mu/x-x*x)

    settings = dict(xtol=np.nextafter(0.,1.), rtol=4*np.finfo(float).eps)
    xb = float(brentq(lapse, 2*mu, 3*mu, **settings))
    xc = float(brentq(lapse, 1/math.sqrt(3), 1., **settings))
    rb,rc = xb/hubble,xc/hubble
    return {"rb": rb, "rc": rc, "nariai_parameter": nariai,
            "root_f_residual": lapse(xb),
            "cosmological_root_f_residual": lapse(xc),
            "mass_identity_residual": 2*mass - (rb - hubble**2 * rb**3)}


def bridge_amplitude(radius, throat, hubble):
    """A(R) in f=(R-rb)A(R), finite and nonzero at the throat."""
    radius = np.asarray(radius, dtype=float)
    return (1 - hubble**2 * (radius**2 + radius*throat + throat**2))/radius


def vacuum_slice_state(radius, mass, hubble, velocity, velocity_radius, radial_prime):
    """Bounded slice geometry for one prescribed velocity. K has no wall delta."""
    radius = np.asarray(radius, dtype=float)
    velocity = np.broadcast_to(np.asarray(velocity, dtype=float), radius.shape)
    velocity_radius = np.broadcast_to(np.asarray(velocity_radius, dtype=float), radius.shape)
    radial_prime = np.broadcast_to(np.asarray(radial_prime, dtype=float), radius.shape)
    if np.any(radius <= 0):
        _fail("areal radius must stay positive away from the regular centre")
    curvature_radius = mass/radius**2 - hubble**2 * radius + velocity*velocity_radius
    radial_curvature = -velocity_radius
    angular_curvature = -velocity/radius
    theta = velocity_radius + 2*velocity/radius
    shear_squared = (velocity_radius - velocity/radius)**2 / 3
    scalar3 = -4*curvature_radius/radius + 2*(1-radial_prime**2)/radius**2
    angular_derivative = -velocity_radius/radius + velocity/radius**2
    return {
        "radius": radius, "velocity": velocity, "velocity_radius": velocity_radius,
        "Rprime": radial_prime, "Rsecond": curvature_radius,
        "Kll": radial_curvature, "Ktheta": angular_curvature, "theta": theta,
        "shear_squared": shear_squared, "scalar3": scalar3,
        "hamiltonian": scalar3 + 4*radial_curvature*angular_curvature
        + 2*angular_curvature**2 - 6*hubble**2,
        "momentum": 2*radial_prime/radius*(radial_curvature-angular_curvature)
        - 2*radial_prime*angular_derivative,
        "misner_sharp": radius*(1+velocity**2 - radial_prime**2)/2,
        "target_mass": mass+hubble**2*radius**3/2,
        "f": metric_lapse(radius, mass, hubble),
    }


def _unwrap(state, scalar):
    if not scalar:
        return state
    return {key: float(np.asarray(value).reshape(-1)[0]) for key, value in state.items()}


def flat_cap_state(radius, hubble):
    """Flat expanding de Sitter cap, including the regular centre."""
    scalar = np.isscalar(radius)
    radius = np.atleast_1d(np.asarray(radius, dtype=float))
    if np.any(radius < 0) or not math.isfinite(hubble):
        _fail("cap radius and Hubble value must be finite and nonnegative")
    velocity = hubble*radius
    state = {
        "radius": radius, "velocity": velocity, "velocity_radius": np.full(radius.shape, hubble),
        "Rprime": np.ones(radius.shape), "Rsecond": np.zeros(radius.shape),
        "Kll": np.full(radius.shape, -hubble), "Ktheta": np.full(radius.shape, -hubble),
        "theta": np.full(radius.shape, 3*hubble), "shear_squared": np.zeros(radius.shape),
        "scalar3": np.zeros(radius.shape), "hamiltonian": np.zeros(radius.shape),
        "momentum": np.zeros(radius.shape), "misner_sharp": hubble**2 * radius**3 / 2,
        "f": 1 - hubble**2 * radius**2,
    }
    return _unwrap(state, scalar)


def positive_sheet_velocity(radius, mass, hubble):
    radius = np.asarray(radius, dtype=float)
    squared = 2*mass/radius + hubble**2 * radius**2
    if np.any(squared <= 0):
        _fail("expanding sheet velocity is not real")
    velocity = np.sqrt(squared)
    return velocity, (hubble**2 * radius - mass/radius**2)/velocity


def positive_sheet_state(radius, mass, hubble):
    scalar = np.isscalar(radius)
    velocity, velocity_radius = positive_sheet_velocity(radius, mass, hubble)
    state = vacuum_slice_state(radius, mass, hubble, velocity, velocity_radius,
                               np.ones(np.shape(velocity)))
    return _unwrap(state, scalar)


def taper_velocity(radius, stop, start, wall_velocity):
    span = start-stop
    if span <= 0:
        _fail("taper must run from its stop to the wall")
    step, step_radius = smooth_step((np.asarray(radius, dtype=float)-stop)/span)
    return wall_velocity*step, wall_velocity*step_radius/span


def transition_velocity(radius, low, high, mass, hubble):
    span = high-low
    if span <= 0:
        _fail("transition interval must increase in areal radius")
    step, step_radius = smooth_step((np.asarray(radius, dtype=float)-low)/span)
    sheet, sheet_radius = positive_sheet_velocity(radius, mass, hubble)
    return step*sheet, step_radius/span*sheet + step*sheet_radius


def signed_slope(radius, mass, hubble, velocity, orientation):
    if orientation not in (-1, 1):
        _fail("slice orientation must be -1 or +1")
    inside = metric_lapse(radius, mass, hubble) + np.asarray(velocity, dtype=float)**2
    if np.any(inside < 0):
        _fail("prescribed velocity leaves the spacelike slice")
    return orientation*np.sqrt(inside)


def bridge_chart(y, throat, hubble, mass):
    """Regular throat chart R=rb+y^2. The supplied mass is not replaced."""
    scalar = np.isscalar(y)
    y = np.atleast_1d(np.asarray(y, dtype=float))
    radius = throat + y**2
    amplitude = bridge_amplitude(radius, throat, hubble)
    if np.any(amplitude <= 0) or not np.all(np.isfinite(amplitude)):
        _fail("bridge amplitude left the positive chart")
    slope = y*np.sqrt(amplitude)
    state = vacuum_slice_state(radius, mass, hubble, np.zeros(y.shape), np.zeros(y.shape), slope)
    state["y"] = y
    state["amplitude"] = amplitude
    state["dl_dy"] = 2/np.sqrt(amplitude)
    return _unwrap(state, scalar)


def tension_wall(radius, mass_minus, hubble_minus, mass_plus, hubble_plus, tension):
    """One static-chart tension wall and its flat initial slice boost."""
    if min(radius, tension) <= 0 or mass_minus < 0 or mass_plus <= 0:
        _fail("wall radius, tension and exterior mass must be positive")
    f_minus = float(metric_lapse(radius, mass_minus, hubble_minus))
    f_plus = float(metric_lapse(radius, mass_plus, hubble_plus))
    f_minus_radius = float(metric_lapse_radius(radius, mass_minus, hubble_minus))
    f_plus_radius = float(metric_lapse_radius(radius, mass_plus, hubble_plus))
    numerator = f_minus-f_plus + (tension*radius)**2
    beta_minus = numerator/(2*tension*radius)
    beta_plus = beta_minus-tension*radius
    numerator_radius = f_minus_radius-f_plus_radius + 2*tension**2 * radius
    beta_radius = (numerator_radius*radius-numerator)/(2*tension*radius**2)
    motion = beta_minus**2 - f_minus
    if motion <= 0 or beta_minus == 0 or beta_plus == 0:
        _fail("wall worldtube is not timelike with nonzero signed normals")
    areal_rate = math.sqrt(motion)
    areal_acceleration = beta_minus*beta_radius - f_minus_radius/2
    temporal_minus = (areal_acceleration+f_minus_radius/2)/beta_minus
    temporal_plus = (areal_acceleration+f_plus_radius/2)/beta_plus
    angular_minus = beta_minus/radius
    angular_plus = beta_plus/radius
    trace_minus = temporal_minus+2*angular_minus
    trace_plus = temporal_plus+2*angular_plus
    mean_trace = (trace_minus+trace_plus)/2
    predicted = (mass_minus + (hubble_minus**2-hubble_plus**2-tension**2)*radius**3/2
                 + tension*radius**2 * beta_minus)
    sigma = tension/(4*math.pi)
    normal_force = sigma*mean_trace + (-3*hubble_plus**2 + 3*hubble_minus**2)/KAPPA
    velocity_minus = math.sqrt(2*mass_minus/radius + hubble_minus**2 * radius**2)
    boost = (beta_minus+areal_rate)/(1+velocity_minus)
    if boost <= 0:
        _fail("wall boost is not a positive Lorentz factor")
    gamma = (boost+1/boost)/2
    gamma_velocity = (boost-1/boost)/2
    speed = gamma_velocity/gamma
    radial_plus = 1-tension*radius*gamma
    velocity_plus = velocity_minus+tension*radius*gamma_velocity
    return {
        "radius": radius, "mass_minus": mass_minus, "mass_plus": mass_plus,
        "hubble_minus": hubble_minus, "hubble_plus": hubble_plus, "tension": tension,
        "f_minus": f_minus, "f_plus": f_plus, "beta_minus": beta_minus,
        "beta_plus": beta_plus, "Rdot": areal_rate, "Rddot": areal_acceleration,
        "worldtube_Ktau_minus": temporal_minus, "worldtube_Ktau_plus": temporal_plus,
        "worldtube_Ktheta_minus": angular_minus, "worldtube_Ktheta_plus": angular_plus,
        "angular_jump": angular_plus-angular_minus, "temporal_jump": temporal_plus-temporal_minus,
        "trace_jump": trace_plus-trace_minus, "mean_trace": mean_trace,
        "mass_plus_predicted": predicted, "mass_map_residual": predicted-mass_plus,
        "normal_force_residual": normal_force, "sigma": sigma,
        "V_minus": velocity_minus, "V_plus": velocity_plus,
        "Rprime_minus": 1.0, "Rprime_plus": radial_plus,
        "gamma": gamma, "gamma_v": gamma_velocity, "v": speed,
        "lorentz_residual": gamma**2-gamma_velocity**2-1,
        "minus_rate_residual": areal_rate-gamma*(velocity_minus+speed),
        "minus_beta_residual": beta_minus-gamma*(1+speed*velocity_minus),
        "plus_rate_residual": areal_rate-(gamma*velocity_plus+gamma_velocity*radial_plus),
        "plus_beta_residual": beta_plus-(gamma*radial_plus+gamma_velocity*velocity_plus),
        "energy": sigma*gamma, "momentum": sigma*gamma_velocity,
    }


def gauss_integral(integrand, lower, upper, order, counter):
    if order not in (32, 64):
        _fail("nested-slice quadrature orders are 32 and 64")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    sample = (upper+lower)/2 + (upper-lower)/2 * nodes
    counter.add(sample.size)
    values = np.asarray(integrand(sample), dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if values.shape != (sample.size, values.shape[1]):
        _fail("integrand must return one row per sample")
    return np.atleast_1d((upper-lower)/2 * weights @ values)


def _columns(proper, volume, theta, theta_squared, shear, scalar3):
    return np.column_stack((proper, volume, theta, theta_squared, shear, scalar3))


def _watch(state, name, length_scale):
    for key in ("hamiltonian", "momentum"):
        worst = float(np.max(np.abs(state[key])))*length_scale**2
        _require(f"{name} {key}", worst)
    mass_gap = state["misner_sharp"]-(state.get("target_mass", 0))
    # target mass is checked by callers that know the region; Misner-Sharp is
    # compared here only when the state carries it explicitly.
    if "target_mass" in state:
        _require(f"{name} Misner-Sharp", float(np.max(np.abs(mass_gap)))/length_scale, 1)
    for key in ("Kll", "Ktheta"):
        if not np.all(np.isfinite(state[key])):
            _fail(f"{name} {key} is not bounded")


def _radius_pack(radius, mass, hubble, velocity, velocity_radius, orientation, length_scale, name):
    slope = signed_slope(radius, mass, hubble, velocity, orientation)
    state = vacuum_slice_state(radius, mass, hubble, velocity, velocity_radius, slope)
    _watch(state, name, length_scale)
    measure = 4*math.pi*radius**2 / state["Rprime"]
    return _columns(1/state["Rprime"], measure, measure*state["theta"],
                    measure*state["theta"]**2, measure*state["shear_squared"],
                    measure*state["scalar3"])


def _integrate_pair(integrand, lower, upper, counter):
    coarse = gauss_integral(integrand, lower, upper, 32, counter)
    fine = gauss_integral(integrand, lower, upper, 64, counter)
    return coarse, fine


def _quadrature_gap(coarse, fine, name):
    gap = float(abs(coarse-fine))
    limit = QUADRATURE_FACTOR*max(1.0, abs(float(fine)))
    if gap > limit:
        _fail(f"{name} Gauss32/64 difference {gap} exceeds indicator {limit}")
    return gap


def _positive_concave_lapse(lower, upper, mass, hubble, name, counter):
    sample = np.array([lower, (lower+upper)/2, upper], dtype=float)
    counter.add(sample.size)
    second = metric_lapse_second(sample, mass, hubble)
    if np.any(second >= 0):
        _fail(f"{name} lapse is not strictly concave")
    ends = metric_lapse(np.array([lower, upper]), mass, hubble)
    counter.add(2)
    if float(np.min(ends)) <= 0:
        _fail(f"{name} endpoint lapse is not positive, so the taper is not outside its horizon")


def nested_slice_record(length_scale=1.0, *, include_partition=False):
    """Scalar milestones and moments through Rcut. The AE end remains analytic."""
    started = time.process_time()
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        _fail("length scale must be finite and positive")
    counter = EvaluationBudget()
    host = kottler_horizons(HOST_MASS*ell, INTERIOR_HUBBLE/ell, counter)
    exterior = kottler_horizons(OUTER_MASS*ell, EXTERIOR_HUBBLE, counter)
    _require("host root lapse", host["root_f_residual"])
    _require("host cosmological lapse", host["cosmological_root_f_residual"])
    _require("host mass identity", host["mass_identity_residual"]/ell)
    _require("exterior mass identity", exterior["mass_identity_residual"]/ell)
    if not 0 < host["rb"] < HOST_STATIC_LEFT*ell < HOST_STATIC_RIGHT*ell < host["rc"] < OUTER_WALL*ell:
        _fail("host roots do not separate the static bridge from the expanding patch")
    if exterior["rb"] != 2*OUTER_MASS*ell or not exterior["rb"] < OUTER_TAPER_STOP*ell < OUTER_CUT*ell:
        _fail("exterior root is not the Schwarzschild throat inside the retained cut")

    inner = tension_wall(INNER_WALL*ell, 0.0, INTERIOR_HUBBLE/ell, HOST_MASS*ell,
                         INTERIOR_HUBBLE/ell, TENSION/ell)
    outer = tension_wall(OUTER_WALL*ell, HOST_MASS*ell, INTERIOR_HUBBLE/ell, OUTER_MASS*ell,
                         EXTERIOR_HUBBLE, TENSION/ell)
    counter.add(2)
    for name, wall in (("inner", inner), ("outer", outer)):
        scaled = {
            "mass_map_residual": wall["mass_map_residual"]/ell,
            "normal_force_residual": wall["normal_force_residual"]*ell**2,
            "lorentz_residual": wall["lorentz_residual"],
            "minus_rate_residual": wall["minus_rate_residual"],
            "minus_beta_residual": wall["minus_beta_residual"],
            "plus_rate_residual": wall["plus_rate_residual"],
            "plus_beta_residual": wall["plus_beta_residual"],
        }
        for key, value in scaled.items():
            _require(f"{name} wall {key}", value)
        _require(f"{name} angular jump", (wall["angular_jump"]+wall["tension"])*ell)
        _require(f"{name} temporal jump", (wall["temporal_jump"]+wall["tension"])*ell)
        _require(f"{name} trace jump", (wall["trace_jump"]+3*wall["tension"])*ell)
        if wall["Rprime_plus"] >= 0 or abs(wall["v"]) >= 1:
            _fail(f"{name} wall does not have the outward negative slope")
    _require("inner mean trace", inner["mean_trace"]*ell)
    _require("outer mean trace", (outer["mean_trace"]+0.75/ell)*ell)

    regions = []
    coarse_regions = []
    regional_comparisons = []

    def remember(name, lower, upper, integrand):
        coarse, fine = _integrate_pair(integrand, lower, upper, counter)
        for index, label in enumerate(("proper", "volume", "theta", "theta_squared", "shear", "scalar")):
            gap = _quadrature_gap(coarse[index], fine[index], f"{name} {label}")
            regional_comparisons.append((gap/max(1.,abs(float(fine[index]))),name+"."+label))
        regions.append((name, fine))
        coarse_regions.append((name,coarse))

    def cap_integrand(radius):
        radius = np.asarray(radius, dtype=float)
        measure = 4*math.pi*radius**2
        theta = np.full(radius.shape, 3/ell)
        return _columns(np.ones(radius.shape), measure, measure*theta, measure*theta**2,
                        np.zeros(radius.shape), np.zeros(radius.shape))

    def host_left(radius):
        velocity, velocity_radius = taper_velocity(radius, HOST_STATIC_LEFT*ell, INNER_WALL*ell,
                                                   inner["V_plus"])
        return _radius_pack(radius, HOST_MASS*ell, 1/ell, velocity, velocity_radius, -1, ell,
                            "host left taper")

    def host_bridge(y):
        state = bridge_chart(y, host["rb"], 1/ell, HOST_MASS*ell)
        _watch(state, "host bridge", ell)
        measure = 4*math.pi*state["radius"]**2 * state["dl_dy"]
        return _columns(state["dl_dy"], measure, measure*state["theta"], measure*state["theta"]**2,
                        measure*state["shear_squared"], measure*state["scalar3"])

    def host_transition(radius):
        velocity, velocity_radius = transition_velocity(
            radius, HOST_STATIC_RIGHT*ell, HOST_EXPANDING_START*ell, HOST_MASS*ell, 1/ell)
        return _radius_pack(radius, HOST_MASS*ell, 1/ell, velocity, velocity_radius, 1, ell,
                            "host right transition")

    def host_sheet(radius):
        velocity, velocity_radius = positive_sheet_velocity(radius, HOST_MASS*ell, 1/ell)
        state = vacuum_slice_state(radius, HOST_MASS*ell, 1/ell, velocity, velocity_radius,
                                   np.ones(np.shape(radius)))
        _watch(state, "host expanding sheet", ell)
        _require("host sheet unit slope identity",
                 float(np.max(np.abs(state["f"]+state["velocity"]**2-1))))
        measure = 4*math.pi*np.asarray(radius, dtype=float)**2
        return _columns(np.ones(np.shape(radius)), measure, measure*state["theta"],
                        measure*state["theta"]**2, measure*state["shear_squared"],
                        measure*state["scalar3"])

    def outer_taper(radius):
        velocity, velocity_radius = taper_velocity(radius, OUTER_TAPER_STOP*ell, OUTER_WALL*ell,
                                                   outer["V_plus"])
        return _radius_pack(radius, OUTER_MASS*ell, 0.0, velocity, velocity_radius, -1, ell,
                            "outer taper")

    def outer_bridge(y):
        state = bridge_chart(y, exterior["rb"], 0.0, OUTER_MASS*ell)
        _watch(state, "outer bridge", ell)
        measure = 4*math.pi*state["radius"]**2 * state["dl_dy"]
        return _columns(state["dl_dy"], measure, measure*state["theta"], measure*state["theta"]**2,
                        measure*state["shear_squared"], measure*state["scalar3"])

    _positive_concave_lapse(HOST_STATIC_LEFT*ell, INNER_WALL*ell, HOST_MASS*ell, 1/ell,
                            "host left taper", counter)
    _positive_concave_lapse(HOST_STATIC_RIGHT*ell, HOST_EXPANDING_START*ell, HOST_MASS*ell, 1/ell,
                            "host right transition", counter)
    _positive_concave_lapse(OUTER_TAPER_STOP*ell, OUTER_WALL*ell, OUTER_MASS*ell, 0.0,
                            "outer taper", counter)
    remember("cap", 0.0, INNER_WALL*ell, cap_integrand)
    remember("host_left", INNER_WALL*ell, HOST_STATIC_LEFT*ell, host_left)
    y_host_left = -math.sqrt(HOST_STATIC_LEFT*ell-host["rb"])
    y_host_right = math.sqrt(HOST_STATIC_RIGHT*ell-host["rb"])
    remember("host_bridge_left", y_host_left, 0.0, host_bridge)
    remember("host_bridge_right", 0.0, y_host_right, host_bridge)
    remember("host_transition", HOST_STATIC_RIGHT*ell, HOST_EXPANDING_START*ell, host_transition)
    remember("host_sheet_to_horizon", HOST_EXPANDING_START*ell, host["rc"], host_sheet)
    remember("host_sheet_to_wall", host["rc"], OUTER_WALL*ell, host_sheet)
    remember("outer_taper", OUTER_WALL*ell, OUTER_TAPER_STOP*ell, outer_taper)
    y_outer_left = -math.sqrt(OUTER_TAPER_STOP*ell-exterior["rb"])
    y_outer_right = math.sqrt(OUTER_CUT*ell-exterior["rb"])
    remember("outer_bridge_left", y_outer_left, 0.0, outer_bridge)
    remember("outer_bridge_right", 0.0, y_outer_right, outer_bridge)

    totals = np.sum([values for _, values in regions], axis=0)
    cap_exact = np.array([
        INNER_WALL*ell, 4*math.pi*(INNER_WALL*ell)**3/3,
        3/ell*4*math.pi*(INNER_WALL*ell)**3/3, 9/ell**2 * 4*math.pi*(INNER_WALL*ell)**3/3,
        0.0, 0.0])
    cap_quadrature = dict(regions)["cap"]
    for index, label in enumerate(("proper", "volume", "theta", "theta_squared")):
        _quadrature_gap(cap_quadrature[index], cap_exact[index], f"exact cap {label}")
    totals += cap_exact-cap_quadrature
    coarse_totals = np.sum([values for _,values in coarse_regions],axis=0)
    coarse_totals += cap_exact-dict(coarse_regions)["cap"]

    def wall_curvature(wall):
        return -16*math.pi*wall["radius"]*(wall["Rprime_plus"]-wall["Rprime_minus"])

    def wall_energy(wall):
        return 4*math.pi*wall["radius"]**2 * wall["energy"]

    def wall_stress(wall):
        gamma = wall["gamma"]
        return 4*math.pi*wall["radius"]**2 * wall["sigma"]*(gamma**2-3)/gamma

    bulk_volume, bulk_theta, bulk_theta2, bulk_shear, bulk_scalar = (
        float(totals[1]), float(totals[2]), float(totals[3]), float(totals[4]), float(totals[5]))
    host_volume = bulk_volume-float(cap_exact[1])-float(dict(regions)["outer_taper"][1]
                  + dict(regions)["outer_bridge_left"][1]+dict(regions)["outer_bridge_right"][1])
    exterior_volume = bulk_volume-float(cap_exact[1])-host_volume
    if min(float(cap_exact[1]), host_volume, exterior_volume) <= 0:
        _fail("a nested region has nonpositive volume")
    rho = 3/(KAPPA*ell**2)
    cap_energy = rho*float(cap_exact[1])
    host_energy = rho*host_volume
    inner_curvature = wall_curvature(inner)
    outer_curvature = wall_curvature(outer)
    curvature = bulk_scalar+inner_curvature+outer_curvature
    energy = cap_energy+host_energy+wall_energy(inner)+wall_energy(outer)
    stress = -3*(cap_energy+host_energy)+wall_stress(inner)+wall_stress(outer)
    theta_moment = bulk_theta
    jay = 2*bulk_theta2/3-2*bulk_shear
    diagnostic = scalar_average_sources(
        bulk_volume/ell**3, theta_moment/ell**2, bulk_theta2/ell, bulk_shear/ell,
        curvature/ell, energy/ell, stress/ell)
    _require("averaged initial Hamiltonian", diagnostic["hamiltonian_residual"])
    _require("averaged scalar-fitting identity", diagnostic["scalar_fitting_acceleration_residual"])
    coarse_host_volume = (coarse_totals[1]-cap_exact[1]
                         -sum(dict(coarse_regions)[name][1] for name in
                              ("outer_taper","outer_bridge_left","outer_bridge_right")))
    coarse_energy = rho*(cap_exact[1]+coarse_host_volume)+wall_energy(inner)+wall_energy(outer)
    coarse_stress = -3*rho*(cap_exact[1]+coarse_host_volume)+wall_stress(inner)+wall_stress(outer)
    coarse_diagnostic = scalar_average_sources(
        coarse_totals[1]/ell**3,coarse_totals[2]/ell**2,coarse_totals[3]/ell,
        coarse_totals[4]/ell,(coarse_totals[5]+inner_curvature+outer_curvature)/ell,
        coarse_energy/ell,coarse_stress/ell)

    lengths = {name: float(values[0]) for name, values in regions}
    lengths["cap"] = float(cap_exact[0])
    cumulative = 0.0
    milestone_l = {}
    for name in ("cap", "host_left", "host_bridge_left", "host_bridge_right", "host_transition",
                 "host_sheet_to_horizon", "host_sheet_to_wall", "outer_taper",
                 "outer_bridge_left", "outer_bridge_right"):
        cumulative += lengths[name]
        milestone_l[name] = cumulative

    def milestone(name, radius, proper, side):
        counter.add(1)
        return {"name": name, "R_over_ell": radius/ell, "l_over_ell": proper/ell, "side": side}

    centre = flat_cap_state(0.0, 1/ell)
    sheet_horizon = positive_sheet_state(host["rc"], HOST_MASS*ell, 1/ell)
    _require("horizon sheet", sheet_horizon["f"]+sheet_horizon["velocity"]**2-1)
    _require("horizon slope", sheet_horizon["Rprime"]-1)
    throat = bridge_chart(0.0, host["rb"], 1/ell, HOST_MASS*ell)
    _require("host throat slope", throat["Rprime"])
    outer_throat = bridge_chart(0.0, exterior["rb"], 0.0, OUTER_MASS*ell)
    cut_y = math.sqrt(OUTER_CUT*ell-exterior["rb"])
    cut = bridge_chart(cut_y, exterior["rb"], 0.0, OUTER_MASS*ell)
    reference = schwarzschild_bridge(math.acosh(math.sqrt(OUTER_CUT/(2*OUTER_MASS))), OUTER_MASS*ell)
    reference_radius = float(np.asarray(reference["R"]).reshape(-1)[0])
    reference_slope = float(np.asarray(reference["Rprime"]).reshape(-1)[0])
    counter.add(6)
    _require("retained cut matches Schwarzschild", (cut["radius"]-reference_radius)/ell)
    _require("retained cut slope", cut["Rprime"]-reference_slope)
    if time.process_time()-started > CPU_SECONDS:
        _fail("nested-slice CPU budget exhausted")
    record = {
        "schema": "nsc-atlas-nested-slice-v1",
        "assumptions": {
            "ell_positive": True, "c": 1, "G": 1, "k_ell": TENSION, "xi": "1/6",
            "phi": 0.0, "scalar_field_solves": 0, "evolution_steps": 0,
            "inner_hubble_ell": 1.0, "host_hubble_ell": 1.0, "exterior_hubble_ell": 0.0,
            "host_mass_over_ell": HOST_MASS, "exterior_mass_over_ell": OUTER_MASS,
            "hubble_is_supplied": True, "cap_recorded_from_exact_ball": True,
        },
        "horizons": {
            "host_rb_over_ell": host["rb"]/ell, "host_rc_over_ell": host["rc"]/ell,
            "host_root_f_residual": host["root_f_residual"],
            "host_cosmological_root_f_residual": host["cosmological_root_f_residual"],
            "host_mass_identity_residual_over_ell": host["mass_identity_residual"]/ell,
            "host_nariai_parameter": host["nariai_parameter"],
            "exterior_rb_over_ell": exterior["rb"]/ell,
            "exterior_root_is_two_mass": True,
        },
        "walls": {
            "inner": _wall_report(inner, ell),
            "outer": _wall_report(outer, ell),
        },
        "milestones": [
            milestone("centre", 0.0, 0.0, "cap"),
            milestone("inner_wall", inner["radius"], milestone_l["cap"], "both"),
            milestone("host_minimum", host["rb"], milestone_l["host_bridge_left"], "bridge"),
            milestone("host_transition_start", HOST_STATIC_RIGHT*ell,
                      milestone_l["host_bridge_right"], "bridge"),
            milestone("host_transition_end", HOST_EXPANDING_START*ell,
                      milestone_l["host_transition"], "sheet"),
            milestone("cosmological_horizon", host["rc"], milestone_l["host_sheet_to_horizon"], "sheet"),
            milestone("outer_wall", outer["radius"], milestone_l["host_sheet_to_wall"], "both"),
            milestone("outer_minimum", exterior["rb"], milestone_l["outer_bridge_left"], "bridge"),
            milestone("retained_cut", OUTER_CUT*ell, milestone_l["outer_bridge_right"], "bridge"),
        ],
        "endpoint_scalars": {
            "centre_theta_ell": centre["theta"]*ell,
            "centre_misner_sharp_over_ell": centre["misner_sharp"]/ell,
            "host_throat_Rprime": throat["Rprime"],
            "host_throat_dl_dy_over_sqrt_ell": throat["dl_dy"]/math.sqrt(ell),
            "horizon_Rprime": sheet_horizon["Rprime"],
            "horizon_velocity": sheet_horizon["velocity"],
            "outer_throat_Rprime": outer_throat["Rprime"],
            "cut_Rprime": cut["Rprime"],
            "cut_Kll_ell": cut["Kll"]*ell,
            "cut_Ktheta_ell": cut["Ktheta"]*ell,
        },
        "moments": {
            "volume_over_ell3": bulk_volume/ell**3,
            "I_theta_over_ell2": theta_moment/ell**2,
            "J_over_ell": jay/ell,
            "intrinsic_curvature_over_ell": curvature/ell,
            "bulk_scalar_over_ell": bulk_scalar/ell,
            "inner_wall_curvature_over_ell": inner_curvature/ell,
            "outer_wall_curvature_over_ell": outer_curvature/ell,
            "cap_energy_over_ell": cap_energy/ell,
            "host_energy_over_ell": host_energy/ell,
            "exterior_energy_over_ell": 0.0,
            "inner_wall_energy_over_ell": wall_energy(inner)/ell,
            "outer_wall_energy_over_ell": wall_energy(outer)/ell,
            "cap_stress_trace_over_ell": -3*cap_energy/ell,
            "host_stress_trace_over_ell": -3*host_energy/ell,
            "inner_wall_stress_trace_over_ell": wall_stress(inner)/ell,
            "outer_wall_stress_trace_over_ell": wall_stress(outer)/ell,
        },
        "average_diagnostic": {
            "clock": "Gaussian-normal initial slice with N=1 and Ndot=0",
            "role": "finite initial-slice diagnostic only",
            "supplied_hubble_is_not_a_derived_lambda": True,
            "physical_clock_and_history_open": True,
            "scalar_outputs_in_ell_units": {
                key: float(value) for key, value in diagnostic.items()},
        },
        "quadrature": {
            "orders": [32, 64],
            "difference_limit": "1e-7*max(1, abs(value64))",
            "differences_are_rigorous_bounds": False,
            "cap_replaced_by_exact_ball_after_agreement": True,
            "largest_regional_scaled_difference": max(regional_comparisons)[0],
            "largest_regional_scaled_difference_at": max(regional_comparisons)[1],
            "gauss64_minus_gauss32_in_ell_units": {
                "proper_distance_over_ell": float((totals[0]-coarse_totals[0])/ell),
                "volume_over_ell3": float((totals[1]-coarse_totals[1])/ell**3),
                "I_theta_over_ell2": float((totals[2]-coarse_totals[2])/ell**2),
                "J_over_ell": float((2*(totals[3]-coarse_totals[3])/3
                                      -2*(totals[4]-coarse_totals[4]))/ell),
                "intrinsic_curvature_over_ell": float((totals[5]-coarse_totals[5])/ell),
            },
            "scalar_diagnostic_64_minus_32": {
                key:float(value-coarse_diagnostic[key]) for key,value in diagnostic.items()},
        },
        "geometric_evaluations": counter.count,
        "cpu_seconds": time.process_time()-started,
        "analytic_end": {
            "retained_R_over_ell": OUTER_CUT,
            "continuation": "exact Schwarzschild right sheet, V=0, K=0, H=0, R to infinity",
            "finite_samples_do_not_infer_infinity": True,
        },
        "limitations": [
            "H is supplied regional initial data, not a derived cosmological Lambda or fractal dimension",
            "No formation history and no common three-dimensional sibling parent",
            "This is not a shared two-wall spacetime: static T_l=V/f is singular at horizons",
            "Kottler Kruskal, a right outgoing Eddington-Finkelstein chart, and matched worldtube clocks are not constructed",
            "No weak-shell well-posedness or nonsingular horizon-passage claim",
            "Quadrature differences are indicators, not rigorous bounds",
            "Beyond Rcut the right AE end is the analytic Schwarzschild sheet; enlarging that cut changes volume only",
        ],
    }
    if include_partition:
        # Expose already integrated scalar moments, never a new quadrature or
        # stored field grid. Keep the historical record's default shape.
        def moment_rows(integrated_regions):
            rows = {}
            for name, integrated in integrated_regions:
                part = cap_exact if name == 'cap' else integrated
                V, I, I2, shear, R3 = (float(part[i]) for i in range(1, 6))
                E = 0. if name.startswith('outer_') else rho*V
                rows[name] = [V/ell**3, I/ell**2, I2/ell,
                              shear/ell, R3/ell, E/ell, -3*E/ell]
            return rows
        surfaces = {}
        for name, wall in (('inner', inner), ('outer', outer)):
            surfaces[name] = dict(curvature=wall_curvature(wall)/ell,
                energy=wall_energy(wall)/ell, stress_trace=wall_stress(wall)/ell,
                volume_flux=4*math.pi*(wall['radius']/ell)**2*wall['v'])
        record['partition'] = dict(
            bulk_regions_in_ell_units=moment_rows(regions),
            coarse_bulk_regions_in_ell_units=moment_rows(coarse_regions),
            surfaces_in_ell_units=surfaces,
            moment_order=['volume', 'theta', 'theta_squared', 'shear_squared',
                          'curvature', 'energy', 'stress_trace'],
            clock='Existing initial common normal frame; no common future observer clock inferred.',
            additional_quadrature_values=0)
    return record


def _wall_report(wall, ell):
    return {
        "R_over_ell": wall["radius"]/ell,
        "beta_minus": wall["beta_minus"], "beta_plus": wall["beta_plus"],
        "Rdot": wall["Rdot"], "gamma": wall["gamma"], "v": wall["v"],
        "Rprime_minus": wall["Rprime_minus"], "Rprime_plus": wall["Rprime_plus"],
        "V_minus": wall["V_minus"], "V_plus": wall["V_plus"],
        "energy_ell": wall["energy"]*ell, "momentum_ell": wall["momentum"]*ell,
        "mean_trace_ell": wall["mean_trace"]*ell,
        "angular_jump_ell": wall["angular_jump"]*ell,
        "temporal_jump_ell": wall["temporal_jump"]*ell,
        "trace_jump_ell": wall["trace_jump"]*ell,
        "mass_map_residual_over_ell": wall["mass_map_residual"]/ell,
        "normal_force_residual_ell2": wall["normal_force_residual"]*ell**2,
        "lorentz_residual": wall["lorentz_residual"],
    }
