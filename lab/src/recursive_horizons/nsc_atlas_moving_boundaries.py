"""Finite imported moving-wall control and the radiation flux obstruction.

The owning note records admission, conventions, source equations and limits.
This module does not attach either spherical candidate to the BL parent.
"""

from dataclasses import dataclass
import math

import numpy as np


SCHEMA = "nsc-atlas-moving-boundaries-v1"


@dataclass(frozen=True)
class VacuumWallParameters:
    """Geometric units: H,k have inverse-length units, mass has length units."""

    H: float = 1.0
    k: float = 2.0
    mass: float = 0.4

    def __post_init__(self):
        if any(not math.isfinite(x) or x <= 0 for x in (self.H, self.k, self.mass)):
            raise ValueError("H, k and Schwarzschild mass must be finite and positive")

    @property
    def sigma(self):
        return self.k / (4.0 * math.pi)


def perfect_fluid_flux(rho, pressure, speed, normal_sign=1):
    """T(U,n) in the fluid orthonormal frame, with n outward for sign +1.

    Speed is the shell's radial speed relative to that fluid, not dR/dtau.
    A perfect vacuum has rho+pressure=0, independent of that speed.
    """
    rho, pressure, speed = np.broadcast_arrays(
        np.asarray(rho, dtype=float), np.asarray(pressure, dtype=float),
        np.asarray(speed, dtype=float),
    )
    if normal_sign not in (-1, 1):
        raise ValueError("normal_sign must be +1 or -1")
    if not all(np.all(np.isfinite(x)) for x in (rho, pressure, speed)):
        raise ValueError("fluid data must be finite")
    if np.any(rho < 0) or np.any(np.abs(speed) >= 1):
        raise ValueError("rho must be nonnegative and the shell timelike")
    return normal_sign * (rho + pressure) * speed / (1.0 - speed**2)


def _wall_algebra(radius, parameters):
    R = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(R)) or np.any(R <= 0):
        raise ValueError("areal radii must be finite and positive")
    H, k, M = parameters.H, parameters.k, parameters.mass
    bm = M / (k * R**2) + (k*k - H*H) * R / (2*k)
    bp = M / (k * R**2) - (k*k + H*H) * R / (2*k)
    fm, fp = 1 - (H*R)**2, 1 - 2*M/R
    q = bm**2 - fm
    bm_prime = -2*M/(k*R**3) + (k*k-H*H)/(2*k)
    acceleration = bm * bm_prime + H*H*R
    return R, bm, bp, fm, fp, q, acceleration


def vacuum_wall_state(radius, parameters=VacuumWallParameters(), expanding=True):
    """Signed local shell data, restricted to a smooth chart segment.

    Turning points, zero beta and coordinate horizons require other adapters;
    rejecting these inputs does not assert a physical singularity there.
    """
    R, bm, bp, fm, fp, q, acceleration = _wall_algebra(radius, parameters)
    guard = 1e-12
    if np.any(q <= guard):
        raise ValueError("non-real motion or turning-point domain exit")
    if np.any(np.abs(bm) <= guard) or np.any(np.abs(bp) <= guard):
        raise ValueError("zero signed normal component: chart domain exit")
    if np.any(np.abs(fm) <= guard) or np.any(np.abs(fp) <= guard):
        raise ValueError("static-coordinate horizon: chart domain exit")
    if not isinstance(expanding, (bool, np.bool_)):
        raise ValueError("expanding must select one time-oriented branch")
    speed = np.sqrt(q) * (1 if expanding else -1)
    ktm = (acceleration - parameters.H**2*R) / bm
    ktp = (acceleration + parameters.mass/R**2) / bp
    sigma = (bm-bp)/(4*math.pi*R)
    pressure = (ktp-ktm + (bp-bm)/R)/(8*math.pi)
    signed_residual = (
        np.sign(bm)*np.sqrt(q+fm) - np.sign(bp)*np.sqrt(q+fp)
        - parameters.k*R
    )
    return {
        "R": R, "Rdot": speed, "Rddot": acceleration,
        "beta_minus": bm, "beta_plus": bp, "f_minus": fm, "f_plus": fp,
        "epsilon_minus": np.sign(bm), "epsilon_plus": np.sign(bp),
        "Ktau_minus": ktm, "Ktau_plus": ktp,
        "sigma_required": sigma, "pressure_required": pressure,
        "sigma_wall": np.full_like(R, parameters.sigma),
        "angular_residual": signed_residual,
        "pressure_residual": pressure + parameters.sigma,
    }


def proper_time_samples(radii, parameters=VacuumWallParameters(), order=32):
    """Proper duration since first radius for this finite expanding segment."""
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 2 or np.any(np.diff(radii) <= 0):
        raise ValueError("radii must be a strictly increasing one-dimensional grid")
    if not isinstance(order, int) or order < 4 or order > 64:
        raise ValueError("finite quadrature order must be between 4 and 64")
    vacuum_wall_state(radii, parameters)
    # q=a/R^4+b/R+c R^2-1 has one positive stationary radius and it is
    # the global minimum. Nodes alone can miss a narrow forbidden interval.
    H, k, M = parameters.H, parameters.k, parameters.mass
    a = M*M/(k*k)
    b = M*(k*k-H*H)/(k*k)
    c = (k*k+H*H)**2/(4*k*k)
    discriminant = math.sqrt(b*b+32*a*c)
    y = (b+discriminant)/(4*c) if b >= 0 else 8*a/(discriminant-b)
    minimum_radius = y**(1/3)
    if radii[0] <= minimum_radius <= radii[-1]:
        if _wall_algebra(minimum_radius, parameters)[5] <= 1e-12:
            raise ValueError("turning or forbidden interval between sample radii")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    durations = [0.0]
    for lower, upper in zip(radii[:-1], radii[1:]):
        radius = (upper+lower)/2 + (upper-lower)/2*nodes
        state = vacuum_wall_state(radius, parameters)
        # Sign/patch crossings between samples are outside this adapter.
        edges = vacuum_wall_state(np.array([lower, upper]), parameters)
        for key in ("beta_minus", "beta_plus", "f_minus", "f_plus"):
            if edges[key][0] * edges[key][1] <= 0:
                raise ValueError("patch or orientation crossing requires another adapter")
        duration = (upper-lower)/2 * np.dot(weights, 1/state["Rdot"])
        durations.append(durations[-1] + float(duration))
    return np.asarray(durations)


def _metric_acceleration_check(radius, parameters, side):
    """Direct Christoffel acceleration; differentiate tangent with a 5-point stencil."""
    R = np.asarray(radius, dtype=float)
    step = 2e-4 / parameters.H

    def tangent(r):
        state = vacuum_wall_state(r, parameters)
        return np.stack((state["beta_"+side]/state["f_"+side], state["Rdot"]), axis=-1)

    state = vacuum_wall_state(R, parameters)
    u = tangent(R)
    derivative = (tangent(R-2*step) - 8*tangent(R-step)
                  + 8*tangent(R+step) - tangent(R+2*step))/(12*step)
    f = state["f_"+side]
    fp = -2*parameters.H**2*R if side == "minus" else 2*parameters.mass/R**2
    ar = derivative[..., 1]*u[..., 1] + f*fp*u[..., 0]**2/2 - fp*u[..., 1]**2/(2*f)
    at = derivative[..., 0]*u[..., 1] + fp/f*u[..., 0]*u[..., 1]
    ktau = -u[..., 1]*at + u[..., 0]*ar
    norm = -f*u[..., 0]**2 + u[..., 1]**2/f
    return ktau, norm


def _maxabs(value):
    return float(np.max(np.abs(value)))


def moving_boundary_record(length_scale=1.0):
    """One fixed dimensionless candidate; no fitting or parameter sweep."""
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length_scale must be finite and positive")
    parameters = VacuumWallParameters(H=1/ell, k=2/ell, mass=0.4*ell)
    radii = ell*np.linspace(1.2, 1.8, 33)
    state = vacuum_wall_state(radii, parameters)
    tau32 = proper_time_samples(radii, parameters, 32)
    tau64 = proper_time_samples(radii, parameters, 64)
    ktm, um_norm = _metric_acceleration_check(radii, parameters, "minus")
    ktp, up_norm = _metric_acceleration_check(radii, parameters, "plus")
    pressure_metric = (ktp-ktm+(state["beta_plus"]-state["beta_minus"])/radii)/(8*math.pi)
    velocities = np.linspace(-0.8, 0.8, 9)
    rho = 3/(8*math.pi*ell**2)
    pressure = rho/3
    gamma = 1/np.sqrt(1-velocities**2)
    U = np.stack((gamma, gamma*velocities), axis=-1)
    n = np.stack((gamma*velocities, gamma), axis=-1)
    tensor = np.diag([rho, pressure])
    contracted = np.einsum("ni,ij,nj->n", U, tensor, n)
    flux = perfect_fluid_flux(rho, pressure, velocities)
    H, k, M = parameters.H, parameters.k, parameters.mass
    expanded_q = M*M/(k*k*radii**4) + M*(k*k-H*H)/(k*k*radii) \
        + (k*k+H*H)**2*radii**2/(4*k*k)-1
    samples = {
        "R_over_ell": (radii/ell).tolist(), "tau_over_ell": (tau32/ell).tolist(),
        "Rdot": state["Rdot"].tolist(), "ell_Rddot": (ell*state["Rddot"]).tolist(),
        "ell_sigma_required": (ell*state["sigma_required"]).tolist(),
        "ell_pressure_required": (ell*state["pressure_required"]).tolist(),
        "ell_sigma_wall": (ell*state["sigma_wall"]).tolist(),
        "ell_pressure_residual": (ell*state["pressure_residual"]).tolist(),
    }
    for key in ("beta_minus", "beta_plus", "f_minus", "f_plus", "epsilon_minus",
                "epsilon_plus", "angular_residual"):
        samples[key] = state[key].tolist()
    return {
        "schema": SCHEMA,
        "claim_type": "application_of_imported_conservation_and_vacuum_wall_dynamics",
        "admission": {
            "owner": "lab/docs/nsc-atlas-moving-boundaries.md",
            "borrowed": ["gr-qc/0507063v2 Eq8-10", "0711.2198v1 Eq2 and Sec3",
                         "gr-qc/0508093v2 Eq12-13"],
            "missing": "whether noncomoving motion alone repairs the v2 radiation wall",
            "decision": "motion alone cannot repair it; a separate vacuum wall moves",
            "unlocks": "source and nonspherical parent-boundary compatibility",
            "cpu_seconds_max": 60, "quadrature_evaluations": 32*(32+64),
            "stop": "nonzero radiation flux or failed unsquared/metric junction check",
            "priority_claim": False,
        },
        "radiation_flux_screen": {
            "status": "MOTION_ALONE_INCOMPATIBLE_WITH_FIXED_SOURCE",
            "assumptions": ["timelike spherical shell", "positive-density perfect radiation",
                            "constant pure tension", "vacuum Schwarzschild exterior",
                            "no additional surface fields, non-fluid layer or exterior flux"],
            "convention": "U=gamma(u+v e), n=gamma(v u+e), [F]=F_plus-F_minus",
            "equation": "F_minus=(rho+p)*v/(1-v^2); F_plus=0; [F]=0 implies v=0",
            "normal_reversal_changes_conclusion": False,
            "necessary_not_sufficient": True,
            "speed_relative_to_fluid": velocities.tolist(),
            "F_minus_over_rho_plus_p": (flux/(rho+pressure)).tolist(),
            "jump_over_rho_plus_p": (-flux/(rho+pressure)).tolist(),
            "existing_comoving_obstruction": "nsc-cosmological-atlas-v2/junction.json",
            "redirects": ["specified non-vacuum exterior flux", "specified non-pure-wall stress",
                          "different bulk such as de Sitter vacuum"],
        },
        "vacuum_wall_control": {
            "status": "IMPORTED_LOCAL_MOVING_SOLUTION",
            "theory": {"gravity": "four-dimensional Einstein gravity", "H_ell": 1.0,
                       "k_ell": 2.0, "M_over_ell": 0.4, "k_definition": "4*pi*sigma",
                       "wall_law": "positive constant sigma, pressure=-sigma",
                       "fixed_within_test": True, "interior": "de Sitter vacuum",
                       "exterior": "Schwarzschild vacuum"},
            "domain": {"R_over_ell": [1.2, 1.8], "sample_count": 33,
                       "expanding": True, "proper_duration_over_ell": float(tau32[-1]/ell),
                       "excluded": ["turning points", "zero beta", "static-coordinate horizons"],
                       "exclusions_are_physical_singularities": False},
            "normal": "single normal from minus bulk toward plus bulk; beta=n^a partial_a R",
            "orientation": {"epsilon_minus": 1, "epsilon_plus": -1,
                            "meaning": "normal toward decreasing Schwarzschild areal radius"},
            "clock": {"recorded": "proper shell duration from R/ell=1.2",
                      "de_sitter_static_t": "spacelike on this segment",
                      "schwarzschild_patch_t": "decreases along chosen shell tangent",
                      "cross_universe_age_established": False},
            "samples": samples, "embedded_in_shared_parent": False,
            "child_formation_established": False, "stability_established": False,
            "sibling_signalling_established": False,
        },
        "independent_checks": {
            "lorentz_contraction_max_ell2_error": ell**2*_maxabs(contracted-flux),
            "vacuum_flux_max_ell2": ell**2*_maxabs(perfect_fluid_flux(rho, -rho, velocities)),
            "unsquared_angular_max_error": _maxabs(state["angular_residual"]),
            "pressure_max_ell_error": ell*_maxabs(state["pressure_residual"]),
            "metric_pressure_max_ell_error": ell*_maxabs(pressure_metric+parameters.sigma),
            "metric_Ktau_max_ell_error": ell*max(_maxabs(ktm-state["Ktau_minus"]),
                                                _maxabs(ktp-state["Ktau_plus"])),
            "tangent_normalization_max_error": max(_maxabs(um_norm+1), _maxabs(up_norm+1)),
            "expanded_potential_max_error": _maxabs(expanded_q-state["Rdot"]**2),
            "proper_time_32_64_max_error_over_ell": _maxabs(tau32-tau64)/ell,
            "unsigned_angular_min_error": float(np.min(np.abs(
                np.abs(state["beta_minus"])-np.abs(state["beta_plus"])-k*radii))),
            "quadrature_comparison_is_rigorous_error_bound": False,
        },
        "source_reading": [
            {"id": "gr-qc/0507063v2", "url": "https://arxiv.org/html/gr-qc/0507063v2",
             "depth": "Section 2 Eq3,4,8-10; Section 3.2 Eq24-28",
             "use": "general tensor conservation and Israel stresses; not static-fluid Eq26 for FLRW"},
            {"id": "0711.2198v1", "url": "https://arxiv.org/html/0711.2198v1",
             "depth": "Sections 2-3 including Eq2 and epsilon interpretation",
             "use": "fixed-tension de Sitter/Schwarzschild dynamics and formation limits"},
            {"id": "gr-qc/0508093v2", "url": "https://arxiv.org/html/gr-qc/0508093v2",
             "depth": "Sections II.1-II.3; Eq12-13 and chart-sign discussion",
             "use": "independent signed equation; its coordinate-time conventions differ by patch"},
        ],
        "limits": ["no new-to-world theorem or priority claim", "no formation history or regular origin",
                   "no perturbative stability test", "not a junction to the nonspherical BL parent",
                   "no inherited-constant, antimatter, dark-sector or observer-age prediction"],
    }
