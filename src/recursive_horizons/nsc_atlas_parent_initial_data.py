"""Borrowed Brill--Lindquist initial data for a nonspherical shared parent.

This is one vacuum, time-symmetric spatial slice, not an evolved spacetime.
Punctures represent asymptotically flat ends, not child formation events or
traversable mouths. All lengths and geometric masses use G=c=1 units.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

SOURCE = {
    "url": "https://arxiv.org/html/2308.10906v1",
    "reading_depth": "Section III, equations (17)-(21), initial-data and MOTS paragraphs",
    "original": "Brill and Lindquist, Phys. Rev. 131, 471 (1963)",
    "original_url": "https://doi.org/10.1103/PhysRev.131.471",
    "original_reading_depth": "Bibliographic attribution; original article not read in this adapter review",
    "borrowed": "Conformal vacuum initial data and outer/puncture-end ADM masses",
}


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be finite and positive")
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return number


@dataclass(frozen=True)
class BrillLindquistData:
    """Positive bare masses at distinct punctures in an isotropic chart."""

    masses: tuple[float, ...]
    positions: tuple[tuple[float, float, float], ...]

    def __post_init__(self):
        masses = tuple(_positive(m, "bare mass") for m in self.masses)
        positions = np.asarray(self.positions, dtype=float)
        if not masses or positions.shape != (len(masses), 3):
            raise ValueError("one three-dimensional position is needed per mass")
        if not np.isfinite(positions).all():
            raise ValueError("puncture positions must be finite")
        if len({tuple(p) for p in positions}) != len(masses):
            raise ValueError("puncture positions must be distinct")
        object.__setattr__(self, "masses", masses)
        object.__setattr__(self, "positions", tuple(tuple(p) for p in positions.tolist()))


def two_punctures(mass_a=1.0, mass_b=1.5, separation=12.0):
    """Place two ends at z=+-d/2; d is coordinate, not proper separation."""
    distance = _positive(separation, "coordinate separation")
    return BrillLindquistData((mass_a, mass_b),
                             ((0.0, 0.0, -distance / 2), (0.0, 0.0, distance / 2)))


def _chart_sample(points, data, exclusion_radius):
    points = np.asarray(points, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 3 or not np.isfinite(points).all():
        raise ValueError("points must have finite coordinates and shape (..., 3)")
    if isinstance(exclusion_radius, (bool, np.bool_)):
        raise ValueError("exclusion radius must be finite and nonnegative")
    excluded = float(exclusion_radius)
    if not math.isfinite(excluded) or excluded < 0:
        raise ValueError("exclusion radius must be finite and nonnegative")
    delta = points[..., None, :] - np.asarray(data.positions)
    radius = np.linalg.norm(delta, axis=-1)
    if np.any(radius <= excluded):
        raise ValueError("sample lies on a puncture or within its excluded coordinate ball")
    return delta, radius


def conformal_factor(points, masses, positions, *, exclusion_radius=0.0):
    """Lightweight vectorized psi evaluation for a finite masked plot grid."""
    data = BrillLindquistData(tuple(masses), tuple(positions))
    _, radius = _chart_sample(points, data, exclusion_radius)
    return 1.0 + np.sum(np.asarray(data.masses) / (2 * radius), axis=-1)


def conformal_fields(points, data, *, exclusion_radius=0.0):
    """Return psi, its analytic flat derivatives, h, and R[h] off punctures.

    Points have shape (..., 3). The excluded closed balls use coordinate
    radii; they are sampling guards and are not horizons or physical cores.
    R[h]=-8 psi^-5 Delta_flat psi; the trace may retain roundoff.
    """
    delta, radius = _chart_sample(points, data, exclusion_radius)
    half_mass = np.asarray(data.masses) / 2
    psi = 1.0 + np.sum(half_mass / radius, axis=-1)
    gradient = -np.sum((half_mass / radius**3)[..., None] * delta, axis=-2)
    outer = delta[..., :, None] * delta[..., None, :]
    hessian = np.sum(half_mass[..., None, None] * (
        3 * outer / radius[..., None, None]**5
        - np.eye(3) / radius[..., None, None]**3), axis=-3)
    laplacian = np.trace(hessian, axis1=-2, axis2=-1)
    return {
        "psi": psi,
        "gradient": gradient,
        "hessian": hessian,
        "flat_laplacian": laplacian,
        "scalar_curvature": -8.0 * laplacian / psi**5,
        "metric": psi[..., None, None]**4 * np.eye(3),
    }


def adm_masses(data):
    """Asymptotic-end charges, not quasi-local/apparent-horizon masses."""
    masses = np.asarray(data.masses)
    positions = np.asarray(data.positions)
    distance = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=-1)
    np.fill_diagonal(distance, np.inf)
    puncture_ends = masses * (1.0 + np.sum(masses[None, :] / (2 * distance), axis=1))
    outer = float(masses.sum())
    return {
        "outer_end": outer,
        "puncture_ends": puncture_ends.tolist(),
        "outer_minus_sum_puncture_ends": float(outer - puncture_ends.sum()),
        "energy_convention": "Difference of asymptotic-end ADM charges; not a horizon-mass binding-energy measurement",
    }


def sphere_distortion(data, centre_index, radius, *, samples=17):
    """Exact angular data on a coordinate sphere in the two-puncture slice.

    cos(theta)=1 points toward the other puncture. The baseline replaces
    the companion term by its value at the selected puncture. The angular
    variation and quadrupole are retained without an acceptance threshold.
    They are chart-specific spatial diagnostics, not a measured 4D tide.
    """
    if len(data.masses) != 2 or isinstance(centre_index, bool) or centre_index not in (0, 1):
        raise ValueError("sphere distortion needs two punctures and centre index 0 or 1")
    if isinstance(samples, bool) or not isinstance(samples, (int, np.integer)) or not 3 <= samples <= 4097:
        raise ValueError("samples must be an integer between 3 and 4097")
    radius = _positive(radius, "coordinate sphere radius")
    other = 1 - centre_index
    distance = float(np.linalg.norm(np.asarray(data.positions[other]) - data.positions[centre_index]))
    if radius >= distance:
        raise ValueError("local sphere radius must be smaller than puncture separation")
    self_mass, other_mass = data.masses[centre_index], data.masses[other]
    cosines = np.linspace(-1.0, 1.0, samples)
    companion_distance = np.sqrt(distance**2 + radius**2 - 2 * distance * radius * cosines)
    psi = 1 + self_mass / (2 * radius) + other_mass / (2 * companion_distance)
    monopole = 1 + self_mass / (2 * radius) + other_mass / (2 * distance)
    first = other_mass * distance * radius / (2 * companion_distance**3)
    second = 3 * other_mass * distance**2 * radius**2 / (2 * companion_distance**5)
    lap_log = (1 - cosines**2) * (second / psi - (first / psi)**2) - 2 * cosines * first / psi
    # For q_AB=r^2 psi^4 Omega_AB, K=(1-2 Delta_S2 log psi)/(r^2 psi^4).
    curvature_times_r_squared = (1 - 2 * lap_log) / psi**4
    relative = psi / monopole - 1
    return {
        "centre_index": int(centre_index),
        "coordinate_radius": radius,
        "coordinate_radius_over_separation": radius / distance,
        "monopole_psi": monopole,
        "cos_theta": cosines.tolist(),
        "psi": psi.tolist(),
        "relative_psi_error": relative.tolist(),
        "relative_angular_metric_error": ((psi / monopole)**4 - 1).tolist(),
        "gaussian_curvature_times_coordinate_radius_squared": curvature_times_r_squared.tolist(),
        "maximum_absolute_relative_psi_error": float(np.max(np.abs(relative))),
        "quadrupole_coefficient_in_psi": other_mass * radius**2 / (2 * distance**3),
        "quadrupole_over_monopole": other_mass * radius**2 / (2 * distance**3 * monopole),
        "is_exactly_round": False,
        "horizon_status": "Not computed; this coordinate sphere is not identified as a horizon or junction",
        "diagnostic_scope": "Spatial chart angular distortion; no lapse, observer, or four-dimensional tidal tensor supplied",
    }


def parent_record():
    """Deterministic, JSON-safe borrowed parent data for the atlas producer."""
    data = two_punctures()
    points = [[0.0, 0.0, 0.0], [2.0, 0.0, -6.0], [0.0, 1.5, 6.0], [3.0, 4.0, 2.0]]
    fields = conformal_fields(points, data)
    return {
        "schema": "nsc-atlas-parent-initial-data-v1",
        "geometry": "Brill-Lindquist time-symmetric vacuum initial slice",
        "source": dict(SOURCE),
        "units": "G=c=1; lengths and geometric masses relative to one arbitrary length unit",
        "inputs": {"bare_masses": list(data.masses), "positions": [list(p) for p in data.positions], "coordinate_separation": 12.0},
        "dimensionless_parameters": {"mass_b_over_mass_a": 1.5, "separation_over_mass_a": 12.0},
        "equations": {"metric": "h_ij=psi^4 delta_ij", "extrinsic_curvature": "K_ij=0", "psi": "1+sum_i m_i/(2 |x-x_i|)", "constraint": "R[h]=-8 psi^-5 Delta_flat psi=0 away from punctures"},
        "spatial_topology": "R^3 minus two punctures; three asymptotically flat ends",
        "globally_spherical": False,
        "axisymmetric": True,
        "adm_masses": adm_masses(data),
        "sample_controls": {"points": points, "psi": fields["psi"].tolist(), "scalar_curvature": fields["scalar_curvature"].tolist()},
        "local_sphere_diagnostics": [sphere_distortion(data, i, data.masses[i]) for i in range(2)],
        "clock": {"proper_time_elapsed": None, "lapse": None, "status": "Initial slice only; K=0 does not specify lapse or elapsed time"},
        "horizons": {"computed": False, "common_horizon_excluded": False, "note": "No MOTS search performed; coordinate punctures and probe spheres do not locate horizons"},
        "causal_channels": [],
        "child_formation_established": False,
        "direct_spherical_child_attachment_established": False,
        "missing_conditions": ["Four-dimensional evolution and observer choice", "Actual horizon and candidate interface geometry", "Compatible induced metric, extrinsic curvature and physical source for any child attachment", "A physical transition identifying an expanding child instead of an existing asymptotic end"],
        "admission": {
            "kind": "Imported integration control, not a new scientific result",
            "borrowed": "Brill-Lindquist equations and ADM charges, arXiv:2308.10906v1 (17)-(21)",
            "missing": "The v1 diagram had no consistent common spatial parent for its two local holes",
            "decision": "Supply one common nonspherical spatial metric; display the finite angular correction lost by separate spherical approximations",
            "unlocks": "A named initial-data candidate for a later interface/domain test; no direct child gluing is admitted",
            "finite_domain": "One time-symmetric slice, two positive masses, finite coordinate separation, four off-puncture points and two 17-point local-sphere profiles",
            "budget": "At most one CPU minute for adapter evaluation and focused tests; no PDE evolution or horizon solver",
            "stop_or_redirect": "Reject puncture samples or intersecting local expansion domains; any failed vacuum/unit control blocks this adapter; exact spherical gluing requires a separate compatibility calculation",
            "independent_verification": "Finite-difference derivatives, asymptotic flux/inverted-end mass checks, sphere mean-value identity and unit rescaling in focused tests",
            "novelty_claim": False,
        },
    }
