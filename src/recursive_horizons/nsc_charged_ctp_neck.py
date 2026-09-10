"""Charged magnetic angular CTP source on the stored Bronnikov neck.

This is the new source owner.  It transports the scale-inherited Gaussian
occupation through the existing horizon frame and actual exterior reflection,
then performs the already adopted fourth-order adiabatic angular subtraction.
The compact Wilsonian local terms are appended from the locked scale branch;
the nonlocal compact remainder is reported separately.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import exp, log, pi, sqrt

import numpy as np
from scipy.special import expit

from .nsc_angular_stress import (
    _unitary_step,
    adiabatic_coefficients,
    curvature_squared_tensor,
    physical_source,
    profile_jets,
    static_cylinder_density,
)
from .nsc_horizon_source import conformal_stress, pg_source
from .nsc_unruh_state import ParentDirac, frequency_grid, horizon_frame


@dataclass(frozen=True)
class ChargedCTPNeckConfig:
    magnetic_flux: int
    omega: float
    horizon_rho: float
    surface_gravity: float
    angular_levels: int = 8
    points_per_frequency_interval: int = 12
    frequency_edges: tuple[float, ...] = (0.0, 0.25, 1.0, 4.0, 16.0)
    evolution_steps: int = 5000
    phase_cutoff: float = 4.0
    horizon_offset: float = 1e-10
    reflection_tolerance: float = 2e-10

    def __post_init__(self):
        if (isinstance(self.magnetic_flux, bool)
                or not isinstance(self.magnetic_flux, int)
                or self.magnetic_flux == 0):
            raise ValueError("nonzero integer magnetic flux required")
        if not all(np.isfinite((self.omega, self.horizon_rho,
                                self.surface_gravity, self.phase_cutoff,
                                self.horizon_offset, self.reflection_tolerance))):
            raise ValueError("finite source parameters required")
        if min(self.omega, self.surface_gravity, self.phase_cutoff,
               self.horizon_offset, self.reflection_tolerance) <= 0:
            raise ValueError("positive source scales required")
        if self.angular_levels < 1 or self.points_per_frequency_interval < 4:
            raise ValueError("resolved angular and frequency grids required")
        if self.evolution_steps < 100 or len(self.frequency_edges) < 2:
            raise ValueError("resolved evolution and frequency intervals required")
        if self.frequency_edges[0] != 0 or any(
            right <= left for left, right in zip(
                self.frequency_edges, self.frequency_edges[1:]
            )
        ):
            raise ValueError("ordered frequency edges beginning at zero required")


@dataclass(frozen=True)
class ChargedCTPNeckSource:
    config: dict
    massive_angular_source: dict
    lll_source: dict
    compact_local_source: dict
    totals: dict
    diagnostics: dict

    def to_dict(self):
        value = asdict(self)
        value["config"]["frequency_edges"] = list(value["config"]["frequency_edges"])
        return value


def _magnetic_levels(flux: int, maximum: int):
    n = np.arange(1, maximum + 1, dtype=float)
    masses = np.sqrt(n * (n + abs(flux)))
    degeneracies = 2.0 * (abs(flux) + 2.0 * n)
    return masses, degeneracies


def _reflection_data(background, frequencies, masses, phase_cutoff, tolerance):
    reflection = np.zeros((len(masses), len(frequencies)), complex)
    transmission = np.ones((len(masses), len(frequencies)), float)
    worst = 0.0
    solves = 0
    for row, mass in enumerate(masses):
        for column, frequency in enumerate(frequencies):
            if frequency > phase_cutoff:
                continue
            value = background.reflection(
                float(frequency), float(mass), tolerance=tolerance
            )
            reflection[row, column] = value["reflection"]
            transmission[row, column] = value["transmission"]
            worst = max(worst, value["current_defect"])
            solves += 1
    return reflection, transmission, worst, solves


def _massive_angular_source(config: ChargedCTPNeckConfig):
    background = ParentDirac(config.horizon_rho, config.surface_gravity)
    frequencies, weights = frequency_grid(
        config.frequency_edges, config.points_per_frequency_interval
    )
    masses_1d, degeneracies = _magnetic_levels(
        config.magnetic_flux, config.angular_levels
    )
    masses = masses_1d[:, None]
    reflection, transmission, current_defect, solves = _reflection_data(
        background,
        frequencies,
        masses_1d,
        config.phase_cutoff,
        config.reflection_tolerance,
    )
    outgoing_occupation = expit(
        -2.0 * pi * frequencies / config.surface_gravity
    )
    incoming_occupation = expit(
        -2.0 * pi * frequencies / (config.omega * config.surface_gravity)
    )
    covariance_trace = 1.0 + transmission * (
        incoming_occupation[None, :] - outgoing_occupation[None, :]
    )

    radius_horizon_squared = 1.0 + background.horizon_rho**2
    signed_offset = (
        -radius_horizon_squared * config.horizon_offset
        + background.horizon_rho * radius_horizon_squared
        * config.horizon_offset**2
    )
    tortoise = background.near_tortoise(signed_offset)
    distance = sqrt(2.0 * config.horizon_offset / config.surface_gravity)
    positive = np.empty_like(reflection)
    negative = np.empty_like(reflection)
    incoming_positive = np.empty_like(reflection)
    incoming_negative = np.empty_like(reflection)
    minimum_covariance_eigenvalue = 1.0
    maximum_covariance_eigenvalue = 0.0
    for row, mass in enumerate(masses_1d):
        for column, frequency in enumerate(frequencies):
            frame = horizon_frame(
                float(frequency), config.surface_gravity, float(mass),
                distance, tortoise, timelike=True,
            )
            occupied = frame @ np.array([
                sqrt(1.0 - outgoing_occupation[column]),
                -1j * sqrt(outgoing_occupation[column])
                * reflection[row, column],
            ])
            incoming = frame @ np.array([
                0.0,
                sqrt(incoming_occupation[column]
                     * transmission[row, column]),
            ])
            covariance = (
                np.outer(occupied, occupied.conj())
                + np.outer(incoming, incoming.conj())
            )
            eigenvalues = np.linalg.eigvalsh(covariance)
            minimum_covariance_eigenvalue = min(
                minimum_covariance_eigenvalue, float(eigenvalues[0])
            )
            maximum_covariance_eigenvalue = max(
                maximum_covariance_eigenvalue, float(eigenvalues[-1])
            )
            positive[row, column], negative[row, column] = occupied
            incoming_positive[row, column], incoming_negative[row, column] = incoming

    positive *= np.exp(-1j * pi / 4.0)
    negative *= np.exp(1j * pi / 4.0)
    incoming_positive *= np.exp(-1j * pi / 4.0)
    incoming_negative *= np.exp(1j * pi / 4.0)

    target = pi / 2.0
    begin = log(config.horizon_offset)
    end = log(background.horizon_q - target)
    step = (end - begin) / config.evolution_steps
    a = (3.0 - 2.0 * sqrt(3.0)) / 12.0
    b = (3.0 + 2.0 * sqrt(3.0)) / 12.0
    c1 = 0.5 - sqrt(3.0) / 6.0
    c2 = 0.5 + sqrt(3.0) / 6.0

    def generator(value):
        delta = exp(value)
        conformal = background.interior_W(delta)
        return (
            masses * delta / sqrt(conformal),
            -frequencies[None, :] * delta / conformal,
        )

    for index in range(config.evolution_steps):
        value = begin + index * step
        y1, z1 = generator(value + c1 * step)
        y2, z2 = generator(value + c2 * step)
        for hy, hz in (
            (b * y1 + a * y2, b * z1 + a * z2),
            (a * y1 + b * y2, a * z1 + b * z2),
        ):
            positive, negative = _unitary_step(
                positive, negative, hy, hz, step
            )
            incoming_positive, incoming_negative = _unitary_step(
                incoming_positive, incoming_negative, hy, hz, step
            )

    canonical_positive = positive * np.exp(1j * pi / 4.0)
    canonical_negative = negative * np.exp(-1j * pi / 4.0)
    canonical_incoming_positive = incoming_positive * np.exp(1j * pi / 4.0)
    canonical_incoming_negative = incoming_negative * np.exp(-1j * pi / 4.0)
    jets = profile_jets(target)
    conformal_scale = sqrt(jets[0])
    momentum = -frequencies[None, :] / conformal_scale
    energy = np.sqrt(masses**2 + momentum**2)
    cosine = masses / np.sqrt(2.0 * energy * (energy - momentum))
    sine = np.sqrt((energy - momentum) / (2.0 * energy))
    beta = cosine * canonical_positive + sine * canonical_negative
    alpha = -sine * canonical_positive + cosine * canonical_negative
    beta_in = (
        cosine * canonical_incoming_positive
        + sine * canonical_incoming_negative
    )
    alpha_in = (
        -sine * canonical_incoming_positive
        + cosine * canonical_incoming_negative
    )
    occupation = np.abs(beta)**2 + np.abs(beta_in)**2
    tangent = -2.0 * np.real(
        beta.conj() * alpha + beta_in.conj() * alpha_in
    )
    energy2, pressure2, energy4, pressure4 = adiabatic_coefficients(
        momentum / masses, 1.0, jets
    )
    factor = degeneracies / (4.0 * pi**2 * conformal_scale)
    excess = 2.0 * energy * occupation + energy * (1.0 - covariance_trace)
    pressure = (
        momentum**2 / energy**2 * excess
        + momentum * masses / energy * tangent
    )
    energy_remainders = factor * (
        (excess - energy2 / masses - energy4 / masses**3) @ weights
    )
    pressure_remainders = factor * (
        (pressure - pressure2 / masses - pressure4 / masses**3) @ weights
    )
    h_energy, h_parallel = curvature_squared_tensor(jets)
    harmonic = float(__import__("mpmath").euler) / 2.0
    bar_rho = (
        static_cylinder_density()
        + harmonic * h_energy / (480.0 * pi**2)
        + float(np.sum(energy_remainders))
    )
    bar_parallel = (
        -static_cylinder_density()
        - harmonic * h_parallel / (480.0 * pi**2)
        + float(np.sum(pressure_remainders))
    )
    w0, w1, w2, w3, w4 = jets
    scalar2 = -w2
    box_scalar = -w0 * w4 - w1 * w3
    trace = (
        -(scalar2 - 2.0)**2 / 60.0
        - 11.0 * scalar2 / 90.0
        - box_scalar / 30.0
    ) / (16.0 * pi**2)
    row = {
        "q": target,
        "rho_coordinate": 0.0,
        "sphere_radius": 1.0,
        "W_jets": jets.tolist(),
        "bar_rho": bar_rho,
        "bar_p_parallel": bar_parallel,
        "bar_trace": trace,
        "bar_p_sphere": (bar_rho - bar_parallel - trace) / 2.0,
    }
    source = physical_source(row)
    power = float(np.sum(
        degeneracies / pi * (
            transmission @ (
                weights * frequencies
                * (outgoing_occupation - incoming_occupation)
            )
        )
    ))
    norm_defect = float(np.max(np.abs(
        np.abs(positive)**2 + np.abs(negative)**2
        + np.abs(incoming_positive)**2 + np.abs(incoming_negative)**2
        - covariance_trace
    )))
    return {
        "source": source,
        "parent_Killing_power": power,
        "angular_masses": masses_1d.tolist(),
        "angular_degeneracies": degeneracies.tolist(),
        "energy_remainders": energy_remainders.tolist(),
        "pressure_remainders": pressure_remainders.tolist(),
        "maximum_scattering_current_defect": current_defect,
        "mode_norm_defect": norm_defect,
        "minimum_initial_covariance_eigenvalue": minimum_covariance_eigenvalue,
        "maximum_initial_covariance_eigenvalue": maximum_covariance_eigenvalue,
        "reflection_solves": solves,
    }


def charged_ctp_neck_source(
    config: ChargedCTPNeckConfig,
    *,
    compact_weyl_coefficient: float,
) -> ChargedCTPNeckSource:
    if not np.isfinite(compact_weyl_coefficient):
        raise ValueError("finite compact Weyl coefficient required")
    massive = _massive_angular_source(config)
    background = ParentDirac(config.horizon_rho, config.surface_gravity)
    outgoing = abs(config.magnetic_flux) * config.surface_gravity**2 / (48.0 * pi)
    incoming = config.omega**2 * outgoing
    metric_a = 1.0 - 3.0 * pi / 2.0
    metric_beta = sqrt(1.0 - metric_a)
    t_uu, t_uv, t_vv = conformal_stress(
        metric_a, 6.0, -3.0 * pi, outgoing, incoming,
        central_charge=abs(config.magnetic_flux),
    )
    lll_pg = pg_source(metric_a, metric_beta, 1.0, t_uu, t_uv, t_vv)
    area = 4.0 * pi
    coordinate = np.array([
        [lll_pg["T_tau_tau_2D"], lll_pg["T_tau_rho_2D"]],
        [lll_pg["T_tau_rho_2D"], lll_pg["T_rho_rho_2D"]],
    ]) / area
    pg_coframe = np.array([[1.0, 0.0], [metric_beta, 1.0]])
    inverse_pg = np.linalg.inv(pg_coframe)
    pg_frame = inverse_pg.T @ coordinate @ inverse_pg
    interior_f = -metric_a
    child_from_pg = np.array([
        [metric_beta / sqrt(interior_f), -1.0 / sqrt(interior_f)],
        [-1.0 / sqrt(interior_f), metric_beta / sqrt(interior_f)],
    ])
    inverse_child = np.linalg.inv(child_from_pg)
    child = inverse_child.T @ pg_frame @ inverse_child
    lll_child = {
        "rho": float(child[0, 0]),
        "T01": float(child[0, 1]),
        "p_parallel": float(child[1, 1]),
        # The reduced LLL state term has no independent sphere pressure;
        # gauge and gravitational local terms are accounted below.
        "p_sphere_state": 0.0,
    }
    lll_child["radial_null_plus"] = (
        lll_child["rho"] + lll_child["p_parallel"]
        + 2.0 * lll_child["T01"]
    )
    lll_child["radial_null_minus"] = (
        lll_child["rho"] + lll_child["p_parallel"]
        - 2.0 * lll_child["T01"]
    )
    if abs(
        lll_child["T01"]
        + lll_pg["parent_Killing_power"] / (4.0 * pi * interior_f)
    ) > 2e-13:
        raise ArithmeticError("LLL child-frame flux and Killing power differ")
    unit_weyl_null = -48.0 + 112.0 / 3.0 + 16.0 * pi
    compact_null = compact_weyl_coefficient * unit_weyl_null
    massive_flux = -massive["parent_Killing_power"] / (4.0 * pi * interior_f)
    massive_plus = massive["source"]["radial_null"] + 2.0 * massive_flux
    massive_minus = massive["source"]["radial_null"] - 2.0 * massive_flux
    compact_child = {
        "rho": -48.0 * compact_weyl_coefficient,
        "T01": 0.0,
        "p_parallel": (112.0 / 3.0 + 16.0 * pi) * compact_weyl_coefficient,
        "p_sphere": (-128.0 / 3.0 - 8.0 * pi) * compact_weyl_coefficient,
    }
    total_rho = (
        massive["source"]["rho"] + lll_child["rho"]
        + compact_child["rho"]
    )
    total_parallel = (
        massive["source"]["p_parallel"] + lll_child["p_parallel"]
        + compact_child["p_parallel"]
    )
    total_sphere = (
        massive["source"]["p_sphere"] + lll_child["p_sphere_state"]
        + compact_child["p_sphere"]
    )
    total_flux = massive_flux + lll_child["T01"]
    total_plus = total_rho + total_parallel + 2.0 * total_flux
    total_minus = total_rho + total_parallel - 2.0 * total_flux
    total_power = massive["parent_Killing_power"] + lll_pg["parent_Killing_power"]
    return ChargedCTPNeckSource(
        config={
            **asdict(config),
            "magnetic_matching_scale": sqrt(abs(config.magnetic_flux)),
        },
        massive_angular_source=massive,
        lll_source={
            "PG_coordinate": {key: float(value) for key, value in lll_pg.items()},
            "child_frame": lll_child,
        },
        compact_local_source={
            "C_Weyl": float(compact_weyl_coefficient),
            "unit_radial_null": unit_weyl_null,
            "radial_null": compact_null,
            "Killing_power": 0.0,
            "vacuum_radial_null": 0.0,
            "gauge_radial_null": 0.0,
            "child_frame": compact_child,
        },
        totals={
            "rho": float(total_rho),
            "T01": float(total_flux),
            "p_parallel": float(total_parallel),
            "p_sphere": float(total_sphere),
            "radial_null_plus": float(total_plus),
            "radial_null_minus": float(total_minus),
            "parent_Killing_power": float(total_power),
            "both_radial_null_components_negative": bool(
                total_plus < 0 and total_minus < 0
            ),
        },
        diagnostics={
            "horizon_q": float(background.horizon_q),
            "massive_child_T01": float(massive_flux),
            "massive_radial_null_plus": float(massive_plus),
            "massive_radial_null_minus": float(massive_minus),
            "retained_four_metric_projections_returned": True,
            "compact_nonlocal_remainder_included": False,
            "compact_excited_covariance_included": False,
            "full_owner_complete": False,
        },
    )
