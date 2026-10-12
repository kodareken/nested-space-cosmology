"""Coupled response of the existing radial RN step.

The adapter differentiates the owners it imports. It does not replace the
discrete product rule, add a pump or a reset, or advance a trajectory.
The production march step remains the cap 0.001 r_m. Half of that step is
0.0005 r_m. Half the full CFL is a larger number and is not this step.

Width and target-energy initials are algebraic on the truncated window.
The central-frequency preparation derivative is a centered difference of
``prepare_radial_packet``, including its phase alignment, at two step sizes.
"""
from __future__ import annotations

import importlib
import math
import os

for _thread_name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_name, "1")

import numpy as np

from . import nsc_rn_observables as observables
from . import nsc_rn_pg as pg
from . import nsc_rn_source as source


API_VERSION = "nsc-rn-response-v1"
BASELINE_MASS_OVER_RM = 1.04
BASELINE_ENERGY_TARGET = 1.0e-3
BASELINE_RANK = 9
BASELINE_WIDTH_OVER_RM = source.PACKET_WIDTH_OVER_RM
BASELINE_CENTER_OVER_RM = source.PACKET_CENTER_OVER_RM
BASELINE_OUTER_OVER_RM = 32.0
STEP_CAP_OVER_RM = 0.001
PRODUCTION_HALF_STEP_OVER_RM = 0.5 * STEP_CAP_OVER_RM
RECONSTRUCTION_WEIGHT = 1.0
SHELL_KAPPA = source.SHELL_KAPPA
SHELL_MULTIPLICITY = source.SHELL_MULTIPLICITY
# Reviewer figure for half the full CFL. It is not the march step.
REJECTED_HALF_CFL_OVER_RM = 0.008742
FREQUENCY_PREPARE_POINT_CAP = 128
OUTER_STATIC_PROTOCOL = "outer_static"
LEGACY_EXCISION_LAPSE = "legacy_fixed_excision_lapse"
TRANSPORTED_NORMAL_PROTOCOL = "transported_normal"
# Field RK4 stages start from the original state. Fractions multiply dt
# and the previous stage rate. A separate observer samples each metric.
FIELD_RK4_STAGE_FRACTIONS = (0.0, 0.5, 0.5, 1.0)
HELD_FORECAST = {
    "width_fraction": 0.05,
    "frequency_fraction": 0.05,
    "fixed_initial_incoming_killing_energy": True,
    "status": "not_run",
    "missing_piece": (
        "integrated capture, content, and delay under the finite +5% width "
        "and +5% frequency controls, compared with this tangent"
    ),
}


def baseline():
    """Declared numerical family point. No scattering solve and no march."""
    shell, insertions = observables.shell_factor(1, SHELL_MULTIPLICITY)
    return {
        "api_version": API_VERSION,
        "mass_over_rm": BASELINE_MASS_OVER_RM,
        "epsilon": BASELINE_ENERGY_TARGET,
        "kappa": SHELL_KAPPA,
        "multiplicity": int(shell),
        "shell_insertions": int(insertions),
        "reconstruction_weight": RECONSTRUCTION_WEIGHT,
        "frequency_samples": BASELINE_RANK,
        "phi_rank": 1,
        "width_over_rm": BASELINE_WIDTH_OVER_RM,
        "center_over_rm": BASELINE_CENTER_OVER_RM,
        "outer_over_rm": BASELINE_OUTER_OVER_RM,
        "step_cap_over_rm": STEP_CAP_OVER_RM,
        "production_half_step_over_rm": PRODUCTION_HALF_STEP_OVER_RM,
        "rejected_half_cfl_over_rm": REJECTED_HALF_CFL_OVER_RM,
        "killing_frequency_metadata_is_spectrum": False,
        "held_forecast": dict(HELD_FORECAST),
        "new_dynamical_relation": False,
    }


def production_step(cfl_dt, r_m):
    """March step from the parent driver. The cap wins on the baseline pair."""
    driver = importlib.import_module("derive_nsc_rn_parent")
    choice = dict(driver.matched_step(cfl_dt, r_m))
    choice["rejected_half_cfl_over_rm"] = REJECTED_HALF_CFL_OVER_RM
    choice["production_half_step"] = 0.5 * float(choice["step_cap"])
    choice["half_cfl_is_the_march_step"] = False
    return choice


def _column(phi):
    array = np.asarray(phi, dtype=complex)
    if array.ndim == 2:
        array = array[:, :, None]
    if array.ndim != 3 or array.shape[0] != 2 or array.shape[2] != 1:
        raise ValueError("response tangent expects the rank-one column (2, points, 1)")
    return array


def weighted_pairing(left, right, weights):
    """Owner Gram off-diagonal ⟨left, right⟩. Multiplicity is not included."""
    stacked = np.concatenate([_column(left), _column(right)], axis=2)
    return source.weighted_gram(stacked, np.asarray(weights, dtype=float))[0, 1]


def real_overlap(left, right, weights):
    return float(np.real(weighted_pairing(left, right, weights)))


def filling_derivative(packet_energy, mass, epsilon, newton_constant, *,
                       d_epsilon, d_packet_energy, multiplicity=SHELL_MULTIPLICITY):
    """ν = ε M / (G_N 4 E), so dν/ν = dε/ε − dE/E. The 4 is already inside ν."""
    energy = float(packet_energy)
    mass = float(mass)
    epsilon = float(epsilon)
    newton = float(newton_constant)
    multiplicity = int(multiplicity)
    if min(energy, mass, epsilon, newton) <= 0.0 or multiplicity != SHELL_MULTIPLICITY:
        raise ValueError("filling derivative needs the positive κ=1 shell data")
    filling = epsilon * mass / (newton * multiplicity * energy)
    relative = float(d_epsilon) / epsilon - float(d_packet_energy) / energy
    return {
        "filling": float(filling),
        "d_filling": float(filling * relative),
        "relative": relative,
        "multiplicity": multiplicity,
        "multiplicity_applied_once": True,
        "formula": "nu = epsilon * M / (G_N * 4 * E_packet)",
    }


def _zero_clocks():
    return {"t": 0.0, "killing": 0.0, "normal": 0.0}


def _zero_stocks():
    return {"excision_flux": 0.0, "outer_flux": 0.0, "sat_debit": 0.0}


def _zero_probability():
    return {"remaining": 0.0, "excision_outflow": 0.0, "outer_outflow": 0.0}


def state_direction(dphi, d_occupations, *, d_inner_mass=0.0, d_cached_covariance=None):
    """Initial small direction. Clock and stock tangents start at zero."""
    phi = _column(dphi)
    occupations = np.asarray(d_occupations, dtype=float).reshape(-1)
    if occupations.shape != (1,):
        raise ValueError("rank-one filling tangent must have one occupation")
    return {
        "dphi": phi,
        "d_occupations": occupations,
        "d_inner_mass": float(d_inner_mass),
        "d_clocks": _zero_clocks(),
        "d_stocks": _zero_stocks(),
        "d_probability_stocks": _zero_probability(),
        "d_cached_covariance": None if d_cached_covariance is None else np.array(
            d_cached_covariance, dtype=complex, copy=True
        ),
        "d_electric_current": 0.0,
        "d_magnetic_flux": 0.0,
        "d_kappa": 0.0,
        "d_killing_frequency_metadata": 0.0,
        "metadata_moves_the_spectrum": False,
        "dN_out_constraint": 0.0,
    }


def energy_initial_tangent(prepared, d_epsilon=1.0):
    """Target ε at fixed width and frequency. The column does not move."""
    phi = _column(prepared["phi"])
    energy = float(prepared["packet_mode_energy"])
    epsilon = float(prepared["energy_target"])
    if abs(float(prepared["actual_energy_target"]) - epsilon) > 1e-12:
        raise ValueError("actual energy target is not the requested ε")
    derived = filling_derivative(
        energy, prepared["background"]["mass"], epsilon, prepared["G_N"],
        d_epsilon=d_epsilon, d_packet_energy=0.0,
    )
    if abs(derived["filling"] - float(prepared["packet_filling"])) > 1e-12 * max(1.0, derived["filling"]):
        raise ValueError("prepared filling does not match ν = ε M / (G_N 4 E)")
    zeros = np.zeros_like(phi)
    direction = state_direction(zeros, np.array([derived["d_filling"]]))
    car = np.asarray(prepared["car"]["physical"], dtype=complex)
    scale = derived["relative"]
    return {
        "kind": "energy",
        "direction": direction,
        "d_packet_energy": 0.0,
        "d_actual_energy_target": float(d_epsilon),
        "d_sample_energies": np.zeros_like(np.asarray(prepared["spectrum"]["requested"], dtype=float)),
        "sample_nodes_move": False,
        "modes_move": False,
        "d_car": scale * car,
        "multiplicity_in_car": False,
        "real_overlap": 0.0,
        "formula": "dphi = 0, dnu/nu = depsilon/epsilon",
    }


def _truncated_window(omega, spectral_width, rank):
    energies, truncated, span = source._positive_frequency_window(omega, spectral_width, rank)
    if not truncated:
        raise ValueError("algebraic width tangent requires the truncated positive window")
    return energies, float(span)


def quadratic_energy_tangent(mode, dmode, lapse, shift, radius, weights, derivative, *, kappa=SHELL_KAPPA):
    """Derivative of mode_quadratic_form through the owner bilinears."""
    mode = np.asarray(mode, dtype=complex)
    dmode = np.asarray(dmode, dtype=complex)
    if mode.shape != dmode.shape or mode.ndim != 2 or mode.shape[0] != 2:
        raise ValueError("quadratic tangent expects a characteristic pair")
    radial = np.ones_like(np.asarray(radius, dtype=float))
    base = source.mode_quadratic_form(
        mode, lapse, shift, radial, radius, weights, derivative, kappa=kappa,
    )
    speed_plus, speed_minus = source._speeds_from_fields(
        np.asarray(lapse, dtype=float), np.asarray(shift, dtype=float), radial,
    )
    d_plus = (
        source._hermitian_flux(mode[0], dmode[0], derivative)
        + source._hermitian_flux(dmode[0], mode[0], derivative)
    )
    d_minus = (
        source._hermitian_flux(mode[1], dmode[1], derivative)
        + source._hermitian_flux(dmode[1], mode[1], derivative)
    )
    d_angular = (
        source._angular_bilinear(mode, dmode) + source._angular_bilinear(dmode, mode)
    )
    density = (
        speed_plus * d_plus + speed_minus * d_minus
        + (float(kappa) * np.asarray(lapse, dtype=float) / np.asarray(radius, dtype=float)) * d_angular
    )
    return {
        "energy": float(base["energy"]),
        "d_energy": float(np.real(np.sum(np.asarray(weights, dtype=float) * density))),
    }


def width_initial_tangent(prepared, d_width=1.0):
    """Truncated-window width derivative at fixed ε and fixed frequency nodes."""
    spectral = np.asarray(prepared["spectral_phi"], dtype=complex)
    if spectral.ndim != 3 or spectral.shape[0] != 2:
        raise ValueError("spectral_phi must have shape (2, points, rank)")
    rank = int(spectral.shape[2])
    omega = float(prepared["radial_moments"]["requested_omega"])
    width = float(prepared["radial_moments"]["requested_width"])
    speed = float(prepared["spectrum"]["ingoing_group_velocity"])
    if width <= 0.0 or omega <= 0.0 or speed >= 0.0:
        raise ValueError("width tangent needs a positive width and an ingoing center")
    spectral_width = abs(speed) / width
    if abs(spectral_width - float(prepared["spectrum"]["spectral_width"])) > 1e-10 * spectral_width:
        raise ValueError("stored spectral width is not |a-| / requested width")
    energies, span = _truncated_window(omega, spectral_width, rank)
    stored = np.asarray(prepared["spectrum"]["requested"], dtype=float)
    if stored.shape != energies.shape or np.max(np.abs(stored - energies)) > 1e-10 * max(1.0, omega):
        raise ValueError("stored sample energies are not the truncated window")
    amplitudes = source._gaussian_amplitudes(energies, omega, spectral_width)
    weights = np.asarray(prepared["background"]["weights"], dtype=float)
    superposition = spectral @ amplitudes
    norm = math.sqrt(float(np.real(np.sum(weights * np.sum(np.abs(superposition) ** 2, axis=0)))))
    packet = superposition / norm
    stored_phi = _column(prepared["phi"])[:, :, 0]
    if np.max(np.abs(packet - stored_phi)) > 1e-8:
        raise ValueError("stored column is not the normalized truncated Gaussian packet")
    xi = (energies - omega) / spectral_width
    d_spectral_width = -spectral_width / width * float(d_width)
    d_xi = -xi / spectral_width * d_spectral_width
    d_amplitudes = amplitudes * (-0.5 * xi * d_xi)
    d_superposition = spectral @ d_amplitudes
    d_norm = real_overlap(superposition[:, :, None], d_superposition[:, :, None], weights) / norm
    eta = d_superposition / norm
    d_packet = eta - packet * (d_norm / norm)
    overlap = real_overlap(packet[:, :, None], d_packet[:, :, None], weights)
    background = prepared["background"]
    energy = quadratic_energy_tangent(
        packet, d_packet, background["lapse"], background["shift"], background["radius"],
        weights, background["derivative"],
    )
    if abs(energy["energy"] - float(prepared["packet_mode_energy"])) > 1e-8 * max(1.0, energy["energy"]):
        raise ValueError("packet energy does not match mode_quadratic_form")
    derived = filling_derivative(
        energy["energy"], background["mass"], prepared["energy_target"], prepared["G_N"],
        d_epsilon=0.0, d_packet_energy=energy["d_energy"],
    )
    alpha = amplitudes / norm
    d_alpha = d_amplitudes / norm - amplitudes * d_norm / norm**2
    coefficient = float(prepared["packet_filling"]) * np.outer(alpha, np.conj(alpha))
    d_coefficient = derived["d_filling"] * np.outer(alpha, np.conj(alpha)) + float(prepared["packet_filling"]) * (
        np.outer(d_alpha, np.conj(alpha)) + np.outer(alpha, np.conj(d_alpha))
    )
    gram = source.weighted_gram(spectral, weights)
    root = source._gram_sqrt(gram)
    d_car = root @ d_coefficient @ root
    d_car = 0.5 * (d_car + d_car.conj().T)
    direction = state_direction(d_packet[:, :, None], np.array([derived["d_filling"]]), d_cached_covariance=d_car)
    return {
        "kind": "width",
        "direction": direction,
        "d_packet_energy": energy["d_energy"],
        "d_actual_energy_target": 0.0,
        "d_sample_energies": np.zeros(rank),
        "sample_nodes_move": False,
        "modes_move": False,
        "d_span": span / width * float(d_width),
        "d_car": d_car,
        "d_coefficient": d_coefficient,
        "coefficient": coefficient,
        "real_overlap": overlap,
        "multiplicity_in_car": False,
        "multiplicity_in_gram": False,
        "formula": "dnu/nu = -dE/E at fixed epsilon; Re<phi, dphi> = 0",
    }


def frequency_preparation_derivative(step, **prepare_kwargs):
    """Centered prepare_radial_packet derivative at ``step`` and ``step/2``.

    Phase alignment is the alignment already inside the owner. Grids above
    the point cap are refused so this call is not a production campaign.
    """
    step = float(step)
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("frequency step must be positive")
    if prepare_kwargs.get("radii") is None or prepare_kwargs.get("omega") is None:
        raise ValueError("frequency derivative needs an explicit omega and caller radii")
    if int(np.asarray(prepare_kwargs["radii"]).shape[0]) > FREQUENCY_PREPARE_POINT_CAP:
        raise ValueError("frequency preparation derivative refuses a production-sized grid")
    omega = float(prepare_kwargs["omega"])

    def prepare_at(value):
        call = dict(prepare_kwargs)
        call["omega"] = float(value)
        return source.prepare_radial_packet(**call)

    samples = {
        "plus": prepare_at(omega + step),
        "minus": prepare_at(omega - step),
        "plus_half": prepare_at(omega + 0.5 * step),
        "minus_half": prepare_at(omega - 0.5 * step),
    }

    def center_ingoing(prepared):
        grid = np.asarray(prepared["background"]["grid"], dtype=float)
        center = float(prepared["radial_moments"]["requested_center"])
        index = int(np.argmin(np.abs(grid - center)))
        value = complex(prepared["phi"][1, index, 0])
        return index, value

    phases = {}
    for name, prepared in samples.items():
        index, value = center_ingoing(prepared)
        phases[name] = {
            "index": index,
            "value": value,
            "real_positive": bool(value.real > 0.0 and abs(value.imag) <= 1e-8 * max(1.0, abs(value))),
        }
    dphi = (samples["plus"]["phi"] - samples["minus"]["phi"]) / (2.0 * step)
    dphi_half = (samples["plus_half"]["phi"] - samples["minus_half"]["phi"]) / step
    d_filling = (samples["plus"]["packet_filling"] - samples["minus"]["packet_filling"]) / (2.0 * step)
    d_filling_half = (
        samples["plus_half"]["packet_filling"] - samples["minus_half"]["packet_filling"]
    ) / step
    d_energy = (samples["plus"]["packet_mode_energy"] - samples["minus"]["packet_mode_energy"]) / (2.0 * step)
    energies_plus = np.asarray(samples["plus"]["spectrum"]["requested"], dtype=float)
    energies_minus = np.asarray(samples["minus"]["spectrum"]["requested"], dtype=float)
    d_energies = (energies_plus - energies_minus) / (2.0 * step)
    energies = 0.5 * (energies_plus + energies_minus)
    filling_base = 0.5 * (samples["plus"]["packet_filling"] + samples["minus"]["packet_filling"])
    energy_base = 0.5 * (samples["plus"]["packet_mode_energy"] + samples["minus"]["packet_mode_energy"])
    direction = state_direction(dphi, np.array([d_filling]))
    return {
        "kind": "frequency",
        "direction": direction,
        "dphi_half_step": dphi_half,
        "d_filling": float(d_filling),
        "d_filling_half_step": float(d_filling_half),
        "d_packet_energy": float(d_energy),
        "d_sample_energies": d_energies,
        "sample_energy_over_frequency": energies / omega,
        "filling_base": float(filling_base),
        "energy_base": float(energy_base),
        "relative_filling_plus_relative_energy": float(
            d_filling / filling_base + d_energy / energy_base
        ),
        "phi_step_discrepancy": float(np.max(np.abs(dphi - dphi_half))),
        "phi_signal": float(np.max(np.abs(dphi_half))),
        "phase_alignment": phases,
        "phase_aligned": all(item["real_positive"] for item in phases.values()),
        "modes_move": True,
        "sample_nodes_move": True,
        "variational_ode": False,
        "scattering_campaign": False,
        "step": step,
        "half_step": 0.5 * step,
    }


def stress_tangent(grid, phi, occupations, dphi, d_occupations, *, kappa=SHELL_KAPPA):
    """Unit-lapse shell stress. Multiplicity 4 multiplies the bilinears once."""
    multiplicity = source.shell_multiplicity(kappa)
    column = _column(phi)[:, :, 0]
    direction = _column(dphi)[:, :, 0]
    nu = float(np.asarray(occupations, dtype=float).reshape(-1)[0])
    dnu = float(np.asarray(d_occupations, dtype=float).reshape(-1)[0])
    derivative = grid["derivative"]
    radius = np.asarray(grid["radius"], dtype=float)

    def herm(left, right):
        return source._hermitian_flux(left, right, derivative)

    plus = herm(column[0], column[0])
    minus = herm(column[1], column[1])
    d_plus = herm(column[0], direction[0]) + herm(direction[0], column[0])
    d_minus = herm(column[1], direction[1]) + herm(direction[1], column[1])
    angular = source._angular_bilinear(column, column)
    d_angular = (
        source._angular_bilinear(column, direction) + source._angular_bilinear(direction, column)
    )
    flux_plus = nu * np.real(plus)
    flux_minus = nu * np.real(minus)
    angular_scalar = nu * np.real(angular)
    d_flux_plus = dnu * np.real(plus) + nu * np.real(d_plus)
    d_flux_minus = dnu * np.real(minus) + nu * np.real(d_minus)
    d_angular_scalar = dnu * np.real(angular) + nu * np.real(d_angular)
    kappa = float(kappa)
    force_n = multiplicity * ((flux_plus - flux_minus) + (kappa / radius) * angular_scalar)
    force_beta = multiplicity * (-(flux_plus + flux_minus))
    force_r = multiplicity * (-(kappa / radius**2) * angular_scalar)
    d_force_n = multiplicity * ((d_flux_plus - d_flux_minus) + (kappa / radius) * d_angular_scalar)
    d_force_beta = multiplicity * (-(d_flux_plus + d_flux_minus))
    d_force_r = multiplicity * (-(kappa / radius**2) * d_angular_scalar)
    base = pg.matter_densities(phi, grid, occupations, kappa=int(kappa))
    gap = max(
        float(np.max(np.abs(force_n - base["eps_N"]))),
        float(np.max(np.abs(force_beta - base["eps_beta"]))),
        float(np.max(np.abs(force_r - base["force_r_unit"]))),
    )
    scale = 1.0 + float(np.max(np.abs(base["eps_N"])))
    if gap > 1e-8 * scale:
        raise RuntimeError("stress tangent does not reproduce matter_densities")
    return {
        "force_n": np.asarray(base["eps_N"], dtype=float),
        "force_beta": np.asarray(base["eps_beta"], dtype=float),
        "force_r_unit": np.asarray(base["force_r_unit"], dtype=float),
        "d_force_n": np.asarray(d_force_n, dtype=float),
        "d_force_beta": np.asarray(d_force_beta, dtype=float),
        "d_force_r_unit": np.asarray(d_force_r, dtype=float),
        "multiplicity": multiplicity,
        "multiplicity_applied_once": True,
        "reconstruction_weight": RECONSTRUCTION_WEIGHT,
        "cached_at_unit_lapse": True,
    }


def reconstruction_tangent(grid, inner_mass, force_n, force_beta, force_r,
                           d_inner_mass, d_force_n, d_force_beta, d_force_r, *,
                           weight=RECONSTRUCTION_WEIGHT):
    """Differentiated trapezoid. Weight stays 1. The second v F_N stays in the flux."""
    if float(weight) != RECONSTRUCTION_WEIGHT:
        raise ValueError("reconstruction tangent uses weight 1; the shell factor is already in F")
    geometry = pg.reconstruct_geometry(
        grid, inner_mass, force_n, force_beta, -np.asarray(force_r, dtype=float) * grid["radius"],
        weight=weight, force_r_unit=force_r,
    )
    geometry = dict(geometry)
    geometry["eps_beta_bare"] = np.asarray(force_beta, dtype=float)
    d_mass, d_lapse, d_shift = pg._linearised_rates(
        grid, geometry, d_force_n, d_force_beta, d_inner_mass, weight,
    )
    radius = np.asarray(grid["radius"], dtype=float)
    velocity = np.asarray(geometry["shift_factor"], dtype=float)
    d_velocity = d_mass / (radius * velocity)
    force_n = np.asarray(force_n, dtype=float)
    force_beta = np.asarray(force_beta, dtype=float)
    force_r = np.asarray(force_r, dtype=float)
    d_force_n = np.asarray(d_force_n, dtype=float)
    d_force_beta = np.asarray(d_force_beta, dtype=float)
    d_force_r = np.asarray(d_force_r, dtype=float)
    bracket = (
        2.0 * velocity * force_n
        + (1.0 + velocity * velocity) * force_beta
        + velocity * radius * force_r
    )
    d_bracket = (
        2.0 * d_velocity * force_n + 2.0 * velocity * d_force_n
        + 2.0 * velocity * d_velocity * force_beta + (1.0 + velocity * velocity) * d_force_beta
        + d_velocity * radius * force_r + velocity * radius * d_force_r
    )
    flux = pg.NEWTON_G * float(weight) * geometry["lapse"] * bracket
    d_flux = pg.NEWTON_G * float(weight) * (d_lapse * bracket + geometry["lapse"] * d_bracket)
    if float(np.max(np.abs(flux - geometry["q_equation_mass_rate"]))) > 1e-8 * (
        1.0 + float(np.max(np.abs(geometry["q_equation_mass_rate"])))
    ):
        raise RuntimeError("flux tangent base missed the second v F_N copy")
    if abs(float(d_lapse[-1])) > 1e-12:
        raise RuntimeError("outer lapse tangent is fixed by N(r_out) = 1")
    return {
        "mass": np.asarray(geometry["mass"], dtype=float),
        "lapse": np.asarray(geometry["lapse"], dtype=float),
        "shift": np.asarray(geometry["shift"], dtype=float),
        "shift_factor": velocity,
        "q_equation_mass_rate": np.asarray(geometry["q_equation_mass_rate"], dtype=float),
        "d_mass": np.asarray(d_mass, dtype=float),
        "d_lapse": np.asarray(d_lapse, dtype=float),
        "d_shift": np.asarray(d_shift, dtype=float),
        "d_shift_factor": np.asarray(d_velocity, dtype=float),
        "d_q_equation_mass_rate": np.asarray(d_flux, dtype=float),
        "dN_out": float(d_lapse[-1]),
        "reconstruction_weight": float(weight),
        "second_v_fn_coefficient": 2.0,
        "force_r_enters_mass_profile": False,
    }


def clock_semantics():
    """Labels for the two base scalars. Neither scalar moves the metric.

    ``killing`` is the static proper clock at the fixed outer node.
    ``normal`` is the legacy integral of the excision lapse. It is not a
    transported observer and its radial curve is spacelike between the
    horizons. Transported normal proper time belongs to another owner.
    """
    return {
        "killing": {
            "protocol": OUTER_STATIC_PROTOCOL,
            "proper_clock": True,
            "node": "outer",
            "formula": "sqrt(N[-1]**2 - shift[-1]**2)",
            "valid_while_gtt_positive": True,
            "moves_the_field": False,
        },
        "normal": {
            "protocol": LEGACY_EXCISION_LAPSE,
            "proper_clock": False,
            "node": "excision",
            "formula": "N[0]",
            "observer_radius_transported": False,
            "curve_is_spacelike_between_horizons": True,
            "moves_the_field": False,
            "replay_tangent_retained": True,
        },
        "admitted_proper_protocols": (OUTER_STATIC_PROTOCOL,),
        "transported_normal_protocol": TRANSPORTED_NORMAL_PROTOCOL,
        "transported_normal_available": False,
    }


def _clock_rates(lapse, shift, d_lapse, d_shift):
    """Replay tangents of the base clocks. The numbers are the PG scalars."""
    metric_tt = float(lapse[-1] ** 2 - shift[-1] ** 2)
    if metric_tt <= 0.0:
        raise ValueError("Killing clock rate is not real on this slice")
    killing = math.sqrt(metric_tt)
    d_killing = float(lapse[-1] * d_lapse[-1] - shift[-1] * d_shift[-1]) / killing
    return {
        "t": 1.0,
        "killing": killing,
        "normal": float(lapse[0]),
        "d_t": 0.0,
        "d_killing": d_killing,
        "d_normal": float(d_lapse[0]),
        "semantics": clock_semantics(),
    }


def dirac_tangent(grid, phi, dphi, lapse, d_lapse, shift, d_shift, *, kappa=SHELL_KAPPA):
    """Characteristic SBP plus outer SAT, differentiated with pg._split."""
    array = _column(phi)
    direction = _column(dphi)
    radius = np.asarray(grid["radius"], dtype=float)
    alpha = np.asarray(lapse, dtype=float) * float(kappa) / radius
    d_alpha = np.asarray(d_lapse, dtype=float) * float(kappa) / radius
    plus = array[0, :, 0]
    minus = array[1, :, 0]
    d_plus = direction[0, :, 0]
    d_minus = direction[1, :, 0]

    def d_split(coeff, d_coeff, component, d_component):
        return pg._split(grid, d_coeff, component) + pg._split(grid, coeff, d_component)

    d_plus_rate = (
        -d_split(lapse, d_lapse, plus, d_plus)
        + d_split(shift, d_shift, plus, d_plus)
        - d_alpha * minus - alpha * d_minus
    )
    d_minus_rate = (
        d_split(lapse, d_lapse, minus, d_minus)
        + d_split(shift, d_shift, minus, d_minus)
        + d_alpha * plus + alpha * d_plus
    )
    incoming = float(-shift[-1] - lapse[-1])
    d_incoming = float(-d_shift[-1] - d_lapse[-1])
    weight = float(grid["weights"][-1])
    d_minus_rate[-1] += (d_incoming / weight) * minus[-1] + (incoming / weight) * d_minus[-1]
    rates = np.zeros_like(array)
    rates[0, :, 0] = d_plus_rate
    rates[1, :, 0] = d_minus_rate
    sat = float(-incoming * np.abs(minus[-1]) ** 2)
    d_sat = float(
        -d_incoming * np.abs(minus[-1]) ** 2
        - incoming * 2.0 * np.real(np.conj(minus[-1]) * d_minus[-1])
    )
    return rates, sat, d_sat


def probability_boundary_tangent(phi, dphi, lapse, d_lapse, shift, d_shift):
    """Probability current. The shell factor is not in this stock."""
    array = _column(phi)
    direction = _column(dphi)
    outgoing = -np.asarray(shift, dtype=float) + np.asarray(lapse, dtype=float)
    ingoing = -np.asarray(shift, dtype=float) - np.asarray(lapse, dtype=float)
    d_outgoing = -np.asarray(d_shift, dtype=float) + np.asarray(d_lapse, dtype=float)
    d_ingoing = -np.asarray(d_shift, dtype=float) - np.asarray(d_lapse, dtype=float)
    plus = np.abs(array[0, :, 0]) ** 2
    minus = np.abs(array[1, :, 0]) ** 2
    d_plus = 2.0 * np.real(np.conj(array[0, :, 0]) * direction[0, :, 0])
    d_minus = 2.0 * np.real(np.conj(array[1, :, 0]) * direction[1, :, 0])
    current = outgoing * plus + ingoing * minus
    d_current = d_outgoing * plus + outgoing * d_plus + d_ingoing * minus + ingoing * d_minus
    return {
        "excision_outflow": float(-current[0]),
        "outer_outflow": float(current[-1]),
        "d_excision_outflow": float(-d_current[0]),
        "d_outer_outflow": float(d_current[-1]),
        "multiplicity_included": False,
    }


def stage_tangent(state, grid, direction):
    """Directional derivative of one stage_rates evaluation."""
    if int(state["kappa"]) != SHELL_KAPPA or int(state["multiplicity"]) != SHELL_MULTIPLICITY:
        raise ValueError("stage tangent is the closed κ=1 shell")
    phi = _column(state["phi"])
    dphi = _column(direction["dphi"])
    if dphi.shape != phi.shape:
        raise ValueError("dphi does not match phi")
    occupations = np.asarray(state["occupations"], dtype=float).reshape(-1)
    d_occupations = np.asarray(direction["d_occupations"], dtype=float).reshape(-1)
    rates = pg.stage_rates(state, grid)
    stress = stress_tangent(grid, phi, occupations, dphi, d_occupations, kappa=state["kappa"])
    geometry = reconstruction_tangent(
        grid, state["inner_mass"], stress["force_n"], stress["force_beta"], stress["force_r_unit"],
        direction["d_inner_mass"], stress["d_force_n"], stress["d_force_beta"], stress["d_force_r_unit"],
    )
    if np.max(np.abs(geometry["lapse"] - rates["lapse"])) > 1e-8:
        raise RuntimeError("reconstruction tangent base does not match stage_rates")
    d_phi_rate, sat_base, d_sat = dirac_tangent(
        grid, phi, dphi, geometry["lapse"], geometry["d_lapse"], geometry["shift"], geometry["d_shift"],
        kappa=state["kappa"],
    )
    if np.max(np.abs(d_phi_rate * 0.0 + rates["phi_rate"] - pg.dirac_rates(
        grid, phi, geometry["lapse"], geometry["shift"], kappa=state["kappa"],
    )[0])) > 1e-8:
        raise RuntimeError("Dirac base does not match the reconstructed metric")
    if abs(sat_base - float(rates["sat_debit_rate"])) > 1e-8 * (1.0 + abs(float(rates["sat_debit_rate"]))):
        raise RuntimeError("SAT base does not match stage_rates")
    boundary = probability_boundary_tangent(
        phi, dphi, geometry["lapse"], geometry["d_lapse"], geometry["shift"], geometry["d_shift"],
    )
    clocks = _clock_rates(geometry["lapse"], geometry["shift"], geometry["d_lapse"], geometry["d_shift"])
    if abs(clocks["killing"] - float(rates["killing_rate"])) > 1e-10:
        raise RuntimeError("Killing rate tangent base does not match stage_rates")
    if abs(clocks["normal"] - float(rates["normal_rate"])) > 1e-12:
        raise RuntimeError("legacy excision lapse tangent base does not match stage_rates")
    weights = grid["weights"]
    coordinate = coordinate_energy_tangent(
        weights, geometry["lapse"], geometry["d_lapse"], geometry["shift"], geometry["d_shift"],
        stress["force_n"], stress["d_force_n"], stress["force_beta"], stress["d_force_beta"],
    )
    return {
        "rates": rates,
        "dphi_rate": d_phi_rate,
        "d_inner_mass_rate": float(geometry["d_q_equation_mass_rate"][0]),
        "d_outer_flux_rate": float(geometry["d_q_equation_mass_rate"][-1]),
        "d_excision_flux_rate": float(geometry["d_q_equation_mass_rate"][0]),
        "d_sat_debit_rate": d_sat,
        "d_probability_excision_outflow_rate": boundary["d_excision_outflow"],
        "d_probability_outer_outflow_rate": boundary["d_outer_outflow"],
        "d_t_rate": clocks["d_t"],
        "d_killing_rate": clocks["d_killing"],
        "d_normal_rate": clocks["d_normal"],
        "d_lapse": geometry["d_lapse"],
        "d_shift": geometry["d_shift"],
        "d_mass": geometry["d_mass"],
        "dN_out": geometry["dN_out"],
        "d_force_n": stress["d_force_n"],
        "d_force_beta": stress["d_force_beta"],
        "d_force_r_unit": stress["d_force_r_unit"],
        "lapse": geometry["lapse"],
        "shift": geometry["shift"],
        "mass": geometry["mass"],
        "force_n": stress["force_n"],
        "force_beta": stress["force_beta"],
        "coordinate_energy": coordinate["energy"],
        "d_coordinate_energy": coordinate["d_energy"],
        "shell_weight": RECONSTRUCTION_WEIGHT,
        "state_multiplicity": int(state["multiplicity"]),
        "multiplicity_in_dirac": False,
        "multiplicity_in_probability": False,
        "second_v_fn_coefficient": 2.0,
        "dt_held_fixed": True,
        "clock_semantics": clocks["semantics"],
        "clocks_move_the_metric": False,
    }


def _rk4_combine(first, second, third, fourth):
    return (first + 2.0 * second + 2.0 * third + fourth) / 6.0


def rk4_tangent(state, grid, dt, direction):
    """Four differentiated reconstructions. dt is held fixed."""
    step = float(dt)
    current = stage_tangent(state, grid, direction)
    if step > float(current["rates"]["cfl_dt"]) * (1.0 + 1e-12):
        raise pg.ChartAdmissionError("step exceeds the full coupled CFL")
    bundles = [current]
    # Same arguments as rk4_step: each stage starts from the original state.
    for fraction in (0.5, 0.5, 1.0):
        previous = bundles[-1]
        moved = pg._advance(state, previous["rates"], fraction * step, grid)
        moved_direction = dict(direction)
        moved_direction["dphi"] = direction["dphi"] + fraction * step * previous["dphi_rate"]
        moved_direction["d_inner_mass"] = (
            float(direction["d_inner_mass"]) + fraction * step * previous["d_inner_mass_rate"]
        )
        bundles.append(stage_tangent(moved, grid, moved_direction))
    first, second, third, fourth = bundles
    dphi = direction["dphi"] + step * _rk4_combine(
        first["dphi_rate"], second["dphi_rate"], third["dphi_rate"], fourth["dphi_rate"],
    )
    d_inner = float(direction["d_inner_mass"]) + step * _rk4_combine(
        first["d_inner_mass_rate"], second["d_inner_mass_rate"],
        third["d_inner_mass_rate"], fourth["d_inner_mass_rate"],
    )
    d_clocks = {}
    for name, rate_name in (("t", "d_t_rate"), ("killing", "d_killing_rate"), ("normal", "d_normal_rate")):
        d_clocks[name] = float(direction["d_clocks"][name]) + step * _rk4_combine(
            first[rate_name], second[rate_name], third[rate_name], fourth[rate_name],
        )
    d_stocks = {
        "excision_flux": float(direction["d_stocks"]["excision_flux"]) + step * _rk4_combine(
            first["d_excision_flux_rate"], second["d_excision_flux_rate"],
            third["d_excision_flux_rate"], fourth["d_excision_flux_rate"],
        ),
        "outer_flux": float(direction["d_stocks"]["outer_flux"]) + step * _rk4_combine(
            first["d_outer_flux_rate"], second["d_outer_flux_rate"],
            third["d_outer_flux_rate"], fourth["d_outer_flux_rate"],
        ),
        "sat_debit": float(direction["d_stocks"]["sat_debit"]) + step * _rk4_combine(
            first["d_sat_debit_rate"], second["d_sat_debit_rate"],
            third["d_sat_debit_rate"], fourth["d_sat_debit_rate"],
        ),
    }
    phi = state["phi"] + (step / 6.0) * (
        first["rates"]["phi_rate"] + 2.0 * second["rates"]["phi_rate"]
        + 2.0 * third["rates"]["phi_rate"] + fourth["rates"]["phi_rate"]
    )
    d_remaining = 2.0 * real_overlap(phi, dphi, grid["weights"])
    d_probability = {
        "remaining": d_remaining,
        "excision_outflow": float(direction["d_probability_stocks"]["excision_outflow"]) + step * _rk4_combine(
            first["d_probability_excision_outflow_rate"], second["d_probability_excision_outflow_rate"],
            third["d_probability_excision_outflow_rate"], fourth["d_probability_excision_outflow_rate"],
        ),
        "outer_outflow": float(direction["d_probability_stocks"]["outer_outflow"]) + step * _rk4_combine(
            first["d_probability_outer_outflow_rate"], second["d_probability_outer_outflow_rate"],
            third["d_probability_outer_outflow_rate"], fourth["d_probability_outer_outflow_rate"],
        ),
    }
    d_closure = (
        d_probability["remaining"] + d_probability["excision_outflow"]
        + d_probability["outer_outflow"] + 2.0 * d_stocks["sat_debit"]
    )
    cached = direction.get("d_cached_covariance")
    return {
        "arm": "coupled",
        "dt": step,
        "dt_held_fixed": True,
        "dphi": dphi,
        "phi": phi,
        "d_inner_mass": d_inner,
        "d_clocks": d_clocks,
        "d_stocks": d_stocks,
        "d_probability_stocks": d_probability,
        "d_occupations": np.array(direction["d_occupations"], dtype=float, copy=True),
        "d_cached_covariance": None if cached is None else np.array(cached, dtype=complex, copy=True),
        "cached_covariance_refreshed": False,
        "d_closure": d_closure,
        "dN_out": float(first["dN_out"]),
        "stages": 4,
        "source_mass_clocks_and_stocks": True,
        "clock_semantics": clock_semantics(),
        "clocks_move_the_metric": False,
        "field_rk4_stage_fractions": FIELD_RK4_STAGE_FRACTIONS,
    }


def frozen_arm_tangent(state, grid, dt, direction):
    """Dirac RK4 on the initial sourced metric, including that metric's tangent.

    The lapse and shift do not advance. Their parameter derivatives are the
    initial reconstruction derivatives and stay present at every stage.
    """
    initial = stage_tangent(state, grid, direction)
    lapse = initial["lapse"]
    shift = initial["shift"]
    d_lapse = initial["d_lapse"]
    d_shift = initial["d_shift"]
    step = float(dt)
    cfl = observables.full_cfl_dt(lapse, shift, grid["spacing"])
    if step > cfl * (1.0 + 1e-12):
        raise ValueError("fixed Dirac step exceeds the frozen-metric CFL")
    kappa = state["kappa"]
    phi0 = _column(state["phi"])
    dphi0 = _column(direction["dphi"])

    def sample(phi_now, dphi_now):
        phi_rate, sat = pg.dirac_rates(grid, phi_now, lapse, shift, kappa=kappa)
        d_phi_rate, sat_check, d_sat = dirac_tangent(
            grid, phi_now, dphi_now, lapse, d_lapse, shift, d_shift, kappa=kappa,
        )
        if abs(sat - sat_check) > 1e-10 * (1.0 + abs(sat)):
            raise RuntimeError("frozen SAT base drifted from dirac_rates")
        boundary = probability_boundary_tangent(phi_now, dphi_now, lapse, d_lapse, shift, d_shift)
        return phi_rate, d_phi_rate, sat, d_sat, boundary

    k1, dk1, s1, ds1, b1 = sample(phi0, dphi0)
    k2, dk2, s2, ds2, b2 = sample(phi0 + 0.5 * step * k1, dphi0 + 0.5 * step * dk1)
    k3, dk3, s3, ds3, b3 = sample(phi0 + 0.5 * step * k2, dphi0 + 0.5 * step * dk2)
    k4, dk4, s4, ds4, b4 = sample(phi0 + step * k3, dphi0 + step * dk3)
    dphi = dphi0 + (step / 6.0) * (dk1 + 2.0 * dk2 + 2.0 * dk3 + dk4)
    phi = phi0 + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    clocks = _clock_rates(lapse, shift, d_lapse, d_shift)
    d_sat = (ds1 + 2.0 * ds2 + 2.0 * ds3 + ds4) / 6.0
    d_excision = (
        b1["d_excision_outflow"] + 2.0 * b2["d_excision_outflow"]
        + 2.0 * b3["d_excision_outflow"] + b4["d_excision_outflow"]
    ) / 6.0
    d_outer = (
        b1["d_outer_outflow"] + 2.0 * b2["d_outer_outflow"]
        + 2.0 * b3["d_outer_outflow"] + b4["d_outer_outflow"]
    ) / 6.0
    d_probability = {
        "remaining": 2.0 * real_overlap(phi, dphi, grid["weights"]),
        "excision_outflow": float(direction["d_probability_stocks"]["excision_outflow"]) + step * d_excision,
        "outer_outflow": float(direction["d_probability_stocks"]["outer_outflow"]) + step * d_outer,
    }
    d_stocks = {
        "excision_flux": float(direction["d_stocks"]["excision_flux"]),
        "outer_flux": float(direction["d_stocks"]["outer_flux"]),
        "sat_debit": float(direction["d_stocks"]["sat_debit"]) + step * d_sat,
    }
    d_clocks = {
        "t": float(direction["d_clocks"]["t"]),
        "killing": float(direction["d_clocks"]["killing"]) + step * clocks["d_killing"],
        "normal": float(direction["d_clocks"]["normal"]) + step * clocks["d_normal"],
    }
    return {
        "arm": "frozen",
        "dt": step,
        "dphi": dphi,
        "phi": phi,
        "d_inner_mass": float(direction["d_inner_mass"]),
        "inner_mass_advanced": False,
        "d_clocks": d_clocks,
        "d_stocks": d_stocks,
        "d_probability_stocks": d_probability,
        "lapse": lapse,
        "shift": shift,
        "d_lapse": d_lapse,
        "d_shift": d_shift,
        "metric_held_fixed_in_time": True,
        "metric_variation_included": True,
        "variation_is_initial_reconstruction": True,
        "mass_flux_advanced": False,
        "dN_out": float(d_lapse[-1]),
        "stages": 4,
        "clock_semantics": clocks["semantics"],
        "clocks_move_the_metric": False,
    }


def declared_proper_clock(stage_bundle, protocol=OUTER_STATIC_PROTOCOL):
    """Proper rate stored on one stage. The legacy normal scalar is refused."""
    if protocol == LEGACY_EXCISION_LAPSE or protocol == "normal":
        raise ValueError("legacy fixed excision lapse integral is not a proper clock")
    if protocol == TRANSPORTED_NORMAL_PROTOCOL:
        raise ValueError("transported normal proper time is not supplied by this adapter")
    if protocol != OUTER_STATIC_PROTOCOL:
        raise ValueError("proper-time comparison needs the outer static protocol")
    semantics = stage_bundle["clock_semantics"]["killing"]
    if semantics["proper_clock"] is not True or semantics["protocol"] != OUTER_STATIC_PROTOCOL:
        raise ValueError("outer static Killing clock is not marked proper on this stage")
    return {
        "protocol": OUTER_STATIC_PROTOCOL,
        "proper_clock": True,
        "node": "outer",
        "rate": float(stage_bundle["rates"]["killing_rate"]),
        "d_rate": float(stage_bundle["d_killing_rate"]),
        "legacy_normal_rate": float(stage_bundle["rates"]["normal_rate"]),
        "legacy_normal_is_proper": False,
    }


def matched_proper_derivative(d_observable_t, observable_rate, d_clock, clock_rate, *, protocol):
    """Chain rule dO_τ = dO_t − Odot dτ / τdot for one declared proper protocol.

    Admitted protocol: ``outer_static``, the fixed outer Killing clock.
    The legacy excision lapse integral is not admitted. Transported normal
    proper time is a separate owner and is not evaluated here.
    """
    if protocol == LEGACY_EXCISION_LAPSE or protocol == "normal":
        raise ValueError("legacy fixed excision lapse integral is not a proper clock")
    if protocol == TRANSPORTED_NORMAL_PROTOCOL:
        raise ValueError("transported normal proper time is not supplied by this adapter")
    if protocol != OUTER_STATIC_PROTOCOL:
        raise ValueError("proper-time comparison needs the outer static protocol")
    rate = float(clock_rate)
    if rate == 0.0 or not math.isfinite(rate):
        raise ValueError("matched proper time needs a nonzero clock rate")
    return float(d_observable_t) - float(observable_rate) * float(d_clock) / rate


def observer_metric_stage(stage_bundle, grid):
    """Reconstructed lapse and shift one field stage exposes to an observer.

    The observer samples these profiles on the grid radius. It does not read
    the legacy normal scalar, and this adapter does not advance an observer
    radius. Field RK4 uses ``FIELD_RK4_STAGE_FRACTIONS`` from the original
    state, with weights 1, 2, 2, 1.
    """
    radius = np.asarray(grid["radius"], dtype=float)
    lapse = np.asarray(stage_bundle["lapse"], dtype=float)
    shift = np.asarray(stage_bundle["shift"], dtype=float)
    d_lapse = np.asarray(stage_bundle["d_lapse"], dtype=float)
    d_shift = np.asarray(stage_bundle["d_shift"], dtype=float)
    if lapse.shape != radius.shape or shift.shape != radius.shape:
        raise ValueError("observer metric stage does not match the grid")
    if d_lapse.shape != radius.shape or d_shift.shape != radius.shape:
        raise ValueError("observer metric tangent does not match the grid")
    return {
        "radius": radius,
        "lapse": lapse,
        "shift": shift,
        "d_lapse": d_lapse,
        "d_shift": d_shift,
        "excision_index": 0,
        "outer_index": -1,
        "field_rk4_stage_fractions": FIELD_RK4_STAGE_FRACTIONS,
        "field_rk4_weights": (1.0, 2.0, 2.0, 1.0),
        "legacy_normal_not_a_sample": True,
        "clocks_move_the_metric": False,
        "transported_clock_computed": False,
    }


def current_car_tangent(phi, dphi, weights, filling, d_filling):
    """Current readout ν G. The cached load-time CAR is not this number."""
    gram = source.weighted_gram(_column(phi), weights)[0, 0]
    d_gram = 2.0 * real_overlap(phi, dphi, weights)
    value = float(filling) * gram
    d_value = float(d_filling) * gram + float(filling) * d_gram
    return {
        "gram": complex(gram),
        "d_gram": d_gram,
        "value": complex(value),
        "d_value": complex(d_value),
        "multiplicity_included": False,
        "cached_covariance_used": False,
    }


def coordinate_energy_tangent(weights, lapse, d_lapse, shift, d_shift, force_n, d_force_n,
                              force_beta, d_force_beta):
    """E = Σ w (N F_N + β F_β) on the reconstructed slice."""
    weights = np.asarray(weights, dtype=float)
    density = np.asarray(lapse, dtype=float) * np.asarray(force_n, dtype=float) + (
        np.asarray(shift, dtype=float) * np.asarray(force_beta, dtype=float)
    )
    d_density = (
        np.asarray(d_lapse, dtype=float) * np.asarray(force_n, dtype=float)
        + np.asarray(lapse, dtype=float) * np.asarray(d_force_n, dtype=float)
        + np.asarray(d_shift, dtype=float) * np.asarray(force_beta, dtype=float)
        + np.asarray(shift, dtype=float) * np.asarray(d_force_beta, dtype=float)
    )
    return {
        "energy": float(np.sum(weights * density)),
        "d_energy": float(np.sum(weights * d_density)),
    }


def trapping_horizon_tangent(radius, lapse, shift, d_lapse, d_shift):
    """Fixed-index derivative of the outermost outgoing-speed zero.

    A finite control that changes the crossing node is outside this formula.
    """
    outgoing, _ingoing = observables.chart_speeds(lapse, shift)
    d_outgoing = -np.asarray(d_shift, dtype=float) + np.asarray(d_lapse, dtype=float)
    radius = np.asarray(radius, dtype=float)
    crossing = None
    for index in range(1, radius.size):
        left, right = float(outgoing[index - 1]), float(outgoing[index])
        if left < 0.0 <= right and right != left:
            span = right - left
            fraction = -left / span
            d_left = float(d_outgoing[index - 1])
            d_right = float(d_outgoing[index])
            d_fraction = (-d_left * span + left * (d_right - d_left)) / span**2
            crossing = {
                "radius": float(radius[index - 1] + fraction * (radius[index] - radius[index - 1])),
                "d_radius": float(d_fraction * (radius[index] - radius[index - 1])),
                "index": index,
                "fixed_index": True,
            }
    if crossing is None:
        raise ValueError("trapping horizon tangent has no fixed-index crossing")
    return crossing


def displace_state(state, direction, scale):
    """Affine probe used by centered differences. Not a time step."""
    scale = float(scale)
    phi = _column(state["phi"]) + scale * _column(direction["dphi"])
    occupations = np.asarray(state["occupations"], dtype=float) + scale * np.asarray(
        direction["d_occupations"], dtype=float,
    )
    if np.any(occupations < 0.0) or np.any(occupations >= 1.0):
        raise ValueError("displaced filling left [0, 1)")
    inner = float(state["inner_mass"]) + scale * float(direction["d_inner_mass"])
    clocks = {
        name: float(state["clocks"][name]) + scale * float(direction["d_clocks"][name])
        for name in ("t", "killing", "normal")
    }
    stocks = {
        name: float(state["stocks"][name]) + scale * float(direction["d_stocks"][name])
        for name in ("excision_flux", "outer_flux", "sat_debit")
    }
    probability = {
        name: float(state["probability_stocks"][name]) + scale * float(direction["d_probability_stocks"][name])
        for name in ("remaining", "excision_outflow", "outer_outflow")
    }
    moved = pg._copy_state(state, phi, inner, clocks, stocks, probability)
    moved["occupations"] = np.array(occupations, dtype=float, copy=True)
    moved["occupation_matrix"] = None
    return moved


def baseline_window(model=None):
    """Truncation of the nine-sample window at μ = 1.04. No radial ODE."""
    if model is None:
        model = pg.action_model(BASELINE_MASS_OVER_RM, 1.0)
    request = source.packet_request(source.magnetic_radius(model))
    chart = source._point_chart(model, request["center"])
    speed = float(chart["speed_minus"])
    width = float(request["width"])
    omega = float(request["omega"])
    spectral_width = abs(speed) / width
    energies, truncated, span = source._positive_frequency_window(omega, spectral_width, BASELINE_RANK)
    wider, wider_truncated, wider_span = source._positive_frequency_window(
        omega, abs(speed) / (1.05 * width), BASELINE_RANK,
    )
    faster, faster_truncated, faster_span = source._positive_frequency_window(
        1.05 * omega, spectral_width, BASELINE_RANK,
    )
    ell = np.linspace(-1.0, 1.0, BASELINE_RANK)
    return {
        "ingoing_speed": speed,
        "spectral_width": spectral_width,
        "truncated": bool(truncated),
        "span": float(span),
        "energies": energies,
        "energies_match_omega_rule": bool(np.max(np.abs(energies - omega * (1.0 + 0.9 * ell))) < 1e-12),
        "plus_five_percent_width_truncated": bool(wider_truncated),
        "plus_five_percent_width_span": float(wider_span),
        "plus_five_percent_frequency_truncated": bool(faster_truncated),
        "plus_five_percent_frequency_span": float(faster_span),
        "width_nodes_fixed": True,
        "frequency_nodes_move": True,
        "held_forecast": "not_run",
        "width_shifted_energies": wider,
        "frequency_shifted_energies": faster,
    }
