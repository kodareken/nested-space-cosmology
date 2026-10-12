"""Coordinate Hamiltonian ledger and normal-observer stress ledger.

The owned generator remains
∫(L C + β D) dx = gravity energy + M Tr(C H), M = 4κ once.
Nodal FL, FQ, Fβ are partial derivatives of that trace, not stresses.
Divide by the node spacing before building a density.

Coordinate matter integrand, a = L/Q:

    e = L FL + β Fβ = M (a K + m S − β J)
    ∂t e + ∂x(Φ_shift + Φ_proper + Φ_cross) = FQ Q̇
    Φ_shift = −β e
    Φ_proper = M a² J
    Φ_cross = −M a β K = β Q FQ

Φ_cross is the shift/conformal-pressure cross term. L and β are held fixed,
so the coordinate work that remains is FQ Q̇, including Q̇ = ∂x(β Q).
That is not observer pressure work.

Normal-observer shell, N = r L, q = r Q, V = 4π q r²:

    nodal normal energy = FL / r = V ρ dx
    ρ = FL / (4π r^4 Q dx)
    p_r = −FQ / (4π r^4 L dx)
    p_⊥ = (L FL + Q FQ) / (8π r^4 L Q dx)
    j = −Fβ / (4π r^4 Q² dx)
    ∂t(Vρ) + ∂x[V (N j/q − β ρ)] = −N V (p_r K_r + 2 p_⊥ K_⊥) − V (j/q) N_x
    K_r = (q̇ − ∂x(β q)) / (N q),   K_⊥ = (ṙ − β r_x) / (N r)

L FL = N (FL/r) is lapse-weighted shell energy, not ρ.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

from .nsc_conformal_adm_source import physical_forces_from_conformal
from .nsc_spherical_coupling import (
    CauchyRate,
    CauchyState,
    _combine,
    build_system,
    field_energy,
    geometric_rates,
    gravity_energy,
    hamilton_constraint,
    rates,
    shift_constraint,
    source_from_columns,
    total_energy,
)
from .nsc_spherical_feedback_action import feedback_F, feedback_V, feedback_Z, partial_F
from .nsc_spherical_galerkin_coupling import (
    blank_state,
    build_grid,
    compose_fine_hamiltonian,
    manufactured_columns,
    prolong_geometry,
    rk4_step,
)

SCHEMA = "NSC-REGIONAL-ENERGY-EXCHANGE-v1"
RENEWAL = False
CONTINUUM_LIMIT_CLAIMED = False
FROZEN_B_EMBEDDING_CLAIMED = False
OLD_CIRCULATION_TRANSFERRED = False
VACUUM_BRANCH_INCLUDED = False
ACCEPTED_TRAJECTORY = False

# Smooth band-limited windows. Not retuned against the v5 residual.
TOL_SMOOTH_WINDOW = 1e-5
TOL_SUM = 1e-8
TOL_WARD = 1e-8
TOL_UNPROJECTED = 1e-7
TOL_PROPER = 5e-5

_LAB_ROOT = Path(__file__).resolve().parents[2]
V5_NPZ = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-refinement-v5.npz"
DIAGNOSTIC_JSON = _LAB_ROOT / "results" / "development" / "nsc-regional-energy-exchange-v1.json"
_MODULE_PATH = Path(__file__).resolve()

GLOBAL_ENERGY_MEANING = (
    "Global energy is the constrained coordinate Hamiltonian: "
    "gravity_energy + M Tr(C H) = dx * sum(L C_tot + beta D_tot). "
    "A value near zero means that smearing of the constraints is small. "
    "It is not the normal-observer energy and it does not say the matter integral vanishes."
)
REGIONAL_ENERGY_MEANING = (
    "Regional coordinate energy is that Hamiltonian density on a smooth window. "
    "Matter and gravity are each large; their sum differs from the windowed "
    "constraint smearing by the summation-by-parts term of 2 L (D F)/Q. "
    "A partition of unity reproduces the global Hamiltonian. "
    "The normal-observer shell energy is the separate nodal sum FL/r."
)
# Coordinate-window figures from the draft this correction replaces.
# They are compared, not used as new targets.
_RETAINED_COORDINATE = {
    "matter_closure_error": 9.864223882161127e-10,
    "transported_closure_error": 9.864782324342514e-10,
    "total_closure_error": -3.073359300387324e-08,
    "shift_transport": -0.18691093233529965,
    "proper_normal_flux": 1.3559256963322972,
    "shift_pressure_cross": -0.3244579548048109,
    "coordinate_metric_work": 0.10814440302646627,
}


def _coordinate_reference_gap(smooth):
    """How far the corrected run sits from the retained coordinate draft."""
    current = {
        "matter_closure_error": smooth["matter_closure_error"],
        "transported_closure_error": smooth["transported_closure_error"],
        "total_closure_error": smooth["total_closure_error"],
        "shift_transport": smooth["channels"]["shift_transport"],
        "proper_normal_flux": smooth["channels"]["proper_normal_flux"],
        "shift_pressure_cross": smooth["channels"]["shift_pressure_cross"],
        "coordinate_metric_work": smooth["channels"]["coordinate_metric_work"],
    }
    gap = {key: float(current[key] - value) for key, value in _RETAINED_COORDINATE.items()}
    return {
        "retained": _RETAINED_COORDINATE,
        "absolute_gap": gap,
        "max_abs_gap": float(max(abs(value) for value in gap.values())),
        "matches_previous_coordinate_draft": bool(max(abs(value) for value in gap.values()) < 1e-9),
    }


CORRECTION = (
    "An earlier draft called L*FL observer energy, FQ/r radial pressure, "
    "the chain-rule radius derivative angular pressure, and beta*Q*FQ a "
    "mode-channel current. Those are Hamiltonian derivatives. "
    "Observer energy is FL/r. Stresses are rho, p_r, p_perp, j below. "
    "beta*Q*FQ is the shift/pressure cross term in the coordinate flux. "
    "FQ*Qdot is coordinate metric work and stays nonzero when K_r = K_perp = 0."
)


def _sha256(path):
    digest = hashlib.sha256()
    file_path = Path(path)
    if not file_path.is_file():
        return None
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _kernels(phi0, phi1, momentum):
    """Per-column real kernels. Sum with occupations matches column_moments."""
    momentum0 = momentum @ phi0
    momentum1 = momentum @ phi1
    kinetic = np.conjugate(phi0) * (-1j * momentum1) + np.conjugate(phi1) * (1j * momentum0)
    sigma = np.conjugate(phi0) * phi1 + np.conjugate(phi1) * phi0
    current = np.conjugate(phi0) * momentum0 + np.conjugate(phi1) * momentum1
    return kinetic.real, sigma.real, current.real


def matter_ledger(system, state):
    """Coordinate integrand, its Noether flux, and normal-observer stresses.

    Hamiltonian samples are summed directly. Stresses are densities: the nodal
    derivative is divided by dx. Normal shell energy FL/r already includes that
    factor through V ρ dx.
    """
    source = source_from_columns(system, state)
    kinetic, sigma, current = _kernels(state.phi0, state.phi1, system.momentum)
    weights = np.asarray(system.occupations, dtype=float)
    if kinetic.shape[1] != weights.size:
        raise ValueError("column count and occupations differ")
    multiplicity = float(system.multiplicity)
    length = system.length_density
    radial = state.Q
    shift = system.shift
    conformal = length / radial
    mass = system.kappa * length
    weighted_k = kinetic @ weights
    weighted_s = sigma @ weights
    weighted_j = current @ weights
    # Owner kernels, used for every summed balance. The column sum is a check.
    kernel_gap = {
        "K": float(np.max(np.abs(weighted_k - source["K"]))),
        "S": float(np.max(np.abs(weighted_s - np.real(source["S1"])))),
        "J": float(np.max(np.abs(weighted_j - source["Pmom"]))),
    }
    lapse_weighted = length * source["force_L"]
    hamiltonian = lapse_weighted + shift * source["force_beta"]
    # Fβ = −M J, so β Fβ = −M β J and e = L FL + β Fβ.
    flux_shift = -shift * hamiltonian
    flux_proper = multiplicity * (conformal ** 2) * source["Pmom"]
    flux_shift_pressure = shift * radial * source["force_Q"]
    stresses = observer_stresses(system, state, source)
    physical = physical_forces_from_conformal(
        source["force_L"], source["force_Q"], source["force_beta"],
        length, radial, state.r,
    )
    # Chain rule for Hamiltonian derivatives. Not a stress identity.
    ward = (
        (state.r * length) * physical["N"]
        + (state.r * radial) * physical["q"]
        + state.r * physical["r"]
    )
    columns = []
    for index, weight in enumerate(weights):
        column_e = multiplicity * weight * (
            conformal * kinetic[:, index] + mass * sigma[:, index] - shift * current[:, index]
        )
        columns.append({
            "occupation": float(weight),
            "shift": -shift * column_e,
            "proper": multiplicity * weight * (conformal ** 2) * current[:, index],
            "shift_pressure": -multiplicity * weight * conformal * shift * kinetic[:, index],
            "lapse_weighted": multiplicity * weight * (
                conformal * kinetic[:, index] + mass * sigma[:, index]
            ),
            "hamiltonian": column_e,
        })
    return {
        "source": source,
        "kernel_gap": kernel_gap,
        "lapse_weighted_integrand": lapse_weighted,
        "normal_energy_nodal": stresses["normal_energy_nodal"],
        "hamiltonian_integrand": hamiltonian,
        "flux_shift": flux_shift,
        "flux_proper": flux_proper,
        "flux_shift_pressure": flux_shift_pressure,
        "flux": flux_shift + flux_proper + flux_shift_pressure,
        "coordinate_conjugate_Q": source["force_Q"],
        "stresses": stresses,
        "hamiltonian_chain_ward": ward,
        "columns": columns,
        "multiplicity": multiplicity,
    }


def observer_stresses(system, state, source=None):
    """Proper-volume densities from nodal Hamiltonian derivatives.

    rho, p_r, p_perp, j include division by dx. normal_energy_nodal = FL/r
    equals V rho dx and is the shell energy sample.
    """
    if source is None:
        source = source_from_columns(system, state)
    force_l = source["force_L"]
    force_q = source["force_Q"]
    force_beta = source["force_beta"]
    radius = state.r
    radial = state.Q
    length = system.length_density
    spacing = system.dx
    sphere = 4.0 * np.pi * radius ** 4
    normal = force_l / radius
    volume = 4.0 * np.pi * (radius * radial) * radius ** 2
    return {
        "normal_energy_nodal": normal,
        "volume_per_dx": volume,
        "rho": force_l / (sphere * radial * spacing),
        "p_r": -force_q / (sphere * length * spacing),
        "p_perp": (length * force_l + radial * force_q) / (
            2.0 * sphere * length * radial * spacing
        ),
        "j": -force_beta / (sphere * radial ** 2 * spacing),
        "shell_gap": normal - volume * (force_l / (sphere * radial * spacing)) * spacing,
    }


def _geometry_pieces(system, state):
    derivative = system.derivative
    radius_x = derivative @ state.r
    chi_x = derivative @ state.chi
    momentum_x = derivative @ state.p_Q
    lapse_x = derivative @ system.length_density
    force_r, force_chi = partial_F(state.r, system.A, system.C_W)
    force_chi = float(force_chi)
    curvature = feedback_F(state.r, state.chi, system.A, system.C_W)
    curvature_x = derivative @ curvature
    product_gap = curvature_x - (force_r * radius_x + force_chi * chi_x)
    vertical = float(feedback_Z(system.A))
    potential = feedback_V(state.r, state.chi, system.A, system.C_W, system.C_F, system.flux)
    pi = state.p_r - (force_r / force_chi) * state.p_chi
    kinetic = (
        state.p_Q * state.p_chi / (2 * force_chi)
        + pi ** 2 / (4 * vertical * state.Q)
        + vertical * radius_x ** 2 / state.Q
        - state.Q * potential
    )
    density = (
        system.length_density * kinetic
        + 2 * curvature_x * lapse_x / state.Q
        + system.shift * (
            state.p_r * radius_x
            + state.p_chi * chi_x
            - state.Q * momentum_x
        )
    )
    radial_coeff = (
        2 * vertical * system.length_density * radius_x / state.Q
        + 2 * force_r * lapse_x / state.Q
        + system.shift * state.p_r
    )
    angular_coeff = 2 * force_chi * lapse_x / state.Q + system.shift * state.p_chi
    momentum_coeff = -state.Q * system.shift
    return {
        "density": density,
        "radius_x": radius_x,
        "chi_x": chi_x,
        "momentum_x": momentum_x,
        "curvature_x": curvature_x,
        "product_gap": product_gap,
        "radial_coeff": radial_coeff,
        "angular_coeff": angular_coeff,
        "momentum_coeff": momentum_coeff,
        "vertical": vertical,
        "force_chi": force_chi,
    }


def geometric_ledger(system, state, rate):
    """Canonical flux of the SBP density along an existing rate.

    Shift transport is the flux evaluated on the shift pieces of the velocities.
    The remainder is the proper-plus-force flux. Pressure work is not put into
    either flux; it is force_Q * Qdot, cancelled against the matter work.
    """
    pieces = _geometry_pieces(system, state)
    shift = system.shift
    shift_radius = shift * pieces["radius_x"]
    shift_chi = shift * pieces["chi_x"]
    shift_momentum = shift * pieces["momentum_x"]
    flux_shift = (
        pieces["radial_coeff"] * shift_radius
        + pieces["angular_coeff"] * shift_chi
        + pieces["momentum_coeff"] * shift_momentum
    )
    flux = (
        pieces["radial_coeff"] * rate.r
        + pieces["angular_coeff"] * rate.chi
        + pieces["momentum_coeff"] * rate.p_Q
    )
    return {
        "density": pieces["density"],
        "nodal_energy": system.dx * pieces["density"],
        "flux_density": flux,
        "nodal_flux": system.dx * flux,
        "shift_flux_density": flux_shift,
        "nodal_shift_flux": system.dx * flux_shift,
        "proper_flux_density": flux - flux_shift,
        "nodal_proper_flux": system.dx * (flux - flux_shift),
        "product_gap": pieces["product_gap"],
        "quasilocal_density": 2 * system.length_density * pieces["curvature_x"] / state.Q,
    }


def coordinate_metric_work(ledger, rate):
    """FQ Q̇. Includes Q̇ = ∂x(β Q). Not the observer pressure work."""
    return ledger["coordinate_conjugate_Q"] * rate.Q


def proper_balance_terms(system, state, rate, ledger=None):
    """Normal-observer flux, expansion, and pressure work along one rate.

    K_r and K_perp use ∂x(β q) and β r_x from the same derivative as the state.
    """
    if ledger is None:
        ledger = matter_ledger(system, state)
    stresses = ledger["stresses"]
    radius = state.r
    radial = state.Q
    length = system.length_density
    shift = system.shift
    lapse = radius * length
    radial_metric = radius * radial
    derivative = system.derivative
    radius_x = derivative @ radius
    expansion_r = (
        (rate.r * radial + radius * rate.Q) - derivative @ (shift * radial_metric)
    ) / (lapse * radial_metric)
    expansion_perp = (rate.r - shift * radius_x) / (lapse * radius)
    volume = stresses["volume_per_dx"]
    flux = volume * (
        lapse * stresses["j"] / radial_metric - shift * stresses["rho"]
    ) * system.dx
    pressure_work = -lapse * volume * (
        stresses["p_r"] * expansion_r + 2.0 * stresses["p_perp"] * expansion_perp
    ) * system.dx
    lapse_gradient_work = -volume * (stresses["j"] / radial_metric) * (derivative @ lapse) * system.dx
    return {
        "flux_nodal": flux,
        "K_r": expansion_r,
        "K_perp": expansion_perp,
        "proper_pressure_work": pressure_work,
        "momentum_lapse_work": lapse_gradient_work,
        "source_nodal": pressure_work + lapse_gradient_work,
    }


def proper_pointwise_slope(system, state, rate, step=1e-6):
    def energy(sample):
        return matter_ledger(system, sample)["normal_energy_nodal"]

    forward = energy(_combine(state, rate, step))
    backward = energy(_combine(state, rate, -step))
    return (forward - backward) / (2 * step)


def proper_pointwise_residual(system, state, rate, step=1e-6):
    ledger = matter_ledger(system, state)
    terms = proper_balance_terms(system, state, rate, ledger)
    slope = proper_pointwise_slope(system, state, rate, step)
    return slope + (system.derivative @ terms["flux_nodal"]) - terms["source_nodal"], terms


def fixed_window(length, count, centers, half_width):
    """Smooth partition bumps on a uniform periodic grid. Sum is identically one."""
    coordinate = np.arange(count, dtype=float) * (float(length) / int(count))
    windows = []
    for center in centers:
        delta = (coordinate - float(center) + 0.5 * length) % length - 0.5 * length
        values = np.zeros(count, dtype=float)
        mask = np.abs(delta) < half_width
        values[mask] = 0.5 * (1.0 + np.cos(np.pi * delta[mask] / half_width))
        windows.append(values)
    total = np.sum(windows, axis=0)
    return coordinate, windows, total


def _smear(weights, values):
    return float(np.sum(np.asarray(weights) * np.asarray(values)))


def matter_window_prediction(system, state, rate, window):
    ledger = matter_ledger(system, state)
    derivative = system.derivative @ window
    work = coordinate_metric_work(ledger, rate)
    channels = {
        "shift_transport": _smear(derivative, ledger["flux_shift"]),
        "proper_normal_flux": _smear(derivative, ledger["flux_proper"]),
        "shift_pressure_cross": _smear(derivative, ledger["flux_shift_pressure"]),
        "coordinate_metric_work": _smear(window, work),
    }
    channels["prediction"] = float(sum(channels.values()))
    channels["mode_columns"] = [
        {
            "occupation": column["occupation"],
            "shift_transport": _smear(derivative, column["shift"]),
            "proper_normal_flux": _smear(derivative, column["proper"]),
            "shift_pressure_cross": _smear(derivative, column["shift_pressure"]),
            "lapse_weighted_energy": _smear(window, column["lapse_weighted"]),
            "hamiltonian_energy": _smear(window, column["hamiltonian"]),
        }
        for column in ledger["columns"]
    ]
    return ledger, channels


def matter_pointwise_slope(system, state, rate, step=1e-6):
    def energy(sample):
        return matter_ledger(system, sample)["hamiltonian_integrand"]

    forward = energy(_combine(state, rate, step))
    backward = energy(_combine(state, rate, -step))
    return (forward - backward) / (2 * step)


def matter_window_slope(system, state, rate, window, step=1e-6):
    return float(np.sum(window * matter_pointwise_slope(system, state, rate, step)))


def transported_matter_slope(system, state, rate, window, step=1e-6):
    """Window Lie-dragged by ∂t W = β ∂x W, the shift flow of the matter flux."""
    motion = system.shift * (system.derivative @ window)

    def energy(sample):
        return matter_ledger(system, sample)["hamiltonian_integrand"]

    forward = energy(_combine(state, rate, step))
    backward = energy(_combine(state, rate, -step))
    return float(
        np.sum((window + step * motion) * forward - (window - step * motion) * backward) / (2 * step)
    )


def total_window_prediction(system, state, rate, window):
    matter, channels = matter_window_prediction(system, state, rate, window)
    geometry = geometric_ledger(system, state, rate)
    derivative = system.derivative @ window
    flux = matter["flux"] - geometry["nodal_flux"]
    prediction = _smear(derivative, flux)
    measured = geometry["nodal_energy"] + matter["hamiltonian_integrand"]
    constraint = _constraint_density(system, state, matter["source"])
    quasilocal = system.dx * geometry["quasilocal_density"]
    return {
        "matter_channels": channels,
        "measured_energy": _smear(window, measured),
        "matter_energy": _smear(window, matter["hamiltonian_integrand"]),
        "lapse_weighted_energy": _smear(window, matter["lapse_weighted_integrand"]),
        "normal_energy": _smear(window, matter["normal_energy_nodal"]),
        "gravity_energy": _smear(window, geometry["nodal_energy"]),
        "constraint_smearing": _smear(window, system.dx * constraint),
        "quasilocal_boundary": -_smear(derivative, quasilocal),
        "prediction": prediction,
        "geometric_shift_flux": _smear(derivative, geometry["nodal_shift_flux"]),
        "geometric_proper_flux": _smear(derivative, geometry["nodal_proper_flux"]),
        "coordinate_metric_work": channels["coordinate_metric_work"],
        "product_gap_max": float(np.max(np.abs(geometry["product_gap"]))),
        "ward_max": float(np.max(np.abs(matter["hamiltonian_chain_ward"]))),
        "kernel_gap": matter["kernel_gap"],
        "shell_gap_max": float(np.max(np.abs(matter["stresses"]["shell_gap"]))),
    }


def total_pointwise_slope(system, state, rate, step=1e-6):
    def energy(sample):
        geometry = _geometry_pieces(system, sample)["density"]
        matter = matter_ledger(system, sample)["hamiltonian_integrand"]
        return system.dx * geometry + matter

    forward = energy(_combine(state, rate, step))
    backward = energy(_combine(state, rate, -step))
    return (forward - backward) / (2 * step)


def total_window_slope(system, state, rate, window, step=1e-6):
    return float(np.sum(window * total_pointwise_slope(system, state, rate, step)))


def _constraint_density(system, state, source):
    return (
        system.length_density * (hamilton_constraint(system, state) + source["force_L"] / system.dx)
        + system.shift * (shift_constraint(system, state) + source["force_beta"] / system.dx)
    )


def band_limited_state(points=64):
    """Low harmonic geometry and AP spinors. Independent of the physical packet."""
    system, state = build_system(points)
    state = state.copy()
    coordinate = system.xi
    length = system.length
    state.Q = state.Q * (1.0 + 0.04 * np.cos(2 * np.pi * coordinate / length))
    state.r = 4.0 + 0.1 * np.cos(2 * np.pi * coordinate / length)
    state.chi = 0.02 * np.sin(2 * np.pi * coordinate / length)
    state.p_Q = 0.03 * np.cos(2 * np.pi * coordinate / length)
    state.p_r = 0.02 * np.sin(2 * np.pi * coordinate / length)
    state.p_chi = 0.01 * np.cos(4 * np.pi * coordinate / length)
    wavenumbers = np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]) * (2 * np.pi / length)
    phi0 = np.zeros((points, 6), dtype=complex)
    phi1 = np.zeros((points, 6), dtype=complex)
    envelope = 0.6 + 0.4 * np.cos(2 * np.pi * coordinate / length)
    for index, wavenumber in enumerate(wavenumbers):
        phase = np.exp(1j * wavenumber * coordinate) / np.sqrt(points)
        phi0[:, index] = envelope * phase / np.sqrt(2)
        phi1[:, index] = (1j if index % 2 == 0 else -1j) * envelope * phase / np.sqrt(2)
        norm = np.sqrt(np.sum(np.abs(phi0[:, index]) ** 2 + np.abs(phi1[:, index]) ** 2))
        phi0[:, index] /= norm
        phi1[:, index] /= norm
    state.phi0 = phi0
    state.phi1 = phi1
    return system, state


def packet_alias_residual(points):
    """Physical packet, geometry frozen. Residual is matter-flux aliasing."""
    system, state = build_system(points)
    source = source_from_columns(system, state)
    rate = CauchyRate(
        np.zeros(points), np.zeros(points), np.zeros(points),
        np.zeros(points), np.zeros(points), np.zeros(points),
        -1j * source["image0"], -1j * source["image1"],
        source["force_L"], source["force_Q"], source["force_beta"], 0.0,
    )
    ledger = matter_ledger(system, state)
    step = 1e-7

    def energy(sample):
        return matter_ledger(system, sample)["hamiltonian_integrand"]

    forward = energy(_combine(state, rate, step))
    backward = energy(_combine(state, rate, -step))
    slope = (forward - backward) / (2 * step)
    residual = slope + (system.derivative @ ledger["flux"])
    return {
        "points": int(points),
        "max_abs": float(np.max(np.abs(residual))),
        "l2": float(np.linalg.norm(residual)),
    }


def smooth_controls(points=64):
    """Band-limited closure. Each omitted channel must account for its own size."""
    system, state = band_limited_state(points)
    rate = rates(system, state, include_matter_force=True)
    _coordinate, windows, partition = fixed_window(system.length, system.points, (1.0, 3.0, 5.0, 7.0), 2.0)
    window = windows[0]
    ledger, channels = matter_window_prediction(system, state, rate, window)
    slope = matter_window_slope(system, state, rate, window)
    closure = slope - channels["prediction"]
    normal_only = rate.Q - (system.derivative @ (system.shift * state.Q))
    work_normal = ledger["coordinate_conjugate_Q"] * normal_only
    without_shift_in_work = (
        channels["shift_transport"]
        + channels["proper_normal_flux"]
        + channels["shift_pressure_cross"]
        + _smear(window, work_normal)
    )
    transported = transported_matter_slope(system, state, rate, window)
    transported_prediction = (
        channels["proper_normal_flux"]
        + channels["shift_pressure_cross"]
        + channels["coordinate_metric_work"]
    )
    proper_residual, proper_terms = proper_pointwise_residual(system, state, rate)
    zero_momenta = state.copy()
    zero_momenta.p_Q = np.zeros_like(state.p_Q)
    zero_momenta.p_r = np.zeros_like(state.p_r)
    zero_momenta.p_chi = np.zeros_like(state.p_chi)
    zero_rate = rates(system, zero_momenta, include_matter_force=True)
    zero_terms = proper_balance_terms(system, zero_momenta, zero_rate)
    total_prediction = total_window_prediction(system, state, rate, window)
    total_slope = total_window_slope(system, state, rate, window)
    measured = geometric_ledger(system, state, rate)["nodal_energy"] + ledger["hamiltonian_integrand"]
    partition_gap = float(np.max(np.abs(partition - 1.0)))
    regional_sum = 0.0
    for piece in windows:
        regional_sum += _smear(piece, measured)
    gravity_gap = abs(gravity_energy(system, state) - float(np.sum(system.dx * _geometry_pieces(system, state)["density"])))
    field_gap = abs(field_energy(system, state) - float(np.sum(ledger["hamiltonian_integrand"])))
    global_gap = abs(total_energy(system, state) - float(np.sum(measured)))
    omitted = {
        "shift_pressure_cross": abs(closure + channels["shift_pressure_cross"]),
        "shift_transport": abs(closure + channels["shift_transport"]),
        "coordinate_metric_work": abs(closure + channels["coordinate_metric_work"]),
        "shift_piece_of_Qdot": abs(slope - without_shift_in_work),
    }
    return {
        "points": int(points),
        "matter_closure_error": float(closure),
        "transported_closure_error": float(transported - transported_prediction),
        "total_closure_error": float(total_slope - total_prediction["prediction"]),
        "proper_closure_max": float(np.max(np.abs(proper_residual))),
        "shell_gap_max": total_prediction["shell_gap_max"],
        "zero_expansion": {
            "K_r_max": float(np.max(np.abs(zero_terms["K_r"]))),
            "K_perp_max": float(np.max(np.abs(zero_terms["K_perp"]))),
            "coordinate_work_l1": float(np.sum(np.abs(coordinate_metric_work(
                matter_ledger(system, zero_momenta), zero_rate
            )))),
            "proper_pressure_work_l1": float(np.sum(np.abs(zero_terms["proper_pressure_work"]))),
        },
        "channels": {key: float(channels[key]) for key in (
            "shift_transport", "proper_normal_flux", "shift_pressure_cross", "coordinate_metric_work",
        )},
        "omitted_channel_error": omitted,
        "gravity_sum_error": float(gravity_gap),
        "field_sum_error": float(field_gap),
        "global_sum_error": float(global_gap),
        "partition_gap": partition_gap,
        "partition_energy_gap": float(regional_sum - total_energy(system, state)),
        "quasilocal_gap": float(
            total_prediction["measured_energy"]
            - total_prediction["constraint_smearing"]
            - total_prediction["quasilocal_boundary"]
        ),
        "product_gap_max": total_prediction["product_gap_max"],
        "ward_max": total_prediction["ward_max"],
        "kernel_gap_max": float(max(ledger["kernel_gap"].values())),
        "geometric_shift_flux": total_prediction["geometric_shift_flux"],
        "geometric_proper_flux": total_prediction["geometric_proper_flux"],
        "matter_energy": total_prediction["matter_energy"],
        "gravity_energy": total_prediction["gravity_energy"],
        "measured_energy": total_prediction["measured_energy"],
        "constraint_smearing": total_prediction["constraint_smearing"],
        "global_energy": float(total_energy(system, state)),
    }


def galerkin_leakage(fermions=32, quadrature=128):
    """Unprojected fine identity versus the lifted rate the stepper actually uses."""
    grid = build_grid(fermions, quadrature=quadrature)
    phi0, phi1 = manufactured_columns(fermions, seed=3)
    state = blank_state(grid, phi0, phi1)
    coordinate = grid.xi_g
    state.r = 4.0 + 0.05 * np.cos(2 * np.pi * coordinate / grid.length)
    state.Q = state.Q * (1.0 + 0.02 * np.cos(2 * np.pi * coordinate / grid.length))
    coarse, bundle = compose_fine_hamiltonian(grid, state, include_matter_force=True)
    fine = bundle["fine_state"]
    source = bundle["source"]
    geometric = geometric_rates(grid.fine, fine)
    unprojected_momentum = geometric[3] - source["force_Q"] / grid.dx_q
    unprojected = CauchyRate(
        geometric[0], geometric[1], geometric[2], unprojected_momentum,
        geometric[4], geometric[5], -1j * source["image0"], -1j * source["image1"],
        source["force_L"], source["force_Q"], source["force_beta"], 0.0,
    )
    lifted = CauchyRate(
        prolong_geometry(grid, coarse.Q),
        prolong_geometry(grid, coarse.r),
        prolong_geometry(grid, coarse.chi),
        prolong_geometry(grid, coarse.p_Q),
        prolong_geometry(grid, coarse.p_r),
        prolong_geometry(grid, coarse.p_chi),
        grid.U_f @ coarse.phi0,
        grid.U_f @ coarse.phi1,
        source["force_L"], source["force_Q"], source["force_beta"], 0.0,
    )
    window = 0.5 * (1.0 + np.cos(2 * np.pi * grid.fine.xi / grid.length))
    prediction = total_window_prediction(grid.fine, fine, unprojected, window)
    unprojected_slope = total_window_slope(grid.fine, fine, unprojected, window)
    lifted_slope = total_window_slope(grid.fine, fine, lifted, window)
    unprojected_proper, _ = proper_pointwise_residual(grid.fine, fine, unprojected)
    lifted_proper, _ = proper_pointwise_residual(grid.fine, fine, lifted)
    return {
        "nf": int(fermions),
        "nq": int(quadrature),
        "unprojected_closure_error": float(unprojected_slope - prediction["prediction"]),
        "lifted_closure_error": float(lifted_slope - prediction["prediction"]),
        "projection_leakage": float(lifted_slope - unprojected_slope),
        "proper_unprojected_max": float(np.max(np.abs(unprojected_proper))),
        "proper_lifted_max": float(np.max(np.abs(lifted_proper))),
        "proper_projection_leakage": float(np.max(np.abs(lifted_proper - unprojected_proper))),
        "qdot_difference_max": float(np.max(np.abs(unprojected.Q - lifted.Q))),
        "phi_difference_max": float(np.max(np.abs(unprojected.phi0 - lifted.phi0))),
        "uses_manufactured_columns": True,
        "frozen_link_copied": False,
    }


def primitives_pass(smooth, leakage, alias_coarse, alias_fine):
    channel_errors = smooth["omitted_channel_error"]
    channels_distinct = all(
        channel_errors[name] > max(1e-3, 0.5 * abs(smooth["channels"][name]))
        for name in ("shift_pressure_cross", "shift_transport", "coordinate_metric_work")
    )
    zero = smooth["zero_expansion"]
    return bool(
        abs(smooth["matter_closure_error"]) < TOL_SMOOTH_WINDOW
        and abs(smooth["transported_closure_error"]) < TOL_SMOOTH_WINDOW
        and abs(smooth["total_closure_error"]) < TOL_SMOOTH_WINDOW
        and smooth["proper_closure_max"] < TOL_PROPER
        and smooth["shell_gap_max"] < TOL_SUM
        and zero["K_r_max"] < 1e-5
        and zero["K_perp_max"] < 1e-8
        and zero["coordinate_work_l1"] > 100.0 * max(zero["proper_pressure_work_l1"], 1e-12)
        and abs(smooth["global_sum_error"]) < TOL_SUM
        and abs(smooth["partition_energy_gap"]) < TOL_SUM
        and abs(smooth["field_sum_error"]) < TOL_SUM
        and abs(smooth["gravity_sum_error"]) < TOL_SUM
        and abs(smooth["quasilocal_gap"]) < TOL_SMOOTH_WINDOW
        and smooth["ward_max"] < TOL_WARD
        and smooth["kernel_gap_max"] < TOL_SUM
        and channels_distinct
        and channel_errors["shift_piece_of_Qdot"] > 1e-3
        and abs(leakage["unprojected_closure_error"]) < TOL_UNPROJECTED
        and leakage["proper_unprojected_max"] < TOL_PROPER
        and leakage["proper_lifted_max"] > 10.0 * max(leakage["proper_unprojected_max"], 1e-12)
        and abs(leakage["projection_leakage"]) > 10 * max(abs(leakage["unprojected_closure_error"]), 1e-12)
        and alias_fine["max_abs"] < 0.5 * alias_coarse["max_abs"]
    )


def _window_record(system, state, rate, window, name, matter_slope=None, total_slope=None):
    if matter_slope is None:
        matter_slope = matter_window_slope(system, state, rate, window)
    _ledger, channels = matter_window_prediction(system, state, rate, window)
    total = total_window_prediction(system, state, rate, window)
    if total_slope is None:
        total_slope = total_window_slope(system, state, rate, window)
    transported = transported_matter_slope(system, state, rate, window)
    transported_prediction = (
        channels["proper_normal_flux"]
        + channels["shift_pressure_cross"]
        + channels["coordinate_metric_work"]
    )
    proper_terms = proper_balance_terms(system, state, rate, _ledger)
    return {
        "name": name,
        "matter_energy": total["matter_energy"],
        "lapse_weighted_energy": total["lapse_weighted_energy"],
        "normal_energy": total["normal_energy"],
        "gravity_energy": total["gravity_energy"],
        "measured_energy": total["measured_energy"],
        "constraint_smearing": total["constraint_smearing"],
        "quasilocal_boundary": total["quasilocal_boundary"],
        "matter_channels": {key: channels[key] for key in (
            "shift_transport", "proper_normal_flux", "shift_pressure_cross", "coordinate_metric_work",
        )},
        "proper_pressure_work": _smear(window, proper_terms["proper_pressure_work"]),
        "momentum_lapse_work": _smear(window, proper_terms["momentum_lapse_work"]),
        "K_r_max": float(np.max(np.abs(proper_terms["K_r"]))),
        "K_perp_max": float(np.max(np.abs(proper_terms["K_perp"]))),
        "mode_columns": channels["mode_columns"],
        "geometric_shift_flux": total["geometric_shift_flux"],
        "geometric_proper_flux": total["geometric_proper_flux"],
        "matter_closure_error": float(matter_slope - channels["prediction"]),
        "transported_matter_closure_error": float(transported - transported_prediction),
        "total_closure_error": float(total_slope - total["prediction"]),
        "ward_max": total["ward_max"],
        "product_gap_max": total["product_gap_max"],
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if np.isfinite(number) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def load_nf512_initial(path=None):
    """Read the saved diagnostic state. Does not rewrite the npz."""
    source = V5_NPZ if path is None else Path(path)
    with np.load(source, allow_pickle=False) as data:
        fermions = int(data["nf512_nf"])
        quadrature = int(data["nf512_nq"])
        length = float(data["nf512_length"])
        state = CauchyState(
            Q=np.array(data["nf512_geometry_Q"], dtype=float, copy=True),
            r=np.array(data["nf512_geometry_r"], dtype=float, copy=True),
            chi=np.array(data["nf512_geometry_chi"], dtype=float, copy=True),
            p_Q=np.array(data["nf512_geometry_p_Q"], dtype=float, copy=True),
            p_r=np.array(data["nf512_geometry_p_r"], dtype=float, copy=True),
            p_chi=np.array(data["nf512_geometry_p_chi"], dtype=float, copy=True),
            phi0=np.array(data["nf512_columns_phi0"], copy=True),
            phi1=np.array(data["nf512_columns_phi1"], copy=True),
        )
        recorded_radius = np.array(data["nf512_quadrature_radius"], dtype=float, copy=True)
    grid = build_grid(fermions, quadrature=quadrature, length=length)
    if state.r.shape != (grid.ng,) or state.phi0.shape[0] != grid.nf:
        raise ValueError("saved nf512 arrays do not match the owned grid")
    return grid, state, recorded_radius, _sha256(source)


def _fine_rates(grid, state):
    coarse, bundle = compose_fine_hamiltonian(grid, state, include_matter_force=True)
    fine = bundle["fine_state"]
    source = bundle["source"]
    geometric = geometric_rates(grid.fine, fine)
    unprojected = CauchyRate(
        geometric[0], geometric[1], geometric[2],
        geometric[3] - source["force_Q"] / grid.dx_q,
        geometric[4], geometric[5],
        -1j * source["image0"], -1j * source["image1"],
        source["force_L"], source["force_Q"], source["force_beta"],
        float(np.sum(source["force_Q"] * geometric[0])),
    )
    lifted = CauchyRate(
        prolong_geometry(grid, coarse.Q),
        prolong_geometry(grid, coarse.r),
        prolong_geometry(grid, coarse.chi),
        prolong_geometry(grid, coarse.p_Q),
        prolong_geometry(grid, coarse.p_r),
        prolong_geometry(grid, coarse.p_chi),
        grid.U_f @ coarse.phi0,
        grid.U_f @ coarse.phi1,
        source["force_L"], source["force_Q"], source["force_beta"],
        float(coarse.fieldwork_power),
    )
    return fine, unprojected, lifted


def _limits():
    return {
        "renewal": RENEWAL,
        "continuum_limit_claimed": CONTINUUM_LIMIT_CLAIMED,
        "frozen_B_embedding_claimed": FROZEN_B_EMBEDDING_CLAIMED,
        "old_circulation_transferred": OLD_CIRCULATION_TRANSFERRED,
        "vacuum_branch_included": VACUUM_BRANCH_INCLUDED,
        "accepted_trajectory": ACCEPTED_TRAJECTORY,
        "new_solver": False,
        "static_gauge": "L and beta are not evolved",
        "multiplicity": "M=4 kappa once, inside the existing nodal forces",
        "half_density": "canonical u = r sqrt(q) psi, sampled as existing ell2 columns",
        "observer_energy": "nodal FL/r = V rho dx, not L*FL",
        "coordinate_work": "FQ*Qdot, including shift drag of Q",
        "proper_pressure_work": "-N V (p_r K_r + 2 p_perp K_perp) dx",
    }


def run_diagnostic(path=None, steps=2, step_size=5e-4):
    """Primitives first. The saved nf512 state is evolved only if they pass.

    Two RK4 steps are a diagnostic of the regional ledger, not a renewal and
    not a rerun of the recorded T=0.005 window.
    """
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    smooth = smooth_controls(64)
    leakage = galerkin_leakage(32, 128)
    alias_coarse = packet_alias_residual(32)
    alias_fine = packet_alias_residual(64)
    passed = primitives_pass(smooth, leakage, alias_coarse, alias_fine)
    record = {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC",
        "correction": CORRECTION,
        "global_energy_meaning": GLOBAL_ENERGY_MEANING,
        "regional_energy_meaning": REGIONAL_ENERGY_MEANING,
        "coordinate_reference_retained": _coordinate_reference_gap(smooth),
        "limits": _limits(),
        "primitives": {
            "smooth": smooth,
            "galerkin_leakage": leakage,
            "packet_alias_n32": alias_coarse,
            "packet_alias_n64": alias_fine,
            "pass": passed,
        },
        "nf512": None,
        "evolution": "not_run",
        "v5_npz_sha256": None,
        "module_sha256": _sha256(_MODULE_PATH),
    }
    if not passed:
        record["evolution"] = "skipped_primitives_failed"
        _write(path, record)
        return record
    grid, state, recorded_radius, digest = load_nf512_initial()
    record["v5_npz_sha256"] = digest
    fine, unprojected, lifted = _fine_rates(grid, state)
    radius_gap = float(np.max(np.abs(fine.r - recorded_radius)))
    _coordinate, windows, partition = fixed_window(grid.length, grid.nq, (1.0, 3.0, 5.0, 7.0), 2.0)
    names = ("packet_0_2", "packet_2_4", "packet_4_6", "complement_6_8")

    def snapshot(sample, label):
        sample_fine, sample_unprojected, sample_lifted = _fine_rates(grid, sample)
        lifted_matter_slope = matter_pointwise_slope(grid.fine, sample_fine, sample_lifted)
        lifted_total_slope = total_pointwise_slope(grid.fine, sample_fine, sample_lifted)
        unprojected_total_slope = total_pointwise_slope(grid.fine, sample_fine, sample_unprojected)
        lifted_proper, _lifted_proper_terms = proper_pointwise_residual(
            grid.fine, sample_fine, sample_lifted
        )
        unprojected_proper, _unprojected_proper_terms = proper_pointwise_residual(
            grid.fine, sample_fine, sample_unprojected
        )
        rows = []
        unprojected_errors = []
        proper_errors = []
        unprojected_proper_errors = []
        for window, name in zip(windows, names):
            row = _window_record(
                grid.fine, sample_fine, sample_lifted, window, name,
                matter_slope=float(np.sum(window * lifted_matter_slope)),
                total_slope=float(np.sum(window * lifted_total_slope)),
            )
            prediction = total_window_prediction(grid.fine, sample_fine, sample_unprojected, window)
            error = float(np.sum(window * unprojected_total_slope) - prediction["prediction"])
            proper_error = float(np.sum(window * lifted_proper))
            unprojected_proper_error = float(np.sum(window * unprojected_proper))
            unprojected_errors.append(error)
            proper_errors.append(proper_error)
            unprojected_proper_errors.append(unprojected_proper_error)
            row["unprojected_total_closure_error"] = error
            row["proper_closure_error"] = proper_error
            row["unprojected_proper_closure_error"] = unprojected_proper_error
            rows.append(row)
        measured = float(sum(row["measured_energy"] for row in rows))
        matter = float(sum(row["matter_energy"] for row in rows))
        gravity = float(sum(row["gravity_energy"] for row in rows))
        normal = float(sum(row["normal_energy"] for row in rows))
        return {
            "label": label,
            "windows": rows,
            "partition_gap": float(np.max(np.abs(partition - 1.0))),
            "sum_measured_energy": measured,
            "sum_matter_energy": matter,
            "sum_normal_energy": normal,
            "sum_gravity_energy": gravity,
            "global_energy": float(total_energy(grid.fine, sample_fine)),
            "sum_minus_global": float(measured - total_energy(grid.fine, sample_fine)),
            "field_energy": float(field_energy(grid.fine, sample_fine)),
            "gravity_energy": float(gravity_energy(grid.fine, sample_fine)),
            "positive_r": bool(np.min(sample_fine.r) > 0.0),
            "positive_Q": bool(np.min(sample_fine.Q) > 0.0),
            "r_min": float(np.min(sample_fine.r)),
            "Q_min": float(np.min(sample_fine.Q)),
            "lifted_projection_gap_max": float(max(abs(row["total_closure_error"]) for row in rows)),
            "unprojected_closure_max": float(max(abs(value) for value in unprojected_errors)),
            "proper_projection_gap_max": float(max(abs(value) for value in proper_errors)),
            "unprojected_proper_closure_max": float(max(abs(value) for value in unprojected_proper_errors)),
        }

    initial = snapshot(state, "initial")
    current = state
    for _step in range(int(steps)):
        current = rk4_step(grid, current, float(step_size))
    final = snapshot(current, "after_short_steps")
    record["nf512"] = {
        "nf": int(grid.nf),
        "ng": int(grid.ng),
        "nq": int(grid.nq),
        "quadrature_radius_gap": radius_gap,
        "steps": int(steps),
        "dt": float(step_size),
        "duration": float(steps * step_size),
        "source_label": "saved v5 nf512 initial state",
        "full_T_0_005_not_rerun": True,
        "initial": initial,
        "final": final,
        "global_energy_drift": float(final["global_energy"] - initial["global_energy"]),
        "regional_measured_drift": [
            float(after["measured_energy"] - before["measured_energy"])
            for before, after in zip(initial["windows"], final["windows"])
        ],
    }
    record["evolution"] = "short_diagnostic"
    _write(path, record)
    return record


def _write(path, record):
    destination = DIAGNOSTIC_JSON if path is None else Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
    os.replace(temporary, destination)
    record["path"] = str(destination)


if __name__ == "__main__":
    result = run_diagnostic()
    print(json.dumps({
        "pass": result["primitives"]["pass"],
        "evolution": result["evolution"],
        "path": result.get("path"),
    }, indent=2))
