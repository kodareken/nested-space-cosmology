"""Neutral positive-frequency packet and its common-action source.

The vacuum geometry is the reference ingoing chart

    ds² = N² dt² − (dr + β dt)² − r² dΩ²,

with r_m = sqrt(P²) = sqrt(magnetic_r2) taken from nsc_rn_reference.
Mass is an independent parameter. The approved ratios mass/r_m are
1.002, 1.01 and 1.04. Vacuum N = 1 is that reference result. The source
evaluator still accepts a coupled lapse.

The production state is one normalized pure packet column, shape
(2, points, 1), with filling ν. Its coefficient matrix on the scattering
modes is kept, and the per-channel CAR is G^{1/2} C G^{1/2}. The shell
factor 4 multiplies the stress once and does not enter that CAR.

The dimensionless target is ε = G_N (4 ν E_packet) / M. Angular-integrated
forces stay ∂E/∂(N, β, q, r). F_N and F_β do not depend on the values of
N or β, and F_r/N is the bare angular force. Orthonormal densities divide
by 4π q r², with j = −F_β/(4π q r²) and normal-observer T_{01} = −j.

Modes are the future-horizon regular κ=1 solutions at E>0, continued onto
every node of the caller grid, including excision between the horizons.
No second grid is built and no node is resampled when radii, weights and
the derivative are supplied.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from . import nsc_rn_reference as reference
from .nsc_covariant_operator import SIGMA1, SIGMA2


SHELL_KAPPA = reference.KAPPA
SHELL_MULTIPLICITY = reference.MULTIPLICITY
ENERGY_TARGETS = (1.0e-3, 1.0e-2)
MASS_OVER_RM = (1.002, 1.01, 1.04)
PACKET_CENTER_OVER_RM = 12.0
PACKET_WIDTH_OVER_RM = 2.0
PACKET_OMEGA_TIMES_RM = 1.0
FOUR_PI = 4.0 * np.pi
# ψ_Pauli = U χ, U† σ₂ U = diag(1, −1), U† σ₁ U = σ₂.
PAULI_FROM_CHARACTERISTIC = np.array(
    [[1.0, 1.0], [1j, -1j]], dtype=complex
) / np.sqrt(2.0)
SHIFT_FORM_AREAL = "dr + beta dt"
SHIFT_FORM_DECLARED = "q_adm dx + beta_one_form dt"


def shell_multiplicity(kappa):
    """κ = 1 shell, four copies once. The factor is the reference multiplicity."""
    if isinstance(kappa, bool) or not isinstance(kappa, (int, np.integer)):
        raise ValueError("integer angular label required")
    if int(kappa) != SHELL_KAPPA:
        raise ValueError(
            "neutral source owns the κ=1 shell only; multiplicity 4 is not "
            "reapplied, and κ=0 is not a mode"
        )
    if SHELL_MULTIPLICITY != 4:
        raise RuntimeError("reference shell multiplicity is not 4")
    return SHELL_MULTIPLICITY


def characteristic_transform():
    """Return U such that U† σ₂ U = diag(1, −1)."""
    basis = PAULI_FROM_CHARACTERISTIC
    diagonal = basis.conj().T @ SIGMA2 @ basis
    mixing = basis.conj().T @ SIGMA1 @ basis
    if not np.allclose(diagonal, np.diag([1.0, -1.0]), atol=1e-12):
        raise RuntimeError("characteristic basis does not diagonalize sigma2")
    if not np.allclose(mixing, SIGMA2, atol=1e-12):
        raise RuntimeError("angular sigma1 did not map to sigma2")
    return basis.copy()


def radial_sbp(points, left, right):
    """Uniform diagonal-norm SBP derivative and quadrature. Nonperiodic."""
    if isinstance(points, bool) or not isinstance(points, (int, np.integer)):
        raise ValueError("point count must be an integer")
    points = int(points)
    left = float(left)
    right = float(right)
    if points < 8 or not np.isfinite(left) or not np.isfinite(right) or not left < right:
        raise ValueError("need at least eight ordered SBP nodes")
    grid = np.linspace(left, right, points)
    step = float(grid[1] - grid[0])
    weights = np.full(points, step)
    weights[0] = weights[-1] = step / 2.0
    derivative = np.zeros((points, points))
    for index in range(1, points - 1):
        derivative[index, index - 1] = -0.5 / step
        derivative[index, index + 1] = 0.5 / step
    derivative[0, 0] = -1.0 / step
    derivative[0, 1] = 1.0 / step
    derivative[-1, -2] = -1.0 / step
    derivative[-1, -1] = 1.0 / step
    weighted = weights[:, None] * derivative
    boundary = weighted + weighted.T
    expected = np.zeros((points, points))
    expected[0, 0] = -1.0
    expected[-1, -1] = 1.0
    closure = float(np.max(np.abs(boundary - expected)))
    return {
        "grid": grid,
        "step": step,
        "weights": weights,
        "derivative": derivative,
        "closure_residual": closure,
        "periodic": False,
    }


def weighted_gram(phi, weights):
    """Per-channel weighted Gram. Multiplicity is not included."""
    return reference.weighted_gram(phi, weights)


def action_reference(mass, *, A=None, C_F=None, flux=None):
    """Reference RN at this mass. A, C_F and flux come from the action.

    The default triple is the audited source action. Passing the PG triple
    (same A, integer flux 1, and the C_F that keeps P²) returns the same
    G_N and the same P². Mass is not rescaled as a charge-free horizon.
    """
    supplied = (A is not None, C_F is not None, flux is not None)
    if not any(supplied):
        A = reference.AUDITED_A
        C_F = reference.AUDITED_C_F
        flux = reference.AUDITED_FLUX
    elif not all(supplied):
        raise ValueError("A, C_F, and flux must be supplied together")
    model = reference.RNReference.from_action(A, C_F, flux, mass)
    if model.V4 != 0.0:
        raise ValueError("V4 is 0 on this branch")
    if model.horizons["r_plus"] is None or model.extremal:
        raise ValueError("packet branch needs a non-extremal magnetic RN hole")
    return model


def magnetic_radius(model):
    """r_m = sqrt(P²). P² is the reference magnetic_r2, not the mass."""
    return float(np.sqrt(model.magnetic_r2))


def packet_request(r_m, *, center=None, width=None, omega=None):
    """Default packet on the magnetic radius: center 12 r_m, width 2 r_m, ω r_m = 1."""
    r_m = float(r_m)
    if not np.isfinite(r_m) or r_m <= 0.0:
        raise ValueError("positive magnetic radius required")
    requested = {
        "r_m": r_m,
        "center": PACKET_CENTER_OVER_RM * r_m if center is None else float(center),
        "width": PACKET_WIDTH_OVER_RM * r_m if width is None else float(width),
        "omega": PACKET_OMEGA_TIMES_RM / r_m if omega is None else float(omega),
    }
    if min(requested["width"], requested["omega"], requested["center"]) <= 0.0:
        raise ValueError("packet center, width, and frequency must be positive")
    if not all(np.isfinite(value) for value in requested.values()):
        raise ValueError("packet parameters must be finite")
    return requested


def vacuum_chart(model, radius):
    """Sample the reference areal chart. q_ADM = 1, so the shift is dr + β dt.

    Conformal Q = q_ADM/r equals 1/r on this chart. That gauge is not
    substituted into the force derivatives.
    """
    radius = np.asarray(radius, dtype=float)
    if radius.ndim != 1 or radius.size < 1 or not np.isfinite(radius).all():
        raise ValueError("areal radius samples must be finite")
    horizons = model.horizons
    if horizons["r_minus"] is None or horizons["r_plus"] is None:
        raise ValueError("static modes need a non-extremal horizon pair")
    outer = float(horizons["r_plus"])
    inner = float(horizons["r_minus"])
    if np.min(radius) <= inner:
        raise ValueError("static modes stay strictly above the inner horizon")
    jets = model.pg_metric_jets(radius)
    outgoing, ingoing = model.characteristic_speeds(radius)
    beta_r = np.asarray(model.beta_derivative(radius), dtype=float)
    lapse = np.asarray(jets["N"], dtype=float)
    shift = np.asarray(jets["beta"], dtype=float)
    q_adm = np.ones_like(radius)
    return {
        "radius": radius.copy(),
        "lapse": lapse,
        "shift": shift,
        "q_adm": q_adm,
        "radial_density": q_adm,
        "conformal_Q": q_adm / radius,
        "redshift": np.asarray(jets["f"], dtype=float),
        "speed_plus": np.asarray(outgoing, dtype=float),
        "speed_minus": np.asarray(ingoing, dtype=float),
        "speed_plus_derivative": -beta_r,
        "speed_minus_derivative": -beta_r,
        "outer_horizon": outer,
        "inner_horizon": float(model.horizons["r_minus"]),
        "G_N": float(model.G_N),
        "P2": float(model.magnetic_r2),
        "r_m": magnetic_radius(model),
        "mass": float(model.mass),
        "magnetic_flux": float(model.flux),
        "magnetic_flux_fixed": True,
        "lapse_derived": True,
        "shift_form": SHIFT_FORM_AREAL,
        "gauge_fix_imposed": False,
        "reference": model,
    }


def _point_chart(model, radius):
    chart = vacuum_chart(model, np.asarray([float(radius)], dtype=float))
    return {name: (float(value[0]) if isinstance(value, np.ndarray) else value)
            for name, value in chart.items() if name != "reference"}


def _speeds_from_fields(lapse, shift, radial_density):
    """β is the one-form coefficient in q_ADM dx + β dt.

    On the areal chart q_ADM = 1 and this is dr + β dt. A different q_ADM
    is the caller's declared conversion; β^x = β/q_ADM is not used here.
    """
    return (lapse - shift) / radial_density, -(lapse + shift) / radial_density


def _validate_metric(lapse, shift, radial_density, radius, weights, derivative):
    fields = [
        np.asarray(lapse, dtype=float), np.asarray(shift, dtype=float),
        np.asarray(radial_density, dtype=float), np.asarray(radius, dtype=float),
        np.asarray(weights, dtype=float),
    ]
    if any(field.ndim != 1 or not np.isfinite(field).all() for field in fields):
        raise ValueError("metric samples and quadrature weights must be finite vectors")
    if len({field.shape for field in fields}) != 1:
        raise ValueError("metric samples and quadrature weights differ in length")
    lapse, shift, radial_density, radius, weights = fields
    if np.min(lapse) <= 0.0 or np.min(radial_density) <= 0.0 or np.min(radius) <= 0.0:
        raise ValueError("lapse, radial density, and areal radius must be positive")
    if np.min(weights) <= 0.0:
        raise ValueError("SBP quadrature weights must be positive")
    derivative = np.asarray(derivative, dtype=float)
    points = lapse.shape[0]
    if derivative.shape != (points, points) or not np.isfinite(derivative).all():
        raise ValueError("SBP derivative must be a finite square matrix on this grid")
    return lapse, shift, radial_density, radius, weights, derivative


def _hermitian_flux(left, right, derivative):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    return (np.conj(left) * (derivative @ right) - np.conj(derivative @ left) * right) / (2j)


def _angular_bilinear(left, right):
    return -1j * np.conj(left[0]) * right[1] + 1j * np.conj(left[1]) * right[0]


def mode_quadratic_form(mode, lapse, shift, radial_density, radius, weights, derivative, kappa=SHELL_KAPPA):
    """Killing matrix element of one characteristic spinor, per angular channel."""
    shell_multiplicity(kappa)
    lapse, shift, radial_density, radius, weights, derivative = _validate_metric(
        lapse, shift, radial_density, radius, weights, derivative
    )
    mode = np.asarray(mode, dtype=complex)
    if mode.shape != (2, lapse.shape[0]):
        raise ValueError("characteristic mode must have shape (2, points)")
    speed_plus, speed_minus = _speeds_from_fields(lapse, shift, radial_density)
    kinetic = (
        speed_plus * _hermitian_flux(mode[0], mode[0], derivative)
        + speed_minus * _hermitian_flux(mode[1], mode[1], derivative)
    )
    angular = (kappa * lapse / radius) * _angular_bilinear(mode, mode)
    density = np.real(kinetic + angular)
    return {
        "density": density,
        "energy": float(np.sum(weights * density)),
        "probability_current": speed_plus * np.abs(mode[0])**2 + speed_minus * np.abs(mode[1])**2,
    }


def _pair_integrand(left, right, lapse, shift, radial_density, radius, derivative, kappa):
    speed_plus, speed_minus = _speeds_from_fields(lapse, shift, radial_density)
    flux_plus = _hermitian_flux(left[0], right[0], derivative)
    flux_minus = _hermitian_flux(left[1], right[1], derivative)
    angular_scalar = _angular_bilinear(left, right)
    density = (
        speed_plus * flux_plus + speed_minus * flux_minus
        + (kappa * lapse / radius) * angular_scalar
    )
    return density, flux_plus, flux_minus, angular_scalar


def _rhs_characteristic(model, radius, plus, minus, energy):
    chart = _point_chart(model, radius)
    angular = chart["lapse"] * float(np.asarray(reference.angular_coupling(radius)))
    if abs(chart["speed_plus"]) < 1e-14 or abs(chart["speed_minus"]) < 1e-14:
        raise RuntimeError("characteristic speed vanished; step off the horizon before dividing")
    d_plus = (
        1j * energy * plus - angular * minus - 0.5 * chart["speed_plus_derivative"] * plus
    ) / chart["speed_plus"]
    d_minus = (
        1j * energy * minus + angular * plus - 0.5 * chart["speed_minus_derivative"] * minus
    ) / chart["speed_minus"]
    return d_plus, d_minus


def _horizon_jet(model, energy):
    """Value and radial derivative of the future-horizon regular spinor at r+.

    a+ vanishes at r+. Differentiating a+ χ+' = (iE − a+'/2) χ+ − α χ− and
    then setting a+ = 0 fixes χ+' without dividing by the outgoing speed.
    """
    shell_multiplicity(SHELL_KAPPA)
    energy = float(energy)
    if energy <= 0.0 or not np.isfinite(energy):
        raise ValueError("exterior Killing frequency must be positive")
    horizon = float(model.horizons["r_plus"])
    chart = _point_chart(model, horizon)
    if abs(chart["speed_plus"]) > 1e-8:
        raise RuntimeError("future horizon is not a root of the outgoing speed")
    angular = chart["lapse"] * float(np.asarray(reference.angular_coupling(horizon)))
    angular_derivative = -chart["lapse"] * float(SHELL_KAPPA) / horizon**2
    a_plus_derivative = chart["speed_plus_derivative"]
    a_minus = chart["speed_minus"]
    a_minus_derivative = chart["speed_minus_derivative"]
    if a_minus == 0.0 or a_plus_derivative == 0.0:
        raise ValueError("horizon characteristic slope vanished")
    minus = 1.0 + 0.0j
    denominator = 1j * energy - 0.5 * a_plus_derivative
    if denominator == 0.0:
        raise ValueError("horizon constraint is singular at this frequency")
    plus = angular * minus / denominator
    minus_derivative = (
        angular * plus + (1j * energy - 0.5 * a_minus_derivative) * minus
    ) / a_minus
    a_second = -float(model.beta_second_derivative(horizon))
    plus_denominator = 1.5 * a_plus_derivative - 1j * energy
    if plus_denominator == 0.0:
        raise ValueError("horizon derivative constraint is singular")
    plus_derivative = (
        -0.5 * a_second * plus
        - angular_derivative * minus
        - angular * minus_derivative
    ) / plus_denominator
    return {
        "radius": horizon,
        "plus": plus,
        "minus": minus,
        "plus_derivative": plus_derivative,
        "minus_derivative": minus_derivative,
        "plus_over_minus": plus / minus,
        "chart": chart,
    }


def future_horizon_initial(model, energy, radius=None):
    """Future-horizon regular root at r+. An off-horizon radius is refused."""
    jet = _horizon_jet(model, energy)
    if radius is not None and abs(float(radius) - jet["radius"]) > 1e-8 * jet["radius"]:
        raise ValueError("the regular root is imposed at the future horizon, not off it")
    return np.array([jet["plus"], jet["minus"]], dtype=complex), jet["chart"]


def _spinor_taylor(jet, radius):
    delta = float(radius) - jet["radius"]
    return (
        jet["plus"] + delta * jet["plus_derivative"],
        jet["minus"] + delta * jet["minus_derivative"],
    )


def _continuation_gap(horizon, inner, grid):
    gap = 1.0e-4 * horizon
    caps = [gap, 0.25 * (horizon - inner)]
    if grid[0] < horizon:
        caps.append(0.25 * (horizon - float(grid[0])))
    if grid[-1] > horizon:
        caps.append(0.25 * (float(grid[-1]) - horizon))
    gap = float(min(caps))
    if not np.isfinite(gap) or gap <= 0.0:
        raise ValueError("no room to step off the future horizon")
    return gap


def _integrate_branch(model, energy, jet, grid, indices, t_start, t_end, mode, rtol, atol):
    if indices.size == 0:
        return 0
    plus0, minus0 = _spinor_taylor(jet, t_start)

    def rhs(radius, state):
        plus = state[0] + 1j * state[1]
        minus = state[2] + 1j * state[3]
        d_plus, d_minus = _rhs_characteristic(model, radius, plus, minus, energy)
        return [float(np.real(d_plus)), float(np.imag(d_plus)),
                float(np.real(d_minus)), float(np.imag(d_minus))]

    solution = solve_ivp(
        rhs, (float(t_start), float(t_end)),
        [plus0.real, plus0.imag, minus0.real, minus0.imag],
        method="DOP853", rtol=rtol, atol=atol, dense_output=True,
    )
    if not solution.success or solution.sol is None:
        raise RuntimeError(solution.message)
    samples = solution.sol(grid[indices])
    mode[0, indices] = samples[0] + 1j * samples[1]
    mode[1, indices] = samples[2] + 1j * samples[3]
    return int(solution.nfev)


def scattering_mode(grid, energy, model, *, rtol=1.0e-7, atol=1.0e-9):
    """E > 0 mode from the future horizon onto every supplied node.

    The grid may start at an excision radius between the horizons. Nodes
    inside the outer horizon are kept. The integration steps off r+ and
    does not project with an absorbing matrix.
    """
    shell_multiplicity(SHELL_KAPPA)
    grid = np.asarray(grid, dtype=float)
    energy = float(energy)
    if grid.ndim != 1 or grid.size < 8 or np.any(np.diff(grid) <= 0.0):
        raise ValueError("increasing radial grid required")
    inner = float(model.horizons["r_minus"])
    if grid[0] <= inner:
        raise ValueError("grid must stay strictly above the inner horizon")
    jet = _horizon_jet(model, energy)
    horizon = jet["radius"]
    gap = _continuation_gap(horizon, inner, grid)
    mode = np.empty((2, grid.size), dtype=complex)
    near = np.abs(grid - horizon) <= gap
    for index in np.flatnonzero(near):
        plus, minus = _spinor_taylor(jet, grid[index])
        mode[0, index] = plus
        mode[1, index] = minus
    outward = np.flatnonzero(grid > horizon + gap)
    inward = np.flatnonzero(grid < horizon - gap)
    evaluations = _integrate_branch(
        model, energy, jet, grid, outward, horizon + gap, float(grid[-1]), mode, rtol, atol,
    )
    evaluations += _integrate_branch(
        model, energy, jet, grid, inward, horizon - gap, float(grid[0]), mode, rtol, atol,
    )
    filled = int(np.count_nonzero(near) + outward.size + inward.size)
    if filled != grid.size or not np.isfinite(mode).all():
        raise RuntimeError("radial nodes were dropped during horizon continuation")
    initial = np.array([jet["plus"], jet["minus"]], dtype=complex)
    return {
        "mode": mode,
        "energy": energy,
        "evaluations": int(evaluations),
        "horizon_initial": initial,
        "horizon_radius": horizon,
        "continuation_gap": gap,
        "interior_nodes": int(np.count_nonzero(grid < horizon)),
        "nodes_dropped": 0,
        "node_count": int(grid.size),
        "absorbing_matrix_projection": False,
    }


def _normalize_mode(mode, weights):
    norm = np.sqrt(np.real(np.sum(weights * np.sum(np.abs(mode)**2, axis=0))))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("scattering mode has vanishing weighted norm")
    return mode / norm, float(norm)


def _phase_align_ingoing(mode, grid, center):
    """Constant phase so the ingoing characteristic is real and positive at the center."""
    center_index = int(np.argmin(np.abs(grid - center)))
    reference_value = mode[1, center_index]
    if reference_value == 0.0:
        raise ValueError("ingoing characteristic vanishes at the packet center")
    return mode * np.exp(-1j * np.angle(reference_value)), center_index


def _gaussian_amplitudes(energies, omega, spectral_width):
    energies = np.asarray(energies, dtype=float)
    if np.any(energies <= 0.0):
        raise ValueError("Gaussian packet is supported on positive Killing frequencies only")
    if spectral_width <= 0.0:
        raise ValueError("positive spectral width required")
    return np.exp(-0.25 * ((energies - omega) / spectral_width) ** 2)


def select_packet_filling(packet_energy, mass, target, newton_constant, multiplicity=SHELL_MULTIPLICITY):
    """ν = ε M / (G_N 4 E_packet), so ε = G_N (4 ν E_packet) / M."""
    packet_energy = float(packet_energy)
    mass = float(mass)
    target = float(target)
    newton = float(newton_constant)
    if packet_energy <= 0.0 or mass <= 0.0 or newton <= 0.0:
        raise ValueError("positive packet energy, mass, and Newton constant required")
    if target not in ENERGY_TARGETS:
        raise ValueError("packet target ε must be 1e-3 or 1e-2")
    filling = target * mass / (newton * multiplicity * packet_energy)
    if not np.isfinite(filling) or filling <= 0.0 or filling >= 1.0:
        raise ValueError("packet filling left the open fermionic interval (0, 1)")
    actual = newton * multiplicity * filling * packet_energy / mass
    return float(filling), float(actual)


def radial_moments(profile, radius, weights):
    density = np.sum(np.abs(profile)**2, axis=0)
    total = float(np.sum(weights * density))
    if total <= 0.0:
        raise ValueError("packet density vanishes")
    center = float(np.sum(weights * radius * density) / total)
    variance = float(np.sum(weights * (radius - center) ** 2 * density) / total)
    return {"center": center, "width_rms": float(np.sqrt(max(variance, 0.0))), "norm": total}


def _gram_sqrt(gram):
    hermitian = 0.5 * (gram + gram.conj().T)
    eigenvalues, vectors = np.linalg.eigh(hermitian)
    if np.min(eigenvalues) <= 0.0:
        raise ValueError("mode Gram is not positive definite")
    return (vectors * np.sqrt(eigenvalues)) @ vectors.conj().T


def physical_car(gram, coefficient_matrix):
    """Per-channel CAR G^{1/2} C G^{1/2}. Multiplicity is not applied."""
    root = _gram_sqrt(gram)
    physical = root @ coefficient_matrix @ root
    return 0.5 * (physical + physical.conj().T)


def ode_residual(model, mode, grid, energy):
    """Relative RMS of the cleared characteristic ODE on the supplied nodes.

    The form a χ' − (right-hand side without division) stays finite at r+.
    Interior nodes are included. This is not a continuum spectral theorem.
    """
    mode = np.asarray(mode, dtype=complex)
    grid = np.asarray(grid, dtype=float)
    if mode.shape != (2, grid.size):
        raise ValueError("mode shape does not match the grid")
    chart = vacuum_chart(model, grid)
    energy = float(energy)
    angular = chart["lapse"] * np.asarray(reference.angular_coupling(grid), dtype=float)
    span = grid[2:] - grid[:-2]
    numerical = (mode[:, 2:] - mode[:, :-2]) / span
    sample = slice(1, -1)
    plus = mode[0, sample]
    minus = mode[1, sample]
    speed_plus = chart["speed_plus"][sample]
    speed_minus = chart["speed_minus"][sample]
    raw_plus = (
        (1j * energy - 0.5 * chart["speed_plus_derivative"][sample]) * plus
        - angular[sample] * minus
    )
    raw_minus = (
        angular[sample] * plus
        + (1j * energy - 0.5 * chart["speed_minus_derivative"][sample]) * minus
    )
    residual_plus = speed_plus * numerical[0] - raw_plus
    residual_minus = speed_minus * numerical[1] - raw_minus
    scale = max(abs(energy), 1.0)
    scale_plus = np.abs(raw_plus) + np.abs(speed_plus * plus) * scale + 1e-30
    scale_minus = np.abs(raw_minus) + np.abs(speed_minus * minus) * scale + 1e-30
    numerator = float(np.sum(span * (np.abs(residual_plus)**2 + np.abs(residual_minus)**2)))
    denominator = float(np.sum(span * (scale_plus**2 + scale_minus**2)))
    if denominator <= 0.0:
        raise ValueError("ODE residual scale vanished")
    return float(np.sqrt(numerator / denominator))


def _positive_frequency_window(omega, spectral_width, rank):
    """Symmetric samples that stay positive. Truncation below zero is recorded."""
    full_span = 2.0
    positive_span = 0.9 * float(omega) / float(spectral_width)
    span = min(full_span, positive_span)
    if span <= 0.0:
        raise ValueError("ingoing spectral window has no positive frequencies")
    offsets = np.linspace(-span, span, rank)
    energies = float(omega) + float(spectral_width) * offsets
    if np.min(energies) <= 0.0:
        raise ValueError("spectral window crossed zero frequency")
    return energies, span < full_span - 1e-12, float(span)


def packet_resolution_scales(model, *, rank=5):
    """Lengths for h <= min(λ_min/24, Δ/16, σ/24).

    σ is the requested packet width 2 r_m. λ_min is 2π min(|a+|, |a-|)/E_max
    at the packet center, using the positive-frequency window of this chart.
    """
    shell_multiplicity(SHELL_KAPPA)
    r_m = magnetic_radius(model)
    request = packet_request(r_m)
    chart = _point_chart(model, request["center"])
    sigma = float(request["width"])
    spectral_width = abs(float(chart["speed_minus"])) / sigma
    energies, truncated, span = _positive_frequency_window(request["omega"], spectral_width, int(rank))
    e_max = float(np.max(energies))
    speed = min(abs(float(chart["speed_plus"])), abs(float(chart["speed_minus"])))
    if e_max <= 0.0 or speed <= 0.0:
        raise ValueError("packet frequency window did not produce a wavelength")
    return {
        "lambda_min": float(2.0 * np.pi * speed / e_max),
        "sigma": sigma,
        "center": float(request["center"]),
        "width": sigma,
        "omega": float(request["omega"]),
        "e_max": e_max,
        "spectral_width": float(spectral_width),
        "positive_span_in_sigmas": span,
        "window_truncated_below_zero": bool(truncated),
        "wavelength": "2 pi min(|a+|, |a-|) / E_max at the packet center",
    }


def _packet_reference(mass_over_rm, A, C_F, flux, model):
    """Audited action, or the caller's RNReference when its G_N is that action."""
    mass_over_rm = float(mass_over_rm)
    if mass_over_rm not in MASS_OVER_RM:
        raise ValueError("mass/r_m must be 1.002, 1.01, or 1.04")
    if model is None:
        probe = action_reference(1.0, A=A, C_F=C_F, flux=flux)
        r_m = magnetic_radius(probe)
        return action_reference(mass_over_rm * r_m, A=A, C_F=C_F, flux=flux)
    if not isinstance(model, reference.RNReference):
        raise TypeError("model must be an RNReference")
    if A is not None and abs(float(A) - model.A) > 1e-12 * max(1.0, model.A):
        raise ValueError("A does not match the supplied reference")
    if C_F is not None and abs(float(C_F) - model.C_F) > 1e-12 * max(1.0, model.C_F):
        raise ValueError("C_F does not match the supplied reference")
    if flux is not None and abs(float(flux) - model.flux) > 1e-12 * max(1.0, abs(model.flux)):
        raise ValueError("flux does not match the supplied reference")
    if model.V4 != 0.0 or model.extremal or model.horizons["r_plus"] is None:
        raise ValueError("packet branch needs a non-extremal magnetic RN hole with V4 = 0")
    r_m = magnetic_radius(model)
    if abs(model.mass / r_m - mass_over_rm) > 1e-8:
        raise ValueError("supplied reference mass/r_m does not match mass_over_rm")
    return model


def _packet_grid(model, request, points, left, right, radii, weights, derivative):
    """Caller nodes, or one SBP grid from the midpoint excision. Never a second sample."""
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    excision = 0.5 * (r_minus + r_plus)
    if radii is None:
        if weights is not None or derivative is not None:
            raise ValueError("weights and derivative require caller radii")
        if points is None:
            points = 81
        if left is None:
            left = excision
        if right is None:
            right = request["center"] + 4.0 * request["width"]
        pack = radial_sbp(points, left, right)
        return pack["grid"], pack["weights"], pack["derivative"], "constructed_sbp"
    if weights is None or derivative is None:
        raise ValueError("caller radii require weights and derivative; refusing to resample")
    grid = np.array(np.asarray(radii, dtype=float), dtype=float, copy=True)
    if points is not None and int(points) != grid.size:
        raise ValueError("points does not match the caller grid; refusing to resample")
    if left is not None and abs(float(left) - float(grid[0])) > 1e-9 * max(1.0, abs(float(grid[0]))):
        raise ValueError("left does not match the caller grid; refusing to resample")
    if right is not None and abs(float(right) - float(grid[-1])) > 1e-9 * max(1.0, abs(float(grid[-1]))):
        raise ValueError("right does not match the caller grid; refusing to resample")
    weight_row = np.array(np.asarray(weights, dtype=float), dtype=float, copy=True)
    differencing = np.array(np.asarray(derivative, dtype=float), dtype=float, copy=True)
    if grid.ndim != 1 or grid.size < 8 or np.any(np.diff(grid) <= 0.0):
        raise ValueError("caller radii must be a strictly increasing grid of at least eight nodes")
    if weight_row.shape != grid.shape or np.any(weight_row <= 0.0) or not np.isfinite(weight_row).all():
        raise ValueError("caller weights must be positive and match the radii")
    if differencing.shape != (grid.size, grid.size) or not np.isfinite(differencing).all():
        raise ValueError("caller derivative must be a finite square matrix on these radii")
    return grid, weight_row, differencing, "caller"


def prepare_radial_packet(*, mass_over_rm=1.01, energy_target=1.0e-3, points=None, rank=5,
                          A=None, C_F=None, flux=None, model=None, center=None, width=None,
                          omega=None, left=None, right=None, radii=None, weights=None,
                          derivative=None, rtol=1.0e-7, atol=1.0e-9):
    """Pure ingoing packet column for nsc_rn_pg.make_state.

    phi has shape (2, points, 1) and occupations is [ν]. Spectral modes and
    the coefficient matrix C are retained beside that column. Pass the PG
    radii, weights and derivative to sample on that grid; nothing is resampled.
    The measured rms width is an indicator, not a fit to the requested width.
    """
    if isinstance(rank, bool) or not isinstance(rank, (int, np.integer)) or int(rank) < 3:
        raise ValueError("at least three positive-frequency samples are required")
    rank = int(rank)
    model = _packet_reference(mass_over_rm, A, C_F, flux, model)
    r_m = magnetic_radius(model)
    mass = float(model.mass)
    request = packet_request(r_m, center=center, width=width, omega=omega)
    grid, weights, derivative, binding = _packet_grid(
        model, request, points, left, right, radii, weights, derivative,
    )
    outer = float(model.horizons["r_plus"])
    inner = float(model.horizons["r_minus"])
    if not inner < float(grid[0]) < request["center"] < float(grid[-1]):
        raise ValueError("packet center must lie strictly inside a grid above the inner horizon")
    if request["center"] <= outer:
        raise ValueError("packet center must lie outside the outer horizon")
    chart = vacuum_chart(model, grid)
    ingoing_speed = float(np.interp(request["center"], grid, chart["speed_minus"]))
    if ingoing_speed >= 0.0:
        raise ValueError("packet center is not on the ingoing exterior characteristic")
    spectral_width = abs(ingoing_speed) / request["width"]
    requested_energies, truncated, span = _positive_frequency_window(
        request["omega"], spectral_width, rank,
    )
    modes = []
    horizon_ratios = []
    evaluations = 0
    for energy in requested_energies:
        solved = scattering_mode(grid, float(energy), model)
        normalized, _raw_norm = _normalize_mode(solved["mode"], weights)
        aligned, _center_index = _phase_align_ingoing(normalized, grid, request["center"])
        modes.append(aligned)
        horizon_ratios.append(solved["horizon_initial"][0] / solved["horizon_initial"][1])
        evaluations += solved["evaluations"]
    spectral_phi = np.stack(modes, axis=-1)
    gram = weighted_gram(spectral_phi, weights)
    per_mode_energy = np.array([
        mode_quadratic_form(
            spectral_phi[:, :, index], chart["lapse"], chart["shift"], chart["q_adm"],
            chart["radius"], weights, derivative,
        )["energy"]
        for index in range(rank)
    ])
    if np.min(per_mode_energy) <= 0.0:
        raise ValueError("measured scattering-mode Killing energy is not positive")
    amplitudes = _gaussian_amplitudes(requested_energies, request["omega"], spectral_width)
    superposition = spectral_phi @ amplitudes
    packet_norm = float(np.sqrt(np.real(np.sum(
        weights * np.sum(np.abs(superposition)**2, axis=0)
    ))))
    packet = superposition / packet_norm
    measured = mode_quadratic_form(
        packet, chart["lapse"], chart["shift"], chart["q_adm"],
        chart["radius"], weights, derivative,
    )
    filling, actual_target = select_packet_filling(
        measured["energy"], mass, energy_target, model.G_N, SHELL_MULTIPLICITY,
    )
    alpha = amplitudes / packet_norm
    coefficient_matrix = filling * np.outer(alpha, np.conj(alpha))
    car = physical_car(gram, coefficient_matrix)
    moments = radial_moments(packet, chart["radius"], weights)
    center_index = int(np.argmin(np.abs(grid - request["center"])))
    current = measured["probability_current"]
    tail_count = min(4, int(grid.size))
    production_phi = packet[:, :, None]
    production_occupations = np.array([filling])
    return {
        "schema": "NSC-RN-NEUTRAL-PACKET-v2",
        "phi": production_phi,
        "occupations": production_occupations,
        "packet_filling": filling,
        "occupation_matrix": np.array([[filling]], dtype=complex),
        "spectral_phi": spectral_phi,
        "spectral_coefficient_matrix": coefficient_matrix,
        "normalization": {
            "packet_weighted_norm": moments["norm"],
            "gram": gram,
            "production_rank": 1,
        },
        "car": {
            "per_channel": True,
            "multiplicity_in_car": False,
            "gram": gram,
            "coefficient_matrix": coefficient_matrix,
            "physical": car,
            "physical_eigenvalues": np.linalg.eigvalsh(car),
            "packet_filling": filling,
        },
        "killing_energy": SHELL_MULTIPLICITY * filling * measured["energy"],
        "sector_killing_energy": filling * measured["energy"],
        "packet_mode_energy": measured["energy"],
        "energy_target": float(energy_target),
        "actual_energy_target": actual_target,
        "energy_target_formula": "G_N * (4 * nu * E_packet) / M",
        "G_N": float(model.G_N),
        "angular_covariance": np.eye(SHELL_MULTIPLICITY) * filling,
        "spectrum": {
            "requested": requested_energies.copy(),
            "quadratic_form": per_mode_energy,
            "packet": measured["energy"],
            "spectral_width": spectral_width,
            "spectral_width_source": "abs(ingoing group velocity) / requested width",
            "ingoing_group_velocity": ingoing_speed,
            "window_truncated_below_zero": truncated,
            "positive_span_in_sigmas": span,
            "negative_frequency_contamination": 0.0,
        },
        "near_horizon_tail": {
            "radius": grid[:tail_count].copy(),
            "plus": spectral_phi[0, :tail_count, :].copy(),
            "minus": spectral_phi[1, :tail_count, :].copy(),
            "horizon_plus_over_minus": np.asarray(horizon_ratios),
            "packet_plus": packet[0, :tail_count].copy(),
            "packet_minus": packet[1, :tail_count].copy(),
        },
        "radial_moments": {
            **moments,
            "requested_center": request["center"],
            "requested_width": request["width"],
            "requested_omega": request["omega"],
            "width_is_fit": False,
        },
        "inward_flux": {
            "center_probability_current": float(current[center_index]),
            "integrated_negative_current": float(np.sum(weights * np.minimum(current, 0.0))),
            "inward_at_center": bool(current[center_index] < 0.0),
        },
        "background": {
            "mass": mass,
            "mass_over_rm": float(mass_over_rm),
            "r_m": r_m,
            "G_N": float(model.G_N),
            "P2": float(model.magnetic_r2),
            "magnetic_flux": float(model.flux),
            "magnetic_flux_fixed": True,
            "outer_horizon": chart["outer_horizon"],
            "inner_horizon": chart["inner_horizon"],
            "lapse_derived": True,
            "q_adm": chart["q_adm"],
            "conformal_Q": chart["conformal_Q"],
            "shift_form": SHIFT_FORM_AREAL,
            "gauge_fix_imposed": False,
            "grid": grid,
            "weights": weights,
            "derivative": derivative,
            "grid_binding": binding,
            "lapse": chart["lapse"],
            "shift": chart["shift"],
            "radial_density": chart["q_adm"],
            "radius": chart["radius"],
            "speed_plus": chart["speed_plus"],
            "speed_minus": chart["speed_minus"],
            "model": model,
        },
        "kappa": SHELL_KAPPA,
        "multiplicity": SHELL_MULTIPLICITY,
        "solver_evaluations": evaluations,
        "absorbing_matrix_projection": False,
        "periodic_identification": False,
        "pg_make_state": pg_make_state_arguments.__name__,
    }


def pg_make_state_arguments(prepared):
    """Arguments accepted by nsc_rn_pg.make_state: one packet column and [ν]."""
    phi = np.asarray(prepared["phi"], dtype=complex)
    occupations = np.asarray(prepared["occupations"], dtype=float)
    if phi.ndim != 3 or phi.shape[0] != 2 or phi.shape[2] != 1:
        raise ValueError("production phi must be the normalized packet column (2, points, 1)")
    if occupations.shape != (1,) or not np.isfinite(occupations[0]):
        raise ValueError("production occupation must be the scalar filling")
    if occupations[0] <= 0.0 or occupations[0] >= 1.0:
        raise ValueError("packet filling must lie in (0, 1)")
    return {"phi": phi, "occupations": occupations}


def mode_count_indicator(*, mass_over_rm=1.01, energy_target=1.0e-3, points=49, ranks=(3, 5)):
    """Cheap rank comparison. Widths are indicators, not a fit to 2 r_m."""
    rows = []
    for rank in ranks:
        prepared = prepare_radial_packet(
            mass_over_rm=mass_over_rm, energy_target=energy_target,
            points=points, rank=int(rank),
        )
        rows.append({
            "rank": int(rank),
            "center": prepared["radial_moments"]["center"],
            "width_rms": prepared["radial_moments"]["width_rms"],
            "requested_width": prepared["radial_moments"]["requested_width"],
            "inward_at_center": prepared["inward_flux"]["inward_at_center"],
            "center_probability_current": prepared["inward_flux"]["center_probability_current"],
            "actual_energy_target": prepared["actual_energy_target"],
        })
    return {"width_is_fit": False, "rows": rows}


def source_evaluate(phi, occupations, lapse, shift, radial_density, radius,
                    weights, derivative, *, kappa=SHELL_KAPPA,
                    occupation_matrix=None, shift_form=None):
    """Dirac stress. Positional arguments match nsc_rn_pg.matter_densities.

    `forces` are angular-integrated ∂E/∂(N, β, q, r), multiplicity 4 once.
    `physical` holds orthonormal ρ, j, p_r, p_Ω. j = −F_β/(4π q r²).
    β is the one-form in q dx + β dt. Pass shift_form only when q ≠ 1.
    """
    multiplicity = shell_multiplicity(kappa)
    lapse, shift, radial_density, radius, weights, derivative = _validate_metric(
        lapse, shift, radial_density, radius, weights, derivative
    )
    if shift_form is None:
        shift_form = SHIFT_FORM_AREAL if np.allclose(radial_density, 1.0) else SHIFT_FORM_DECLARED
    if shift_form not in (SHIFT_FORM_AREAL, SHIFT_FORM_DECLARED):
        raise ValueError("shift form must be the areal chart or the declared one-form conversion")
    if shift_form == SHIFT_FORM_AREAL and np.max(np.abs(radial_density - 1.0)) > 1e-12:
        raise ValueError("areal shift form dr+beta dt requires q_adm = 1")
    phi = np.asarray(phi, dtype=complex)
    if phi.ndim != 3 or phi.shape[0] != 2 or phi.shape[1] != lapse.shape[0]:
        raise ValueError("phi must have shape (2, points, rank)")
    rank = phi.shape[2]
    if occupation_matrix is None:
        occupations = np.asarray(occupations, dtype=float)
        if occupations.shape != (rank,) or not np.isfinite(occupations).all():
            raise ValueError("one finite occupation per mode required")
        if np.min(occupations) < 0.0 or np.max(occupations) >= 1.0:
            raise ValueError("diagonal occupations must lie in [0, 1)")
        coefficients = np.diag(occupations.astype(complex))
        coherent_rank_one = rank == 1
    else:
        coefficients = np.asarray(occupation_matrix, dtype=complex)
        if coefficients.shape != (rank, rank) or not np.isfinite(coefficients).all():
            raise ValueError("occupation matrix must be finite and rank by rank")
        if not np.allclose(coefficients, coefficients.conj().T, atol=1e-10):
            raise ValueError("occupation matrix must be Hermitian")
        spectrum = np.linalg.eigvalsh(0.5 * (coefficients + coefficients.conj().T))
        if np.min(spectrum) < -1e-10 or np.max(spectrum) >= 1.0:
            raise ValueError("occupation eigenvalues must lie in [0, 1)")
        coherent_rank_one = False
    points = lapse.shape[0]
    kinetic_plus = np.zeros(points, dtype=complex)
    kinetic_minus = np.zeros(points, dtype=complex)
    angular_scalar = np.zeros(points, dtype=complex)
    probability = np.zeros(points)
    speed_plus, speed_minus = _speeds_from_fields(lapse, shift, radial_density)
    for column in range(rank):
        for row in range(rank):
            amplitude = coefficients[row, column]
            if amplitude == 0.0:
                continue
            _density, flux_plus, flux_minus, angular = _pair_integrand(
                phi[:, :, row], phi[:, :, column], lapse, shift, radial_density,
                radius, derivative, kappa,
            )
            kinetic_plus += amplitude * flux_plus
            kinetic_minus += amplitude * flux_minus
            angular_scalar += amplitude * angular
            probability += np.real(amplitude * (
                speed_plus * np.conj(phi[0, :, row]) * phi[0, :, column]
                + speed_minus * np.conj(phi[1, :, row]) * phi[1, :, column]
            ))
    flux_plus = np.real(kinetic_plus)
    flux_minus = np.real(kinetic_minus)
    angular_scalar = np.real(angular_scalar)
    force_n = multiplicity * (
        (flux_plus - flux_minus) / radial_density + (kappa / radius) * angular_scalar
    )
    force_beta = multiplicity * (-(flux_plus + flux_minus) / radial_density)
    force_q = multiplicity * (
        -(speed_plus * flux_plus + speed_minus * flux_minus) / radial_density
    )
    force_r = multiplicity * (-(kappa * lapse / radius**2) * angular_scalar)
    angular_force = multiplicity * (kappa / radius) * angular_scalar
    coordinate_energy = lapse * force_n + shift * force_beta
    volume = radial_density * radius**2
    angular_rho = force_n / volume
    angular_momentum = force_beta / volume
    physical_rho = force_n / (FOUR_PI * volume)
    physical_j = -force_beta / (FOUR_PI * volume)
    physical_pr = (force_n - angular_force) / (FOUR_PI * volume)
    physical_pa = -force_r / (2.0 * lapse * radial_density * radius * FOUR_PI)
    conformal_q = radial_density / radius
    return {
        "schema": "NSC-RN-NEUTRAL-SOURCE-v2",
        "kappa": int(kappa),
        "multiplicity": multiplicity,
        "multiplicity_applied_once": True,
        "gauge_fix_imposed": False,
        "q_adm": radial_density,
        "conformal_Q": conformal_q,
        "shift_form": shift_form,
        "neutral": True,
        "electric_current": np.zeros(points),
        "magnetic_flux_fixed": True,
        "weyl_variation_included": False,
        "rho": physical_rho,
        "current": physical_j,
        "radial_pressure": physical_pr,
        "angular_pressure": physical_pa,
        "coordinate_energy": coordinate_energy,
        "probability_current": multiplicity * probability,
        "forces": {"N": force_n, "beta": force_beta, "q": force_q, "r": force_r},
        "angular_integrated": {
            "force_N": force_n,
            "force_beta": force_beta,
            "force_q": force_q,
            "force_r": force_r,
            "rho": angular_rho,
            "radial_momentum": angular_momentum,
            "includes_four_pi": False,
        },
        "physical": {
            "rho": physical_rho,
            "j": physical_j,
            "radial_pressure": physical_pr,
            "angular_pressure": physical_pa,
            "four_pi": FOUR_PI,
            "momentum_sign": "j = -F_beta / (4π q r²)",
            "pg_combination": "4π G r² (rho - v j)",
            "do_not_multiply_four_pi_again": True,
            "do_not_flip_j": True,
        },
        "killing_energy": float(np.sum(weights * coordinate_energy)),
        "sector_killing_energy": float(np.sum(weights * coordinate_energy) / multiplicity),
        "normal_energy": float(np.sum(weights * force_n)),
        "massless_trace_residual": physical_rho - physical_pr - 2.0 * physical_pa,
        "dense_covariance": False,
        "coherent_rank_one_column": coherent_rank_one,
    }


def _split_derivative(coeff, component, derivative):
    return 0.5 * (coeff * (derivative @ component) + derivative @ (coeff * component))


def _weighted_pairing(left, right, weights):
    return np.sum(weights * np.sum(np.conj(left) * right, axis=0))


def _require_diagnostic_column(phi, radial_density):
    """Energy diagnostics own one column on the areal chart.

    The source holder may still store several spectral columns. These
    ledgers do not choose one of them silently, and they do not accept a
    radial density other than q_adm = 1.
    """
    array = np.asarray(phi)
    if array.ndim != 3 or array.shape[0] != 2 or array.shape[2] != 1:
        raise ValueError(
            "energy diagnostics require shape (2, points, 1); "
            "the source holder may still store several columns"
        )
    if radial_density is not None and float(np.max(np.abs(np.asarray(radial_density, dtype=float) - 1.0))) > 1e-12:
        raise ValueError("energy diagnostics require q_adm = 1")


def _polarized_force_n(phi, eta, lapse, radius, derivative, occupations, kappa):
    """Directional derivative of F_N along eta. Shell 4 and ν enter once."""
    _require_diagnostic_column(phi, None)
    column = np.asarray(phi[:, :, 0], dtype=complex)
    direction = np.asarray(eta[:, :, 0], dtype=complex)
    plus_flux = 2.0 * np.real(_hermitian_flux(column[0], direction[0], derivative))
    minus_flux = 2.0 * np.real(_hermitian_flux(column[1], direction[1], derivative))
    angular = np.real(
        _angular_bilinear(column, direction) + _angular_bilinear(direction, column)
    )
    density = plus_flux - minus_flux + (float(kappa) / radius) * angular
    return float(SHELL_MULTIPLICITY * occupations[0]) * density


def energy_account(phi, occupations, lapse, shift, radius, weights, derivative, *,
                   phi_rate=None, lapse_t=None, shift_t=None, lapse_r=None, shift_r=None,
                   kappa=SHELL_KAPPA, radial_density=None):
    """Coordinate and normal energies, boundary flux, SAT energy, and metric power.

    The evolved column is one rank. H is defined by the bulk characteristic
    rate, bulk = -i H. The SBP identity is H†W - WH = i B A with
    A = diag(N-β, -N-β) and B = diag(-1, 0, ..., 1). Shell 4 and the
    occupation multiply that pairing once. The raw probability SAT stock is
    not an energy.
    """
    if radial_density is None:
        radial_density = np.ones_like(np.asarray(lapse, dtype=float))
    _require_diagnostic_column(phi, radial_density)
    evaluated = source_evaluate(
        phi, occupations, lapse, shift, radial_density, radius, weights, derivative, kappa=kappa,
    )
    phi = np.asarray(phi, dtype=complex)
    occupations = np.asarray(occupations, dtype=float)
    weights = np.asarray(weights, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    radius = np.asarray(radius, dtype=float)
    derivative = np.asarray(derivative, dtype=float)
    force_n = evaluated["forces"]["N"]
    force_beta = evaluated["forces"]["beta"]
    force_r_unit = evaluated["forces"]["r"] / lapse
    speed_plus, speed_minus = _speeds_from_fields(lapse, shift, radial_density)
    plus = phi[0, :, 0]
    minus = phi[1, :, 0]
    alpha = float(kappa) * lapse / radius
    bulk = np.zeros_like(phi)
    bulk[0, :, 0] = (
        -_split_derivative(lapse, plus, derivative) + _split_derivative(shift, plus, derivative)
        - alpha * minus
    )
    bulk[1, :, 0] = (
        _split_derivative(lapse, minus, derivative) + _split_derivative(shift, minus, derivative)
        + alpha * plus
    )
    incoming = float(speed_minus[-1])
    sat = np.zeros_like(phi)
    sat[1, -1, 0] = (incoming / weights[-1]) * minus[-1]
    probability_sat_raw = float(-incoming * np.abs(minus[-1]) ** 2)
    column = phi[:, :, 0]
    hamiltonian = 1j * bulk[:, :, 0]
    channel_energy = float(np.real(_weighted_pairing(column, hamiltonian, weights)))
    shell = float(evaluated["multiplicity"] * occupations[0])
    full_rate = bulk if phi_rate is None else np.asarray(phi_rate, dtype=complex)
    if full_rate.shape != phi.shape:
        raise ValueError("phi_rate must have the column shape")
    full = full_rate[:, :, 0]
    bulk_column = bulk[:, :, 0]
    sat_column = sat[:, :, 0]

    def boundary_scalar(other):
        value = -speed_plus[0] * np.conj(column[0, 0]) * other[0, 0]
        value += -speed_minus[0] * np.conj(column[1, 0]) * other[1, 0]
        value += speed_plus[-1] * np.conj(column[0, -1]) * other[0, -1]
        value += speed_minus[-1] * np.conj(column[1, -1]) * other[1, -1]
        return value

    bulk_pairing = shell * 2.0 * np.real(_weighted_pairing(hamiltonian, bulk_column, weights))
    sat_energy_rate = shell * 2.0 * np.real(_weighted_pairing(hamiltonian, sat_column, weights))
    phi_motion_rate = shell * (
        2.0 * np.real(_weighted_pairing(hamiltonian, full, weights))
        + np.imag(boundary_scalar(full))
    )
    boundary_energy_flux = shell * float(np.imag(boundary_scalar(bulk_column)))
    boundary_re_A_H = shell * float(np.real(boundary_scalar(hamiltonian)))
    if lapse_t is None:
        lapse_t = np.zeros_like(lapse)
    if shift_t is None:
        shift_t = np.zeros_like(shift)
    metric_power = float(np.sum(weights * (force_n * lapse_t + force_beta * shift_t)))
    continuity_residual = None
    if lapse_r is not None and shift_r is not None and phi_rate is not None:
        # F_N does not depend on N or β at fixed r and q_adm = 1.
        force_n_t = _polarized_force_n(
            phi, full_rate, lapse, radius, derivative, occupations, kappa,
        )
        flux = -lapse * force_beta - shift * force_n
        divergence = derivative @ flux
        right_hand = (
            np.asarray(shift_r, dtype=float) * (force_n + radius * force_r_unit)
            - shift * force_r_unit
            + np.asarray(lapse_r, dtype=float) * force_beta
        )
        continuity_residual = float(np.max(np.abs(force_n_t + divergence - right_hand)))
    return {
        "E_coordinate": float(np.sum(weights * evaluated["coordinate_energy"])),
        "E_normal": float(evaluated["normal_energy"]),
        "channel_energy": channel_energy,
        "shell_occupation_factor": shell,
        "multiplicity_applied_once": True,
        "e_H_density": evaluated["coordinate_energy"],
        "F_N": force_n,
        "F_beta": force_beta,
        "F_r_unit": force_r_unit,
        "j_int": -force_beta,
        "pr_int": force_n + radius * force_r_unit,
        "pt_int": -0.5 * radius * force_r_unit,
        "metric_power": metric_power,
        "boundary_energy_flux": boundary_energy_flux,
        "boundary_re_A_H": boundary_re_A_H,
        "bulk_pairing_rate": float(np.real(bulk_pairing)),
        "sat_energy_rate": float(np.real(sat_energy_rate)),
        "phi_motion_rate": float(np.real(phi_motion_rate)),
        "predicted_coordinate_rate": float(np.real(phi_motion_rate) + metric_power),
        "probability_sat_raw": probability_sat_raw,
        "probability_sat_is_energy": False,
        "normal_continuity_residual": continuity_residual,
        "field_energy_not_added_to_mass": True,
    }


def normal_energy_ledger(phi, occupations, lapse, shift, radius, weights, derivative, *,
                         phi_rate, lapse_t, shift_t, lapse_r, shift_r, kappa=SHELL_KAPPA,
                         radial_density=None):
    """Attribute the normal-energy balance without forcing it to zero.

    J_N = -N F_β - β F_N. Pressure work is β_r (F_N + r F_r/N) - β F_r/N
    and lapse work is N_r F_β. The SAT piece is the polarized F_N along the
    outer incoming SAT rate, not the raw probability coefficient.
    """
    if radial_density is None:
        radial_density = np.ones_like(np.asarray(lapse, dtype=float))
    _require_diagnostic_column(phi, radial_density)
    evaluated = source_evaluate(
        phi, occupations, lapse, shift, radial_density, radius, weights, derivative, kappa=kappa,
    )
    phi = np.asarray(phi, dtype=complex)
    occupations = np.asarray(occupations, dtype=float)
    weights = np.asarray(weights, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    radius = np.asarray(radius, dtype=float)
    derivative = np.asarray(derivative, dtype=float)
    phi_rate = np.asarray(phi_rate, dtype=complex)
    if phi_rate.shape != phi.shape:
        raise ValueError("phi_rate must have shape (2, points, 1)")
    lapse_r = np.asarray(lapse_r, dtype=float)
    shift_r = np.asarray(shift_r, dtype=float)
    lapse_t = np.asarray(lapse_t, dtype=float)
    shift_t = np.asarray(shift_t, dtype=float)
    if not np.all(np.isfinite(lapse_t)) or not np.all(np.isfinite(shift_t)):
        raise ValueError("lapse and shift rates must be finite; they do not enter F_N")
    force_n = evaluated["forces"]["N"]
    force_beta = evaluated["forces"]["beta"]
    force_r_unit = evaluated["forces"]["r"] / lapse
    # F_N is independent of N and β on this fixed-radius, q_adm = 1 slice.
    force_n_t = _polarized_force_n(
        phi, phi_rate, lapse, radius, derivative, occupations, kappa,
    )
    speed_minus = _speeds_from_fields(lapse, shift, radial_density)[1]
    sat = np.zeros_like(phi)
    sat[1, -1, 0] = (float(speed_minus[-1]) / weights[-1]) * phi[1, -1, 0]
    sat_density = _polarized_force_n(phi, sat, lapse, radius, derivative, occupations, kappa)
    flux = -lapse * force_beta - shift * force_n
    divergence = derivative @ flux
    pressure = shift_r * (force_n + radius * force_r_unit) - shift * force_r_unit
    lapse_work = lapse_r * force_beta
    full_residual = force_n_t + divergence - pressure - lapse_work
    bulk_residual = full_residual - sat_density
    commutator_beta = derivative @ (shift * force_n) - (
        shift_r * force_n + shift * (derivative @ force_n)
    )
    commutator_lapse = derivative @ (lapse * force_beta) - (
        lapse_r * force_beta + lapse * (derivative @ force_beta)
    )
    dE_N = float(np.sum(weights * force_n_t))
    sat_energy = float(np.sum(weights * sat_density))
    pressure_integral = float(np.sum(weights * pressure))
    lapse_integral = float(np.sum(weights * lapse_work))
    global_balance = float(flux[0] - flux[-1] + pressure_integral + lapse_integral + sat_energy)
    interior = slice(1, -1)
    full_index = int(np.argmax(np.abs(full_residual)))
    bulk_index = int(np.argmax(np.abs(bulk_residual)))
    return {
        "E_normal": float(np.sum(weights * force_n)),
        "dE_normal": dE_N,
        "J_N_left": float(flux[0]),
        "J_N_right": float(flux[-1]),
        "pressure_integral": pressure_integral,
        "lapse_integral": lapse_integral,
        "sat_normal_energy": sat_energy,
        "global_residual": dE_N - global_balance,
        "full_residual_max": float(np.max(np.abs(full_residual))),
        "full_residual_index": full_index,
        "full_residual_radius": float(radius[full_index]),
        "full_residual_value": float(full_residual[full_index]),
        "full_residual_at_outer": float(full_residual[-1]),
        "outer_is_endpoint": full_index == radius.size - 1,
        "sat_density_at_outer": float(sat_density[-1]),
        "bulk_residual_max_outside_stencil": float(np.max(np.abs(bulk_residual[interior]))),
        "sbp_nodes_excluded_total": 2,
        "bulk_residual_max": float(np.max(np.abs(bulk_residual))),
        "bulk_residual_index": bulk_index,
        "bulk_residual_radius": float(radius[bulk_index]),
        "product_rule_beta_max": float(np.max(np.abs(commutator_beta))),
        "product_rule_beta_index": int(np.argmax(np.abs(commutator_beta))),
        "product_rule_lapse_max": float(np.max(np.abs(commutator_lapse))),
        "shell_occupation_factor": float(SHELL_MULTIPLICITY * occupations[0]),
        "probability_sat_used_as_energy": False,
        "F_N_time_derivative": "bilinear_polarization",
        "lapse_and_shift_rates_enter_F_N": False,
        "endpoint_stencil": "SBP endpoints, index 0 and the last node",
    }
