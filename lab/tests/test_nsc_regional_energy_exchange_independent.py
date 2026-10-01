"""Identity controls for observer stress and the two energy accounts.

The oracle is the variation derived in this file. It does not import the
regional exchange module and it does not evolve a self-gravitating solution.

Signature is +---. The line element is

    ds² = N² dt² − q² (dx + β dt)² − r² dΩ²,

with physical fields N = r L and q = r Q. The future unit normal is the
vector n = N⁻¹ (∂_t − β ∂_x). Its metric dual is N dt. The 1-form
(dt − β dx)/N is not that normal. The spatial leg is q (dx + β dt).

V = 4π q r² is the proper spatial volume per unit coordinate x. FL, FQ and
Fβ in the formulas below are continuum densities: nodal partials of
M Tr(C H) divided by the node spacing. The conformal owner stores those
partials before the division.

Raised orthonormal flux j = T^{0̂1̂} and the pressures are

    ρ = FL / (4π r⁴ Q) = F_N / V,
    p_r = −FQ / (4π r⁴ L) = −F_q / (4π N r²),
    p_⊥ = (L FL + Q FQ) / (8π r⁴ L Q) = −F_r / (8π N q r),
    j = −Fβ / (4π r⁴ Q²) = −Fβ / (4π q² r²).

The ADM constraint note's symbol j is the lowered component T_{0̂1̂} = −j.
Positive expansion along n is

    K_r = (q_t − β q_x − q β_x) / (N q),
    K_⊥ = (r_t − β r_x) / (N r),

and ∇_μ n^μ = K_r + 2 K_⊥. The normal-energy balance

    ∂_t(V ρ) + ∂_x[V (N j / q − β ρ)]
        = −N V (p_r K_r + 2 p_⊥ K_⊥) − V (j / q) N_x

is exactly (W_t − β W_x) / N, where W_t and W_x are the two diffeomorphism
Ward residuals of (F_N, F_q, F_r, F_β). It holds when those residuals vanish.
It is not an extra dynamical law.

The coordinate Hamiltonian density e = L FL + β Fβ = N V ρ − β q V j
integrates to M Tr(C H). Its chart power is FL L_t + FQ Q_t + Fβ β_t.
That power is not the proper-pressure term −N V (p_r K_r + 2 p_⊥ K_⊥).
Vanishing expansion still allows Q_t = (β Q)_x, and then FQ Q_t can be
nonzero while the pressure term is zero. The normal-energy integral changes
by the lapse-gradient term instead. A fixed window keeps the shift advection
−β ρ in its flux. A window dragged by ∂_t W = β ∂_x W, so its points move at
dx/dt = −β, keeps the orthonormal flux 4π r² N j. A state rate that is not
the Hamiltonian commutator adds Tr(Ċ H) to the canonical ledger and a local
chart-flux residual. Neither addition is proper pressure work.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_conformal_adm_source import (
    conformal_source,
    fixed_gaussian_covariance,
)
from recursive_horizons.nsc_covariant_operator import smooth_metric


FOUR_PI = 4.0 * np.pi


def _periodic_operators(points=65, length=4.0):
    """Integer Fourier modes on an odd periodic grid. Differentiation is exact there."""
    if points % 2 != 1:
        raise ValueError("odd point count keeps a single zero mode")
    coordinate = (np.arange(points) - points // 2) * (length / points)
    wave = 2.0 * np.pi * (np.arange(points) - points // 2) / length
    phase = np.exp(1j * coordinate[:, None] * wave[None, :]) / np.sqrt(points)

    def coefficients(values):
        return phase.conj().T @ np.asarray(values, dtype=complex)

    def synthesize(coeff):
        return phase @ np.asarray(coeff, dtype=complex)

    def differentiate(values):
        return synthesize(1j * wave * coefficients(values))

    def momentum(values):
        return synthesize(wave * coefficients(values))

    def project(values, max_mode):
        coeff = coefficients(values)
        mode = np.arange(points) - points // 2
        coeff = np.where(np.abs(mode) <= max_mode, coeff, 0.0)
        return synthesize(coeff)

    return {
        "x": coordinate,
        "length": float(length),
        "dx": float(length / points),
        "wave": wave,
        "differentiate": differentiate,
        "momentum": momentum,
        "project": project,
        "synthesize": synthesize,
        "coefficients": coefficients,
    }


def _real_derivative(grid, values):
    return np.real(grid["differentiate"](values))


def _observer_from_physical(lapse, radial, radius, force_n, force_q, force_r, force_beta):
    """Map continuum force densities to (ρ, j, p_r, p_⊥). j is raised."""
    volume = FOUR_PI * radial * radius**2
    density = force_n / volume
    flux = -force_beta / (FOUR_PI * radial**2 * radius**2)
    radial_pressure = -force_q / (FOUR_PI * lapse * radius**2)
    angular_pressure = -force_r / (8.0 * np.pi * lapse * radial * radius)
    return density, flux, radial_pressure, angular_pressure, volume


def _observer_from_conformal(length_density, radial_density, radius, force_l, force_q, force_beta):
    """Parent dictionary. Inputs are continuum densities, already divided by dx."""
    denominator = FOUR_PI * radius**4
    density = force_l / (denominator * radial_density)
    radial_pressure = -force_q / (denominator * length_density)
    angular_pressure = (
        (length_density * force_l + radial_density * force_q)
        / (2.0 * denominator * length_density * radial_density)
    )
    flux = -force_beta / (denominator * radial_density**2)
    return density, flux, radial_pressure, angular_pressure


def _expansion(lapse, radial, radius, shift, radial_t, radius_t, radial_x, radius_x, shift_x):
    radial_curvature = (radial_t - shift * radial_x - radial * shift_x) / (lapse * radial)
    angular_curvature = (radius_t - shift * radius_x) / (lapse * radius)
    return radial_curvature, angular_curvature


def _normal_balance_terms(lapse, radial, radius, shift, density, flux, radial_pressure,
                          angular_pressure, radial_curvature, angular_curvature, lapse_x):
    volume = FOUR_PI * radial * radius**2
    coordinate_flux = volume * (lapse * flux / radial - shift * density)
    proper_work = -lapse * volume * (
        radial_pressure * radial_curvature + 2.0 * angular_pressure * angular_curvature
    )
    lapse_gradient = -volume * (flux / radial) * lapse_x
    return coordinate_flux, proper_work, lapse_gradient


def _ward_residuals(grid, lapse, radial, radius, shift, force_n, force_q, force_r, force_beta,
                    lapse_t, radial_t, radius_t, shift_t, force_n_t, force_q_t, force_r_t,
                    force_beta_t):
    """Continuum Ward residuals W_t and W_x. Zero when the forces are a scalar density."""
    derivative = lambda values: _real_derivative(grid, values)
    energy = lapse * force_n + shift * force_beta
    current = (
        -lapse * shift * force_n
        - (shift**2 + lapse**2 / radial**2) * force_beta
        + radial * shift * force_q
    )
    radial_current = radial * force_q - shift * force_beta
    energy_t = lapse_t * force_n + lapse * force_n_t + shift_t * force_beta + shift * force_beta_t
    source_t = (
        force_n * lapse_t + force_q * radial_t + force_r * radius_t + force_beta * shift_t
    )
    source_x = (
        force_n * derivative(lapse)
        + force_q * derivative(radial)
        + force_r * derivative(radius)
        + force_beta * derivative(shift)
    )
    time_residual = energy_t + derivative(current) - source_t
    space_residual = force_beta_t + derivative(radial_current) - source_x
    return time_residual, space_residual, energy, current


def _trig_background(grid):
    """Low trigonometric metric and force profiles. Smooth enough for this grid."""
    angle = 2.0 * np.pi * grid["x"] / grid["length"]
    length_density = 1.2 + 0.15 * np.cos(angle)
    radial_density = 0.85 + 0.12 * np.sin(angle)
    radius = 2.0 + 0.25 * np.sin(angle)
    shift = 0.35 * np.cos(angle)
    force_l = 0.30 + 0.05 * np.sin(angle) + 0.02 * np.cos(2.0 * angle)
    force_q = -0.12 + 0.03 * np.cos(angle)
    force_beta = 0.04 * np.sin(angle) - 0.02 * np.cos(2.0 * angle)
    length_t = 0.02 * np.sin(angle)
    radial_t = 0.03 * np.cos(2.0 * angle)
    radius_t = -0.015 * np.sin(angle)
    shift_t = 0.01 * np.cos(angle)
    force_l_t = 0.004 * np.cos(angle)
    force_q_t = -0.003 * np.sin(2.0 * angle)
    force_beta_t = 0.002 * np.cos(angle)
    return {
        "L": length_density,
        "Q": radial_density,
        "r": radius,
        "beta": shift,
        "FL": force_l,
        "FQ": force_q,
        "FB": force_beta,
        "Lt": length_t,
        "Qt": radial_t,
        "rt": radius_t,
        "bt": shift_t,
        "FLt": force_l_t,
        "FQt": force_q_t,
        "FBt": force_beta_t,
    }


def _physical_from_conformal_state(state):
    """Chain rule at fixed state. F_r is not an independent conformal partial."""
    radius = state["r"]
    lapse = radius * state["L"]
    radial = radius * state["Q"]
    force_n = state["FL"] / radius
    force_q = state["FQ"] / radius
    force_beta = state["FB"]
    force_r = -(state["L"] * state["FL"] + state["Q"] * state["FQ"]) / radius
    force_n_t = state["FLt"] / radius - state["FL"] * state["rt"] / radius**2
    force_q_t = state["FQt"] / radius - state["FQ"] * state["rt"] / radius**2
    force_beta_t = state["FBt"]
    amplitude_t = (
        state["Lt"] * state["FL"] + state["L"] * state["FLt"]
        + state["Qt"] * state["FQ"] + state["Q"] * state["FQt"]
    )
    radius_force_t = -amplitude_t / radius - force_r * state["rt"] / radius
    lapse_t = state["rt"] * state["L"] + radius * state["Lt"]
    radial_t = state["rt"] * state["Q"] + radius * state["Qt"]
    return {
        "N": lapse,
        "q": radial,
        "r": radius,
        "beta": state["beta"],
        "FN": force_n,
        "Fq": force_q,
        "Fr": force_r,
        "FB": force_beta,
        "Nt": lapse_t,
        "qt": radial_t,
        "rt": state["rt"],
        "bt": state["bt"],
        "FNt": force_n_t,
        "Fqt": force_q_t,
        "Frt": radius_force_t,
        "FBt": force_beta_t,
    }


def _dirac_images(grid, first, second, length_density, radial_density, shift, kappa=1.0):
    """One block of H = σ2 {L/Q, P}/2 + σ1 κ L − {β, P}/2."""
    momentum = grid["momentum"]

    def anticommutator(values, component):
        return 0.5 * (values * momentum(component) + momentum(values * component))

    slope = length_density / radial_density
    mass = kappa * length_density
    image0 = (
        -1j * anticommutator(slope, second)
        + mass * second
        - anticommutator(shift, first)
    )
    image1 = (
        1j * anticommutator(slope, first)
        + mass * first
        - anticommutator(shift, second)
    )
    return image0, image1


def _dirac_kernels(grid, first, second):
    momentum = grid["momentum"]
    momentum0 = momentum(first)
    momentum1 = momentum(second)
    kinetic = np.conjugate(first) * (-1j * momentum1) + np.conjugate(second) * (1j * momentum0)
    sigma = np.conjugate(first) * second + np.conjugate(second) * first
    current = np.conjugate(first) * momentum0 + np.conjugate(second) * momentum1
    return kinetic.real, sigma.real, current.real


def _deterministic_spinor(grid):
    points = len(grid["x"])
    first = np.zeros(points, dtype=complex)
    second = np.zeros(points, dtype=complex)
    center = points // 2
    first_coeff = np.zeros(points, dtype=complex)
    second_coeff = np.zeros(points, dtype=complex)
    first_coeff[center + 1] = 0.4
    first_coeff[center - 1] = 0.2
    second_coeff[center + 2] = 0.3
    second_coeff[center] = 0.25
    first = grid["synthesize"](first_coeff)
    second = grid["synthesize"](second_coeff)
    norm = np.sqrt(np.sum((np.abs(first) ** 2 + np.abs(second) ** 2) * grid["dx"]))
    return first / norm, second / norm


def _chart_samples(grid, first, second, length_density, radial_density, shift, kappa=1.0):
    """Densities whose integral is M Tr(C H), M = 4κ, plus the chart flux."""
    multiplicity = 4.0 * kappa
    kinetic, sigma, current = _dirac_kernels(grid, first, second)
    slope = length_density / radial_density
    mass = kappa * length_density
    energy = multiplicity * (slope * kinetic + mass * sigma - shift * current)
    force_l = multiplicity * (kinetic / radial_density + kappa * sigma)
    force_q = -multiplicity * length_density * kinetic / radial_density**2
    force_beta = -multiplicity * current
    flux = (
        -shift * energy
        + multiplicity * slope**2 * current
        - multiplicity * slope * shift * kinetic
    )
    return {
        "energy": energy,
        "flux": flux,
        "FL": force_l,
        "FQ": force_q,
        "FB": force_beta,
        "K": kinetic,
        "S": sigma,
        "J": current,
    }


def _advance_spinor(grid, first, second, length_density, radial_density, shift,
                    length_t, radial_t, shift_t, step, kappa=1.0):
    """One Heun step of i ∂_t ψ = H ψ with an explicit chart velocity."""
    def images_at(sample0, sample1, length, radial, shift_value):
        return _dirac_images(grid, sample0, sample1, length, radial, shift_value, kappa)

    image0, image1 = images_at(first, second, length_density, radial_density, shift)
    midpoint0 = first - 1j * (step / 2.0) * image0
    midpoint1 = second - 1j * (step / 2.0) * image1
    image0, image1 = images_at(
        midpoint0,
        midpoint1,
        length_density + (step / 2.0) * length_t,
        radial_density + (step / 2.0) * radial_t,
        shift + (step / 2.0) * shift_t,
    )
    return first - 1j * step * image0, second - 1j * step * image1


def _integral(grid, values):
    return float(np.sum(np.real(values)) * grid["dx"])


def test_normal_vector_and_stress_dictionary_are_the_parent_map():
    """Raised flux, chain rule, and the unit normal dual. Lowered j is the other sign."""
    time, radius_coord = sp.symbols("t x", real=True)
    lapse, radial, radius, shift = (
        sp.Function(name)(time, radius_coord) for name in ("N", "q", "r", "beta")
    )
    density, lowered, radial_pressure, angular_pressure = (
        sp.Function(name)(time, radius_coord) for name in ("rho", "j_down", "pr", "pperp")
    )
    metric = sp.Matrix([
        [lapse**2 - radial**2 * shift**2, -radial**2 * shift],
        [-radial**2 * shift, -radial**2],
    ])
    inverse = sp.simplify(metric.inv())
    vector = sp.Matrix([1 / lapse, -shift / lapse])
    lowered_normal = sp.simplify(metric * vector)
    assert lowered_normal == sp.Matrix([lapse, 0])
    claimed_form = sp.Matrix([1 / lapse, -shift / lapse])
    claimed_vector = sp.simplify(inverse * claimed_form)
    assert sp.simplify(claimed_vector - vector) != sp.Matrix([0, 0])

    coframe = sp.Matrix([[lapse, 0], [radial * shift, radial]])
    lowered_stress = sp.Matrix([[density, lowered], [lowered, radial_pressure]])
    coordinate_stress = sp.simplify(coframe.T * lowered_stress * coframe)
    volume = 4 * sp.pi * lapse * radial * radius**2

    def metric_force(field):
        contraction = sum(
            coordinate_stress[row, column] * sp.diff(inverse[row, column], field)
            for row in range(2) for column in range(2)
        )
        return sp.factor(sp.simplify(-volume / 2 * contraction))

    assert sp.simplify(metric_force(lapse) - 4 * sp.pi * radial * radius**2 * density) == 0
    assert sp.simplify(metric_force(shift) - 4 * sp.pi * radial**2 * radius**2 * lowered) == 0
    assert sp.simplify(metric_force(radial) + 4 * sp.pi * lapse * radius**2 * radial_pressure) == 0
    angular_force = sp.simplify(-volume / 2 * (4 * angular_pressure / radius))
    assert sp.simplify(angular_force + 8 * sp.pi * lapse * radial * radius * angular_pressure) == 0
    eta = sp.diag(1, -1)
    raised = sp.simplify(eta * lowered_stress * eta)
    assert raised[0, 1] == -lowered

    grid = _periodic_operators()
    state = _physical_from_conformal_state(_trig_background(grid))
    parent = _observer_from_conformal(
        state["N"] / state["r"], state["q"] / state["r"], state["r"],
        state["FN"] * state["r"], state["Fq"] * state["r"], state["FB"],
    )
    physical = _observer_from_physical(
        state["N"], state["q"], state["r"], state["FN"], state["Fq"], state["Fr"], state["FB"],
    )
    for parent_value, physical_value in zip(parent, physical[:4]):
        assert np.max(np.abs(parent_value - physical_value)) < 1e-12


def test_normal_balance_residual_equals_the_weighted_ward_residuals():
    """The parent PDE is (W_t − β W_x)/N, including the static zero-shift reduction."""
    time, space = sp.symbols("t x", real=True)
    lapse, radial, radius, shift = (
        sp.Function(name)(time, space) for name in ("N", "q", "r", "beta")
    )
    force_n, force_q, force_r, force_beta = (
        sp.Function(name)(time, space) for name in ("FN", "Fq", "Fr", "Fb")
    )
    fields = (lapse, radial, radius, shift)
    forces = (force_n, force_q, force_r, force_beta)
    energy = lapse * force_n + shift * force_beta
    current = (
        -lapse * shift * force_n
        - (shift**2 + lapse**2 / radial**2) * force_beta
        + radial * shift * force_q
    )
    time_ward = (
        sp.diff(energy, time) + sp.diff(current, space)
        - sum(force * sp.diff(field, time) for force, field in zip(forces, fields))
    )
    space_ward = (
        sp.diff(force_beta, time) + sp.diff(radial * force_q - shift * force_beta, space)
        - sum(force * sp.diff(field, space) for force, field in zip(forces, fields))
    )
    volume = 4 * sp.pi * radial * radius**2
    density = force_n / volume
    flux = -force_beta / (4 * sp.pi * radial**2 * radius**2)
    radial_pressure = -force_q / (4 * sp.pi * lapse * radius**2)
    angular_pressure = -force_r / (8 * sp.pi * lapse * radial * radius)
    radial_curvature = (
        sp.diff(radial, time) - shift * sp.diff(radial, space) - radial * sp.diff(shift, space)
    ) / (lapse * radial)
    angular_curvature = (
        sp.diff(radius, time) - shift * sp.diff(radius, space)
    ) / (lapse * radius)
    coordinate_flux = volume * (lapse * flux / radial - shift * density)
    balance = (
        sp.diff(volume * density, time) + sp.diff(coordinate_flux, space)
        + lapse * volume * (
            radial_pressure * radial_curvature + 2 * angular_pressure * angular_curvature
        )
        + volume * (flux / radial) * sp.diff(lapse, space)
    )
    assert sp.simplify(sp.expand(balance - time_ward / lapse + shift * space_ward / lapse)) == 0

    divergence = (
        sp.diff(radial * radius**2, time) - sp.diff(shift * radial * radius**2, space)
    ) / (lapse * radial * radius**2)
    assert sp.simplify(sp.expand(divergence - radial_curvature - 2 * angular_curvature)) == 0

    static_space = sp.symbols("y", real=True)
    static_n, static_q, static_r = (sp.Function(name)(static_space) for name in ("N", "q", "r"))
    static_fn, static_fq, static_fr = (
        sp.Function(name)(static_space) for name in ("FN", "Fq", "Fr")
    )
    static_density = static_fn / (4 * sp.pi * static_q * static_r**2)
    static_pr = -static_fq / (4 * sp.pi * static_n * static_r**2)
    static_pp = -static_fr / (8 * sp.pi * static_n * static_q * static_r)
    static_ward = (
        sp.diff(static_q * static_fq, static_space)
        - static_fn * sp.diff(static_n, static_space)
        - static_fq * sp.diff(static_q, static_space)
        - static_fr * sp.diff(static_r, static_space)
    )
    hydrostatic = (
        sp.diff(static_pr, static_space)
        + sp.diff(static_n, static_space) / static_n * (static_density + static_pr)
        + 2 * sp.diff(static_r, static_space) / static_r * (static_pr - static_pp)
    )
    assert sp.simplify(static_ward - (-4 * sp.pi * static_n * static_q * static_r**2) * hydrostatic) == 0

    grid = _periodic_operators()
    state = _physical_from_conformal_state(_trig_background(grid))
    derivative = lambda values: _real_derivative(grid, values)
    time_residual, space_residual, _energy, _current = _ward_residuals(
        grid, state["N"], state["q"], state["r"], state["beta"],
        state["FN"], state["Fq"], state["Fr"], state["FB"],
        state["Nt"], state["qt"], state["rt"], state["bt"],
        state["FNt"], state["Fqt"], state["Frt"], state["FBt"],
    )
    density, flux, radial_pressure, angular_pressure, _volume = _observer_from_physical(
        state["N"], state["q"], state["r"], state["FN"], state["Fq"], state["Fr"], state["FB"],
    )
    radial_curvature, angular_curvature = _expansion(
        state["N"], state["q"], state["r"], state["beta"],
        state["qt"], state["rt"], derivative(state["q"]), derivative(state["r"]),
        derivative(state["beta"]),
    )
    coordinate_flux, proper_work, lapse_gradient = _normal_balance_terms(
        state["N"], state["q"], state["r"], state["beta"],
        density, flux, radial_pressure, angular_pressure,
        radial_curvature, angular_curvature, derivative(state["N"]),
    )
    balance_residual = (
        state["FNt"] + derivative(coordinate_flux) - proper_work - lapse_gradient
    )
    weighted = (time_residual - state["beta"] * space_residual) / state["N"]
    assert np.max(np.abs(balance_residual - weighted)) < 1e-9

    lowered_flux = -flux
    lowered_coordinate, _work, lowered_gradient = _normal_balance_terms(
        state["N"], state["q"], state["r"], state["beta"],
        density, lowered_flux, radial_pressure, angular_pressure,
        radial_curvature, angular_curvature, derivative(state["N"]),
    )
    lowered_residual = (
        state["FNt"] + derivative(lowered_coordinate) - proper_work - lowered_gradient
    )
    assert np.max(np.abs(lowered_residual)) > 1e-3

    theta = (
        state["qt"] * state["r"]**2 + 2.0 * state["q"] * state["r"] * state["rt"]
        - derivative(state["beta"] * state["q"] * state["r"]**2)
    ) / (state["N"] * state["q"] * state["r"]**2)
    assert np.max(np.abs(theta - radial_curvature - 2.0 * angular_curvature)) < 1e-11


def test_fixed_window_keeps_shift_advection_and_moving_window_does_not():
    """A fixed weight sees V(Nj/q − βρ). Dragging the weight at dx/dt = −β removes −βρ."""
    grid = _periodic_operators()
    state = _physical_from_conformal_state(_trig_background(grid))
    derivative = lambda values: _real_derivative(grid, values)
    density, flux, radial_pressure, angular_pressure, volume = _observer_from_physical(
        state["N"], state["q"], state["r"], state["FN"], state["Fq"], state["Fr"], state["FB"],
    )
    radial_curvature, angular_curvature = _expansion(
        state["N"], state["q"], state["r"], state["beta"],
        state["qt"], state["rt"], derivative(state["q"]), derivative(state["r"]),
        derivative(state["beta"]),
    )
    coordinate_flux, proper_work, lapse_gradient = _normal_balance_terms(
        state["N"], state["q"], state["r"], state["beta"],
        density, flux, radial_pressure, angular_pressure,
        radial_curvature, angular_curvature, derivative(state["N"]),
    )
    source = proper_work + lapse_gradient
    angle = 2.0 * np.pi * grid["x"] / grid["length"]
    window = 1.0 + 0.8 * np.cos(angle - 1.1)
    window_x = derivative(window)
    normal_sample = volume * density
    fixed_boundary = _integral(grid, window_x * coordinate_flux)
    fixed_volume = _integral(grid, window * source)
    fixed_from_divergence = _integral(grid, window * (-derivative(coordinate_flux) + source))
    assert abs(fixed_from_divergence - (fixed_boundary + fixed_volume)) < 1e-9

    comoving_flux = volume * state["N"] * flux / state["q"]
    shift_advection = state["beta"] * normal_sample
    assert np.max(np.abs(coordinate_flux - (comoving_flux - shift_advection))) < 1e-12
    assert np.max(np.abs(shift_advection)) > 1e-2
    dragged = _integral(grid, window_x * (coordinate_flux + shift_advection))
    orthonormal = _integral(grid, window_x * comoving_flux)
    assert abs(dragged - orthonormal) < 1e-9
    assert abs(dragged - fixed_boundary) > 1e-2
    assert np.max(np.abs(comoving_flux - FOUR_PI * state["r"]**2 * state["N"] * flux)) < 1e-12
    # ∂_t W = β ∂_x W is the scalar transport for point velocity dx/dt = −β.
    dragged_rate = _integral(grid, shift_advection * window_x) + fixed_from_divergence
    expected_dragged = orthonormal + fixed_volume
    assert abs(dragged_rate - expected_dragged) < 1e-9


def test_expansion_free_shift_leaves_coordinate_power_but_not_pressure_work():
    """Q_t = (β Q)_x and r_t = β r_x make K vanish. FQ Q_t is still chart power."""
    grid = _periodic_operators()
    angle = 2.0 * np.pi * grid["x"] / grid["length"]
    length_density = 1.2 + 0.15 * np.cos(angle)
    radial_density = 0.85 + 0.12 * np.sin(angle)
    radius = 2.0 + 0.25 * np.sin(angle)
    shift = 0.35 * np.cos(angle)
    derivative = lambda values: _real_derivative(grid, values)
    radial_t = derivative(shift * radial_density)
    radius_t = shift * derivative(radius)
    lapse = radius * length_density
    radial = radius * radial_density
    radial_velocity = radius_t * radial_density + radius * radial_t
    radial_curvature, angular_curvature = _expansion(
        lapse, radial, radius, shift, radial_velocity, radius_t,
        derivative(radial), derivative(radius), derivative(shift),
    )
    assert np.max(np.abs(radial_curvature)) < 1e-11
    assert np.max(np.abs(angular_curvature)) < 1e-11
    assert np.max(np.abs(radial_t - derivative(shift) * radial_density - shift * derivative(radial_density))) < 1e-11

    first, second = _deterministic_spinor(grid)
    samples = _chart_samples(grid, first, second, length_density, radial_density, shift)
    coordinate_power = samples["FQ"] * radial_t
    _density, flux, radial_pressure, angular_pressure, volume = _observer_from_physical(
        lapse, radial, radius,
        samples["FL"] / radius, samples["FQ"] / radius, 
        -(length_density * samples["FL"] + radial_density * samples["FQ"]) / radius,
        samples["FB"],
    )
    _coordinate_flux, proper_work, lapse_gradient = _normal_balance_terms(
        lapse, radial, radius, shift, _density, flux, radial_pressure, angular_pressure,
        radial_curvature, angular_curvature, derivative(lapse),
    )
    assert np.max(np.abs(proper_work)) < 1e-10
    assert np.max(np.abs(coordinate_power)) > 0.5
    assert abs(_integral(grid, coordinate_power)) > 1.0
    assert abs(_integral(grid, lapse_gradient)) > 0.1
    assert abs(_integral(grid, coordinate_power) - _integral(grid, lapse_gradient)) > 1.0

    step = 1e-6
    zeros = np.zeros_like(radial_density)

    def rates(direction):
        evolved0, evolved1 = _advance_spinor(
            grid, first, second, length_density, radial_density, shift,
            zeros, radial_t, zeros, direction * step,
        )
        moved_q = radial_density + direction * step * radial_t
        moved_r = radius + direction * step * radius_t
        evolved = _chart_samples(grid, evolved0, evolved1, length_density, moved_q, shift)
        lapse_now = moved_r * length_density
        radial_now = moved_r * moved_q
        density_now, _flux_now, _pr, _pp, volume_now = _observer_from_physical(
            lapse_now, radial_now, moved_r,
            evolved["FL"] / moved_r, evolved["FQ"] / moved_r,
            -(length_density * evolved["FL"] + moved_q * evolved["FQ"]) / moved_r,
            evolved["FB"],
        )
        return _integral(grid, evolved["energy"]), _integral(grid, volume_now * density_now)

    forward = rates(1.0)
    backward = rates(-1.0)
    hamiltonian_rate = (forward[0] - backward[0]) / (2.0 * step)
    normal_rate = (forward[1] - backward[1]) / (2.0 * step)
    assert hamiltonian_rate == pytest.approx(_integral(grid, coordinate_power), abs=1e-6)
    assert normal_rate == pytest.approx(_integral(grid, lapse_gradient), abs=1e-6)
    assert proper_work.sum() * grid["dx"] == pytest.approx(0.0, abs=1e-8)


def test_chart_flux_uses_every_explicit_hamiltonian_velocity():
    """∂_t e + ∂_x Φ = FL L_t + FQ Q_t + Fβ β_t. Dropping L_t or β_t leaves their power."""
    grid = _periodic_operators()
    angle = 2.0 * np.pi * grid["x"] / grid["length"]
    length_density = 1.15 + 0.1 * np.cos(angle)
    radial_density = 0.9 + 0.08 * np.sin(angle)
    shift = 0.22 * np.cos(angle) + 0.04 * np.sin(2.0 * angle)
    length_t = 0.02 * np.sin(angle)
    radial_t = 0.05 * np.cos(angle)
    shift_t = -0.03 * np.cos(2.0 * angle)
    first, second = _deterministic_spinor(grid)
    samples = _chart_samples(grid, first, second, length_density, radial_density, shift)
    step = 1e-6

    def density(direction):
        evolved0, evolved1 = _advance_spinor(
            grid, first, second, length_density, radial_density, shift,
            length_t, radial_t, shift_t, direction * step,
        )
        return _chart_samples(
            grid, evolved0, evolved1,
            length_density + direction * step * length_t,
            radial_density + direction * step * radial_t,
            shift + direction * step * shift_t,
        )["energy"]

    energy_t = (density(1.0) - density(-1.0)) / (2.0 * step)
    full_power = (
        samples["FL"] * length_t + samples["FQ"] * radial_t + samples["FB"] * shift_t
    )
    residual = energy_t + _real_derivative(grid, samples["flux"]) - full_power
    radial_only = energy_t + _real_derivative(grid, samples["flux"]) - samples["FQ"] * radial_t
    explicit_gauge = samples["FL"] * length_t + samples["FB"] * shift_t
    assert np.max(np.abs(residual)) < 1e-6
    assert np.max(np.abs(radial_only - explicit_gauge)) < 1e-6
    assert np.max(np.abs(explicit_gauge)) > 1e-3


def test_projected_spinor_rate_exchanges_locally_beyond_chart_work():
    """A mode cutoff on ψ̇ leaves a local flux residual and almost no global source."""
    grid = _periodic_operators()
    angle = 2.0 * np.pi * grid["x"] / grid["length"]
    length_density = 1.15 + 0.1 * np.cos(angle)
    radial_density = 0.9 + 0.08 * np.sin(angle)
    shift = 0.22 * np.cos(angle)
    first, second = _deterministic_spinor(grid)
    samples = _chart_samples(grid, first, second, length_density, radial_density, shift)
    image0, image1 = _dirac_images(grid, first, second, length_density, radial_density, shift)
    full0, full1 = -1j * image0, -1j * image1
    projected0 = grid["project"](full0, 1)
    projected1 = grid["project"](full1, 1)
    assert np.max(np.abs(projected0 - full0)) > 1e-3

    def energy_rate(rate0, rate1):
        step = 1e-7
        forward = _chart_samples(
            grid, first + step * rate0, second + step * rate1,
            length_density, radial_density, shift,
        )["energy"]
        backward = _chart_samples(
            grid, first - step * rate0, second - step * rate1,
            length_density, radial_density, shift,
        )["energy"]
        return (forward - backward) / (2.0 * step)

    divergence = _real_derivative(grid, samples["flux"])
    full_residual = energy_rate(full0, full1) + divergence
    projected_residual = energy_rate(projected0, projected1) + divergence
    assert np.max(np.abs(full_residual)) < 1e-5
    assert np.max(np.abs(projected_residual)) > 1.0
    assert abs(_integral(grid, projected_residual)) < 1e-6


def test_commutator_defect_is_canonical_power_outside_the_metric_variation():
    """Tr(Ċ H) survives when Ċ is not −i[H, C]. It is not Tr(C ∂_t H)."""
    hamiltonian = np.array([[0.4, 0.2 - 0.1j], [0.2 + 0.1j, -0.3]], dtype=complex)
    hamiltonian_t = np.array([[0.05, 0.1], [0.1, -0.02]], dtype=complex)
    defect = np.array([[0.0, 0.15], [0.15, 0.0]], dtype=complex)
    state_vector = np.array([0.6, 0.8j], dtype=complex)
    state_vector /= np.linalg.norm(state_vector)
    covariance = np.outer(state_vector, state_vector.conj())
    commutator = -1j * (hamiltonian @ covariance - covariance @ hamiltonian)
    assert np.allclose(commutator, commutator.conj().T)
    rate = commutator + defect
    geometric = float(np.trace(covariance @ hamiltonian_t).real)
    canonical_defect = float(np.trace(rate @ hamiltonian).real)
    assert abs(float(np.trace(commutator @ hamiltonian).real)) < 1e-12
    assert abs(canonical_defect) > 1e-3
    step = 1e-6

    def energy(direction):
        moved = covariance + direction * step * rate
        moved = 0.5 * (moved + moved.conj().T)
        changed = hamiltonian + direction * step * hamiltonian_t
        return float(np.trace(moved @ changed).real)

    finite_difference = (energy(1.0) - energy(-1.0)) / (2.0 * step)
    assert finite_difference == pytest.approx(geometric + canonical_defect, abs=1e-8)
    assert finite_difference != pytest.approx(geometric, abs=1e-4)


def test_owner_nodal_partials_use_this_dictionary_without_another_factor():
    """conformal_source partials match the trace and the parent conformal map."""
    metric = smooth_metric(12, general=True)
    angle = 2.0 * np.pi * metric.x / metric.length
    shift = 0.08 * np.cos(angle)
    state = fixed_gaussian_covariance(2 * metric.points)
    source = conformal_source(metric, shift, state, 1)
    spacing = source["spacing"]
    radius = metric.sphere_radius
    length_density = metric.lapse / radius
    radial_density = metric.radial_scale / radius
    force_l = source["nodal"]["L"]
    force_q = source["nodal"]["Q"]
    force_beta = source["nodal"]["beta"]
    force_n = source["nodal"]["N"]
    force_radial = source["nodal"]["q"]
    force_radius = source["nodal"]["r"]
    assert np.max(np.abs(force_n - force_l / radius)) < 1e-12
    assert np.max(np.abs(force_radial - force_q / radius)) < 1e-12
    assert np.max(np.abs(
        force_radius + (length_density * force_l + radial_density * force_q) / radius
    )) < 1e-12
    parent = _observer_from_conformal(
        length_density, radial_density, radius,
        force_l / spacing, force_q / spacing, force_beta / spacing,
    )
    physical = _observer_from_physical(
        metric.lapse, metric.radial_scale, radius,
        force_n / spacing, force_radial / spacing, force_radius / spacing, force_beta / spacing,
    )
    for parent_value, physical_value in zip(parent, physical[:4]):
        assert np.max(np.abs(parent_value - physical_value)) < 1e-12

    hamiltonian_samples = length_density * force_l + shift * force_beta
    assert sum(hamiltonian_samples) == pytest.approx(source["energy"], abs=1e-9)
    normal_samples = force_n
    assert abs(float(np.sum(normal_samples) - source["energy"])) > 1e-3

    def energy(changed, changed_shift):
        return conformal_source(changed, changed_shift, state, 1)["energy"]

    node = 4
    step = 1e-6

    def central(field):
        samples = []
        for sign in (-1.0, 1.0):
            if field == "L":
                values = length_density.copy()
                values[node] += sign * step
                changed = replace(metric, lapse=values * radius, radial_scale=radial_density * radius)
                beta = shift
            elif field == "Q":
                values = radial_density.copy()
                values[node] += sign * step
                changed = replace(metric, lapse=length_density * radius, radial_scale=values * radius)
                beta = shift
            elif field == "beta":
                changed = metric
                beta = np.array(shift, dtype=float, copy=True)
                beta[node] += sign * step
            elif field == "r":
                values = radius.copy()
                values[node] += sign * step
                changed = replace(metric, sphere_radius=values)
                beta = shift
            else:
                raise AssertionError(field)
            samples.append(energy(changed, beta))
        return (samples[1] - samples[0]) / (2.0 * step)

    for field in ("L", "Q", "beta", "r"):
        assert central(field) == pytest.approx(source["nodal"][field][node], rel=1e-6, abs=1e-7)
