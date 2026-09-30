"""Collocation constraint algebra on resolved smooth data.

These checks use the discrete Hamiltonian vector field. They do not load a
saved trajectory and they do not certify a later Galerkin discretization.
The mode-31 case records raw-collocation aliasing: the same physical wave
fails on 64 points and passes once that wave sits below Nyquist.
"""
from dataclasses import replace

import numpy as np

from recursive_horizons import nsc_spherical_coupling as coupling


_ABS = 1e-9
_FOURIER = 1e-11


def _gauge(system, lapse, shift, occupations=None):
    payload = {
        "length_density": np.asarray(lapse, dtype=float),
        "shift": np.asarray(shift, dtype=float),
    }
    if occupations is not None:
        payload["occupations"] = np.asarray(occupations, dtype=float)
    return replace(system, **payload)


def _zeros(system):
    return np.zeros(system.points)


def _manufactured(system, bare, q_variation=0.0):
    xi = system.xi
    length = system.length
    q0 = float(bare.Q[0])
    state = bare.copy()
    state.Q = q0 * (
        1.0
        + q_variation * (0.04 * np.cos(2 * np.pi * xi / length) + 0.02 * np.sin(4 * np.pi * xi / length))
    )
    state.r = 2.0 + 0.05 * np.cos(2 * np.pi * xi / length)
    state.chi = 0.03 * np.sin(2 * np.pi * xi / length)
    state.p_Q = 0.02 * np.cos(4 * np.pi * xi / length)
    state.p_r = 0.01 * np.sin(2 * np.pi * xi / length)
    state.p_chi = 0.012 * np.cos(2 * np.pi * xi / length)
    return state


def _ap_columns(system):
    xi = system.xi
    length = system.length
    columns0 = []
    columns1 = []
    occupations = []
    for harmonic, upper, lower, weight in ((0, 0.6, 0.8, 0.5), (1, -0.3, 0.4, 0.25)):
        wavenumber = 2 * np.pi * (harmonic + 0.5) / length
        wave = np.exp(1j * wavenumber * xi) / np.sqrt(system.points)
        first = upper * wave
        second = 1j * lower * wave
        norm = np.sqrt(np.vdot(first, first) + np.vdot(second, second))
        columns0.append(first / norm)
        columns1.append(second / norm)
        occupations.append(weight)
    return np.column_stack(columns0), np.column_stack(columns1), np.asarray(occupations, dtype=float)


def _images(system, lapse, shift, state):
    return coupling.apply_dirac(
        state.phi0,
        state.phi1,
        np.asarray(lapse, dtype=float),
        state.Q,
        np.asarray(shift, dtype=float),
        system.kappa,
        system.momentum,
    )


def _matter_bracket(system, lapse_n, lapse_m, state, occupations):
    shift = _zeros(system)
    left = _images(system, lapse_n, shift, state)
    right = _images(system, lapse_m, shift, state)
    overlap = np.sum(
        np.conjugate(left[0]) * right[0] + np.conjugate(left[1]) * right[1],
        axis=0,
    )
    # One angular weight. Occupations stay inside the Gaussian sum.
    return 2.0 * system.multiplicity * np.sum(occupations * np.imag(overlap))


def _geometric_bracket(system, state, lapse_n, lapse_m, occupations=None, matter_force=False):
    shift = _zeros(system)
    left = list(coupling.geometric_rates(_gauge(system, lapse_n, shift), state))
    right = list(coupling.geometric_rates(_gauge(system, lapse_m, shift), state))
    if matter_force:
        source_n = coupling.source_from_columns(
            _gauge(system, lapse_n, shift, occupations), state,
        )
        source_m = coupling.source_from_columns(
            _gauge(system, lapse_m, shift, occupations), state,
        )
        left[3] = left[3] - source_n["force_Q"] / system.dx
        right[3] = right[3] - source_m["force_Q"] / system.dx
    q_left, p_left = np.stack(left[:3]), np.stack(left[3:])
    q_right, p_right = np.stack(right[:3]), np.stack(right[3:])
    return system.dx * np.sum(q_left * p_right - p_left * q_right)


def _positive_smearings(system):
    xi = system.xi
    length = system.length
    lapse_n = 1.0 + 0.2 * np.cos(2 * np.pi * xi / length)
    lapse_m = (
        1.15
        + 0.08 * np.sin(2 * np.pi * xi / length)
        + 0.03 * np.cos(4 * np.pi * xi / length)
    )
    return lapse_n, lapse_m


def _shift_vector(system, lapse_n, lapse_m, radial_density):
    derivative = system.derivative
    return (
        lapse_n * (derivative @ lapse_m) - lapse_m * (derivative @ lapse_n)
    ) / radial_density ** 2


def test_periodic_and_antiperiodic_fourier_symbols():
    points = 32
    length = coupling.PERIOD
    derivative = coupling.periodic_derivative(points, length)
    xi = np.arange(points) * (length / points)
    for mode in (0, 1, 2, -3, 5):
        wave = np.exp(2j * np.pi * mode * xi / length)
        wavenumber = 2 * np.pi * mode / length
        mismatch = np.max(np.abs(derivative @ wave - 1j * wavenumber * wave))
        assert mismatch < _FOURIER
    nyquist = np.array([(-1.0) ** index for index in range(points)])
    assert np.max(np.abs(derivative @ nyquist)) < _FOURIER

    momentum, _metric = coupling.antiperiodic_momentum(points, length)
    for harmonic in (0, 1, -1, 2):
        wavenumber = 2 * np.pi * (harmonic + 0.5) / length
        wave = np.exp(1j * wavenumber * xi)
        mismatch = np.max(np.abs(momentum @ wave - wavenumber * wave))
        assert mismatch < _FOURIER


def test_low_band_geometry_matter_and_combined_brackets_use_one_angular_weight():
    system, bare = coupling.build_system(64)
    assert system.multiplicity == 4 * int(system.kappa)
    assert system.multiplicity == 4
    lapse_n, lapse_m = _positive_smearings(system)
    assert np.min(lapse_n) > 0.0 and np.min(lapse_m) > 0.0
    phi0, phi1, occupations = _ap_columns(system)

    for variation in (0.0, 0.15):
        state = _manufactured(system, bare, variation)
        state.phi0 = phi0
        state.phi1 = phi1
        generated = _shift_vector(system, lapse_n, lapse_m, state.Q)
        geometry = _geometric_bracket(system, state, lapse_n, lapse_m)
        diffeomorphism = system.dx * np.sum(
            generated * coupling.shift_constraint(system, state)
        )
        assert abs(geometry - diffeomorphism) < _ABS

        matter = _matter_bracket(system, lapse_n, lapse_m, state, occupations)
        source = coupling.source_from_columns(
            _gauge(system, lapse_n, _zeros(system), occupations),
            state,
        )
        matter_constraint = float(np.sum(generated * source["force_beta"]))
        assert abs(matter - matter_constraint) < _ABS
        doubled = matter * system.multiplicity
        # A second factor of M moves the bracket by the whole matter scale.
        assert abs(doubled - matter_constraint) > 0.5 * max(abs(matter_constraint), 1.0)

        combined = _geometric_bracket(
            system, state, lapse_n, lapse_m, occupations, matter_force=True,
        )
        assert abs((combined + matter) - (diffeomorphism + matter_constraint)) < _ABS


def test_beta_only_lie_laws_give_c_weight_one_and_d_weight_two():
    system, bare = coupling.build_system(64)
    xi = system.xi
    length = system.length
    derivative = system.derivative
    q0 = float(bare.Q[0])
    state = bare.copy()
    state.Q = q0 * (
        1.0 + 0.05 * np.cos(2 * np.pi * xi / length) + 0.02 * np.sin(4 * np.pi * xi / length)
    )
    state.r = 2.0 + 0.08 * np.cos(2 * np.pi * xi / length) + 0.03 * np.sin(4 * np.pi * xi / length)
    state.chi = 0.04 * np.sin(2 * np.pi * xi / length) + 0.01 * np.cos(6 * np.pi * xi / length)
    state.p_Q = 0.02 * np.cos(4 * np.pi * xi / length)
    state.p_r = 0.015 * np.sin(2 * np.pi * xi / length)
    state.p_chi = 0.01 * np.cos(2 * np.pi * xi / length)
    beta = 0.3 + 0.1 * np.cos(2 * np.pi * xi / length) + 0.05 * np.sin(4 * np.pi * xi / length)
    alpha = 0.2 + 0.07 * np.sin(4 * np.pi * xi / length)
    lapse = 1.0 + 0.2 * np.sin(2 * np.pi * xi / length) + 0.04 * np.cos(4 * np.pi * xi / length)
    zero = _zeros(system)

    transport = coupling.geometric_rates(_gauge(system, zero, beta), state)
    expected = (
        derivative @ (beta * state.Q),
        beta * (derivative @ state.r),
        beta * (derivative @ state.chi),
        beta * (derivative @ state.p_Q),
        derivative @ (beta * state.p_r),
        derivative @ (beta * state.p_chi),
    )
    for actual, target in zip(transport, expected):
        assert np.max(np.abs(actual - target)) < _ABS

    hamilton_flow = coupling.geometric_rates(_gauge(system, lapse, zero), state)
    diffeomorphism_flow = transport
    bracket = system.dx * np.sum(
        np.stack(hamilton_flow[:3]) * np.stack(diffeomorphism_flow[3:])
        - np.stack(hamilton_flow[3:]) * np.stack(diffeomorphism_flow[:3])
    )
    density = coupling.hamilton_constraint(system, state)
    smeared = system.dx * np.sum((-beta * (derivative @ lapse)) * density)
    assert abs(bracket - smeared) < _ABS

    other = coupling.geometric_rates(_gauge(system, zero, alpha), state)
    algebra = system.dx * np.sum(
        np.stack(transport[:3]) * np.stack(other[3:])
        - np.stack(transport[3:]) * np.stack(other[:3])
    )
    commutator = beta * (derivative @ alpha) - alpha * (derivative @ beta)
    momentum = coupling.shift_constraint(system, state)
    assert abs(algebra - system.dx * np.sum(commutator * momentum)) < _ABS

    step = 1e-4

    def advance(dt):
        moved = state.copy()
        moved.Q = state.Q + dt * transport[0]
        moved.r = state.r + dt * transport[1]
        moved.chi = state.chi + dt * transport[2]
        moved.p_Q = state.p_Q + dt * transport[3]
        moved.p_r = state.p_r + dt * transport[4]
        moved.p_chi = state.p_chi + dt * transport[5]
        return moved

    hamilton_rate = (
        coupling.hamilton_constraint(system, advance(step))
        - coupling.hamilton_constraint(system, advance(-step))
    ) / (2 * step)
    momentum_rate = (
        coupling.shift_constraint(system, advance(step))
        - coupling.shift_constraint(system, advance(-step))
    ) / (2 * step)
    weight_one = derivative @ (beta * density)
    weight_two = derivative @ (beta * momentum) + (derivative @ beta) * momentum
    assert np.max(np.abs(hamilton_rate - weight_one)) < 1e-6
    assert np.max(np.abs(hamilton_rate - beta * (derivative @ density))) > 1e-2
    assert np.max(np.abs(momentum_rate - weight_two)) < 1e-8
    assert np.max(np.abs(momentum_rate - (derivative @ (beta * momentum)))) > 1e-4


def _mode_bracket_gap(points, mode):
    system, bare = coupling.build_system(points)
    xi = system.xi
    length = system.length
    state = bare.copy()
    state.Q = np.full(system.points, float(bare.Q[0]))
    state.r = 2.0 + 0.05 * np.cos(2 * np.pi * mode * xi / length)
    state.chi = 0.02 * np.sin(2 * np.pi * min(mode, 3) * xi / length)
    state.p_Q = 0.01 * np.cos(2 * np.pi * xi / length)
    state.p_r = 0.02 * np.sin(2 * np.pi * mode * xi / length)
    state.p_chi = 0.01 * np.cos(2 * np.pi * xi / length)
    state.phi0 = np.zeros((system.points, 1), dtype=complex)
    state.phi1 = np.zeros((system.points, 1), dtype=complex)
    wave = np.exp(2j * np.pi * mode * xi / length)
    wavenumber = 2 * np.pi * mode / length
    symbol_error = np.max(np.abs(system.derivative @ wave - 1j * wavenumber * wave))
    lapse_n = 1.0 + 0.2 * np.cos(2 * np.pi * xi / length)
    lapse_m = 1.1 + 0.15 * np.sin(4 * np.pi * xi / length)
    generated = _shift_vector(system, lapse_n, lapse_m, state.Q)
    bracket = _geometric_bracket(system, state, lapse_n, lapse_m)
    predicted = system.dx * np.sum(generated * coupling.shift_constraint(system, state))
    return symbol_error, bracket - predicted


def test_raw_collocation_mode31_aliases_on_n64_and_resolves_on_n128():
    """Same physical mode 31. N=64 raw products alias; N=128 resolves them.

    This is a property of the current collocation products, not a trajectory
    requirement and not a failure demanded of a Galerkin replacement.
    """
    coarse_symbol, coarse_gap = _mode_bracket_gap(64, 31)
    fine_symbol, fine_gap = _mode_bracket_gap(128, 31)
    assert coarse_symbol < _FOURIER
    assert fine_symbol < _FOURIER
    assert abs(coarse_gap) > 1.0
    assert abs(fine_gap) < _ABS
