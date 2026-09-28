"""Synthetic KS-envelope controls; no authenticated field or 7601-grid runs."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeIncoming, axial_center, computational_z_grid,
    continuum_speed_distance, envelope_z_derivative, evolve_ks_source_envelope,
    frame_potential_identity, physical_incoming_interval,
    trigonometric_polynomial, usual_axial_support,
)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_retarded_radial_response import ks_generator
from recursive_horizons.nsc_transmitting_dirac_domain import S2, S3


MASS = np.pi / 2
ANGULAR = np.sqrt(5)
RHO_UP = 1.03
SOLVER = dict(rtol=2e-10, atol=2e-12, max_step=0.002)


def tiny_source(energies=(0.55, -0.4), weights=(1.1, 0.8)):
    blocks = np.array([
        [[0.45, 0.06 + 0.03j], [0.06 - 0.03j, 0.5]],
        [[0.55, -0.04j], [0.04j, 0.4]],
    ], complex)
    return FixedSourcePreparation.from_signed_blocks(
        np.asarray(energies, float), np.asarray(weights, float), blocks)


def one_source(energy=0.7):
    return FixedSourcePreparation.from_signed_blocks(
        np.array([energy]), np.array([1.2]), np.array([[[0.8]]]))


def A_up(nsrc, seed=2):
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(2, nsrc)) + 1j * rng.normal(size=(2, nsrc))) * 0.2


def zero_family():
    return LocalIncomingFamily(np.zeros((2, 8)))


def bump_family(amplitude=0.002):
    center = axial_center()
    w = LocalAxialFunction((0.05, -0.02, 0.01), center)
    U = LocalAxialFunction((0.02, -0.005), center)
    return CompatibleIncomingMetric(
        (float(amplitude),), (CompatibleRadiusDirection(w, U, 0.007, 0.03),))


def evolve(family, source, initial, z=None, target=None, tangents='zero', **extra):
    z = computational_z_grid(16) if z is None else z
    interval = physical_incoming_interval()
    if target is None:
        target = np.linspace(interval[0], interval[1], 7)
    return evolve_ks_source_envelope(
        source, initial, family, z, target, MASS, ANGULAR, RHO_UP,
        axial_support=usual_axial_support(), tangents=tangents, **SOLVER, **extra)


def test_frame_map_cancels_radius_and_matches_ks_potential():
    report = frame_potential_identity(1.01, MASS, ANGULAR)
    assert report['radius_cancellation'] < 3e-15
    assert report['potential_residual'] < 3e-15
    assert report['depends_only_on_beta_a']
    assert report['a'] == pytest.approx(np.sqrt(report['beta'] ** 2 - 1))


def test_homogeneous_envelope_agrees_with_ks_generator_radial_ode():
    source = tiny_source()
    initial = A_up(4)
    result = evolve(zero_family(), source, initial, tangents='zero')
    phase = np.exp(1j * result.source_energies * result.z[:, None])
    recovered = result.columns * phase[:, None, :]
    assert float(np.max(np.abs(recovered - recovered[:1]))) < 3e-11

    def rhs(rho, state):
        field = state.reshape(2, 4)
        G = ks_generator(rho, source.energies, MASS, ANGULAR)
        return np.einsum('sab,bs->as', G, field).ravel()

    radial = solve_ivp(
        rhs, (RHO_UP, 1.), initial.ravel(), method='DOP853',
        rtol=SOLVER['rtol'], atol=SOLVER['atol'], max_step=SOLVER['max_step'])
    assert radial.success
    np.testing.assert_allclose(
        recovered[0], radial.y[:, -1].reshape(2, 4), atol=3e-11, rtol=0)
    result.require_history(zero_family())
    result.require_history(zero_family().metric())
    assert result.diagnostics['extra_r_z_term'] is False
    assert result.coverage['PG_discretization_copied'] is False
    assert not hasattr(result, 'node_fields')
    assert not hasattr(result, 'history')


def test_fft_envelope_derivative_matches_analytic_plane_wave():
    z = computational_z_grid(32)
    wave = 2 * np.pi * np.fft.fftfreq(len(z), d=float(z[1] - z[0]))
    assert wave.shape == z.shape
    for mode in (1, 2, -3):
        field = np.exp(1j * wave[mode] * (z - z[0]))
        observed = envelope_z_derivative(field, z)
        np.testing.assert_allclose(observed, 1j * wave[mode] * field, atol=1e-11, rtol=1e-13)
        rebuilt = trigonometric_polynomial(field, z, z)
        np.testing.assert_allclose(rebuilt, field, atol=1e-11, rtol=1e-13)
        deriv = trigonometric_polynomial(field, z, z, derivative=True)
        np.testing.assert_allclose(deriv, observed, atol=1e-11, rtol=1e-13)
        mid = 0.5 * (z[4] + z[5])
        exact = np.exp(1j * wave[mode] * (mid - z[0]))
        assert abs(trigonometric_polynomial(field, z, np.array([mid]))[0] - exact) < 1e-11
    spinor = np.stack((field, -0.3j * field))
    dspinor = envelope_z_derivative(spinor, z)
    np.testing.assert_allclose(dspinor[0], 1j * wave[mode] * field, atol=1e-11, rtol=1e-13)


def test_nonzero_radius_finite_difference_includes_F_z_tangent():
    source = one_source()
    initial = A_up(1, seed=5)
    z = computational_z_grid(32)
    interval = physical_incoming_interval()
    target = np.linspace(interval[0], interval[1], 11)
    h = 1e-4
    base = evolve(bump_family(0.002), source, initial, z=z, target=target, tangents='all')
    plus = evolve(bump_family(0.002 + h), source, initial, z=z, target=target, tangents='zero')
    minus = evolve(bump_family(0.002 - h), source, initial, z=z, target=target, tangents='zero')
    dF = (plus.columns - minus.columns) / (2 * h)
    dFz = (plus.axial_columns - minus.axial_columns) / (2 * h)
    np.testing.assert_allclose(dF, base.column_tangents[0], atol=3e-8, rtol=0)
    np.testing.assert_allclose(dFz, base.axial_tangents[0], atol=3e-8, rtol=0)
    interior = slice(1, -1)
    step = float(target[1] - target[0])
    fd_z = (base.columns[2:] - base.columns[:-2]) / (2 * step)
    np.testing.assert_allclose(fd_z, base.axial_columns[interior], atol=2e-6, rtol=0)
    phase = np.exp(-1j * base.source_energies * base.z[:, None])
    recovered_X = base.columns / phase[:, None, :]
    fd_X = (recovered_X[2:] - recovered_X[:-2]) / (2 * step)
    reconstructed = phase[interior, None, :] * (
        fd_X - 1j * base.source_energies * recovered_X[interior])
    np.testing.assert_allclose(reconstructed, base.axial_columns[interior], atol=2e-6, rtol=0)


def test_whole_grid_discrete_operator_is_anti_hermitian_and_preserves_norm():
    z = computational_z_grid(16)
    rho = 1.015
    source = one_source()
    family = bump_family(0.003)
    T1 = chart_coordinates(1.)[0]
    shift = chart_coordinates(rho)[0] - T1
    r_ref = np.sqrt(1 + rho * rho)
    r_g = r_ref + family.amplitudes[0] * np.array(
        [family.directions[0].value(shift, float(point)) for point in z])
    beta = float(geometry(rho)[0])
    axial = np.sqrt(beta * beta - 1)
    G = ks_generator(rho, source.energies, MASS, ANGULAR)

    def apply(field):
        freq = envelope_z_derivative(field, z)
        term_z = np.einsum('ab,bsz->asz', S3, freq) / (axial * axial)
        term_G = np.einsum('sab,bsz->asz', G, field)
        term_r = (1j * ANGULAR / axial) * (1 / r_g - 1 / r_ref) * np.einsum('ab,bsz->asz', S2, field)
        return term_z + term_G + term_r

    n = 2 * len(z)
    matrix = np.zeros((n, n), complex)
    for index in range(n):
        unit = np.zeros(n, complex)
        unit[index] = 1
        matrix[:, index] = apply(unit.reshape(2, 1, len(z))).ravel()
    assert float(np.max(np.abs(matrix + matrix.conj().T))) < 3e-14
    rng = np.random.default_rng(3)
    sample = rng.normal(size=(2, 1, len(z))) + 1j * rng.normal(size=(2, 1, len(z)))
    image = apply(sample)
    pairing = np.vdot(sample.ravel(), image.ravel())
    assert abs(pairing.real) < 3e-12
    result = evolve(family, source, A_up(1, seed=8), z=z, tangents='zero')
    assert result.diagnostics['computational_norm_residual'] < 3e-11


def test_source_frame_constant_initial_phase_weights_and_coherence():
    source = tiny_source()
    assert source.energies.tolist() == [0.55, 0.55, -0.4, -0.4]
    initial = A_up(4, seed=9)
    initial[:, 1] = initial[:, 0]
    result = evolve(zero_family(), source, initial, tangents='zero')
    np.testing.assert_allclose(result.columns[:, :, 0], result.columns[:, :, 1], atol=3e-12, rtol=0)
    phase = np.exp(1j * result.source_energies * result.z[:, None])
    recovered = result.columns * phase[:, None, :]
    assert float(np.max(np.abs(recovered - recovered[:1]))) < 3e-11
    np.testing.assert_allclose(
        result.weighted_columns, result.columns * result.column_weights, atol=0, rtol=0)
    np.testing.assert_allclose(
        result.weighted_axial_columns, result.axial_columns * result.column_weights, atol=0, rtol=0)
    kernel = result.covariance()
    assert float(np.max(np.abs(kernel - kernel.transpose(2, 3, 0, 1).conj()))) < 3e-14
    assert result.covariance_tangent().shape[0] == 0
    assert result.coverage['source_frequency_controls_spatial_mesh'] is False
    assert result.coverage['quadrature_applied_once'] is True


def test_padding_preparation_radius_and_history_guards():
    source = one_source()
    initial = A_up(1)
    family = bump_family(0.002)
    z16 = computational_z_grid(16)
    z32 = computational_z_grid(32)
    interval = physical_incoming_interval()
    target = np.linspace(interval[0], interval[1], 5)
    result = evolve(family, source, initial, z=z16, target=target, tangents='all')
    result.require_history(family, z_grid=z16, rho_up=RHO_UP, solver_options=SOLVER)
    clone = CompatibleIncomingMetric(
        family.amplitudes,
        (CompatibleRadiusDirection(
            lambda point, order, _w=result.binding.w_samples[0], _z=z16:
            0. if order else float(_w[np.argmin(np.abs(_z - point))]),
            lambda point, order, _U=result.binding.U_samples[0], _z=z16:
            0. if order else float(_U[np.argmin(np.abs(_z - point))]),
            family.directions[0].inner_radius, family.directions[0].outer_radius),))
    result.require_history(clone)
    with pytest.raises(ValueError, match='generating amplitudes'):
        result.require_history(bump_family(0.003))
    warped = CompatibleIncomingMetric(
        family.amplitudes,
        (CompatibleRadiusDirection(
            LocalAxialFunction((0.08, 0.01), axial_center()),
            LocalAxialFunction((0.,), axial_center()),
            0.007, 0.03),))
    with pytest.raises(ValueError, match='evaluated w/U profiles'):
        result.require_history(warped)
    with pytest.raises(ValueError, match='computational z-grid'):
        result.require_history(family, z_grid=z32)
    with pytest.raises(ValueError, match='rho endpoints'):
        result.require_history(family, rho_up=1.04)
    with pytest.raises(ValueError, match='solver options'):
        result.require_history(family, solver_options=dict(rtol=1e-8, atol=1e-10, max_step=0.01))
    altered = result.binding.w_samples.copy()
    altered[0, 4] += 0.02
    with pytest.raises(ValueError, match='fingerprint'):
        replace(result, binding=replace(result.binding, w_samples=altered))
    other = evolve(family, source, initial, z=z32, target=target, tangents='zero')
    assert other.fixed_preparation_digest == result.fixed_preparation_digest
    assert other.binding.fingerprint != result.binding.fingerprint
    D = continuum_speed_distance(RHO_UP)
    assert result.diagnostics['padding'][0] > 2 * D
    with pytest.raises(ValueError, match='rho_up'):
        evolve_ks_source_envelope(
            source, initial, family, z16, target, MASS, ANGULAR, 1.02,
            axial_support=usual_axial_support(), **SOLVER)
    with pytest.raises(ValueError, match='padding'):
        evolve_ks_source_envelope(
            source, initial, family, z16, target, MASS, ANGULAR, 1.05,
            axial_support=usual_axial_support(), **SOLVER)
    with pytest.raises(ValueError, match='radius'):
        evolve(LocalIncomingFamily(np.array([[100.] + [0.] * 7, [0.] * 8])), source, initial)
    wide = CompatibleIncomingMetric(
        family.amplitudes, (CompatibleRadiusDirection(
            family.directions[0].w, family.directions[0].U, 0.007, 0.04),))
    with pytest.raises(ValueError, match='normal support'):
        evolve(wide, source, initial)
    with pytest.raises(ValueError, match='physical interval'):
        evolve(family, source, initial, target=np.array([interval[0] - 0.02]))
    with pytest.raises(TypeError, match='FixedSourcePreparation'):
        evolve(family, np.eye(1), initial)
    with pytest.raises(ValueError, match='node_fields'):
        evolve(family, source, initial, node_fields=np.ones(3))
    with pytest.raises(ValueError, match='HistoryBinding'):
        evolve(family, source, initial, history=family)
    with pytest.raises(TypeError, match='LocalIncomingFamily or CompatibleIncomingMetric'):
        evolve_ks_source_envelope(
            source, initial, object(), z16, target, MASS, ANGULAR, RHO_UP,
            axial_support=usual_axial_support(), **SOLVER)
    assert isinstance(result, KSEnvelopeIncoming)
    assert result.coverage['physical_confidence_from_tolerances'] is False
    assert result.diagnostics['rtol'] == SOLVER['rtol']
    assert result.diagnostics['atol'] == SOLVER['atol']
    assert result.diagnostics['max_step'] == SOLVER['max_step']
