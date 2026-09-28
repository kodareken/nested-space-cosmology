"""Focused checks of the reference/difference coordinate change."""
from dataclasses import replace

import numpy as np
import pytest

from recursive_horizons import nsc_ks_source_envelope as E
from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric, CompatibleRadiusDirection
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_difference_envelope import difference_rhs, evolve_ks_difference_envelope
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction


def inputs(amplitude=0.002):
    source = FixedSourcePreparation.from_signed_blocks(
        np.array([0.7]), np.array([1.2]), np.array([[[0.6, 0.04j], [-0.04j, 0.4]]]))
    initial = np.array([[0.6, 0.3j], [0.2j, 0.7]], complex)
    w = LocalAxialFunction((0., 1.), E.axial_center())
    U = LocalAxialFunction((0.2, -0.1), E.axial_center())
    family = CompatibleIncomingMetric((amplitude,), (CompatibleRadiusDirection(w, U, .007, .03),))
    return source, initial, family


def evolve(amplitude=0.002, owner=evolve_ks_difference_envelope, tangents='all'):
    return owner(*inputs(amplitude), E.computational_z_grid(32),
                 np.linspace(*E.physical_incoming_interval(), 9), np.pi / 2, np.sqrt(5), 1.03,
                 axial_support=E.usual_axial_support(), rtol=2e-12, atol=2e-14,
                 max_step=.001, tangents=tangents)


def test_sum_equation_is_the_original_operator_and_retarded_variation():
    rng = np.random.default_rng(11)
    z = E.computational_z_grid(16)
    wave = 2 * np.pi * np.fft.fftfreq(len(z), d=z[1] - z[0])
    A = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    D = .01 * (rng.normal(size=(2, 2, 16)) + 1j * rng.normal(size=(2, 2, 16)))
    Y = .01 * (rng.normal(size=(1, 2, 2, 16)) + 1j * rng.normal(size=(1, 2, 2, 16)))
    rho, angular, mass = 1.015, np.sqrt(5), np.pi / 2
    axial = E.reference_chart(rho)[2]
    r_ref = np.sqrt(1 + rho**2)
    basis = .01 * np.sin(2 * np.pi * np.arange(16) / 16)[None]
    delta_r = .2 * basis[0]
    G = E.ks_generator(rho, np.array([.7, -.4]), mass, angular)
    Adot, Ddot, Ydot = difference_rhs(A, D, Y, G, axial, angular, r_ref, delta_r, basis, wave)
    total = A[:, :, None] + D

    def original(field):
        return (np.einsum('ab,...bsz->...asz', E.S3, E.envelope_z_derivative(field, z)) / axial**2
                + np.einsum('sab,...bsz->...asz', G, field)
                + (1j * angular / axial) * (1 / (r_ref + delta_r) - 1 / r_ref)
                * np.einsum('ab,...bsz->...asz', E.S2, field))

    np.testing.assert_allclose(Adot[:, :, None] + Ddot, original(total), rtol=0, atol=5e-13)
    dL = -(1j * angular / axial) * basis / (r_ref + delta_r)**2
    expected = original(Y) + dL[:, None, None, :] * np.einsum('ab,bsz->asz', E.S2, total)[None]
    np.testing.assert_allclose(Ydot, expected, rtol=0, atol=5e-13)


def test_zero_history_has_exactly_zero_difference_and_no_fft_background_derivative():
    result = evolve(0., tangents='zero')
    assert np.count_nonzero(result.envelope_difference) == 0
    assert np.count_nonzero(result.envelope_difference_z) == 0
    np.testing.assert_allclose(result.axial_columns, -1j * result.source_energies * result.columns,
                               rtol=0, atol=2e-16)
    assert result.diagnostics['stress_drift_subtracted'] is False
    assert result.coverage['physical_source_error'] == 'OPEN'


def test_nonzero_evolution_matches_total_coordinates_and_preserves_fixed_source():
    result = evolve()
    old = evolve(owner=E.evolve_ks_source_envelope)
    assert result.fixed_preparation_digest == old.fixed_preparation_digest
    for name in ('columns', 'axial_columns', 'column_tangents', 'axial_tangents'):
        np.testing.assert_allclose(getattr(result, name), getattr(old, name), atol=3e-11, rtol=0)
    assert np.max(abs(result.envelope_difference)) > 1e-9
    with pytest.raises(ValueError, match='reconstruct'):
        replace(result, reference_amplitudes=result.reference_amplitudes + 1e-5)


def test_retarded_tangent_matches_changed_history_and_axial_derivative():
    h = 1e-4
    base, plus, minus = evolve(), evolve(.002 + h, tangents='zero'), evolve(.002 - h, tangents='zero')
    for name, derivative in (('columns', 'column_tangents'), ('axial_columns', 'axial_tangents')):
        fd = (getattr(plus, name) - getattr(minus, name)) / (2 * h)
        np.testing.assert_allclose(fd, getattr(base, derivative)[0], atol=3e-9, rtol=0)
    np.testing.assert_allclose((plus.covariance() - minus.covariance()) / (2 * h),
                               base.covariance_tangent()[0], atol=3e-9, rtol=0)
    assert base.fixed_preparation_digest == plus.fixed_preparation_digest == minus.fixed_preparation_digest
