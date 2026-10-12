"""Opposite-angular KS signed-state map; synthetic only, no source campaigns."""
import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_signed_state import (
    SOURCE_COMPLEMENT_TOLERANCE, actual_radius_envelope_rhs,
    negative_angular_partner, s3_conjugate, source_complement_residual,
)
from recursive_horizons.nsc_ks_source_envelope import (
    axial_center, computational_z_grid, envelope_z_derivative,
    physical_incoming_interval, usual_axial_support,
)
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_retarded_radial_response import ks_generator
from recursive_horizons.nsc_transmitting_dirac_domain import S2, S3


CONFIG = {'surface_gravity': 0.24, 'omega': 3.9}
MASS = np.pi / 2
ANGULAR = np.sqrt(5)
RHO_UP = 1.03
SOLVER = dict(rtol=2e-12, atol=2e-14, max_step=0.001)
ENERGY = 0.7
WEIGHT = 1.2


def horizon_blocks(energies):
    return np.array([
        source_covariance(float(energy), CONFIG['surface_gravity'], CONFIG['omega'], MASS)
        for energy in np.asarray(energies, float)])


def signed_sources():
    energies = np.array([ENERGY])
    weights = np.array([WEIGHT])
    positive = FixedSourcePreparation.from_signed_blocks(
        energies, weights, horizon_blocks(energies))
    negative = FixedSourcePreparation.from_signed_blocks(
        -energies, weights, horizon_blocks(-energies))
    return positive, negative


def A_up(nsrc, seed=5):
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(2, nsrc)) + 1j * rng.normal(size=(2, nsrc))) * 0.25


def bump_family(amplitude=0.002):
    w = LocalAxialFunction((0., 1.), axial_center())
    U = LocalAxialFunction((0.2, -0.1), axial_center())
    return CompatibleIncomingMetric(
        (float(amplitude),), (CompatibleRadiusDirection(w, U, 0.007, 0.03),))


def evolve(source, initial, angular, amplitude=0.002, tangents='all'):
    return evolve_ks_difference_envelope(
        source, initial, bump_family(amplitude), computational_z_grid(16),
        np.linspace(*physical_incoming_interval(), 5), MASS, angular, RHO_UP,
        axial_support=usual_axial_support(), tangents=tangents, **SOLVER)


def test_independent_negative_evolution_agrees_in_F_Fz_and_tangents():
    positive_source, negative_source = signed_sources()
    initial = A_up(positive_source.covariance.shape[0])
    prepared = evolve(positive_source, initial, ANGULAR)
    partner = negative_angular_partner(prepared, negative_source)
    independent = evolve(
        negative_source, s3_conjugate(initial), -ANGULAR)
    for name in (
            'columns', 'axial_columns', 'reference_amplitudes',
            'envelope_difference', 'envelope_difference_z'):
        np.testing.assert_allclose(
            getattr(partner, name), getattr(independent, name), atol=3e-11, rtol=0,
            err_msg=name)
    np.testing.assert_allclose(
        partner.column_tangents, independent.column_tangents, atol=2e-10, rtol=0,
        err_msg='column_tangents')
    # dF_z compares two independent Fourier-differentiated tangent ODEs.
    np.testing.assert_allclose(
        partner.axial_tangents, independent.axial_tangents, atol=3e-8, rtol=0,
        err_msg='axial_tangents')
    np.testing.assert_array_equal(partner.source_energies, -prepared.source_energies)
    np.testing.assert_array_equal(partner.column_weights, prepared.column_weights)
    assert partner.angular == pytest.approx(-ANGULAR)
    assert partner.binding.fingerprint == prepared.binding.fingerprint
    assert partner.fixed_preparation_digest == independent.fixed_preparation_digest
    assert partner.diagnostics['energy_folding_factor'] == 1
    assert np.max(np.abs(prepared.envelope_difference)) > 1e-9
    assert np.max(np.abs(partner.envelope_difference)) > 1e-9


def test_full_pde_antiunitary_identity_includes_z_derivative():
    rng = np.random.default_rng(19)
    z = computational_z_grid(16)
    nz = len(z)
    energies = np.array([0.55, 1.1])
    X = rng.normal(size=(nz, 2, 2)) + 1j * rng.normal(size=(nz, 2, 2))
    Xz = np.moveaxis(envelope_z_derivative(np.moveaxis(X, 0, -1), z), -1, 0)
    rho = 1.015
    r_ref = np.sqrt(1. + rho * rho)
    radius = r_ref + 0.03 * np.sin(2 * np.pi * np.arange(nz) / nz)
    rhs = actual_radius_envelope_rhs(X, Xz, rho, energies, MASS, ANGULAR, radius)
    mapped = actual_radius_envelope_rhs(
        s3_conjugate(X), s3_conjugate(Xz), rho, -energies, MASS, -ANGULAR, radius)
    np.testing.assert_allclose(mapped, s3_conjugate(rhs), atol=1e-13, rtol=0)
    G = ks_generator(rho, energies, MASS, ANGULAR)
    a = np.sqrt(float(geometry(rho)[0]) ** 2 - 1.)
    homogeneous = np.einsum('ab,zbs->zas', S3, Xz) / a ** 2 + np.einsum(
        'sab,zbs->zas', G, X)
    np.testing.assert_allclose(
        actual_radius_envelope_rhs(X, Xz, rho, energies, MASS, ANGULAR, r_ref),
        homogeneous, atol=1e-13, rtol=0)
    extra = (1j * ANGULAR / a) * (1. / radius - 1. / r_ref)[:, None, None] * np.einsum(
        'ab,zbs->zas', S2, X)
    np.testing.assert_allclose(rhs, homogeneous + extra, atol=1e-13, rtol=0)


def test_weights_and_coherence_are_retained_once():
    positive_source, negative_source = signed_sources()
    prepared = evolve(positive_source, A_up(3), ANGULAR, tangents='zero')
    partner = negative_angular_partner(prepared, negative_source)
    np.testing.assert_array_equal(partner.column_weights, prepared.column_weights)
    np.testing.assert_array_equal(partner.column_weights, negative_source.column_weights)
    np.testing.assert_array_equal(partner.source_covariance, negative_source.covariance)
    assert abs(partner.source_covariance[0, 1]) > 0
    np.testing.assert_array_equal(
        partner.weighted_columns, s3_conjugate(prepared.weighted_columns))
    assert partner.coverage['quadrature_applied_once'] is True
    assert partner.coverage['angular_multiplicity_added'] is False


def test_invalid_energy_and_weights_are_rejected():
    positive_source, negative_source = signed_sources()
    prepared = evolve(positive_source, A_up(3), ANGULAR, amplitude=0., tangents='zero')
    wrong_energy = FixedSourcePreparation(
        negative_source.covariance, negative_source.column_weights, prepared.source_energies)
    with pytest.raises(ValueError, match='energies'):
        negative_angular_partner(prepared, wrong_energy)
    wrong_weights = FixedSourcePreparation(
        negative_source.covariance, 2. * negative_source.column_weights,
        negative_source.energies)
    with pytest.raises(ValueError, match='weights'):
        negative_angular_partner(prepared, wrong_weights)
    with pytest.raises(TypeError, match='FixedSourcePreparation'):
        negative_angular_partner(prepared, negative_source.covariance)
    with pytest.raises(TypeError, match='KSDifferenceIncoming'):
        negative_angular_partner(object(), negative_source)
    bad = np.array(negative_source.covariance, copy=True)
    bad[0, 0] = bad[0, 0] + 1e-6
    with pytest.raises(ValueError, match='source law'):
        negative_angular_partner(
            prepared, FixedSourcePreparation(
                bad, negative_source.column_weights, negative_source.energies))


def test_source_complement_is_not_the_evolved_field_isometry():
    positive_source, negative_source = signed_sources()
    nsrc = positive_source.covariance.shape[0]
    C_neg = np.array(negative_source.covariance, copy=True)
    C_neg[0, 0] = C_neg[0, 0] + 1e-15
    supplied = FixedSourcePreparation(
        C_neg, negative_source.column_weights, negative_source.energies)
    assert source_complement_residual(
        positive_source.covariance, supplied.covariance) < SOURCE_COMPLEMENT_TOLERANCE
    prepared = evolve(positive_source, A_up(nsrc), ANGULAR)
    partner = negative_angular_partner(prepared, supplied)
    np.testing.assert_array_equal(partner.column_tangents,s3_conjugate(prepared.column_tangents))
    np.testing.assert_array_equal(partner.axial_tangents,s3_conjugate(prepared.axial_tangents))
    np.testing.assert_array_equal(partner.source_covariance, supplied.covariance)
    complement = np.eye(nsrc) - prepared.source_covariance.conj()
    assert not np.array_equal(partner.source_covariance, complement)
    assert partner.diagnostics['covariance_replaced_by_rounded_complement'] is False
    assert partner.diagnostics['field_isometry_assumed'] is False
    assert partner.diagnostics['source_complement_residual'] == pytest.approx(
        source_complement_residual(prepared.source_covariance, supplied.covariance))
    point = prepared.weighted_columns[2]
    gram = point @ point.conj().T
    assert np.max(np.abs(gram - np.eye(2))) > 1e-3
    density = point @ prepared.source_covariance @ point.conj().T
    mapped_point = partner.weighted_columns[2]
    mapped_density = mapped_point @ partner.source_covariance @ mapped_point.conj().T
    naive = np.eye(2) - density
    naive_conj = np.eye(2) - density.conj()
    assert np.max(np.abs(mapped_density - naive)) > 1e-3
    assert np.max(np.abs(mapped_density - naive_conj)) > 1e-3
    nz = len(prepared.z)
    identity = np.eye(nz * 2, dtype=complex).reshape(nz, 2, nz, 2)
    field_naive = identity - prepared.covariance()
    assert np.max(np.abs(partner.covariance() - field_naive)) > 1e-3
    F = partner.weighted_columns.reshape(nz * 2, nsrc)
    expected = (F @ partner.source_covariance @ F.conj().T).reshape(nz, 2, nz, 2)
    np.testing.assert_allclose(partner.covariance(), expected, atol=2e-15, rtol=0)
