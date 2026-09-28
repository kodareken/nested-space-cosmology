"""Synthetic KS local-assembly controls; no field generators or archived science."""
from dataclasses import dataclass

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, source_column_matter, surface_geometry_response,
)
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_local_constraints import (
    _REFERENCE_CACHE, assemble_local_incoming, homogeneous_reference_amplitudes,
    stationary_computational_reference,
)
from recursive_horizons.nsc_local_incoming_constraints import (
    assess_local_residual, local_error_budget,
)
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeIncoming, computational_z_grid, evolve_ks_source_envelope,
    physical_incoming_interval, usual_axial_support,
)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily
from recursive_horizons.nsc_pg_ks_metric_pullback import ks_to_pg, reference_chart
from recursive_horizons.nsc_retarded_radial_response import ks_generator
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain


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


def A_up(nsrc, seed=2):
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(2, nsrc)) + 1j * rng.normal(size=(2, nsrc))) * 0.2


def bump_family(amplitude=0.002):
    center = chart_coordinates(1.)[1] + 0.15
    w = LocalAxialFunction((0.05, -0.02, 0.01), center)
    U = LocalAxialFunction((0.02, -0.005), center)
    return CompatibleIncomingMetric(
        (float(amplitude),), (CompatibleRadiusDirection(w, U, 0.007, 0.03),))


def zero_family():
    return LocalIncomingFamily(np.zeros((2, 8)))


def evolve(family, source, initial, z=None, target=None, tangents='all',
           mass=MASS, angular=ANGULAR):
    z = computational_z_grid(16) if z is None else z
    interval = physical_incoming_interval()
    if target is None:
        target = np.linspace(interval[0], interval[1], 7)
    return evolve_ks_source_envelope(
        source, initial, family, z, target, mass, angular, RHO_UP,
        axial_support=usual_axial_support(), tangents=tangents, **SOLVER)


def coefficients_for():
    _, _, axial, _, _ = reference_chart(1.)
    r = np.sqrt(2.)
    N, beta, q, r = ks_to_pg(1., 1., 0., axial, r)
    domain = TransmittingDiracSeamDomain(
        lapse=float(N), radial_scale=float(q), shift=float(beta), radius=float(r))
    return {
        'a': domain.induced_axial_scale,
        'r': float(r),
        'Hr': -0.4,
        'A0': 0.05,
        'A1': 0.02,
        'd': 0.3,
        'e': -0.05,
        'C': 0.11,
        'D': np.array([0., 0.01, 0.02, 0.0, 0.0]),
        'F': np.array([16.5, 0.1, 0.0]),
    }


def channel(index, mass, angular, sign, copies=2, degeneracy=12):
    return {
        'index': index,
        'compact_mass': mass,
        'angular_eigenvalue': abs(angular),
        'angular_sign': sign,
        'copy_count': copies,
        'degeneracy': degeneracy,
    }


def z_interval(prepared):
    z = np.asarray(prepared.z, float)
    return (float(z[0]), float(z[-1]))


def test_homogeneous_reference_matches_radial_ode_and_is_cached_across_histories():
    source = tiny_source()
    initial = A_up(4)
    zero = evolve(zero_family(), source, initial, tangents='zero')
    bump = evolve(bump_family(), source, initial, tangents='zero')
    assert zero.fixed_preparation_digest == bump.fixed_preparation_digest
    _REFERENCE_CACHE.clear()
    Aref = homogeneous_reference_amplitudes(zero)
    assert len(_REFERENCE_CACHE) == 1

    def rhs(rho, state):
        field = state.reshape(2, 4)
        G = ks_generator(rho, source.energies, MASS, ANGULAR)
        return np.einsum('sab,bs->as', G, field).ravel()

    radial = solve_ivp(
        rhs, (RHO_UP, 1.), initial.ravel(), method='DOP853',
        rtol=SOLVER['rtol'], atol=SOLVER['atol'], max_step=SOLVER['max_step'])
    assert radial.success
    np.testing.assert_allclose(Aref, radial.y[:, -1].reshape(2, 4), atol=3e-11, rtol=0)
    again = homogeneous_reference_amplitudes(bump)
    assert again is Aref
    assert len(_REFERENCE_CACHE) == 1
    reference, reference_z = stationary_computational_reference(zero)
    phase = np.exp(-1j * zero.source_energies * zero.z[:, None])
    np.testing.assert_allclose(reference, Aref[None] * phase[:, None, :], atol=0, rtol=0)
    np.testing.assert_allclose(reference_z, -1j * zero.source_energies * reference, atol=0, rtol=0)
    coeff = coefficients_for()
    ch = {'14_1': channel(14, MASS, ANGULAR, 1)}
    budget = local_error_budget()
    interval = z_interval(zero)
    zero_res = assemble_local_incoming(
        {'14_1': zero}, zero_family(), np.zeros(2), coeff, ch, interval, budget)
    bump_res = assemble_local_incoming(
        {'14_1': bump}, bump_family(), np.zeros(2), coeff, ch, interval, budget)
    np.testing.assert_allclose(
        zero_res['family_reference_matter']['14_1'],
        bump_res['family_reference_matter']['14_1'], atol=3e-14, rtol=0)
    assert len(_REFERENCE_CACHE) == 1
    assert zero_res['scope']['reference_tangent'] == 0
    assert zero_res['physical_constraint_status'] == 'OPEN'


def test_pde_derivative_assembly_matches_raw_weighted_pair():
    prepared = evolve(bump_family(), tiny_source(), A_up(4))
    coeff = coefficients_for()
    ch = channel(14, MASS, ANGULAR, 1, copies=2, degeneracy=12)
    result = assemble_local_incoming(
        {'14_1': prepared}, bump_family(), np.array([0.31, -0.17]), coeff,
        {'14_1': ch}, z_interval(prepared), local_error_budget())
    F, Fz = prepared.weighted_columns, prepared.weighted_axial_columns
    dF, dFz = prepared.weighted_column_tangents, prepared.weighted_axial_tangents
    current = source_column_matter(
        F, Fz, prepared.source_covariance, dF, dFz,
        mass=MASS, angular=ANGULAR, axial_scale=coeff['a'], radius=coeff['r'],
        multiplicity=2 * 12 / 2)
    reference, reference_z = stationary_computational_reference(prepared)
    zeros = np.zeros_like(dF)
    ref = source_column_matter(
        reference * prepared.column_weights, reference_z * prepared.column_weights,
        prepared.source_covariance, zeros, zeros,
        mass=MASS, angular=ANGULAR, axial_scale=coeff['a'], radius=coeff['r'],
        multiplicity=12.)
    np.testing.assert_allclose(
        result['family_current_matter']['14_1'], current['action_gradient'], atol=3e-14, rtol=0)
    np.testing.assert_allclose(
        result['family_corrections']['14_1'],
        current['action_gradient'] - ref['action_gradient'], atol=3e-14, rtol=0)
    frozen = source_column_matter(
        F, -1j * prepared.source_energies * F, prepared.source_covariance,
        dF, -1j * prepared.source_energies * dF,
        mass=MASS, angular=ANGULAR, axial_scale=coeff['a'], radius=coeff['r'],
        multiplicity=12.)
    assert np.max(np.abs(current['action_gradient'] - frozen['action_gradient'])) > 1e-8
    assert result['scope']['pde_axial_derivatives'] is True
    assert result['family_records']['14_1']['multiplicity'] == 12.


def test_history_jacobian_matches_state_finite_difference():
    source = tiny_source()
    initial = A_up(4)
    z = computational_z_grid(32)
    interval = physical_incoming_interval()
    target = np.linspace(interval[0], interval[1], 7)
    h = 1e-4
    center = evolve(bump_family(0.002), source, initial, z=z, target=target, tangents='all')
    plus = evolve(bump_family(0.002 + h), source, initial, z=z, target=target, tangents='zero')
    minus = evolve(bump_family(0.002 - h), source, initial, z=z, target=target, tangents='zero')
    coeff = coefficients_for()
    ch = {'14_1': channel(14, MASS, ANGULAR, 1)}
    baseline = np.array([0.31, -0.17])
    bounds = z_interval(center)
    budget = local_error_budget()
    assembled = assemble_local_incoming(
        {'14_1': center}, bump_family(0.002), baseline, coeff, ch, bounds, budget)
    plus_res = assemble_local_incoming(
        {'14_1': plus}, bump_family(0.002 + h), baseline, coeff, ch, bounds, budget)
    minus_res = assemble_local_incoming(
        {'14_1': minus}, bump_family(0.002 - h), baseline, coeff, ch, bounds, budget)
    fd = (plus_res['action_gradient'] - minus_res['action_gradient']) / (2 * h)
    np.testing.assert_allclose(fd, assembled['history_jacobian'][0], atol=3e-8, rtol=0)
    assert assembled['scope']['reference_state_tangent_subtracted'] is False
    assert assembled['history_jacobian'].shape == (1, len(center.z), 2)


def test_two_family_baseline_and_geometry_counted_once():
    source = tiny_source()
    initial = A_up(4)
    family = bump_family()
    a = evolve(family, source, initial, mass=MASS, angular=ANGULAR)
    b = evolve(family, source, initial, mass=0.8, angular=-1.)
    coeff = coefficients_for()
    baseline = np.array([0.31, -0.17])
    interval = z_interval(a)
    budget = local_error_budget()
    channels = {
        '14_1': channel(14, MASS, ANGULAR, 1, copies=2, degeneracy=12),
        '7_m': channel(7, 0.8, 1., -1, copies=1, degeneracy=4),
    }
    both = assemble_local_incoming(
        {'14_1': a, '7_m': b}, family, baseline, coeff, channels, interval, budget)
    only_a = assemble_local_incoming(
        {'14_1': a}, family, baseline, coeff, {'14_1': channels['14_1']}, interval, budget)
    only_b = assemble_local_incoming(
        {'7_m': b}, family, baseline, coeff, {'7_m': channels['7_m']}, interval, budget)
    geometry = both['local_reference_gradient_change']
    np.testing.assert_allclose(
        both['action_gradient'] - only_a['action_gradient'] - only_b['action_gradient']
        + baseline + geometry, 0, atol=3e-14, rtol=0)
    np.testing.assert_allclose(
        both['family_corrections']['14_1'] + both['family_corrections']['7_m'] + baseline + geometry,
        both['action_gradient'], atol=3e-14, rtol=0)
    assert both['scope']['baseline_and_geometry_counted_once'] is True
    assert both['family_records']['14_1']['multiplicity'] == 12.
    assert both['family_records']['7_m']['multiplicity'] == 2.
    slots, tangent = compatible_history_slots(family, a.z, 1)
    geometric = surface_geometry_response(slots, tangent, coeff)
    np.testing.assert_allclose(
        both['local_reference_gradient_change'], geometric['action_gradient_change'],
        atol=3e-14, rtol=0)


def test_coherence_weights_and_signed_multiplicity():
    prepared = evolve(bump_family(), tiny_source(), A_up(4))
    coeff = coefficients_for()
    result = assemble_local_incoming(
        {'14_1': prepared}, bump_family(), np.zeros(2), coeff,
        {'14_1': channel(14, MASS, ANGULAR, 1, copies=2, degeneracy=12)},
        z_interval(prepared), local_error_budget())
    F, Fz = prepared.weighted_columns, prepared.weighted_axial_columns
    dF, dFz = prepared.weighted_column_tangents, prepared.weighted_axial_tangents
    C = prepared.source_covariance
    params = dict(mass=MASS, angular=ANGULAR, axial_scale=coeff['a'], radius=coeff['r'],
                  multiplicity=12.)
    coherent = source_column_matter(F, Fz, C, dF, dFz, **params)
    diagonal = source_column_matter(F, Fz, np.diag(np.diag(C)), dF, dFz, **params)
    assert np.max(np.abs(coherent['action_gradient'] - diagonal['action_gradient'])) > 1e-8
    weights = prepared.column_weights
    doubled = source_column_matter(F * weights, Fz * weights, C, dF * weights, dFz * weights, **params)
    assert np.max(np.abs(coherent['action_gradient'] - doubled['action_gradient'])) > 1e-4
    flipped = source_column_matter(F, Fz, C, dF, dFz, mass=MASS, angular=-ANGULAR,
                                   axial_scale=coeff['a'], radius=coeff['r'], multiplicity=12.)
    assert np.max(np.abs(coherent['action_gradient'] - flipped['action_gradient'])) > 1e-8
    assert result['family_records']['14_1']['n_actual_angular_signs'] == 2
    assert result['scope']['finite_energy_samples_are_complete_spectrum'] is False
    assert result['physical_constraint_status'] == 'OPEN'


def test_mismatched_prep_meshes_and_legacy_rejection():
    family = bump_family()
    source = tiny_source()
    initial = A_up(4)
    prepared = evolve(family, source, initial)
    other = evolve(family, tiny_source(weights=(1.0, 0.9)), initial, z=computational_z_grid(32),
                   mass=MASS, angular=-ANGULAR)
    coeff = coefficients_for()
    interval = z_interval(prepared)
    budget = local_error_budget()
    with pytest.raises(ValueError, match='preparation digest'):
        assemble_local_incoming(
            {'14_1': prepared}, family, np.zeros(2), coeff,
            {'14_1': channel(14, MASS + 0.2, ANGULAR, 1)}, interval, budget)
    with pytest.raises(ValueError, match='preparation digest'):
        assemble_local_incoming(
            {'14_1': prepared}, family, np.zeros(2), coeff,
            {'14_1': channel(14, MASS, ANGULAR, -1)}, interval, budget)
    bad = dict(coeff)
    bad['a'] = coeff['a'] + 0.01
    with pytest.raises(ValueError, match='intrinsic'):
        assemble_local_incoming(
            {'14_1': prepared}, family, np.zeros(2), bad,
            {'14_1': channel(14, MASS, ANGULAR, 1)}, interval, budget)
    with pytest.raises(ValueError, match='generating amplitudes'):
        assemble_local_incoming(
            {'14_1': prepared}, bump_family(0.003), np.zeros(2), coeff,
            {'14_1': channel(14, MASS, ANGULAR, 1)}, interval, budget)
    with pytest.raises(TypeError, match='KSEnvelopeIncoming'):
        assemble_local_incoming(
            {'14_1': prepared.source_covariance}, family, np.zeros(2), coeff,
            {'14_1': channel(14, MASS, ANGULAR, 1)}, interval, budget)
    with pytest.raises(ValueError, match='frozen incoming C0'):
        assemble_local_incoming(
            {'14_1': {'incoming_C0': np.eye(2), 'stress_approximant': np.zeros(4)}},
            family, np.zeros(2), coeff, {'14_1': channel(14, MASS, ANGULAR, 1)}, interval, budget)
    with pytest.raises(ValueError, match='node_fields'):
        assemble_local_incoming(
            {'14_1': {'node_fields': np.ones(3), 'history': True}},
            family, np.zeros(2), coeff, {'14_1': channel(14, MASS, ANGULAR, 1)}, interval, budget)
    independent = assemble_local_incoming(
        {'14_1': prepared, '14_other': other}, family, np.zeros(2), coeff,
        {'14_1': channel(14, MASS, ANGULAR, 1),
         '14_other': channel(14, MASS, ANGULAR, -1)}, interval, budget)
    assert independent['family_records']['14_1']['source_digest'] != independent['family_records']['14_other']['source_digest']
    assert independent['family_records']['14_1']['computational_count'] == 16
    assert independent['family_records']['14_other']['computational_count'] == 32
    assert independent['scope']['computational_meshes_need_not_match'] is True
    wrong = {**channel(14, MASS, ANGULAR, -1), 'expected_source_digest': source.digest}
    with pytest.raises(ValueError, match='declared family source'):
        assemble_local_incoming(
            {'14_other': other}, family, np.zeros(2), coeff, {'14_other': wrong}, interval, budget)
    with pytest.raises(ValueError, match='duplicate family'):
        assemble_local_incoming(
            [('14_1', prepared), ('14_1', other)], family, np.zeros(2), coeff,
            {'14_1': channel(14, MASS, ANGULAR, 1)}, interval, budget)
    with pytest.raises(ValueError, match='mismatched incoming z'):
        moved = evolve(family, tiny_source(weights=(1.0, 0.9)), initial,
                       target=np.linspace(physical_incoming_interval()[0],
                                          physical_incoming_interval()[1], 5),
                       mass=MASS, angular=-ANGULAR)
        assemble_local_incoming(
            {'14_1': prepared, '14_other': moved}, family, np.zeros(2), coeff,
            {'14_1': channel(14, MASS, ANGULAR, 1),
             '14_other': channel(14, MASS, ANGULAR, -1)}, interval, budget)


def test_missing_budget_stays_open():
    prepared = evolve(bump_family(), tiny_source(), A_up(4), tangents='zero')
    coeff = coefficients_for()
    ch = {'14_1': channel(14, MASS, ANGULAR, 1)}
    interval = z_interval(prepared)
    with pytest.raises(ValueError, match='error_budget'):
        assemble_local_incoming(
            {'14_1': prepared}, bump_family(), np.zeros(2), coeff, ch, interval, None)
    with pytest.raises(ValueError, match='missing mandatory'):
        assemble_local_incoming(
            {'14_1': prepared}, bump_family(), np.zeros(2), coeff, ch, interval,
            {'baseline_covered_regions': None})
    result = assemble_local_incoming(
        {'14_1': prepared}, bump_family(), np.zeros(2), coeff, ch, interval, local_error_budget())
    assert result['physical_constraint_status'] == 'OPEN'
    assert result['physical_EXISTENCE_certificate'] is False
    assert result['error_budget']['baseline_low_subgap'] is None
    assessed = assess_local_residual(result)
    assert assessed['conditional_numeric_tolerance'] is False
    assert assessed['physical_EXISTENCE_certificate'] is False
    assert assessed['physical_constraint_status'] == 'OPEN'
    assert assessed['missing_bounds_are_not_zero'] is True
    zero = {
        'z': result['z'],
        'action_gradient': np.zeros_like(result['action_gradient']),
        'interval': interval,
        'error_budget': local_error_budget(),
    }
    assert assess_local_residual(zero)['physical_constraint_status'] == 'OPEN'
    closed = local_error_budget(
        baseline_covered_regions=(0., 0.),
        baseline_low_subgap=(0., 0.),
        changed_history_finite_energy=(0., 0.),
        changed_history_tail=(0., 0.),
        preparation_field_axial_derivative=(0., 0.),
        coefficient_arithmetic=(0., 0.),
        between_node_remainder=(0., 0.),
    )
    numeric = assess_local_residual(zero, error_budget=closed)
    assert numeric['conditional_numeric_tolerance'] is True
    assert numeric['physical_EXISTENCE_certificate'] is False
    assert numeric['physical_constraint_status'] == 'OPEN'


def test_subclass_isinstance_and_ell0_documentation():
    prepared = evolve(bump_family(), tiny_source(), A_up(4), mass=1.1, angular=0.)

    @dataclass(frozen=True)
    class DifferenceEnvelopeIncoming(KSEnvelopeIncoming):
        pass

    child = DifferenceEnvelopeIncoming(
        **{field: getattr(prepared, field)
           for field in KSEnvelopeIncoming.__dataclass_fields__})
    assert isinstance(child, KSEnvelopeIncoming)
    result = assemble_local_incoming(
        {'0_1': child}, bump_family(), np.array([0.2, 0.1]), coefficients_for(),
        {'0_1': channel(0, 1.1, 0., 1, copies=3, degeneracy=4)},
        z_interval(child), local_error_budget())
    record = result['family_records']['0_1']
    assert record['ell0_pure_radius_identity'] is True
    assert record['n_actual_angular_signs'] == 1
    assert record['multiplicity'] == 12.
    assert 'identically zero' in record['ell0_documentation']
    assert result['physical_constraint_status'] == 'OPEN'
