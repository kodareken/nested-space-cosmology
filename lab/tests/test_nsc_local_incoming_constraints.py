"""Synthetic local-assembly controls on the real exact-phase owner; no field generators."""
import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, source_column_matter, surface_geometry_response,
)
from recursive_horizons.nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric, FixedSourcePreparation,
)
from recursive_horizons.nsc_exact_phase_prepared_state import (
    prepare_exact_phase_incoming,
)
from recursive_horizons.nsc_cf4_prepared_state import prepare_cf4_incoming
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_local_incoming_constraints import (
    assemble_local_incoming, assess_local_residual, local_error_budget,
    stationary_computational_reference,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
from recursive_horizons.nsc_transmitting_resolvent import profile


def signed_source(weights=(1.2, 0.7)):
    return FixedSourcePreparation.from_signed_blocks(
        np.array([0.4, -0.55]),
        np.array(weights, dtype=float),
        np.array([
            [[0.45, 0.07 + 0.04j], [0.07 - 0.04j, 0.6]],
            [[0.50, 0.02 - 0.03j], [0.02 + 0.03j, 0.35]],
        ]),
    )


def delayed_bump(center, halfwidth):
    def value(z, order):
        if order not in (0, 1, 2, 3):
            raise ValueError('owned incoming derivative order 0..3 required')
        u = (z - center) / halfwidth
        if abs(u) >= 1:
            return 0.
        b, bp = profile(u)
        if order == 0:
            return b
        if order == 1:
            return bp / halfwidth
        d = 1 - u * u
        g1 = -2 * u / d ** 2
        g2 = -2 / d ** 2 - 8 * u * u / d ** 3
        if order == 2:
            return b * (g1 * g1 + g2) / halfwidth ** 2
        g3 = -24 * u / d ** 3 - 48 * u ** 3 / d ** 4
        return b * (g1 ** 3 + 3 * g1 * g2 + g3) / halfwidth ** 3
    return value


def delayed_family(amplitude=0.003, inner=0.05, outer=0.25, shift=0.025, halfwidth=0.018):
    x = np.linspace(-2., 1.2, 17)
    i8 = int(np.argmin(np.abs(x - 0.8)))
    center = chart_coordinates(float(x[i8]))[1] + shift
    directions = (CompatibleRadiusDirection(
        delayed_bump(center, halfwidth), lambda z, n: 0., inner, outer),)
    return CompatibleIncomingMetric((float(amplitude),), directions)


def initial_columns(owner, nsrc):
    x = owner.x
    n = len(x)
    rng = np.linspace(0., 1., n)
    columns = np.zeros((2 * n, nsrc), dtype=complex)
    for k in range(nsrc):
        columns[:n, k] = (0.2 + 0.03 * k) * np.exp(-rng * rng) + 0.04j * np.sin((k + 1) * rng)
        columns[n:, k] = (0.15 - 0.02 * k) * np.exp(-0.3 * rng * rng) + 0.05j * np.cos((k + 2) * rng)
        columns[n - 1, k] += 0.12 + 0.03j * (k + 1)
        columns[2 * n - 1, k] += 0.09 - 0.02j * (k + 1)
    return columns


def setup(mass=1.3, angular=2., amplitude=0.003, ndir=1, source=None):
    x = np.linspace(-2., 1.2, 17)
    owner = FourthOrderModePropagator(x, mass, angular)
    source = signed_source() if source is None else source
    phi = initial_columns(owner, source.covariance.shape[0])
    times = np.linspace(0., 0.04, 5)
    family = delayed_family(amplitude)
    provider = family if ndir else AmplitudeOnlyMetric(family)
    return owner, provider, family, times, phi, source


def test_cf4_prepared_state_integrates_without_dropping_its_gauss_binding():
    owner, provider, family, times, phi, source = setup()
    prepared = prepare_cf4_incoming(owner, times, provider, phi, source)
    result = assemble_local_incoming(
        {'14_1': prepared}, provider, np.zeros(2), coefficients_for(prepared),
        {'14_1': {**channel(14, 1.3, 2., 1), 'expected_source_digest': source.digest}},
        (prepared.state.z[0], prepared.state.z[-1]), local_error_budget())
    assert np.isfinite(result['action_gradient']).all()
    assert result['history_jacobian'].shape == (1, len(times), 2)
    assert result['physical_constraint_status'] == 'OPEN'
    prepared.require_history(provider, times, owner.x)


def coefficients_for(prepared):
    N, beta, q, r = prepared.state.history.rho1_metric
    domain = TransmittingDiracSeamDomain(lapse=N, radial_scale=q, shift=beta, radius=r)
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
        'angular_eigenvalue': angular,
        'angular_sign': sign,
        'copy_count': copies,
        'degeneracy': degeneracy,
    }


def z_interval(prepared):
    z = np.asarray(prepared.state.z, float)
    return (float(z[0]), float(z[-1]))


def prepare_family(mass=1.3, angular=2., amplitude=0.003, ndir=1, source=None):
    owner, provider, family, times, phi, source = setup(
        mass, angular, amplitude, ndir, source)
    prepared = prepare_exact_phase_incoming(owner, times, provider, phi, source)
    return prepared, provider, family


def test_pde_derivative_assembly_matches_raw_weighted_pair():
    prepared, provider, _ = prepare_family()
    coeff = coefficients_for(prepared)
    ch = channel(14, 1.3, 2., 1, copies=2, degeneracy=12)
    result = assemble_local_incoming(
        {'14_1': prepared}, provider, np.array([0.31, -0.17]), coeff,
        {'14_1': ch}, z_interval(prepared), local_error_budget())
    F, Fz = prepared.weighted_columns, prepared.weighted_axial_columns
    dF, dFz = prepared.weighted_column_tangents, prepared.weighted_axial_tangents
    current = source_column_matter(
        F, Fz, prepared.state.source_covariance, dF, dFz,
        mass=1.3, angular=2., axial_scale=coeff['a'], radius=coeff['r'],
        multiplicity=2 * 12 / 2)
    reference, reference_z = stationary_computational_reference(prepared)
    zeros = np.zeros_like(dF)
    ref = source_column_matter(
        reference * prepared.state.column_weights,
        reference_z * prepared.state.column_weights,
        prepared.state.source_covariance, zeros, zeros,
        mass=1.3, angular=2., axial_scale=coeff['a'], radius=coeff['r'],
        multiplicity=12.)
    np.testing.assert_allclose(
        result['family_current_matter']['14_1'], current['action_gradient'], atol=3e-14, rtol=0)
    np.testing.assert_allclose(
        result['family_corrections']['14_1'],
        current['action_gradient'] - ref['action_gradient'], atol=3e-14, rtol=0)
    energy = prepared.state.source_energies
    frozen = source_column_matter(
        F, -1j * energy * F, prepared.state.source_covariance,
        dF, -1j * energy * dF,
        mass=1.3, angular=2., axial_scale=coeff['a'], radius=coeff['r'],
        multiplicity=12.)
    assert np.max(np.abs(current['action_gradient'] - frozen['action_gradient'])) > 1e-6
    assert result['scope']['pde_axial_derivatives'] is True
    assert result['scope']['interpolated_axial_derivatives'] is False
    assert result['family_records']['14_1']['multiplicity'] == 12.


def test_history_jacobian_matches_state_finite_difference():
    center, provider, _ = prepare_family(amplitude=0.003, ndir=1)
    plus, plus_provider, _ = prepare_family(amplitude=0.003 + 1e-4, ndir=0)
    minus, minus_provider, _ = prepare_family(amplitude=0.003 - 1e-4, ndir=0)
    coeff = coefficients_for(center)
    ch = {'14_1': channel(14, 1.3, 2., 1)}
    baseline = np.array([0.31, -0.17])
    interval = z_interval(center)
    budget = local_error_budget()
    assembled = assemble_local_incoming(
        {'14_1': center}, provider, baseline, coeff, ch, interval, budget)
    plus_res = assemble_local_incoming(
        {'14_1': plus}, plus_provider, baseline, coeff, ch, interval, budget)
    minus_res = assemble_local_incoming(
        {'14_1': minus}, minus_provider, baseline, coeff, ch, interval, budget)
    fd = (plus_res['action_gradient'] - minus_res['action_gradient']) / (2e-4)
    np.testing.assert_allclose(fd, assembled['history_jacobian'][0], atol=3e-8, rtol=0)
    assert assembled['scope']['reference_tangent'] == 0
    assert assembled['scope']['reference_state_tangent_subtracted'] is False


def test_two_family_baseline_and_geometry_counted_once():
    a, provider, _ = prepare_family(mass=1.3, angular=2.)
    b, provider_b, _ = prepare_family(mass=0.8, angular=-1.)
    np.testing.assert_allclose(provider.amplitudes, provider_b.amplitudes, atol=0, rtol=0)
    coeff = coefficients_for(a)
    baseline = np.array([0.31, -0.17])
    interval = z_interval(a)
    budget = local_error_budget()
    channels = {
        '14_1': channel(14, 1.3, 2., 1, copies=2, degeneracy=12),
        '7_m': channel(7, 0.8, 1., -1, copies=1, degeneracy=4),
    }
    both = assemble_local_incoming(
        {'14_1': a, '7_m': b}, provider, baseline, coeff, channels, interval, budget)
    only_a = assemble_local_incoming(
        {'14_1': a}, provider, baseline, coeff, {'14_1': channels['14_1']}, interval, budget)
    only_b = assemble_local_incoming(
        {'7_m': b}, provider, baseline, coeff, {'7_m': channels['7_m']}, interval, budget)
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
    slots, tangent = compatible_history_slots(provider, a.state.z, 1)
    geometric = surface_geometry_response(slots, tangent, coeff)
    np.testing.assert_allclose(
        both['local_reference_gradient_change'], geometric['action_gradient_change'],
        atol=3e-14, rtol=0)
    np.testing.assert_allclose(only_a['local_reference_gradient_change'],
                               only_b['local_reference_gradient_change'], atol=3e-14, rtol=0)


def test_coherence_weights_and_signed_multiplicity():
    prepared, provider, _ = prepare_family()
    coeff = coefficients_for(prepared)
    ch = channel(14, 1.3, 2., 1, copies=2, degeneracy=12)
    result = assemble_local_incoming(
        {'14_1': prepared}, provider, np.zeros(2), coeff,
        {'14_1': ch}, z_interval(prepared), local_error_budget())
    F, Fz = prepared.weighted_columns, prepared.weighted_axial_columns
    dF, dFz = prepared.weighted_column_tangents, prepared.weighted_axial_tangents
    C = prepared.state.source_covariance
    params = dict(mass=1.3, angular=2., axial_scale=coeff['a'], radius=coeff['r'], multiplicity=12.)
    coherent = source_column_matter(F, Fz, C, dF, dFz, **params)
    diagonal = source_column_matter(F, Fz, np.diag(np.diag(C)), dF, dFz, **params)
    assert np.max(np.abs(coherent['action_gradient'] - diagonal['action_gradient'])) > 1e-8
    weights = prepared.state.column_weights
    doubled = source_column_matter(F * weights, Fz * weights, C, dF * weights, dFz * weights, **params)
    assert np.max(np.abs(coherent['action_gradient'] - doubled['action_gradient'])) > 1e-4
    flipped = source_column_matter(F, Fz, C, dF, dFz, mass=1.3, angular=-2.,
                                   axial_scale=coeff['a'], radius=coeff['r'], multiplicity=12.)
    assert np.max(np.abs(coherent['action_gradient'] - flipped['action_gradient'])) > 1e-8
    assert result['family_records']['14_1']['n_actual_angular_signs'] == 2
    assert result['scope']['finite_energy_samples_are_complete_spectrum'] is False
    assert result['physical_constraint_status'] == 'OPEN'


def test_ell0_pure_radius_identity_is_documented_not_inferred():
    prepared, provider, _ = prepare_family(mass=1.1, angular=0.)
    coeff = coefficients_for(prepared)
    result = assemble_local_incoming(
        {'0_1': prepared}, provider, np.array([0.2, 0.1]), coeff,
        {'0_1': channel(0, 1.1, 0., 1, copies=3, degeneracy=4)},
        z_interval(prepared), local_error_budget())
    record = result['family_records']['0_1']
    assert record['ell0_pure_radius_identity'] is True
    assert record['n_actual_angular_signs'] == 1
    assert record['multiplicity'] == 12.
    assert 'identically zero' in record['ell0_documentation']
    assert 'baseline retained' in record['ell0_documentation']
    assert result['scope']['ell0_pure_radius_identities']['0_1'] == record['ell0_documentation']
    # The identity is documented from ell=0, never from a small numerical tangent.
    assert result['family_matter_tangents']['0_1'].shape[0] == 1


def test_intrinsic_and_preparation_binding():
    prepared, provider, family = prepare_family()
    coeff = coefficients_for(prepared)
    interval = z_interval(prepared)
    budget = local_error_budget()
    with pytest.raises(ValueError, match='preparation digest'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), coeff,
            {'14_1': channel(14, 1.7, 2., 1)}, interval, budget)
    with pytest.raises(ValueError, match='preparation digest'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), coeff,
            {'14_1': channel(14, 1.3, 2., -1)}, interval, budget)
    bad = dict(coeff)
    bad['a'] = coeff['a'] + 0.01
    with pytest.raises(ValueError, match='intrinsic'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), bad,
            {'14_1': channel(14, 1.3, 2., 1)}, interval, budget)
    other = delayed_family(0.003, shift=0.030)
    with pytest.raises(ValueError, match='actual sampled'):
        assemble_local_incoming(
            {'14_1': prepared}, other, np.zeros(2), coeff,
            {'14_1': channel(14, 1.3, 2., 1)}, interval, budget)
    local = LocalIncomingFamily(np.zeros((2, 8)))
    with pytest.raises(ValueError, match='bound to its generating metric family'):
        assemble_local_incoming(
            {'14_1': prepared}, local, np.zeros(2), coeff,
            {'14_1': channel(14, 1.3, 2., 1)}, interval, budget)
    with pytest.raises(TypeError, match='ExactPhaseIncoming'):
        assemble_local_incoming(
            {'14_1': prepared.state}, provider, np.zeros(2), coeff,
            {'14_1': channel(14, 1.3, 2., 1)}, interval, budget)
    with pytest.raises(ValueError, match='frozen incoming C0'):
        assemble_local_incoming(
            {'14_1': {'incoming_C0': np.eye(2), 'stress_approximant': np.zeros(4)}},
            provider, np.zeros(2), coeff,
            {'14_1': channel(14, 1.3, 2., 1)}, interval, budget)


def test_duplicate_and_coverage_rejection():
    a, provider, _ = prepare_family()
    b, _, _ = prepare_family(mass=0.8, angular=-1.)
    coeff = coefficients_for(a)
    interval = z_interval(a)
    budget = local_error_budget()
    channels = {
        '14_1': channel(14, 1.3, 2., 1),
        '7_m': channel(7, 0.8, 1., -1, copies=1, degeneracy=4),
    }
    with pytest.raises(ValueError, match='duplicate family'):
        assemble_local_incoming(
            [('14_1', a), ('14_1', b)], provider, np.zeros(2), coeff,
            channels, interval, budget)
    with pytest.raises(ValueError, match='duplicate family'):
        assemble_local_incoming(
            {'14_1': a, '14_again': a}, provider, np.zeros(2), coeff,
            {'14_1': channel(14, 1.3, 2., 1), '14_again': channel(14, 1.3, 2., 1)},
            interval, budget)
    with pytest.raises(ValueError, match='sampled family coverage'):
        assemble_local_incoming(
            {'14_1': a, '7_m': b}, provider, np.zeros(2), coeff, channels, interval, budget,
            coverage={'sampled_family_ids': ['14_1']})
    with pytest.raises(ValueError, match='complete spectrum'):
        assemble_local_incoming(
            {'14_1': a}, provider, np.zeros(2), coeff, {'14_1': channels['14_1']},
            interval, budget, coverage={'complete_spectrum': True})
    other, _, _ = prepare_family(mass=1.3, angular=-2., source=signed_source(weights=(1.1, 0.8)))
    # Separate signs may use distinct energy quadratures; fixed preparation
    # is checked against each actual family's declaration, not another sign.
    independent = assemble_local_incoming(
        {'14_1': a, '14_other': other}, provider, np.zeros(2), coeff,
        {'14_1': channel(14, 1.3, 2., 1),
         '14_other': channel(14, 1.3, 2., -1)}, interval, budget)
    assert independent['family_records']['14_1']['source_digest'] != independent['family_records']['14_other']['source_digest']
    wrong = {**channel(14, 1.3, 2., -1), 'expected_source_digest': a.state.history.source_digest}
    with pytest.raises(ValueError, match='declared family source'):
        assemble_local_incoming({'14_other': other}, provider, np.zeros(2), coeff,
                                {'14_other': wrong}, interval, budget)
    times = np.linspace(0., 0.05, 6)
    owner, _, family, _, phi, source = setup()
    moved = prepare_exact_phase_incoming(owner, times, family, phi, source)
    with pytest.raises(ValueError, match='mismatched times'):
        assemble_local_incoming(
            {'14_1': a, '7_m': moved}, provider, np.zeros(2), coeff,
            {'14_1': channels['14_1'],
             '7_m': channel(7, 1.3, 2., 1)},
            interval, budget)


def test_missing_budget_interval_and_endpoints_prevent_false_pass():
    prepared, provider, _ = prepare_family()
    coeff = coefficients_for(prepared)
    ch = {'14_1': channel(14, 1.3, 2., 1)}
    interval = z_interval(prepared)
    with pytest.raises(ValueError, match='error_budget'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), coeff, ch, interval, None)
    with pytest.raises(ValueError, match='missing mandatory'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), coeff, ch, interval,
            {'baseline_covered_regions': None})
    with pytest.raises(ValueError, match='interval'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), coeff, ch, None, local_error_budget())
    with pytest.raises(ValueError, match='positive-length'):
        assemble_local_incoming(
            {'14_1': prepared}, provider, np.zeros(2), coeff, ch,
            (interval[0], interval[0]), local_error_budget())
    result = assemble_local_incoming(
        {'14_1': prepared}, provider, np.zeros(2), coeff, ch, interval, local_error_budget())
    assert result['physical_constraint_status'] == 'OPEN'
    assert result['physical_EXISTENCE_certificate'] is False
    zero = {
        'z': result['z'],
        'action_gradient': np.zeros_like(result['action_gradient']),
        'interval': interval,
        'error_budget': local_error_budget(),
    }
    assessed = assess_local_residual(zero)
    assert assessed['conditional_numeric_tolerance'] is False
    assert assessed['physical_EXISTENCE_certificate'] is False
    assert assessed['physical_constraint_status'] == 'OPEN'
    assert assessed['missing_bounds_are_not_zero'] is True
    with pytest.raises(ValueError, match='endpoints'):
        assess_local_residual(result, interval=(interval[0] - 1., interval[0] - 0.5))
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
    S1 = chart_coordinates(1.)[1]
    with pytest.raises(ValueError, match='endpoints'):
        assess_local_residual(result, interval=(S1 + 0.12, S1 + 0.18))


def test_partial_one_family_control_remains_open():
    prepared, provider, _ = prepare_family()
    result = assemble_local_incoming(
        {'14_1': prepared}, provider, np.array([0.31, -0.17]),
        coefficients_for(prepared), {'14_1': channel(14, 1.3, 2., 1)},
        z_interval(prepared), local_error_budget(),
        coverage={'sampled_family_ids': ['14_1'], 'missing_family_ids': 'OPEN'})
    assert result['scope']['sampled_family_ids'] == ('14_1',)
    assert result['scope']['missing_family_ids'] == 'OPEN'
    assert result['scope']['energy_coverage'] == 'OPEN'
    assert result['scope']['all_retained_ids_cannot_certify_energy_coverage'] is True
    assert result['physical_constraint_status'] == 'OPEN'
    assessed = assess_local_residual(result)
    assert assessed['conditional_numeric_tolerance'] is False
    assert assessed['physical_EXISTENCE_certificate'] is False
