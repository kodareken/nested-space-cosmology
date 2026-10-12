"""Synthetic batched KS constraint controls; no inventory campaign or generators."""
import numpy as np
import pytest
from dataclasses import replace

from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, source_column_matter, surface_geometry_response,
)
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_batched_constraints import (
    KSConstraintAccumulator, bind_upstream_batch, map_negative_batch,
    assemble_batched_incoming, signed_incoming_pair,
)
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_local_constraints import assemble_local_incoming
from recursive_horizons.nsc_ks_source_envelope import (
    computational_z_grid, physical_incoming_interval, usual_axial_support,
)
from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from recursive_horizons.nsc_local_incoming_constraints import (
    assess_local_residual, local_error_budget,
)
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from test_nsc_ks_local_constraints import bump_family, channel, coefficients_for


CONFIG = {'surface_gravity': 0.24, 'omega': 3.9}
MASS = np.pi / 2
ANGULAR = np.sqrt(5)
RHO_UP = 1.03
SOLVER = dict(rtol=2e-10, atol=2e-12, max_step=0.002)
BASELINE = np.array([0.31, -0.17])


def make_panel(name='synthetic', group=14, sign=1, mass=MASS, angular=ANGULAR,
               energies=(0.4, 0.8, 2.0), weights=(0.1, 0.2, 0.3), seed=33):
    energies = np.asarray(energies, float)
    rng = np.random.default_rng(seed)
    columns = rng.normal(size=(len(energies), 2, 3)) + 1j * rng.normal(size=(len(energies), 2, 3))
    blocks = np.array([
        source_covariance(float(energy), CONFIG['surface_gravity'], CONFIG['omega'], mass)
        for energy in energies])
    return ReferenceSourcePanel(
        name, group, sign, mass, abs(angular), energies, np.asarray(weights, float),
        blocks, columns, {'fixture': True})


def bind(panel):
    return bind_upstream_batch(panel, rho_up=RHO_UP, **SOLVER)


def test_original_rows_cannot_mislabel_or_split_a_coherent_fiber():
    batch=bind(make_panel())
    with pytest.raises(ValueError,match='three complete'):
        replace(batch,rows=(0,2))
    with pytest.raises(ValueError,match='integer original'):
        replace(batch,rows=(.5,3.5))
    with pytest.raises(ValueError,match='angular label'):
        replace(batch,angular_sign=-1)
    covariance=np.array(batch.source.covariance,copy=True)
    covariance[0,3]=covariance[3,0]=1e-8
    source=FixedSourcePreparation(covariance,batch.source.column_weights,batch.source.energies)
    with pytest.raises(ValueError,match='cross-energy coherences'):
        replace(batch,source=source)


def evolve_batch(batch, family, z=None, target=None, tangents='all'):
    z = computational_z_grid(16) if z is None else z
    if target is None:
        target = np.linspace(*physical_incoming_interval(), 5)
    return evolve_ks_difference_envelope(
        batch.source, batch.initial_columns, family, z, target, batch.mass, batch.angular,
        batch.rho_up, axial_support=usual_axial_support(), tangents=tangents, **SOLVER)


def negative_source(panel):
    partner = panel.negative_partner(CONFIG)
    return FixedSourcePreparation.from_signed_blocks(
        partner.energies, partner.weights, partner.covariance)


def interval_of(prepared):
    z = np.asarray(prepared.z, float)
    return (float(z[0]), float(z[-1]))


def assemble(entries, family, baseline=BASELINE, coverage=None):
    first = entries[0][1] if len(entries[0]) == 4 else entries[0][0]
    return assemble_batched_incoming(
        entries, family, baseline, coefficients_for(), interval_of(first),
        local_error_budget(), coverage=coverage)


def test_split_vs_unsplit_coherent_source_contraction():
    panel = make_panel()
    family = bump_family()
    whole_batch = bind(panel)
    whole = evolve_batch(whole_batch, family)
    coeff = coefficients_for()
    ch = channel(14, MASS, ANGULAR, 1)
    params = dict(mass=MASS, angular=ANGULAR, axial_scale=coeff['a'], radius=coeff['r'],
                  multiplicity=12.)
    F, Fz = whole.weighted_columns, whole.weighted_axial_columns
    dF, dFz = whole.weighted_column_tangents, whole.weighted_axial_tangents
    C = whole.source_covariance
    full = source_column_matter(F, Fz, C, dF, dFz, **params)
    head = source_column_matter(
        F[..., :6], Fz[..., :6], C[:6, :6], dF[..., :6], dFz[..., :6], **params)
    tail = source_column_matter(
        F[..., 6:], Fz[..., 6:], C[6:, 6:], dF[..., 6:], dFz[..., 6:], **params)
    np.testing.assert_allclose(
        head['action_gradient'] + tail['action_gradient'], full['action_gradient'],
        atol=3e-14, rtol=0)
    assert abs(C[0, 1]) > 0
    pieces = list(panel.batches(2))
    split_batches = [bind(piece) for piece in pieces]
    split_states = [evolve_batch(batch, family) for batch in split_batches]
    unsplit = assemble([(whole, whole_batch, ch)], family)
    split = assemble(list(zip(split_states, split_batches, [ch, ch])), family)
    np.testing.assert_allclose(split['action_gradient'], unsplit['action_gradient'], atol=3e-11, rtol=0)
    np.testing.assert_allclose(split['history_jacobian'], unsplit['history_jacobian'], atol=3e-11, rtol=0)
    existing = assemble_local_incoming(
        {'14_1': whole}, family, BASELINE, coeff, {'14_1': ch}, interval_of(whole),
        local_error_budget())
    np.testing.assert_allclose(unsplit['action_gradient'], existing['action_gradient'], atol=3e-14, rtol=0)
    assert unsplit['scope']['dense_all_energy_covariance_allocated'] is False
    assert split['physical_constraint_status'] == 'OPEN'
    assert split['scope']['sampled_family_ids']==unsplit['scope']['sampled_family_ids']
    assert len(split['scope']['sampled_batch_ids'])==2
    assert len(split['scope']['sampled_family_ids'])==1


def test_nonzero_history_finite_difference_tangent():
    panel = make_panel(energies=(0.4, 0.8), weights=(0.1, 0.2))
    pieces = list(panel.batches(1))
    batches = [bind(piece) for piece in pieces]
    h = 1e-4
    ch = channel(14, MASS, ANGULAR, 1)
    z = computational_z_grid(16)
    target = np.linspace(*physical_incoming_interval(), 5)

    def run(amplitude, tangents):
        family = bump_family(amplitude)
        states = [evolve_batch(batch, family, z=z, target=target, tangents=tangents)
                  for batch in batches]
        return assemble(list(zip(states, batches, [ch, ch])), family)

    center, plus, minus = run(0.002, 'all'), run(0.002 + h, 'zero'), run(0.002 - h, 'zero')
    fd = (plus['action_gradient'] - minus['action_gradient']) / (2 * h)
    np.testing.assert_allclose(fd, center['history_jacobian'][0], atol=3e-8, rtol=0)
    assert center['scope']['reference_state_tangent_subtracted'] is False
    assert center['history_jacobian'].shape == (1, 5, 2)


def test_baseline_and_geometry_counted_once():
    family = bump_family()
    a_panel = make_panel(name='a', group=14, seed=33)
    b_panel = make_panel(name='b', group=7, sign=-1, mass=0.8, angular=1.,
                         energies=(0.5,), weights=(0.2,), seed=8)
    a_batch, b_batch = bind(a_panel), bind(b_panel)
    a_state, b_state = evolve_batch(a_batch, family), evolve_batch(b_batch, family)
    ch_a, ch_b = channel(14, MASS, ANGULAR, 1, copies=2, degeneracy=12), channel(
        7, 0.8, 1., -1, copies=1, degeneracy=4)
    both = assemble([(a_state, a_batch, ch_a), (b_state, b_batch, ch_b)], family)
    only_a = assemble([(a_state, a_batch, ch_a)], family)
    only_b = assemble([(b_state, b_batch, ch_b)], family)
    geometry = both['local_reference_gradient_change']
    np.testing.assert_allclose(
        both['action_gradient'] - only_a['action_gradient'] - only_b['action_gradient']
        + BASELINE + geometry, 0, atol=3e-14, rtol=0)
    np.testing.assert_allclose(
        sum(both['batch_corrections'].values()) + BASELINE + geometry, both['action_gradient'],
        atol=3e-14, rtol=0)
    assert both['scope']['baseline_and_geometry_counted_once'] is True
    assert both['scope']['baseline_and_geometry_added_per_batch'] is False
    assert both['family_records']['14_1_E+1']['multiplicity'] == 12.
    assert both['family_records']['7_-1_E+1']['multiplicity'] == 2.
    slots, tangent = compatible_history_slots(family, a_state.z, 1)
    geometric = surface_geometry_response(slots, tangent, coefficients_for())
    np.testing.assert_allclose(
        both['local_reference_gradient_change'], geometric['action_gradient_change'],
        atol=3e-14, rtol=0)


def test_overlapping_rows_and_changed_prep_profile_rejection():
    panel = make_panel()
    family = bump_family()
    pieces = list(panel.batches(2))
    batch_a, batch_b = bind(pieces[0]), bind(pieces[1])
    state_a = evolve_batch(batch_a, family)
    state_b = evolve_batch(batch_b, family, z=computational_z_grid(32))
    ch = channel(14, MASS, ANGULAR, 1)
    acc = KSConstraintAccumulator(
        family, BASELINE, coefficients_for(), interval_of(state_a), local_error_budget())
    acc.add(state_a, batch_a, ch, batch_id='a')
    with pytest.raises(ValueError, match='overlapping original panel rows'):
        acc.add(state_a, batch_a, ch, batch_id='again')
    overlap_panel = ReferenceSourcePanel(
        panel.name + '/overlap', panel.group, panel.angular_sign, panel.mass,
        panel.angular_magnitude, panel.energies[1:], panel.weights[1:],
        panel.covariance[1:], panel.amplitudes_at_one[1:],
        {'parent_panel': panel.name, 'rows': [1, 3], 'source': panel.provenance})
    overlap_batch = bind(overlap_panel)
    overlap_state = evolve_batch(overlap_batch, family, tangents='zero')
    with pytest.raises(ValueError, match='overlapping original panel rows'):
        acc.add(overlap_state, overlap_batch, ch, batch_id='overlap')
    with pytest.raises(ValueError, match='preparation digest'):
        acc.add(state_b, batch_b, channel(14, MASS + 0.2, ANGULAR, 1), batch_id='mass')
    other = bind(make_panel(name='other', seed=99))
    other_state = evolve_batch(other, family, tangents='zero')
    with pytest.raises(ValueError, match='prepared source differs'):
        acc.add(other_state, batch_b, ch, batch_id='source')
    diff_A = ReferenceSourcePanel(
        'diffA', panel.group, panel.angular_sign, panel.mass, panel.angular_magnitude,
        pieces[1].energies, pieces[1].weights, pieces[1].covariance,
        pieces[1].amplitudes_at_one * 1.7,
        {'parent_panel': 'diffA', 'rows': [0, 1], 'source': {'fixture': True}})
    diff_state = evolve_batch(bind(diff_A), family, tangents='zero')
    with pytest.raises(ValueError, match='A_up/preparation digest'):
        acc.add(diff_state, batch_b, ch, batch_id='prep')
    with pytest.raises(ValueError, match='generating amplitudes'):
        KSConstraintAccumulator(
            bump_family(0.003), BASELINE, coefficients_for(), interval_of(state_a),
            local_error_budget()).add(state_a, batch_a, ch)
    with pytest.raises(TypeError, match='KSDifferenceIncoming'):
        acc.add(object(), batch_b, ch)
    mixed = assemble(
        [('fine', state_a, batch_a, ch), ('coarse', state_b, batch_b, ch)], family)
    assert mixed['batch_records']['fine']['computational_count'] == 16
    assert mixed['batch_records']['coarse']['computational_count'] == 32
    assert mixed['scope']['computational_meshes_need_not_match'] is True
    with pytest.raises(ValueError, match='LLL'):
        bind(make_panel(name='lll', group=0, angular=0.))


def test_sign_inventory_and_weights_preserved():
    panel = make_panel(energies=(0.7,), weights=(1.2,))
    family = bump_family()
    batch = bind(panel)
    prepared = evolve_batch(batch, family)
    source_neg = negative_source(panel)
    prepared, batch, partner, neg_batch = signed_incoming_pair(prepared, batch, source_neg)
    np.testing.assert_array_equal(partner.column_weights, prepared.column_weights)
    np.testing.assert_array_equal(partner.column_weights, source_neg.column_weights)
    np.testing.assert_array_equal(neg_batch.source.column_weights, batch.source.column_weights)
    np.testing.assert_array_equal(partner.source_covariance, source_neg.covariance)
    assert abs(partner.source_covariance[0, 1]) > 0
    assert partner.diagnostics['energy_folding_factor'] == 1
    assert partner.diagnostics['field_isometry_assumed'] is False
    assert neg_batch.energy_sign == -1
    assert neg_batch.angular_sign == -1
    ch_p, ch_m = channel(14, MASS, ANGULAR, 1), channel(14, MASS, ANGULAR, -1)
    both = assemble([(prepared, batch, ch_p), (partner, neg_batch, ch_m)], family)
    positive = assemble([(prepared, batch, ch_p)], family)
    pos_corr = next(iter(positive['batch_corrections'].values()))
    both_corr = sum(both['batch_corrections'].values())
    assert np.max(np.abs(both_corr - 2 * pos_corr)) > 1e-8
    assert np.max(np.abs(both['action_gradient'] - 2 * positive['action_gradient'])) > 1e-8
    assert both['scope']['positive_result_doubled'] is False
    assert both['scope']['energy_folding_factor'] == 1
    assert both['batch_records'][next(i for i, rec in both['batch_records'].items()
                                      if rec['energy_sign'] < 0)]['energy_folding_factor'] == 1
    gram = prepared.weighted_columns[2] @ prepared.weighted_columns[2].conj().T
    assert np.max(np.abs(gram - np.eye(2))) > 1e-3
    mapped = map_negative_batch(batch, source_neg)
    np.testing.assert_array_equal(mapped.initial_columns, partner.initial_columns)
    assert mapped.continuation['history_independent'] is True
    assert batch.continuation['history_independent'] is True


def test_missing_error_budget_stays_open():
    panel = make_panel(name='ell0', group=13, angular=0., energies=(0.5,), weights=(0.2,))
    family = bump_family()
    batch = bind(panel)
    prepared = evolve_batch(batch, family, tangents='zero')
    ch = channel(13, MASS, 0., 1, copies=3, degeneracy=4)
    coeff = coefficients_for()
    interval = interval_of(prepared)
    with pytest.raises(ValueError, match='error_budget'):
        KSConstraintAccumulator(family, np.zeros(2), coeff, interval, None)
    with pytest.raises(ValueError, match='missing mandatory'):
        assemble_batched_incoming(
            [(prepared, batch, ch)], family, np.zeros(2), coeff, interval,
            {'baseline_covered_regions': None})
    result = assemble([(prepared, batch, ch)], family, baseline=np.zeros(2))
    record = next(iter(result['batch_records'].values()))
    assert record['ell0_pure_radius_identity'] is True
    assert record['n_actual_angular_signs'] == 1
    assert record['multiplicity'] == 12.
    assert 'identically zero' in record['ell0_documentation']
    assert result['physical_constraint_status'] == 'OPEN'
    assert result['physical_EXISTENCE_certificate'] is False
    assert result['error_budget']['baseline_low_subgap'] is None
    assert result['error_budget']['changed_history_finite_energy'] is None
    assert result['scope']['missing_family_ids'] == 'OPEN'
    assert result['scope']['missing_energy_ranges'] == 'OPEN'
    assert result['scope']['energy_coverage'] == 'OPEN'
    assessed = assess_local_residual(result)
    assert assessed['conditional_numeric_tolerance'] is False
    assert assessed['physical_EXISTENCE_certificate'] is False
    assert assessed['physical_constraint_status'] == 'OPEN'
    assert assessed['missing_bounds_are_not_zero'] is True
    zero = {
        'z': result['z'], 'action_gradient': np.zeros_like(result['action_gradient']),
        'interval': interval, 'error_budget': local_error_budget(),
    }
    assert assess_local_residual(zero)['physical_constraint_status'] == 'OPEN'


def test_exact_zero_radius_groups_keep_baseline_and_reject_double_counting():
    family=bump_family();batch=bind(make_panel(energies=(.7,),weights=(.2,)))
    state=evolve_batch(batch,family,tangents='zero')
    acc=KSConstraintAccumulator(family,BASELINE,coefficients_for(),interval_of(state),local_error_budget())
    acc.add(state,batch,channel(14,MASS,ANGULAR,1))
    for group,mass in ((0,0.),(13,np.pi/2),(23,np.pi)):
        acc.declare_radius_invariant_channel(channel(group,mass,0.,1,copies=1,degeneracy=4))
    with pytest.raises(ValueError,match='already included'):
        acc.declare_radius_invariant_channel(channel(13,np.pi/2,0.,1))
    with pytest.raises(ValueError,match='ell=0 groups'):
        acc.declare_radius_invariant_channel(channel(14,MASS,ANGULAR,1))
    result=acc.finalize()
    np.testing.assert_array_equal(result['baseline_action_gradient'],np.broadcast_to(BASELINE,(len(state.z),2)))
    assert set(result['scope']['analytic_zero_response_groups'])=={0,13,23}
    assert result['scope']['analytic_zero_response_groups'][13]['all_source_energies']
    assert result['physical_constraint_status']=='OPEN'
