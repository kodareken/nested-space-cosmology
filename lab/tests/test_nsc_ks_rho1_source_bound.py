"""Coverage, preservation and validator controls for the rho=1 source bound."""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from recursive_horizons.nsc_ks_upstream_error import packed_upper_fraction
from recursive_horizons.nsc_ks_rho1_source_bound import (
    ELL0_GROUPS,
    HISTORICAL_COARSE_AGGREGATE,
    REMAINING_RECONSTRUCTION_PANELS,
    REMAINING_RECONSTRUCTION_ROWS,
    V5_ACTION_COVERED_PANELS,
    V5_ACTION_COVERED_ROWS,
    attach_reconstruction_bound,
    audit_row_threshold_partition,
    assert_coherence_retained,
    assert_ell0_baseline_retained,
    build_rho1_coverage,
    classify_rho1_windows,
    combine_directed_n_beta,
    coverage_identity,
    inspect_authenticated_sample,
    load_authenticated_inventory_payload,
    pair_negative_energy,
    physical_source_error_for_upstream,
    reject_double_counting,
    reject_historical_coarse,
    restore_directed_upper,
    restore_rational,
    row_identities,
    serialize_directed_upper,
    serialize_rational,
    signed_channel_multiplicity,
    validate_rho1_source_bound_record,
)
from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_retarded_radial_response import ks_generator
from recursive_horizons.nsc_transmitting_dirac_domain import S3
import derive_nsc_ks_rho1_source_bound_v1 as D


CONFIG = {'surface_gravity': .24, 'omega': 3.9}


def synthetic_panel(*, group=14, sign=1, angular=np.sqrt(5), mass=1.5, name='synthetic'):
    energies = np.array([.4, .8, 2.])
    rng = np.random.default_rng(33)
    columns = rng.normal(size=(3, 2, 3)) + 1j * rng.normal(size=(3, 2, 3))
    covariance = np.array([
        source_covariance(energy, CONFIG['surface_gravity'], CONFIG['omega'], mass)
        for energy in energies])
    return ReferenceSourcePanel(
        name, group, sign, mass, angular, energies, np.array([.1, .2, .3]),
        covariance, columns, {'fixture': True, 'method': 'accepted physical low-mode reconstruction'})


def inventory_and_v5():
    inventory = json.loads((ROOT / 'results/development/nsc-ks-source-inventory.json').read_text())
    v5 = json.loads((ROOT / 'results/development/nsc-incoming-source-update-v5.json').read_text())
    decision = json.loads((ROOT / 'results/development/nsc-ks-low-subgap-decision-v2.json').read_text())
    channels = json.loads(
        (ROOT / 'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    return inventory, v5, decision, channels


def contribution(identity, *, n=None, beta=None, method='synthetic-directed'):
    return {
        'identity': identity,
        'N': n,
        'beta': beta,
        'method': method,
        'coherence_retained': True,
        'weights_retained': True,
        'minus_energy_retained': True,
        'degeneracy_retained': True,
    }


def test_covariance_coherence_is_retained_and_independent_max_rejected():
    panel = synthetic_panel()
    report = assert_coherence_retained(panel.covariance, require_offdiag=True)
    assert report['coherence_retained']
    assert report['offdiag_max'] > 0
    assert abs(panel.covariance[0, 0, 1]) > 0
    pieces = list(panel.batches(2))
    np.testing.assert_array_equal(
        np.concatenate([piece.covariance for piece in pieces]), panel.covariance)
    assert abs(pieces[0].covariance[0, 0, 1]) > 0
    with pytest.raises(ValueError, match='forbidden'):
        attach_reconstruction_bound({'panel': panel.name}, Fraction(1, 10 ** 12),
                                    'independent_column_max')
    with pytest.raises(ValueError, match='forbidden'):
        attach_reconstruction_bound({'panel': panel.name}, Fraction(1, 10 ** 12),
                                    'drop_coherence')
    diagonal = np.zeros_like(panel.covariance)
    for index in range(3):
        diagonal[:, index, index] = panel.covariance[:, index, index]
    with pytest.raises(ValueError, match='off-diagonal'):
        assert_coherence_retained(diagonal, require_offdiag=True)


def test_signed_multiplicity_counts_each_sign_once():
    _inventory, _v5, _decision, channels = inventory_and_v5()
    per_sign, signs = signed_channel_multiplicity(channels[14])
    assert signs == 2
    assert per_sign == channels[14]['copy_count'] * channels[14]['degeneracy'] / 2
    assert per_sign * signs == channels[14]['copy_count'] * channels[14]['degeneracy']
    ell0_per_sign, ell0_signs = signed_channel_multiplicity(channels[13])
    assert ell0_signs == 1
    assert ell0_per_sign == channels[13]['copy_count'] * channels[13]['degeneracy']
    folded = channels[14]['copy_count'] * channels[14]['degeneracy'] * 2
    assert folded != per_sign * signs


def test_ell0_zero_response_retains_baseline():
    assert ELL0_GROUPS == (0, 13, 23)
    assert assert_ell0_baseline_retained(13, baseline_retained=True) == 'exact_zero'
    with pytest.raises(ValueError, match='baseline'):
        assert_ell0_baseline_retained(13, baseline_retained=False)
    with pytest.raises(ValueError, match='no synthetic'):
        assert_ell0_baseline_retained(0, baseline_retained=True, panel='lll/fake')
    with pytest.raises(ValueError, match='nontrivial'):
        assert_ell0_baseline_retained(14, baseline_retained=True, radius_response='exact_zero')
    panel = synthetic_panel(group=13, sign=1, angular=0., mass=np.pi / 2, name='group13/low')
    partner = pair_negative_energy(panel, CONFIG)
    assert partner.angular_sign == 1
    assert panel.angular_magnitude == 0


def test_energy_pairing_is_antiunitary_not_folding():
    positive = synthetic_panel()
    negative = pair_negative_energy(positive, CONFIG)
    assert negative.angular_sign == -1
    assert negative.provenance['energy_folding_factor'] == 1
    np.testing.assert_array_equal(negative.weights, positive.weights)
    np.testing.assert_array_equal(negative.energies, -positive.energies)
    np.testing.assert_allclose(
        negative.covariance, np.eye(3) - positive.covariance.conj(), atol=2e-16, rtol=0)
    for rho in (1., 1.015, 1.03):
        plus = ks_generator(rho, positive.energies, positive.mass, positive.angular)
        minus = ks_generator(rho, negative.energies, negative.mass, negative.angular)
        np.testing.assert_array_equal(minus, S3 @ plus.conj() @ S3)
    with pytest.raises(ValueError, match='forbidden'):
        attach_reconstruction_bound({}, Fraction(1, 10 ** 12), 'drop_minus_energy')


def test_malformed_coverage_is_rejected():
    inventory, v5, decision, channels = inventory_and_v5()
    windows = deepcopy(decision['coverage']['windows'])
    build_rho1_coverage(windows, inventory_record=inventory, v5_coverage=v5['coverage'],
                        channels=channels)
    duplicate = deepcopy(windows)
    extra = deepcopy(next(row for row in duplicate if row.get('panel') == 'low/1_1'))
    extra['panel'] = 'low/1_-1'
    extra['sign'] = -1
    duplicate.append(extra)
    with pytest.raises(ValueError, match='duplicate'):
        build_rho1_coverage(duplicate, inventory_record=inventory, v5_coverage=v5['coverage'],
                            channels=channels)
    joined = deepcopy(windows)
    target = next(row for row in joined if row.get('panel') == 'low/1_1')
    target['v5_joined'] = [16, 'infinity']
    with pytest.raises(ValueError, match='joined'):
        build_rho1_coverage(joined, inventory_record=inventory, v5_coverage=v5['coverage'],
                            channels=channels)
    ell0 = deepcopy(windows)
    target = next(row for row in ell0 if row.get('panel') == 'group13/low')
    target['sign'] = -1
    with pytest.raises(ValueError, match='plus sign'):
        build_rho1_coverage(ell0, inventory_record=inventory, v5_coverage=v5['coverage'],
                            channels=channels)
    lll = deepcopy(windows)
    target = next(row for row in lll if row.get('group') == 0)
    target['panel'] = 'lll/fake'
    with pytest.raises(ValueError, match='synthetic mode panel'):
        build_rho1_coverage(lll, inventory_record=inventory, v5_coverage=v5['coverage'],
                            channels=channels)


def test_double_count_rejection():
    reject_double_counting(['low/1_1', 'subgap8/14_1', 'group14/low32_1'])
    with pytest.raises(ValueError, match='duplicate'):
        reject_double_counting(['low/1_1', 'low/1_1'])
    with pytest.raises(ValueError, match='replace'):
        reject_double_counting(['low_ref/32_1', 'low/32_1'])
    keys = row_identities(14, 1, 'group14/low32_1', 2)
    assert keys[0] == coverage_identity(14, 1, 'group14/low32_1', 0)
    with pytest.raises(ValueError, match='double-counted'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20)),
             contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20))],
            declared_universe_size=2)


def test_lossless_directed_serialization():
    value = Fraction(3, 7)
    assert restore_rational(serialize_rational(value)) == value
    packed = serialize_directed_upper(value)
    assert packed_upper_fraction(packed) >= value
    assert restore_directed_upper(packed).is_finite()
    assert serialize_directed_upper(None) is None
    assert restore_directed_upper(None) is None
    missing = attach_reconstruction_bound({'panel': 'synthetic'}, None, 'unevaluated')
    assert missing['reconstruction_error_upper'] is None


def test_historical_coarse_aggregate_is_rejected():
    with pytest.raises(ValueError, match='4.998e-9'):
        reject_historical_coarse(HISTORICAL_COARSE_AGGREGATE)
    with pytest.raises(ValueError, match='forbidden'):
        attach_reconstruction_bound({}, HISTORICAL_COARSE_AGGREGATE,
                                    'historical_coarse_4.998e-9')


def test_missing_bounds_keep_directed_n_beta_open():
    result = combine_directed_n_beta(
        [contribution((14, 1, 'synthetic', 0)), contribution((14, 1, 'synthetic', 1))],
        declared_universe_size=2)
    assert result['N'] is None and result['beta'] is None
    assert result['status'].startswith('OPEN')
    closed = combine_directed_n_beta(
        [contribution((14, 1, 'synthetic', 0), n=Fraction(1, 10), beta=Fraction(1, 20)),
         contribution((14, -1, 'synthetic', 0), n=Fraction(1, 10), beta=Fraction(1, 20))],
        declared_universe_size=2)
    assert closed['N'] is None and closed['beta'] is None
    assert closed['status'].startswith('OPEN')
    assert closed['finite_inputs_do_not_close_v1']
    with pytest.raises(ValueError, match='declared universe'):
        combine_directed_n_beta(
            [contribution((14, 1, 'synthetic', 0), n=Fraction(1, 10), beta=Fraction(1, 20))],
            declared_universe_size=REMAINING_RECONSTRUCTION_PANELS)
    feed = physical_source_error_for_upstream(None)
    assert feed['nine_component_feed']['combined_upstream_bound_N_beta'] is None
    assert feed['nine_component_feed']['combined_status'].startswith('OPEN')
    with pytest.raises(ValueError, match='no authenticated physical'):
        physical_source_error_for_upstream((8.916253820835393e-12,
                                            3.882373567206641e-18))
    with pytest.raises(ValueError, match='no authenticated physical'):
        physical_source_error_for_upstream((4.997e-9, 4.997e-9))


def test_authenticated_coverage_counts_match_v5_and_inventory():
    inventory, v5, decision, channels = inventory_and_v5()
    coverage = build_rho1_coverage(
        decision['coverage']['windows'], inventory_record=inventory,
        v5_coverage=v5['coverage'], channels=channels)
    assert coverage['remaining_reconstruction_panels'] == REMAINING_RECONSTRUCTION_PANELS
    assert coverage['remaining_reconstruction_rows'] == REMAINING_RECONSTRUCTION_ROWS
    assert coverage['v5_action_covered_panels'] == V5_ACTION_COVERED_PANELS
    assert coverage['v5_action_covered_rows'] == V5_ACTION_COVERED_ROWS
    assert coverage['inventory_positive_panels'] == inventory['positive_panels']
    assert coverage['inventory_positive_rows'] == inventory['positive_energy_rows']
    assert coverage['remaining_by_kind'] == {
        'low_mode_reconstruction': 64, 'subgap_rows': 38}
    assert coverage['remaining_physical_rho1_source_error_bound'] is None
    assert coverage['directed_N_beta_contribution_bound']['N'] is None
    assert not coverage['double_counted']
    D.ensure_inventory_payload(inventory)
    arrays, _meta, _digest = load_authenticated_inventory_payload(ROOT, inventory)
    partition = audit_row_threshold_partition(arrays, coverage)
    assert partition['rows_below_joined_in_remaining_panels'] == 27828
    assert partition['rows_at_or_above_joined_in_remaining_panels'] == 1528
    assert partition['straddling_remaining_panels'] == 59
    assert partition['physical_column_accuracy_universe_panels'] == 164
    assert partition['physical_column_accuracy_universe_rows'] == 49372
    assert partition['remaining_panel_campaign_is_not_complete_physical_source_universe']
    assert sha256((ROOT / 'results/development/nsc-incoming-source-update-v5.json').read_bytes()
                  ).hexdigest() == json.loads(
        (ROOT / 'results/development/nsc-ks-low-subgap-decision-v2.json').read_text()
    )['input_hashes']['results/development/nsc-incoming-source-update-v5.json']


def test_classifier_keeps_reconstruction_unbounded():
    panel = synthetic_panel()
    classified = classify_rho1_windows(
        {'radius_response_exactly_zero_groups': [0, 13, 23]},
        [panel],
        {'14': {'joined': [16, 'infinity'], 'actual_signs': [1, -1]}})
    recon = next(row for row in classified['windows'] if row.get('panel') == panel.name)
    assert recon['source_reconstruction'] == 'unbounded'
    assert recon['error_bound'] is None
    assert recon['baseline_action_region'] == 'not_covered'
    assert classified['remaining_low_subgap_source_error_bound'] is None


def test_small_authenticated_real_sample_preserves_fibers():
    inventory, v5, decision, channels = inventory_and_v5()
    D.ensure_inventory_payload(inventory)
    arrays, meta, digest = load_authenticated_inventory_payload(ROOT, inventory)
    assert digest == inventory['payload']['sha256']
    coverage = build_rho1_coverage(
        decision['coverage']['windows'], inventory_record=inventory,
        v5_coverage=v5['coverage'], channels=channels)
    sample = inspect_authenticated_sample(
        arrays, meta, coverage['remaining_windows'], meta['config'],
        panel_names=('group14/low32_1', 'group13/low'), rows_per_panel=2)
    assert sample['inspected_sample_panels'] == 2
    assert sample['inspected_sample_rows'] == 3
    assert sample['reconstruction_bounds_evaluated_panels'] == 0
    assert not sample['full_remaining_campaign_run']
    by_name = {row['panel']: row for row in sample['panels']}
    assert by_name['group14/low32_1']['offdiag_max'] > 0
    assert by_name['group14/low32_1']['coherence_retained']
    assert by_name['group14/low32_1']['energy_folding_factor'] == 1
    assert by_name['group13/low']['radius_history_response'] == 'exact_zero'
    assert by_name['group13/low']['baseline_retained']
    assert by_name['group14/low32_1']['reconstruction_error_upper'] is None
    remaining_names = {item['panel'] for item in coverage['remaining_windows']}
    assert 'group14/low32_1' in remaining_names
    assert 'group13/low' in remaining_names
    with np.load(ROOT / inventory['payload']['path'], allow_pickle=False) as payload:
        assert len(payload['group14/low32_1/energies']) == 96
        assert len(payload['group13/low/energies']) == 192


def test_recorded_register_stays_open_and_authenticates_inputs():
    inventory, _v5, _decision, _channels = inventory_and_v5()
    D.ensure_inventory_payload(inventory)
    record = D.calculate()
    D.check(record)
    validate_rho1_source_bound_record(record)
    assert record['coverage']['remaining_reconstruction_panels'] == 102
    assert record['coverage']['remaining_reconstruction_rows'] == 29356
    assert record['sample']['inspected_sample_rows'] == 3
    assert record['remaining_physical_rho1_source_error_bound'] is None
    assert record['directed_N_beta_contribution_bound']['N'] is None
    assert record['directed_N_beta_contribution_bound']['beta'] is None
    assert record['physical_EXISTENCE_certificate'] is False
    assert record['physical_NONEXISTENCE_certificate'] is False
    assert not record['historical_coarse_4p998e_minus9_is_current_remaining_bound']
    if D.OUTPUT.is_file():
        saved = json.loads(D.OUTPUT.read_text())
        assert D.stable(record) == D.stable(saved)
        D.check(saved)
        for path, expected in {**saved['source_hashes'], **saved['input_hashes']}.items():
            assert sha256((ROOT / path).read_bytes()).hexdigest() == expected

    forged = deepcopy(record)
    forged['status'] = 'PASS'
    with pytest.raises(ValueError, match='must remain OPEN'):
        validate_rho1_source_bound_record(forged)
    forged = deepcopy(record)
    forged['physical_NONEXISTENCE_certificate'] = True
    with pytest.raises(ValueError, match='cannot issue'):
        validate_rho1_source_bound_record(forged)
    forged = deepcopy(record)
    forged['coverage']['directed_N_beta_contribution_bound']['N'] = 0.0
    with pytest.raises(ValueError, match='no finite directed'):
        validate_rho1_source_bound_record(forged)
    forged = deepcopy(record)
    forged['upstream_physical_slot']['combined_upstream_bound_N_beta'] = [0.0, 0.0]
    with pytest.raises(ValueError, match='must remain OPEN'):
        validate_rho1_source_bound_record(forged)


def test_runtime_memory_probe_is_cross_platform():
    value, units, source = D.peak_rss()
    assert value is None or value > 0
    assert units in {'bytes', 'kilobytes', 'unavailable'}
    assert isinstance(source, str) and source
