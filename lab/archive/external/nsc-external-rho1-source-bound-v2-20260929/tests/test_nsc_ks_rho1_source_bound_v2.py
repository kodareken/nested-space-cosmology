"""Split, construction-trace, Lipschitz and validator controls for rho=1 v2."""
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

from recursive_horizons.nsc_ks_rho1_source_bound import (
    HISTORICAL_COARSE_AGGREGATE,
    INVENTORY_POSITIVE_PANELS,
    INVENTORY_POSITIVE_ROWS,
    REMAINING_RECONSTRUCTION_PANELS,
    SCHEMA as SCHEMA_V1,
    STRADDLING_REMAINING_PANELS,
    V5_ACTION_COVERED_PANELS,
    assert_coherence_retained,
    build_rho1_coverage,
    load_authenticated_inventory_payload,
    pair_negative_energy,
)
from recursive_horizons.nsc_ks_rho1_source_bound_v2 import (
    MISSING_PRIMITIVE,
    SCHEMA,
    SIDE_ABOVE_UNCOVERED,
    SIDE_BELOW,
    SIDE_V5_ACTION,
    SPLIT_WINDOW_COUNT,
    V1_RECORD,
    assert_v1_immutable,
    attach_reconstruction_bound,
    combine_directed_n_beta,
    correlated_budget_feeds,
    coverage_identity,
    pauli_generator_matches_ks,
    physical_source_error_for_upstream,
    reject_action_as_column_substitution,
    reject_forbidden_bound_method,
    reject_straddling_misclassification,
    rho1_chart,
    row_propagation_coefficients,
    split_straddling_windows,
    validate_rho1_source_bound_v2_record,
    window_identity,
)
from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from recursive_horizons.nsc_ks_upstream_error import packed_upper_fraction, restore_packed_upper
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_transmitting_dirac_domain import S3
import derive_nsc_ks_rho1_source_bound_v1 as V1
import derive_nsc_ks_rho1_source_bound_v2 as D


CONFIG = {'surface_gravity': .24, 'omega': 3.9}


def inventory_and_v5():
    inventory = json.loads((ROOT / 'results/development/nsc-ks-source-inventory.json').read_text())
    v5 = json.loads((ROOT / 'results/development/nsc-incoming-source-update-v5.json').read_text())
    decision = json.loads((ROOT / 'results/development/nsc-ks-low-subgap-decision-v2.json').read_text())
    channels = json.loads(
        (ROOT / 'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    return inventory, v5, decision, channels


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


def contribution(identity, *, n=None, beta=None, method='synthetic-directed',
                 payload=None, **flags):
    row = {
        'identity': identity,
        'N': n,
        'beta': beta,
        'method': method,
        'coherence_retained': True,
        'weights_retained': True,
        'minus_energy_retained': True,
        'degeneracy_retained': True,
        'independent_column_max_used': False,
        'enclosure_payload_sha256': payload,
    }
    row.update(flags)
    return row


def coverage_bundle():
    inventory, v5, decision, channels = inventory_and_v5()
    D.ensure_inventory_payload(inventory)
    arrays, meta, digest = load_authenticated_inventory_payload(ROOT, inventory)
    coverage = build_rho1_coverage(
        decision['coverage']['windows'], inventory_record=inventory,
        v5_coverage=v5['coverage'], channels=channels)
    return inventory, v5, decision, channels, arrays, meta, digest, coverage


def test_v1_bytes_remain_immutable():
    digest = assert_v1_immutable(ROOT)
    saved = json.loads((ROOT / V1_RECORD).read_text())
    assert saved['schema'] == SCHEMA_V1
    assert digest == sha256((ROOT / V1_RECORD).read_bytes()).hexdigest()
    assert saved['remaining_physical_rho1_source_error_bound'] is None
    V1.ensure_inventory_payload(
        json.loads((ROOT / 'results/development/nsc-ks-source-inventory.json').read_text()))
    v1 = V1.calculate()
    assert v1['schema'] == SCHEMA_V1
    assert V1.stable(v1) == V1.stable(saved)


def test_pauli_generator_two_norm_encloses_ks_generator():
    numeric, lower, upper = pauli_generator_matches_ks(1.5, 1.5707963267948966, np.sqrt(5))
    assert lower <= numeric <= upper
    chart = rho1_chart()
    assert float(restore_packed_upper({'mantissa': '1', 'exponent': 0})) == 1.0
    assert float(chart['r']) > 1.4
    assert float(chart['a']) > 0.8


def test_lipschitz_keeps_coherent_frobenius_and_minus_energy():
    panel = synthetic_panel()
    report = assert_coherence_retained(panel.covariance, require_offdiag=True)
    assert report['coherence_retained']
    assert report['offdiag_max'] > 0
    partner = pair_negative_energy(panel, CONFIG)
    np.testing.assert_array_equal(partner.weights, panel.weights)
    np.testing.assert_allclose(
        partner.amplitudes_at_one, S3 @ panel.amplitudes_at_one.conj(), atol=0, rtol=0)
    ones = np.full((2, 3), 1 / np.sqrt(6), dtype=complex)
    coeffs = row_propagation_coefficients(
        2.0, 2 * np.pi, ones, mass=0.0, angular=0.0, multiplicity=1,
        chart=rho1_chart())
    assert coeffs['minus_energy_factor'] == 2
    assert coeffs['weights_retained']
    assert not coeffs['independent_column_max_used']
    # ||A||_F=1, sqrt(w/2pi)=1, mu=1, m=ell=0 -> linear_N = 4 G / a, linear_beta = 4 G
    lin_n = packed_upper_fraction(coeffs['linear_N_upper'])
    lin_b = packed_upper_fraction(coeffs['linear_beta_upper'])
    assert lin_n > 0 and lin_b > 0


def test_action_bound_cannot_substitute_for_column_error():
    v5 = json.loads((ROOT / 'results/development/nsc-incoming-source-update-v5.json').read_text())
    action = v5['partial_error_budget']['covered_spectral_regions_action_error_upper']
    with pytest.raises(ValueError, match='not a rho=1 column'):
        reject_action_as_column_substitution(action[0], action)
    with pytest.raises(ValueError, match='forbidden'):
        attach_reconstruction_bound({}, Fraction(1, 10 ** 12), 'v5_action_as_column',
                                    v5_action_bound=action)
    with pytest.raises(ValueError, match='forbidden'):
        attach_reconstruction_bound({}, Fraction(1, 10 ** 12),
                                    'action_bound_as_column_reconstruction')
    with pytest.raises(ValueError, match='forbidden'):
        reject_forbidden_bound_method('ode_tolerance_as_bound')
    with pytest.raises(ValueError, match='forbidden'):
        reject_forbidden_bound_method('svd_amplification_as_bound')
    with pytest.raises(ValueError, match='4.998e-9'):
        attach_reconstruction_bound({}, HISTORICAL_COARSE_AGGREGATE,
                                    'historical_coarse_4.998e-9')


def test_split_covers_164_and_splits_59():
    _inventory, _v5, _decision, _channels, arrays, _meta, _digest, coverage = coverage_bundle()
    split = split_straddling_windows(arrays, coverage)
    reject_straddling_misclassification(arrays, coverage, split)
    assert split['straddling_remaining_panels'] == STRADDLING_REMAINING_PANELS
    assert split['split_window_count'] == SPLIT_WINDOW_COUNT
    assert split['below_joined_rows'] == 27828
    assert split['at_or_above_joined_uncovered_rows'] == 1528
    assert split['v5_action_rows'] == 20016
    assert (split['below_joined_rows'] + split['at_or_above_joined_uncovered_rows']
            + split['v5_action_rows'] == INVENTORY_POSITIVE_ROWS)
    parents = [w for w in split['windows'] if w.get('straddling_parent')]
    assert len(parents) == 118
    panels = {w['panel'] for w in split['windows']}
    assert len(panels) == INVENTORY_POSITIVE_PANELS
    assert coverage['remaining_reconstruction_panels'] == REMAINING_RECONSTRUCTION_PANELS
    assert coverage['v5_action_covered_panels'] == V5_ACTION_COVERED_PANELS
    below = next(w for w in split['windows']
                 if w['panel'] == 'low/1_1' and w['side'] == SIDE_BELOW)
    key = window_identity(1, 1, 'low/1_1', SIDE_BELOW, below['n_rows'])
    assert key[3] == SIDE_BELOW


def test_straddling_misclassification_is_rejected():
    _inventory, _v5, _decision, _channels, arrays, _meta, _digest, coverage = coverage_bundle()
    split = split_straddling_windows(arrays, coverage)
    forged = deepcopy(split)
    target = next(w for w in forged['windows']
                  if w['panel'] == 'low/1_1' and w['side'] == SIDE_BELOW)
    other = next(w for w in forged['windows']
                 if w['panel'] == 'low/1_1' and w['side'] == SIDE_ABOVE_UNCOVERED)
    target['n_rows'] = target['n_rows'] + other['n_rows']
    target['row_indices'] = target['row_indices'] + other['row_indices']
    forged['windows'] = [w for w in forged['windows'] if w is not other]
    with pytest.raises(ValueError, match='not split'):
        reject_straddling_misclassification(arrays, coverage, forged)
    promoted = deepcopy(split)
    below = next(w for w in promoted['windows']
                 if w['panel'] == 'low/1_1' and w['side'] == SIDE_BELOW)
    below['side'] = SIDE_V5_ACTION
    with pytest.raises(ValueError, match='V5 action'):
        reject_straddling_misclassification(arrays, coverage, promoted)


def test_dropped_coherence_weights_minus_energy_and_duplicates_rejected():
    keys = [coverage_identity(14, 1, 'synthetic', 0),
            coverage_identity(14, 1, 'synthetic', 1)]
    payload = sha256(b'enclosure').hexdigest()
    with pytest.raises(ValueError, match='coherence'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload,
                          coherence_retained=False),
             contribution(keys[1], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload)],
            declared_universe_size=2)
    with pytest.raises(ValueError, match='weights'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload,
                          weights_retained=False),
             contribution(keys[1], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload)],
            declared_universe_size=2)
    with pytest.raises(ValueError, match='minus-energy'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload,
                          minus_energy_retained=False),
             contribution(keys[1], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload)],
            declared_universe_size=2)
    with pytest.raises(ValueError, match='double-counted'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload),
             contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload)],
            declared_universe_size=2)
    with pytest.raises(ValueError, match='independent maxima'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload,
                          independent_column_max_used=True),
             contribution(keys[1], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload)],
            declared_universe_size=2)
    with pytest.raises(ValueError, match='zero-filled'):
        combine_directed_n_beta(
            [contribution(keys[0], n=0, beta=0),
             contribution(keys[1], n=0, beta=0)],
            declared_universe_size=2)
    with pytest.raises(ValueError, match='forbidden'):
        combine_directed_n_beta(
            [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20),
                          method='drop_minus_energy', payload=payload),
             contribution(keys[1], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload)],
            declared_universe_size=2)


def test_combinator_closes_only_with_enclosure_payloads():
    keys = [coverage_identity(14, 1, 'synthetic', 0),
            coverage_identity(14, -1, 'synthetic', 0)]
    payload = sha256(b'authenticated-synthetic-enclosure').hexdigest()
    open_result = combine_directed_n_beta(
        [contribution(keys[0]), contribution(keys[1])],
        declared_universe_size=2, universe_identities=keys)
    assert open_result['N'] is None and open_result['beta'] is None
    assert open_result['status'].startswith('OPEN')
    closed = combine_directed_n_beta(
        [contribution(keys[0], n=Fraction(1, 10), beta=Fraction(1, 20), payload=payload),
         contribution(keys[1], n=Fraction(2, 10), beta=Fraction(3, 20), payload=payload)],
        declared_universe_size=2, universe_identities=keys)
    assert closed['reconstruction_enclosures_bound'] == 2
    assert packed_upper_fraction(closed['N']) >= Fraction(3, 10)
    assert packed_upper_fraction(closed['beta']) >= Fraction(1, 5)
    assert closed['gate_2_closed'] is False
    with pytest.raises(ValueError, match='enclosure payload'):
        physical_source_error_for_upstream((Fraction(1, 10), Fraction(1, 20)))
    feeds = correlated_budget_feeds()
    assert feeds['baseline_low_subgap_N_beta'] is None
    assert feeds['nine_component_feed']['combined_status'].startswith('OPEN')
    assert feeds['covered_regions_not_fed_from_column_error']
    bound = physical_source_error_for_upstream(
        (Fraction(1, 10), Fraction(1, 20)),
        enclosure_payload_sha256=payload,
        universe_digest=closed['identity_digest'])
    assert bound['physical_upstream_N_beta'] is not None
    assert bound['gate_2_closed'] is False


def test_authenticated_universe_trace_and_open_production_record():
    inventory, _v5, _decision, _channels, arrays, meta, digest, coverage = coverage_bundle()
    assert digest == inventory['payload']['sha256']
    split = split_straddling_windows(arrays, coverage)
    from recursive_horizons.nsc_ks_rho1_source_bound_v2 import trace_physical_rho1_universe
    channels = json.loads(
        (ROOT / 'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    trace = trace_physical_rho1_universe(arrays, meta, split, channels)
    assert trace['n_rows'] == INVENTORY_POSITIVE_ROWS
    assert trace['n_panels'] == INVENTORY_POSITIVE_PANELS
    assert trace['reconstruction_enclosures_evaluated_rows'] == 0
    assert trace['propagation_coefficients_evaluated_rows'] == INVENTORY_POSITIVE_ROWS
    assert trace['construction_census']['low_mode_reconstruction']['panels'] == 64
    assert trace['construction_census']['subgap_rows']['rows'] == 304
    assert trace['construction_census']['middle_mode_columns']['panels'] == 62
    assert trace['construction_census']['low_mode_reconstruction']['straddling_panels'] == 59
    assert packed_upper_fraction(trace['propagation_sums']['linear_N_upper']) > 0
    assert meta['panels']['group13/low']['provenance']['method'].startswith('accepted physical')
    assert meta['panels']['group13/low']['group'] == 13
    record = D.calculate()
    D.check({k: v for k, v in record.items() if k != '_payload_bytes'})
    validate_rho1_source_bound_v2_record(
        {k: v for k, v in record.items() if k != '_payload_bytes'})
    public = {k: v for k, v in record.items() if k != '_payload_bytes'}
    assert public['schema'] == SCHEMA
    assert public['status'].startswith('OPEN')
    assert public['missing_primitive'] == MISSING_PRIMITIVE
    assert public['coverage']['straddling_remaining_panels'] == 59
    assert public['coverage']['split_window_count'] == SPLIT_WINDOW_COUNT
    assert public['directed_N_beta_contribution_bound']['N'] is None
    assert public['gate_2_closed'] is False
    assert public['resource_forecast']['generators_launched'] == 0
    assert public['sample']['reconstruction_bounds_evaluated_rows'] == 0
    assert public['propagation']['propagation_coefficients_evaluated_rows'] == 49372
    if D.OUTPUT.is_file():
        saved = json.loads(D.OUTPUT.read_text())
        assert D.stable(public) == D.stable(saved)
        D.check(saved)
    forged = deepcopy(public)
    forged['status'] = 'PASS'
    with pytest.raises(ValueError, match='must remain OPEN'):
        validate_rho1_source_bound_v2_record(forged)
    forged = deepcopy(public)
    forged['gate_2_closed'] = True
    with pytest.raises(ValueError, match='Gate 2'):
        validate_rho1_source_bound_v2_record(forged)
    forged = deepcopy(public)
    forged['directed_N_beta_contribution_bound']['N'] = {
        'mantissa': '0', 'exponent': 0}
    with pytest.raises(ValueError, match='no finite directed'):
        validate_rho1_source_bound_v2_record(forged)
    forged = deepcopy(public)
    forged['reused_successors'] = dict(public['reused_successors'])
    forged['directed_N_beta_contribution_bound'] = dict(
        public['directed_N_beta_contribution_bound'])
    forged['upstream_physical_slot'] = deepcopy(public['upstream_physical_slot'])
    forged['upstream_physical_slot']['physical_upstream_N_beta'] = ['0', '0']
    forged['upstream_physical_slot']['nine_component_feed']['combined_status'] = 'PASS'
    with pytest.raises(ValueError, match='must remain OPEN'):
        validate_rho1_source_bound_v2_record(forged)


def test_runtime_memory_probe_is_cross_platform():
    value, units, source = D.peak_rss()
    assert value is None or value > 0
    assert units in {'bytes', 'kilobytes', 'unavailable'}
    assert isinstance(source, str) and source
