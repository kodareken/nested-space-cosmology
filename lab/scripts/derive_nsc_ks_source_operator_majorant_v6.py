#!/usr/bin/env python3
"""Insert the 18 approved actual18 dyadics into the accepted v4 aggregate.

The negative source identity is a signed-zero equivalence. Both raw digests are
kept. Quadrature weights are applied once. No ODE or direct-window campaign is run.
"""
import argparse, hashlib, json, sys
from pathlib import Path
from flint import arb, ctx
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from derive_nsc_ks_source_operator_majorant import digest, read_bound, source_error_moments
from derive_nsc_ks_source_operator_majorant_v4 import require_same_slice
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_operator_majorant import history_source_insertion
from recursive_horizons.nsc_threshold18_replay import (
    ADJUDICATION, CHECKPOINT_RELATIVE, IMPORT_RELATIVE, IMPORT_TREE_SHA256, ROWS,
    import_tree_sha256, replay_threshold18)
PRIOR = 'results/development/nsc-ks-source-operator-majorant-v4.json'
OUTPUT = 'results/development/nsc-ks-source-operator-majorant-v6.json'
V4_SHA256 = '4b2462987f7784ee77b58687521226f274400ac6b0924bdf493428037aabbdf1'
PRIOR_COVERED, PRIOR_REMAINING = 1042, 126
ADDED_SIGNED = 18
COVERED, REMAINING = 1060, 108
OWNERS = (
    'scripts/derive_nsc_ks_source_operator_majorant_v6.py',
    'scripts/derive_nsc_ks_source_operator_majorant_v4.py',
    'scripts/derive_nsc_ks_source_operator_majorant.py',
    'src/recursive_horizons/nsc_threshold18_replay.py',
    'src/recursive_horizons/nsc_direct_vacuum_source.py',
    'src/recursive_horizons/nsc_source_occupation_enclosure.py',
    'src/recursive_horizons/nsc_ks_signed_state.py',
    'src/recursive_horizons/nsc_ks_source_operator_majorant.py',
    'tests/test_nsc_ks_source_operator_majorant_v6.py')


def accepted_v4():
    raw = (ROOT / PRIOR).read_bytes()
    if hashlib.sha256(raw).hexdigest() != V4_SHA256:
        raise ValueError('accepted v4 source aggregate changed')
    record = read_bound(PRIOR)
    coverage = record.get('source_coverage') if isinstance(record, dict) else None
    if (record.get('schema') != 'NSC-KS-SOURCE-OPERATOR-MAJORANT-v4'
            or record.get('family') != [14, 1] or record.get('physical_local_gate') != 'OPEN'
            or record.get('physical_upstream_budget_component') is not None
            or record.get('source_weights_applied_once') is not True
            or record.get('source_quadrature_error_included') is not False
            or record.get('numerical_field_error_included') is not False
            or record.get('value_assembly_arithmetic_included') is not False
            or record.get('negative_bounds_supplied_separately') is not True
            or record.get('old_middle_pair_replaced_not_added') is not True
            or not isinstance(coverage, dict)
            or coverage.get('covered_signed_rows') != PRIOR_COVERED
            or coverage.get('remaining_signed_rows_in_family') != PRIOR_REMAINING
            or coverage.get('all_source_families') is not False
            or coverage.get('complete_middle_window_replayed') is not True
            or coverage.get('complete_direct_source_window_replayed') is not True
            or coverage.get('v3_signed_rows_reused') != 868
            or coverage.get('direct_signed_rows_added') != 174):
        raise ValueError('accepted v4 family14_1 coverage required')
    return record


def prior_rows(record):
    coverage = record.get('source_coverage') if isinstance(record, dict) else None
    rows = []
    for row in (coverage or {}).get('rows') or []:
        if not isinstance(row, dict) or row.get('covariance_error_upper') is None:
            raise ValueError('missing v4 source bound')
        if row.get('energy_sign') not in (-1, 1) or row.get('bound_owner') is None:
            raise ValueError('missing v4 source bound')
        rows.append({**row, 'angular_sign': row['energy_sign'], 'owner': row['bound_owner'],
                     'epsilon': row['covariance_error_upper']})
    if len(rows) != PRIOR_COVERED:
        raise ValueError('accepted v4 family14_1 coverage required')
    return rows


def threshold_rows(window):
    """Unweighted imported dyadics identified by the local raw digests."""
    if not isinstance(window, dict):
        raise ValueError('replayed threshold18 equivalence required')
    if (window.get('mode') != 'check' or window.get('new_rows_captured') != 0
            or window.get('ode_rerun') is not False or window.get('uses_imported_dyadic') is not True):
        raise ValueError('imported dyadic replay required')
    if window.get('weight_applied_to_error') is not False:
        raise ValueError('source weight already applied')
    equivalence = window.get('equivalence')
    if (not isinstance(equivalence, dict) or equivalence.get('relation') != 'signed_zero_only'
            or equivalence.get('nonzero_numeric_differences') != 0
            or equivalence.get('positive_digests_match') is not True
            or equivalence.get('original_raw_source_digest') == equivalence.get('local_raw_source_digest')
            or equivalence.get('batch_signed_zeros') != 16):
        raise ValueError('signed-zero equivalence required')
    if (window.get('coverage_complete') is not True or window.get('completed_positive_rows') != 9
            or window.get('completed_signed_rows') != ADDED_SIGNED
            or window.get('physical_local_gate') != 'OPEN'
            or window.get('physical_upstream_budget_component') is not None
            or window.get('negative_endpoint_adjudication') != ADJUDICATION
            or window.get('negative_density_endpoint') != '(nx,-ny,-nz)'):
        raise ValueError('complete actual18 census required')
    rows = window.get('rows')
    if not isinstance(rows, list) or len(rows) != 9:
        raise ValueError('complete actual18 census required')
    if [item.get('row') for item in rows] != list(ROWS):
        raise ValueError('actual18 rows changed')
    result = []
    for row in rows:
        if row.get('weight_applied_to_error') is not False or row.get('panel') != 'group14/low32_1':
            raise ValueError('source weight already applied')
        comparisons = row.get('signed_comparisons')
        if (not isinstance(comparisons, list) or len(comparisons) != 2
                or [item.get('energy_sign') for item in comparisons] != [1, -1]
                or [item.get('angular_sign') for item in comparisons] != [1, -1]):
            raise ValueError('separate signed source bounds required')
        if (comparisons[0].get('local_raw_source_digest') == comparisons[1].get('local_raw_source_digest')
                or comparisons[0].get('local_raw_preparation_digest') == comparisons[1].get('local_raw_preparation_digest')):
            raise ValueError('duplicate source identity')
        if (comparisons[1].get('original_raw_source_digest') != equivalence['original_raw_source_digest']
                or comparisons[1].get('local_raw_source_digest') != equivalence['local_raw_source_digest']
                or comparisons[0].get('original_raw_source_digest') != comparisons[0].get('local_raw_source_digest')):
            raise ValueError('signed-zero equivalence required')
        for item in comparisons:
            if item.get('signed_total_upper') is None:
                raise ValueError('missing signed total')
            if item.get('weight_applied_to_error') is not False:
                raise ValueError('source weight already applied')
            result.append({
                'panel': row['panel'], 'row': row['row'], 'energy_sign': item['energy_sign'],
                'angular_sign': item['angular_sign'], 'energy_hex': item['energy_hex'],
                'owner': row['checkpoint'], 'epsilon': item['signed_total_upper'],
                'source_digest': item['local_raw_source_digest'],
                'preparation_digest': item['local_raw_preparation_digest']})
    if len(result) != ADDED_SIGNED:
        raise ValueError('actual18 row census changed')
    return result


def add_threshold_pairs(rows, added):
    keys = {(row['panel'], row['row'], row['energy_sign']) for row in rows}
    if len(keys) != len(rows) or len(keys) != PRIOR_COVERED:
        raise ValueError('duplicate bounded source row')
    for item in added:
        key = (item['panel'], item['row'], item['energy_sign'])
        if key in keys:
            raise ValueError('duplicate actual18 source row')
        if item.get('epsilon') is None:
            raise ValueError('missing signed total')
        keys.add(key)
    if len(added) != ADDED_SIGNED or len(keys) != COVERED:
        raise ValueError('source census changed')
    return list(rows) + list(added)


def bind_import(inputs, window):
    if not isinstance(inputs, dict) or window.get('import_tree_sha256') != IMPORT_TREE_SHA256:
        raise ValueError('dependency binding required')
    root = ROOT / IMPORT_RELATIVE
    bound = {}
    files = sorted(path for path in root.rglob('*') if path.is_file())
    if len(files) != 25:
        raise ValueError('imported threshold18 tree changed')
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        if relative in inputs or relative in bound:
            raise ValueError('threshold18 dependency overlaps an existing binding')
        bound[relative] = digest(relative)
    records = {item['checkpoint'] for item in window.get('rows') or []}
    witnesses = {checkpoint[:-len('record.json')] + 'witness.npz' for checkpoint in records}
    if len(records) != 9 or not records <= set(bound) or not witnesses <= set(bound):
        raise ValueError('unbound threshold18 witness')
    if not all(relative.startswith(CHECKPOINT_RELATIVE + '/rows/') for relative in records):
        raise ValueError('unbound threshold18 witness')
    inputs.update(bound)
    return inputs


def calculate():
    if import_tree_sha256(ROOT) != IMPORT_TREE_SHA256:
        raise ValueError('imported threshold18 tree changed')
    prior = accepted_v4()
    rows = prior_rows(prior)
    archive = RetainedUpstreamArchive(ROOT)
    rho, entries = require_same_slice(prior, archive)
    window = replay_threshold18(ROOT, archive)
    added = threshold_rows(window)
    rows = add_threshold_pairs(rows, added)
    family, identity = require_original_history(ROOT)
    if identity != prior.get('profile_identity'):
        raise ValueError('profile identity changed')
    channel = archive.meta['channels'][14]
    if channel['copy_count'] * channel['degeneracy'] / 2 != 12:
        raise ValueError('multiplicity changed')
    with ctx.workprec(192):
        moment, energy, covered, missing = source_error_moments(archive, entries, rows)
        if len(covered) != COVERED or missing != REMAINING:
            raise ValueError('source census changed')
        bound = history_source_insertion(
            family, rho, moment, energy,
            mass=arb(channel['compact_mass']).union(arb.pi() / 2),
            absolute_angular=arb(channel['angular_eigenvalue']).union(arb(5).sqrt()),
            multiplicity=12)
    return {
        'schema': 'NSC-KS-SOURCE-OPERATOR-MAJORANT-v6',
        'status': 'continuous source-error bound for 1060 signed rows; full upstream component OPEN',
        'profile_identity': identity, 'rho_up_hex': rho.hex(), 'family': [14, 1],
        'multiplicity_per_signed_family': 12,
        'weighted_error_moment': exact_upper(moment),
        'weighted_energy_error_moment': exact_upper(energy),
        'continuous_action_error_upper': bound,
        'negative_source_equivalence': window['equivalence'],
        'threshold18_closure': window['closure'],
        'source_coverage': {
            'covered_signed_rows': COVERED, 'remaining_signed_rows_in_family': REMAINING,
            'rows': covered, 'all_source_families': False,
            'complete_middle_window_replayed': True,
            'complete_direct_source_window_replayed': True,
            'v4_signed_rows_reused': PRIOR_COVERED,
            'threshold18_signed_rows_added': ADDED_SIGNED,
            'threshold18_dyadic_source': 'imported per-row signed_total_upper',
            'ode_rerun': False, 'failed_rows_entered_as_zero': False},
        'sign_audit': ADJUDICATION, 'negative_bloch_convention': '(nx,-ny,-nz)',
        'covariance_complement_included': True, 'negative_bounds_supplied_separately': True,
        'source_weights_applied_once': True, 'source_quadrature_error_included': False,
        'numerical_field_error_included': False, 'value_assembly_arithmetic_included': False,
        'physical_upstream_budget_component': None, 'physical_local_gate': 'OPEN',
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': bind_import(dict(prior['input_hashes']), window)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = calculate()
    encoded = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    if args.record:
        publish_exclusive_file(ROOT, OUTPUT, encoded)
    elif (ROOT / OUTPUT).read_bytes() != encoded:
        raise ValueError('source insertion or dependencies changed')
    action = result['continuous_action_error_upper']
    coverage = result['source_coverage']
    equivalence = result['negative_source_equivalence']
    print(json.dumps({
        'covered_signed_rows': coverage['covered_signed_rows'],
        'remaining_signed_rows_in_family': coverage['remaining_signed_rows_in_family'],
        'N': float(restored_upper(action['N'])),
        'beta': float(restored_upper(action['beta'])),
        'N_mantissa': action['N']['mantissa'], 'N_exponent': action['N']['exponent'],
        'beta_mantissa': action['beta']['mantissa'], 'beta_exponent': action['beta']['exponent'],
        'original_raw_source_digest': equivalence['original_raw_source_digest'],
        'local_raw_source_digest': equivalence['local_raw_source_digest'],
        'v6_sha256': hashlib.sha256(encoded).hexdigest()}, indent=2))
