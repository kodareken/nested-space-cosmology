#!/usr/bin/env python3
"""Insert replayed subgap witnesses into the accepted v4 source aggregate.

Reads the accepted v4 record and replays saved subgap traces. It does not
solve, retune a rejected row, apply a quadrature weight twice, or close the
local gate. Failed and missing rows stay out of the sum.
"""
import argparse, json, sys
from hashlib import sha256
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
from recursive_horizons.nsc_subgap_source_window import (
    EXPECTED_POSITIVE_ROWS, FAILURE_SCHEMA, PANEL_ROWS, SCHEMA, checkpoint_relative,
    cover, failure_relative, settings_for)
WINDOW = 'results/development/nsc-subgap-source-window-v1'
PRIOR = 'results/development/nsc-ks-source-operator-majorant-v4.json'
OUTPUT = 'results/development/nsc-ks-source-operator-majorant-v5.json'
V4_SHA256 = '4b2462987f7784ee77b58687521226f274400ac6b0924bdf493428037aabbdf1'
PRIOR_COVERED, PRIOR_REMAINING = 1042, 126
SIGN_AUDIT = 'c82b600b'
OWNERS = (
    'scripts/derive_nsc_ks_source_operator_majorant_v5.py',
    'scripts/derive_nsc_ks_source_operator_majorant_v4.py',
    'scripts/derive_nsc_ks_source_operator_majorant.py',
    'scripts/derive_nsc_subgap_source_window.py',
    'src/recursive_horizons/nsc_subgap_source_window.py',
    'src/recursive_horizons/nsc_subgap_row_witness.py',
    'src/recursive_horizons/nsc_subgap_upstream_covariance.py',
    'src/recursive_horizons/nsc_ks_source_operator_majorant.py',
    'tests/test_nsc_ks_source_operator_majorant_v5.py')


def _expected_keys():
    return [(panel, row) for panel, rows in PANEL_ROWS for row in rows]


def accepted_v4():
    """The published v4 file, unchanged, with its declared dependencies."""
    raw = (ROOT / PRIOR).read_bytes()
    if sha256(raw).hexdigest() != V4_SHA256:
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


def _dyadic(value, label):
    if not isinstance(value, dict) or not isinstance(value.get('mantissa'), str):
        raise ValueError(label)
    try:
        upper = restored_upper(value)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(label) from error
    if not upper.is_finite() or not upper >= 0:
        raise ValueError(label)
    return value


def subgap_rows(window):
    """Signed pairs from a replay. Rejected rows are not given a zero bound."""
    if EXPECTED_POSITIVE_ROWS != 54 or SCHEMA != 'NSC-SUBGAP-SOURCE-WINDOW-v1':
        raise ValueError('subgap window census changed')
    expected = _expected_keys()
    if len(expected) != 54 or len(set(expected)) != 54:
        raise ValueError('subgap window census changed')
    if not isinstance(window, dict):
        raise ValueError('replayed subgap window required')
    if (window.get('mode') != 'check' or window.get('new_rows_captured') != 0
            or window.get('replay_uses_saved_phase_and_bloch_witness') is not True
            or window.get('preparation_jost_rerun_on_replay') is not False):
        raise ValueError('replayed subgap window required')
    if window.get('quadrature_weight_applied') is not False or window.get('covariance_errors_unweighted') is not True:
        raise ValueError('subgap source weight already applied')
    if (window.get('schema') != SCHEMA or window.get('group') != 14 or window.get('angular_sign') != 1
            or window.get('failed_rows_entered_as_zero') is not False
            or window.get('missing_rows_entered_as_zero') is not False
            or window.get('finite_window_closes_gate') is not False
            or window.get('physical_upstream_budget_component') is not None
            or window.get('physical_local_gate') != 'OPEN'
            or window.get('negative_bounds_supplied_separately') is not True
            or window.get('source_columns_replaced') is not False
            or window.get('all_source_families') is not False):
        raise ValueError('replayed subgap window required')
    rows, failed, missing = window.get('rows'), window.get('failed'), window.get('missing')
    if not all(isinstance(item, list) for item in (rows, failed, missing)):
        raise ValueError('subgap census does not add up')
    if (len(rows) + len(failed) + len(missing) != 54
            or window.get('completed_positive_rows') != len(rows)
            or window.get('completed_signed_rows') != 2 * len(rows)
            or window.get('failed_positive_rows') != len(failed)
            or window.get('missing_positive_rows') != len(missing)):
        raise ValueError('subgap census does not add up')
    seen = {(item.get('panel'), item.get('row')) for item in rows}
    failed_keys = {(item.get('panel'), item.get('row')) for item in failed}
    missing_keys = {(item.get('panel'), item.get('row')) for item in missing}
    if seen & failed_keys or seen & missing_keys or failed_keys & missing_keys:
        raise ValueError('a row cannot be both witnessed and failed')
    if seen | failed_keys | missing_keys != set(expected):
        raise ValueError('subgap census does not add up')
    result = []
    for row in rows:
        panel, index = row.get('panel'), row.get('row')
        if (panel, index) not in set(expected) or type(index) is not int:
            raise ValueError('subgap window panel changed')
        settings = settings_for(panel, index)
        if (row.get('transport_profile') != settings['transport_profile']
                or row.get('tube') != settings['tube']
                or row.get('defect_subdivisions') != settings['defect_subdivisions']
                or row.get('quadrature_weight_applied') is not False
                or row.get('negative_error_copied_from_positive') is not False
                or row.get('physical_local_gate') != 'OPEN'
                or row.get('physical_upstream_budget_component') is not None
                or row.get('checkpoint') != checkpoint_relative(panel, index)):
            raise ValueError('subgap witness does not match the established control')
        comparisons = row.get('signed_comparisons')
        fiber, negative = row.get('energy_hex'), row.get('negative_energy_hex')
        if (not isinstance(comparisons, list) or len(comparisons) != 2
                or [item.get('energy_sign') for item in comparisons] != [1, -1]
                or [item.get('angular_sign') for item in comparisons] != [1, -1]):
            raise ValueError('separate signed source bounds required')
        try:
            energy = float.fromhex(fiber)
        except (TypeError, ValueError) as error:
            raise ValueError('signed subgap energy mismatch') from error
        if not 0 < energy < float.fromhex('0x1.921fb54442d18p+0') or negative != (-energy).hex():
            raise ValueError('signed subgap energy mismatch')
        if (comparisons[0].get('source_digest') == comparisons[1].get('source_digest')
                or comparisons[0].get('preparation_digest') == comparisons[1].get('preparation_digest')):
            raise ValueError('duplicate source identity')
        if (comparisons[0].get('covariance_error_upper') != row.get('positive_covariance_error_upper')
                or comparisons[1].get('covariance_error_upper') != row.get('negative_covariance_error_upper')):
            raise ValueError('signed subgap bound mismatch')
        owner = str(Path(WINDOW) / row['checkpoint'] / 'record.json')
        for sign, item, expected_energy in ((1, comparisons[0], fiber), (-1, comparisons[1], negative)):
            if (item.get('energy_hex') != expected_energy or item.get('source_digest') is None
                    or item.get('preparation_digest') is None):
                raise ValueError('signed subgap energy mismatch')
            result.append({
                'panel': panel, 'row': index, 'energy_sign': sign, 'angular_sign': sign,
                'energy_hex': expected_energy, 'owner': owner,
                'epsilon': _dyadic(item.get('covariance_error_upper'), 'nonnegative covariance bound required'),
                'source_digest': item['source_digest'],
                'preparation_digest': item['preparation_digest']})
    for item in failed:
        if (item.get('entered_as_zero') is not False or item.get('covariance_bound_used_in_aggregate') is not False
                or item.get('retuned') is not False or item.get('physical_local_gate') != 'OPEN'):
            raise ValueError('a failed row must not be entered as zero')
    if len(result) != 2 * len(rows):
        raise ValueError('subgap census does not add up')
    return result


def add_subgap_pairs(rows, added):
    keys = {(row['panel'], row['row'], row['energy_sign']) for row in rows}
    if len(keys) != len(rows) or len(keys) != PRIOR_COVERED:
        raise ValueError('duplicate bounded source row')
    for item in added:
        key = (item['panel'], item['row'], item['energy_sign'])
        if key in keys:
            raise ValueError('duplicate subgap source row')
        if item.get('epsilon') is None:
            raise ValueError('missing subgap source bound')
        keys.add(key)
    if len(added) % 2 or len(keys) != PRIOR_COVERED + len(added):
        raise ValueError('source census changed')
    return list(rows) + list(added)


def _record_bytes(path):
    try:
        return (ROOT / path).read_bytes()
    except OSError as error:
        raise ValueError('missing subgap window file') from error


def require_saved_witness(record):
    if (record.get('negative_complement_applied') is not True
            or record.get('negative_error_copied_from_positive') is not False
            or record.get('quadrature_weight_applied') is not False
            or record.get('method_too_loose') is not False
            or record.get('covariance_error_unweighted') is not True
            or record.get('physical_local_gate') != 'OPEN'):
        raise ValueError('covariance complement required')
    return record


def bind_subgap_inputs(inputs, window):
    if not isinstance(inputs, dict) or not isinstance(window, dict):
        raise ValueError('dependency binding required')
    bound = {}
    freeze_path = f'{WINDOW}/proof-freeze.json'
    freeze_bytes = _record_bytes(freeze_path)
    try:
        freeze = json.loads(freeze_bytes)
    except json.JSONDecodeError as error:
        raise ValueError('missing subgap window file') from error
    proof = freeze.get('sha256') if isinstance(freeze, dict) else None
    if not isinstance(proof, str) or len(proof) != 64:
        raise ValueError('frozen proof code changed')
    bound[freeze_path] = sha256(freeze_bytes).hexdigest()
    for row in window.get('rows') or []:
        record_path = str(Path(WINDOW) / row['checkpoint'] / 'record.json')
        witness_path = str(Path(WINDOW) / row['checkpoint'] / 'witness.npz')
        record_bytes, witness_bytes = _record_bytes(record_path), _record_bytes(witness_path)
        try:
            record = json.loads(record_bytes)
        except json.JSONDecodeError as error:
            raise ValueError('missing subgap window file') from error
        require_saved_witness(record)
        if record.get('proof_sha256') != proof:
            raise ValueError('frozen proof code changed')
        payload = record.get('payload') if isinstance(record, dict) else None
        witness_hash = sha256(witness_bytes).hexdigest()
        if (not isinstance(payload, dict) or payload.get('path') != f"{row['checkpoint']}/witness.npz"
                or payload.get('sha256') != witness_hash or payload.get('bytes') != len(witness_bytes)):
            raise ValueError('payload hash changed')
        for path, raw in ((record_path, record_bytes), (witness_path, witness_bytes)):
            if path in inputs or path in bound:
                raise ValueError('subgap window dependency overlaps an existing binding')
            bound[path] = sha256(raw).hexdigest()
    for item in window.get('failed') or []:
        path = str(Path(WINDOW) / failure_relative(item['panel'], item['row']) / 'failure.json')
        raw = _record_bytes(path)
        try:
            record = json.loads(raw)
        except json.JSONDecodeError as error:
            raise ValueError('missing subgap window file') from error
        if (record.get('schema') != FAILURE_SCHEMA or record.get('entered_as_zero') is not False
                or record.get('covariance_bound_used_in_aggregate') is not False
                or record.get('retuned') is not False or record.get('proof_sha256') != proof):
            raise ValueError('a failed row must not be entered as zero')
        if path in inputs or path in bound:
            raise ValueError('subgap window dependency overlaps an existing binding')
        bound[path] = sha256(raw).hexdigest()
    if freeze_path in inputs:
        raise ValueError('subgap window dependency overlaps an existing binding')
    actual = {path.relative_to(ROOT).as_posix()
              for path in (ROOT / WINDOW).rglob('*') if path.is_file()}
    if actual != set(bound):
        raise ValueError('unbound subgap window file')
    inputs.update(bound)
    return inputs


def calculate():
    prior = accepted_v4()
    rows = prior_rows(prior)
    archive = RetainedUpstreamArchive(ROOT)
    rho, entries = require_same_slice(prior, archive)
    window = cover(ROOT, ROOT / WINDOW, mode='check')
    added = subgap_rows(window)
    rows = add_subgap_pairs(rows, added)
    family, identity = require_original_history(ROOT)
    if identity != prior.get('profile_identity'):
        raise ValueError('profile identity changed')
    channel = archive.meta['channels'][14]
    if channel['copy_count'] * channel['degeneracy'] / 2 != 12:
        raise ValueError('multiplicity changed')
    covered_count = PRIOR_COVERED + len(added)
    remaining = PRIOR_REMAINING - len(added)
    with ctx.workprec(192):
        moment, energy, covered, missing = source_error_moments(archive, entries, rows)
        if len(covered) != covered_count or missing != remaining:
            raise ValueError('source census changed')
        bound = history_source_insertion(
            family, rho, moment, energy,
            mass=arb(channel['compact_mass']).union(arb.pi() / 2),
            absolute_angular=arb(channel['angular_eigenvalue']).union(arb(5).sqrt()),
            multiplicity=12)
    witnessed = window['completed_positive_rows']
    failed = window['failed_positive_rows']
    missing_rows = window['missing_positive_rows']
    if witnessed + failed + missing_rows != 54:
        raise ValueError('subgap census does not add up')
    return {
        'schema': 'NSC-KS-SOURCE-OPERATOR-MAJORANT-v5',
        'status': (
            f'continuous source-error bound for {covered_count} signed rows; '
            f'{witnessed} new positive subgap rows witnessed, {failed} rejected, '
            f'{missing_rows} not witnessed; full upstream component OPEN'),
        'profile_identity': identity, 'rho_up_hex': rho.hex(), 'family': [14, 1],
        'multiplicity_per_signed_family': 12,
        'weighted_error_moment': exact_upper(moment),
        'weighted_energy_error_moment': exact_upper(energy),
        'continuous_action_error_upper': bound,
        'source_coverage': {
            'covered_signed_rows': covered_count,
            'remaining_signed_rows_in_family': remaining,
            'rows': covered, 'all_source_families': False,
            'complete_middle_window_replayed': True,
            'complete_direct_source_window_replayed': True,
            'v4_signed_rows_reused': PRIOR_COVERED,
            'subgap_signed_rows_added': len(added),
            'subgap_positive_rows_witnessed': witnessed,
            'subgap_positive_rows_failed': failed,
            'subgap_positive_rows_missing': missing_rows,
            'subgap_window_complete': window['coverage_complete'] is True,
            'failed_rows_entered_as_zero': False,
            'above_mass_low32_rows_32_to_40_included': False},
        'sign_audit': SIGN_AUDIT,
        'negative_bloch_convention': '(nx,-ny,-nz)',
        'covariance_complement_included': True,
        'negative_bounds_supplied_separately': True,
        'source_weights_applied_once': True,
        'source_quadrature_error_included': False,
        'numerical_field_error_included': False,
        'value_assembly_arithmetic_included': False,
        'physical_upstream_budget_component': None,
        'physical_local_gate': 'OPEN',
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': bind_subgap_inputs(dict(prior['input_hashes']), window)}


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
    print(json.dumps({
        'status': result['status'],
        'N': float(restored_upper(action['N'])),
        'beta': float(restored_upper(action['beta'])),
        'N_mantissa': action['N']['mantissa'], 'N_exponent': action['N']['exponent'],
        'beta_mantissa': action['beta']['mantissa'], 'beta_exponent': action['beta']['exponent'],
        'covered_signed_rows': coverage['covered_signed_rows'],
        'remaining_signed_rows_in_family': coverage['remaining_signed_rows_in_family'],
        'subgap_signed_rows_added': coverage['subgap_signed_rows_added'],
        'subgap_positive_rows_failed': coverage['subgap_positive_rows_failed'],
        'subgap_positive_rows_missing': coverage['subgap_positive_rows_missing'],
    }, indent=2))
