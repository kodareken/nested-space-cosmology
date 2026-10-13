"""Actual v4/v5 receipt helpers retained by the v7 majorant owner.

These replay local authenticated files only; no source preparation is performed.
Original helper bytes are identified by the explicit 329ef809 producer context.
"""
import argparse, json, sys
from hashlib import sha256
from pathlib import Path
from recursive_horizons.nsc_posix import ROOT
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_subgap_source_window import EXPECTED_POSITIVE_ROWS, FAILURE_SCHEMA, PANEL_ROWS, SCHEMA, checkpoint_relative, cover, failure_relative, settings_for
WINDOW = 'results/development/nsc-subgap-source-window-v1'

def _expected_keys():
    return [(panel, row) for panel, rows in PANEL_ROWS for row in rows]

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
    if window.get('mode') != 'check' or window.get('new_rows_captured') != 0 or window.get('replay_uses_saved_phase_and_bloch_witness') is not True or (window.get('preparation_jost_rerun_on_replay') is not False):
        raise ValueError('replayed subgap window required')
    if window.get('quadrature_weight_applied') is not False or window.get('covariance_errors_unweighted') is not True:
        raise ValueError('subgap source weight already applied')
    if window.get('schema') != SCHEMA or window.get('group') != 14 or window.get('angular_sign') != 1 or (window.get('failed_rows_entered_as_zero') is not False) or (window.get('missing_rows_entered_as_zero') is not False) or (window.get('finite_window_closes_gate') is not False) or (window.get('physical_upstream_budget_component') is not None) or (window.get('physical_local_gate') != 'OPEN') or (window.get('negative_bounds_supplied_separately') is not True) or (window.get('source_columns_replaced') is not False) or (window.get('all_source_families') is not False):
        raise ValueError('replayed subgap window required')
    rows, failed, missing = (window.get('rows'), window.get('failed'), window.get('missing'))
    if not all((isinstance(item, list) for item in (rows, failed, missing))):
        raise ValueError('subgap census does not add up')
    if len(rows) + len(failed) + len(missing) != 54 or window.get('completed_positive_rows') != len(rows) or window.get('completed_signed_rows') != 2 * len(rows) or (window.get('failed_positive_rows') != len(failed)) or (window.get('missing_positive_rows') != len(missing)):
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
        panel, index = (row.get('panel'), row.get('row'))
        if (panel, index) not in set(expected) or type(index) is not int:
            raise ValueError('subgap window panel changed')
        settings = settings_for(panel, index)
        if row.get('transport_profile') != settings['transport_profile'] or row.get('tube') != settings['tube'] or row.get('defect_subdivisions') != settings['defect_subdivisions'] or (row.get('quadrature_weight_applied') is not False) or (row.get('negative_error_copied_from_positive') is not False) or (row.get('physical_local_gate') != 'OPEN') or (row.get('physical_upstream_budget_component') is not None) or (row.get('checkpoint') != checkpoint_relative(panel, index)):
            raise ValueError('subgap witness does not match the established control')
        comparisons = row.get('signed_comparisons')
        fiber, negative = (row.get('energy_hex'), row.get('negative_energy_hex'))
        if not isinstance(comparisons, list) or len(comparisons) != 2 or [item.get('energy_sign') for item in comparisons] != [1, -1] or ([item.get('angular_sign') for item in comparisons] != [1, -1]):
            raise ValueError('separate signed source bounds required')
        try:
            energy = float.fromhex(fiber)
        except (TypeError, ValueError) as error:
            raise ValueError('signed subgap energy mismatch') from error
        if not 0 < energy < float.fromhex('0x1.921fb54442d18p+0') or negative != (-energy).hex():
            raise ValueError('signed subgap energy mismatch')
        if comparisons[0].get('source_digest') == comparisons[1].get('source_digest') or comparisons[0].get('preparation_digest') == comparisons[1].get('preparation_digest'):
            raise ValueError('duplicate source identity')
        if comparisons[0].get('covariance_error_upper') != row.get('positive_covariance_error_upper') or comparisons[1].get('covariance_error_upper') != row.get('negative_covariance_error_upper'):
            raise ValueError('signed subgap bound mismatch')
        owner = (Path(WINDOW) / row['checkpoint'] / 'record.json').as_posix()
        for sign, item, expected_energy in ((1, comparisons[0], fiber), (-1, comparisons[1], negative)):
            if item.get('energy_hex') != expected_energy or item.get('source_digest') is None or item.get('preparation_digest') is None:
                raise ValueError('signed subgap energy mismatch')
            result.append({'panel': panel, 'row': index, 'energy_sign': sign, 'angular_sign': sign, 'energy_hex': expected_energy, 'owner': owner, 'epsilon': _dyadic(item.get('covariance_error_upper'), 'nonnegative covariance bound required'), 'source_digest': item['source_digest'], 'preparation_digest': item['preparation_digest']})
    for item in failed:
        if item.get('entered_as_zero') is not False or item.get('covariance_bound_used_in_aggregate') is not False or item.get('retuned') is not False or (item.get('physical_local_gate') != 'OPEN'):
            raise ValueError('a failed row must not be entered as zero')
    if len(result) != 2 * len(rows):
        raise ValueError('subgap census does not add up')
    return result

def _record_bytes(path):
    try:
        return (ROOT / path).read_bytes()
    except OSError as error:
        raise ValueError('missing subgap window file') from error

def require_saved_witness(record):
    if record.get('negative_complement_applied') is not True or record.get('negative_error_copied_from_positive') is not False or record.get('quadrature_weight_applied') is not False or (record.get('method_too_loose') is not False) or (record.get('covariance_error_unweighted') is not True) or (record.get('physical_local_gate') != 'OPEN'):
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
        record_path = (Path(WINDOW) / row['checkpoint'] / 'record.json').as_posix()
        witness_path = (Path(WINDOW) / row['checkpoint'] / 'witness.npz').as_posix()
        record_bytes, witness_bytes = (_record_bytes(record_path), _record_bytes(witness_path))
        try:
            record = json.loads(record_bytes)
        except json.JSONDecodeError as error:
            raise ValueError('missing subgap window file') from error
        require_saved_witness(record)
        if record.get('proof_sha256') != proof:
            raise ValueError('frozen proof code changed')
        payload = record.get('payload') if isinstance(record, dict) else None
        witness_hash = sha256(witness_bytes).hexdigest()
        if not isinstance(payload, dict) or payload.get('path') != f"{row['checkpoint']}/witness.npz" or payload.get('sha256') != witness_hash or (payload.get('bytes') != len(witness_bytes)):
            raise ValueError('payload hash changed')
        for path, raw in ((record_path, record_bytes), (witness_path, witness_bytes)):
            if path in inputs or path in bound:
                raise ValueError('subgap window dependency overlaps an existing binding')
            bound[path] = sha256(raw).hexdigest()
    for item in window.get('failed') or []:
        path = (Path(WINDOW) / failure_relative(item['panel'], item['row']) / 'failure.json').as_posix()
        raw = _record_bytes(path)
        try:
            record = json.loads(raw)
        except json.JSONDecodeError as error:
            raise ValueError('missing subgap window file') from error
        if record.get('schema') != FAILURE_SCHEMA or record.get('entered_as_zero') is not False or record.get('covariance_bound_used_in_aggregate') is not False or (record.get('retuned') is not False) or (record.get('proof_sha256') != proof):
            raise ValueError('a failed row must not be entered as zero')
        if path in inputs or path in bound:
            raise ValueError('subgap window dependency overlaps an existing binding')
        bound[path] = sha256(raw).hexdigest()
    if freeze_path in inputs:
        raise ValueError('subgap window dependency overlaps an existing binding')
    actual = {path.relative_to(ROOT).as_posix() for path in (ROOT / WINDOW).rglob('*') if path.is_file()}
    if actual != set(bound):
        raise ValueError('unbound subgap window file')
    inputs.update(bound)
    return inputs

def require_same_slice(record, archive):
    rho_hex = record.get('rho_up_hex') if isinstance(record, dict) else None
    if not isinstance(rho_hex, str):
        raise ValueError('same upstream slice required')
    rho = float.fromhex(rho_hex)
    if rho.hex() != rho_hex:
        raise ValueError('same upstream slice required')
    entries = archive.family_entries((14, 1))
    if any((batch.rho_up != rho for batch, _ in entries)):
        raise ValueError('upstream archive slice changed')
    return (rho, entries)
