"""Replay of one archived group14/low16_1 row and its signed partner."""
from copy import deepcopy
from hashlib import sha256
from io import BytesIO
from pathlib import Path
import runpy
import time

import numpy as np
import pytest
from flint import arb
from scipy.integrate._ivp.common import OdeSolution
from scipy.integrate._ivp.rk import Dop853DenseOutput

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_massive_jost_mixed_transport import (
    subgap_mixed_phase_transport_bound)
from recursive_horizons.nsc_massive_jost_transport_bound import original_dense_solution
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_subgap_row_witness import (
    BLOCH_Y_START,
    PHASE_COLUMNS,
    load_witness_arrays,
    phase_prefix_from_dense,
    replay_mode,
    require_exact_source,
    require_same_phase_segments,
    retain_dense_solution,
    run_from_prefix,
    segments_of_prefix,
    signed_source_errors,
    validate_phase_prefix,
)
from recursive_horizons.nsc_transmitting_dirac_domain import S3


DRIVER = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'scripts/derive_nsc_subgap_row15_upstream.py'))
HORIZON = 2.0
INNER = HORIZON + 1.01e-4


def _step(y0, y1, theta0, theta1, correction=0.0):
    corrections = np.zeros(6, dtype=np.float64)
    corrections[0] = correction
    return (float(y0), float(y1), float(theta0), float(theta1), corrections)


def _dense(steps):
    interpolants = []
    for y0, y1, theta0, theta1, corrections in steps:
        coefficients = np.zeros((7, 1), dtype=np.complex128)
        coefficients[0, 0] = theta1 - theta0
        coefficients[1:, 0] = corrections
        y_old = np.array([theta0], dtype=np.complex128)
        interpolants.append(Dop853DenseOutput(y0, y1, y_old, coefficients))
    nodes = np.empty(len(steps) + 1, dtype=np.float64)
    nodes[0] = steps[0][0]
    nodes[1:] = [step[1] for step in steps]
    return OdeSolution(nodes, interpolants)


def _joined_prefix():
    return np.array([
        [1.0, 0.5, 0.25, 0.5, 0.25, 1e-6, 0, 0, 0, 0, 0],
        [0.5, 0.0, 0.5, 0.75, 0.25, -1e-6, 0, 0, 0, 0, 0],
    ], dtype=np.float64)


def _arrays(raw):
    with np.load(BytesIO(raw), allow_pickle=False) as saved:
        return {name: saved[name].copy() for name in saved.files}


def _rehash(record, arrays):
    raw = deterministic_npz_bytes(arrays)
    forged = deepcopy(record)
    forged['payload'] = {
        'path': record['payload']['path'],
        'sha256': sha256(raw).hexdigest(),
        'bytes': len(raw),
    }
    return forged, raw


def _upper(item):
    return float(restored_upper(item))


def test_witness_does_not_bind_the_solvers():
    import recursive_horizons.nsc_subgap_row_witness as witness
    assert 'solve_jost' not in witness.__dict__
    assert 'capture_bloch' not in witness.__dict__


def test_prefix_rejects_broken_joins_f0_and_shape():
    good = _joined_prefix()
    assert validate_phase_prefix(good).shape == (2, 11)
    assert list(PHASE_COLUMNS)[4] == 'F0'
    broken_f0 = good.copy()
    broken_f0[0, 4] += 1e-12
    with pytest.raises(ValueError, match=r'theta0\+F0'):
        validate_phase_prefix(broken_f0)
    broken_radius = good.copy()
    broken_radius[1, 0] += 1e-12
    broken_radius[1, 1] += 1e-12
    with pytest.raises(ValueError, match='identical radii'):
        validate_phase_prefix(broken_radius)
    broken_angle = good.copy()
    broken_angle[1, 2] += 1e-12
    broken_angle[1, 3] += 1e-12
    with pytest.raises(ValueError, match='identical angles'):
        validate_phase_prefix(broken_angle)
    outward = good.copy()
    outward[0, 0], outward[0, 1] = outward[0, 1], outward[0, 0]
    with pytest.raises(ValueError, match='step inward'):
        validate_phase_prefix(outward)
    with pytest.raises(ValueError, match='11 columns'):
        validate_phase_prefix(good[:, :10])
    missing = good.copy()
    missing[0, 2] = np.nan
    with pytest.raises(ValueError, match='finite'):
        validate_phase_prefix(missing)


def test_reconstructed_prefix_matches_original_segments_including_a_discarded_neighbor():
    kept = (
        _step(0.0, -1.0, 0.25, 0.5, 2**-10),
        _step(-1.0, -2.0, 0.5, 0.75, -2**-12),
    )
    discarded = _step(-2.0, -20.0, 0.75, 1.0, 2**-8)
    dense = _dense((*kept, discarded))
    prefix = phase_prefix_from_dense(dense, HORIZON, INNER)
    assert prefix.shape == (2, 11)
    require_same_phase_segments(
        segments_of_prefix(prefix, HORIZON, INNER),
        segments_of_prefix(phase_prefix_from_dense(_dense(kept), HORIZON, INNER), HORIZON, INNER))
    run, recovered = run_from_prefix(prefix)
    assert run.sol.__code__.co_freevars == ('dense',)
    assert len(run.sol.__closure__) == 1
    assert original_dense_solution(run) is recovered
    assert isinstance(recovered, OdeSolution)
    assert all(isinstance(item, Dop853DenseOutput) for item in recovered.interpolants)
    assert retain_dense_solution(recovered).__code__.co_freevars == ('dense',)


def test_mode_label_mismatch_is_rejected_before_transport():
    prefix = validate_phase_prefix(_joined_prefix())
    background = type('Background', (), {'horizon_rho': 1.0})()
    mode = replay_mode(
        prefix, energy=1.0, mass=np.pi / 2, angular=np.sqrt(5.0),
        background=background, outer_radius=60.0)
    with pytest.raises(ValueError, match='energy interval does not contain the Jost label'):
        subgap_mixed_phase_transport_bound(
            0.25, np.pi / 2, np.sqrt(5.0), background, mode=mode,
            inner_radius=1.2, degree=8, metric_terms=8, bits=80)


def test_negative_complement_is_not_copied_from_the_positive_bound():
    columns = np.array([[1, 0, 0], [0, 2, 0]], complex)
    covariance = np.diag([1.0, 0.0, 0.0])
    proof = {
        'bits': 160,
        'endpoint': (arb(0), arb(0), arb(1)),
        'bloch_error': arb(0),
    }
    positive, negative = signed_source_errors(
        columns, covariance, S3 @ columns.conj(), np.eye(3) - covariance.conj(), proof)
    assert positive == 0
    assert negative >= 3
    assert negative < arb('3.00000000000001')


def test_source_array_mismatch_is_rejected_without_a_proof():
    columns = np.zeros((2, 3), dtype=np.complex128)
    covariance = np.zeros((3, 3), dtype=np.complex128)
    payload = {
        'positive_columns': columns,
        'negative_columns': columns,
        'positive_covariance': covariance,
        'negative_covariance': covariance,
    }
    changed = columns.copy()
    changed[0, 0] = 1
    payload['negative_columns'] = changed
    with pytest.raises(ValueError, match='archived upstream source changed in witness: negative_columns'):
        require_exact_source(payload, columns, columns, covariance, covariance)


def test_malformed_payload_inventory_is_rejected():
    arrays = {
        'phase_prefix': _joined_prefix(),
        'bloch_trace': np.array([[BLOCH_Y_START, -17.0, *([0.0] * 24)]], dtype=np.float64),
        'positive_columns': np.zeros((2, 3), dtype=np.complex128),
        'negative_columns': np.zeros((2, 3), dtype=np.complex128),
        'positive_covariance': np.zeros((3, 3), dtype=np.complex128),
        'negative_covariance': np.zeros((3, 3), dtype=np.complex128),
        'outer_radius': np.array([60.0], dtype=np.float64),
        'phase_nfev': np.array([1], dtype=np.int64),
        'bloch_nfev': np.array([1], dtype=np.int64),
    }
    # The synthetic Bloch row above is only an inventory fixture; join checks need a real step.
    arrays['bloch_trace'][0, 1] = BLOCH_Y_START + 0.1
    loaded = load_witness_arrays(arrays)
    assert loaded['phase_prefix'].shape[1] == 11
    missing = dict(arrays)
    del missing['negative_covariance']
    with pytest.raises(ValueError, match='complete row witness payload required'):
        load_witness_arrays(missing)
    bad_shape = dict(arrays)
    bad_shape['positive_columns'] = np.zeros((2, 2), dtype=np.complex128)
    with pytest.raises(ValueError, match='malformed row witness payload: positive_columns'):
        load_witness_arrays(bad_shape)


@pytest.fixture(scope='module')
def evidence():
    audit = {}
    started, wall = time.process_time(), time.perf_counter()
    record, raw = DRIVER['calculate'](audit=audit)
    audit['record'] = record
    audit['raw'] = raw
    audit['cpu_seconds'] = time.process_time() - started
    audit['wall_seconds'] = time.perf_counter() - wall
    return audit


def test_saved_prefix_matches_the_one_captured_row(evidence):
    assert evidence['prefix'].shape == (843, 11)
    assert len(evidence['original_segments']) == 843
    require_same_phase_segments(
        evidence['original_segments'],
        segments_of_prefix(evidence['prefix'], evidence['horizon_rho'], evidence['inner_radius']))
    record = evidence['record']
    assert record['phase_cells'] == 843
    assert record['bloch_cells'] == 671
    assert record['rate_signs'] == {
        'contracting': 820, 'expanding': 23, 'zero': 0, 'straddling': 0}
    assert record['source_panel'] == 'group14/low16_1'
    assert record['source_row'] == 15
    assert record['source_signs'] == [1, -1]
    assert record['angular_signs'] == [1, -1]
    phase = _upper(record['phase_error_inner_upper'])
    positive = _upper(record['positive_covariance_error_upper'])
    negative = _upper(record['negative_covariance_error_upper'])
    assert phase == pytest.approx(1.814146078e-9, abs=5e-18), phase
    assert positive == pytest.approx(1.047331588e-10, abs=5e-20), positive
    assert negative == pytest.approx(1.046705871e-10, abs=5e-20), negative
    assert record['positive_covariance_error_upper'] != record['negative_covariance_error_upper']


def test_initial_errors_are_nonzero_and_the_gate_stays_open(evidence):
    record = evidence['record']
    assert _upper(record['initializer_phase_error_upper']) > 0
    assert _upper(record['initial_bloch_error_upper']) > 0
    assert _upper(record['phase_error_inner_upper']) > 0
    assert record['physical_local_gate'] == 'OPEN'
    assert record['n_beta_aggregate'] is None
    assert record['physical_upstream_budget_component'] is None
    assert record['physical_rho1_source_error'] is None
    assert record['horizon_sewing_error'] is None
    assert record['frame_error_included_in_initial_bloch'] is True
    assert record['separate_frame_error_addition_required'] is False
    assert record['covariance_error_unweighted'] is True
    assert record['quadrature_weight_applied'] is False
    assert record['negative_error_copied_from_positive'] is False
    assert record['negative_complement_applied'] is True
    assert record['source_columns_replaced'] is False
    assert record['all_source_families'] is False
    assert record['amplitude_component_retained'] is False
    assert 'cpu_seconds' not in record
    assert 'wall_seconds' not in record
    assert evidence['cpu_seconds'] > 0
    assert evidence['wall_seconds'] > 0


def test_negative_partner_comes_from_the_archive(evidence):
    binding = DRIVER['load_binding']()
    arrays = _arrays(evidence['raw'])
    assert np.array_equal(arrays['positive_columns'], binding['positive_columns'])
    assert np.array_equal(arrays['negative_columns'], binding['negative_columns'])
    assert np.array_equal(arrays['positive_covariance'], binding['positive_covariance'])
    assert np.array_equal(arrays['negative_covariance'], binding['negative_covariance'])
    assert binding['negative'].angular == -binding['positive'].angular
    assert binding['negative'].energy_sign == -1
    assert not np.array_equal(binding['positive_columns'], binding['negative_columns'])
    assert evidence['positive_angular'] == binding['positive'].angular
    assert evidence['negative_angular'] == binding['negative'].angular


def test_saved_payload_rejects_hash_shape_joins_and_labels(evidence):
    record, raw = evidence['record'], evidence['raw']
    arrays = _arrays(raw)
    stale = deterministic_npz_bytes({**arrays, 'phase_nfev': arrays['phase_nfev'] + 1})
    with pytest.raises(ValueError, match='payload hash'):
        DRIVER['calculate'](saved=(record, stale))
    missing = dict(arrays)
    del missing['bloch_trace']
    forged, mutated = _rehash(record, missing)
    with pytest.raises(ValueError, match='complete row witness payload required'):
        DRIVER['calculate'](saved=(forged, mutated))
    short = dict(arrays)
    short['phase_prefix'] = arrays['phase_prefix'][:, :10].copy()
    forged, mutated = _rehash(record, short)
    with pytest.raises(ValueError, match='11 columns'):
        DRIVER['calculate'](saved=(forged, mutated))
    broken_f0 = dict(arrays)
    broken_f0['phase_prefix'] = arrays['phase_prefix'].copy()
    broken_f0['phase_prefix'][0, 4] += 1e-12
    forged, mutated = _rehash(record, broken_f0)
    with pytest.raises(ValueError, match=r'theta0\+F0'):
        DRIVER['calculate'](saved=(forged, mutated))
    broken_join = dict(arrays)
    broken_join['phase_prefix'] = arrays['phase_prefix'].copy()
    broken_join['phase_prefix'][1, 0] += 1e-12
    broken_join['phase_prefix'][1, 1] += 1e-12
    forged, mutated = _rehash(record, broken_join)
    with pytest.raises(ValueError, match='identical radii'):
        DRIVER['calculate'](saved=(forged, mutated))
    broken_trace = dict(arrays)
    broken_trace['bloch_trace'] = arrays['bloch_trace'].copy()
    broken_trace['bloch_trace'][1, 0] += 1e-8
    forged, mutated = _rehash(record, broken_trace)
    with pytest.raises(ValueError, match='do not join'):
        DRIVER['calculate'](saved=(forged, mutated))
    shifted = dict(arrays)
    shifted['bloch_trace'] = arrays['bloch_trace'].copy()
    shifted['bloch_trace'][0, 0] += 1e-12
    forged, mutated = _rehash(record, shifted)
    with pytest.raises(ValueError, match='initializer'):
        DRIVER['calculate'](saved=(forged, mutated))
    changed_source = dict(arrays)
    changed_source['positive_columns'] = arrays['positive_columns'].copy()
    changed_source['positive_columns'][0, 0] += 1e-12
    forged, mutated = _rehash(record, changed_source)
    with pytest.raises(ValueError, match='archived upstream source'):
        DRIVER['calculate'](saved=(forged, mutated))
    relabeled = deepcopy(record)
    relabeled['positive_energy_hex'] = '0x1p+0'
    with pytest.raises(ValueError, match='source label mismatch'):
        DRIVER['calculate'](saved=(relabeled, raw))
    blank = deepcopy(record)
    del blank['payload']
    with pytest.raises(ValueError, match='complete row witness payload required'):
        DRIVER['calculate'](saved=(blank, raw))


def test_replay_checks_the_saved_witness_without_solving(evidence, monkeypatch, tmp_path):
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError('replay called a solver')
    monkeypatch.setitem(DRIVER['calculate'].__globals__, 'solve_jost', fail_if_called)
    monkeypatch.setitem(DRIVER['calculate'].__globals__, 'capture_bloch', fail_if_called)
    root = Path(tmp_path).resolve()
    (root / 'results/development/artifacts').mkdir(parents=True)
    DRIVER['publish_witness'](root, evidence['record'], evidence['raw'])
    checked = DRIVER['verify_saved'](root)
    assert checked == evidence['record']
    assert (root / DRIVER['OUTPUT']).is_file()
    assert (root / DRIVER['PAYLOAD']).is_file()
    with pytest.raises(FileExistsError, match='already exists'):
        DRIVER['publish_witness'](root, evidence['record'], evidence['raw'])


def test_changed_phase_coefficient_does_not_replay_as_the_original(evidence):
    record, raw = evidence['record'], evidence['raw']
    arrays = _arrays(raw)
    arrays['phase_prefix'] = arrays['phase_prefix'].copy()
    arrays['phase_prefix'][0, 5] += 1e-6
    forged, mutated = _rehash(record, arrays)
    assert forged['payload']['sha256'] != record['payload']['sha256']
    try:
        fresh, fresh_raw = DRIVER['calculate'](saved=(forged, mutated))
    except (ValueError, ArithmeticError) as exc:
        text = str(exc)
        assert (
            'replay differs' in text or 'tube' in text or 'defect' in text
            or 'unresolved' in text or 'not inside' in text)
    else:
        assert fresh != record
        assert fresh_raw != raw
