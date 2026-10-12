"""Identity, overlap, null and weight checks for the actual18 v6 insertion."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import runpy
import numpy as np
import pytest
from flint import ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_threshold18_replay import IMPORT_TREE_SHA256, import_tree_sha256, require_signed_zero_only
ROOT = Path(__file__).resolve().parents[1]
HISTORIC = (
    ROOT / 'results/development/nsc-ks-source-operator-majorant-v4.json',
    ROOT / 'results/development/nsc-threshold18-import-v1',
    ROOT / 'src/recursive_horizons/nsc_direct_vacuum_source.py',
    ROOT / 'scripts/derive_nsc_ks_source_operator_majorant_v5.py',
)


def _snapshot(path):
    path = Path(path)
    digest = sha256()
    if path.is_file():
        digest.update(path.read_bytes())
        return digest.hexdigest()
    for item in sorted(path.rglob('*'), key=lambda candidate: candidate.relative_to(path).as_posix()):
        digest.update(item.relative_to(path).as_posix().encode())
        digest.update(b'\0')
        if item.is_file():
            digest.update(item.read_bytes())
        digest.update(b'\0')
    return digest.hexdigest()


@pytest.fixture(scope='module', autouse=True)
def historic_bytes_stay_unchanged():
    before = {str(path): _snapshot(path) for path in HISTORIC}
    yield
    assert {str(path): _snapshot(path) for path in HISTORIC} == before
    assert import_tree_sha256(ROOT) == IMPORT_TREE_SHA256


@pytest.fixture(scope='module')
def data():
    driver = runpy.run_path(str(ROOT / 'scripts/derive_nsc_ks_source_operator_majorant_v6.py'))
    archive = RetainedUpstreamArchive(ROOT)
    window = driver['replay_threshold18'](ROOT, archive)
    prior = driver['accepted_v4']()
    result = driver['calculate']()
    return driver, archive, window, prior, result


def test_signed_zero_identity_is_explicit_and_rows_stay_disjoint(data):
    driver, _archive, window, prior, result = data
    equivalence = result['negative_source_equivalence']
    assert equivalence['original_raw_source_digest'] == '63e52bb2986c13bd72d91f3965a6f48b3d11380a6d29ae8326658c4fb2bc8d2c'
    assert equivalence['local_raw_source_digest'] == '641b187a3c37aa47a5841b559163eba5f78ec3ca2c79b535dc21dae8db657bdf'
    assert equivalence['original_raw_preparation_digest'] == 'af92bd35174e8d76f368c0aa3f031d630cff06aa128fa28c1e06c4392624ed81'
    assert equivalence['local_raw_preparation_digest'] == '91a903bab1a082cb8a01e8f28c03f745d6f3240b882442d55df3a87db4d1eccf'
    assert equivalence['nonzero_numeric_differences'] == 0
    assert equivalence['batch_signed_zeros'] == 16
    assert equivalence['propagator_resolution']['blob'] == '31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227'
    assert equivalence['propagator_resolution']['sha256'] == '5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94'
    assert equivalence['propagator_resolution']['head_sha256'] != equivalence['propagator_resolution']['sha256']
    assert result['threshold18_closure']['source_hashes_checked'] == 67
    assert result['threshold18_closure']['archive_digests_checked'] == 64
    coverage = result['source_coverage']
    assert coverage['covered_signed_rows'] == 1060
    assert coverage['remaining_signed_rows_in_family'] == 108
    assert coverage['v4_signed_rows_reused'] == 1042
    assert coverage['threshold18_signed_rows_added'] == 18
    assert coverage['ode_rerun'] is False
    assert result['source_weights_applied_once'] is True
    keys = [(row['panel'], row['row'], row['energy_sign']) for row in coverage['rows']]
    assert len(keys) == len(set(keys)) == 1060
    added = {('group14/low32_1', index, sign) for index in range(32, 41) for sign in (1, -1)}
    assert added <= set(keys)
    base = {(row['panel'], row['row'], row['energy_sign']) for row in driver['prior_rows'](prior)}
    assert base.isdisjoint(added) and len(base) == 1042
    by_key = {(row['panel'], row['row'], row['energy_sign']): row for row in coverage['rows']}
    sample = by_key[('group14/low32_1', 32, -1)]
    assert sample['source_digest'] == equivalence['local_raw_source_digest']
    assert sample['source_digest'] != equivalence['original_raw_source_digest']
    assert sample['measure_upper'] != sample['covariance_error_upper']
    touched = np.array([[0.0]], float)
    touched[0, 0] = 1e-12
    with pytest.raises(ValueError, match='numerical source arrays differ'):
        require_signed_zero_only(touched, np.zeros((1, 1)), 'mutated')
    assert 'capture_direct_vacuum' not in Path(driver['__file__']).read_text()
    assert 'capture_direct_vacuum' not in (ROOT / 'src/recursive_horizons/nsc_threshold18_replay.py').read_text()


def test_overlap_null_and_weight_are_refused(data):
    driver, _archive, window, prior, _result = data
    base, added = driver['prior_rows'](prior), driver['threshold_rows'](window)
    with pytest.raises(ValueError, match='duplicate'):
        driver['add_threshold_pairs'](base, added + added[:1])
    with pytest.raises(ValueError, match='duplicate'):
        driver['add_threshold_pairs'](base + added[:1], added)
    missing = deepcopy(window)
    missing['rows'][0]['signed_comparisons'][1]['signed_total_upper'] = None
    with pytest.raises(ValueError, match='missing'):
        driver['threshold_rows'](missing)
    weighted = deepcopy(window)
    weighted['weight_applied_to_error'] = True
    with pytest.raises(ValueError, match='weight'):
        driver['threshold_rows'](weighted)
    row_weighted = deepcopy(window)
    row_weighted['rows'][0]['signed_comparisons'][0]['weight_applied_to_error'] = True
    with pytest.raises(ValueError, match='weight'):
        driver['threshold_rows'](row_weighted)
    copied = deepcopy(window)
    copied['rows'][0]['signed_comparisons'][1]['local_raw_source_digest'] = (
        copied['rows'][0]['signed_comparisons'][0]['local_raw_source_digest'])
    with pytest.raises(ValueError, match='duplicate'):
        driver['threshold_rows'](copied)
    dropped = deepcopy(window)
    dropped['equivalence']['original_raw_source_digest'] = dropped['equivalence']['local_raw_source_digest']
    with pytest.raises(ValueError, match='signed-zero'):
        driver['threshold_rows'](dropped)
    short = deepcopy(window)
    short['rows'] = short['rows'][:-1]
    short['completed_positive_rows'] = 8
    short['completed_signed_rows'] = 16
    with pytest.raises(ValueError, match='census'):
        driver['threshold_rows'](short)
