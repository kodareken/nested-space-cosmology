"""Accepted v6 plus the disjoint replayed subgap rows. The family census closes; the gate stays open."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import runpy
import pytest
from flint import ctx
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
ROOT = Path(__file__).resolve().parents[1]
HISTORIC = (
    ROOT / 'results/development/nsc-ks-source-operator-majorant-v4.json',
    ROOT / 'results/development/nsc-ks-source-operator-majorant-v5.json',
    ROOT / 'results/development/nsc-ks-source-operator-majorant-v6.json',
    ROOT / 'results/development/nsc-subgap-source-window-v1',
    ROOT / 'results/development/nsc-threshold18-import-v1',
    ROOT / 'results/development/nsc-direct-source-window-v1',
    ROOT / '.source-history',
    ROOT / 'src/recursive_horizons/nsc_ks_energy_propagator.py',
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
        if item.is_symlink():
            digest.update(b'link')
        elif item.is_file():
            digest.update(item.read_bytes())
        else:
            digest.update(b'dir')
        digest.update(b'\0')
    return digest.hexdigest()


@pytest.fixture(scope='module', autouse=True)
def historic_bytes_stay_unchanged():
    before = {str(path): _snapshot(path) for path in HISTORIC}
    yield
    assert {str(path): _snapshot(path) for path in HISTORIC} == before


@pytest.fixture(scope='module')
def data():
    driver = runpy.run_path(str(ROOT / 'scripts/derive_nsc_ks_source_operator_majorant_v7.py'))
    result = driver['calculate']()
    encoded = driver['encoded_result'](result)
    output = ROOT / driver['OUTPUT']
    if output.exists():
        assert output.read_bytes() == encoded
    else:
        driver['publish_exclusive_file'](ROOT, driver['OUTPUT'], encoded)
    prior = json.loads((ROOT / driver['PRIOR']).read_text())
    v5 = json.loads((ROOT / driver['V5']).read_text())
    archive = RetainedUpstreamArchive(ROOT)
    return driver, result, prior, v5, archive, encoded


def test_family_census_joins_disjoint_rows_once(data):
    driver, result, prior, v5, _archive, encoded = data
    coverage = result['source_coverage']
    assert coverage['covered_signed_rows'] == 1168
    assert coverage['remaining_signed_rows_in_family'] == 0
    assert coverage['missing_keys'] == []
    assert coverage['v6_signed_rows_reused'] == 1060
    assert coverage['v4_signed_rows_inside_v6'] == 1042
    assert coverage['threshold18_signed_rows_retained'] == 18
    assert coverage['subgap_signed_rows_added'] == 108
    assert coverage['subgap_positive_rows_witnessed'] == 54
    assert coverage['subgap_positive_rows_failed'] == 0
    assert coverage['subgap_positive_rows_missing'] == 0
    assert coverage['v5_covered_total_not_added'] is True
    assert coverage['subgap_proof_replay_repeated'] is False
    assert coverage['ode_rerun'] is False
    assert coverage['failed_rows_entered_as_zero'] is False
    assert coverage['family14_1_signed_rows_complete'] is True
    assert coverage['all_source_families'] is False
    assert coverage['subgap_tree_sha256'] == driver['SUBGAP_TREE_SHA256']
    assert driver['subgap_tree_sha256']() == driver['SUBGAP_TREE_SHA256']
    assert sha256(encoded).hexdigest() == sha256((ROOT / driver['OUTPUT']).read_bytes()).hexdigest()
    assert result['physical_local_gate'] == 'OPEN'
    assert result['physical_upstream_budget_component'] is None
    assert result['source_weights_applied_once'] is True
    assert result['source_quadrature_error_included'] is False
    assert result['numerical_field_error_included'] is False
    assert result['value_assembly_arithmetic_included'] is False
    assert result['negative_bounds_supplied_separately'] is True
    assert 'ENCLOSED' not in result['status']
    assert coverage['covered_signed_rows'] != (
        prior['source_coverage']['covered_signed_rows'] + v5['source_coverage']['covered_signed_rows'])
    keys = [(row['panel'], row['row'], row['energy_sign']) for row in coverage['rows']]
    assert len(keys) == len(set(keys)) == 1168
    v6_keys = {(row['panel'], row['row'], row['energy_sign']) for row in prior['source_coverage']['rows']}
    v5_keys = {(row['panel'], row['row'], row['energy_sign']) for row in v5['source_coverage']['rows']}
    present = set(keys)
    assert v6_keys <= present and len(v6_keys) == 1060
    added = present - v6_keys
    assert added == v5_keys - v6_keys and len(added) == 108
    assert v6_keys - v5_keys == {('group14/low32_1', row, sign) for row in range(32, 41) for sign in (1, -1)}
    assert ('group14/low32_1', 32, 1) in v6_keys and ('group14/low32_1', 32, 1) not in added
    assert ('group14/low16_1', 1, 1) in added and ('subgap8/14_1', 7, -1) in added
    equivalence = result['negative_source_equivalence']
    assert equivalence['relation'] == 'signed_zero_only'
    assert equivalence['nonzero_numeric_differences'] == 0
    assert equivalence['original_raw_source_digest'] == '63e52bb2986c13bd72d91f3965a6f48b3d11380a6d29ae8326658c4fb2bc8d2c'
    assert equivalence['local_raw_source_digest'] == '641b187a3c37aa47a5841b559163eba5f78ec3ca2c79b535dc21dae8db657bdf'
    assert equivalence['original_raw_source_digest'] != equivalence['local_raw_source_digest']
    assert equivalence['original_raw_preparation_digest'] == 'af92bd35174e8d76f368c0aa3f031d630cff06aa128fa28c1e06c4392624ed81'
    assert equivalence['local_raw_preparation_digest'] == '91a903bab1a082cb8a01e8f28c03f745d6f3240b882442d55df3a87db4d1eccf'
    resolution = equivalence['propagator_resolution']
    assert resolution['blob'] == '31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227'
    assert resolution['sha256'] == '5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94'
    assert resolution['sha256'] != resolution['head_sha256']
    assert result['threshold18_closure']['historical_objects'] == [resolution]
    assert result['threshold18_closure']['source_hashes_checked'] == 67
    assert result['threshold18_closure']['archive_digests_checked'] == 64
    assert float(restored_upper(result['continuous_action_error_upper']['N'])) > float(
        restored_upper(prior['continuous_action_error_upper']['N']))
    assert float(restored_upper(result['continuous_action_error_upper']['beta'])) > float(
        restored_upper(prior['continuous_action_error_upper']['beta']))
    assert float(restored_upper(result['continuous_action_error_upper']['N'])) > float(
        restored_upper(v5['continuous_action_error_upper']['N']))
    source = (ROOT / 'scripts/derive_nsc_ks_source_operator_majorant_v7.py').read_text()
    assert 'solve_jost' not in source and 'capture_row' not in source and "cover(" not in source
    assert "mode='resume'" not in source


def test_overlap_weight_and_identity_mutations_are_refused(data):
    driver, result, prior, _v5, archive, _encoded = data
    base = driver['prior_rows'](prior)
    with pytest.raises(ValueError, match='duplicate'):
        driver['add_disjoint_rows'](base, base[:1])
    with pytest.raises(ValueError, match='duplicate'):
        driver['add_disjoint_rows'](base + base[:1], [])
    rho, entries = driver['require_same_slice'](prior, archive)
    assert rho.hex() == prior['rho_up_hex'] == result['rho_up_hex']
    equivalence = result['negative_source_equivalence']
    rebound, rebound_count = driver['rebind_signed_zero_spelling'](base, entries)
    assert rebound_count == 530
    assert sum(row['source_digest'] == equivalence['local_raw_source_digest'] for row in rebound) == 0
    covered_rows = []
    for row in result['source_coverage']['rows']:
        covered_rows.append({**row, 'angular_sign': row['energy_sign'], 'owner': row['bound_owner'],
                             'epsilon': row['covariance_error_upper']})
    with ctx.workprec(192):
        driver['refuse_mutated_identity'](archive, entries, covered_rows)
        _moment, _energy, covered, missing = driver['source_error_moments'](archive, entries, rebound)
    assert len(covered) == 1060 and missing == 108
    window = {
        'schema': driver['SCHEMA'], 'mode': 'resume', 'new_rows_captured': 1,
        'replay_uses_saved_phase_and_bloch_witness': True,
        'preparation_jost_rerun_on_replay': False, 'quadrature_weight_applied': False,
        'covariance_errors_unweighted': True, 'group': 14, 'angular_sign': 1,
        'failed_rows_entered_as_zero': False, 'missing_rows_entered_as_zero': False,
        'finite_window_closes_gate': False, 'physical_upstream_budget_component': None,
        'physical_local_gate': 'OPEN', 'negative_bounds_supplied_separately': True,
        'source_columns_replaced': False, 'all_source_families': False,
        'rows': [], 'failed': [], 'missing': [],
        'completed_positive_rows': 0, 'completed_signed_rows': 0,
        'failed_positive_rows': 0, 'missing_positive_rows': 0}
    with pytest.raises(ValueError, match='replayed'):
        driver['subgap_rows'](window)
    weighted = deepcopy(window)
    weighted['mode'] = 'check'
    weighted['new_rows_captured'] = 0
    weighted['quadrature_weight_applied'] = True
    with pytest.raises(ValueError, match='weight'):
        driver['subgap_rows'](weighted)
    copied = {'negative_complement_applied': True, 'negative_error_copied_from_positive': True,
              'quadrature_weight_applied': False, 'method_too_loose': False,
              'covariance_error_unweighted': True, 'physical_local_gate': 'OPEN'}
    with pytest.raises(ValueError, match='covariance complement'):
        driver['require_saved_witness'](copied)
