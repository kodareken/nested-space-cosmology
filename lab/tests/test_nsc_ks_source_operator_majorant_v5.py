"""Accepted v4 plus replayed disjoint subgap rows. Partial coverage stays open."""
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
    ROOT / 'results/development/nsc-direct-source-window-v1',
    ROOT / 'results/development/nsc-subgap-row15-upstream-v1.json',
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
    """Use the record that --check already recomputed. Do not replay the proofs."""
    driver = runpy.run_path(str(ROOT / 'scripts/derive_nsc_ks_source_operator_majorant_v5.py'))
    result = json.loads((ROOT / driver['OUTPUT']).read_text())
    archive = RetainedUpstreamArchive(ROOT)
    prior = driver['accepted_v4']()
    return driver, result, archive, prior


def test_new_rows_stay_disjoint_and_partial_counts_match(data):
    driver, result, archive, prior = data
    assert sha256((ROOT / driver['PRIOR']).read_bytes()).hexdigest() == driver['V4_SHA256']
    coverage = result['source_coverage']
    added = coverage['subgap_signed_rows_added']
    witnessed = coverage['subgap_positive_rows_witnessed']
    failed = coverage['subgap_positive_rows_failed']
    missing = coverage['subgap_positive_rows_missing']
    assert added == 2 * witnessed
    assert witnessed + failed + missing == 54
    assert witnessed >= 1
    assert coverage['covered_signed_rows'] == 1042 + added
    assert coverage['remaining_signed_rows_in_family'] == 126 - added
    assert coverage['failed_rows_entered_as_zero'] is False
    assert coverage['above_mass_low32_rows_32_to_40_included'] is False
    assert result['physical_local_gate'] == 'OPEN'
    assert result['physical_upstream_budget_component'] is None
    assert result['source_weights_applied_once'] is True
    assert result['source_quadrature_error_included'] is False
    assert result['covariance_complement_included'] is True
    assert result['negative_bloch_convention'] == '(nx,-ny,-nz)'
    assert result['sign_audit'] == 'c82b600b'
    assert result['status'].startswith('continuous source-error bound')
    assert 'ENCLOSED' not in result['status']
    keys = [(row['panel'], row['row'], row['energy_sign']) for row in coverage['rows']]
    assert len(keys) == len(set(keys)) == coverage['covered_signed_rows']
    present = set(keys)
    assert ('group14/low16_1', 1, 1) in present
    assert ('group14/low16_1', 1, -1) in present
    assert ('group14/low32_1', 32, 1) not in present
    assert ('group14/low32_1', 40, -1) not in present
    base_rows = driver['prior_rows'](prior)
    base = {(row['panel'], row['row'], row['energy_sign']) for row in base_rows}
    new = present - base
    assert base.isdisjoint(new) and len(base) == 1042 and len(new) == added
    assert float(restored_upper(result['continuous_action_error_upper']['N'])) > float(
        restored_upper(prior['continuous_action_error_upper']['N']))
    assert float(restored_upper(result['continuous_action_error_upper']['beta'])) > float(
        restored_upper(prior['continuous_action_error_upper']['beta']))
    source = (ROOT / 'scripts/derive_nsc_ks_source_operator_majorant_v5.py').read_text()
    assert "mode='check'" in source and "mode='resume'" not in source and 'solve_jost' not in source
    rho, entries = driver['require_same_slice'](prior, archive)
    assert rho.hex() == prior['rho_up_hex']
    with ctx.workprec(192):
        _moment, _energy, covered, missing_family = driver['source_error_moments'](
            archive, entries, driver['add_subgap_pairs'](base_rows, []))
    assert len(covered) == 1042 and missing_family == 126
    with pytest.raises(ValueError, match='duplicate'):
        driver['add_subgap_pairs'](base_rows, base_rows[:1])


def test_mutated_identity_weight_or_failure_cannot_enter(data):
    driver, result, archive, prior = data
    base = driver['prior_rows'](prior)
    with pytest.raises(ValueError, match='duplicate'):
        driver['add_subgap_pairs'](base, [base[-1]])
    rho, entries = driver['require_same_slice'](prior, archive)
    joined = driver['add_subgap_pairs'](base, [])
    altered = deepcopy(joined)
    altered[-1]['preparation_digest'] = '0' * 64
    with pytest.raises(ValueError, match='original preparation'):
        driver['source_error_moments'](archive, entries, altered)
    altered = deepcopy(joined)
    altered[-1]['epsilon'] = {'exponent': 0, 'mantissa': '-1'}
    with pytest.raises(ValueError, match='nonnegative'):
        driver['source_error_moments'](archive, entries, altered)
    with pytest.raises(ValueError, match='duplicate'):
        driver['source_error_moments'](archive, entries, joined + [joined[-1]])
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
    omitted = deepcopy(copied)
    omitted['negative_error_copied_from_positive'] = False
    omitted['negative_complement_applied'] = False
    with pytest.raises(ValueError, match='covariance complement'):
        driver['require_saved_witness'](omitted)
    assert result['source_coverage']['v4_signed_rows_reused'] == 1042
