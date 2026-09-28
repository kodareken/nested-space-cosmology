"""The current low/subgap decision reuses completed middle-window bounds."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_low_subgap_decision_v2 as D


def test_successor_reuses_group12_and_group32_and_keeps_low_columns_open():
    value = D.calculate()
    assert value['reused_successors']['group12_order24_lapse_error_upper'] < 3e-11
    assert max(value['reused_successors']['group32_combined_action_error_upper']) < 3e-11
    assert value['remaining_uncovered_reconstruction_window_count'] > 0
    assert value['remaining_uncovered_reconstruction_rows'] > 0
    assert value['remaining_low_subgap_source_error_bound'] is None
    assert not value['historical_coarse_4p998e_minus9_is_current_remaining_bound']
    assert not value['redo_group12_or_group32_order24_middle_bounds']
