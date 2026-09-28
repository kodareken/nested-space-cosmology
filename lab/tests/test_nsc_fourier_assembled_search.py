"""Assembled Fourier search stays short of EXISTENCE and does not rewrite the class."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def test_assembled_fourier_search_refuses_evolve():
    search = json.loads(
        (ROOT / 'results/development/nsc-ks-fourier-assembled-search.json').read_text())
    leftover = np.asarray(search['all_node_unclipped_leftover'], float)
    assert search['acts_on_assembled_residual'] is True
    assert search['not_geometry_leftover_alone'] is True
    assert search['no_invented_matter_columns'] is True
    assert search['assembled_transfer_faithful'] is False
    assert search['n8']['slot_transfer_relative_residual'] > 1e-2
    assert search['evolve_authorized'] is False
    assert search['families_evolved'] == 0
    assert search['unclipped_reaches_existence'] is False
    assert np.all(leftover > 3e-11)
    assert search['named_search_stall'] == 'finite_wu_geometry_image_cannot_reach_existence'
    assert search['scope']['s5_not_next'] is True
    assert search['scope']['no_class_rewrite'] is True
    forbidden = set(search['forbidden_identities'])
    assert '2d2588c3c9b1c83fe1a3f291314b8030300c4f14d3b12c24eadd2664eb20231b' in forbidden
    assert '949a188181bd9ba0df43bc8ea00db2c89ce35a458017052b7af928059f819680' in forbidden


def test_gate_stays_open_and_paper_blocked():
    status = json.loads(
        (ROOT / 'results/development/nsc-local-incoming-gate-status-after-assembled-search.json').read_text())
    bounds = json.loads(
        (ROOT / 'results/development/nsc-ks-deciding-bounds-after-assembled-search.json').read_text())
    blocked = json.loads(
        (ROOT / 'results/development/nsc-027-blocked-after-assembled-search.json').read_text())
    assert status['physical_local_gate'] == 'OPEN'
    assert status['physical_EXISTENCE_certificate'] is False
    assert status['physical_NONEXISTENCE_certificate'] is False
    assert status['s5_retracted'] is True
    assert status['named_search_stall'] == 'finite_wu_geometry_image_cannot_reach_existence'
    assert 's^5' not in (status['next_named_owner'] or '')
    assert 'no class rewrite' in status['next_named_owner']
    assert status['families_evolved'] == 0
    assert bounds['bounds_spent'] is False
    assert bounds['changed_history_UV_tail'] is None
    assert bounds['field_error_bound'] is None
    assert bounds['baseline_low_subgap'] is None
    assert blocked['manuscript_027_written'] is False
    assert blocked['arxiv_metadata_prepared'] is False
    assert blocked['public_checkout_unused'] is True
