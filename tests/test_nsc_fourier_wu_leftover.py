"""Fourier leftover, flux currents and 0.27.0 stay short of a certificate."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def test_fourier_leftover_refuses_evolve():
    leftover = json.loads(
        (ROOT / 'results/development/nsc-ks-fourier-wu-leftover.json').read_text())
    refusal = json.loads(
        (ROOT / 'results/development/nsc-ks-fourier-wu-evolve-refusal.json').read_text())
    unclipped = np.asarray(leftover['all_node_unclipped_leftover'], float)
    assert leftover['evolve_authorized'] is False
    assert leftover['representation'] == 'Fourier-on-I times owned axial plateau'
    assert leftover['selected_mode_count'] in (8, 16)
    assert leftover['all_node_condition_number'] < 1e16
    assert refusal['families_evolved'] == 0
    assert np.all(unclipped > 3e-11)
    assert all(row['control_max'] > 3e-11 for row in leftover['truncated_svd_all_node'])


def test_flux_first_integral_is_not_locked():
    flux = json.loads(
        (ROOT / 'results/development/nsc-ks-flux-first-integral.json').read_text())
    assert flux['leftover_v_not_reused'] is True
    assert flux['L_E_not_reevaluated'] is True
    assert flux['physical_NONEXISTENCE_certificate'] is False
    assert flux['walk_integral_sign_stable'] is False
    assert flux['unenclosed_edge_tail'] is None
    assert flux['edge_can_eat_gap'] is True
    assert flux['walk_integral_min'] < 0 < flux['walk_integral_max']


def test_gate_stays_open_and_paper_blocked():
    status = json.loads(
        (ROOT / 'results/development/nsc-local-incoming-gate-status-after-fourier.json').read_text())
    bounds = json.loads(
        (ROOT / 'results/development/nsc-ks-deciding-bounds-after-fourier.json').read_text())
    blocked = json.loads(
        (ROOT / 'results/development/nsc-027-blocked-after-fourier.json').read_text())
    assert status['physical_local_gate'] == 'OPEN'
    assert status['physical_EXISTENCE_certificate'] is False
    assert status['physical_NONEXISTENCE_certificate'] is False
    assert status['named_search_stall'] == 'declared_wu_image_cannot_reach_existence'
    assert bounds['bounds_spent'] is False
    assert bounds['changed_history_UV_tail'] is None
    assert blocked['manuscript_027_written'] is False
    assert blocked['arxiv_metadata_prepared'] is False
    assert blocked['public_checkout_unused'] is True
