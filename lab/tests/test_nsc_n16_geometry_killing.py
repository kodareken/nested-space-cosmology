"""L(E) and n=32 geometry leftover refuse a certificate and an evolve."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def test_geometry_killing_current_is_not_sign_stable():
    killing = json.loads(
        (ROOT / 'results/development/nsc-ks-n16-geometry-killing-current.json').read_text())
    assert killing['walk_sign_stable'] is False
    assert killing['physical_NONEXISTENCE_certificate'] is False
    assert killing['L_pieces']['baseline']['sign_stable'] is True
    assert killing['L_pieces']['assembled']['sign_stable'] is False
    assert killing['walk_L_min'] < 0 < killing['walk_L_max']


def test_n32_geometry_leftover_refuses_evolve():
    leftover = json.loads(
        (ROOT / 'results/development/nsc-ks-n32-geometry-leftover.json').read_text())
    refusal = json.loads(
        (ROOT / 'results/development/nsc-ks-n32-evolve-refusal.json').read_text())
    unclipped = np.asarray(leftover['all_node_unclipped_leftover'], float)
    assert leftover['evolve_authorized'] is False
    assert leftover['all_node_rank_or_conditioning'] is True
    assert refusal['families_evolved'] == 0
    assert np.all(unclipped > 3e-11)
    truncated = leftover['truncated_svd_all_node']
    assert all(row['control_max'] > 3e-11 for row in truncated)


def test_deciding_bounds_stay_none_after_n32():
    bounds = json.loads(
        (ROOT / 'results/development/nsc-ks-deciding-bounds-after-n32.json').read_text())
    blocked = json.loads(
        (ROOT / 'results/development/nsc-027-blocked-after-n32.json').read_text())
    assert bounds['bounds_spent'] is False
    assert bounds['changed_history_UV_tail'] is None
    assert blocked['manuscript_027_written'] is False
    assert blocked['arxiv_metadata_prepared'] is False
