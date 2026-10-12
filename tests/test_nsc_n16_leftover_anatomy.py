"""Leftover anatomy refuses a sign-stable necessary relation and a 1e-4 evolve."""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def _module(name):
    scripts = str(ROOT / 'scripts')
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    return importlib.import_module(name)


def test_classify_splits_collocation_svd_and_necessary_relation():
    anatomy = _module('derive_nsc_ks_n16_leftover_anatomy')
    truncated = [
        {'control_max': 0.05, 'kept_modes': 17},
        {'control_max': 0.004, 'kept_modes': 25},
        {'control_max': 0.0005, 'kept_modes': 30},
    ]
    well_conditioned = np.geomspace(1.0, 1e-3, 32)
    assert anatomy.classify(1e-12, 6.3e-4, 3e-12, truncated, well_conditioned) == (
        'all_node_collocation')
    tiny_modes = np.array([1.0] + [1e-9] * 31)
    assert anatomy.classify(6.26e-4, 6.34e-4, 3e-12, truncated, tiny_modes) == (
        'truncated_svd')
    no_tiny = np.geomspace(1.0, 1e-4, 32)
    far_truncated = [
        {'control_max': 0.05, 'kept_modes': 32},
        {'control_max': 0.02, 'kept_modes': 32},
        {'control_max': 0.01, 'kept_modes': 32},
    ]
    assert anatomy.classify(6.26e-4, 6.34e-4, 3e-12, far_truncated, no_tiny) == (
        'necessary_relation_candidate')
    assert anatomy.classify(6.26e-4, 6.34e-4, 6.26e-4, far_truncated, no_tiny) == (
        'freed_jet')


def test_recorded_identification_refuses_nonexistence_and_evolve():
    identification = json.loads(
        (ROOT / 'results/development/nsc-ks-n16-leftover-identification.json').read_text())
    proposal = json.loads(
        (ROOT / 'results/development/nsc-ks-n16-freed-jet-proposal.json').read_text())
    leftover = np.asarray(proposal['unclipped_all_node_leftover'], float)
    assert identification['walk_projection_sign_stable'] is False
    assert identification['physical_NONEXISTENCE_certificate'] is False
    assert identification['evolve_authorized'] is False
    assert proposal['evolve_authorized'] is False
    assert np.all(leftover > 3e-11)


def test_deciding_bounds_stay_none_against_a_1e4_leftover():
    bounds = json.loads(
        (ROOT / 'results/development/nsc-ks-deciding-bounds-after-leftover.json').read_text())
    assert bounds['bounds_spent'] is False
    assert bounds['changed_history_UV_tail'] is None
    assert bounds['field_error_bound'] is None
    assert bounds['baseline_low_subgap'] is None
    assert bounds['full_between_node_remainder'] is None
    assert bounds['geometry_between_node_remainder'] == [1e-13, 1e-13]
