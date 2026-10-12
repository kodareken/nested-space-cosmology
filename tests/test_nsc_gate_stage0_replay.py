#!/usr/bin/env python3
"""Dirac-free provenance checks for the truncated trail and the ruler."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n32_truncated_iterate as N32
import derive_nsc_ks_coupled_newton_n64_high_modes as N64
import derive_nsc_ks_gate_ruler_drift as Ruler

ITERATE6 = 'b4b53f0bba7290cfde14af741c82fda9af1b49a13f7a3b207a7230f39ea106e9'
TRAIL = (
    'nsc-ks-coupled-newton-n32-truncated-trial',
    'nsc-ks-coupled-newton-n32-truncated-iterate2',
    'nsc-ks-coupled-newton-n32-truncated-iterate3',
    'nsc-ks-coupled-newton-n32-truncated-iterate4',
    'nsc-ks-coupled-newton-n32-truncated-iterate5',
    'nsc-ks-coupled-newton-n32-truncated-iterate6',
)


def test_name_rejects_path_characters():
    for name in ('../etc', 'iterate6.json', 'Iterate', 'a/b', ''):
        try:
            N32.paths(name)
        except ValueError:
            continue
        raise AssertionError('path-like name was accepted: '+name)
    directory, output, proposal = N32.paths('iterate6')
    assert directory.name == 'nsc-ks-coupled-newton-n32-truncated-iterate6'
    assert output.parent == proposal.parent


def test_campaign_best_regenerates_iterate6_identity():
    parent = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate5.json'
    best = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate3.json'
    proposal = N32.build_proposal(parent, 1e6, 1.0, best)
    assert proposal['candidate_profile_identity'] == ITERATE6
    assert proposal['campaign_best_constraint_maxima'] == json.loads(best.read_text())['constraint_maxima']
    assert proposal['immediate_parent_constraint_maxima'] == json.loads(parent.read_text())['constraint_maxima']
    saved = json.loads((ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6-proposal.json').read_text())
    assert np.array_equal(
        np.asarray(proposal['candidate_coefficients'], float),
        np.asarray(saved['candidate_coefficients'], float))
    assert proposal['kept_modes'] == saved['kept_modes']
    assert proposal['step_norm'] == saved['step_norm']


def test_truncated_trail_payloads_match_their_records():
    for name in TRAIL:
        directory = ROOT/'results/development/artifacts'/name
        records = sorted(directory.glob('family-*.json'))
        assert len(records) == 60
        first = json.loads(records[0].read_text())
        payload = ROOT/first['payload']['path']
        assert N32.digest(payload) == first['payload']['sha256']
        for path in records:
            record = json.loads(path.read_text())
            assert (ROOT/record['payload']['path']).stat().st_size == record['payload']['bytes']


def test_iterate6_replay_keeps_the_measured_residual():
    result = N32.replay_name('iterate6', family_key=(1, -1))
    saved = json.loads((ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json').read_text())
    assert result['replay_dirac_evolution'] is False
    assert result['constraint_maxima'] == saved['constraint_maxima']
    assert result['profile_identity'] == ITERATE6


def test_n64_column_replay_keeps_the_rank_table():
    result = N64.replay(family_key=(1, -1))
    saved = json.loads(N64.OUTPUT.read_text())
    assert result['replay_dirac_evolution'] is False
    assert result['ranks'] == saved['ranks']
    assert result['profile_identity'] == ITERATE6


def test_ruler_drift_is_the_saved_matter_sums():
    result = Ruler.compute()
    assert result['families_evolved'] == 0
    assert result['physical_EXISTENCE_certificate'] is False
    assert result['drift_subtracted_as_physics'] is False
    n_drift, beta_drift = result['sum_difference_max_abs']
    assert 5.784e-5 < n_drift < 5.785e-5
    assert 1.284e-10 < beta_drift < 1.285e-10
    assert result['per_family_relative']['median'][0] > 1e-5
    assert result['per_family_relative']['maximum'][0] > 1e-3
    saved = json.loads(Ruler.OUTPUT.read_text())
    assert result == saved
