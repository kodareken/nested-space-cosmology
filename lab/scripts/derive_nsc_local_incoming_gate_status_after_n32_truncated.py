#!/usr/bin/env python3
"""Best measured n=32 truncated residual. Still OPEN."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
BEST = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate3.json'
WORSE = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate4.json'
OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-n32-truncated.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_n32_truncated.py',
    'docs/nsc-local-incoming-gate-status-after-n32-truncated.md',
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def compute():
    best = json.loads(BEST.read_text())
    worse = json.loads(WORSE.read_text())
    measured = np.asarray(best['constraint_maxima'], float)
    if not best['residual_improved_on_all_components']:
        raise ValueError('best truncated iterate must improve on both components')
    if worse['residual_improved_on_all_components']:
        raise ValueError('a later both-component improvement would replace this best point')
    if np.any(measured <= EXISTENCE_TOLERANCE):
        raise ValueError('a component at 3e-11 would require an enclosed-epsilon certificate')
    gap = measured - EXISTENCE_TOLERANCE
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-N32-TRUNCATED-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: best measured n=32 truncated residual is still above 3e-11'),
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'history_class': 'local source-fixed pure-radius compatible histories',
        'interval': 'S(1)+[.12,.18]',
        'profile_identity': best['profile_identity'],
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'best_constraint_maxima': measured.tolist(),
        'best_prediction_minus_measured': best['prediction_minus_measured'],
        'later_step_improved_both_components': False,
        'later_constraint_maxima': worse['constraint_maxima'],
        'gap_above_existence_tolerance': gap.tolist(),
        'families_evolved_at_best': best['families_evolved'] if 'families_evolved' in best else best['completed_positive_families'],
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_search_stall': 'finite_wu_n32_truncated_image_above_existence',
        'next_named_owner': 'same declared (w,U) class on I; no class rewrite',
        's5_retracted': True,
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'prediction_match_is_not_existence': True,
            'no_class_rewrite': True,
            'metric_timestep': False,
        },
        'source_hashes': {path: digest(ROOT/path) for path in OWNERS},
        'input_hashes': {
            str(BEST.relative_to(ROOT)): digest(BEST),
            str(WORSE.relative_to(ROOT)): digest(WORSE),
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('n=32 truncated gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=32 truncated gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'best_constraint_maxima': result['best_constraint_maxima'],
        'best_prediction_minus_measured': result['best_prediction_minus_measured'],
        'physical_local_gate': result['physical_local_gate'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
