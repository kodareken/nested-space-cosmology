#!/usr/bin/env python3
"""Audit current evidence with certificate-v2 arithmetic; record only a closed result."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_local_gate_certificate_v2 import build_certificate

OUTPUT = ROOT / 'results/development/nsc-local-incoming-gate-certificate-v2.json'
REANCHOR = ROOT / 'results/development/nsc-ks-gate-reanchor-v2.json'
BUDGET = ROOT / 'results/development/nsc-ks-gate-budget-v4.json'


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def current_open_audit():
    reanchor = json.loads(REANCHOR.read_text())
    budget = json.loads(BUDGET.read_text())
    best = next(row for row in reanchor['rows']
                if row['profile_identity'] == reanchor['best_measured_profile_identity'])
    components = {name: row['bound'] for name, row in budget['components'].items()}
    mutations = {name: False for name in (
        'drop-delta-C', 'drop-coherence', 'double-weights', 'k=-E',
        'drop-minus-E-sector')}
    return build_certificate(
        mode='existence', profile_identity=best['profile_identity'],
        source_identity=reanchor['source_identity'], interval=reanchor['interval'],
        state_law=reanchor['state_law'], history_class=reanchor['class'],
        residual_upper=best['nodal_maxima'], components=components,
        mutations=mutations, numerical_settings={'certificate_numerics': None},
        coverage={'full_retained_source': reanchor['full_retained_source_covered'],
                  'continuous_interval': False, 'node_count': reanchor['node_count']},
        provenance={'input_hashes': {
            str(REANCHOR.relative_to(ROOT)): digest(REANCHOR),
            str(BUDGET.relative_to(ROOT)): digest(BUDGET)}})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--audit-open', action='store_true')
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        value = json.loads(OUTPUT.read_text())
    else:
        value = current_open_audit()
        if args.record:
            if value['verdict'] == 'OPEN':
                raise ValueError('certificate v2 refuses to record an OPEN result')
            publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
                (json.dumps(value, indent=2, sort_keys=True) + '\n').encode())
    print(json.dumps({
        'schema': value['schema'], 'verdict': value['verdict'],
        'open_is_not_done': value['open_is_not_done'],
        'decision': value['decision']}, indent=2))
