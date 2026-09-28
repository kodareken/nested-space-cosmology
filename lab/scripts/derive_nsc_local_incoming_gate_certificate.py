#!/usr/bin/env python3
"""Local incoming-gate certificate. EXISTENCE stays false while any epsilon is missing."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-certificate.json'
REANCHOR = ROOT/'results/development/nsc-ks-gate-reanchor.json'
FEASIBILITY = ROOT/'results/development/nsc-ks-gate-budget-feasibility.json'
STALL = ROOT/'results/development/nsc-ks-gate-search-stall.json'
TOLERANCE = 3e-11
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_certificate.py',
    'docs/nsc-local-incoming-gate-certificate.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    reanchor = json.loads(REANCHOR.read_text())
    feasibility = json.loads(FEASIBILITY.read_text())
    stall = json.loads(STALL.read_text())
    best = next(row for row in reanchor['rows'] if row['name'] == reanchor['campaign_best_name'])
    maxima = best['value_only_constraint_maxima']
    epsilons = {
        'field_space_time': None,
        'changed_history_UV_tail': None,
        'baseline_low_subgap': None,
        'upstream': None,
        'energy_interpolation': None,
        'covered_regions': None,
        'phase_value': None,
        'between_node': None,
        'arithmetic': None,
    }
    if any(value is None for value in epsilons.values()):
        existence = False
    else:
        existence = all(abs(maxima[i]) + sum(epsilons.values()) <= TOLERANCE for i in (0, 1))
    if existence:
        raise ValueError('this owner refuses EXISTENCE while the search stall stands')
    if stall['physical_NONEXISTENCE_certificate']:
        raise ValueError('a numerical-floor stall is not NON-EXISTENCE')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-CERTIFICATE-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: certificate blocked; merit and epsilons miss 3e-11',
        'class': 'delta r = chi(s)[s w(z) + s^3 U(z)/6]',
        'interval': 'S(1)+[0.12,0.18]',
        'state_law': 'C_Sigma[g]=U_g C_up U_g†',
        'profile_identity': best['profile_identity'],
        'nodal_constraint_maxima': maxima,
        'merit': best['value_only_merit'],
        'epsilons': epsilons,
        'tolerance': TOLERANCE,
        'verdict': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'mutations_cannot_pass': ['drop-delta-C', 'drop-coherence', 'double-weights', 'k=-E', 'drop-minus-E-sector'],
        'named_gap': feasibility['named_gap'],
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(REANCHOR.relative_to(ROOT)): digest(REANCHOR),
            str(FEASIBILITY.relative_to(ROOT)): digest(FEASIBILITY),
            str(STALL.relative_to(ROOT)): digest(STALL),
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
            raise FileExistsError('certificate register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('certificate replay differs')
    print(json.dumps({
        'status': result['status'],
        'verdict': result['verdict'],
        'merit': result['merit'],
        'physical_EXISTENCE_certificate': result['physical_EXISTENCE_certificate'],
    }, indent=2))
