#!/usr/bin/env python3
"""Small exact endpoint contraction record; no archived calculation reruns."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.nsc_ks_linear_uv_contractions import contraction_identity

OUTPUT = ROOT / 'results/development/nsc-ks-linear-uv-contractions.json'
OWNERS = ('src/recursive_horizons/nsc_ks_linear_uv_contractions.py',
          'docs/nsc-ks-linear-uv-contractions.md', 'scripts/derive_nsc_ks_linear_uv_contractions.py')
INPUTS = ('src/recursive_horizons/nsc_incoming_fixed_transfer.py',
          'results/development/nsc-incoming-fixed-transfer.json',
          'src/recursive_horizons/nsc_incoming_surface_integrability.py',
          'results/development/nsc-incoming-surface-integrability.json')


def record():
    digest = lambda p: sha256((ROOT / p).read_bytes()).hexdigest()
    return {'schema': 'NSC-KS-LINEAR-UV-CONTRACTIONS-v1', 'accountable_author': 'Douglas Ek',
            'status': 'PASS: exact linear endpoint coefficients; finite-history local gate OPEN',
            'result': contraction_identity(), 'source_hashes': {p: digest(p) for p in OWNERS},
            'input_hashes': {p: digest(p) for p in INPUTS},
            'reproducer': 'python3 scripts/derive_nsc_ks_linear_uv_contractions.py --check'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--record', action='store_true')
    group.add_argument('--check', action='store_true')
    result = record()
    if parser.parse_args().record:
        if OUTPUT.exists():
            raise FileExistsError('existing contraction record; use --check')
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    elif json.loads(OUTPUT.read_text()) != result:
        raise ValueError('linear contraction identity or owner binding changed')
    print(json.dumps({'status': result['status'], 'residuals': result['result']['residuals']}, indent=2))
