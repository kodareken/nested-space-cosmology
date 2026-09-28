#!/usr/bin/env python3
"""0.27.0 stays blocked after the assembled Chebyshev leftover floor."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_local_incoming_gate_status_after_assembled_chebyshev_gn as G

OUTPUT = ROOT/'results/development/nsc-027-blocked-after-assembled-chebyshev-gn.json'
OWNERS = (
    'scripts/derive_nsc_027_blocked_after_assembled_chebyshev_gn.py',
    'docs/nsc-027-blocked-after-assembled-chebyshev-gn.md',
)
PUBLIC = Path('/Users/admin/Documents/nested-space-cosmology')


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    status = json.loads(G.OUTPUT.read_text())
    for path, expected in {**status['source_hashes'], **status['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('0.27.0 Chebyshev-GN block input changed: '+path)
    if status['physical_EXISTENCE_certificate'] or status['physical_NONEXISTENCE_certificate']:
        raise ValueError('a closed certificate would require the public manuscript owners')
    if status['physical_local_gate'] != 'OPEN':
        raise ValueError('this owner records only the OPEN block')
    return {
        'schema': 'NSC-027-BLOCKED-AFTER-ASSEMBLED-CHEBYSHEV-GN-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: 0.27.0, GitHub release and arXiv stay blocked',
        'physical_local_gate': 'OPEN',
        'manuscript_027_written': False,
        'public_release_027_cut': False,
        'arxiv_metadata_prepared': False,
        'public_checkout_unused': True,
        'public_paper_exists': (PUBLIC/'paper/nested-space-cosmology.md').is_file(),
        'named_search_stall': status['named_search_stall'],
        'scope': {
            'open_is_not_done': True,
            'open_is_not_a_paper': True,
            'public_checkout_not_edited_by_this_owner': True,
            'no_class_rewrite': True,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(G.OUTPUT.relative_to(ROOT)): digest(G.OUTPUT),
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
            raise FileExistsError('0.27.0 Chebyshev-GN block exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('0.27.0 Chebyshev-GN block replay differs')
    print(json.dumps({
        'status': result['status'],
        'manuscript_027_written': result['manuscript_027_written'],
        'named_search_stall': result['named_search_stall'],
        'arxiv_metadata_prepared': result['arxiv_metadata_prepared'],
    }, indent=2))
