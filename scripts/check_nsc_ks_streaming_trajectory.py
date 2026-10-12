#!/usr/bin/env python3
"""Verify streaming KS trajectory capture; do not run a field/source control."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.nsc_ks_streaming_trajectory import (
    CHUNK_ARRAYS, INTENDED_LATER_CONTROL, MAX_SEGMENTS_PER_CHUNK, SCHEMA,
    StreamSummary, stream_ks_trajectory,
)

OWNED = (
    'src/recursive_horizons/nsc_ks_streaming_trajectory.py',
    'tests/test_nsc_ks_streaming_trajectory.py',
    'docs/nsc-ks-streaming-trajectory.md',
    'scripts/check_nsc_ks_streaming_trajectory.py',
)


def manifest_schema():
    return {
        'schema': SCHEMA,
        'complete': 'boolean; true only after a successful close',
        'max_segments_per_chunk': MAX_SEGMENTS_PER_CHUNK,
        'accepted_steps': 'sum of complete-chunk segment counts',
        'chunk_count': 'number of monotonically numbered NPZ files',
        'state_size': 'complex endpoint length',
        'rho_start': 'first stored node',
        'rho_end': 'last stored node',
        'chunks': [{
            'index': '0-based integer',
            'filename': 'chunk-000000.npz',
            'sha256': 'hex digest of the deterministic NPZ bytes',
            'bytes': 'file length',
            'segment_count': '1..4',
            'rho_start': 'first node of this chunk',
            'rho_end': 'last node of this chunk',
        }],
        'chunk_arrays': list(CHUNK_ARRAYS),
        'public_api': {
            'stream_ks_trajectory': '(*evolve_args, on_segment=callback, **evolve_options) '
                                    '-> (prepared, StreamSummary)',
            'StreamSummary': list(StreamSummary.__dataclass_fields__),
            'on_segment': 'called once per accepted step; streamer retains no segment list',
        },
        'intended_later_control_not_executed': INTENDED_LATER_CONTROL,
        'physical_local_gate': 'OPEN',
        'owned': list(OWNED),
    }


def check():
    env = os.environ.copy()
    env['PYTHONPATH'] = str(ROOT / 'src') + os.pathsep + env.get('PYTHONPATH', '')
    subprocess.check_call(
        [sys.executable, '-m', 'pytest', '-q', 'tests/test_nsc_ks_streaming_trajectory.py'],
        cwd=ROOT, env=env)
    assert callable(stream_ks_trajectory)
    return manifest_schema()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--run', action='store_true')
    args = parser.parse_args()
    if args.run:
        raise SystemExit('this owner does not execute the N=2048 field control')
    result = check()
    print(json.dumps(
        {'status': 'streaming capture self-check passed; physical gate OPEN',
         'schema': result['schema'],
         'public_api': result['public_api'],
         'intended_later_control_not_executed': result['intended_later_control_not_executed']},
        indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
