#!/usr/bin/env python3
"""Validate the reference defect on the authenticated family-14 cell 122."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_reference_residual import reference_segment_defect
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
import validate_nsc_ks_current_field_pilot_v3 as capture_owner

OUTPUT = 'results/development/nsc-ks-reference-residual-pilot-v1.json'
OWNERS = (
    'src/recursive_horizons/nsc_ks_reference_residual.py',
    'src/recursive_horizons/nsc_ks_upstream_error.py',
    'src/recursive_horizons/nsc_ks_ball_geometry.py',
    'src/recursive_horizons/nsc_ks_ball_trajectory.py',
    'src/recursive_horizons/nsc_ks_trajectory.py',
    'scripts/derive_nsc_ks_reference_residual_pilot.py',
    'scripts/validate_nsc_ks_current_field_pilot_v3.py',
    'tests/test_nsc_ks_reference_residual.py',
)


def calculate():
    capture, arrays = capture_owner.load_capture()
    segment, _, _, _, _, weights, energies = capture_owner.segment_from_capture(capture,arrays)
    bound = reference_segment_defect(segment, energies, weights,
                                    capture['mass'],capture['angular'],bits=120)
    dependencies = [str(capture_owner.CAPTURE.relative_to(ROOT)),capture['payload']['path']]
    return {
        'schema': 'NSC-KS-REFERENCE-RESIDUAL-PILOT-v1',
        'status': 'ENCLOSED: numerical homogeneous reference defect on one saved cell',
        'family': capture['family'], 'cell_index': capture['cell_index'],
        'profile_identity': capture['profile_identity'],
        'source_preparation_digest': capture['source_preparation_digest'],
        'bounds': bound, 'all_history_cells': False,
        'physical_source_error_included': False, 'physical_local_gate': 'OPEN',
        'source_hashes': {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        'input_hashes': {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies},
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=calculate()
    if args.record:
        publish_exclusive_file(ROOT,OUTPUT,(json.dumps(result,sort_keys=True,indent=2)+'\n').encode())
    elif json.loads((ROOT/OUTPUT).read_text()) != result:
        raise ValueError('reference residual replay differs')
    print(json.dumps({'status':result['status'],'cell_index':result['cell_index'],
        'integral_upper_display':float(restored_upper(result['bounds']['weighted_row_residual_integral'])),
        'physical_local_gate':'OPEN'},indent=2))
