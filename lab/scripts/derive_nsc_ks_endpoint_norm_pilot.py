#!/usr/bin/env python3
"""Enclose endpoint norms on the existing captured original-source cell."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_endpoint_contraction import endpoint_difference_norms
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
import validate_nsc_ks_current_field_pilot_v3 as capture_owner

OUTPUT='results/development/nsc-ks-endpoint-norm-pilot-v1.json'
OWNERS=(
    'src/recursive_horizons/nsc_ks_endpoint_contraction.py',
    'src/recursive_horizons/nsc_ks_finite_matter_error.py',
    'src/recursive_horizons/nsc_ks_difference_error.py',
    'src/recursive_horizons/nsc_ks_ball_trajectory.py',
    'src/recursive_horizons/nsc_ks_reference_error.py',
    'scripts/derive_nsc_ks_endpoint_norm_pilot.py',
    'scripts/validate_nsc_ks_current_field_pilot_v3.py',
    'tests/test_nsc_ks_endpoint_contraction.py',
    'docs/nsc-ks-endpoint-contraction.md',
)


def calculate():
    capture,arrays=capture_owner.load_capture()
    segment,grid,_,_,period,weights,energies=capture_owner.segment_from_capture(capture,arrays)
    size=2*len(weights)
    difference=segment.end[size:size*(1+len(grid))].reshape(2,len(weights),len(grid))
    norms=endpoint_difference_norms(difference,weights,energies,period,bits=120)
    dependencies=[str(capture_owner.CAPTURE.relative_to(ROOT)),capture['payload']['path']]
    return {
        'schema':'NSC-KS-ENDPOINT-NORM-PILOT-v1',
        'status':'ENCLOSED: continuous spatial norms of captured difference interpolant',
        'profile_identity':capture['profile_identity'],
        'source_preparation_digest':capture['source_preparation_digest'],
        'family':capture['family'],'cell_index':capture['cell_index'],
        'rho_endpoint':segment.rho_end.hex(),
        'bounds':norms,'all_history_cells':False,
        'incoming_constraint_bound':None,'physical_local_gate':'OPEN',
        'source_hashes':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        'input_hashes':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies},
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args();value=calculate()
    if args.record:
        publish_exclusive_file(ROOT,OUTPUT,(json.dumps(value,sort_keys=True,indent=2)+'\n').encode())
    elif json.loads((ROOT/OUTPUT).read_text()) != value:
        raise ValueError('endpoint norm replay differs')
    print(json.dumps({'status':value['status'],
        'display_norms':{k:float(restored_upper(value['bounds'][k])) for k in
                         ('difference_norm','difference_axial_norm')},
        'physical_local_gate':'OPEN'},indent=2))
