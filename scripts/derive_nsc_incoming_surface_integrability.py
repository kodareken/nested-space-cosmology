#!/usr/bin/env python3
"""Authenticate exact incoming reference contractions; no numerical generator."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_surface_integrability import certificate

OUTPUT=ROOT/'results/development/nsc-incoming-surface-integrability.json'
SOURCES=('src/recursive_horizons/nsc_incoming_surface_integrability.py',
         'scripts/derive_nsc_incoming_surface_integrability.py',
         'tests/test_nsc_incoming_surface_integrability.py','docs/nsc-incoming-surface-integrability.md')
INPUTS=('results/development/nsc-reference-band-bulk.json',
        'src/recursive_horizons/nsc_spatial_reference_symbol.py',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
        'src/recursive_horizons/nsc_incoming_reference_response.py',
        'src/recursive_horizons/nsc_reference_band_action.py',
        'results/development/nsc-mode-resolved-cauchy-state.json')


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def make_record():
    uv=json.loads((ROOT/INPUTS[0]).read_text())
    for field in ('source_hashes','input_hashes'):
        for path,expected in uv[field].items():
            if digest(path)!=expected:raise ValueError('imported UV certificate changed: '+path)
    channels=json.loads((ROOT/INPUTS[-1]).read_text())['channels']
    if len(channels)!=33 or any(c['compact_mass']==c['angular_eigenvalue']==0 for c in channels[1:]):
        raise ValueError('32 gapped groups and separately owned LLL required')
    result=certificate(uv)
    return {'schema':'NSC-INCOMING-SURFACE-INTEGRABILITY-v1','accountable_author':'Douglas Ek',
            'status':'PASS: included lapse/shift reference changes integrable for finite compatible incoming jets',
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'result':result,'exact_residuals':'all matrix, trace and retained first-Weyl residuals exactly zero',
            'tolerance':0,'computation':'exact2x2 coefficient algebra; imported existing UV rows; no reference quadrature or probe',
            'reproducer':'python3 scripts/derive_nsc_incoming_surface_integrability.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write and OUTPUT.exists():raise FileExistsError('existing identity record is not overwritten')
    data=make_record()
    if args.write:OUTPUT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    elif json.loads(OUTPUT.read_text())!=data:raise ValueError('exact identity replay changed')
    print(json.dumps({k:data[k] for k in ('status','exact_residuals','tolerance')},indent=2))


if __name__=='__main__':main()
