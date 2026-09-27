#!/usr/bin/env python3
"""Record/check the bounded exact controls; never close the physical gate."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_nested_qualities import exact_controls

OUTPUT=ROOT/'results/development/nsc-nested-qualities-v1.json'
OWNERS=('src/recursive_horizons/nsc_nested_qualities.py',
        'scripts/check_nsc_nested_qualities.py','tests/test_nsc_nested_qualities.py',
        'docs/nsc-nested-qualities.md')
REUSED=('results/development/nsc-local-observer-correspondence.json',
        'results/nsc-3-regulated-recursion.json')

def record():
    return {'schema':'NSC-NESTED-QUALITIES-v1','status':'PASS_FINITE_OPERATOR_CONTROLS',
            'claim':'finite-window inheritance and nested local response',
            'physical_local_gate':'OPEN','cosmological_identification':False,
            'eternity_claimed':False,'new_physical_laws_claimed':False,
            'controls':exact_controls(),
            'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
            'reused_record_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in REUSED}}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    a=p.parse_args()
    if a.check and not OUTPUT.is_file():raise FileNotFoundError(OUTPUT)
    value=record()
    if a.record:
        with OUTPUT.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')
    elif value!=json.loads(OUTPUT.read_text()):raise ValueError('nested-quality evidence differs')
    print(json.dumps({'status':value['status'],'exact_controls':len(value['controls']['checks']),
                      'physical_local_gate':'OPEN'},indent=2))

if __name__=='__main__':main()
