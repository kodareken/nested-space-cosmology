#!/usr/bin/env python3
"""One geometry-matched endpoint certificate; no propagation or source solve."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_matched_horizon_initializer import (
    matched_ratio_identity, matched_initializer_bound, represent_initializer,
)
from recursive_horizons.nsc_incoming_source_quadrature_bound import unpack,float_enclosure
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_hi

OUTPUT=ROOT/'results/development/nsc-incoming-matched-horizon-initializer.json'
SOURCES=('src/recursive_horizons/nsc_incoming_matched_horizon_initializer.py',
         'scripts/derive_nsc_incoming_matched_horizon_initializer.py',
         'tests/test_nsc_incoming_matched_horizon_initializer.py',
         'docs/nsc-incoming-matched-horizon-initializer.md')
INPUTS=('results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json',
        'results/development/nsc-incoming-horizon-endpoint-bound.json',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
        'src/recursive_horizons/nsc_incoming_middle_bound.py',
        'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
        'src/recursive_horizons/nsc_incoming_state_moments.py',
        'src/recursive_horizons/nsc_paired_horizon_preparation.py',
        'src/recursive_horizons/nsc_unruh_state.py')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signature():return{p:digest(p) for p in(*SOURCES,*INPUTS)}


def prepare():
    before=signature()
    channel=json.loads((ROOT/INPUTS[0]).read_text())['channels'][32]
    config=json.loads((ROOT/INPUTS[1]).read_text())['scattering_provenance']['config']
    data={'signature':before,'channel':channel,'config':config,
          'identity':matched_ratio_identity(),'certificate':matched_initializer_bound(channel,config),
          'point_representations':[represent_initializer(channel,config,E,s) for s in(1,-1) for E in(24,28,32)]}
    if signature()!=before:raise ValueError('matched initializer input changed during evaluation')
    raw=json.dumps(data,sort_keys=True,allow_nan=False).encode();h=sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-matched-horizon-initializer.{h}.json'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-addressed collision')
    if not path.exists():path.write_bytes(raw)
    return path


def make_record(path):
    data=json.loads(path.read_text())
    if data['signature']!=signature() or data['identity']!=matched_ratio_identity():raise ValueError('matched initializer owner/input/identity changed')
    cert=data['certificate']
    def check_intervals(value):
        if isinstance(value,dict):
            if 'binary_interval'in value and float_enclosure(unpack(value['binary_interval']))!=value['float_enclosure']:
                raise ValueError('initializer interval serialization differs')
            for item in value.values():check_intervals(item)
        elif isinstance(value,list):
            for item in value:check_intervals(item)
    with _precision(cert['precision']):
        check_intervals(cert);check_intervals(data['point_representations'])
        total=unpack(cert['total_lapse_error_upper']['binary_interval'])
        parts=sum((unpack(row['lapse_error_upper']['binary_interval']) for row in cert['per_sign']))
        if total._mpi_!=parts._mpi_:raise ValueError('angular-pair bound sum changed')
        passed=_hi(total)<=cert['contextual_remaining_lapse_budget']
        if passed!=cert['mathematical_endpoint_budget_pass']:raise ValueError('matched endpoint gate changed')
    return {'schema':'NSC-GEOMETRY-MATCHED-HORIZON-INITIALIZER-v1','accountable_author':'Douglas Ek',
            'status':('PASS: geometry-matched mathematical vacuum endpoint; propagation OPEN' if passed else 'OPEN: matched initializer enclosure exceeds budget'),
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'payload':{'path':str(path.relative_to(ROOT)),'sha256':digest(path)},
            'identity':data['identity'],'certificate':cert,'point_representations':data['point_representations'],
            'scope':{'new_numerical_approximation_of_same_affine_vacuum':True,
                     'archived_initializers_or_sources_changed':False,'source_kappa_occupation_changed':False,
                     'ODE_source_metric_solves':0,'physical_IV_or_constraint_root':False,
                     'uniform_numeric_error_from_six_samples':False,'push_or_PDF':False},
            'next_connection':'validated unitary propagation of interval initial data on group32 E24–32, with exact-profile coordinate conversion; source quadrature and positive thermal bound remain separate',
            'reproducer':'python3 scripts/derive_nsc_incoming_matched_horizon_initializer.py --check'}


def main():
    p=argparse.ArgumentParser();m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=p.parse_args()
    if args.write:
        if OUTPUT.exists():raise FileExistsError('matched initializer record is not overwritten')
        result=make_record(prepare());OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        old=json.loads(OUTPUT.read_text());path=ROOT/old['payload']['path']
        if digest(path)!=old['payload']['sha256']:raise ValueError('matched initializer artifact changed')
        result=make_record(path)
        if result!=old:raise ValueError('matched initializer replay differs')
    print(json.dumps({'status':result['status'],'lapse_error_upper':result['certificate']['total_lapse_error_upper']['float_enclosure']},indent=2))


if __name__=='__main__':main()
