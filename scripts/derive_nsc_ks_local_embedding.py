#!/usr/bin/env python3
"""Verify the continuum local embedding of the two saved KS numerical cells."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_energy_propagator_control as O
from recursive_horizons.nsc_ks_local_embedding import continuum_local_embedding
from recursive_horizons.nsc_ks_profile_identity import profile_identity

OUTPUT=ROOT/'results/development/nsc-ks-local-embedding.json'
FINE=ROOT/'results/development/nsc-ks-fine-trajectory.json'
RADIUS=ROOT/'results/development/nsc-ks-radius-bound.json'
OWNED=('scripts/derive_nsc_ks_local_embedding.py','src/recursive_horizons/nsc_ks_local_embedding.py',
       'docs/nsc-ks-local-embedding.md','src/recursive_horizons/nsc_ks_ball_geometry.py',
       'src/recursive_horizons/nsc_ks_profile_identity.py')


def digest(path):return O.digest(path)


def compute():
    fine=json.loads(FINE.read_text());radius=json.loads(RADIUS.read_text())
    for record in (fine,radius):
        for category in ('source_hashes','input_hashes'):
            for path,expected in record.get(category,{}).items():
                if digest(path)!=expected:raise ValueError('embedding input dependency changed: '+path)
    if digest(fine['payload']['path'])!=fine['payload']['sha256']:raise ValueError('fine trajectory payload changed')
    with np.load(ROOT/fine['payload']['path'],allow_pickle=False) as f:meta=json.loads(f['metadata_json'].tobytes())
    family=O.C.R.C.P.family(radius['control']['history_amplitude'])
    if meta['history_identity']!=profile_identity(family):raise ValueError('radius certificate and fine history differ')
    positive=Q(radius['bounds']['radius_lower']['exact_rational'])>0
    amin=Q(radius['bounds']['axial_lower']['exact_rational'])
    kwargs=dict(rho_up=meta['rho_up'],rho_sigma=1,axial_lower=amin,radius_positive=positive)
    fine_case=continuum_local_embedding(**kwargs,period_left=meta['period_origin'],period_length=meta['period_length'],
        physical_interval=meta['physical_interval'],axial_support=meta['axial_support'])
    old,oldmeta,_=O.inputs();saved=json.loads(O.OUTPUT.read_text())
    if digest(O.PAYLOAD)!=saved['payload']['sha256']:raise ValueError('saved operator changed')
    with np.load(O.PAYLOAD,allow_pickle=False) as f:op=O.restore({k:f[k] for k in f.files},oldmeta)
    if op.digest!=saved['control']['operator_digest']:raise ValueError('operator binding changed')
    if tuple(op.binding.amplitudes)!=tuple(family.amplitudes):raise ValueError('operator and radius amplitude differ')
    grid=op.binding.computational_z
    length=Q(float(grid[1]))-Q(float(grid[0]));length*=len(grid)
    operator_case=continuum_local_embedding(**kwargs,period_left=float(grid[0]),period_length=length,
        physical_interval=(float(op.z[0]),float(op.z[-1])),axial_support=op.binding.axial_support)
    paths=(FINE,RADIUS,ROOT/fine['payload']['path'],O.OUTPUT,O.PAYLOAD)
    return {'schema':'NSC-KS-LOCAL-EMBEDDING-CONTROL-v1','accountable_author':'Douglas Ek',
        'status':'PASS: continuum local domain embedding; numerical/source and physical gates separate',
        'cases':{'fine_field_cell':fine_case,'retained_operator_cell':operator_case},
        'history_identity':meta['history_identity'],'axial_profile_identity':meta['axial_profile_identity'],
        'source_law':'C_Sigma[g]=U_g C_up U_gdagger; same original upstream covariance',
        'source_hashes':{p:digest(p) for p in OWNED},
        'input_hashes':{str(p.relative_to(ROOT)):digest(p) for p in paths},
        'scope':{'local_interval':'S(1)+[.12,.18]','physical_local_gate':'OPEN',
            'finite_Fourier_causality_certified':False,'field_and_source_error_bound':None,
            'metric_timestep':False,'global_parent_child_matching_required':False,
            'new_field_or_source_runs':0,'physical_periodic_boundary_condition':False},
        'reproducer':'python3 scripts/derive_nsc_ks_local_embedding.py --check'}


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();result=compute()
    if args.record:
        if OUTPUT.exists():raise FileExistsError('embedding record exists; use --check')
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif result!=json.loads(OUTPUT.read_text()):raise ValueError('embedding replay differs')
    print(json.dumps({'status':result['status'],'two_D_buffer_margins':{
        k:[float(Q(v)) for v in row['two_distance_buffer_margin']] for k,row in result['cases'].items()}},indent=2))
