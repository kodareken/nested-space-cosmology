#!/usr/bin/env python3
"""Enclose one original rho=1 source covariance without replacing its bytes."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import acb,arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_metric_horizon_frame import metric_horizon_frame
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_subgap_source_covariance import (
    BlochSource,initial_bloch,capture_bloch,validate_bloch,original_covariance_error,
    original_covariance_distance_bounds,
)

OUTPUT='results/development/nsc-subgap-source-covariance-v1.json'
PAYLOAD='results/development/artifacts/nsc-subgap-source-covariance-v1.npz'
HORIZON='results/development/nsc-metric-horizon-frame-v1.json'
INVENTORY='results/development/nsc-ks-source-inventory.json'
OWNERS=('scripts/derive_nsc_subgap_source_covariance.py',
        'src/recursive_horizons/nsc_subgap_source_covariance.py',
        'tests/test_nsc_subgap_source_covariance.py','docs/nsc-subgap-source-covariance.md')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def ball(value):
    return restored_upper(value['lower']).union(restored_upper(value['upper']))


def calculate():
    horizon=json.loads((ROOT/HORIZON).read_text())
    inventory=json.loads((ROOT/INVENTORY).read_text())
    for group in ('source_hashes','input_hashes'):
        for path,expected in horizon[group].items():
            if digest(path)!=expected:raise ValueError('horizon dependency changed: '+path)
    payload=inventory['payload']
    if digest(payload['path'])!=payload['sha256']:
        raise ValueError('original inventory payload changed')
    panel='group14/low16_1';row=0
    with np.load(ROOT/payload['path'],allow_pickle=False) as data:
        metadata=json.loads(data['metadata_json'].tobytes());config=metadata['config']
        energy=float(data[panel+'/energies'][row])
        F=data[panel+'/amplitudes_at_one'][row].copy();C=data[panel+'/covariance'][row].copy()
        labels=metadata['panels'][panel]
    if horizon['source_panel']!=panel or horizon['source_row']!=row or energy.hex()!=horizon['energy_hex']:
        raise ValueError('source/horizon energy or row mismatch')
    if labels['mass']!=np.pi/2 or labels['angular_magnitude']!=np.sqrt(5) or labels['angular_sign']!=1:
        raise ValueError('original massive group14 channel required')
    traces={};cases=[]
    with ctx.workprec(192):
        mass,angular=ball(horizon['mass']),ball(horizon['angular'])
        frame=metric_horizon_frame(energy,mass,angular,config['horizon_rho'],bits=192)
        R=acb(ball(horizon['affine_reflection_direct_frame']['real']),
              ball(horizon['affine_reflection_direct_frame']['imag']))
        y_start=float(ball(horizon['compact_distance']).log().mid())
        y_target=(frame.q-3*arb.pi()/4).log();y_end=float(y_target.mid())
        initial=initial_bloch(frame,R,config['surface_gravity'],y_start)
        model=BlochSource(frame.q,energy,mass,angular,bits=192)
        for name,step in (('control',.05),('fine',.025)):
            trace,nfev=capture_bloch(model,initial,y_start,y_end,max_step=step)
            proof=validate_bloch(model,initial,trace,y_target,degree=12)
            error=original_covariance_error(F,C,proof)
            traces[name]=trace
            cases.append({'name':name,'max_step_hex':step.hex(),'nfev':nfev,'cells':proof['cells'],
                'initial_bloch_error_upper':exact_upper(proof['initial_error']),
                'normalized_defect_integral_upper':exact_upper(proof['defect_integral']),
                'endpoint_bridge_upper':exact_upper(proof['endpoint_bridge']),
                'bloch_error_upper':exact_upper(proof['bloch_error']),
                'rho1_covariance_operator_error_upper':exact_upper(error)})
        coherence_dropped=original_covariance_distance_bounds(F,np.diag(np.diag(C)),proof)['lower']
        if not coherence_dropped>arb('.1'):
            raise ArithmeticError('coherence mutation did not separate from the source')
        raw=deterministic_npz_bytes(traces)
        record={'schema':'NSC-SUBGAP-SOURCE-COVARIANCE-v1',
            'status':'ENCLOSED: one original subgap source covariance at rho=1; full preparation budget OPEN',
            'source_panel':panel,'source_row':row,'group':14,'energy_hex':energy.hex(),
            'angular_sign':1,'energy_sign':1,'target_rho_hex':float(1).hex(),
            'source_kappa_hex':float(config['surface_gravity']).hex(),
            'start_log_delta_hex':y_start.hex(),'end_log_delta_hex':y_end.hex(),
            'bits':192,'metric_terms':48,'defect_degree':12,'cases':cases,
            'trajectory_columns':['y0','y1','n_start[3]','n_end[3]','F1..F6[3]'],
            'coherence_dropped_operator_difference_lower':exact_upper(coherence_dropped),
            'source_columns_replaced':False,'source_covariance_replaced':False,
            'reference_preparation_only':True,'changed_history_state_evaluated':False,
            'energy_quadrature_error_included':False,'all_source_families':False,
            'numerical_continuation_to_rho_up_error':None,
            'physical_upstream_budget_component':None,'physical_local_gate':'OPEN',
            'payload':{'path':PAYLOAD,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)},
            'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
            'input_hashes':{p:digest(p) for p in (HORIZON,INVENTORY,payload['path'])}}
    return record,raw


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();record,raw=calculate()
    encoded=(json.dumps(record,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:
        if (ROOT/OUTPUT).exists() or (ROOT/PAYLOAD).exists():
            raise FileExistsError('source covariance record already exists; use --check')
        publish_exclusive_file(ROOT,PAYLOAD,raw);publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded or (ROOT/PAYLOAD).read_bytes()!=raw:
        raise ValueError('source covariance proof, witness or dependencies changed')
    print(json.dumps({'status':record['status'],'cases':[
        {'name':c['name'],'cells':c['cells'],
         'source_covariance_error_upper':float(restored_upper(c['rho1_covariance_operator_error_upper']))}
         for c in record['cases']],'physical_upstream_budget_component':None},indent=2))
