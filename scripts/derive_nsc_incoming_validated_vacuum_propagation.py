#!/usr/bin/env python3
"""Two fixed E24 field propagations; independent resumable residual receipts."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile

import mpmath as mp
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_validated_vacuum_propagation import propagate_sign,replay_propagation
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision
from recursive_horizons.nsc_incoming_source_quadrature_bound import float_enclosure

OUTPUT=ROOT/'results/development/nsc-incoming-validated-vacuum-propagation.json'
SOURCES=('src/recursive_horizons/nsc_incoming_validated_vacuum_propagation.py',
         'scripts/derive_nsc_incoming_validated_vacuum_propagation.py',
         'tests/test_nsc_incoming_validated_vacuum_propagation.py',
         'docs/nsc-incoming-validated-vacuum-propagation.md')
INPUTS=('results/development/nsc-incoming-matched-horizon-initializer.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json',
        'src/recursive_horizons/nsc_incoming_matched_horizon_initializer.py',
        'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
        'src/recursive_horizons/nsc_incoming_middle_bound.py',
        'src/recursive_horizons/nsc_incoming_state_moments.py')


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signature():return{p:digest(p) for p in(*SOURCES,*INPUTS)}


def cache():
    key=sha256(json.dumps(signature(),sort_keys=True).encode()).hexdigest()
    path=Path(tempfile.gettempdir())/('nsc-validated-E24-'+key[:20]);path.mkdir(exist_ok=True)
    return path,key


def initial_inputs():
    initial=json.loads((ROOT/INPUTS[0]).read_text())
    for field in('source_hashes','input_hashes'):
        for p,h in initial[field].items():
            if digest(p)!=h:raise ValueError('matched initializer owner or input changed')
    if digest(initial['payload']['path'])!=initial['payload']['sha256']:raise ValueError('initializer artifact changed')
    channel=json.loads((ROOT/INPUTS[1]).read_text())['channels'][32]
    config=json.loads((ROOT/INPUTS[2]).read_text())['scattering_provenance']['config']
    return channel,config,initial


def read_sign(spec):
    if digest(spec['path'])!=spec['sha256']:raise ValueError('propagation artifact bytes changed')
    data=json.loads((ROOT/spec['path']).read_text())
    if data['signature']!=signature() or data['propagation']['binding']['sign']!=spec['sign']:
        raise ValueError('propagation input/producer/sign changed')
    checked=replay_propagation(data['propagation'])
    if checked!=data['checked']:raise ValueError('continuous propagation residual replay differs')
    return checked


def prepare_sign(sign):
    directory,key=cache();pointer=directory/f'sign{sign}.receipt.json'
    if pointer.exists():
        spec=json.loads(pointer.read_text());read_sign(spec)
        print(json.dumps({'resumed_completed_sign':sign}),flush=True)
        return spec
    channel,config,initial=initial_inputs();before=signature()
    print(json.dumps({'started_sign':sign,'checkpoint':str(directory/f'sign{sign}.checkpoint.json')}),flush=True)
    with threadpool_limits(limits=1):
        result=propagate_sign(channel,config,sign,initial,checkpoint=directory/f'sign{sign}.checkpoint.json',
                              checkpoint_key=key,progress=lambda r:print(json.dumps(r),flush=True))
        checked=replay_propagation(result)
    if signature()!=before:raise ValueError('propagation owner/input changed during fixed run')
    raw=json.dumps({'schema':'NSC-VALIDATED-VACUUM-SIGN-v1','signature':before,'propagation':result,'checked':checked},
                   sort_keys=True,allow_nan=False).encode();h=sha256(raw).hexdigest()
    relative=f'results/development/artifacts/nsc-incoming-validated-vacuum-sign{sign}.{h}.json'
    path=ROOT/relative
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-addressed collision')
    if not path.exists():path.write_bytes(raw)
    spec={'sign':sign,'path':relative,'sha256':h,'bytes':len(raw)}
    pointer.write_text(json.dumps(spec,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'completed_sign':sign,'checked':checked,'artifact':relative}),flush=True)
    return spec


def make_record(specs):
    if [s['sign'] for s in specs]!=[1,-1]:raise ValueError('both actual angular signs once required')
    signs={str(s['sign']):read_sign(s) for s in specs}
    with _precision(80):
        total=sum((mp.iv.mpf(v['pointwise_weighted_lapse_kernel_error_upper']) for v in signs.values()),mp.iv.mpf(0))
        upper=float_enclosure(total)[1]
    passed=upper<=1e-13
    return {'schema':'NSC-INCOMING-VALIDATED-VACUUM-PROPAGATION-v1','accountable_author':'Douglas Ek',
            'status':('PASS: E24 both-sign vacuum field propagation certified; spectral source OPEN' if passed else
                      'OPEN: fixed propagation pointwise lapse-kernel bound exceeds1e-13'),
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'sign_payloads':specs,'signs':signs,'group':32,'energy':24,
            'pointwise_weighted_lapse_kernel_error_upper':upper,'pilot_pointwise_tolerance':1e-13,
            'numerical_gate_pass':passed,
            'method':{'Taylor_degree':36,'maximum_real_step':'1/32','analytic_radius':'1/8',
                      'Cauchy_sample_radius':'1/16','directed_DFT_nodes':128,'interval_precision':80,
                      'sinc_series_degree':36,'separate_column_error_radius':True,'maximum_CPU_workers':2,
                      'exact_geometry_disks_nonzero':True,'continuous_residuals_enclosed':True,
                      'coefficient_step_centering_rounding_included':True},
            'scope':{'physical_state':'same affine limiting vacuum; certified new initializer and unitary fixed-profile field equation',
                     'Hamiltonian_mass_angular':'unchanged stored binary channel labels',
                     'thermal_source_kappa_changed':False,'archived_modes_normalized_or_changed':False,
                     'source_quadrature_or_window_integral_claimed':False,'full_source_accuracy_claimed':False,
                     'physical_IV_or_stationarity_root':False,'metric_timestep':False,'push_or_PDF':False},
            'next_connection':'extend validated field error to an energy panel with an independent spectral integration bound; one E24 result does not certify E24–32',
            'reproducer':'python3 scripts/derive_nsc_incoming_validated_vacuum_propagation.py --check'}


def main():
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true')
    p.add_argument('--workers',type=int,default=2);args=p.parse_args()
    if args.workers not in(1,2):raise ValueError('at most2 CPU workers required')
    if args.prepare:
        if OUTPUT.exists():raise FileExistsError('completed propagation receipt is not overwritten')
        print('resumable checkpoints: '+str(cache()[0]),flush=True)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:specs=list(pool.map(prepare_sign,(1,-1)))
        result=make_record(specs);OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        old=json.loads(OUTPUT.read_text());result=make_record(old['sign_payloads'])
        if old!=result:raise ValueError('propagation record replay differs')
    print(json.dumps({'status':result['status'],'pointwise_weighted_lapse_kernel_error_upper':result['pointwise_weighted_lapse_kernel_error_upper']},indent=2))


if __name__=='__main__':main()
