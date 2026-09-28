#!/usr/bin/env python3
"""One energy-polynomial step, no full-panel propagation or source integral."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_energy_panel import preflight,panel_identities,replay_step
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision

OUTPUT=ROOT/'results/development/nsc-incoming-energy-panel-preflight.json'
SOURCES=('src/recursive_horizons/nsc_incoming_energy_panel.py',
         'scripts/derive_nsc_incoming_energy_panel_preflight.py',
         'tests/test_nsc_incoming_energy_panel.py',
         'docs/nsc-incoming-energy-panel-preflight.md')
INPUTS=('results/development/nsc-incoming-matched-horizon-initializer.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json',
        'results/development/nsc-incoming-validated-vacuum-propagation.json',
        'src/recursive_horizons/nsc_incoming_matched_horizon_initializer.py',
        'src/recursive_horizons/nsc_incoming_validated_vacuum_propagation.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
        'src/recursive_horizons/nsc_incoming_middle_bound.py',
        'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py')


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signature():return{p:digest(p) for p in(*SOURCES,*INPUTS)}


def prepare():
    before=signature();matched=json.loads((ROOT/INPUTS[0]).read_text())
    for field in('source_hashes','input_hashes'):
        for p,h in matched[field].items():
            if digest(p)!=h:raise ValueError('matched initializer input/owner changed')
    if digest(matched['payload']['path'])!=matched['payload']['sha256']:raise ValueError('initializer artifact changed')
    c=json.loads((ROOT/INPUTS[1]).read_text())['channels'][32]
    cfg=json.loads((ROOT/INPUTS[2]).read_text())['scattering_provenance']['config']
    result=preflight(c,cfg,matched)
    if signature()!=before:raise ValueError('preflight input/producer changed')
    raw=json.dumps({'schema':'NSC-ENERGY-PANEL-PREFLIGHT-ARTIFACT-v1','signature':before,'preflight':result},sort_keys=True,allow_nan=False).encode()
    h=sha256(raw).hexdigest();path=ROOT/f'results/development/artifacts/nsc-incoming-energy-panel-preflight.{h}.json'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-addressed collision')
    if not path.exists():path.write_bytes(raw)
    return path


def make_record(path):
    data=json.loads(path.read_text())
    if data['signature']!=signature():raise ValueError('energy preflight producer/input changed')
    result=data['preflight']
    if result['identities']!=panel_identities():raise ValueError('energy identities changed')
    checks={}
    with _precision(80):
        for sign,case in result['cases'].items():
            replay=replay_step(case['step'])
            if any(case['step'][k]!=v for k,v in replay.items()):raise ValueError('uniform energy residual replay changed')
            checks[sign]={'initial_interpolation_error_upper':case['initial_interpolant']['column_interpolation_error_upper'],
                          'initial_DCT_rounding_upper':case['initial_interpolant']['DCT_coefficient_radius_l1_upper'],
                          'inherited_projector_error_separate':case['inherited_mathematical_projector_error_upper'],
                          'step_elapsed_seconds':case['step']['elapsed_seconds'],**replay}
    passed=all(row['local_gate_pass'] for row in checks.values())
    return {'schema':'NSC-INCOMING-ENERGY-PANEL-PREFLIGHT-v1','accountable_author':'Douglas Ek',
            'status':('PASS: initial Cheb24 and one-step energy residual; full panel not run' if passed else
                      'OPEN: fixed Cheb24 one-step residual exceeds1e-20; full panel not run'),
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'payload':{'path':str(path.relative_to(ROOT)),'sha256':digest(path)},
            'identities':result['identities'],'signs':checks,
            'shared_geometry_seconds':result['shared_geometry']['elapsed_seconds'],
            'elapsed_seconds':result['elapsed_seconds'],'future_work':result['future_work'],
            'initial_interpolation_source':'https://raw.githubusercontent.com/chebfun/ATAP/development/chap8.m (Theorem8.2 Eq8.3)',
            'scope':{'full_panel_run':False,'physical_source_integral':False,
                     'one_time_step_each_sign':True,'shared_geometry_preparations':1,
                     'physical_state_scales_or_source_kappa_changed':False,
                     'endpoint_column_renormalized':False,'archived_modes_modified':False,
                     'metric_or_physical_IV':False,'push_or_PDF':False},
            'next_decision':'review full residual and measured cost before authorizing any whole-panel run',
            'reproducer':'python3 scripts/derive_nsc_incoming_energy_panel_preflight.py --check'}


def main():
    p=argparse.ArgumentParser();m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=p.parse_args()
    if args.write:
        if OUTPUT.exists():raise FileExistsError('energy preflight record is not overwritten')
        result=make_record(prepare());OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        old=json.loads(OUTPUT.read_text());path=ROOT/old['payload']['path']
        if digest(path)!=old['payload']['sha256']:raise ValueError('energy preflight artifact changed')
        result=make_record(path)
        if result!=old:raise ValueError('energy preflight replay differs')
    print(json.dumps({'status':result['status'],'signs':result['signs'],'elapsed_seconds':result['elapsed_seconds'],
                      'future_work':result['future_work']},indent=2))


if __name__=='__main__':main()
