#!/usr/bin/env python3
"""Record/replay finite-band source-phase transport with original inventory labels."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_massive_jost_transport_bound import (
    original_dense_solution, phase_transport_segments, subgap_phase_transport_bound)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap

OUTPUT = 'results/development/nsc-massive-jost-transport-bound-v1.json'
PAYLOAD = 'results/development/artifacts/nsc-massive-jost-transport-bound-v1.npz'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
OWNERS = (
    'scripts/derive_nsc_massive_jost_transport_bound.py',
    'src/recursive_horizons/nsc_massive_jost_transport_bound.py',
    'tests/test_nsc_massive_jost_transport_bound.py',
    'docs/nsc-massive-jost-transport-bound.md',
)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def calculate():
    source_record = json.loads((ROOT/INVENTORY).read_text())
    source_payload = source_record['payload']
    if digest(source_payload['path']) != source_payload['sha256']:
        raise ValueError('original source inventory payload changed')
    panel = 'group14/low16_1'
    with np.load(ROOT/source_payload['path'], allow_pickle=False) as archive:
        meta = json.loads(archive['metadata_json'].tobytes())
        energy = float(archive[panel+'/energies'][0])
        info, config = meta['panels'][panel], meta['config']
    mass, angular = float(info['mass']), float(info['angular_magnitude'])
    if info['group'] != 14 or info['angular_sign'] != 1 or mass != np.pi/2 or angular != np.sqrt(5.):
        raise ValueError('recorded group-14 pilot channel required')
    prep = PairedHorizonSeedMap(config['horizon_rho'], config['surface_gravity'],
        config['omega'], config['horizon_offset'], config['scattering_tolerance'], config['outer_floor'])
    solver = {'order':8, 'radial_collar':1e-10, 'rtol':2e-13, 'atol':2e-15}
    mode = solve_jost(prep.background,energy,mass,angular,**solver)
    near_radius = prep.horizon_rho+1.01e-4
    segments = phase_transport_segments(original_dense_solution(mode.run), prep.horizon_rho, near_radius)
    trace = np.array([[s.y_start,s.y_end,s.theta_start,s.theta_end,*s.corrections] for s in segments])
    raw = deterministic_npz_bytes({'phase_segments':trace})
    cases = []
    with ctx.workprec(192):
        m = arb(mass).union(arb.pi()/2)
        ell = arb(angular).union(arb(5).sqrt())
        for name, floor, degree, terms, subdivisions in (
            ('outer_control',8.,12,12,1),
            ('near_single_interval',near_radius,16,48,1),
            ('near_four_subintervals',near_radius,16,48,4),
        ):
            proof = subgap_phase_transport_bound(energy,m,ell,prep.background,
                mode=mode,inner_radius=floor,degree=degree,metric_terms=terms,
                defect_subdivisions=subdivisions)
            cases.append({'name':name,'requested_inner_radius_hex':float(floor).hex(),'proof':proof})
    result = {
        'schema':'NSC-MASSIVE-JOST-TRANSPORT-BOUND-v1',
        'status':'ENCLOSED: exterior finite-band phase transport; horizon sewing and physical source error OPEN',
        'source_panel':panel,'source_row':0,'family':[14,1],
        'energy_hex':energy.hex(),'mass_hex':mass.hex(),'angular_hex':angular.hex(),
        'horizon_coordinate_origin_hex':float(prep.horizon_rho).hex(),
        'solver_settings':solver, 'solver_outer_radius':float(mode.outer_radius),
        'captured_trace_columns':['y_start','y_end','theta_start','theta_end','F1','F2','F3','F4','F5','F6'],
        'cases':cases,
        'same_trajectory_for_all_enclosures':True,
        'source_inventory_unchanged':True,'physical_source_columns_replaced':False,
        'full_source_coverage':False,'horizon_sewing_error':None,
        'physical_rho1_source_error':None,'physical_local_gate':'OPEN',
        'payload':{'path':PAYLOAD,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)},
        'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
        'input_hashes':{p:digest(p) for p in (INVENTORY,source_payload['path'])},
    }
    return result,raw


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result,raw=calculate()
    if args.record:
        if (ROOT/OUTPUT).exists() or (ROOT/PAYLOAD).exists():
            raise FileExistsError('transport record exists; use --check')
        publish_exclusive_file(ROOT,PAYLOAD,raw)
        publish_exclusive_file(ROOT,OUTPUT,(json.dumps(result,sort_keys=True,indent=2)+'\n').encode())
    else:
        if json.loads((ROOT/OUTPUT).read_text())!=result or (ROOT/PAYLOAD).read_bytes()!=raw:
            raise ValueError('transport proof, trace or dependencies changed')
    print(json.dumps({'status':result['status'],'source_panel':result['source_panel'],
        'cases':[{'name':case['name'],'cells':case['proof']['cells'],
          'inner_offset_lower':float(restored_upper(case['proof']['inner_offset_certified_lower'])),
          'phase_error_upper':float(restored_upper(case['proof']['phase_error_inner_upper']))}
          for case in result['cases']], 'physical_rho1_source_error':None},indent=2))
