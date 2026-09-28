#!/usr/bin/env python3
"""Record/replay the full-profile horizon frame and one reflection enclosure."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_metric_horizon_frame import metric_horizon_frame, reflection_from_phase

OUTPUT = 'results/development/nsc-metric-horizon-frame-v1.json'
INPUT = 'results/development/nsc-massive-jost-transport-bound-v1.json'
OWNERS = ('scripts/derive_nsc_metric_horizon_frame.py',
          'src/recursive_horizons/nsc_metric_horizon_frame.py',
          'tests/test_nsc_metric_horizon_frame.py',
          'docs/nsc-metric-horizon-frame.md')


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def interval(value):
    return {'lower':exact_upper(value.lower()), 'upper':exact_upper(value.upper())}


def complex_interval(value):
    return {'real':interval(value.real), 'imag':interval(value.imag)}


def calculate():
    source = json.loads((ROOT/INPUT).read_text())
    for group in ('input_hashes','source_hashes'):
        for path, expected in source[group].items():
            if digest(path) != expected:
                raise ValueError('source phase dependency changed: '+path)
    payload = source['payload']
    if digest(payload['path']) != payload['sha256']:
        raise ValueError('phase trace payload changed')
    if source['source_panel'] != 'group14/low16_1' or source['source_row'] != 0:
        raise ValueError('representative original low-energy row required')
    with np.load(ROOT/payload['path'],allow_pickle=False) as archive:
        trace = archive['phase_segments']
    case = next(c for c in source['cases'] if c['name']=='near_four_subintervals')
    if len(trace) != case['proof']['cells'] or trace.shape[1] != 10 or not np.isfinite(trace).all():
        raise ValueError('phase trace coverage mismatch')
    with ctx.workprec(192):
        E = float.fromhex(source['energy_hex'])
        m = arb(float.fromhex(source['mass_hex'])).union(arb.pi()/2)
        ell = arb(float.fromhex(source['angular_hex'])).union(arb(5).sqrt())
        h = float.fromhex(source['horizon_coordinate_origin_hex'])
        rho = arb(h)+arb(float(trace[-1,1])).exp()
        error = restored_upper(case['proof']['phase_error_inner_upper'])
        theta = arb(float(trace[-1,3]))+arb(0,error)
        frame = metric_horizon_frame(E,m,ell,h,order=16,bits=192)
        reflection, distance, tail = reflection_from_phase(frame,rho,theta)
        exterior, _ = frame.evaluate(distance)
        interior, _ = frame.evaluate(distance,interior=True)
        record = {
            'schema':'NSC-METRIC-HORIZON-FRAME-v1',
            'status':'ENCLOSED: exact-profile horizon series and representative affine reflection',
            'source_panel':source['source_panel'],'source_row':source['source_row'],
            'energy_hex':source['energy_hex'], 'mass':interval(m),'angular':interval(ell),
            'order':16,'bits':192,'analytic_t_radius':'0.1',
            'exact_horizon_q':interval(frame.q),'exact_metric_slope':interval(frame.slope),
            'analytic_generator_majorant_upper':exact_upper(frame.generator_majorant),
            'radial_coordinate':interval(rho),'compact_distance':interval(distance),
            'phase':interval(theta),'phase_error_upper':exact_upper(error),
            'frame_column_l1_tail_upper':exact_upper(tail),
            'exterior_frame':[[complex_interval(v) for v in row] for row in exterior],
            'interior_frame_same_compact_distance':[[complex_interval(v) for v in row] for row in interior],
            'affine_reflection_direct_frame':complex_interval(reflection),
            'reflection_modulus':interval(abs(reflection)),
            'source_inventory_unchanged':True,'physical_source_columns_replaced':False,
            'source_phase_convention_matched':False,'interior_transport_error':None,
            'physical_rho1_source_error':None,'full_source_coverage':False,
            'physical_local_gate':'OPEN',
            'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
            'input_hashes':{INPUT:digest(INPUT),payload['path']:digest(payload['path'])},
        }
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true')
    mode.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.record:
        publish_exclusive_file(ROOT,OUTPUT,(json.dumps(result,sort_keys=True,indent=2)+'\n').encode())
    elif json.loads((ROOT/OUTPUT).read_text()) != result:
        raise ValueError('horizon frame proof or dependency changed')
    print(json.dumps({'status':result['status'],
        'frame_tail_upper':float(restored_upper(result['frame_column_l1_tail_upper'])),
        'physical_rho1_source_error':None,'physical_local_gate':'OPEN'},indent=2))
