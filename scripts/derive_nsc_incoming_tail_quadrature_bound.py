#!/usr/bin/env python3
"""One numerical-tail certificate against unchanged archived group22 values."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_tail_quadrature_bound import prepare_integral, replay_integral, directed_thermal_tail
from recursive_horizons.nsc_incoming_source_quadrature_bound import pack, unpack, float_enclosure
from recursive_horizons.nsc_incoming_middle_bound import _parameters
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-tail-quadrature-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_tail_quadrature_bound.py',
           'tests/test_nsc_incoming_tail_quadrature_bound.py',
           'scripts/derive_nsc_incoming_tail_quadrature_bound.py',
           'docs/nsc-incoming-tail-quadrature-bound.md')
INPUTS = ('results/development/nsc-incoming-source-tail.json',
          'results/development/nsc-incoming-vacuum-tail-bound.json',
          'results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-compact-matched-restart.json',
          'src/recursive_horizons/nsc_incoming_source_tail.py',
          'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_compact_ctp_neck.py')


def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())
def signature(): return {p:digest(p) for p in (*SOURCES,*INPUTS)}


def inputs():
    for path in INPUTS[:4]:
        record = read(path)
        for field in ('source_hashes','input_hashes'):
            for p,expected in record.get(field,{}).items():
                if digest(p) != expected: raise ValueError('tail input changed: '+p)
        if 'payload' in record and digest(record['payload']['path']) != record['payload']['sha256']:
            raise ValueError('tail payload changed')
    tail,vacuum,inventory,config = [read(p) for p in INPUTS[:4]]
    row = tail['groups']['22']; physical = vacuum['groups'][21]
    if row['lower'] != 160. or physical['group'] != 22 or physical['lower'] != 160.:
        raise ValueError('unchanged group22 endpoint160 required')
    return inventory['channels'][22],config['scattering_provenance']['config'],row,physical


def prepare():
    channel,config,_,_ = inputs(); before = signature()
    integral = prepare_integral(channel)
    thermal = directed_thermal_tail(channel,config)
    if signature() != before: raise ValueError('producer changed during pilot')
    payload = {'signature':before,'integral':integral,'thermal':thermal,
               'versions':{'mpmath':mp.__version__,'sympy':sp.__version__}}
    raw = (json.dumps(payload,sort_keys=True,indent=2)+'\n').encode(); h = sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-tail-quadrature-bound.{h}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    channel,config,archived,physical = inputs(); payload = json.loads(path.read_text())
    if payload['signature'] != signature() or payload['versions'] != {'mpmath':mp.__version__,'sympy':sp.__version__}:
        raise ValueError('pilot signature or numerical library changed')
    integral = replay_integral(payload['integral']); thermal = directed_thermal_tail(channel,config)
    if thermal != payload['thermal']: raise ValueError('outward thermal replay changed')
    old = archived['tail_by_order']['13']
    with _precision(80):
        intervals = [unpack(v) for v in integral['integral_interval']]
        error = [mp.iv.mpf(_hi(abs(value-mp.iv.mpf(stored)))) for value,stored in zip(intervals,old)]
        a,r,_,_,_ = _parameters(channel)
        numerical = [4*mp.iv.pi*a*r*r*error[0],4*mp.iv.pi*a*a*r*r*error[2]]
        physical_stress = [mp.iv.mpf(v) for v in physical['vacuum_tail_error_upper']]
        thermal_stress = [unpack(v) for v in thermal['intervals']]
        combined = [4*mp.iv.pi*a*r*r*(error[0]+physical_stress[0]+thermal_stress[0]),
                    4*mp.iv.pi*a*a*r*r*(error[2]+physical_stress[2]+thermal_stress[2])]
        zero_safe = lambda v: 0. if _hi(v) == 0 else _up_float(v)
        numerical_upper = [zero_safe(v) for v in numerical]
        combined_upper = [zero_safe(v) for v in combined]
        numeric_stress = [zero_safe(v) for v in error]
        display = [float_enclosure(v) for v in intervals]
    return {'schema':'NSC-INCOMING-TAIL-QUADRATURE-BOUND-v1','accountable_author':'Douglas Ek',
            'status':('PASS' if max(combined_upper)<=3e-11 else 'OPEN')+': group22 archived order16 infinite tail numerical and physical budget; full source OPEN',
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'payload':{'path':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest()},
            'group':22,'lower_energy':160.,'riccati_order':16,'reference_order':4,
            'kernel_order':['rho','p_parallel','T01','p_perp'],
            'archived_order13_source':old,'archived_source_unchanged':True,
            'certified_integral':integral,'integral_float_enclosure':display,
            'archived_numerical_stress_error_upper':numeric_stress,
            'archived_numerical_action_error_upper':numerical_upper,
            'imported_vacuum_physical_stress_error_upper':physical['vacuum_tail_error_upper'],
            'thermal_tail':thermal,'combined_action_error_upper':combined_upper,
            'stationarity_tolerance':3e-11,
            'endpoint_certificate':payload['integral']['endpoint_certificate'],
            'analytic_adapter':payload['integral']['analytic_adapter'],
            'residuals':{'exact_endpoint_identities':0.,'saved_gauss_weight_replay':0.,'thermal_replay':0.},
            'verification_tolerances':{'exact_endpoint_identities':0.,'saved_gauss_weight_replay':0.,'thermal_replay':0.},
            'scope':{'old_Taylor_regenerated':False,'radial_or_mode_solves':False,'source_recipe_changed':False,
                     'new_source_value_composed':False,'groups_calculated':1,'complete32_group_tail_certified':False,
                     'finite_offset_low_modes_certified':False,'constraints_solved':False,
                     'metric_evolution':False,'Gamma_rest_assigned':False,'push_or_PDF':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_tail_quadrature_bound.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true'); mode.add_argument('--check',action='store_true'); args=parser.parse_args()
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('pilot is not overwritten')
        result = replay(prepare()); OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        old = json.loads(OUTPUT.read_text())
        if digest(old['payload']['path']) != old['payload']['sha256']: raise ValueError('pilot artifact changed')
        result = replay(ROOT/old['payload']['path'])
        if result != old: raise ValueError('tail numerical record differs from replay')
    print(json.dumps({k:result[k] for k in ('status','archived_numerical_action_error_upper','combined_action_error_upper','residuals')},indent=2))


if __name__ == '__main__': main()
