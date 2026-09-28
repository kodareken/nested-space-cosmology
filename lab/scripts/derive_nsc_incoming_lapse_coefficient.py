#!/usr/bin/env python3
"""Closed lapse-response coefficient; no spectral or history generator."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_lapse_coefficient import reference_identity,closed_lapse_coefficient,local_density_stencil

OUTPUT=ROOT/'results/development/nsc-incoming-lapse-coefficient.json'
SOURCES=('src/recursive_horizons/nsc_incoming_lapse_coefficient.py','scripts/derive_nsc_incoming_lapse_coefficient.py',
         'tests/test_nsc_incoming_lapse_coefficient.py','docs/nsc-incoming-lapse-coefficient.md')
INPUTS=('results/development/nsc-incoming-constraint-germ.json','results/development/nsc-incoming-local-constraints.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py','src/recursive_horizons/nsc_incoming_local_constraints.py',
        'src/recursive_horizons/nsc_light_restoration_action.py','src/recursive_horizons/nsc_magnetic_light_reference.py',
        'src/recursive_horizons/nsc_spherical_local_history.py','src/recursive_horizons/nsc_spatial_reference_symbol.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def make_record():
    records=[json.loads((ROOT/p).read_text()) for p in INPUTS[:3]]
    for record in records:
        for key in('source_hashes','input_hashes'):
            for p,h in record.get(key,{}).items():
                if sha(p)!=h:raise ValueError('coefficient input changed: '+p)
        if 'payload' in record and sha(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('coefficient control artifact changed')
    germ,local,inventory=records
    result=closed_lapse_coefficient(inventory['channels'],local['locked_inputs'])
    coarse=local_density_stencil(local['locked_inputs'],1.)
    fine=local_density_stencil(local['locked_inputs'],.5)
    distance=lambda value,bounds:max(bounds[0]-value,value-bounds[1],0.)
    residuals={'density_stencil_steps':max(abs(coarse[k]-fine[k]) for k in coarse),
               'local_closed_vs_density':max(distance(fine[k],bounds) for k,bounds in result['local_coefficient_intervals'].items()),
               'analytic_zero_local_channels':max(abs(fine[k]) for k in result['zero_channels'])}
    control_payload=json.loads((ROOT/germ['payload']['path']).read_text())['controls']['u']
    u=control_payload['probe'][0]['delta']
    reference=control_payload['references']['24']['reference_action_gradient_change'][0]/u
    residuals['reference_closed_vs_previous_integral']=distance(reference,result['reference_coefficient_interval'])
    old=germ['numerical_coefficients']['c_u']
    old_distance=distance(old,result['total_coefficient_interval'])
    residuals['old_probe_unscaled_action_difference']=old_distance*abs(u)
    if max(residuals.values())>3e-12 or not result['strictly_positive']:
        raise ArithmeticError('closed coefficient control or positivity failed')
    return {'schema':'NSC-INCOMING-CLOSED-LAPSE-COEFFICIENT-v1','accountable_author':'Douglas Ek',
            'status':'PASS: closed positive lapse-response coefficient; physical root OPEN',
            'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
            'identity':reference_identity(),'coefficient':result,
            'independent_density_stencils':{'step1':coarse,'step_half':fine},
            'previous_contour_coefficient':old,'old_coefficient_distance_to_interval':old_distance,
            'old_probe_interpretation':'dividing the small action response by0.01 amplifies its arithmetic/contour indicator; no new physical source correction',
            'residuals':residuals,'numerical_control_tolerance':3e-12,
            'scope':{'coefficient_positive_enclosure':True,'source_accuracy_required_for_this_coefficient':False,
                     'new_physical_parameter_or_term':False,'other_response_coefficients_exactly_enclosed':False,
                     'physical_initial_data_selected':False,'constraints_solved':False,
                     'metric_evolution':False,'Gamma_rest_assigned':False,'push_or_PDF':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_lapse_coefficient.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write and OUTPUT.exists():raise FileExistsError('existing coefficient record is not overwritten')
    result=make_record()
    if args.write:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(json.loads(OUTPUT.read_text()),result)
    print(json.dumps({'status':result['status'],'interval':result['coefficient']['total_coefficient_interval'],
                      'residuals':result['residuals']},indent=2))


if __name__=='__main__':main()
