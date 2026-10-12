#!/usr/bin/env python3
"""Compose authenticated source refinements; no scientific producer rerun."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_source_update import refined_matter_source
from recursive_horizons.nsc_incoming_joint_constraints import joint_constraint_approximant,source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets,IncomingNormalJetChange
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_up_float

OUTPUT=ROOT/'results/development/nsc-incoming-source-update.json'
NAMES=('nsc-incoming-spectral-source','nsc-incoming-subgap-source','nsc-incoming-source-tail',
       'nsc-incoming-middle-order-correction','nsc-incoming-retained-order-correction',
       'nsc-incoming-low-high-source','nsc-incoming-centered-order24','nsc-incoming-retained-order24-bound',
       'nsc-incoming-joint-constraints','nsc-incoming-local-constraints','nsc-incoming-energy-error-budget',
       'nsc-incoming-vacuum-tail-bound','nsc-incoming-middle-bound')
SOURCES=('src/recursive_horizons/nsc_incoming_source_update.py','scripts/derive_nsc_incoming_source_update.py',
         'tests/test_nsc_incoming_source_update.py','docs/nsc-incoming-source-update.md')
OWNERS=('src/recursive_horizons/nsc_incoming_joint_constraints.py',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py','src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def authenticated_records():
    records={}
    for name in NAMES:
        record=json.loads((ROOT/f'results/development/{name}.json').read_text())
        for key in('source_hashes','input_hashes','numerical_owner_and_input_hashes','prepared_dependency_hashes'):
            for path,digest in record.get(key,{}).items():
                if sha(path)!=digest:raise ValueError('source-update dependency changed: '+path)
        payloads=list(record.get('input_payloads',[]))
        if 'payload' in record:payloads.append(record['payload'])
        for payload in payloads:
            if sha(payload['path'])!=payload['sha256']:raise ValueError('source-update artifact changed')
        records[name]=record
    return records


def make_record():
    r=authenticated_records()
    matter=refined_matter_source(*(r[n] for n in NAMES[:6]))
    b14,b12,low=r[NAMES[6]],r[NAMES[7]],r[NAMES[5]]
    for bound,group in ((b14,14),(b12,12)):
        if bound['radial']['group']!=group or bound['energy']['physical_Riccati_order']!=24:
            raise ValueError('matching physical order24 error certificates required')
        if bound['stationarity_tolerance']!=3e-11:raise ValueError('stationarity tolerance changed')
    if low['physical_mode_energy_bound']['physical_Riccati_order']!=24:
        raise ValueError('group22 physical bound must match new source order24')
    old=r['nsc-incoming-joint-constraints'];baseline_domain=incoming_cauchy_jets()
    reference=old['reference_controls']['n24']
    probe=incoming_cauchy_jets(tuple(IncomingNormalJetChange(row['field'],row['normal_order'],row['spatial_order'],row['delta']) for row in reference['changed_normal_entries']))
    zero={**reference,'changed_normal_entries':[],'reference_action_gradient_change':[0.,0.,0.,0.]}
    baseline=joint_constraint_approximant(baseline_domain,matter,r['nsc-incoming-local-constraints']['baseline'],zero)
    changed=joint_constraint_approximant(probe,matter,old['local_probe'],reference)
    residuals=dict(matter['composition_residuals'])
    increment=np.asarray(matter['action_gradient_increment'])[:2]
    for name,new in (('baseline',baseline),('mixed_normal_probe',changed)):
        residuals[name+'_increment_identity']=float(np.max(abs(np.asarray(new['action_gradient_approximant'])-np.asarray(old[name]['action_gradient_approximant'])-increment)))
    if max(residuals.values())>3e-13:raise ArithmeticError('source-to-constraint composition failed')
    with _precision(40):
        a=mp.iv.sqrt(3*mp.iv.pi/2-4)
        thermal14=mp.iv.mpf(r['nsc-incoming-middle-order-correction']['thermal_scattering_stress_error_upper'][0])
        new14=mp.iv.mpf(b14['energy']['lapse_action_error_upper'])+8*mp.iv.pi*a*thermal14
        new12=mp.iv.mpf(b12['conditional_combined_lapse_error_upper'])
        new22=mp.iv.mpf(low['physical_error_budget']['lapse_action_error_upper'])
        remaining=sum((mp.iv.mpf(row['lapse_action_error_upper']) for row in r['nsc-incoming-energy-error-budget']['all_group_order16_energy_bounds'] if row['group'] not in(12,14)),mp.iv.mpf(0))
        tail=r['nsc-incoming-vacuum-tail-bound'];vac=mp.iv.mpf(tail['aggregate']['vacuum_tail_error_upper'][0])
        thermal=mp.iv.mpf(tail['thermal_source_bound']['decimal_upper_bounds'][0])
        tailN=8*mp.iv.pi*a*(vac+thermal)
        partial={'middle14':_up_float(new14),'middle12':_up_float(new12),'low22_32_40':_up_float(new22),'all_group_high_tail':_up_float(tailN)}
        budget={'certified_refined_and_tail_lapse_error_upper':_up_float(new14+new12+new22+tailN),
                'remaining30_middle_lapse_error_upper':_up_float(remaining),
                'all_middle_plus_refined_low_cell_and_tail_lapse_error_upper':_up_float(remaining+new14+new12+new22+tailN),
                'full_source_error_bound':None,'rigorous_source_quadrature_error_included':False,
                'other_low_subgap_errors_bounded':False,'stationarity_tolerance':3e-11,
                'full_source_status':'OPEN'}
    paths=[f'results/development/{n}.json' for n in NAMES]
    return {'schema':'NSC-INCOMING-SOURCE-UPDATE-v1','accountable_author':'Douglas Ek',
            'status':'PASS: disjoint source refinements composed; physical constraints OPEN',
            'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in(*paths,*OWNERS)},
            'matter':matter,'baseline':baseline,'mixed_normal_probe':changed,
            'residuals':residuals,'composition_tolerance':3e-13,
            'certified_lapse_component_budgets':partial,'partial_error_budget':budget,
            'scope':{'old_source_records_overwritten':False,'local_light_geometry_count':1,
                     'same_physical_source_law':True,'physical_initial_data_selected':False,
                     'extended_stationarity':'OPEN','metric_evolution':False,'Gamma_rest_assigned':False,
                     'old_scientific_generators_rerun':False,'push_or_PDF':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_source_update.py --check'}


def encode(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.floating):return float(value)
    raise TypeError(type(value).__name__)


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write and OUTPUT.exists():raise FileExistsError('existing source update is not overwritten')
    result=json.loads(json.dumps(make_record(),default=encode))
    if args.write:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(json.loads(OUTPUT.read_text()),result)
    print(json.dumps({k:result[k] for k in('status','baseline','partial_error_budget','residuals')},indent=2))


if __name__=='__main__':main()
