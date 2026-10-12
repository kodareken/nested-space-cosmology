#!/usr/bin/env python3
"""Compose adjacent incoming source windows without rerunning producers."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_up_float

OUTPUT=ROOT/'results/development/nsc-incoming-source-update-v2.json'
NAMES=('nsc-incoming-source-update','nsc-incoming-window-reuse','nsc-incoming-group32-source','nsc-incoming-energy-error-budget')
SOURCES=('scripts/derive_nsc_incoming_source_update_v2.py','tests/test_nsc_incoming_source_update_v2.py','docs/nsc-incoming-source-update-v2.md')
OWNERS=('src/recursive_horizons/nsc_incoming_joint_constraints.py','src/recursive_horizons/nsc_incoming_cauchy_jets.py','src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def validate_windows(windows):
    expected={'group12_low':(12,'low_vacuum',[32.,40.]),'group22_middle':(22,'middle_correction',[40.,160.]),
              'group32_low':(32,'low_vacuum',[32.,40.]),'group32_middle':(32,'middle_correction',[40.,320.])}
    if set(windows)!=set(expected):raise ValueError('exactly the four explicit windows required')
    # Previously applied windows are included in the overlap check too.
    coverage={12:[[40.,320.]],14:[[16.,160.]],22:[[32.,40.]],32:[]}
    for name,row in windows.items():
        if (row['group'],row['kind'],row['interval'])!=expected[name] or row['source_order']!=24:
            raise ValueError('source window group, interval, kind or numerical order changed')
        if row['interval_binding']['energy']['physical_Riccati_order']!=24:
            raise ValueError('matching order24 physical certificate required')
        coverage[row['group']].append(row['interval'])
    for group,intervals in coverage.items():
        ordered=sorted(intervals)
        if any(a[1]>b[0] for a,b in zip(ordered,ordered[1:])):
            raise ValueError('source windows overlap in positive measure')
    return coverage


def make_record():
    records={}
    for name in NAMES:
        d=json.loads((ROOT/f'results/development/{name}.json').read_text())
        for key in('source_hashes','input_hashes'):
            for p,h in d[key].items():
                if sha(p)!=h:raise ValueError('source-update input changed: '+p)
        if 'payload' in d and sha(d['payload']['path'])!=d['payload']['sha256']:raise ValueError('source artifact changed')
        records[name]=d
    previous=records[NAMES[0]];windows={**records[NAMES[1]]['windows'],**records[NAMES[2]]['windows']}
    coverage=validate_windows(windows);residuals={};increments={}
    for name,row in windows.items():
        delta=np.asarray(row['explicit_source_delta']);increments[name]=delta
        residuals[name+'_difference']=float(np.max(abs(np.asarray(row['updated_source_approximant'])-row['old_source']-delta)))
    addition=sum(increments.values());action=source_action_gradient(incoming_cauchy_jets(),addition)
    independent=sum(np.asarray(row['additive_action_gradient']) for row in windows.values())
    residuals['source_to_action']=float(np.max(abs(action-independent)))
    matter={**previous['matter'],'previous_update_stress_approximant':previous['matter']['stress_approximant'],
            'stress_approximant':(np.asarray(previous['matter']['stress_approximant'])+addition).tolist(),
            'additional_refinement_parts':{k:v.tolist() for k,v in increments.items()},
            'addition_since_previous':addition.tolist(),
            'total_stress_increment':(np.asarray(previous['matter']['total_stress_increment'])+addition).tolist(),
            'action_gradient_increment':(np.asarray(previous['matter']['action_gradient_increment'])+action).tolist(),
            'full_source_error_bound':None}
    cases={}
    for name in('baseline','mixed_normal_probe'):
        old=previous[name];new=np.asarray(old['action_gradient_approximant'])+action[:2]
        case={**old,'action_gradient_approximant':new.tolist(),'force_approximant':(-new).tolist(),
              'matter_baseline_action_gradient':(np.asarray(old['matter_baseline_action_gradient'])+action[:2]).tolist(),
              'maximum_absolute_approximant':float(np.max(abs(new)))}
        residuals[name+'_increment']=float(np.max(abs(new-np.asarray(old['action_gradient_approximant'])-action[:2])))
        cases[name]=case
    if max(residuals.values())>3e-13:raise ArithmeticError('source-window composition failed')
    with _precision(40):
        certified=mp.iv.mpf(previous['partial_error_budget']['certified_refined_and_tail_lapse_error_upper'])
        for row in windows.values():certified+=mp.iv.mpf(row['physical_error_budget']['lapse_action_error_upper'])
        remaining=sum((mp.iv.mpf(row['lapse_action_error_upper']) for row in records[NAMES[3]]['all_group_order16_energy_bounds'] if row['group'] not in(12,14,22,32)),mp.iv.mpf(0))
        budget={'certified_refined_and_tail_lapse_error_upper':_up_float(certified),
                'remaining28_middle_lapse_error_upper':_up_float(remaining),
                'covered_spectral_domains_lapse_error_upper':_up_float(certified+remaining),
                'full_source_error_bound':None,'rigorous_source_quadrature_error_included':False,
                'other_low_subgap_errors_bounded':False,'stationarity_tolerance':3e-11,'full_source_status':'OPEN'}
    paths=[f'results/development/{n}.json' for n in NAMES]
    return {'schema':'NSC-INCOMING-SOURCE-UPDATE-v2','accountable_author':'Douglas Ek',
            'status':'PASS: four further disjoint source updates composed; physical constraints OPEN',
            'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in(*paths,*OWNERS)},
            'previous_record':paths[0],'new_source_windows':list(windows),'refined_coverage':coverage,
            'matter':matter,**cases,'additional_action_gradient':action.tolist(),'residuals':residuals,
            'composition_tolerance':3e-13,'partial_error_budget':budget,
            'scope':{'same_physical_source_law':True,'old_records_overwritten':False,'union_counted_as_extra_source':False,
                     'local_light_geometry_count':1,'new_source_or_radial_generators_run':False,
                     'physical_initial_data_selected':False,'metric_evolution':False,'extended_stationarity':'OPEN',
                     'Gamma_rest_assigned':False,'push_or_PDF':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_source_update_v2.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write and OUTPUT.exists():raise FileExistsError('existing source update is not overwritten')
    data=make_record()
    # JSON roundtrip normalizes integer dictionary keys for exact replay.
    data=json.loads(json.dumps(data))
    if args.write:OUTPUT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(json.loads(OUTPUT.read_text()),data)
    print(json.dumps({k:data[k] for k in('status','additional_action_gradient','baseline','partial_error_budget','residuals')},indent=2))


if __name__=='__main__':main()
