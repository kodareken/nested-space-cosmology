#!/usr/bin/env python3
"""Certify the exact initial reference against the saved onset event."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_reference_initial_bound import reference_initial_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={
 'source':('results/summaries/nsc-atlas-source-cell-bound-v2.json',
           '707bb8d76359edaf2a26228657a219933ebd164edcf89eb87742f9e66bec2bdd'),
 'gauge':('results/summaries/nsc-atlas-reference-gauge-bound-v1.json',
          'e8a8400cec0f6ae0bd831723696ebf1d41c38bdf0b052301ffab837465f1d216'),
 'onset':('results/summaries/nsc-atlas-fixed-enriched-history-onset-v1.json',
          'f6b5995960483d9285a4226e2c7ae1be4e6bd08ad53fdef2a30855edae7be62f')}
OWNERS=('scripts/check_nsc_atlas_reference_initial_bound.py',
        'src/recursive_horizons/nsc_atlas_reference_initial_bound.py',
        'src/recursive_horizons/nsc_atlas_fixed_cap_stage.py',
        'src/recursive_horizons/nsc_atlas_characteristic_background.py',
        'src/recursive_horizons/nsc_atlas_background_profile.py',
        'src/recursive_horizons/nsc_source_binding_graph.py','scripts/run_python.py','pyproject.toml')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new summary path')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for name,expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('freeze changed producer before execution: '+name)
    saved={}
    for key,(name,expected) in INPUTS.items():
        raw=(ROOT/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=expected: raise ValueError('input identity changed: '+name)
        saved[key]=json.loads(raw)
    gauge,onset=saved['gauge'],saved['onset'];prep=onset['preparation']
    if (gauge['status']!='certified_exact_reference_gauge_cell'
            or onset['status']!='completed_numerical_interval'
            or gauge['cell']['left']!=onset['geometry']['x_left']
            or prep['inputs']['quiet_tolerances']!=[0.,0.,0.]
            or prep['diagnostics']['quiet_jet_discrepancies']!=[0.,0.,0.]
            or not prep['enrichment']['onset_piece_unchanged']):
        raise ValueError('initial cell or quiet-prefix evidence differs')
    source=saved['source']
    if (source['cell']['left']!=gauge['cell']['left']
            or any(Q(v)!=0 for v in source['cell']['normalized_exact_binary_powers'][:2])):
        raise ValueError('source cell does not have the original quiet affine initial jets')
    base_name=source['input']['path'];base_raw=(ROOT/base_name).read_bytes()
    if hashlib.sha256(base_raw).hexdigest()!=source['input']['sha256']:
        raise ValueError('original source geometry input differs')
    source_geometry=json.loads(base_raw)['geometry']
    budget=Budget(maximum_evaluations=100000,cpu_seconds=20.)
    report=dict(schema='nsc-atlas-reference-initial-bound-v1',source_commit=commit,
                source_bindings=bindings,inputs={k:dict(path=v[0],sha256=v[1]) for k,v in INPUTS.items()},
                original_geometry_input=dict(path=base_name,sha256=source['input']['sha256']),
                retained_profile_sha256=prep['input_sha256'],status='incomplete')
    try:
        certificate=reference_initial_bound(gauge['certificate'],
            seed=Q(onset['geometry']['reference_radius'][0]),
            retained_wall={k:Q(v) for k,v in onset['initial']['wall'].items()},
            history_offset=Q(prep['clock_offset']),history_onset=Q(prep['onset']['label']),
            cap_knots=tuple(Q(v) for v in onset['geometry']['exact_cap_knots']),
            source_geometry=source_geometry,budget=budget)
        report.update(status='certified_quiet_exact_reference_initial_event',certificate=certificate,
                      retained_initial_source_constraint=onset['initial']['source_constraint'],
                      reconstruction='exact reference plus original zero departure coefficients; future continuous residual not yet certified')
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,time_steps=0,root_searches=0,retries=0,raw_history_saved=False)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','execution','failure')}))
    return 0 if report['status']=='certified_quiet_exact_reference_initial_event' else 1


if __name__=='__main__':
    raise SystemExit(main())
