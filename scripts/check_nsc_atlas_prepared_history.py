#!/usr/bin/env python3
"""One physical prepared-history transfer; summaries and hashes only."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import plain
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget,solve_prepared_field
from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
from recursive_horizons.nsc_atlas_background_profile import build_prepared_background

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('scripts/check_nsc_atlas_prepared_history.py',
    *(f'src/recursive_horizons/nsc_atlas_{name}.py' for name in
      ('prepared_response','cap_field','background_profile','characteristic_background',
       'cap_background_jets','cap_mass_derivatives','field_preparation','field_transfer',
       'causal_reference','causal_observers')))


def array_hash(arrays):
    digest=hashlib.sha256()
    for name,value in sorted(arrays.items()):
        value=np.ascontiguousarray(value,dtype='<f8')
        digest.update(name.encode());digest.update(str(value.shape).encode());digest.update(value.tobytes())
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise ValueError('preserve the existing result; choose a new admitted output')
    if not args.output.parent.is_dir(): raise ValueError('output parent must already exist')
    subprocess.run(['git','diff','--exit-code','HEAD','--',*SOURCES],cwd=ROOT,check=True,capture_output=True)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    budget=PreparedFieldBudget(cpu_seconds=60.,kernel_values_max=10000000)
    report=dict(schema='nsc-atlas-prepared-history-v1',source_commit=commit,
        inputs=dict(field_nodes=257,geometry_knots=257,derivative_order=4,quiet_tolerances=[0.,0.,0.]),
        limits=dict(charged_values=10000000,cpu_seconds=60.),
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES})
    phase='existing prepared source solve';source_completed=False
    try:
        field=solve_prepared_field(nodes=257,geometry_knots=257,budget=budget)
        source_completed=True
        cap=CapLeadingField(field,derivative_order=4)
        report['physical_input_domain']=cap.domain
        phase='native history transfer with zero quiet tolerances'
        result=build_prepared_background(cap)
        background=result['background'];data=background.inputs
        report['transfer']={key:result[key] for key in
            ('onset','diagnostics','charged_evaluations','curvature_method','numerical_error_enclosure','scope')}
        report['input_sha256']=array_hash(data)
        phase='existing-interval midpoint comparison'
        mid=data['labels'][:-1]+np.diff(data['labels'])/2
        direct=cap.profile(mid,second=True);represented=background.profile.jet(mid,budget=budget)
        comparisons={}
        for original,retained in zip(('f','f_prime','f_second'),('F','F_prime','F_second')):
            expected=np.asarray(direct[original]);observed=np.asarray(represented[retained])
            delta=abs(expected-observed);index=int(np.argmax(delta));peak=float(np.max(abs(expected)))
            comparisons[retained]=dict(maximum_absolute=float(delta[index]),label=float(mid[index]),
                direct_peak=peak,peak_relative=float(np.max(delta)/peak) if peak else None)
        budget.check(6*len(mid))
        report.update(status='completed_interpolation_comparison',comparison=comparisons,
            comparison_points=len(mid),comparison_sha256=array_hash({**direct,**represented,'labels':mid}),
            maximum_inverse_bracket_width=cap.max_inverse_bracket_width,
            interpretation='finite comparisons on one prepared numerical field; no continuum error bound, finite-amplitude evolution or observational prediction')
    except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError) as error:
        report.update(status='stopped',failure=dict(phase=phase,reason=str(error)))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_integrations_completed=int(source_completed),source_integrations_attempted=1,
        retries=0,time_steps=0,raw_history_saved=False)
    args.output.write_text(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n',
                           encoding='utf-8',newline='\n')
    print(json.dumps(plain(dict(status=report['status'],execution=report['execution'],
        failure=report.get('failure'),comparison=report.get('comparison')))),flush=True)
    return 0 if report['status']=='completed_interpolation_comparison' else 1


if __name__=='__main__':raise SystemExit(main())
