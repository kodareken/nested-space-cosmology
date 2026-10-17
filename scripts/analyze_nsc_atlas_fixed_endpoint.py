#!/usr/bin/env python3
"""Compare three completed endpoint receipts without rerunning the evolution."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

from recursive_horizons.nsc_atlas_endpoint_comparison import (
    ENDPOINT_SOURCE_X, jordan_radius, refinement, scaled_radius,
)

ROOT=Path(__file__).resolve().parents[1]
GAPS={
    'mass_trace':'connection.physical_temporal_mass_actual',
    'mass_matched':'connection.physical_temporal_mass_matched',
    'wall_flux':'connection.conformal_wall_flux',
}


def compare(records):
    if len(records)!=3 or any(type(r.get('level')) is not int for r in records):
        raise ValueError('three numbered resolution receipts required')
    records=sorted(records,key=lambda r:r['level'])
    if [r['level'] for r in records]!=[0,1,2]:
        raise ValueError('require levels 0, 1 and 2 exactly once')
    for r in records:
        if r['schema']!='nsc-atlas-fixed-endpoint-v1' or r['status']!='completed_numerical_endpoint':
            raise ValueError('all three evolutions must complete their endpoint')
        if (r['execution']['last_accepted_x']!=r['geometry']['x_right'] or
            not math.isclose(r['geometry']['x_right'],ENDPOINT_SOURCE_X,rel_tol=0,abs_tol=2e-15) or
            r['geometry']['x_left']!=.4 or r['geometry']['reference_radius_endpoint']!=.55):
            raise ValueError('receipt does not cover the full P001 interval')
        if len(r['samples'])!=9:
            raise ValueError('nine common sample cuts required')
        if r['execution']['raw_history_saved'] or r['execution']['final_state_saved']:
            raise ValueError('endpoint comparison expects compact receipts')
        if r['endpoint']['frame_formula_needs_check']:
            raise ValueError('unresolved endpoint frame conversion')
        if r['radial_accuracy']!=dict(rtol=2.5e-14,atol=1e-15,adjoint_order=32,
                                    gradient_absolute_tolerance=2e-8):
            raise ValueError('endpoint did not use the confirmed radial accuracy')
    if ([r['geometry']['cap_cells'] for r in records]!=[2,4,8] or
        [r['geometry']['exterior_cells'] for r in records]!=[32,64,128]):
        raise ValueError('unexpected mesh family')
    for group,key in [('input','sha256'),('input','alpha'),('input','source_rows'),
                      ('preparation','input_sha256'),('preparation','original_input_sha256')]:
        if any(r[group][key]!=records[0][group][key] for r in records[1:]):
            raise ValueError('physical input differs: '+group+'.'+key)
    cuts=[s['x'] for s in records[0]['samples']]
    if any([s['x'] for s in r['samples']]!=cuts for r in records[1:]):
        raise ValueError('sample cuts differ between resolutions')
    maxima={name:refinement([r['maximum_absolute'][key] for r in records])
            for name,key in GAPS.items()}
    accumulated={name:refinement([r['transport_comparison'][key] for r in records])
        for name,key in [('mass_endpoint','mass_absolute_gap'),
                         ('wall_scalar_endpoint','scalar_absolute_gap'),
                         ('mass_running_maximum','max_abs_running_mass_gap'),
                         ('wall_scalar_running_maximum','max_abs_running_scalar_gap')]}
    radii=[]
    for r in records:
        endpoint=r['endpoint']
        radius=jordan_radius(endpoint['einstein_radius'],endpoint['c'],r['input']['alpha'])
        if not math.isclose(radius,endpoint['jordan_radius'],rel_tol=4e-15,abs_tol=0):
            raise ValueError('saved Jordan radius differs from invariant conversion')
        radii.append(radius)
    paired=[abs(b-a) for a,b in zip(radii,radii[1:])]
    radius_contracts=paired[1]<paired[0]
    passes=all(r['strictly_decreasing'] for r in maxima.values())
    return dict(schema='nsc-atlas-fixed-endpoint-comparison-v1',
        status='completed_saved_receipt_comparison',
        physical_evolutions=0, all_endpoints_completed=True,
        residual_maxima=maxima, accumulated_transport_gaps=accumulated,
        requested_mass_and_flux_contraction=passes,
        endpoint_radius_pair_differences=paired,
        endpoint_radius_pair_difference_contracts=radius_contracts,
        physical_observable=scaled_radius(radii),
        prediction_status=('finite_amplitude_numerical_estimate_with_resolution_envelope'
                           if passes and radius_contracts else 'unconverged_numerical_estimate'),
        sample_cuts=cuts,
        source_commits=[r['source_commit'] for r in records],
        shared_alpha=records[0]['input']['alpha'],
        radial_accuracy=records[0]['radial_accuracy'],
        numerical_method_note='Nonuniform fixed cap and exterior meshes, with simultaneous time-step refinement; no asymptotic order or rigorous error bound inferred.',
        operator_revision_note='All three successors use the confirmed tighter forward radial solve, the unchanged 32-point adjoint and gradient guard, and the same physical input. Earlier stopped or interrupted executions remain separate evidence.',
        independent_observational_test=False)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,nargs=3,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    inputs=[];records=[]
    for path in args.inputs:
        full=path.resolve()
        if not full.is_relative_to(ROOT):
            raise ValueError('receipt must belong to this repository')
        raw=full.read_bytes();records.append(json.loads(raw))
        inputs.append(dict(path=full.relative_to(ROOT).as_posix(),bytes=len(raw),
                           sha256=hashlib.sha256(raw).hexdigest()))
    result=compare(records)
    result['inputs']=inputs
    result['analysis_source_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    result['analysis_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ('status','requested_mass_and_flux_contraction','prediction_status')}))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
