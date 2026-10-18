"""Locate the active retained history in the three saved P001 cap spaces.

Read-only diagnosis: no field evaluation, mass integration or evolution.
The incoming-front support belongs to the retained background, not to the
unarchived evolved departures. The all-stage maximizing states are unknown.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import time

from diagnose_nsc_atlas_endpoint_gaps import authenticate

ROOT=Path(__file__).resolve().parents[1]
CPU_LIMIT=60.0
BYTE_LIMIT=16000


def support_at_cut(knots, *, radius, central_time, clock_offset, onset_label):
    """F(T+z)-F(T-z) before the centre receives the quiet-onset front.

    The original profile is identically zero at labels <= onset_label.
    With T=central_time+clock_offset < onset, F(T-z)=0 for every z>=0.
    This locates the possible incoming support; it does not infer its size.
    """
    knots=tuple(Fraction(v) for v in knots)
    radius,central_time,clock_offset,onset_label=map(float,
        (radius,central_time,clock_offset,onset_label))
    if (len(knots)<2 or knots[0]!=0 or any(b<=a for a,b in zip(knots,knots[1:]))
            or not all(math.isfinite(v) for v in (radius,central_time,clock_offset,onset_label))
            or not 0<radius<1):
        raise ValueError('finite cap cut and increasing knots from zero required')
    wall=math.atanh(radius)
    if not float(knots[-2])<wall<float(knots[-1]):
        raise ValueError('wall must lie in the last fixed cell')
    history_time=central_time+clock_offset
    front=onset_label-history_time
    if front<=0:
        raise ValueError('this incoming-only diagnosis ends at central arrival')
    active=[i for i,(a,b) in enumerate(zip(knots,knots[1:]))
            if min(float(b),wall)>max(float(a),front)]
    return dict(history_time=history_time,incoming_front_z=front,wall_z=wall,
        last_cell_left=float(knots[-2]),last_cell_full_width=float(knots[-1]-knots[-2]),
        last_cell_clipped_width=wall-float(knots[-2]),
        incoming_supported_width=max(0.,wall-front),
        incoming_supported_cell_indices=active,all_inner_cells_before_onset=front>float(knots[-2]),
        front_to_last_inner_knot=front-float(knots[-2]))


def relative_change(left,right):
    if left==0: return None
    return (right-left)/abs(left)


def build_report(records,binding,*,cpu_seconds=0.):
    rows=[]
    for record in records:
        knots=record['geometry']['exact_cap_knots']
        preparation=record['preparation']
        kwargs=dict(clock_offset=preparation['clock_offset'],onset_label=preparation['onset']['label'])
        cuts=[]
        for sample in record['samples']:
            item=support_at_cut(knots,radius=sample['einstein_radius'],
                central_time=sample['central_proper_clock_T'],**kwargs)
            cuts.append(dict(x=sample['x'],**item))
        end=record['endpoint']
        rows.append(dict(level=record['level'],cap_cells=record['geometry']['cap_cells'],
            exterior_cells=record['geometry']['exterior_cells'],polynomial_degree=5,
            endpoint_support=support_at_cut(knots,radius=end['einstein_radius'],
                central_time=end['central_proper_clock_T'],**kwargs),
            samples_with_supported_field=sum(bool(c['incoming_supported_cell_indices']) for c in cuts),
            samples_with_support_only_in_last_cell=sum(c['incoming_supported_cell_indices']==[len(knots)-2] for c in cuts),
            all_sampled_inner_cells_before_onset=all(c['all_inner_cells_before_onset'] for c in cuts),
            common_cuts=cuts,
            endpoint_projection=end['background_projection'],
            endpoint_gap=end['mass_gap_decomposition'],
            endpoint_physical_mass_gap=end['temporal_residual']['mass_actual'],
            actual_mass_envelope=record['maximum_absolute']['temporal.physical_temporal_mass_actual'],
            integrated_absolute_mass_gap=record['integrated_defect']['absolute']['physical_temporal_mass_actual'],
            integrated_signed_mass_gap=record['integrated_defect']['signed']['physical_temporal_mass_actual']))
    pairs=[]
    for left,right,lr,rr in zip(rows,rows[1:],records,records[1:]):
        lk={Fraction(v) for v in lr['geometry']['exact_cap_knots'][:-1]}
        rk={Fraction(v) for v in rr['geometry']['exact_cap_knots'][:-1]}
        pairs.append(dict(levels=[left['level'],right['level']],old_internal_knots_preserved=lk<=rk,
            clipped_outer_width_ratio=right['endpoint_support']['last_cell_clipped_width']/left['endpoint_support']['last_cell_clipped_width'],
            actual_mass_envelope_relative_change=relative_change(left['actual_mass_envelope'],right['actual_mass_envelope']),
            absolute_integral_relative_change=relative_change(left['integrated_absolute_mass_gap'],right['integrated_absolute_mass_gap'])))
    return dict(schema='nsc-atlas-endpoint-active-resolution-v1',comparison=binding,levels=rows,pairs=pairs,
        scope='Retained-background support and saved residuals in the original spaces; no new grid or evolution.',
        diagnosis=[
            'The active retained incoming layer occupies one unbisected quintic outer cell at every sampled active cut on all three levels.',
            'Doubling the cap-cell count primarily adds cells below that retained front. The outer widths are .80/.73/.66, not a halved-h family.',
            'The spaces change their internal knots and are not nested refinements. Strong wall traces and background time projections need not improve monotonically.',
            'The absolute integrated mass defect also grows. A more densely sampled maximum alone cannot explain that trend.',
            'The known physical residual is alpha times the strong/weak trace difference plus the explicit affine-background remainder. Neither is a radial mass/lapse integration error.',
        ],
        not_established=[
            'The evolved departure coefficients and maximizing stages were not archived. Background support is not proof of departure support.',
            'This identifies an unresolved approximation and invalidates a cap-convergence inference; it does not uniquely attribute the all-stage maximum or exclude another coupling error.',
            'No Richardson order, continuum limit, certified uncertainty or completed Step 2 follows from these three runs.',
        ],
        next_calculation='Check the wall trace and affine-background transport on the same physical active-region space before redesigning any evolution grid.',
        runtime_cpu_seconds=cpu_seconds,counts=dict(evolutions=0,new_grids=0,mass_integrations=0,field_evaluations=0),
        limits=dict(CPU_seconds=CPU_LIMIT,summary_bytes=BYTE_LIMIT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison',default='results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    if not args.output.parent.is_dir(): raise ValueError('output parent must exist')
    started=time.process_time()
    binding,_,records,_=authenticate(ROOT,args.comparison)
    result=build_report(records,binding)
    result['runtime_cpu_seconds']=time.process_time()-started
    if result['runtime_cpu_seconds']>CPU_LIMIT: raise RuntimeError('diagnosis CPU budget exhausted')
    result['producer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>BYTE_LIMIT: raise RuntimeError('diagnosis summary exceeds byte budget')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ('pairs','runtime_cpu_seconds','counts')}))


if __name__=='__main__': main()
