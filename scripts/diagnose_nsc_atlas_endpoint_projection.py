"""Same physical retained background on the three existing P001 cap spaces.

This is an instantaneous counterfactual, not a state recovered from P001.
There is no exterior solve, material-wall match, new mesh or time step.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_limits

from diagnose_nsc_atlas_endpoint_gaps import authenticate
from run_nsc_atlas_fixed_onset import plain,prepare_history
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

ROOT=Path(__file__).resolve().parents[1]
CPU_LIMIT=80.0
VALUE_LIMIT=100_000_000
BYTE_LIMIT=12000


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    if not args.output.parent.is_dir(): raise ValueError('output parent must exist')
    producer=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    budget=Budget(maximum_evaluations=VALUE_LIMIT,cpu_seconds=CPU_LIMIT)
    binding,_,records,_=authenticate(ROOT,'results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json')
    fine=records[2]
    sample=fine['samples'][5]
    if sample['x']!=.45094013016102275: raise ValueError('the admitted common cut changed')
    radius=float(sample['einstein_radius'])
    history_time=float(sample['central_proper_clock_T'])+fine['preparation']['clock_offset']
    alpha=fine['input']['alpha']
    preparation={}
    cases=[]
    with threadpool_limits(limits=1):
        background,offset=prepare_history(budget,preparation,enrich_ordinary_pieces=True)
        if (preparation['input_sha256']!=fine['preparation']['input_sha256']
                or offset!=fine['preparation']['clock_offset']):
            raise ValueError('the regenerated original retained history changed')
        Z=math.atanh(radius)
        jet=background.reduced_jets(history_time,Z,budget=budget)
        incoming=(jet['p']+jet['h'])/radius
        for record in records:
            before=time.process_time()
            basis=FixedCapBasis(tuple(Fraction(v) for v in record['geometry']['exact_cap_knots']),budget=budget)
            zeros=np.zeros(basis.dimension)
            row=fixed_cap_reference(basis,zeros,zeros,Z,alpha=alpha,T_x=1.,radius_rate=0.,
                incoming_w_plus=incoming,budget=budget,background=background,background_time=history_time,
                rtol=2.5e-14,atol=1e-15,background_panel_sides=True,compiled_coefficients=True)
            w=row['wall']
            trace=w['A']**2*math.exp(row['lapse'])*(w['p_weak']*w['h_weak']-w['p_trace']*w['h_trace'])/2
            remainder=row['explicit_time_remainder']
            physical=alpha*(trace+remainder)
            cases.append(dict(level=record['level'],dimension=basis.dimension,
                full_outer_width=float(basis.exact_knots[-1]-basis.exact_knots[-2]),
                clipped_outer_width=Z-float(basis.exact_knots[-2]),
                mass=row['mass'],lapse=row['lapse'],
                A_wall=w['A'],B_wall=w['B'],p_trace=w['p_trace'],h_trace=w['h_trace'],
                p_weak=w['p_weak'],h_weak=w['h_weak'],
                velocity_projection_defect=w['velocity_discrepancy'],
                background_time_wall_projection_gap=row['background']['wall_time_projection_gap'],
                normalized_trace_term=trace,normalized_affine_remainder=remainder,
                physical_trace_term=alpha*trace,physical_affine_remainder=alpha*remainder,
                physical_mass_gap=physical,
                weak_balance_identity_defect=row['remainder_identity_defect'],
                incoming_strong_trace_defect=(w['p_trace']+w['h_trace'])/radius-incoming,
                A_lower_certificate=row['metric_certificate']['A_lower_float_down'],
                CPU_seconds=time.process_time()-before))
            budget.check()
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=producer:
        raise ValueError('producer changed during the diagnostic')
    elapsed=time.process_time()-budget.started
    result=dict(schema='nsc-atlas-endpoint-frozen-projection-v1',comparison=binding,
        source_context_commit=commit,producer_sha256=producer,preparation=preparation,
        fixed_point=dict(source_x=sample['x'],radius=radius,history_time=history_time,
            alpha=alpha,T_x=1.,R_x=0.,incoming_w_plus=incoming,departure_coefficients='identically zero'),
        cases=cases,scope='Counterfactual original background at one common physical cut. Diagnostic rates are not the material-wall dynamics.',
        not_inferred='No evolved departure, actual exterior, trajectory maximum, continuum convergence or error range is recovered.',
        runtime_cpu_seconds=elapsed,charged_values=budget.evaluations,
        counts=dict(original_history_preparations=1,instantaneous_cap_solves=3,new_meshes=0,time_steps=0),
        earlier_format_failure=dict(completed_cap_solves=3,original_history_preparations=1,
            wall_seconds=6.8481225,CPU_seconds=None,
            error='numpy preparation diagnostics rejected by JSON after the three solves; no result values saved',
            correction='reuse existing plain conversion; same state, spaces, equations and parameters'),
        limits=dict(CPU_seconds=CPU_LIMIT,charged_values=VALUE_LIMIT,summary_bytes=BYTE_LIMIT))
    raw=json.dumps(plain(result),indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>BYTE_LIMIT: raise RuntimeError('frozen projection summary exceeds its budget')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(CPU_seconds=elapsed,charged_values=budget.evaluations,
        gaps=[r['physical_mass_gap'] for r in cases],
        remainders=[r['normalized_affine_remainder'] for r in cases])))


if __name__=='__main__': main()
