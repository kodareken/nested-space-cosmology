#!/usr/bin/env python3
"""Audit retained cap/exterior common values and the existing quiet corner.

No profile transfer, radial solve, prefix solve or time step. Terminal
comparisons do not alter either stored field to manufacture a junction.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate, Budget, plain
from record_nsc_atlas_operator_comparison import subset
from validation_source_successors import verify_saved_source
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice

ROOT=Path(__file__).resolve().parents[1]
REFERENCE='results/summaries/nsc-atlas-cap-fixed-reference-v1.json'
SOURCES=('scripts/check_nsc_atlas_cap_initial_junction.py',
         'scripts/diagnose_nsc_atlas_boundary_rows.py',
         'src/recursive_horizons/nsc_atlas_nonlinear_boundary.py',
         'src/recursive_horizons/nsc_atlas_vacuum_departure.py')


def main():
    budget=Budget(maximum_evaluations=10000,cpu_seconds=15.)
    current,_,cache,identity=authenticate(ROOT,budget)
    raw=(ROOT/REFERENCE).read_bytes();reference=json.loads(raw)
    if (reference['schema']!='nsc-atlas-cap-fixed-reference-v1' or
            reference['input']['sha256']!=identity['sha256']):
        raise ValueError('cap reference does not bind the retained exterior payload')
    for name,expected in reference['source_bindings'].items():
        verify_saved_source(name,expected,ROOT)
    cap=reference['result'];selected=reference['selection']
    R,f,alpha=cap['radius'],cap['wall']['f'],selected['alpha']
    q=math.sqrt(alpha/6)*f
    c=R*f*(math.sinh(q)/q if q else 1.)
    terminal=[]
    for case in current['cases']:
        row=case['snapshots'][selected['cut']];wall=row['physical']['wall']
        differences=dict(radius=R-wall['r'],canonical_time=selected['T']-wall['T'],
                         invariant_scalar=c-wall['c'],alpha=alpha-case['alpha'])
        # This is a necessary common-data check, not an existence criterion.
        matches=all(math.isclose(value,0.,rel_tol=0.,abs_tol=2e-12)
                    for value in differences.values())
        terminal.append(dict(amplitude=case['amplitude'],x=row['x'],
                             cap_minus_stored=differences,common_data_match=matches))
        budget.check(256)

    gauge=subset(cache['gauge_rows'],0)
    source=subset(cache['source_rows']['finite'],0)
    x=float(cache['x_nodes'][0])
    if (x!=.4 or gauge['x']!=x or source['x']!=x or
            any(source[k]!=0. for k in ('c','c_x'))):
        raise ValueError('stored source does not have the quiet corner')
    src=vacuum_slice(0.,x,budget=budget)
    vac=vacuum_slice(gauge['W'],x,budget=budget)
    residual={
        'source_r':source['r']-float(src['r']),
        'source_log_C':source['log_C']-float(src['L']),
        'source_r_x':source['r_x']-float(src['r_x']),
        'source_log_C_x':source['log_C_x']-float(src['log_C_x']),
        'gauge_r':gauge['r']-float(vac['r']),
        'gauge_L':gauge['L']-float(vac['L']),
        'gauge_boost':gauge['eta_plus']+.5*math.log(gauge['W_x']),
    }
    for case in current['cases']:
        first=case['snapshots'][0]
        if first['x']!=x or first['physical']['wall']['c']!=0.:
            raise ValueError('stored initial wall is not the shared quiet event')
        for group in first['state'].values():
            for values in group.values():
                array=np.asarray(values,dtype=float)
                budget.check(int(array.size))
                if np.any(array!=0.):
                    raise ValueError('stored initial departures are not quiet')
        for key in ('r','L','T','s'):
            residual[case['amplitude']+'_initial_'+key]=first['physical']['wall'][key]-gauge[key]
    outgoing=dict(r_V=float(vac['r_x']),c_V=0.,log_C_V=float(vac['log_C_x']))
    common=dict(r=gauge['r'],c=0.,eta_plus=gauge['eta_plus'])
    cap_wall=dict(A=1-gauge['r']**2,delta=0.,w_minus=0.)
    first=hybrid_wall_boundary(common,outgoing,cap_wall,
                               dict(C_wall=math.exp(gauge['L']),log_C_s=0.),alpha=alpha)
    xs=float(first['Vdot_plus'])
    Lx=2*gauge['eta_plus_x']-2*float(first['k_plus'])/xs+2*outgoing['log_C_V']
    hybrid=hybrid_wall_boundary(common,outgoing,cap_wall,
                                dict(C_wall=math.exp(gauge['L']),log_C_s=Lx*xs),alpha=alpha)
    rates=dict(r=float(hybrid['r_s'])/xs,L=Lx,T=float(hybrid['T_s'])/xs,s=1/xs,
               c=float(hybrid['c_s'])/xs)
    for key,name in (('r','r_x_total'),('L','L_x_total'),('T','T_x'),('s','s_x')):
        residual['rate_'+key]=rates[key]-gauge[name]
    residual.update(incoming_scalar=float(hybrid['w_plus']),scalar_rate=rates['c'],
        incoming_r_jet=float(hybrid['r_U_plus'])-float(vac['r_Ubar']),
        incoming_scalar_jet=float(hybrid['c_U_plus']),
        incoming_lapse_jet=float(hybrid['log_C_U_plus'])-float(vac['log_C_Ubar']))
    budget.check(1024)
    if not all(math.isfinite(v) and abs(v)<=2e-12 for v in residual.values()):
        raise ValueError('quiet corner algebra does not reproduce the stored vacuum junction')
    report=dict(schema='nsc-atlas-cap-initial-junction-v1',
        terminal_cap=dict(radius=R,T=selected['T'],alpha=alpha,f=f,invariant_scalar=c),
        terminal_comparisons=terminal,
        quiet_corner=dict(x=x,radius=gauge['r'],alpha=alpha,rates=rates,residuals=residual,
                          common_data_match=True,finite_time_solution=False),
        scope='necessary common-data audit and quiet vacuum junction; no evolved nonzero scalar state',
        input={key:identity[key] for key in ('path','sha256','manifest_sha256','authenticated_input_files')},
        cap_reference=dict(path=REFERENCE,sha256=hashlib.sha256(raw).hexdigest()),
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                       terminal_comparisons=2,quiet_hybrid_calls=2,radial_solves=0,prefix_solves=0,time_steps=0),
        limits=dict(charged_values=10000,cpu_seconds=15.))
    budget.check()
    print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__': main()
