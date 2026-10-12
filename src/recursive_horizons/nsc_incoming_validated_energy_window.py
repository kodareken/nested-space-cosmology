"""Full fixed group32 E24–32 panel using the frozen energy preflight method.

The successful preflight step is reused. Geometry coefficients are shared
through locked per-step cache files. No endpoint column is normalized.
"""
import fcntl
import json
from pathlib import Path
from time import process_time

import mpmath as mp
import numpy as np

from .nsc_incoming_energy_panel import shared_geometry, energy_step, replay_step
from .nsc_incoming_matched_horizon_initializer import _geometry
from .nsc_incoming_validated_vacuum_propagation import (
    ExactProfileGenerator, RADIUS, MAX_STEP, point, unpoint, cpoint, uncpoint,
    _upper_float, _ivc,
)
from .nsc_incoming_source_quadrature_bound import (
    pack, unpack, float_enclosure, certify_gauss_brackets, interval_bloch,
)
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range
from .nsc_incoming_state_moments import incoming_group_factor


class WindowObstruction(RuntimeError):
    """Fixed-resolution gate or authorized work budget needs root review."""


def _decode_coefficients(values):
    return [[uncpoint(v) for v in row] for row in values]


def _basis(geometry):
    return [tuple(mp.iv.mpc(unpack(v['real']),unpack(v['imag'])) for v in row)
            for row in geometry['basis_coefficient_enclosures']]


def _write_json(path, value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.writing')
    temporary.write_text(json.dumps(value,sort_keys=True,allow_nan=False))
    temporary.replace(path)


def geometry_at(channel,config,index,center,directory,key):
    """One actual geometry coefficient calculation per shared real step."""
    path=Path(directory)/f'geometry-{index:04d}.json'
    with path.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if path.exists():
            data=json.loads(path.read_text())
            if data['key']!=key or data['index']!=index or data['geometry']['center']!=point(center):
                raise ValueError('shared geometry cache belongs to changed inputs')
        else:
            _,geometry=shared_geometry(channel,config,center)
            data={'key':key,'index':index,'geometry':geometry}
            _write_json(path,data)
    return _basis(data['geometry']),data['geometry']


def uniform_generator_norm(channel,config,center):
    generator=ExactProfileGenerator(channel,config,1,24)
    y=mp.iv.mpc(_range(center-RADIUS,center+RADIUS),_range(-RADIUS,RADIUS))
    D,plus,minus=generator.values(y)
    return mp.iv.sqrt(abs(plus)**2+abs(minus)**2+abs(64*D/24)**2)


def factor_and_vertices(channel):
    a=mp.iv.sqrt(3*mp.iv.pi/2-4);r=mp.iv.sqrt(2)
    m=mp.iv.mpf(channel['compact_mass']);ell=mp.iv.mpf(channel['angular_eigenvalue'])
    exact=channel['copy_count']*channel['degeneracy']/(8*mp.iv.pi**2*r*r*a)
    stored=mp.mpf(incoming_group_factor(channel,(1,-1)))
    factor=_range(min(_lo(exact),stored),max(_hi(exact),stored))
    H=mp.iv.sqrt(m*m+(ell/r)**2+(32/a)**2)
    return a,r,m,ell,factor,(H,32/a,32/a,ell/(2*r))


def sign_error_bounds(channel,column_error,initial_projector):
    a,r,_,_,factor,vertices=factor_and_vertices(channel)
    e=mp.iv.mpf(column_error);normal=(2+e)*e
    projector=mp.iv.mpf(initial_projector)+normal
    stress=[2*factor*8*norm*projector for norm in vertices]
    # The initial pure-projector discrepancy has zero trace. The only
    # vacuum-current error of the raw polynomial is its column norm defect.
    stress[2]=factor*8*vertices[2]*normal
    lapse=4*mp.iv.pi*a*r*r*stress[0]
    return {'uniform_column_error_upper':_upper_float(e),
            'uniform_projector_error_upper':_upper_float(projector),
            'uniform_raw_norm_defect_upper':_upper_float(normal),
            'field_stress_error_upper':[_upper_float(v) for v in stress],
            'field_lapse_error_upper':_upper_float(lapse)}


def propagate_panel_sign(channel,config,sign,preflight,*,directory,key,cpu_budget,progress=None):
    """Reuse first step, then run missing fixed steps with a separate radius."""
    with _precision(80):
        directory=Path(directory);directory.mkdir(exist_ok=True)
        runlock=(directory/f'sign{sign}.run.lock').open('a')
        try:
            fcntl.flock(runlock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            runlock.close();raise RuntimeError('a live worker already owns this sign; no duplicate restart')
        try:
            source=preflight['cases'][str(sign)]
            g=_geometry(channel,config)
            start=unpoint(preflight['shared_geometry']['center'])
            true_start=mp.iv.ln(g['delta0']);true_end=mp.iv.ln(g['qh']-3*mp.iv.pi/4)
            stop=(_lo(true_end)+_hi(true_end))/2
            count=int(mp.ceil((stop-start)/MAX_STEP))
            grid=[start+j*MAX_STEP for j in range(count)]+[stop]
            initial=source['initial_interpolant']['total_initial_column_error_upper']
            initial_projector=_hi(unpack(source['inherited_mathematical_projector_error_upper']['binary_interval']))
            start_error=_upper_float(mp.iv.mpf(preflight['shared_geometry']['phase_gauged_generator_norm_upper'])*abs(true_start-mp.iv.mpf(start)))
            end_error=_upper_float(uniform_generator_norm(channel,config,stop)*abs(true_end-mp.iv.mpf(stop)))
            binding={'key':key,'group':32,'sign':sign,'energy_interval':[24,32],
                     'start':point(start),'stop':point(stop),'steps':count,
                     'time_degree':36,'energy_degree':24,'channel':channel,'config':config}
            first=dict(source['step'],center=point(start),index=0,reused_preflight=True)
            if first['step']!=point(MAX_STEP) or not first['local_gate_pass']:
                raise ValueError('successful exact first preflight step is required for reuse')
            rows=[first];cpu_previous=0.
            checkpoint=directory/f'sign{sign}.checkpoint.json'
            if checkpoint.exists():
                saved=json.loads(checkpoint.read_text())
                if saved['binding']!=binding:raise ValueError('panel checkpoint input binding changed')
                rows=saved['steps'];cpu_previous=saved['cpu_seconds']
                if not rows or rows[0]!=first:raise ValueError('cached panel no longer reuses the certified preflight')
            for j,row in enumerate(rows):
                if row['center']!=point(grid[j]) or row['step']!=point(grid[j+1]-grid[j]):
                    raise ValueError('panel checkpoint has a real-domain gap or overlap')
                replay=replay_step(row)
                if any(row[k]!=v for k,v in replay.items()):raise ValueError('panel checkpoint residual changed')
            coefficients=_decode_coefficients(rows[-1]['column_end_coefficients'])
            accumulated=mp.iv.mpf(initial)+start_error+end_error+sum((mp.iv.mpf(row['local_uniform_column_error_upper']) for row in rows),mp.iv.mpf(0))
            began=process_time()
            def checkpoint_now():
                _write_json(checkpoint,{'binding':binding,'steps':rows,
                                        'cpu_seconds':cpu_previous+process_time()-began})
            for j in range(len(rows),count):
                if cpu_previous+process_time()-began>=cpu_budget:
                    checkpoint_now()
                    raise WindowObstruction(json.dumps({'reason':'fixed CPU planning allowance reached','sign':sign,'completed_steps':len(rows),'checkpoint':str(checkpoint)}))
                basis,geometry=geometry_at(channel,config,j,grid[j],directory,key)
                row=energy_step(basis,geometry,coefficients,channel,sign,grid[j+1]-grid[j])
                row.update(center=point(grid[j]),index=j,reused_preflight=False,
                           geometry_B_real_lower=geometry['geometry_disk']['B_real_lower'],
                           geometry_sin_real_lower=geometry['geometry_disk']['sin_real_lower'])
                rows.append(row);coefficients=_decode_coefficients(row['column_end_coefficients'])
                accumulated+=mp.iv.mpf(row['local_uniform_column_error_upper'])
                if (j+1)%16==0 or j+1==count:
                    checkpoint_now()
                    bounds=sign_error_bounds(channel,_hi(accumulated),initial_projector)
                    if progress:progress({'sign':sign,'completed_steps':len(rows),'total_steps':count,
                                          'cpu_seconds':cpu_previous+process_time()-began,**bounds})
                    if bounds['field_lapse_error_upper']>3e-12:
                        raise WindowObstruction(json.dumps({'reason':'fixed cumulative field enclosure exceeds full window budget','checkpoint':str(checkpoint),**bounds}))
            checkpoint_now()
            return {'binding':binding,'steps':rows,'initial_column_error_upper':initial,
                    'initial_projector_error_upper':_upper_float(mp.iv.mpf(initial_projector)),
                    'start_coordinate_error_upper':start_error,'end_coordinate_error_upper':end_error,
                    'column_end_coefficients':rows[-1]['column_end_coefficients'],
                    'cpu_seconds':cpu_previous+process_time()-began,'preflight_reused_once':True}
        finally:
            runlock.close()


def replay_panel_sign(payload):
    """Uniform field error and endpoint coefficients, no field generation."""
    with _precision(80):
        b=payload['binding'];rows=payload['steps']
        if b['group']!=32 or b['energy_interval']!=[24,32] or b['time_degree']!=36 or b['energy_degree']!=24 or len(rows)!=b['steps']:
            raise ValueError('complete fixed group32 E24–32 propagation required')
        at=unpoint(b['start']);error=mp.iv.mpf(payload['initial_column_error_upper'])+payload['start_coordinate_error_upper']+payload['end_coordinate_error_upper']
        energy_error=mp.iv.mpf(0);max_local=0.
        for j,row in enumerate(rows):
            h=unpoint(row['step'])
            if row['index']!=j or unpoint(row['center'])!=at or not 0<h<=MAX_STEP or row['reused_preflight']!=(j==0):
                raise ValueError('panel domain or reused preflight changed')
            replay=replay_step(row)
            if any(row[k]!=v for k,v in replay.items()):raise ValueError('panel continuous residual replay differs')
            if j and min(row['geometry_B_real_lower'],row['geometry_sin_real_lower'])<=0:raise ValueError('panel geometry disk margin failed')
            error+=mp.iv.mpf(row['local_uniform_column_error_upper']);energy_error+=mp.iv.mpf(row['energy_degree25_residual_upper'])
            max_local=max(max_local,row['local_uniform_column_error_upper']);at+=h
        if at!=unpoint(b['stop']) or payload['column_end_coefficients']!=rows[-1]['column_end_coefficients']:
            raise ValueError('incoming endpoint changed')
        return {'sign':b['sign'],'steps':len(rows),'energy_interval':[24,32],
                'cumulative_energy_degree25_column_error_upper':_upper_float(energy_error),
                'maximum_local_column_error_upper':max_local,
                'column_end_coefficients':payload['column_end_coefficients'],
                'cpu_seconds':payload['cpu_seconds'],**sign_error_bounds(b['channel'],_hi(error),payload['initial_projector_error_upper'])}


def _moment(n):
    return mp.iv.mpf(0) if n%2 else mp.iv.mpf(2)/(1-n*n)


def _energy_moment(n):
    xm=_moment(1) if n==0 else (_moment(n-1)+_moment(n+1))/2
    return 28*_moment(n)+4*xm


def raw_polynomial_integral(channel,sign,coefficients):
    """Exact Chebyshev moments of degree48 covariance and linear vertices."""
    columns=[[_ivc(v) for v in row] for row in _decode_coefficients(coefficients)]
    covariance=[[[mp.iv.mpc(0) for _ in range(2)] for _ in range(2)] for _ in range(49)]
    for i,left in enumerate(columns):
        for j,right in enumerate(columns):
            conjugate=[mp.iv.mpc(v.real,-v.imag) for v in right]
            for n in(i+j,abs(i-j)):
                for a in range(2):
                    for b in range(2):covariance[n][a][b]+=left[a]*conjugate[b]/2
    ordinary=[mp.iv.mpc(0) for _ in range(4)];weighted=[mp.iv.mpc(0) for _ in range(4)]
    for n,C in enumerate(covariance):
        vector=[C[0][0]+C[1][1],C[0][1]+C[1][0],mp.iv.j*(C[0][1]-C[1][0]),C[0][0]-C[1][1]]
        ordinary=[v+4*_moment(n)*c for v,c in zip(ordinary,vector)]
        weighted=[v+4*_energy_moment(n)*c for v,c in zip(weighted,vector)]
    a=mp.iv.sqrt(3*mp.iv.pi/2-4);r=mp.iv.sqrt(2)
    m=mp.iv.mpf(channel['compact_mass']);ell=sign*mp.iv.mpf(channel['angular_eigenvalue'])
    parallel=-weighted[3]/a
    result=[-m*ordinary[1]+ell/r*ordinary[2]+parallel,parallel,
            (weighted[0]-224)/a,ell/(2*r)*ordinary[2]]
    if any(not(_lo(v.imag)<=0<=_hi(v.imag)) for v in result):raise ArithmeticError('polynomial source lost real-axis Hermiticity')
    return [v.real for v in result]


def reference_kernel(channel,sign,E):
    """The same full ad4 Bloch reference, with frozen binary mass/angular."""
    a=mp.iv.sqrt(3*mp.iv.pi/2-4);r=mp.iv.sqrt(2)
    m=mp.iv.mpf(channel['compact_mass']);ell=sign*mp.iv.mpf(channel['angular_eigenvalue'])
    terms=interval_bloch()(3*mp.iv.pi/4,m,ell,E)
    b=[sum((terms[j][i,0] for j in range(5)),mp.iv.mpf(0)) for i in range(3)]
    parallel=-E/a*b[2]
    return[-m*b[0]+ell/r*b[1]+parallel,parallel,mp.iv.mpc(0),ell/(2*r)*b[1]]


def reference_integral(channel,sign,rule):
    """Directed Gauss48 and analytic disk bound for reference alone."""
    nodes=certify_gauss_brackets(rule['roots'],48)
    samples=[];quadrature=[mp.iv.mpf(0) for _ in range(4)]
    for x,w in nodes:
        values=reference_kernel(channel,sign,28+4*x)
        samples.append([pack(v.real) for v in values])
        quadrature=[q+4*w*v.real for q,v in zip(quadrature,values)]
    maxima=[mp.mpf(0)]*4
    # E diskcenter28 radius8 has ReE>=20>|ImE|; every non-LLL
    # energy-square gap has Re>0 and its positive-root branch is analytic.
    for arc in range(32):
        angle=2*mp.iv.pi*_range(mp.mpf(arc)/32,mp.mpf(arc+1)/32)
        E=mp.iv.mpc(28+8*mp.iv.cos(angle),8*mp.iv.sin(angle))
        for j,value in enumerate(reference_kernel(channel,sign,E)):
            maxima[j]=max(maxima[j],_hi(abs(value)))
    if not all(mp.isfinite(v) for v in maxima):raise ArithmeticError('reference disk bound is not finite')
    maximum_floats=[float_enclosure(mp.iv.mpf(v))[1] if v else 0. for v in maxima]
    error=[mp.iv.mpf(16)*mp.iv.mpf(M)*(mp.iv.mpf('.5')**96)/(1-mp.iv.mpf('.5')) for M in maximum_floats]
    return {'sign':sign,'samples':samples,'rule':rule,
            'quadrature_interval':[pack(v) for v in quadrature],
            'circle_maxima_upper':maximum_floats,
            'analytic_error_upper':[_upper_float(v) for v in error],
            'analytic_domain':'E diskcenter28 radius8, ReE>=20>|ImE|, gapped positive-root branch',
            'reference_current_in_polynomial_identity':True}


def source_from_endpoints(channel,sign_payloads,references):
    """Complete direct vacuum approximation, with all arithmetic enclosed."""
    source=[mp.iv.mpf(0) for _ in range(4)];quadrature=[mp.iv.mpf(0) for _ in range(4)]
    per_sign={}
    for sign in(1,-1):
        field=sign_payloads[str(sign)];reference=references[str(sign)]
        polynomial=raw_polynomial_integral(channel,sign,field['column_end_coefficients'])
        nodes=certify_gauss_brackets(reference['rule']['roots'],48)
        values=[sum((4*w*unpack(row[j]) for (_,w),row in zip(nodes,reference['samples'])),mp.iv.mpf(0)) for j in range(4)]
        if [pack(v) for v in values]!=reference['quadrature_interval']:raise ValueError('reference Gauss contraction changed')
        proof=[_upper_float(mp.iv.mpf(16)*mp.iv.mpf(M)*(mp.iv.mpf('.5')**96)/(1-mp.iv.mpf('.5')))
               for M in reference['circle_maxima_upper']]
        if proof!=reference['analytic_error_upper']:raise ValueError('reference analytic remainder replay changed')
        a,r,_,_,factor,_=factor_and_vertices(channel)
        piece=[factor*(v-q) for v,q in zip(polynomial,values)]
        source=[s+p for s,p in zip(source,piece)]
        quadrature=[q+factor*mp.iv.mpf(e) for q,e in zip(quadrature,reference['analytic_error_upper'])]
        per_sign[str(sign)]={'raw_polynomial_integral':[pack(v) for v in polynomial],
                             'reference_quadrature_integral':reference['quadrature_interval'],
                             'direct_vacuum_integral':[pack(v) for v in piece]}
    raw=[float((_lo(v)+_hi(v))/2) for v in source]
    arithmetic=[_upper_float(abs(v-mp.iv.mpf(x))) for v,x in zip(source,raw)]
    approx=list(raw);approx[2]=0.;arithmetic[2]=0.
    return {'direct_vacuum_source':approx,'raw_polynomial_vacuum_source':raw,
            'exact_vacuum_current':0.,
            'current_policy':'true vacuum projector and ad4 reference each have trace1; exact T01=0 used explicitly, raw polynomial current retained without normalization',
            'arithmetic_interval':[pack(v) for v in source],
            'source_arithmetic_error_upper':arithmetic,
            'reference_analytic_quadrature_error_upper':[_upper_float(v) for v in quadrature],
            'per_sign':per_sign,'raw_covariance_degree':48,'weighted_kernel_degree':49,
            'polynomial_integration':'exact Chebyshev moments with directed arithmetic',
            'endpoint_column_normalized':False}
