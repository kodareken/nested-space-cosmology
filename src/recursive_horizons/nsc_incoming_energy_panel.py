"""Uniform-real-energy Chebyshev column preflight, not a source-window result.

The scalar-phase gauge leaves the vacuum projector unchanged. Chebyshev
degree25 leakage remains part of the continuous residual of degree24 data.
"""
from math import ceil
from time import perf_counter

import mpmath as mp
import sympy as sp

from .nsc_incoming_matched_horizon_initializer import _geometry, _coefficients
from .nsc_incoming_validated_vacuum_propagation import (
    DEGREE, RADIUS, MAX_STEP, ExactProfileGenerator, generator_coefficients,
    _cmid, _ivc, _norm, _upper_float, point, unpoint, cpoint,
)
from .nsc_incoming_source_quadrature_bound import pack, unpack, float_enclosure
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range


ENERGY_DEGREE = 24


def panel_identities():
    """Exact phase, Chebyshev multiplication and endpoint-degree checks."""
    x = sp.Symbol('x', real=True)
    cheb = [sp.chebyshevt(j,x) for j in range(ENERGY_DEGREE+2)]
    residuals = [sp.expand(x*cheb[0]-cheb[1])]
    residuals += [sp.expand(x*cheb[j]-(cheb[j-1]+cheb[j+1])/2) for j in range(1,ENERGY_DEGREE+1)]
    z = sp.Symbol('z'); C=sp.Matrix(2,2,sp.symbols('c0:4'))
    phase = z*sp.eye(2)*C-C*z*sp.eye(2)
    if residuals != [0]*25 or phase != sp.zeros(2):
        raise ArithmeticError('energy multiplication or scalar phase identity failed')
    return {'scalar_commutator': [str(v) for v in phase],
            'Chebyshev_multiplication_residuals': [str(v) for v in residuals],
            'phase_gauge': 'G->G+(E/B)I, diag(0,2E/B), same projector',
            'retained_energy_degree': 24, 'residual_energy_degree': 25,
            'unresolved_edge': 'lower component degree25 is 4*(1/B)*c24_lower',
            'raw_covariance_degree': 48, 'linear_vertex_kernel_degree': 49,
            'endpoint_normalization_permitted': False,
            'uniform_norm': 'sup_real_x norm(sum Rk Tk) <= sum norm(Rk)'}


def _complex_report(value):
    return {'real':pack(value.real), 'imag':pack(value.imag)}


def initial_interpolant(channel, config, sign, *, degree=24):
    """Directed Chebyshev-Lobatto DCT of the explicit matched initial column."""
    if degree != ENERGY_DEGREE or sign not in (1,-1):
        raise ValueError('fixeddegree24 and actual angular sign required')
    g=_geometry(channel,config)
    g['mass']=mp.iv.mpf(channel['compact_mass'])
    g['angular']=mp.iv.mpf(channel['angular_eigenvalue'])
    # Closed initial formula only: no propagated complex-energy assertion.
    rho=mp.iv.mpf(8);axis=2*(rho+1/rho);height=2*(rho-1/rho)
    Ebox=mp.iv.mpc(_range(28-_hi(axis),28+_hi(axis)),_range(-_hi(height),_hi(height)))
    c=_coefficients(g,Ebox,sign)
    w=mp.iv.sqrt(g['delta0'])*(c['w0']+g['delta0']*c['w1'])
    s=mp.iv.mpf(_hi(abs(w)));product=s*s
    if _hi(product)>=1:raise ArithmeticError('initial normalization pole not excluded on ellipse')
    root=mp.iv.sqrt(1-product)
    component_M=[product/(root*(1+root)),s/root]
    component_error=[4*M*rho**(-degree)/(rho-1) for M in component_M]
    interpolation=_norm(component_error)
    nodes=[];values=[]
    for j in range(degree+1):
        x=mp.iv.cos(mp.iv.pi*j/degree);E=28+4*x
        c=_coefficients(g,E,sign)
        w=mp.iv.sqrt(g['delta0'])*(c['w0']+g['delta0']*c['w1'])
        normalization=mp.iv.sqrt(1+w.real*w.real+w.imag*w.imag)
        values.append([1/normalization-1,w/normalization])
        nodes.append(pack(x))
    coefficients=[];enclosures=[];rounding=mp.iv.mpf(0)
    for k in range(degree+1):
        coefficient=[mp.iv.mpc(0),mp.iv.mpc(0)]
        for j,v in enumerate(values):
            weight=mp.iv.mpf('0.5') if j in(0,degree) else mp.iv.mpf(1)
            weight*=mp.iv.cos(mp.iv.pi*k*j/degree)*2/degree
            if k in(0,degree):weight/=2
            coefficient=[a+weight*b for a,b in zip(coefficient,v)]
        if k==0:coefficient[0]+=1  # The exact constant vacuum column.
        center=[_cmid(v) for v in coefficient]
        rounding+=_norm([v-_ivc(c) for v,c in zip(coefficient,center)])
        coefficients.append(center);enclosures.append([_complex_report(v) for v in coefficient])
    return coefficients,{'sign':sign,'degree':degree,'ellipse_parameter':8,
                         'energy_rectangle':{'real':[11.75,44.25],'imag':[-15.75,15.75]},
                         'ratio_modulus_upper':_upper_float(s),'normalization_product_upper':_upper_float(product),
                         'component_ellipse_bounds':[_upper_float(v) for v in component_M],
                         'component_interpolation_error_upper':[_upper_float(v) for v in component_error],
                         'column_interpolation_error_upper':_upper_float(interpolation),
                         'DCT_coefficient_radius_l1_upper':_upper_float(rounding),
                         'total_initial_column_error_upper':_upper_float(interpolation+rounding),
                         'nodes':nodes,'coefficient_enclosures':enclosures,
                         'column_coefficients':[[cpoint(v) for v in c] for c in coefficients],
                         'theorem':'ATAP8.2 Eq8.3 applied componentwise to v-e1; combine Euclidean component bounds',
                         'complex_policy':'coefficient-conjugate normalization at same complex E; no conjugation of E'}


def shared_geometry(channel, config, center):
    """Reuse the frozen point generator coefficient proof once for both signs."""
    generator=ExactProfileGenerator(channel,config,1,24)
    disk=generator.disk(center)
    old=generator_coefficients(generator,center,disk)
    m=mp.iv.mpf(channel['compact_mass']);ell=mp.iv.mpf(channel['angular_eigenvalue'])
    if _lo(m)<=0 or _lo(ell)<=0:raise ValueError('this group32 basis factorization requires nonzero frozen labels')
    # u=1/B, fr=r*sqrt(delta/B), f=sqrt(delta/B).
    basis=[(D/24,-(plus+minus)/(2*m),(plus-minus)/(2*mp.iv.j*ell)) for D,plus,minus in old]
    rectangle=mp.iv.mpc(_range(center-RADIUS,center+RADIUS),_range(-RADIUS,RADIUS))
    D,plus,minus=generator.values(rectangle)
    M=mp.iv.sqrt(abs(plus)**2+abs(minus)**2+abs(64*D/24)**2)
    return basis,{'center':point(center),'geometry_disk':disk,
                  'phase_gauged_generator_norm_upper':_upper_float(M),
                  'basis_coefficient_enclosures':[[_complex_report(v) for v in row] for row in basis],
                  'shared_between_signs':True,'geometry_coefficient_preparations':1}


def _apply_energy(g, coefficients, sign, mass, angular):
    """Full degree25 product; the edge is returned even during projected steps."""
    u,fr,f=g;plus=-mass*fr+mp.iv.j*sign*angular*f;minus=-mass*fr-mp.iv.j*sign*angular*f
    size=len(coefficients);out=[]
    zero=mp.mpc(0)
    def lower(k):return _ivc(coefficients[k][1]) if 0<=k<size else mp.iv.mpc(0)
    for k in range(size+1):
        a,b=coefficients[k] if k<size else (zero,zero)
        if k==0:xb=lower(1)/2
        elif k==1:xb=lower(0)+lower(2)/2
        else:xb=(lower(k-1)+lower(k+1))/2
        out.append([minus*b,plus*a+56*u*b+8*u*xb])
    return out


def _sum_panels(values):
    result=[[mp.iv.mpc(0),mp.iv.mpc(0)] for _ in range(ENERGY_DEGREE+2)]
    for value in values:
        result=[[a+b for a,b in zip(row,addition)] for row,addition in zip(result,value)]
    return result


def replay_step(row):
    """Uniform-real-energy continuous residual; degree25 is never discarded."""
    h=mp.iv.mpf(unpoint(row['step']))
    low=mp.iv.mpf(0);high=mp.iv.mpf(0);edge=mp.iv.mpf(0)
    for n,values in enumerate(row['residual_time_energy_norm_upper']):
        factor=h**(n+1)/(n+1)
        edge+=mp.iv.mpf(values[25])*factor
        term=sum((mp.iv.mpf(v) for v in values[:25]),mp.iv.mpf(0))*factor
        if n<DEGREE:low+=term
        else:high+=term
    ratio=h/RADIUS
    tail=mp.iv.mpf(row['uniform_generator_norm_upper'])*mp.iv.mpf(row['column_polynomial_Cheb_l1_upper'])*h*ratio**37/(38*(1-ratio))
    total=low+high+edge+tail+mp.iv.mpf(row['centering_coefficient_l1_error_upper'])
    return {'projected_low_order_arithmetic_upper':_upper_float(low),
            'temporal_polynomial_tail_upper':_upper_float(high),
            'energy_degree25_residual_upper':_upper_float(edge),
            'analytic_temporal_generator_tail_upper':_upper_float(tail),
            'local_uniform_column_error_upper':_upper_float(total),
            'local_uniform_column_tolerance':1e-20,'local_gate_pass':_hi(total)<=mp.mpf('1e-20')}


def energy_step(basis, geometry, initial, channel, sign, step):
    """One degree36 time/degree24 energy polynomial, with complete residual."""
    start=perf_counter();mass=mp.iv.mpf(channel['compact_mass']);angular=mp.iv.mpf(channel['angular_eigenvalue'])
    coefficients=[initial]
    for n in range(DEGREE):
        value=_sum_panels(_apply_energy(basis[j],coefficients[n-j],sign,mass,angular) for j in range(n+1))
        # Projection defines the approximation. Its missing degree25 is
        # explicitly evaluated below in the physical equation's residual.
        coefficients.append([[_cmid(-mp.iv.j*v/(n+1)) for v in row] for row in value[:25]])
    norms=[]
    for n in range(2*DEGREE+1):
        value=_sum_panels(_apply_energy(basis[j],coefficients[n-j],sign,mass,angular)
                          for j in range(max(0,n-DEGREE),min(DEGREE,n)+1))
        value=[[mp.iv.j*v for v in row] for row in value]
        if n<DEGREE:
            for k in range(25):
                value[k]=[v+(n+1)*_ivc(c) for v,c in zip(value[k],coefficients[n+1][k])]
        norms.append([_upper_float(_norm(row)) for row in value])
    h=mp.iv.mpf(step);output=[[mp.iv.mpc(0),mp.iv.mpc(0)] for _ in range(25)]
    for coefficient in reversed(coefficients):
        output=[[v*h+_ivc(c) for v,c in zip(row,addition)] for row,addition in zip(output,coefficient)]
    next_center=[[_cmid(v) for v in row] for row in output]
    centering=sum((_norm([v-_ivc(c) for v,c in zip(row,mid)]) for row,mid in zip(output,next_center)),mp.iv.mpf(0))
    poly_norm=sum((sum((_norm([_ivc(v) for v in row]) for row in coefficient),mp.iv.mpf(0))*h**n
                   for n,coefficient in enumerate(coefficients)),mp.iv.mpf(0))
    row={'sign':sign,'step':point(step),'time_degree':36,'energy_degree':24,
         'residual_time_energy_norm_upper':norms,
         'uniform_generator_norm_upper':geometry['phase_gauged_generator_norm_upper'],
         'column_polynomial_Cheb_l1_upper':_upper_float(poly_norm),
         'centering_coefficient_l1_error_upper':_upper_float(centering),
         'column_end_coefficients':[[cpoint(v) for v in c] for c in next_center]}
    row.update(replay_step(row));row['elapsed_seconds']=perf_counter()-start
    return row


def preflight(channel, config, matched_receipt):
    """Initial DCT/ellipse certificate plus ONE shared-geometry time step."""
    start=perf_counter()
    with _precision(80):
        g=_geometry(channel,config);center=(_lo(mp.iv.ln(g['delta0']))+_hi(mp.iv.ln(g['delta0'])))/2
        began=perf_counter();basis,geometry=shared_geometry(channel,config,center)
        geometry['elapsed_seconds']=perf_counter()-began
        cases={}
        for sign in(1,-1):
            began=perf_counter();initial,certificate=initial_interpolant(channel,config,sign)
            certificate['elapsed_seconds']=perf_counter()-began
            endpoint=[r for r in matched_receipt['certificate']['per_sign'] if r['sign']==sign][0]
            if matched_receipt['certificate']['frozen_source_kappa']!=config['surface_gravity']:
                raise ValueError('source occupation kappa changed')
            row=energy_step(basis,geometry,initial,channel,sign,mp.mpf(1)/32)
            cases[str(sign)]={'initial_interpolant':certificate,'step':row,
                              'inherited_mathematical_projector_error_upper':endpoint['projector_error_upper'],
                              'initial_projector_error_added_to_column_error':False}
        stop=(_lo(mp.iv.ln(g['qh']-3*mp.iv.pi/4))+_hi(mp.iv.ln(g['qh']-3*mp.iv.pi/4)))/2
        steps=int(mp.ceil((stop-center)/MAX_STEP))
        measured=geometry['elapsed_seconds']+sum(v['step']['elapsed_seconds'] for v in cases.values())
        return {'group':32,'energy_interval':[24,32],'shared_geometry':geometry,'cases':cases,
                'identities':panel_identities(),'elapsed_seconds':perf_counter()-start,
                'future_work':{'real_steps_at_current_spacing':steps,'signs':2,
                               'geometry_preparations':steps,
                               'energy_operator_applications_per_sign_step':2035,
                               'padded_energy_vector_products_per_operator':26,
                               'full_upper_vector_product_count':steps*2*2035*26,
                               'linear_runtime_estimate_seconds':steps*measured,
                               'planning_allowance_seconds':2*steps*measured,
                               'runtime_interpretation':'measured-first-step extrapolation and2x planning allowance, not a wall-clock guarantee',
                               'authorized_full_panel':False},
                'scope':{'only_one_time_step':True,'both_angular_signs':True,
                         'frozen_binary_mass_angular_everywhere':True,'source_kappa_changed':False,
                         'whole_panel_field_error_or_source_integral':False,
                         'metric_evolution_or_physical_IV':False,'maximum_CPU_workers':1}}
