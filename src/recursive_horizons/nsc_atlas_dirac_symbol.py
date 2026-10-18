"""Actual radar-metric jets for the existing fourth-order Dirac symbol core.

Only local, one-sided subtraction data are supplied. These formal orders
are not a Gaussian state, a wall transmission theorem or a vacuum stress.
Geometry and analytic field conventions are owned by the existing atlas.
"""
from __future__ import annotations

import numpy as np

from .nsc_atlas_dirac_compact_source import K_WALL
from .nsc_spatial_reference_symbol import (
    INDEX,SymbolJet,SIGMA,reference_projector)

PAULI_CYCLE=(np.eye(2)+1j*SIGMA.sum(axis=0))/2


def _wall_geometry(chart,r):
    d=chart.parent_mass-chart.central_mass
    bm=d/(K_WALL*r*r)+3*r/4
    f=1-2*chart.central_mass/r-r*r
    v=(bm*bm-f).power(.5)
    return bm,bm-K_WALL*r,v


def _ode_series(initial,rhs):
    """Taylor coefficients from y'=rhs(y), evaluated in the same jet algebra."""
    y=SymbolJet.constant(initial)
    for n in range(4):
        y.data[INDEX[n+1,0,0]]=rhs(y).data[INDEX[n,0,0]]/(n+1)
    return y


def _compose(series,increment):
    out=SymbolJet.constant(series.value[:,0,0])
    power=SymbolJet.constant(1.)
    for n in range(1,5):
        power=power*increment
        out+=power*SymbolJet.constant(series.data[INDEX[n,0,0]][:,0,0])
    return out


def _optical_series(chart,R,initial,*,cap,direction):
    radius=_ode_series(R,lambda r:_wall_geometry(chart,r)[2])
    bm,bp,v=_wall_geometry(chart,radius)
    rate=1/((bm if cap else bp)+direction*v)
    optical=SymbolJet.constant(initial)
    for n in range(1,5):
        optical.data[INDEX[n,0,0]]=rate.data[INDEX[n-1,0,0]]/n
    return radius,optical


def radar_metric_jets(chart,U,V,*,region='nested'):
    """Analytic Taylor jets in actual (T,x), to physical order four.

    U=T-x,V=T+x are wall proper-clock labels. At x=0, explicit cap or
    exterior requests select one-sided jets; nested defaults to the cap.
    chart is the existing exact CompactMassChart, including b=0.
    """
    if (np.ndim(U) or np.ndim(V) or not np.isfinite([U,V]).all()
            or not 0<=U<=chart.duration or not 0<=V<=chart.duration):
        raise ValueError('one event inside the actual proper-wall diamond required')
    if region not in ('nested','cap','exterior'):
        raise ValueError('nested, cap or exterior region required')
    cap=U>=V if region=='nested' else region=='cap'
    side='cap' if cap else 'exterior'
    Ru,Rv=chart.radius_at_clock(np.array([U,V]))
    maps_u,maps_v=chart.maps(Ru),chart.maps(Rv)
    field=chart.from_wall_radii(Ru,Rv,region=side)
    dt=SymbolJet.variable(0.,0);dx=SymbolJet.variable(0.,1)
    ru,u=_optical_series(chart,Ru,maps_u[1 if cap else 3],cap=cap,direction=1)
    rv,v=_optical_series(chart,Rv,maps_v[2 if cap else 4],cap=cap,direction=-1)
    ru,u=_compose(ru,dt-dx),_compose(u,dt-dx)
    rv,v=_compose(rv,dt+dx),_compose(v,dt+dx)
    optical=(v-u)/2
    mass=chart.central_mass if cap else chart.parent_mass
    def F(r):
        return 1-2*mass/r-r*r if cap else 1-2*mass/r
    radius_series=_ode_series(float(field['radius']),F)
    r=_compose(radius_series,optical-SymbolJet.constant(optical.value[:,0,0]))
    bmu,bpu,vu=_wall_geometry(chart,ru)
    bmv,bpv,vv=_wall_geometry(chart,rv)
    p=1/((bmu if cap else bpu)+vu)
    q=1/((bmv if cap else bpv)-vv)
    N=(F(r)*p*q).power(.5)
    return dict(lapse=N,radial_length=N,radius=r,shift=SymbolJet.constant(0.),
        potential_per_k=N/r,region=side,clock_u=float(U),clock_v=float(V))


def radar_reference_projector(chart,U,V,momenta,angular_sectors,*,region='nested'):
    """Orders zero through four, reusing the existing Weyl/transport algebra.

    Momentum is conjugate to radar x, not an incoming source frequency.
    Both p and k scale in a joint UV limit. The symbol is a subtraction
    input; pointwise positivity or an exact projector is not asserted.
    """
    p=np.asarray(momenta)
    k=np.asarray(angular_sectors)
    if (p.dtype.kind not in 'iuf' or k.dtype.kind not in 'iuf' or p.ndim!=1 or not len(p)
            or not np.isfinite(p).all() or not np.isfinite(k).all()):
        raise ValueError('finite real momentum vector and angular sectors required')
    try:
        k=np.broadcast_to(k,p.shape)
    except ValueError as error:
        raise ValueError('angular sectors must broadcast to momenta') from error
    if np.any(k==0) or np.any(k!=np.round(k)):
        raise ValueError('nonzero integer angular Dirac sectors required')
    fields=radar_metric_jets(chart,U,V,region=region)
    raw=reference_projector(fields['lapse'],fields['shift'],fields['radial_length'],
        fields['radius'],p,np.zeros_like(p),k)
    orders=PAULI_CYCLE@raw['orders']@PAULI_CYCLE.conj().T
    return dict(orders=orders,residuals=raw['residuals'],gap=raw['gap'],
        metric_jets=fields,formal_order=4,
        scope='local one-sided subtraction symbol; no state or UV remainder supplied')


def second_order_symbol_jump(momentum,potential,potential_xx_jump):
    """Formal c2 jump with all lower and tangential wall jets continuous.

    H=p*sigma2+W*sigma1, J=[W_xx]. This is the jump of the LOCAL
    subtraction symbol, not an imposed jump of a physical covariance.
    """
    raw=[np.asarray(v) for v in (momentum,potential,potential_xx_jump)]
    if any(v.dtype.kind not in 'iuf' for v in raw):
        raise ValueError('finite real momentum and potential jets required')
    p,W,J=np.broadcast_arrays(*raw)
    omega=np.hypot(p,W)
    if (not all(np.isfinite(v).all() for v in (p,W,J)) or np.any(omega==0)):
        raise ValueError('finite jets and a nonzero energy gap required')
    first=5*W*W*p*p
    second=W*p*(2*p*p-3*W*W)
    return (-J/(16*omega**7))[...,None,None]*(
        first[...,None,None]*SIGMA[0]+second[...,None,None]*SIGMA[1])
