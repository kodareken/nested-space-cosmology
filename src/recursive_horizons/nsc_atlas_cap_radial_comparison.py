"""Compose native and smooth-reference cap maps on one frozen backward ray.

This computes a local map-error sum and its derivative, not a new evolution.
Reference solver tolerances remain distinct from an integration-error bound.
"""
import numpy as np

from .nsc_atlas_background_cap import _budget, _gauss_rule
from .nsc_atlas_cap_mass_tangent import _Jet, _poly
from .nsc_atlas_cap_radial_reference import affine_radial_map
from .nsc_atlas_cap_transport_fields import CapTransportFields
from .nsc_atlas_characteristic_background import SMALL_RADIUS
from .nsc_atlas_polar_cap import _scalar


def _compose(first, second):
    """Second map after first; derivatives include both incoming factors."""
    a,b=first,second
    return dict(S=b['S']*a['S'],T=b['S']*a['T']+b['T'],
        S_x=b['S_x']*a['S']+b['S']*a['S_x'],
        T_x=b['S_x']*a['T']+b['S']*a['T_x']+b['T_x'])


def _reference_pieces(panel):
    """Add the fixed formula switch in backward-ray order, including slivers."""
    left=(panel['left'],panel['left_rate']); right=(panel['right'],panel['right_rate'])
    key=lambda p:(p[0],-p[1])
    if key(right)<key(left):
        raise ValueError('reference panel is reversed on the backward ray')
    fixed=(SMALL_RADIUS,0.)
    points=[left]+([fixed] if key(left)<key(fixed)<key(right) else [])+[right]
    result=[]
    for low,high in zip(points[:-1],points[1:]):
        if low[0]<high[0]:
            side='left' if high[0]<=SMALL_RADIUS else 'right'
            if low[0]<SMALL_RADIUS<high[0]:
                raise ValueError('reference piece still crosses the formula switch')
        else:
            midpoint_rate=(low[1]+high[1])/2
            side=('left' if low[0]<SMALL_RADIUS or
                  low[0]==SMALL_RADIUS and midpoint_rate>=0 else 'right')
        result.append(dict(left=low[0],right=high[0],left_rate=low[1],right_rate=high[1],
                           owner=panel['owner'],small_radius_side=side,born=low[0]==high[0]))
    return result


def _native_map(view,panel,order,budget):
    left=_Jet(panel['left'],np.array([panel['left_rate']]))
    right=_Jet(panel['right'],np.array([panel['right_rate']]))
    width=right-left
    if width.v==0:
        if left.v==0:
            a=b=0.
        else:
            row=view.sample(panel['owner'],[float(left.v)],budget=budget)
            a,b=float(row['a'][0]),float(row['b'][0])
        budget.check(4)
        return dict(S=1.,T=0.,S_x=-float(width.d[0])*a,T_x=float(width.d[0])*b)
    mid,h=(left+right)/2,width/2
    nodes,weights,basis=_gauss_rule(order)
    def at(z,names):
        row=view.sample(panel['owner'],z.v,budget=budget)
        return [_Jet(row[name],row[name+'_x'][:,None]+row[name+'_r'][:,None]*z.d) for name in names]
    primary=mid+h*nodes
    a,=at(primary,('a',))
    edge=left+h*(1+nodes)
    ap,bp=at(edge,('a','b'))
    theta=(edge-mid)/h
    Gpoly=sum(h*a[k]*_poly(poly.coef,theta) for k,poly in enumerate(basis))
    G=h*(ap*weights).sum()
    S=(-G).exp(); T=h*((Gpoly-G).exp()*bp*weights).sum()
    budget.check(8*order+4)
    return dict(S=float(S.v),T=float(T.v),S_x=float(S.d[0]),T_x=float(T.d[0]))


def _reference_map(view,panel,budget,rtol,atol):
    result=dict(S=1.,T=0.,S_x=0.,T_x=0.)
    pieces=_reference_pieces(panel); calls=0
    for piece in pieces:
        def coefficients(r):
            if r==0:
                return view.reference_axis_coefficients(budget=budget)
            side=piece['small_radius_side'] if r==SMALL_RADIUS else None
            row=view.reference_sample(piece['owner'],[r],budget=budget,small_radius_side=side)
            return {key:float(row[key][0]) for key in ('a','b','a_x','b_x','a_r','b_r')}
        flow=affine_radial_map(coefficients,**{key:piece[key] for key in
            ('left','right','left_rate','right_rate')},budget=budget,rtol=rtol,atol=atol)
        calls+=flow['coefficient_calls']
        result=_compose(result,flow)
        budget.check(8)
    return result,pieces,calls


def compare_cap_radial_maps(metric,direction,tangent,*,budget,rtol=2e-11,atol=2e-13):
    """Return the independently weighted local error and derivative.

    The supplied tangent is the existing native mass_tangent result on this
    same direction. Native face values, panel ownership and the full tangent
    are checked before accepting the comparison. No total-derivative
    difference is used to construct the local-error derivative.

    The native tangent's domain still applies: a native quadrature sample or
    born panel exactly on a formula switch is rejected. Reference-only pieces
    born at the fixed radial switch use the explicitly named limiting side.
    """
    _budget(budget,metric); before=budget.evaluations
    native_rate=_scalar(tangent['forward'],'native tangent')
    native_roundoff=_scalar(tangent['contraction_roundoff_allowance'],'native roundoff',nonnegative=True)
    view=CapTransportFields(metric,direction,budget=budget)
    panels=view.native_panels(budget=budget)
    if (view.limiter_branches!={k:tuple(v) for k,v in tangent['limiter_branches'].items()}
            or len(panels)!=len(tangent['panels'])):
        raise ValueError('native tangent partition/limiter differs')
    for p,saved in zip(panels,tangent['panels']):
        if any(p[k]!=saved[k] for k in ('left','right','owner','born')):
            raise ValueError('native tangent owner or endpoint differs')
    order=metric['quadrature_order']; _gauss_rule(order)
    q=qx=e=ex=0.; qabs=qxabs=eabs=exabs=0.
    qfaces=np.zeros(len(metric['r_faces'])); efaces=qfaces.copy()
    qxfaces=qfaces.copy(); exfaces=qfaces.copy()
    rows=[]; reference_calls=0
    for panel in panels:
        Q=_native_map(view,panel,order,budget)
        E,pieces,calls=_reference_map(view,panel,budget,rtol,atol)
        reference_calls+=calls
        if not all(np.isfinite(v) for m in (Q,E) for v in m.values()):
            raise ValueError('nonfinite radial affine map')
        eps=(Q['S']-E['S'])*q+(Q['T']-E['T'])
        epsx=(Q['S_x']-E['S_x'])*q+(Q['S']-E['S'])*qx+(Q['T_x']-E['T_x'])
        epsabs=(abs(Q['S'])+abs(E['S']))*qabs+abs(Q['T'])+abs(E['T'])
        epsxabs=((abs(Q['S_x'])+abs(E['S_x']))*qabs+
                 (abs(Q['S'])+abs(E['S']))*qxabs+abs(Q['T_x'])+abs(E['T_x']))
        rows.append(dict(panel=panel,native=Q,reference=E,pieces=pieces,
            entering_mass=q,entering_rate=qx,local_error=eps,local_error_rate=epsx,
            absolute_error_work=epsabs,absolute_rate_work=epsxabs))
        q,qx=Q['S']*q+Q['T'],Q['S_x']*q+Q['S']*qx+Q['T_x']
        e,ex=E['S']*e+E['T'],E['S_x']*e+E['S']*ex+E['T_x']
        qabs,qxabs=(abs(Q['S'])*qabs+abs(Q['T']),
            abs(Q['S_x'])*qabs+abs(Q['S'])*qxabs+abs(Q['T_x']))
        eabs,exabs=(abs(E['S'])*eabs+abs(E['T']),
            abs(E['S_x'])*eabs+abs(E['S'])*exabs+abs(E['T_x']))
        if panel['right_tag'][0]=='cell':
            j=panel['right_tag'][1]
            qfaces[j],efaces[j],qxfaces[j],exfaces[j]=q,e,qx,ex
        budget.check(32)
    if not np.allclose(qfaces,metric['n_faces'],rtol=2e-12,atol=2e-13):
        raise ValueError('native affine composition differs from retained primal faces')
    tangent_allowance=max(2e-13+2e-12*abs(native_rate),native_roundoff)
    if abs(qx-native_rate)>tangent_allowance:
        raise ValueError('native affine composition differs from full mass tangent')
    sigma=1.; sigmax=0.; sigmaabs=1.; sigmaxabs=0.
    error=errorx=errorwork=errorxwork=0.
    for row in reversed(rows):
        error+=sigma*row['local_error']
        errorx+=sigmax*row['local_error']+sigma*row['local_error_rate']
        errorwork+=sigmaabs*row['absolute_error_work']
        errorxwork+=sigmaxabs*row['absolute_error_work']+sigmaabs*row['absolute_rate_work']
        S,Sx=row['reference']['S'],row['reference']['S_x']
        sigma,sigmax=sigma*S,sigmax*S+sigma*Sx
        sigmaabs,sigmaxabs=sigmaabs*abs(S),sigmaxabs*abs(S)+sigmaabs*abs(Sx)
        budget.check(16)
    # Arithmetic check conditional on the computed maps. Absolute recurrences
    # retain cancellation work; 32 operations per panel plus 64 terminal
    # operations exceed every expanded product path in the two telescopes.
    # This does not include coefficient evaluation or radial solver error.
    k=32*len(rows)+64; gamma=k*np.finfo(float).eps/(1-k*np.finfo(float).eps)
    bounds=gamma*np.array([qabs+eabs+errorwork,qxabs+exabs+errorxwork])/(1-gamma)
    defects=np.array([error-(q-e),errorx-(qx-ex)])
    if not np.all(np.isfinite(bounds)) or np.any(abs(defects)>bounds):
        raise ValueError('radial local-error telescope failed its arithmetic check')
    budget.check()
    return dict(native_mass=q,native_rate=qx,reference_mass=e,reference_rate=ex,
        error=error,error_rate=errorx,direct_error=q-e,direct_error_rate=qx-ex,
        telescope_defects=defects,telescope_roundoff_bounds=bounds,
        native_tangent_defect=qx-native_rate,native_tangent_allowance=tangent_allowance,
        native_faces=qfaces,reference_faces=efaces,native_rate_faces=qxfaces,reference_rate_faces=exfaces,
        panels=rows,reference_coefficient_calls=reference_calls,evaluations=budget.evaluations-before,
        requested_rtol=rtol,requested_atol=atol,uniform_error_enclosure=False)
