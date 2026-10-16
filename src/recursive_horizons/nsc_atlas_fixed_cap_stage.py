"""Common fixed-cap, material-wall and exterior departure stage.

All field equations remain in the existing cap, hybrid and exterior owners.
This joins one simultaneous state and keeps its independent wall residuals.
No time step or stability claim is made for the clipped Galerkin operator.
"""
from collections.abc import Mapping
import math

import numpy as np

from .nsc_atlas_cap_fixed_basis import FixedCapBasis
from .nsc_atlas_cap_fixed_hybrid import fixed_cap_hybrid_reference
from .nsc_atlas_coupled_evolution import _finite, _grid
from .nsc_atlas_compatible_evolution import FIELDS, WALL_FIELDS
from .nsc_atlas_departure_gradient import quadratic_departure_gradient
from .nsc_atlas_polar_cap import _scalar
from .nsc_atlas_stage_ssprk3 import _shared_budget
from .nsc_atlas_vacuum_departure import vacuum_slice, departure_slice


def fixed_state(value,exterior_cells,basis):
    """Copy and validate coefficients without calling them cell integrals."""
    if not isinstance(basis,FixedCapBasis): raise ValueError('require a fixed cap basis')
    if type(exterior_cells) is not int or not 2<=exterior_cells<=512:
        raise ValueError('require 2..512 exterior cells')
    if not isinstance(value,Mapping) or set(value)!={'departure','wall','cap'}:
        raise ValueError('fixed state must contain departure, wall and cap')
    result={}
    for group,names,shape in (
            ('departure',FIELDS,(exterior_cells-1,)),('wall',WALL_FIELDS,()),
            ('cap',('chi','canonical_density'),(basis.dimension,))):
        if not isinstance(value[group],Mapping) or set(value[group])!=set(names):
            raise ValueError('invalid fixed '+group+' fields')
        result[group]={key:_finite(value[group][key],group+'.'+key,shape) for key in names}
    return result


def quiet_fixed_state(basis,*,exterior_cells):
    """Zero coefficients and vacuum departures; any background is separate."""
    if not isinstance(basis,FixedCapBasis): raise ValueError('require a fixed cap basis')
    if type(exterior_cells) is not int or not 2<=exterior_cells<=512:
        raise ValueError('require 2..512 exterior cells')
    return dict(departure={k:np.zeros(exterior_cells-1) for k in FIELDS},
                wall={k:0. for k in WALL_FIELDS},
                cap={k:np.zeros(basis.dimension) for k in ('chi','canonical_density')})


def exterior_departure_stage(departure,wall,*,x,eta,source,W,W_x,alpha,tube,budget,
                             departure_gradient=quadratic_departure_gradient,
                             quadrature_order=2):
    """Existing endpoint reconstruction, one gradient and ordered prefixes.

    Endpoint data are supplied by the common stage. The returned outgoing
    jets are numerical prefix outputs, never replaced by required wall jets.
    """
    _shared_budget(budget);before=budget.evaluations
    eta=_grid(eta,'eta',2,512)
    x=_scalar(x,'x');W=_scalar(W,'W',positive=True);Wx=_scalar(W_x,'W_x',positive=True)
    alpha=_scalar(alpha,'alpha',nonnegative=True)
    if not callable(departure_gradient): raise ValueError('departure gradient must be callable')
    if not isinstance(departure,Mapping) or set(departure)!=set(FIELDS):
        raise ValueError('exterior departures must contain r,L,c')
    delta={k:_finite(departure[k],'departure.'+k,(len(eta)-2,)) for k in FIELDS}
    wall={k:_scalar(wall[k],'wall.'+k) for k in FIELDS}
    vacuum=vacuum_slice(eta*W,x,budget=budget)
    src={k:_scalar(source['log_C' if k=='L' else k],'source.'+k) for k in FIELDS}
    base={k:vacuum[k] if k!='c' else np.zeros_like(eta) for k in FIELDS}
    offsets={k:src[k]-base[k][0] for k in FIELDS}
    departures={k:np.r_[0.,delta[k],wall[k]-base[k][-1]-offsets[k]] for k in FIELDS}
    X={k:base[k]+offsets[k]+departures[k] for k in FIELDS}
    for k in FIELDS: X[k][0],X[k][-1]=src[k],wall[k]
    gradient_eta,values=eta.copy(),np.stack([departures[k] for k in FIELDS])
    gradient_eta.flags.writeable=values.flags.writeable=False
    G=_finite(departure_gradient(gradient_eta,values,budget=budget),'departure gradient',values.shape)
    gradients={k:G[i] for i,k in enumerate(FIELDS)}
    jets=dict(r_Ubar=gradients['r']/W,log_C_Ubar=gradients['L']/W,c_Ubar=gradients['c']/W)
    ext=departure_slice(eta*W,X,jets,source,x=x,alpha=alpha,tube=tube,budget=budget,
                        quadrature_order=quadrature_order)
    longitudinal={k:ext[name] for k,name in (('r','r_x'),('L','log_C_x'),('c','c_x'))}
    first={k:ext[name] for k,name in (('r','r_Ubar'),('L','log_C_Ubar'),('c','c_Ubar'))}
    delta_q={k:ext[name] for k,name in (('r','delta_r_x'),('L','delta_log_C_x'),('c','delta_c_x'))}
    source_delta=dict(r=_scalar(source['r_x'],'source.r_x')-1.,
                      L=_scalar(source['log_C_x'],'source.log_C_x'),
                      c=_scalar(source['c_x'],'source.c_x'))
    rates={k:(delta_q[k]-source_delta[k]+eta*(Wx/W)*gradients[k])[1:-1] for k in FIELDS}
    outgoing=dict(r_V=float(longitudinal['r'][-1]),c_V=float(longitudinal['c'][-1]),
                   log_C_V=float(longitudinal['L'][-1]))
    budget.check(3*len(eta))
    return dict(exterior=X,base=base,source_values=src,gradients=gradients,prefix=ext,
                longitudinal=longitudinal,first_jets=first,outgoing=outgoing,rates=rates,
                physical_increment={k:X[k][1:-1]-src[k] for k in FIELDS},
                charged_evaluations=budget.evaluations-before)


def coupled_fixed_stage(state,*,basis,x,eta,source,gauge,alpha,tube,budget,
                         departure_gradient=quadratic_departure_gradient,
                         exterior_quadrature_order=2,rtol=2e-13,atol=2e-15,
                         background=None):
    """One actual exterior/cap/hybrid state with separately measured defects.

    With a retained history, cap coefficients describe departures. Evaluate
    that history at this state's physical cap clock, including its evolved
    departure; the coupled wall supplies the clock rate and total trace.
    """
    _shared_budget(budget);before=budget.evaluations
    eta=_grid(eta,'eta',2,512);state=fixed_state(state,len(eta)-1,basis)
    x=_scalar(x,'x');alpha=_scalar(alpha,'alpha',nonnegative=True)
    if not callable(departure_gradient): raise ValueError('departure gradient must be callable')
    W,Wx,boost,boost_x=(_scalar(gauge[k],'gauge.'+k) for k in ('W','W_x','eta_plus','eta_plus_x'))
    if W<=0 or Wx<=0 or not math.isclose(boost,-.5*math.log(Wx),rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError('require positive W,W_x and their compatible boost')
    reference={k:_scalar(gauge[k],'gauge.'+k) for k in WALL_FIELDS}
    reference_rates={k:_scalar(gauge[name],'gauge.'+name) for k,name in
                     (('r','r_x_total'),('L','L_x_total'),('T','T_x'),('s','s_x'))}
    wall={k:reference[k]+float(state['wall'][k]) for k in WALL_FIELDS}
    radius=wall['r']
    if not 0<radius<1: raise ValueError('fixed cap requires 0<r<1')
    Z=float(np.arctanh(radius));exterior={}
    # Source curvature is independent input, not obtained from a mixed RHS.
    source_constraint=None
    if 'r_xx' in source:
        r,c=(_scalar(source[k],'source.'+k) for k in ('r','c'))
        if r<=0: raise ValueError('source radius must be positive')
        rx,cx,Lx,rxx=(_scalar(source[k],'source.'+k) for k in ('r_x','c_x','log_C_x','r_xx'))
        ux=(cx-c*rx/r)/r
        source_constraint=rxx-Lx*rx+alpha*r*ux*ux/(2*(1+alpha*(c/r)**2/6))
        if not math.isfinite(source_constraint): raise ValueError('nonfinite source constraint')
    def outgoing(common):
        difference=common['r']-radius
        if abs(difference)>4*np.finfo(float).eps*radius:
            raise ValueError('unresolved cap radius coordinate conversion')
        wall.update(common)
        result=exterior_departure_stage(state['departure'],wall,x=x,eta=eta,source=source,
            W=W,W_x=Wx,alpha=alpha,tube=tube,budget=budget,departure_gradient=departure_gradient,
            quadrature_order=exterior_quadrature_order)
        exterior.update(result)
        return result['outgoing']
    connection=fixed_cap_hybrid_reference(basis,state['cap']['chi'],state['cap']['canonical_density'],Z,
        alpha=alpha,exterior_wall=outgoing,C_wall=math.exp(wall['L']),eta_plus=boost,eta_plus_x=boost_x,
        budget=budget,rtol=rtol,atol=atol,background=background,
        background_time=wall['T'] if background is not None else None)
    hybrid_rates=connection['wall_rates'];actual_rates=dict(hybrid_rates,c=connection['cap_c_x'])
    q,p=exterior['longitudinal'],exterior['first_jets']
    hybrid_defect={k:hybrid_rates[k]-q[k][-1]-Wx*p[k][-1] for k in FIELDS}
    actual_defect={k:actual_rates[k]-q[k][-1]-Wx*p[k][-1] for k in FIELDS}
    hybrid=connection['hybrid']
    required={k:float(hybrid[name]) for k,name in
              (('r','r_U_plus'),('L','log_C_U_plus'),('c','c_U_plus'))}
    jet_defect={k:float(p[k][-1])-required[k] for k in FIELDS}
    rates=dict(departure=exterior['rates'],
               wall={k:hybrid_rates[k]-reference_rates[k] for k in WALL_FIELDS},
               cap=connection['cap']['rates'])
    if not all(np.all(np.isfinite(v)) for group in rates.values() for v in group.values()):
        raise ValueError('nonfinite fixed common stage rate')
    budget.check(256)
    return dict(rates=rates,physical_wall=wall,wall_rates=hybrid_rates,shared_wall_rates=actual_rates,
        physical_increment=exterior['physical_increment'],exterior=exterior['exterior'],
        exterior_first_jets=p,exterior_longitudinal=q,exterior_prefix=exterior['prefix'],
        departure_gradients=exterior['gradients'],connection=connection,
        wall_kinematic_defect=hybrid_defect,shared_wall_update_defect=actual_defect,
        wall_incoming_jet_defect=jet_defect,source_constraint=source_constraint,
        coordinate_radius_difference=wall['r']-radius,
        exterior_characteristic_courant_rate=Wx/(W*(eta[1]-eta[0])),
        charged_evaluations=budget.evaluations-before,
        scope='one common finite-dimensional stage; no Galerkin time stability or pure bulk constraint certificate')
