"""Ordered cap reference channels and independent native mass-balance comparison.

No field evolution, source fitting or grid search occurs here. Requested ODE
tolerances and sampled metric margins are not numerical error enclosures.
"""
import math

import numpy as np

from .nsc_atlas_background_cap import _budget,background_cap_profile
from .nsc_atlas_cap_radial_comparison import _reference_pieces,compare_cap_radial_maps
from .nsc_atlas_cap_transport_fields import CapTransportFields
from .nsc_atlas_cap_transport_reference import STATE_KEYS,transport_radial_piece
from .nsc_atlas_characteristic_background import SMALL_RADIUS
from .nsc_atlas_polar_cap import _scalar


def _fields(view,piece,r,budget):
    side=piece['small_radius_side'] if r==SMALL_RADIUS else None
    row=view.reference_sample(piece['owner'],[r],budget=budget,small_radius_side=side)
    return {k:float(np.asarray(v)[0]) for k,v in row.items()}


def _trace(r,n,delta,row,alpha):
    A=1-r*r-2*alpha*n/r
    rho=r*r*(row['Pdiff']/2+A*(row['Pi']**2+row['Phi']**2)/4)
    J=r*r*A*A*math.exp(delta)*row['Pi']*row['Phi']/2
    if not np.all(np.isfinite([A,rho,J])):
        raise ValueError('nonfinite cap reference trace')
    return dict(A=A,rho=rho,J=J)


def _integrate_pieces(view,pieces,*,center_delta,alpha,clock_rate,budget,rtol,atol):
    """Ordered piece chain; only the first/last events are physical endpoints."""
    if not pieces or pieces[0]['left']!=0:
        raise ValueError('cap reference pieces must start at the centre')
    state=dict.fromkeys(STATE_KEYS,0.); state['delta']=center_delta
    jumps=[]; stages=[]; minimum_A=1.; calls=rhs_calls=0
    for i,piece in enumerate(pieces):
        if i:
            previous=pieces[i-1]; r=piece['left']; velocity=piece['left_rate']
            if previous['right']!=r or previous['right_rate']!=velocity:
                raise ValueError('cap reference interface coordinates/rates disagree')
            if r==0:
                if state['n']!=0:
                    raise ValueError('nonregular reference mass at an axis event')
                rho_jump=flux_jump=0.
            else:
                left=_fields(view,previous,r,budget); right=_fields(view,piece,r,budget)
                A=1-r*r-2*alpha*state['n']/r
                minimum_A=min(minimum_A,A)
                energy_jump=sum((right[k]-left[k])*(right[k]+left[k]) for k in ('Pi','Phi'))
                rho_jump=r*r*((right['Pdiff']-left['Pdiff'])/2+A*energy_jump/4)
                product_jump=(right['Pi']-left['Pi'])*right['Phi']+left['Pi']*(right['Phi']-left['Phi'])
                flux_jump=r*r*A*A*math.exp(state['delta'])*product_jump/2
            mass_motion=-velocity*rho_jump; clock_flux=-clock_rate*flux_jump
            state['Z_jump']+=mass_motion+clock_flux
            if not np.all(np.isfinite([rho_jump,flux_jump,state['Z_jump']])):
                raise ValueError('nonfinite cap reference jump')
            jumps.append(dict(r=r,velocity=velocity,left_owner=previous['owner'],right_owner=piece['owner'],
                left_side=previous['small_radius_side'],right_side=piece['small_radius_side'],
                rho_jump=rho_jump,flux_jump=flux_jump,mass_motion=mass_motion,clock_flux=clock_flux,
                delta=state['delta']))
            budget.check(16)
        result=transport_radial_piece(lambda r:_fields(view,piece,r,budget),
            left=piece['left'],right=piece['right'],initial=state,alpha=alpha,
            clock_rate=clock_rate,budget=budget,rtol=rtol,atol=atol)
        state={k:result[k] for k in STATE_KEYS}
        minimum_A=min(minimum_A,result['minimum_sampled_A'])
        calls+=result['coefficient_calls']; rhs_calls+=result['rhs_calls']
        stages.append(dict(piece=piece,state=state.copy()))
    wall_r=pieces[-1]['right']
    if wall_r<=0:
        raise ValueError('cap reference requires a positive wall radius')
    wall_fields=_fields(view,pieces[-1],wall_r,budget)
    wall=_trace(wall_r,state['n'],state['delta'],wall_fields,alpha)
    minimum_A=min(minimum_A,wall['A'])
    for jump in jumps:
        jump['attenuation']=math.exp(jump['delta']-state['delta'])
    jump_sum=math.fsum(j['attenuation']*(j['mass_motion']+j['clock_flux']) for j in jumps)
    budget.check(8+3*len(jumps))
    return dict(**state,wall=wall,wall_fields=wall_fields,jumps=jumps,stages=stages,
        jump_sum_from_lapse=jump_sum,jump_sum_defect=state['Z_jump']-jump_sum,
        minimum_sampled_A=minimum_A,sampled_A_positive=minimum_A>0,
        coefficient_calls=calls,rhs_calls=rhs_calls,uniform_error_enclosure=False)


def reference_transport_cap(metric,direction,*,budget,rtol=2e-11,atol=2e-13):
    """Integrate all reference channels without reading a native mass tangent."""
    _budget(budget,metric); before=budget.evaluations
    view=CapTransportFields(metric,direction,budget=budget)
    center=view.reference_center(budget=budget)
    if (not all(np.isclose(center[k],metric['center'][k],rtol=2e-12,atol=2e-13) for k in center)
            or not np.isclose(center['delta'],metric['delta_faces'][0],rtol=2e-12,atol=2e-13)):
        raise ValueError('cap reference centre differs from retained field/lapse')
    panels=view.native_panels(budget=budget)
    pieces=[piece for panel in panels for piece in _reference_pieces(panel)]
    result=_integrate_pieces(view,pieces,center_delta=float(metric['center']['delta']),
        alpha=float(metric['alpha']),clock_rate=float(direction['T']),budget=budget,rtol=rtol,atol=atol)
    result.update(panels=panels,center=center,requested_rtol=rtol,requested_atol=atol,
                  evaluations=budget.evaluations-before)
    return result


def cap_transport_balance(metric,direction,tangent,*,budget,rtol=2e-11,atol=2e-13):
    """Compare local-map, wall, bulk and interface channels to native B_mass.

    A returned defect is a numerical comparison, not automatically a pass.
    Its acceptance requires the caller's independent error estimate.
    The native wall sample keeps its existing formula-switch exclusions.
    """
    _budget(budget,metric); before=budget.evaluations
    R,alpha=float(metric['radius']),float(metric['alpha'])
    if R==SMALL_RADIUS:
        raise ValueError('full cap balance excludes a wall on the small-radius formula switch')
    reverse=_scalar(tangent['reverse'],'native reverse tangent')
    forward=_scalar(tangent['forward'],'native forward tangent')
    allowance=_scalar(tangent['contraction_roundoff_allowance'],'native contraction allowance',nonnegative=True)
    if abs(reverse-forward)>allowance:
        raise ValueError('native forward/reverse tangent contractions differ')
    view=CapTransportFields(metric,direction,budget=budget)
    owner=len(metric['r_faces'])-2
    native={k:float(np.asarray(v)[0]) for k,v in view.sample(owner,[R],budget=budget).items()}
    maps=compare_cap_radial_maps(metric,direction,tangent,budget=budget,rtol=rtol,atol=atol)
    reference=reference_transport_cap(metric,direction,budget=budget,rtol=rtol,atol=atol)
    if reference['panels']!=[row['panel'] for row in maps['panels']]:
        raise ValueError('cap map and transport partitions differ')
    rx,tx=(float(direction[k]) for k in ('R','T'))
    original=background_cap_profile(metric,np.array([R]),budget=budget)
    if (not all(np.allclose(native[k],original[k],rtol=2e-12,atol=2e-13) for k in ('f','Pi','Phi'))
            or not np.isclose(native['f'],metric['wall']['f'],rtol=2e-12,atol=2e-13)):
        raise ValueError('native mass-balance wall trace differs')
    Nq=float(metric['n_faces'][-1]); deltaq=_scalar(metric['wall']['delta'],'native wall lapse')
    qwall=_trace(R,Nq,deltaq,native,alpha)
    if not np.isclose(qwall['A'],metric['wall']['A'],rtol=2e-12,atol=2e-13):
        raise ValueError('native wall metric differs from radial mass')
    ew=reference['wall_fields']
    if any(native[k]!=ew[k] for k in ('f','Pi','Phi','a')):
        raise ValueError('native and reference wall fields require a common limiting trace')
    # Stable analytic potential gap; shared wall fields make a_E == a_Q.
    density_correction=ew['a']*maps['error']+R*R*ew['implementation_gap']/2
    channels=dict(radial_map_error_rate=maps['error_rate'],
        wall_density_conversion=rx*density_correction,
        wall_flux_conversion=tx*(reference['wall']['J']-qwall['J']),
        field_f=reference['Z_f'],field_Pi=reference['Z_Pi'],field_Phi=reference['Z_Phi'],
        interfaces=reference['Z_jump'])
    total=math.fsum(channels.values())
    native_balance=reverse-rx*qwall['rho']-tx*qwall['J']
    reference_rate=math.fsum(reference[k] for k in ('Z_f','Z_Pi','Z_Phi','Z_jump'))
    reference_rate+=rx*reference['wall']['rho']+tx*reference['wall']['J']
    if not np.all(np.isfinite([*channels.values(),total,native_balance,reference_rate,
                              alpha*native_balance,alpha*(total-native_balance)])):
        raise ValueError('nonfinite full cap balance')
    budget.check(32)
    return dict(channels=channels,total=total,native_balance=native_balance,
        balance_defect=total-native_balance,actual_native_balance=alpha*native_balance,
        actual_balance_defect=alpha*(total-native_balance),native_wall=qwall,
        reference_rate_from_transport=reference_rate,
        reference_rate_defect=reference_rate-maps['reference_rate'],
        reference_mass_defect=reference['n']-maps['reference_mass'],
        maps=maps,reference=reference,evaluations=budget.evaluations-before,
        uniform_error_enclosure=False)
