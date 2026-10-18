"""Clocked causal elimination of geometric variations along the owned ODE.

The full nonlinear reference and original coordinate-time tangent are kept.
Only variations are split by the fixed constant-mode canonical projector;
all fermion variables are retained. This is a finite local-Hamiltonian
response, not the curved real-time completion of the spectral modulus.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nsc_conformal_tangent_evolution import conformal_jvp
from .nsc_spherical_coupling import CauchyState, _combine
from .nsc_spherical_galerkin_coupling import compose_fine_hamiltonian


CANONICAL_GEOMETRY = ('a','b','chi','pi_a','pi_b','p_chi')
CANONICAL_FIELDS = CANONICAL_GEOMETRY + ('phi0','phi1')


@dataclass
class CanonicalVariation:
    a: np.ndarray
    b: np.ndarray
    chi: np.ndarray
    pi_a: np.ndarray
    pi_b: np.ndarray
    p_chi: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray


def combine_variations(*terms):
    """Linear combination without mutating a stage or an input direction."""
    if not terms:
        raise ValueError('at least one canonical variation is required')
    return CanonicalVariation(**{name:sum(
        (weight*getattr(value,name) for weight,value in terms),
        np.zeros_like(getattr(terms[0][1],name))) for name in CANONICAL_FIELDS})


def _positive_coordinates(state):
    if (not np.isfinite(state.r).all() or not np.isfinite(state.Q).all()
            or np.min(state.r)<=0 or np.min(state.Q)<=0):
        raise ValueError('canonical conversion needs positive finite nodal r,Q')


def to_canonical_tangent(state, direction):
    """DS(z)u, with y=(ln r,ln Q,chi,r*p_r,Q*p_Q,p_chi,phi0,phi1)."""
    _positive_coordinates(state)
    return CanonicalVariation(a=direction.r/state.r,b=direction.Q/state.Q,
        chi=direction.chi.copy(),pi_a=state.r*direction.p_r+state.p_r*direction.r,
        pi_b=state.Q*direction.p_Q+state.p_Q*direction.Q,p_chi=direction.p_chi.copy(),
        phi0=direction.phi0.copy(),phi1=direction.phi1.copy())


def to_physical_tangent(state, direction):
    """DS(z)^(-1)v. These arrays are variations, not a new Cauchy state."""
    _positive_coordinates(state)
    return CauchyState(r=state.r*direction.a,Q=state.Q*direction.b,
        chi=direction.chi.copy(),p_r=direction.pi_a/state.r-state.p_r*direction.a,
        p_Q=direction.pi_b/state.Q-state.p_Q*direction.b,p_chi=direction.p_chi.copy(),
        phi0=direction.phi0.copy(),phi1=direction.phi1.copy())


def canonical_jvp(grid,state,direction,*,rate=None,bundle=None,checkpoint=None):
    """D(DS*f)DS^(-1)v, including the second derivative of the coordinate map."""
    if rate is None or bundle is None:
        rate,bundle=compose_fine_hamiltonian(grid,state)
    physical=to_physical_tangent(state,direction)
    differentiated=conformal_jvp(grid,state,physical,bundle=bundle,checkpoint=checkpoint)
    result=to_canonical_tangent(state,differentiated)
    result.a-=rate.r*physical.r/state.r**2
    result.b-=rate.Q*physical.Q/state.Q**2
    result.pi_a+=physical.r*rate.p_r+physical.p_r*rate.r
    result.pi_b+=physical.Q*rate.p_Q+physical.p_Q*rate.Q
    return result


def observer_values(state,row):
    weights=np.asarray(row,dtype=float)
    if weights.shape!=state.r.shape or not np.isfinite(weights).all():
        raise ValueError('observer row must match the geometry nodes')
    radius,scale=float(weights@state.r),float(weights@state.Q)
    if radius<=0 or scale<=0 or not np.isfinite([radius,scale]).all():
        raise ValueError('the prolonged observer must have positive r,Q')
    return radius,scale,radius*scale


def observer_clock_direction(state,physical_direction,row):
    """DN_o u at the prolonged observer, in the original coordinate-time family."""
    radius,scale,lapse=observer_values(state,row)
    return lapse*(float(row@physical_direction.r)/radius+float(row@physical_direction.Q)/scale)


def clocked_canonical_jvp(grid,state,direction,row,*,rate=None,bundle=None,checkpoint=None):
    """[J_y-f_y tensor D ln N_o]v, with derivative taken along baseline T.

    The proper-time variational generator is this operator divided by N_o.
    Its use here follows dtau_o/dT=N_o; it is not the Hamiltonian of H/N_o.
    """
    if rate is None or bundle is None:
        rate,bundle=compose_fine_hamiltonian(grid,state)
    radius,scale,_lapse=observer_values(state,row)
    dlog=(float(row@(state.r*direction.a))/radius
          +float(row@(state.Q*direction.b))/scale)
    tangent_rate=canonical_jvp(grid,state,direction,rate=rate,bundle=bundle,checkpoint=checkpoint)
    return combine_variations((1,tangent_rate),(-dlog,to_canonical_tangent(state,rate)))


def retained_projection(direction):
    """Exact fixed rank-one geometric span; identity on both fermion blocks."""
    values={name:np.full_like(getattr(direction,name),np.mean(getattr(direction,name)))
            for name in CANONICAL_GEOMETRY}
    values.update({name:getattr(direction,name).copy() for name in ('phi0','phi1')})
    return CanonicalVariation(**values)


def environment_projection(direction):
    return combine_variations((1,direction),(-1,retained_projection(direction)))


def causal_split_rates(apply,retained,initial,feedback):
    """Closed Schur split, with two operator actions shared by all three rates."""
    drive_rate=apply(initial)
    coupled_rate=apply(combine_variations((1,feedback),(1,retained)))
    return (retained_projection(combine_variations((1,drive_rate),(1,coupled_rate))),
            environment_projection(drive_rate),environment_projection(coupled_rate))


@dataclass
class ResolutionPacket:
    state: CauchyState
    coordinate_tangent: CauchyState
    retained: CanonicalVariation
    initial_drive: CanonicalVariation
    feedback: CanonicalVariation
    clock: float=0.0
    clock_tangent: float=0.0


def initial_packet(state,tangent):
    canonical=to_canonical_tangent(state,tangent)
    return ResolutionPacket(state.copy(),tangent.copy(),retained_projection(canonical),
        environment_projection(canonical),combine_variations((0,canonical)))


def proper_time_reference(packet,rate,row):
    """Independent fixed-proper-time tangent from original u(T) and delta tau(T)."""
    lapse=observer_values(packet.state,row)[2]
    return combine_variations((1,to_canonical_tangent(packet.state,packet.coordinate_tangent)),
        (-packet.clock_tangent/lapse,to_canonical_tangent(packet.state,rate)))


def reconstructed_tangent(packet):
    return combine_variations((1,packet.retained),(1,packet.initial_drive),(1,packet.feedback))


def fixed_receiver_clock_redshift(dlog_receiver,dlog_emitter,dotlog_receiver,
                                 dotlog_emitter,clock_tangent,receiver_lapse):
    """d log(1+z) at fixed receiver proper time for the existing fixed-x rays.

    In conformal gauge delta T_emitter=delta T_receiver=-delta tau_o/N_o.
    dlog arguments are original fixed-T variations; dotlog are T derivatives.
    """
    values=np.asarray([dlog_receiver,dlog_emitter,dotlog_receiver,dotlog_emitter,
                       clock_tangent,receiver_lapse],dtype=float)
    if not np.isfinite(values).all() or receiver_lapse<=0:
        raise ValueError('fixed receiver clock requires finite data and positive lapse')
    return float(dlog_receiver-dlog_emitter
                 -(dotlog_receiver-dotlog_emitter)*clock_tangent/receiver_lapse)


def resolution_stage(grid,packet,row,*,checkpoint=None):
    rate,bundle=compose_fine_hamiltonian(grid,packet.state)
    original=conformal_jvp(grid,packet.state,packet.coordinate_tangent,
        bundle=bundle,checkpoint=checkpoint)
    def apply(direction):
        return clocked_canonical_jvp(grid,packet.state,direction,row,
            rate=rate,bundle=bundle,checkpoint=checkpoint)
    retained,initial,feedback=causal_split_rates(apply,packet.retained,packet.initial_drive,packet.feedback)
    return ResolutionPacket(rate,original,retained,initial,feedback,
        observer_values(packet.state,row)[2],
        observer_clock_direction(packet.state,packet.coordinate_tangent,row))


def _packet_offset(packet,rate,step):
    return ResolutionPacket(_combine(packet.state,rate.state,step),
        _combine(packet.coordinate_tangent,rate.coordinate_tangent,step),
        combine_variations((1,packet.retained),(step,rate.retained)),
        combine_variations((1,packet.initial_drive),(step,rate.initial_drive)),
        combine_variations((1,packet.feedback),(step,rate.feedback)),
        packet.clock+step*rate.clock,packet.clock_tangent+step*rate.clock_tangent)


def resolution_rk4_step(grid,packet,dt,row,*,checkpoint=None):
    """Advance the original flow, original tangent, clock and causal split together."""
    step=float(dt)
    if not np.isfinite(step) or step<=0:
        raise ValueError('causal resolution requires a finite positive time step')
    rates=[]
    for weight,index in ((0,None),(.5,0),(.5,1),(1,2)):
        if checkpoint is not None:checkpoint()
        stage=packet if index is None else _packet_offset(packet,rates[index],step*weight)
        rates.append(resolution_stage(grid,stage,row,checkpoint=checkpoint))
    weights=(1/6,1/3,1/3,1/6)
    def physical(field):
        states=[getattr(r,field) for r in rates]
        return CauchyState(**{name:sum(w*getattr(s,name) for w,s in zip(weights,states))
                              for name in ('Q','r','chi','p_Q','p_r','p_chi','phi0','phi1')})
    average=ResolutionPacket(physical('state'),physical('coordinate_tangent'),
        *(combine_variations(*[(w,getattr(r,name)) for w,r in zip(weights,rates)])
          for name in ('retained','initial_drive','feedback')),
        sum(w*r.clock for w,r in zip(weights,rates)),
        sum(w*r.clock_tangent for w,r in zip(weights,rates)))
    return _packet_offset(packet,average,step)
