"""Whole-line retarded Dirac response of the eight PG packet observables.

The finite cuts only evaluate the global ingoing/outgoing solutions. They
impose no reflecting wall. This avoids subtraction of large horizon-basis
terms in a high-angular subgap spectral density.
"""
from dataclasses import dataclass
from math import exp,log,sqrt

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_complex_horizon_modes import complex_outgoing_ratio,finite_pg_horizon_basis
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_lorentzian import geometry
from .nsc_pg_massive_modes import MassivePGModeResolution
from .nsc_pg_packet_modes import project_horizon_bases
from .nsc_pg_lll_preparation import radial_packet
from .nsc_transmitting_dirac_domain import I2,S2


def finite_cut_outgoing_ratio(z,mass,angular,background,cut=4.,*,rtol=2e-13,atol=2e-15):
    """Restrict the same outgoing Jost solution to a finite exterior cut.

    Existing R_h data do not stably determine this ratio by subtracting two
    large basis columns. Only the previously unretained infinity-to-cut
    solution is evaluated; no horizon/state generator is run.
    """
    z=complex(z)
    if z.real<=0 or z.imag<=0 or cut<=background.horizon_rho:
        raise ValueError('upper-half-plane frequency and exterior evaluation cut required')
    momentum=np.sqrt(z*z-mass*mass)
    end=max(60.,30*max(1.,abs(angular),mass)/max(abs(momentum),.2),2*cut)
    for _ in range(16):
        initial,derivative,_,_=complex_outgoing_ratio(z,mass,angular,end)
        A=background.A_from_offset(end-background.horizon_rho)
        v1=angular*sqrt(A)/sqrt(1+end*end);v2=mass*sqrt(A)
        residual=abs(A*derivative-1j*(v1+1j*v2)+2j*z*initial-1j*(v1-1j*v2)*initial*initial)
        if residual<3e-13:break
        end*=2
    else:raise ArithmeticError('finite-cut Jost series unresolved')
    def rhs(y,state):
        delta=exp(y);rho=background.horizon_rho+delta;A=background.A_from_offset(delta)
        v1=angular*sqrt(A)/sqrt(1+rho*rho);v2=mass*sqrt(A);r=state[0]
        return [delta/A*(1j*(v1+1j*v2)-2j*z*r+1j*(v1-1j*v2)*r*r)]
    run=solve_ivp(rhs,(log(end-background.horizon_rho),log(cut-background.horizon_rho)),
                  np.array([initial]),method='DOP853',rtol=rtol,atol=atol,max_step=.15)
    if not run.success:raise ArithmeticError(run.message)
    ratio=run.y[0,-1]
    return complex(ratio),{'outer_Riccati_residual':float(residual),'finite_cut_ratio_modulus':float(abs(ratio)),
                          'outer_radius':float(end),'function_evaluations':run.nfev}


def _parts(rho,z,mass,angular):
    beta,bp,_=geometry(rho);v=S2-beta*I2
    M=CommonTimeBulkSplit.local_potential(1.,np.sqrt(1+rho*rho),mass,angular)
    iv=np.linalg.inv(v)
    return iv@(1j*(z*I2-M)+.5*bp*I2),1j*iv


def _inner_probes(rho):
    return np.hstack([float(radial_packet(np.array(rho),kind)[0])*I2 for kind in ('parent','child','child_bulk')])


@dataclass
class RetardedPacketResponse:
    z: complex
    response: np.ndarray
    residuals: dict


def retarded_packets(z,mass,angular,preparation,*,basis=None,projected_basis=None,rtol=2e-13,atol=2e-15):
    """J^dagger(H_D-z)^-1 J, with the original fixed PG domain and probes.

    A supplied basis/projection must be the same analytic horizon solution;
    these optional arguments let the covariance integrator reuse its fields.
    There are six trapped observables and two exterior ones. The reverse
    exterior response to trapped sources remains zero by causal support.
    """
    z=complex(z)
    if z.real<=0 or z.imag<=0:raise ValueError('retarded upper-half-plane frequency required')
    if basis is None:basis=finite_pg_horizon_basis(z,mass,angular,preparation,rtol=rtol,atol=atol)
    if projected_basis is None:
        projected_basis,_=project_horizon_bases(z,mass,angular,basis.interior,basis.exterior,rtol=rtol,atol=atol)
    if basis.frequency!=z or basis.mass!=mass or basis.angular!=angular:
        raise ValueError('incoming basis belongs to a different operator/frequency')
    resolution=MassivePGModeResolution(preparation,rtol=rtol,atol=atol)
    ratio,jost_residuals=finite_cut_outgoing_ratio(z,mass,angular,preparation.background,4.,rtol=rtol,atol=atol)
    right=resolution.frame(4.)[0]@np.array([1.,ratio])
    right=right/np.linalg.norm(right)
    incoming=basis.exterior[:,1]
    scale=np.linalg.norm(incoming)
    if scale==0:raise ArithmeticError('ingoing boundary column vanished')
    left=incoming/scale

    # Retarded trapped source solution, reusing the same zero-inflow rule.
    state=np.zeros(12+36,complex)
    for start,end in ((1.,0.),(0.,-1.)):
        def rhs_inner(rho,y):
            A,B=_parts(rho,z,mass,angular);J=_inner_probes(rho);u=y[:12].reshape(2,6)
            return np.r_[(A@u+B@J).ravel(),(-J.conj().T@u).ravel()]
        run=solve_ivp(rhs_inner,(start,end),state,method='DOP853',rtol=rtol,atol=atol)
        if not run.success:raise ArithmeticError(run.message)
        state=run.y[:,-1]
    inner=state[12:].reshape(6,6)

    # Exterior two-sided source problem on [3,4], with global Weyl data.
    initial=np.r_[I2.ravel(),np.zeros(12,complex)]
    def rhs_outer(rho,y):
        A,B=_parts(rho,z,mass,angular);f=float(radial_packet(np.array(rho),'exterior_bulk')[0])
        Y=y[:4].reshape(2,2);part=y[4:8].reshape(2,2)
        return np.r_[(A@Y).ravel(),(A@part+B*f).ravel(),(f*Y).ravel(),(f*part).ravel()]
    run=solve_ivp(rhs_outer,(3.,4.),initial,method='DOP853',rtol=rtol,atol=atol)
    if not run.success:raise ArithmeticError(run.message)
    Y,part,L,integral=run.y[:,-1].reshape(4,2,2)
    forward=Y@left;forward_scale=np.linalg.norm(forward)
    matching=np.column_stack((forward/forward_scale,-right))
    coefficients=np.linalg.solve(matching,-part)
    c=coefficients[0]/forward_scale
    outer=np.outer(L@left,c)+integral
    endpoint=forward[:,None]*c[None,:]+part-right[:,None]*coefficients[1][None,:]
    G=np.zeros((8,8),complex);G[:6,:6]=inner;G[6:,6:]=outer
    G[:6,6:]=np.outer(projected_basis[:6,1]/scale,c)
    imaginary=(G-G.conj().T)/(2j)
    eig=np.linalg.eigvalsh(imaginary)
    residuals={**jost_residuals,
               'endpoint_matching':float(np.linalg.norm(endpoint)),
               'imaginary_part_lower_violation':max(0.,-float(eig.min())),
               'resolvent_norm':float(np.linalg.norm(G,2)),
               'resolvent_bound_violation':max(0.,float(np.linalg.norm(G,2))-1/z.imag),
               'right_matching_condition':float(np.linalg.cond(matching)),
               'forward_basis_norm':float(forward_scale),
               'reverse_retarded_response':float(np.linalg.norm(G[6:,:6]))}
    return RetardedPacketResponse(z,G,residuals)
