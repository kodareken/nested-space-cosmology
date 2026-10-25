"""Finite radial and quadrupole response on the retained evolving envelope.

The spherical source is replayed unchanged. These are real-linear tangent
equations on that time-dependent source, not stationary eigenvalues or a
nonlinear three-dimensional stability theorem.
"""
from dataclasses import dataclass,asdict
import hashlib
import math
import time

import numpy as np
from numba import njit

from .nsc_atlas_hierarchy_evolution import (
    background,initial_envelope,potential,kinetic_step)
from .nsc_atlas_scale_hierarchy import ScaleHierarchy


def multipole_potential(psi,q,r,dr,A,ell):
    """q=r*eta, Delta_l v=Re(psi* eta)/A. One column per real direction."""
    density=2*np.real(np.conj(psi[:,None])*(q/r[:,None]))
    lower=dr*np.cumsum(r[:,None]**(ell+2)*density,axis=0)
    weighted=r[:,None]**(1-ell)*density
    upper=dr*(np.sum(weighted,axis=0)[None,:]-np.cumsum(weighted,axis=0))
    return -(lower/r[:,None]**(ell+1)+r[:,None]**ell*upper)/(2*A*(2*ell+1))


@njit(cache=True)
def angular_kinetic(q,dr,integrated_inverse_scale_squared,ell):
    """Unitary CN for -1/2 d_rr + ell(ell+1)/(2r^2), q(0)=q(L)=0."""
    n,columns=q.shape;h=integrated_inverse_scale_squared
    off=-1j*h/(4*dr*dr);upper=np.empty(n,dtype=np.complex128)
    result=np.empty_like(q)
    for i in range(n):
        radius=(i+1)*dr
        diagonal=1/(dr*dr)+ell*(ell+1)/(2*radius*radius)
        left=1+.5j*h*diagonal
        if i>0:left-=off*upper[i-1]
        upper[i]=off/left
        for j in range(columns):
            rhs=(1-.5j*h*diagonal)*q[i,j]
            if i>0:rhs-=off*q[i-1,j]+off*result[i-1,j]
            if i<n-1:rhs-=off*q[i+1,j]
            result[i,j]=rhs/left
    for i in range(n-2,-1,-1):
        for j in range(columns):result[i,j]-=upper[i]*result[i+1,j]
    return result


def seed_basis(u,r,dr,A,ambient,ell):
    """Two compact smooth shapes and their phase partners, real Gram=I.

    Radial seeds preserve total envelope mass to first order. Quadrupoles
    do so automatically by angular orthogonality. Angular functions have
    mean square one: Y_l0=sqrt(2l+1)*P_l.
    """
    rho=A*r;z=np.clip((rho-15)/5,0,1);window=1-(6*z**5-15*z**4+10*z**3)
    profiles=[r*rho**ell*np.exp(-rho*rho/(2*w*w))*window for w in (2.,4.)]
    q=np.column_stack((profiles[0],profiles[1],1j*profiles[0],1j*profiles[1]))
    psi=ambient+u/r
    if ell==0:
        direction=r*psi*window
        denominator=dr*np.real(np.vdot(r*psi,direction))
        projections=dr*np.real((r*np.conj(psi))@q)/denominator
        q-=direction[:,None]*projections[None,:]
    gram=dr*np.real(q.conj().T@q);values,vectors=np.linalg.eigh(gram)
    if values[0]<=1e-10*values[-1]:raise RuntimeError('declared initial radial span is ill-conditioned')
    q=q@(vectors/np.sqrt(values)[None,:])
    core_norm=math.sqrt(float(dr*np.sum(abs(u*window)**2)))
    return q,core_norm,float(values[-1]/values[0])


def diagnostics(q,u,r,dr,A,ambient,ell,core_norm,initial_symplectic):
    psi=ambient+u/r;gram=dr*np.real(q.conj().T@q)
    density=2*np.real(np.conj(psi[:,None])*q/r[:,None])
    inside=A*r<=20;radial=r[inside,None]*density[inside]
    contrast=abs(u/r)**2+2*ambient*np.real(u/r)
    density_norm=math.sqrt(float(dr*np.sum((r[inside]*contrast[inside])**2)))
    density_gram=dr*radial.T@radial
    eigen,coeff=np.linalg.eigh(density_gram)
    field_values,field_coefficients=np.linalg.eigh(gram)
    field_direction=field_coefficients[:,-1]
    density_gain=core_norm*math.sqrt(max(0.,float(eigen[-1])))/density_norm
    symmetry=dr*np.imag(q.conj().T@q)
    moments=dr*np.sum(r[inside,None]**4*density[inside],axis=0)
    denominator=dr*np.sum(r[inside]**4*contrast[inside])
    quad=core_norm*np.linalg.norm(moments)/(math.sqrt(5)*abs(denominator)) if ell==2 else None
    return dict(field_norm_gain=math.sqrt(max(0.,float(field_values[-1]))),
        core_density_gain=density_gain,quadrupole_moment_gain=None if quad is None else float(quad),
        density_gain_of_largest_field_direction=core_norm*math.sqrt(max(0.,float(field_direction@density_gram@field_direction)))/density_norm,
        symplectic_gap=float(np.max(abs(symmetry-initial_symplectic))),
        total_mass_first_variation_norm=float(core_norm*dr*np.linalg.norm(np.sum(r[:,None]**2*density,axis=0))) if ell==0 else 0.,
        worst_density_coefficients=coeff[:,-1].tolist())


@dataclass(frozen=True)
class PerturbationRun:
    outer: float=1024.
    points: int=8191
    step: float=.02
    end: float=600.
    late_seed: float=480.


def run_perturbations(config=PerturbationRun(),*,retain=False):
    c=config;cpu=time.process_time();wall=time.perf_counter()
    if not (c.points>=1023 and c.points<=32767 and 128<=c.outer<=2048 and 0<c.step<=.04 and 0<c.end<=600):
        raise ValueError('outside finite perturbation domain')
    steps=round(c.end/c.step);late_step=round(c.late_seed/c.step)
    if abs(steps*c.step-c.end)>1e-10 or abs(late_step*c.step-c.late_seed)>1e-10:raise ValueError('unaligned events')
    model=ScaleHierarchy(1e9,.001,.5);dr=c.outer/(c.points+1)
    r=np.arange(1,c.points+1)*dr;k=np.pi*np.arange(1,c.points+1)/c.outer
    delta,b,prep=initial_envelope(model,r);u=r*delta
    parent=background(model,c.end);times=np.arange(steps+1)*c.step;A,Q=parent.sol(times)
    Am=parent.sol(times[:-1]+c.step/2)[0]
    groups={};initial_bases={};histories={};final_states={}
    for name,ell in (('radial',0),('quadrupole',2),('passive_quadrupole',2)):
        q,norm,condition=seed_basis(u,r,dr,1.,b,ell)
        if name=='passive_quadrupole':q=q[:,:1]
        groups[name]=dict(q=q,ell=ell,core_norm=norm,condition=condition,
            symplectic=dr*np.imag(q.conj().T@q),start=0.)
        initial_bases[name]=q.copy();histories[name]=[]
    initial_u=u.copy();stride=max(1,round(1/c.step))
    for iteration in range(steps):
        if iteration==late_step:
            q,norm,condition=seed_basis(u,r,dr,A[iteration],b,2)
            groups['late_quadrupole']=dict(q=q,ell=2,core_norm=norm,condition=condition,
                symplectic=dr*np.imag(q.conj().T@q),start=float(times[iteration]))
            initial_bases['late_quadrupole']=q.copy();histories['late_quadrupole']=[]
        for half in range(2):
            V,_=potential(u,r,dr,Am[iteration],b);psi=b+u/r;phase=np.exp(-.5j*c.step*V)
            for name,g in groups.items():
                q=g['q']
                dv=0. if name=='passive_quadrupole' else multipole_potential(psi,q,r,dr,Am[iteration],g['ell'])
                g['q']=phase[:,None]*(q-.5j*c.step*(r*psi)[:,None]*dv)
            u=(u+b*r)*phase-b*r
            if half==0:
                interval=Q[iteration+1]-Q[iteration]
                u=kinetic_step(u,k,interval)
                for g in groups.values():
                    if g['ell']==0:g['q']=kinetic_step(g['q'].T,k,interval).T
                    else:g['q']=angular_kinetic(g['q'],dr,interval,g['ell'])
        if (iteration+1)%stride==0 or iteration==steps-1:
            for name,g in groups.items():
                row=diagnostics(g['q'],u,r,dr,A[iteration+1],b,g['ell'],g['core_norm'],g['symplectic'])
                row['tau']=float(times[iteration+1]);histories[name].append(row)
    result={}
    for name,g in groups.items():
        rows=histories[name]
        result[name]=dict(ell=g['ell'],start=g['start'],end=c.end,initial_basis_condition=g['condition'],
            initial_core_field_norm=g['core_norm'],final=rows[-1],
            maximum_field_norm_gain=max(x['field_norm_gain'] for x in rows),
            maximum_core_density_gain=max(x['core_density_gain'] for x in rows),
            maximum_symplectic_gap=max(x['symplectic_gap'] for x in rows),
            maximum_quadrupole_moment_gain=max(x['quadrupole_moment_gain'] for x in rows) if g['ell']==2 else None)
        final_states[name]=g['q']
    digest=hashlib.sha256(np.asarray(np.column_stack((r,u.real,u.imag)),dtype='<f8').tobytes()).hexdigest()
    record=dict(config=asdict(c),responses=result,source_final_field_sha256=digest,
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    if retain:return record,dict(radius=r,initial_u=initial_u,final_u=u,initial_bases=initial_bases,final_states=final_states,histories=histories)
    return record


def nonlinear_radial_check(config,field_file,coefficients,epsilon):
    """Independent centred difference of the full nonlinear radial flow."""
    c=config;cpu=time.process_time();wall=time.perf_counter();data=np.load(field_file)
    r=data['radius'];dr=c.outer/(c.points+1);k=np.pi*np.arange(1,c.points+1)/c.outer
    model=ScaleHierarchy(1e9,.001,.5);_,b,_=initial_envelope(model,r)
    base=data['initial_u'];seed=data['radial_seed']@np.asarray(coefficients)
    reference=data['radial_final']@np.asarray(coefficients)
    # Same normalization used by the response matrices.
    _,core_norm,_=seed_basis(base,r,dr,1.,b,0)
    fields=np.column_stack((base+epsilon*core_norm*seed,base-epsilon*core_norm*seed))
    steps=round(c.end/c.step);times=np.arange(steps+1)*c.step;parent=background(model,c.end)
    A,Q=parent.sol(times);Am=parent.sol(times[:-1]+c.step/2)[0]
    for it in range(steps):
        for half in range(2):
            for j in range(2):
                V,_=potential(fields[:,j],r,dr,Am[it],b)
                fields[:,j]=(fields[:,j]+b*r)*np.exp(-.5j*c.step*V)-b*r
            if half==0:fields=kinetic_step(fields.T,k,Q[it+1]-Q[it]).T
    finite=(fields[:,0]-fields[:,1])/(2*epsilon*core_norm)
    relative=np.linalg.norm(finite-reference)/np.linalg.norm(reference)
    psi=b+data['final_u']/r
    linear_density=2*np.real(np.conj(psi)*reference/r)
    finite_density=(abs(b+fields[:,0]/r)**2-abs(b+fields[:,1]/r)**2)/(2*epsilon*core_norm)
    inside=A[-1]*r<=20
    density_error=np.linalg.norm(r[inside]*(finite_density-linear_density)[inside])/np.linalg.norm(r[inside]*linear_density[inside])
    return dict(epsilon=epsilon,coefficients=list(map(float,coefficients)),relative_field_derivative_error=float(relative),
        relative_core_density_derivative_error=float(density_error),cpu_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall)
