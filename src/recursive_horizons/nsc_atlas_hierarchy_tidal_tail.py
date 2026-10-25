"""Leading tidal-phase matching for the existing quadrupole response.

The background and tangent bulk equations are unchanged. The old finite
Dirichlet experiment remains retained; this separately declared boundary
allows the forced homogeneous-ambient harmonic phase tail.
"""
from dataclasses import asdict
import hashlib
import time
import numpy as np

from .nsc_atlas_hierarchy_evolution import background,initial_envelope,potential,kinetic_step
from .nsc_atlas_scale_hierarchy import ScaleHierarchy
from .nsc_atlas_hierarchy_perturbations import (
    PerturbationRun,angular_kinetic,multipole_potential,seed_basis,diagnostics)


def tidal_kinetic(q,dr,interval,ell=2):
    """Exact rank-one correction of the verified Dirichlet Cayley map.

    q(L)=(r_last/L)^ell q(r_last) continues q proportional to r^-ell.
    This changes one real Hamiltonian diagonal and remains unitary.
    """
    n=q.shape[0];e=np.zeros((n,1),dtype=complex);e[-1,0]=1.
    propagated=angular_kinetic(np.column_stack((q,e)),dr,interval,ell)
    Uq=propagated[:,:-1];inverse_e=(propagated[:,-1]+e[:,0])/2
    gamma=(n/(n+1))**ell
    alpha=-1j*interval*gamma/(4*dr*dr)
    correction=alpha*(q[-1]+Uq[-1])/(1+alpha*inverse_e[-1])
    return Uq-inverse_e[:,None]*correction[None,:]


def run_tidal_tail(config=PerturbationRun(),*,retain=False):
    c=config;cpu=time.process_time();wall=time.perf_counter();model=ScaleHierarchy(1e9,.001,.5)
    dr=c.outer/(c.points+1);r=np.arange(1,c.points+1)*dr;k=np.pi*np.arange(1,c.points+1)/c.outer
    delta,b,prep=initial_envelope(model,r);u=r*delta
    q,Nc,condition=seed_basis(u,r,dr,1.,b,2);S=dr*np.imag(q.conj().T@q)
    steps=round(c.end/c.step);times=np.arange(steps+1)*c.step;parent=background(model,c.end)
    A,Q=parent.sol(times);Am=parent.sol(times[:-1]+c.step/2)[0]
    stride=round(1/c.step);history=[]
    for it in range(steps):
        for half in range(2):
            V,_=potential(u,r,dr,Am[it],b);psi=b+u/r
            dv=multipole_potential(psi,q,r,dr,Am[it],2);phase=np.exp(-.5j*c.step*V)
            q=phase[:,None]*(q-.5j*c.step*(r*psi)[:,None]*dv)
            u=(u+b*r)*phase-b*r
            if half==0:
                interval=Q[it+1]-Q[it]
                u=kinetic_step(u,k,interval);q=tidal_kinetic(q,dr,interval)
        if (it+1)%stride==0 or it==steps-1:
            row=diagnostics(q,u,r,dr,A[it+1],b,2,Nc,S);row['tau']=float(times[it+1]);history.append(row)
    psi=b+u/r;v=multipole_potential(psi,q,r,dr,A[-1],2)
    probe_atlas=2*prep['reservoir_outer_radius'];probe_x=probe_atlas*1e6/A[-1]
    vector=np.array([np.interp(probe_x,r,v[:,j]) for j in range(4)])
    amplitude=1e-6*Nc*np.linalg.norm(vector)
    record=dict(config=asdict(c),source_final_field_sha256=hashlib.sha256(
        np.asarray(np.column_stack((r,u.real,u.imag)),dtype='<f8').tobytes()).hexdigest(),
        initial_basis_condition=condition,final=history[-1],
        maximum_field_norm_gain=max(x['field_norm_gain'] for x in history),
        maximum_core_density_gain=max(x['core_density_gain'] for x in history),
        maximum_symplectic_gap=max(x['symplectic_gap'] for x in history),
        probe_radius_atlas=probe_atlas,potential_RMS_coefficient_per_initial_fractional_field_norm=float(amplitude),
        potential_pointwise_coefficient_per_initial_fractional_field_norm=float(np.sqrt(5)*amplitude),
        boundary='Leading harmonic tidal phase: q proportional to r^-2, implemented as a real self-adjoint outer diagonal correction. No claim of an exact transparent boundary for arbitrary radiative perturbations.',
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    if retain:return record,dict(radius=r,final_u=u,final_q=q,history=history)
    return record
