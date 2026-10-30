"""A same-state massive scalar handoff through a smooth observer horizon.

Two intersecting null traces are restrictions of the retained global parent.
The conditional two-mode correlation fixes the second trace; it is not a
new causal signal from the horizon or a daughter spacetime.
"""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import json
import math
from pathlib import Path
import time

import numpy as np
from numba import njit
from scipy.integrate import quad,simpson

from .nsc_atlas_horizon_response import solve_source
from .nsc_atlas_parent_circulation import horizon_flux,static_patch_energy

ROOT=Path(__file__).resolve().parents[2]
C=np.sqrt([2.,12.])/math.pi


@dataclass(frozen=True)
class Diamond:
    u0: float
    v0: float
    width: float = .125

    @classmethod
    def inherited(cls,parent):
        record=json.loads((ROOT/'results/summaries/nsc-atlas-parent-circulation-v1.json').read_bytes())
        eta=math.atan(math.sinh(parent.H*record['background']['finite_proper_endpoint']))
        return cls(2*eta-math.pi/2,math.pi/2)

    @property
    def eta_corner(self):return (self.u0+self.v0)/2

    @property
    def chi_receiver(self):return (self.v0-self.u0)/2

    def axes(self,n):
        if type(n) is not int or n<8 or n>512:
            raise ValueError('finite 8 <= cells per null leg <= 512 required')
        if self.width<=0 or self.eta_corner+self.width>=math.pi/2:
            raise ValueError('diamond must remain before the conformal endpoint')
        if self.chi_receiver-self.width/2<=0 or self.chi_receiver+self.width/2>=math.pi:
            raise ValueError('diamond must avoid the spherical-coordinate poles')
        return np.linspace(self.u0,self.u0+self.width,n+1),np.linspace(self.v0,self.v0+self.width,n+1)


def geometry(source,u,v):
    u,v=np.broadcast_arrays(u,v);eta=(u+v)/2;chi=(v-u)/2
    a=1/(source.parent.H*np.cos(eta))
    potential=(source.parent.mass_squared/source.parent.H**2-2)/np.cos(eta)**2+2/np.sin(chi)**2
    return eta,chi,a,potential


def inherited_reference(source,u,v):
    """Signed radial fields and first derivatives, evaluated from the parent."""
    eta,chi,a,_=geometry(source,u,v);w=np.arcsinh(np.tan(eta))
    if np.min(w)<0 or np.max(w)>source.endpoint_u:
        raise ValueError('outside the declared parent field interval')
    y=source.modes.sol(w.ravel()).reshape((4,)+eta.shape)
    q=np.stack([y[0],y[2]],axis=-1)*np.cos(eta)[...,None]**1.5
    p=source.parent.H*np.cos(eta)[...,None]**1.5*(
        np.stack([y[1],y[3]],axis=-1)-1.5*np.sin(eta)[...,None]*np.stack([y[0],y[2]],axis=-1))
    sin,cos=np.sin(chi),np.cos(chi)
    radial=np.stack([sin*sin,sin*sin*cos],axis=-1)
    radial_derivative=np.stack([2*sin*cos,2*sin*cos*cos-sin**3],axis=-1)
    psi=a[...,None]*C*q*radial
    eta_derivative=np.tan(eta)[...,None]*psi+a[...,None]**2*C*p*radial
    chi_derivative=a[...,None]*C*q*radial_derivative
    return psi,(eta_derivative-chi_derivative)/2,(eta_derivative+chi_derivative)/2


def inherited_second_leg(source,domain,v):
    """Second trace from the signed horizon trace in the frozen two-mode sector.

    The paired points have equal eta and are spacelike separated. This is
    a global-state correlation, not a signal emitted along the horizon.
    """
    v=np.asarray(v);uh=domain.u0+v-domain.v0
    f,_,_=inherited_reference(source,uh,np.full_like(v,domain.v0))
    eta=(domain.u0+v)/2;chi=(v-domain.u0)/2;hchi=math.pi/2-eta
    ratio1=np.sin(chi)**2/np.sin(hchi)**2
    if np.min(abs(np.cos(hchi)))<1e-8:
        raise ValueError('use the derivative limit at a mode-2 horizon node')
    ratio=np.stack([ratio1,ratio1*np.cos(chi)/np.cos(hchi)],axis=-1)
    return f*ratio,ratio


@njit(cache=True,nogil=True)
def diamond_advance(f,g,potential,area_step):
    """Implicit four-corner midpoint rule; two radial components, no reset."""
    n,m=potential.shape
    scale=max(np.max(np.abs(f[0])),np.max(np.abs(g[0])),np.finfo(np.float64).tiny)
    if np.max(np.abs(f[0]-g[0]))>64*np.finfo(np.float64).eps*scale:
        raise ValueError('the two inherited traces must agree at the corner')
    psi=np.empty((n+1,m+1,f.shape[1]))
    psi[:,0,:]=f
    psi[0,:,:]=g
    for i in range(n):
        for j in range(m):
            c=area_step*potential[i,j]/16
            for k in range(f.shape[1]):
                psi[i+1,j+1,k]=(1-c)/(1+c)*(psi[i+1,j,k]+psi[i,j+1,k])-psi[i,j,k]
    return psi


def fields_and_derivatives(source,u,v,psi,du,dv):
    eta,chi,a,_=geometry(source,u,v)
    r=a*np.sin(chi);cot=1/np.tan(chi)
    F=psi/r[...,None]
    Fu=(du-.5*(np.tan(eta)-cot)[...,None]*psi)/r[...,None]
    Fv=(dv-.5*(np.tan(eta)+cot)[...,None]*psi)/r[...,None]
    return F,Fu,Fv


def currents(source,u,v,psi,du,dv):
    _,chi,a,_=geometry(source,u,v)
    F,Fu,Fv=fields_and_derivatives(source,u,v,psi,du,dv)
    N=4*math.pi/3
    Tuu=N*np.sum(Fu*Fu,axis=-1);Tvv=N*np.sum(Fv*Fv,axis=-1)
    Tuv=N/4*(2/np.sin(chi)**2+a*a*source.parent.mass_squared)*np.sum(F*F,axis=-1)
    Ku=source.parent.H*np.cos(u);Kv=source.parent.H*np.cos(v)
    Ju=a*a*np.sin(chi)**2*(Tuv*Ku+Tvv*Kv)
    Jv=a*a*np.sin(chi)**2*(Tuu*Ku+Tuv*Kv)
    Wu=psi[...,0]*du[...,1]-psi[...,1]*du[...,0]
    Wv=psi[...,0]*dv[...,1]-psi[...,1]*dv[...,0]
    return Ju,Jv,Wv,Wu


def flux_balance(u,v,Ju,Jv):
    terms=np.array([simpson(Ju[-1],x=v),-simpson(Ju[0],x=v),
                    simpson(Jv[:,-1],x=u),-simpson(Jv[:,0],x=u)])
    absolute=abs(float(sum(terms)));scale=float(np.sum(abs(terms)))
    return dict(signed_edges=terms.tolist(),absolute_residual=absolute,
                relative_residual=absolute/scale if scale else 0.)


def discrete_symplectic_balance(psi):
    """Exact diamond Green identity, avoiding differentiated-field cancellation."""
    A,B=psi[...,0],psi[...,1]
    uv=A[:,:-1]*B[:,1:]-B[:,:-1]*A[:,1:]
    vu=A[:-1,:]*B[1:,:]-B[:-1,:]*A[1:,:]
    terms=np.array([np.sum(uv[-1]),-np.sum(uv[0]),np.sum(vu[:,-1]),-np.sum(vu[:,0])])
    absolute=abs(float(sum(terms)));scale=float(np.sum(abs(terms)))
    return dict(signed_edges=terms.tolist(),absolute_residual=absolute,
                relative_residual=absolute/scale if scale else 0.)


def receiver_record(source,domain,psi,du,dv,axis):
    u,v=domain.u0+domain.width,domain.v0+domain.width
    eta,chi,a,_=geometry(source,u,v)
    F,Fu,Fv=fields_and_derivatives(source,u,v,psi,du,dv)
    n=np.asarray(axis[:2]);field=float(np.dot(n,F));p=float(np.dot(n,Fu+Fv))/float(a)
    radial=float(np.dot(n,Fv-Fu))
    angular=float(np.dot(F,F)-field*field)
    rho=(p*p+(radial*radial+angular/math.sin(float(chi))**2)/float(a*a)
         +source.parent.mass_squared*field*field)/2
    q1=F[0]/(C[0]*math.sin(float(chi)))
    p1=(Fu[0]+Fv[0])/(float(a)*C[0]*math.sin(float(chi)))
    phase=math.atan2(-p1,source.parent.reference_frequency*q1)
    return dict(eta=float(eta),chi=float(chi),proper_readout=math.asinh(math.tan(float(eta)))/source.parent.H,
                field_over_epsilon=field,proper_derivative_over_epsilon=p,
                density_over_epsilon_squared_vacuum=rho/source.parent.vacuum,
                inherited_first_mode_phase_mod_2pi=phase)


def evaluate_case(source,domain,n,axis,*,alter_second_leg=False,return_arrays=False):
    start=time.perf_counter();u,v=domain.axes(n);h=domain.width/n
    f=inherited_reference(source,u,np.full_like(u,domain.v0))[0]
    g,ratios=inherited_second_leg(source,domain,v)
    direct_g=inherited_reference(source,np.full_like(v,domain.u0),v)[0]
    relation_error=float(np.max(abs(g-direct_g)))
    if alter_second_leg:
        z=(v-domain.v0)/domain.width;bump=np.zeros_like(z);inside=(z>0)&(z<1)
        bump[inside]=np.exp(4-1/(z[inside]*(1-z[inside])))
        g=g.copy();g[:,0]+=.05*float(np.max(abs(g)))*bump
    U,V=np.meshgrid(u,v,indexing='ij')
    potential=geometry(source,(U[:-1,:-1]+U[1:,1:])/2,(V[:-1,:-1]+V[1:,1:])/2)[3]
    psi=diamond_advance(f,g,potential,h*h)
    exact,exact_du,exact_dv=inherited_reference(source,U,V)
    du,dv=np.gradient(psi,h,h,axis=(0,1),edge_order=2)
    rms=float(np.sqrt(np.mean(exact*exact)))
    error=psi-exact
    Ju,Jv,Su,Sv=currents(source,U,V,psi,du,dv)
    exactJ=currents(source,U,V,exact,exact_du,exact_dv)
    receiver=receiver_record(source,domain,psi[-1,-1],du[-1,-1],dv[-1,-1],axis)
    direct=receiver_record(source,domain,exact[-1,-1],exact_du[-1,-1],exact_dv[-1,-1],axis)
    phase_error=receiver['inherited_first_mode_phase_mod_2pi']-direct['inherited_first_mode_phase_mod_2pi']
    row=dict(cells_per_leg=n,relative_RMS_field_error=float(np.sqrt(np.mean(error*error)))/rms,
             relative_max_field_error=float(np.max(abs(error)))/float(np.max(abs(exact))),
             horizon_normal_derivative_relative_RMS=float(np.sqrt(np.mean((dv[:,0]-exact_dv[:,0])**2)))/float(np.sqrt(np.mean(exact_dv[:,0]**2))),
             inherited_second_trace_max_error=relation_error,
             receiver=receiver,reference_receiver=direct,
             receiver_wrapped_phase_error=math.atan2(math.sin(phase_error),math.cos(phase_error)),
             noether_flux=flux_balance(u,v,Ju,Jv),symplectic_flux=flux_balance(u,v,Su,Sv),
             discrete_symplectic_flux=discrete_symplectic_balance(psi),
             exact_reference_noether_flux=flux_balance(u,v,exactJ[0],exactJ[1]),
             transfer_factors_at_second_leg_end=ratios[-1].tolist(),elapsed_seconds=time.perf_counter()-start)
    if return_arrays:return row,(u,v,psi,exact)
    return row


def calculation():
    source=solve_source(endpoint_u=2.)
    domain=Diamond.inherited(source.parent)
    earlier=json.loads((ROOT/'results/summaries/nsc-atlas-horizon-response-v1.json').read_bytes())
    axis=earlier['angular_transfer']['axis_in_n_coordinates']
    # Compile once; independent grids then run without the Python GIL.
    diamond_advance(np.zeros((2,2)),np.zeros((2,2)),np.ones((1,1)),.01)
    with ThreadPoolExecutor(max_workers=3) as pool:
        rows=list(pool.map(lambda n:evaluate_case(source,domain,n,axis),(64,128,256)))
    fine,arrays=evaluate_case(source,domain,256,axis,return_arrays=True)
    altered,altarrays=evaluate_case(source,domain,256,axis,alter_second_leg=True,return_arrays=True)
    reference=arrays[2];changed=altarrays[2]
    rms=float(np.sqrt(np.mean(reference*reference)))
    one_trace_control=dict(second_leg_change='C-infinity bump, 5% reference amplitude, only on the transverse datum',
        horizon_trace_max_difference=float(np.max(abs(changed[:,0]-reference[:,0]))),
        interior_relative_RMS_change=float(np.sqrt(np.mean((changed[1:,1:]-reference[1:,1:])**2)))/rms,
        receiver_field_change_over_epsilon=altered['receiver']['field_over_epsilon']-fine['receiver']['field_over_epsilon'],
        interpretation='One null trace alone is not unique for the full KG solution space. The frozen two-mode global-state correlation supplies the second trace in the main calculation.')
    u,v,psi,_=arrays;U,V=np.meshgrid(u,v,indexing='ij');h=domain.width/256
    potential=geometry(source,(U[:-1,:-1]+U[1:,1:])/2,(V[:-1,:-1]+V[1:,1:])/2)[3]
    middle=128
    left=diamond_advance(psi[:,0],psi[0,:middle+1],potential[:,:middle],h*h)
    right=diamond_advance(left[:,-1],psi[0,middle:],potential[:,middle:],h*h)
    composed=np.concatenate([left,right[:,1:]],axis=1)
    composition=dict(max_field_difference=float(np.max(abs(composed-psi))),
                     shared_trace_is_first_region_output=True,
                     other_incoming_leg='Unchanged remainder of the same inherited transverse trace',
                     no_new_source_reset=True)
    eta0=domain.eta_corner;eta1=eta0+domain.width/2
    t0=math.asinh(math.tan(eta0))/source.parent.H;t1=math.asinh(math.tan(eta1))/source.parent.H
    outgoing,error=quad(lambda t:horizon_flux(t,source.state(source.parent.H*t),source.parent)[1],t0,t1,
                         epsabs=2e-12,epsrel=2e-12)
    volume_loss=static_patch_energy(t0,source.state(source.parent.H*t0),source.parent)-static_patch_energy(t1,source.state(source.parent.H*t1),source.parent)
    horizon_connection=dict(proper_start=t0,proper_end=t1,
        parent_patch_energy_loss_over_epsilon_squared=volume_loss,
        same_horizon_flux_over_epsilon_squared=outgoing,
        relative_volume_flux_residual=abs(volume_loss-outgoing)/abs(volume_loss),
        quadrature_error=error,
        finite_grid_horizon_flux_relative_difference=abs(-fine['noether_flux']['signed_edges'][3]/outgoing-1))
    orders={key:[math.log(rows[i][key]/rows[i+1][key],2) for i in (0,1)]
            for key in ('relative_RMS_field_error','horizon_normal_derivative_relative_RMS')}
    return dict(schema='nsc-atlas-characteristic-handoff-v1',
        model=dict(vacuum=source.parent.vacuum,H=source.parent.H,mass_squared=source.parent.mass_squared,
                   angular_momentum_l=1,global_S3_eigenvalues=[3,8],no_new_source_normalization=True),
        domain=dict(u0=domain.u0,v0=domain.v0,null_extent=domain.width,
                    eta_corner=domain.eta_corner,chi_receiver=domain.chi_receiver,
                    conformal_future_margin=math.pi/2-domain.eta_corner-domain.width,
                    receiver_angular_axis=axis,window_is_a_declared_finite_domain=True),
        exact_relations=dict(equation='4 Psi_uv+[(m^2-2)sec^2 eta+2csc^2 chi]Psi=0',
            horizon_potential='V|v=pi/2=m^2 sec^2 eta',
            second_trace='g1=f1(u0+v-v0)*sin^2chi/sin^2chiH; g2 adds coschi/coschiH',
            second_trace_is_global_state_correlation=True,spacelike_signal_created=False,
            two_trace_uniqueness='Volterra iterates bounded by (M du dv/4)^n/(n!)^2',
            constant_potential_dispersion='cos(omega*h/2)=((1-c)/(1+c))*cos(k*h/2), c=V*h^2/16'),
        resolutions=rows,observed_orders=orders,one_trace_control=one_trace_control,
        composition=composition,horizon_energy_connection=horizon_connection,
        finite_resolution_indicator=dict(
            receiver_field_fine_middle_difference=abs(rows[2]['receiver']['field_over_epsilon']-rows[1]['receiver']['field_over_epsilon']),
            receiver_density_fine_middle_difference=abs(rows[2]['receiver']['density_over_epsilon_squared_vacuum']-rows[1]['receiver']['density_over_epsilon_squared_vacuum']),
            meaning='Empirical discretization differences, not higher-order gravitational or observational uncertainties'),
        execution=dict(independent_grid_workers=3,numba_nogil_kernel=True),
        scope=['Leading inherited real scalar on the smooth frozen parent; no material wall',
               'A finite exterior receiving region, not a daughter universe',
               'Both traces derived conditionally from the same two-mode parent state; the primordial global preparation is not selected',
               'North-Killing Noether charge can have either sign outside its patch; receiver density uses a timelike comoving observer',
               'Scalar metric feedback, finite-amplitude collapse, regular child continuation and an SI ruler remain open'],
        no_observational_fit=True)
