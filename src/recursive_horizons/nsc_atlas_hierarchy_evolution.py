"""Finite classical envelope evolution of the supplied two-scale parent data.

This is a leading weak-field/multiple-time-scale reduction, not the full
Einstein--Klein--Gordon stepper. Parent expansion and the coherent ambient
field remain dynamical inputs. No absorber removes the reservoir.
"""
from dataclasses import dataclass,asdict
import hashlib
import math
from pathlib import Path
import time

import numpy as np
from scipy.fft import dst,idst
from scipy.integrate import solve_ivp

from .nsc_atlas_scale_hierarchy import ScaleHierarchy
from .nsc_atlas_core_reservoir import compact_bump


@dataclass(frozen=True)
class HierarchyRun:
    frequency: float=1e9
    weak_amplitude: float=.001
    depletion: float=.5
    outer: float=1024.
    points: int=8191
    step: float=.02
    end: float=600.
    gravity: bool=True
    expansion: bool=True


def background(model,end,*,expansion=True):
    s=model.reference;scale=model.frequency*model.weak_amplitude**2
    vacuum=model.vacuum;dust=s.scalar.kappa*model.excitation;C=model.curvature
    def rhs(t,y):
        A=y[0]
        rate=math.sqrt(vacuum*A*A/3+dust/(3*A)-C)/scale if expansion else 0.
        return rate,1/(A*A)
    sol=solve_ivp(rhs,(0,end),(1.,0.),method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def initial_envelope(model,r):
    record=model.solve();w=model.frequency;z=model.weak_amplitude
    normal=math.sqrt(2/model.reference.scalar.kappa)*z*z
    ambient=model.ambient_coefficient/(w*normal)
    R=r/(w*z);core=model.core.field(R)[0]/normal
    b,_=compact_bump(R,model.radius,record['reservoir_outer_radius'])
    return np.asarray(core-model.depletion*ambient*b,dtype=complex),ambient,record


def potential(u,r,dr,A,ambient,gravity=True):
    # u=r*(psi-ambient); compute the contrast without subtracting
    # nearly equal background densities. The discrete shell kernel is
    # symmetric in its radial quadrature weights.
    contrast=abs(u/r)**2+2*ambient*np.real(u/r)
    radial=r*r*contrast
    if not gravity:return np.zeros_like(r),radial
    mass=dr*np.cumsum(radial);tail=dr*np.cumsum(radial/r)
    return -(mass/r+tail[-1]-tail)/(2*A),radial


def kinetic_step(u,k,integrated_inverse_scale_squared):
    return idst(dst(u,type=1,norm='ortho')*np.exp(-.5j*k*k*integrated_inverse_scale_squared),
        type=1,norm='ortho')


def enclosed(r,radial,dr,cut):
    if not 0<cut<r[-1]:raise ValueError('mass cut outside the interior')
    integral=dr*(np.cumsum(radial)-radial/2)
    return float(np.interp(cut,np.r_[0.,r],np.r_[0.,integral]))


def cut_flux(u,r,dr,A,h,ambient,cut,physical=True):
    """Outward current, including the moving comoving cut of fixed R."""
    x=cut/A if physical else cut
    j=max(1,min(len(r)-3,int(x/dr)-1));nodes=r[j-1:j+3];values=u[j-1:j+3]
    # A cubic through four neighbouring radial values, differentiated
    # independently of the enclosed-mass quadrature.
    q=(x-nodes[1])/dr
    weights=(-q*(q-1)*(q-2)/6,(q+1)*(q-1)*(q-2)/2,-(q+1)*q*(q-2)/2,(q+1)*q*(q-1)/6)
    slopes=(-(3*q*q-6*q+2)/6,(3*q*q-4*q-1)/2,-(3*q*q-2*q-2)/2,(3*q*q-1)/6)
    value=sum(weight*v for weight,v in zip(weights,values))
    derivative=sum(weight*v for weight,v in zip(slopes,values))/dr
    psi=ambient+value/x;psir=derivative/x-value/x**2
    current=x*x*np.imag(np.conj(psi)*psir)/(A*A)
    contrast=abs(value/x)**2+2*ambient*np.real(value/x)
    return float(current+(h*x**3*contrast if physical else 0.))


def run_hierarchy(config=HierarchyRun(),*,retain=False):
    c=config;cpu=time.process_time();wall=time.perf_counter()
    if (c.points<1023 or c.points>32767 or not 128<=c.outer<=2048
            or not 0<c.step<=.04 or not 0<c.end<=600):
        raise ValueError('outside declared finite envelope domain')
    steps=round(c.end/c.step)
    if not math.isclose(steps*c.step,c.end,abs_tol=1e-10):raise ValueError('endpoint is not a whole step')
    model=ScaleHierarchy(c.frequency,c.weak_amplitude,c.depletion)
    dr=c.outer/(c.points+1);r=dr*np.arange(1,c.points+1);k=np.pi*np.arange(1,c.points+1)/c.outer
    delta,ambient,initial=initial_envelope(model,r);u=r*delta
    if model.frequency*model.weak_amplitude*initial['reservoir_outer_radius']>.7*c.outer:
        raise ValueError('outer domain does not retain a sufficiently large ambient collar')
    parent=background(model,c.end,expansion=c.expansion)
    times=np.arange(steps+1)*c.step;A,Q=parent.sol(times)
    Am=parent.sol(times[:-1]+c.step/2)[0]
    slow=c.frequency*c.weak_amplitude**2
    def hubble(a):
        return math.sqrt(model.vacuum/3+model.reference.scalar.kappa*model.excitation/(3*a**3)
            -model.curvature/a**2)/slow if c.expansion else 0.
    h=np.array([hubble(a) for a in A])
    def diagnostics(u,a):
        V,radial=potential(u,r,dr,a,ambient,c.gravity)
        T=float(dr*np.sum(k*k*abs(dst(u,type=1,norm='ortho'))**2)/(2*a*a))
        W=float(dr*np.sum(V*radial)/2)
        return T,W,float(dr*np.sum(radial)),radial
    T,W,N,radial=diagnostics(u,1.);E0=T+W;N0=N
    initial_core=enclosed(r,radial,dr,20.)
    cuts=(20.,40.);inside0=np.array([enclosed(r,radial,dr,x) for x in cuts])
    integrated=np.zeros(2);oldflux=np.array([cut_flux(u,r,dr,1.,h[0],ambient,x) for x in cuts])
    def boundary_flux(u,a):
        derivative=(-48*u[-1]+36*u[-2]-16*u[-3]+3*u[-4])/(12*dr)
        return float(ambient*c.outer*np.imag(derivative)/(a*a))
    exterior_flux=0.;old_boundary=boundary_flux(u,1.)
    work=0.;lastwork=h[0]*(2*T+W);lastworktime=0.
    stride=max(1,round(.1/c.step));rows=[];maximum_balance=np.zeros(2);max_norm_gap=0.;max_energy_gap=0.
    initial_u=u.copy()
    for iteration in range(steps):
        # Midpoint scale factor plus the integrated kinetic coefficient
        # gives a symmetric second-order nonautonomous split step.
        V,_=potential(u,r,dr,Am[iteration],ambient,c.gravity)
        u=(u+ambient*r)*np.exp(-.5j*c.step*V)-ambient*r
        u=kinetic_step(u,k,Q[iteration+1]-Q[iteration])
        V,_=potential(u,r,dr,Am[iteration],ambient,c.gravity)
        u=(u+ambient*r)*np.exp(-.5j*c.step*V)-ambient*r
        a=A[iteration+1];hh=h[iteration+1]
        flux=np.array([cut_flux(u,r,dr,a,hh,ambient,x) for x in cuts])
        integrated+=c.step*(oldflux+flux)/2;oldflux=flux
        outgoing=boundary_flux(u,a);exterior_flux+=c.step*(old_boundary+outgoing)/2;old_boundary=outgoing
        if (iteration+1)%stride==0 or iteration==steps-1:
            tau=times[iteration+1];T,W,N,radial=diagnostics(u,a)
            newwork=hh*(2*T+W);work+=(tau-lastworktime)*(lastwork+newwork)/2
            lastwork=newwork;lastworktime=tau
            inside=np.array([enclosed(r,radial,dr,x/a) for x in cuts])
            balance=inside-inside0+integrated;maximum_balance=np.maximum(maximum_balance,abs(balance))
            normgap=N-N0+exterior_flux;energygap=T+W-E0+work
            max_norm_gap=max(max_norm_gap,abs(normgap));max_energy_gap=max(max_energy_gap,abs(energygap))
            cumulative=dr*(np.cumsum(radial)-radial/2);eligible=r<=20/a
            half=float(np.interp(inside[0]/2,np.r_[0.,cumulative[eligible]],np.r_[0.,r[eligible]]))*a if inside[0]>0 else None
            central=(4*u[0]/r[0]-u[1]/r[1])/3+ambient
            rows.append(dict(tau=float(tau),scale_factor=float(a),core_mass_ratio=float(inside[0]/initial_core),
                physical_half_mass_radius=half,central_squared_amplitude=float(abs(central)**2/a**3),
                mass_flux_gaps=balance.tolist(),norm_balance=normgap,energy_work_balance=energygap))
    tail=[row for row in rows if row['tau']>=.8*c.end]
    digest=hashlib.sha256(np.asarray(np.column_stack((r,u.real,u.imag)),dtype='<f8').tobytes()).hexdigest()
    result=dict(config=asdict(c),preparation=initial,background_envelope=ambient,
        initial_envelope_core_geometric_mass=c.weak_amplitude*initial_core/(2*c.frequency),
        initial_GR_core_geometric_mass=initial['core']['mass'],
        initial_envelope_total_excess_over_core=N0/initial_core,
        final_scale_factor=float(A[-1]),final_core_mass_ratio=rows[-1]['core_mass_ratio'],
        late_interval=[.8*c.end,c.end],late_mean_core_mass_ratio=float(np.mean([x['core_mass_ratio'] for x in tail])),
        late_minimum_core_mass_ratio=min(x['core_mass_ratio'] for x in tail),
        late_mean_physical_half_mass_radius=float(np.mean([x['physical_half_mass_radius'] for x in tail])),
        maximum_core_mass_flux_gaps_over_initial_core=(maximum_balance/initial_core).tolist(),
        maximum_total_mass_flux_gap_over_initial_core=max_norm_gap/initial_core,
        maximum_energy_work_gap_over_initial_energy=max_energy_gap/abs(E0),
        outer_net_flux_over_initial_core=exterior_flux/initial_core,
        end_atlas_proper_time=c.end/slow,initial_weak_field_parameter=c.weak_amplitude**2,
        initial_H_over_frequency=model.reference.expansions[1]/c.frequency,
        parent_spatial_curvature_at_initial_reservoir=model.curvature*initial['reservoir_outer_radius']**2,
        final_field_sha256=digest,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    if retain:return result,dict(radius=r,initial_u=initial_u,final_u=u,history=rows)
    return result


def full_equation_comparison(zeta=.08,*,refinement=1):
    """Short carrier-resolved comparison along a weak-field family.

    Keep omega*zeta^2=12.8 and the atlas duration .4 fixed. Halving
    zeta weakens both relativistic and fast-expansion corrections.
    The original full EKG source functions remain unchanged.
    """
    from .nsc_atlas_core_evolution import moving_mesh,evolve_arrays,geometry,interpolate
    cpu=time.process_time();wall=time.perf_counter();omega=12.8/zeta**2
    model=ScaleHierarchy(omega,zeta,.5);prep=model.solve();end=.4
    # The fine comparison has a constant scaled spatial spacing .08
    # and carrier step omega*dt=.1 across the two physical members.
    dr=.08/(omega*zeta)/refinement;dt=.1/omega/refinement
    n=round(1/dr);dr=1/n;r=np.arange(n+1)*dr
    delta=model.core.field(r)[0]
    b,_=compact_bump(r,model.radius,prep['reservoir_outer_radius'])
    delta-=model.depletion*model.ambient_coefficient/omega*b
    Pi=np.zeros_like(r);s=model.reference
    bg=np.array([model.ambient_coefficient/omega,0.,s.cap_geometry['parent_maximum_radius'],s.expansions[1]])
    mesh=moving_mesh(r);a0=bg[2];cuts=np.array([model.radius,2*model.radius])
    steps=round(end/dt)
    out=evolve_arrays(delta,Pi,bg,dr,dt,steps,omega,1.,9/32,.7,cuts,max(1,steps//20),mesh,a0,True,0.)
    g=geometry(out[0],out[1],out[2],dr,omega,1.,mesh,a0)
    mass=interpolate(g[3],model.radius,g[8]);initial=out[-1][0]
    return dict(zeta=zeta,omega=omega,refinement=refinement,end_atlas_time=end,end_slow_time=12.8*end,
        core_mass_ratio=float(mass/initial),initial_core_mass=float(initial),final_core_mass=float(mass),
        maximum_mass_flux_gap_over_initial_core=float(np.max(out[6])/initial),
        minimum_lapse=out[8],minimum_radial_chart=out[7],final_scale_factor=float(out[2][2]/a0),
        initial_support=prep['reservoir_outer_radius'],outer_radius=1.,
        maximum_grid_characteristic_speed=out[9],
        finite_boundary_margin=float(1-prep['reservoir_outer_radius']-out[9]*end),
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
