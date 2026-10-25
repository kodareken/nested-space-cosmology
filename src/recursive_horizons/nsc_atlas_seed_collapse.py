"""Nonlinear compact-seed evolution in the declared earlier parent.

Initial source constraints and reduced evolution are distinct calculations.
No stationary-core profile is used to prepare the field.
"""
from dataclasses import dataclass,asdict
import hashlib
import math
import time

import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
from scipy.fft import dst

from .nsc_atlas_formation_window import FormationBackground
from .nsc_atlas_core_reservoir import compact_bump
from .nsc_atlas_hierarchy_evolution import potential,kinetic_step,enclosed,cut_flux


LENGTH=.0004815623798469198


def boundary_potential(u,r,dr,outer,gravity=True):
    """Symmetric Poisson Green function with V(outer)=0."""
    V,radial=potential(u,r,dr,1.,1.,gravity)
    if gravity:V=V+dr*np.sum(radial)/(2*outer)
    return V,radial


def seed_shape(r):
    """Exact -Laplacian(g*g)/int|grad g|^2, g=(1-r^2)^4_+."""
    r=np.asarray(r,dtype=float);x=np.clip(r,0,2)
    poly=np.polyval([21,336,2310,8736,18944,20864,4096,-9728,-8192,-2048],x)
    deriv=np.polyval([21,294,1722,5292,8320,3584,-7168,-7168,-2048],x)
    value=-(x-2)**8*poly/524288
    slope=-17*x*(x-2)**7*deriv/524288
    return np.where((r>=0)&(r<2),value,0.),np.where((r>=0)&(r<2),slope,0.)


def seed_mass_primitive(r):
    x=np.clip(np.asarray(r,dtype=float),0,2)
    poly=np.polyval([63,1134,8960,40320,111360,187392,178176,92160,20480],x)
    return -x**3*(x-2)**9*poly/31457280


def initial_constraint(amplitude=.05,Ai=.01,L=LENGTH,*,rtol=2e-11):
    b=FormationBackground.inherited();w=b.frequency;D=b.excitation
    q=1/(w*w*Ai*Ai*L*L);gamma=D*L*L/Ai;c=b.curvature*L*L
    def field(r,correction):
        S,Sp=seed_shape(r);W,Wp=compact_bump(r,1.,2.)
        delta=amplitude*S-correction*W
        if np.any(1+delta<=0):raise ValueError('negative initial squared envelope')
        f=np.sqrt(1+delta);fp=(amplitude*Sp-correction*Wp)/(2*f)
        return delta,f,fp
    def integrate(correction):
        def rhs(r,n):
            delta,f,fp=field(r,correction)
            B=1-c*r*r-(gamma*n[0]/r if r else 0.)
            if B<=0:raise RuntimeError('initial radial chart is not positive')
            return (r*r*(delta+q*B*fp*fp),)
        sol=solve_ivp(rhs,(0,2),(0.,),method='DOP853',rtol=rtol,atol=rtol*1e-6,max_step=.025,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        return sol
    first=integrate(0.);J=quad(lambda r:r*r*compact_bump(r,1.,2.)[0],1,2,epsabs=1e-13)[0]
    # Solve in a scaled correction to avoid an absolute-tolerance floor.
    reference=max(abs(first.y[0,-1])/J,1e-30)
    scale=brentq(lambda z:integrate(z*reference).y[0,-1]/reference,0,3,xtol=2e-9)
    correction=scale*reference;sol=integrate(correction)
    r=np.linspace(1e-5,2,1001);n=sol.sol(r)[0];delta,f,fp=field(r,correction)
    B=1-c*r*r-gamma*n/r
    # Independent proper-distance Hamilton system: h=n/r, dr/dxi=sqrt(B).
    def proper_rhs(xi,y):
        r,h=y;delta,f,fp=field(r,correction);B=1-c*r*r-gamma*h
        beta=math.sqrt(B)
        hp=r*(delta+q*B*fp*fp)-(h/r if r else 0.)
        return beta,beta*hp
    def stop(xi,y):return y[0]-2
    stop.terminal=True;stop.direction=1
    proper=solve_ivp(proper_rhs,(0,3),(0.,0.),method='DOP853',rtol=rtol,atol=rtol*1e-6,
        max_step=.025,events=stop)
    if not proper.success or len(proper.t_events[0])!=1:raise RuntimeError('proper-distance check did not reach support edge')
    peak_radius=brentq(lambda r:float(seed_shape(r)[0]),.1,1.5)
    positive_norm=float(sol.sol(peak_radius)[0]);mass_unit=D*L**3/2
    record=dict(amplitude=amplitude,initial_scale=Ai,length=L,frequency=w,
        correction=correction,correction_reference=reference,
        uncorrected_GR_endpoint_norm=float(first.y[0,-1]),
        corrected_GR_endpoint_norm=float(sol.y[0,-1]),
        proper_distance_endpoint_norm=float(2*proper.y_events[0][0,1]),
        initial_envelope_total_norm=-correction*J,
        initial_positive_excess_norm=positive_norm,initial_positive_radius=peak_radius,
        initial_positive_geometric_mass=mass_unit*positive_norm,mass_unit=mass_unit,
        minimum_B=float(np.min(B)),minimum_squared_envelope=float(np.min(f*f)),
        gradient_mass_coefficient=q,curvature_support_parameter=4*c,
        initial_H_over_frequency=float(b.H(1/Ai)/w),
        envelope_source='sqrt(1+amplitude*S-correction*W), real; W is the fixed C4 annular bump on1<r<2. Correction solves the initial GR mass constraint, not observational fitting.',
        limit='Exact spherical turning-point initial constraint on the declared homogeneous parent patch. Envelope momentum differs at O(H/omega); no earlier material-wall solution is asserted.')
    return record,lambda r:field(r,correction)


@dataclass(frozen=True)
class SeedRun:
    amplitude: float=.1
    initial_scale: float=.01
    final_scale: float=1.2
    length: float=LENGTH
    outer: float=8.
    points: int=4095
    log_step: float=.002
    phase_cap: float=.08
    gravity: bool=True


def shell_positions(r,density,targets,central):
    """Positive-density, exact-volume integration with a linear cell profile."""
    x=np.r_[0.,r];n=np.r_[central,density]
    a=x[:-1];b=x[1:];slope=np.diff(n)/np.diff(x)
    cell=n[:-1]*(b**3-a**3)/3+slope*((b**4-a**4)/4-a*(b**3-a**3)/3)
    cum=np.r_[0.,np.cumsum(cell)];out=[]
    for target in targets:
        j=min(len(a)-1,max(0,int(np.searchsorted(cum,target))-1))
        def value(z):return cum[j]+n[j]*(z**3-a[j]**3)/3+slope[j]*((z**4-a[j]**4)/4-a[j]*(z**3-a[j]**3)/3)-target
        # A target at a cell edge may differ by an ulp when a cumulative
        # sum and the primitive are rounded separately. Accept only that
        # explicit endpoint tolerance; a genuine missing bracket still fails.
        tol=32*np.finfo(float).eps*max(abs(target),abs(cum[j+1]),1e-300)
        if abs(value(a[j]))<=tol:out.append(float(a[j]))
        elif abs(value(b[j]))<=tol:out.append(float(b[j]))
        else:out.append(brentq(value,a[j],b[j],xtol=1e-14))
    return out,cum


def evolve_seed(config=SeedRun(),*,retain=False):
    c=config;cpu=time.process_time();wall=time.perf_counter();bg=FormationBackground.inherited()
    if not (0<c.amplitude<=.1 and c.initial_scale==.01 and 1<=c.final_scale<=16
            and 4<=c.outer<=16 and 1023<=c.points<=32767 and 0<c.log_step<=.004 and 0<c.phase_cap<=.1):
        raise ValueError('outside admitted compact-seed evolution')
    prep,initial=initial_constraint(c.amplitude,c.initial_scale,c.length)
    nodes=np.linspace(math.log(c.initial_scale),math.log(c.final_scale),601)
    A=np.exp(nodes)
    dr=c.outer/(c.points+1);r=dr*np.arange(1,c.points+1);k=np.pi*np.arange(1,c.points+1)/c.outer
    delta,f,fp=initial(r);u=np.asarray(r*(f-1),dtype=complex);initial_u=u.copy()
    cuts=(.25,.5,1.,2.)
    def field_potential(u):
        # The field has fixed boundary phase. Its potential must use that
        # same boundary clock: Vg=V_infinity-V_infinity(L). The added rank-one
        # constant leaves spatial forces unchanged and keeps the Green kernel
        # symmetric. It prevents an incompatible driven phase at u(L)=0.
        return boundary_potential(u,r,dr,c.outer,c.gravity)
    def alpha(a):return 1/(bg.frequency*c.length**2*a*a*float(bg.H(1/a)))
    def state(u,a):
        V,radial=field_potential(u)
        T=float(dr*np.sum(k*k*abs(dst(u,type=1,norm='ortho'))**2)/(2*bg.frequency*c.length*c.length*a*a))
        W=float(bg.frequency*bg.excitation*c.length*c.length/a*dr*np.sum(V*radial)/2)
        return V,radial,T,W
    def currents(u,a):return alpha(a)*np.array([cut_flux(u,r,dr,1.,0.,1.,x,False) for x in cuts])
    def boundary(u,a):
        derivative=(-48*u[-1]+36*u[-2]-16*u[-3]+3*u[-4])/(12*dr)
        return alpha(a)*c.outer*float(np.imag(derivative))
    V,radial,T,W=state(u,A[0]);E0=T+W;N0=dr*np.sum(radial)
    positive_norm=prep['initial_positive_excess_norm']
    inside0=np.array([enclosed(r,radial,dr,x) for x in cuts])
    shell_radii=(.05,.1,.25,.5,1.)
    total_initial=abs(1+u/r)**2
    _,cum=shell_positions(r,total_initial,[],1+c.amplitude)
    shell_masses=[float(np.interp(x,np.r_[0.,r],cum)) for x in shell_radii]
    flux=np.zeros(len(cuts));oldcurr=currents(u,A[0]);outflow=0.;oldbound=boundary(u,A[0])
    work=0.;oldwork=2*T+W
    max_flux=0.;max_total=0.;max_work=0.;history=[]
    def diagnostic(u,a,N):
        nonlocal max_flux,max_total,max_work
        V,radial,T,W=state(u,a);newwork=2*T+W
        inside=np.array([enclosed(r,radial,dr,x) for x in cuts]);gaps=inside-inside0+flux
        normgap=dr*np.sum(radial)-N0+outflow;energygap=T+W-E0+work
        max_flux=max(max_flux,float(np.max(abs(gaps)))/positive_norm)
        max_total=max(max_total,abs(float(normgap))/positive_norm)
        max_work=max(max_work,abs(float(energygap))/max(abs(E0),1e-30))
        psi=1+u/r;central=(4*psi[0]-psi[1])/3
        contrast=abs(psi)**2-1;central_contrast=float(abs(central)**2-1)
        crossings=np.flatnonzero(contrast<central_contrast/2)
        j=int(crossings[0]) if len(crossings) else len(r)-1
        half=float(np.interp(central_contrast/2,[contrast[j],contrast[max(0,j-1)]],[r[j],r[max(0,j-1)]])) if j else float(r[0])
        shells,_=shell_positions(r,abs(psi)**2,shell_masses,abs(central)**2)
        return dict(scale_factor=float(a),central_density_contrast=central_contrast,
            physical_half_height_radius=a*c.length*half,
            physical_mass_shell_radii=[a*c.length*x for x in shells],
            comoving_mass_shell_radii=shells,cut_excess_norms=inside.tolist(),
            total_excess_norm=float(dr*np.sum(radial)),potential_max=float(np.max(abs(V))*bg.excitation*c.length*c.length/a),
            mass_flux_gap_over_initial_positive=float(np.max(abs(gaps))/positive_norm),
            energy_work_gap_over_initial_energy=float(energygap/max(abs(E0),1e-30)))
    history.append(diagnostic(u,A[0],nodes[0]))
    N=nodes[0];steps=0;maxphase=0.;minimum_step=c.log_step
    for target in nodes[1:]:
        while target-N>1e-13:
            a=math.exp(N);V,_=field_potential(u)
            rate=bg.frequency*bg.excitation*c.length*c.length/(a*float(bg.H(1/a)))
            phase_rate=rate*float(np.max(abs(V)))
            dn=min(c.log_step,target-N,c.phase_cap/max(phase_rate,1e-30))
            minimum_step=min(minimum_step,dn)
            Ag=np.exp(N+dn/2+np.array([-1.,1.])*dn/(2*math.sqrt(3)));Hg=bg.H(1/Ag)
            Q=dn*np.mean(1/(bg.frequency*c.length**2*Ag**2*Hg))
            P=dn*np.mean(bg.frequency*bg.excitation*c.length**2/(Ag*Hg))
            before=float(np.max(abs(V)))
            u=(u+r)*np.exp(-.5j*P*V)-r
            u=kinetic_step(u,k,Q)
            V,_=field_potential(u)
            maxphase=max(maxphase,P*max(before,float(np.max(abs(V)))))
            u=(u+r)*np.exp(-.5j*P*V)-r
            N+=dn;a=math.exp(N);steps+=1
            now=currents(u,a);flux+=dn*(oldcurr+now)/2;oldcurr=now
            nb=boundary(u,a);outflow+=dn*(oldbound+nb)/2;oldbound=nb
            _,_,T,W=state(u,a);newwork=2*T+W
            work+=dn*(oldwork+newwork)/2;oldwork=newwork
        history.append(diagnostic(u,math.exp(target),target))
    final=history[-1];positions=np.array([h['physical_mass_shell_radii'] for h in history])
    maxima=np.argmax(positions,axis=0)
    digest=hashlib.sha256(np.asarray(np.column_stack((r,u.real,u.imag)),dtype='<f8').tobytes()).hexdigest()
    result=dict(config=asdict(c),preparation=prep,steps=steps,minimum_log_step=minimum_step,
        maximum_potential_phase_per_step=maxphase,
        final_central_density_contrast=final['central_density_contrast'],
        initial_mass_shell_comoving_radii=list(shell_radii),
        conserved_shell_geometric_masses=(np.array(shell_masses)*prep['mass_unit']).tolist(),
        final_physical_half_height_radius=final['physical_half_height_radius'],
        final_physical_mass_shell_radii=final['physical_mass_shell_radii'],
        maximum_physical_mass_shell_radii=np.max(positions,axis=0).tolist(),
        shell_maximum_scale_factors=[history[i]['scale_factor'] for i in maxima],
        final_over_maximum_shell_radius=(positions[-1]/np.max(positions,axis=0)).tolist(),
        final_cut_excess_norms=final['cut_excess_norms'],
        final_total_excess_norm=final['total_excess_norm'],initial_total_excess_norm=float(N0),
        initial_energy=E0,maximum_mass_flux_gap_over_initial_positive=max_flux,
        maximum_total_flux_gap_over_initial_positive=max_total,
        maximum_energy_work_gap_over_initial_energy=max_work,
        outer_net_flux_over_initial_positive=outflow/positive_norm,
        maximum_potential_over_c_squared=max(h['potential_max'] for h in history),
        final_field_sha256=digest,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    if retain:return result,dict(radius=r,initial_u=initial_u,final_u=u,history=history)
    return result
