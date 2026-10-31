"""Smooth spherical continuation of the declared one-mode valve preparation.

The action is unchanged. Spherical symmetry is exact for this NEW occupation,
not an average of the old two-mode state. Synchronous proper time and regular
pole variables are used. A future-trapped neck is not identified with a
regular child, an event horizon, or a received-energy fraction.
"""
import math
import time

import numpy as np
from scipy.integrate import solve_ivp, quad, cumulative_simpson
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer

from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_valve_neck import solve_necks


def canonical_potential(S,chi,sstar):
    """The unchanged two-canonical-field potential, shared by both frames."""
    s=sstar+S/math.sqrt(84)
    a=7*np.exp(-12*s);b=-7.5*np.exp(-14*s);c=25/32*np.exp(-24*s)
    mu2=400*np.exp(-12*s)
    return (a+b+c+.5*mu2*chi**2,
            (-12*a-14*b-24*c-6*mu2*chi**2)/math.sqrt(84),mu2*chi)


def gegenbauer_jets(u, modes):
    u=np.atleast_1d(np.asarray(u,dtype=float))
    if np.any(~np.isfinite(u)) or np.any(abs(u)>1):
        raise ValueError('finite u in [-1,1] required')
    values=[]
    for derivative,factor in enumerate((1.,2.,8.)):
        values.append(np.column_stack([
            factor*eval_gegenbauer(l-derivative,1+derivative,u)
            if l>=derivative else np.zeros_like(u) for l in range(modes)]))
    return values


class SphericalValve:
    """Exact SO(3) reduction of Einstein gravity and the two canonical fields.

    h=q² dchi²+b² sin²chi dOmega², u=cos(chi), q=b+(1-u²)c.
    K_rad=K_tan+(1-u²)d enforces regular metric/curvature poles.
    State rows are b,c,K_tan,d,S,chi,Pi_S,Pi_chi. No shell is present.
    """
    def __init__(self,modes=32,*,epsilon=1.,initial_degree=40,budget_seconds=60.):
        if isinstance(modes,bool) or not isinstance(modes,int) or not 12<=modes<=128:
            raise ValueError('integer modes in [12,128] required')
        self.modes=modes
        self.initial_data=solve_necks(initial_degree,(epsilon,))[0]
        self.parent=self.initial_data.parent;self.H=self.parent.H
        self.kappa=self.parent.kappa;self.epsilon=float(epsilon)
        self.sstar=next(x['s'] for x in VolumeParent().vacuum_roots() if x['kind']=='minimum')
        self.budget_seconds=float(budget_seconds)
        if not math.isfinite(self.budget_seconds) or self.budget_seconds<=0:
            raise ValueError('positive finite compute budget required')
        angles=np.pi*np.arange(1,modes+1)/(modes+1)
        self.u=np.cos(angles)
        Y,Yu,Yuu=gegenbauer_jets(self.u,modes)
        self.inverse=np.linalg.inv(Y);self.D=Yu@self.inverse;self.D2=Yuu@self.inverse
        self.weights=4*np.pi**2/(modes+1)*np.sin(angles)**2
        psi=self.initial_data.values(self.u)
        self.initial=np.zeros((8,modes));self.initial[0]=psi**2/self.H
        z=self.parent.mass_squared/self.H**2
        self.initial[5]=self.epsilon*math.sqrt(24/(self.kappa*(z+3)))*self.u
        self.calls=0;self.started=None

    def potential(self,S,chi):
        return canonical_potential(S,chi,self.sstar)

    def coefficients(self,state):
        state=np.asarray(state,dtype=float)
        if state.shape!=(8,self.modes) or np.any(~np.isfinite(state)):
            raise ValueError('finite 8 by modes state required')
        return state@self.inverse.T

    def jets(self,state,u=None):
        if u is None:
            values=np.asarray(state);first=(self.D@values.T).T;second=(self.D2@values.T).T
        else:
            Y,Yu,Yuu=gegenbauer_jets(u,self.modes);coef=self.coefficients(state)
            values=coef@Y.T;first=coef@Yu.T;second=coef@Yuu.T
        return values,first,second

    def pieces(self,state,u=None):
        u=self.u if u is None else np.atleast_1d(u)
        y,yu,yuu=self.jets(state,None if u is self.u else u)
        b,c,kb,d,S,chi,PS,PC=y
        bu,cu,Su,chi_u=yu[[0,1,4,5]];buu,Suu,chi_uu=yuu[[0,4,5]]
        v=1-u*u;q=b+v*c;qu=bu-2*u*c+v*cu;ka=kb+v*d;K=ka+2*kb
        if np.any(q<=0) or np.any(b<=0):
            raise ValueError('positive spatial metric ended; change foliation before continuing')
        rr=u*b-v*bu
        radial=(-b-3*u*bu+v*buu+(qu/q)*rr)/(q*q*b)
        angular=(b*b+2*b*c+2*u*b*bu+v*(c*c-bu*bu))/(q*q*b*b)
        Rr=-2*radial;Rt=angular-radial
        Rdiff=(b/q*(-2*b*c-u*c*bu-2*c*c-u*b*cu+bu*bu+v*cu*bu)-b*buu-c*c+bu*bu)/(q*q*b*b)
        gradient_square=v*(Su*Su+chi_u*chi_u)/(q*q)
        W,WS,WC=self.potential(S,chi)
        drift=(v*(2*bu/b-qu/q)-3*u)/(q*q)
        lapS=v*Suu/(q*q)+drift*Su;lapC=v*chi_uu/(q*q)+drift*chi_u
        pi_square=PS*PS+PC*PC
        rho=.5*(pi_square+gradient_square)+W
        pr=.5*(pi_square+gradient_square)-W
        pt=.5*(pi_square-gradient_square)-W
        # j is physical energy flux in the increasing-chi unit direction.
        j=np.sqrt(v)/q*(PS*Su+PC*chi_u)
        return dict(y=y,yu=yu,u=u,q=q,b=b,ka=ka,kb=kb,K=K,Rr=Rr,Rt=Rt,Rdiff=Rdiff,
                    grad2=gradient_square,W=W,WS=WS,WC=WC,lapS=lapS,lapC=lapC,
                    rho=rho,pr=pr,pt=pt,j=j,rr=rr,radius=np.sqrt(v)*b)

    def rhs(self,t,flat):
        self.calls+=1
        if self.started is not None and self.calls%50==0 and time.perf_counter()-self.started>self.budget_seconds:
            raise TimeoutError('bounded spherical evolution compute budget exceeded')
        y=np.asarray(flat).reshape(8,self.modes);z=self.pieces(y)
        b,c,kb,d,S,chi,PS,PC=y
        v=1-self.u*self.u;Su,chi_u=z['yu'][[4,5]]
        out=np.array([-kb*b,-kb*c-d*z['q'],z['Rt']+z['K']*kb-self.kappa*z['W'],
                      z['Rdiff']+z['K']*d-self.kappa*(Su*Su+chi_u*chi_u)/z['q']**2,
                      PS,PC,z['K']*PS+z['lapS']-z['WS'],z['K']*PC+z['lapC']-z['WC']])
        # Exact reflection parity; this projects roundoff only, not a physical mode.
        parity=np.array([1,1,1,1,1,-1,1,-1])[:,None]
        return (.5*(out+parity*out[:,::-1])).ravel()

    def evolve(self,endpoint=.3,*,rtol=3e-11,atol=1e-13,max_step=.002):
        if not math.isfinite(endpoint) or not 0<endpoint<=2:
            raise ValueError('finite proper endpoint in (0,2] required')
        self.started=time.perf_counter();self.calls=0
        sol=solve_ivp(self.rhs,(0.,endpoint),self.initial.ravel(),method='DOP853',rtol=rtol,
                      atol=atol,max_step=max_step,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        return sol

    def expansions(self,state,u):
        z=self.pieces(state,u)
        if np.any(abs(z['u'])>=1):raise ValueError('closed spheres require |u|<1')
        curvature=2*z['rr']/(z['q']*z['radius'])
        return -2*z['kb']+curvature,-2*z['kb']-curvature

    def diagnostic(self,t,state,*,validation_order=None):
        order=validation_order or 3*self.modes+7
        angles=np.pi*np.arange(1,order+1)/(order+1)
        u=np.cos(angles);w=4*np.pi**2/(order+1)*np.sin(angles)**2
        z=self.pieces(state,u);measure=z['q']*z['b']**2
        volume=float(w@measure)
        average=lambda f:float(w@(measure*f))/volume
        Ham=z['Rr']+2*z['Rt']+z['K']**2-z['ka']**2-2*z['kb']**2-2*self.kappa*z['rho']
        kbu=z['yu'][2];bu=z['yu'][0]
        geom=2*((u/np.sqrt(1-u*u)-np.sqrt(1-u*u)*bu/z['b'])*(z['ka']-z['kb'])+np.sqrt(1-u*u)*kbu)/z['q']
        Mom=geom-self.kappa*z['j']
        neck=self.pieces(state,[0.]);tp,tm=self.expansions(state,[0.])
        r0=float(neck['radius'][0]);b0=float(neck['b'][0]);kb0=float(neck['kb'][0])
        norm_dr=-(kb0*b0)**2 # r_chi=0 exactly on the reflection neck.
        pressure=average((z['pr']+2*z['pt'])/3);density=average(z['rho'])
        work=float(w@(measure*(z['ka']*z['pr']+2*z['kb']*z['pt'])))
        return dict(proper_time=float(t),modes=self.modes,epsilon=self.epsilon,
            theta_plus_over_H=float(tp[0]/self.H),theta_minus_over_H=float(tm[0]/self.H),
            closed_neck_future_trapped=bool(tp[0]<0 and tm[0]<0),neck_radius=r0,
            neck_radial_scale=float(neck['q'][0]),neck_rho=float(neck['rho'][0]),
            neck_radial_pressure=float(neck['pr'][0]),neck_tangential_pressure=float(neck['pt'][0]),
            neck_timelike_energy_flux=float(neck['j'][0]),
            neck_radius_gradient_norm_squared=norm_dr,
            neck_geometric_mass=0.5*r0*(1-norm_dr),
            volume=volume,half_volume=volume/2,energy=volume*density,half_energy=volume*density/2,
            mean_density=density,mean_isotropic_pressure=pressure,mean_w=pressure/density,
            pressure_work_rate=work,half_pressure_work_rate=work/2,
            Hamiltonian_RMS_over_H_squared=math.sqrt(average(Ham*Ham))/self.H**2,
            momentum_RMS_over_H_squared=math.sqrt(average(Mom*Mom))/self.H**2,
            min_radial_scale=float(min(z['q'])),min_angular_scale=float(min(z['b'])),
            initial_constraint_degree=self.initial_data.degree,validation_order=order,
            scope='Numerical finite continuation of explicitly supplied one-mode occupation; no regular child or inherited preparation selection.')

    def trapped_band(self,state):
        """Locate the first marginal sphere on each side of the neck.

        Returns the connected future-trapped band containing u=0, not an
        event horizon or a new universe. A root is accepted only after both
        expansions are checked in the band and marginality is re-evaluated.
        """
        tp,tm=self.expansions(state,[0.])
        if tp[0]>=0 or tm[0]>=0:return None
        def edge(u):
            a,b=self.expansions(state,[u]);return float(max(a[0],b[0]))
        samples=np.linspace(0.,.999,257)
        values=np.maximum(*self.expansions(state,samples))
        idx=np.flatnonzero(values>=0)
        if not len(idx):return None
        i=int(idx[0]);root=brentq(edge,float(samples[i-1]),float(samples[i]),xtol=2e-13)
        def integrands(u):
            z=self.pieces(state,[u]);measure=4*np.pi*z['q'][0]*z['b'][0]**2*math.sqrt(1-u*u)
            return (measure,measure*z['rho'][0],measure*(z['ka'][0]*z['pr'][0]+2*z['kb'][0]*z['pt'][0]),
                    measure*(z['pr'][0]+2*z['pt'][0])/3)
        values=[2*quad(lambda u:integrands(u)[i],0,root,epsabs=2e-11,epsrel=2e-11)[0] for i in range(4)]
        z=self.pieces(state,[root]);a,b=self.expansions(state,[root]);radius=float(z['radius'][0])
        excitation=values[1]-self.parent.vacuum*values[0]
        excitation_pressure=values[3]+self.parent.vacuum*values[0]
        return dict(u_edge=root,volume=values[0],energy=values[1],pressure_work_rate=values[2],
                    mean_density=values[1]/values[0],mean_isotropic_pressure=values[3]/values[0],
                    mean_w=values[3]/values[1],
                    reference_vacuum_energy=self.parent.vacuum*values[0],
                    excitation_energy_above_reference_vacuum=excitation,
                    excitation_w=excitation_pressure/excitation if abs(excitation)>1e-14 else None,
                    marginal_expansion=float(max(a[0],b[0])),other_expansion=float(min(a[0],b[0])),
                    boundary_radius=radius,boundary_rho=float(z['rho'][0]),
                    # At +u the outward direction opposes increasing chi.
                    outward_energy_flux=float(-z['j'][0]),boundary_radial_scale=float(z['q'][0]),
                    scope='Connected trapped band; moving-boundary sweeping, flux and pressure work must be distinguished.')


def conformal_forcing(y,yt,yx,sstar,kappa):
    """Wave-equation sources for e^(2a)(-dT²+dX²)+R²dOmega².

    Rows are R,a,S,chi. Null coordinates are (T-X)/sqrt(2) and
    (T+X)/sqrt(2); no membrane or incoming boundary state is added.
    """
    R,a,S,chi=y;RT,aT,ST,CT=yt;RX,aX,SX,CX=yx
    if np.min(R)<=0 or np.max(abs(a))>100:
        raise ValueError('positive finite-radius conformal patch ended')
    e2=np.exp(2*a);W,WS,WC=canonical_potential(S,chi,sstar)
    mink_R=RT*RT-RX*RX
    return np.array([-mink_R/R-e2/R+kappa*W*R*e2,
                     mink_R/(R*R)+e2/(R*R)-kappa/2*(ST*ST+CT*CT-SX*SX-CX*CX),
                     -2*(RT*ST-RX*SX)/R-e2*WS,
                     -2*(RT*CT-RX*CX)/R-e2*WC])

def conformal_diagnostics(T,tau,y,yt,yx,yxx,h,H,kappa,sstar):
    """Finite numerical controls; source-inferred curvature is conditional.

    The spatial R_ll uses the conformal slices. R_ll+K_rad*R_proper
    is the spacetime Hessian in their radial unit direction.
    """
    # All jets share the current time slice and shrunken causal patch.
    R,a,S,chi=y;RT,aT,ST,CT=yt;RX,aX,SX,CX=yx
    W,WS,WC=canonical_potential(S,chi,sstar);em=np.exp(-a);em2=em*em
    Ka=-em*aT;Kb=-em*RT/R
    radial=em2*(yxx[0]-aX*RX)/R
    Rr=-2*radial;Rt=(1-em2*RX*RX)/(R*R)-radial
    pi2=em2*(ST*ST+CT*CT);grad2=em2*(SX*SX+CX*CX)
    rho=.5*(pi2+grad2)+W
    ham=Rr+2*Rt+4*Ka*Kb+2*Kb*Kb-2*kappa*rho
    # Differentiate the same centered time jet in space for momentum.
    Rtx=(RT[2:]-RT[:-2])/(2*h)
    Kbx=-em[1:-1]*(Rtx/R[1:-1]-RT[1:-1]*RX[1:-1]/R[1:-1]**2-aX[1:-1]*RT[1:-1]/R[1:-1])
    j=-em2*(ST*SX+CT*CX)
    mom=2*em[1:-1]*RX[1:-1]/R[1:-1]*(Ka[1:-1]-Kb[1:-1])-2*em[1:-1]*Kbx-kappa*j[1:-1]
    Rproper=em*RT;Gamma=em*RX;massratio=(1+Rproper**2-Gamma**2)/(2*R*R)
    ricci=kappa*(-pi2+grad2+4*W)
    mss=em2*(-ST*ST+SX*SX);mcc=em2*(-CT*CT+CX*CX);msc=em2*(-ST*CT+SX*CX)
    rsq=kappa*kappa*(mss*mss+mcc*mcc+2*msc*msc+2*W*(-pi2+grad2)+4*W*W)
    Weyl=massratio-kappa*(pi2-grad2+2*W)/12
    kretsch=48*Weyl*Weyl+2*rsq-ricci*ricci/3
    center=len(R)//2
    F=conformal_forcing(y,yt,yx,sstar,kappa)
    proper_acc=em2*(yxx[0]+F[0]-aT*RT)
    tplus=2*(Rproper+Gamma)/R;tminus=2*(Rproper-Gamma)/R
    weight=np.exp(a)*R*R
    rms=lambda f:math.sqrt(float(np.sum(weight*f*f)/np.sum(weight)))
    # Also report a central half-patch norm, not just a moving full-patch norm.
    sl=slice(len(R)//4,3*len(R)//4+1)
    wc=weight[sl];hc=ham[sl]
    budget=abs(Rr)+2*abs(Rt)+4*abs(Ka*Kb)+2*Kb*Kb+2*kappa*abs(rho)
    return dict(conformal_time=float(T),center_proper_time=float(tau),remaining_points=len(R),
                remaining_half_width=float(h*(len(R)-1)/2),center_radius=float(R[center]),
                center_a=float(a[center]),center_lapse=float(np.exp(a[center])),
                center_S=float(S[center]),center_chi=float(chi[center]),
                center_Pi_S=float(em[center]*ST[center]),center_Pi_chi=float(em[center]*CT[center]),
                center_Rproper=float(Rproper[center]),center_Rproper_acc=float(proper_acc[center]),
                center_Rll_over_R_H2=float(radial[center]/H**2),
                center_K_rad=float(Ka[center]),
                center_spacetime_Hess_ss_over_R_H2=float((radial[center]+Ka[center]*Rproper[center]/R[center])/H**2),
                center_theta_plus_over_H=float(tplus[center]/H),center_theta_minus_over_H=float(tminus[center]/H),
                center_theta_product_over_H2=float(tplus[center]*tminus[center]/H**2),
                center_R_gradient_norm_squared=float(-Rproper[center]**2+Gamma[center]**2),
                center_energy_flux=float(j[center]),
                Hamiltonian_RMS_over_H2=rms(ham)/H**2,
                Hamiltonian_relative_RMS_term_budget=rms(ham)/rms(budget),
                Hamiltonian_relative_RMS_Riemann_scale=math.sqrt(float(np.sum(weight*ham*ham)/np.sum(weight*abs(kretsch)))),
                central_half_Hamiltonian_RMS_over_H2=math.sqrt(float(np.sum(wc*hc*hc)/np.sum(wc)))/H**2,
                momentum_RMS_over_H2=float(np.sqrt(np.mean(mom*mom))/H**2),
                max_abs_Ricci_over_H2=float(np.max(abs(ricci))/H**2),
                max_abs_source_inferred_Kretschmann_over_H4=float(np.max(abs(kretsch))/H**4),
                min_radius=float(np.min(R)),min_lapse=float(np.min(np.exp(a))),
                maximum_metric_log=float(np.max(abs(a))),
                center_rho=float(rho[center]),
                center_radial_pressure=float(.5*(pi2[center]+grad2[center])-W[center]),
                center_tangential_pressure=float(.5*(pi2[center]-grad2[center])-W[center]),
                patch_max_theta_plus_over_H=float(np.max(tplus)/H),
                patch_max_theta_minus_over_H=float(np.max(tminus)/H),
                full_current_patch_future_trapped=bool(np.all(tplus<0)&np.all(tminus<0)),
                scope='Causal shrinking patch; source-inferred curvature conditional on constraints; Rll belongs to conformal slicing; center proper-time invariants can be compared across frames.')


CONFORMAL_CENTER_COLUMNS=(
    'conformal_time','center_proper_time','center_radius','center_S',
    'center_Pi_S','center_theta_plus_over_H','center_Rproper')


def _conformal_points(points):
    if isinstance(points,(bool,np.bool_)) or not isinstance(points,(int,np.integer)) or points<33 or points%2!=1:
        raise ValueError('odd integer points >=33 required')
    return int(points)


def _conformal_evolve(initial,*,H,kappa,sstar,proper_time,endpoint,sample_every=.05,budget_seconds=60.):
    """Centered characteristic diamonds on a shrinking Cauchy domain.

    Every accepted level discards one point at either end. Spatial and
    temporal centered jets use an extra surrounding diamond; no ghost,
    reflecting, vacuum, or independently specified boundary values exist.
    """
    started=time.perf_counter()
    endpoint=float(endpoint);budget_seconds=float(budget_seconds);sample_every=float(sample_every)
    if not math.isfinite(endpoint) or not 0<endpoint<initial['half_width']:
        raise ValueError('positive endpoint strictly inside the initial causal half-width required')
    if not math.isfinite(budget_seconds) or budget_seconds<=0 or not math.isfinite(sample_every) or sample_every<=0:
        raise ValueError('positive finite budget and sample interval required')
    y0,v0,x0,xx0=(np.asarray(initial[k],dtype=float) for k in ('y','yt','yx','yxx'))
    n=_conformal_points(y0.shape[1]);h=float(initial['step'])
    if any(a.shape!=(4,n) or np.any(~np.isfinite(a)) for a in (y0,v0,x0,xx0)) or h<=0:
        raise ValueError('finite four-field uniform Cauchy jets required')
    if np.min(y0[0])<=0:raise ValueError('strictly positive areal radius on the patch required')
    if endpoint<h:raise ValueError('endpoint must include at least one actual diamond step')
    first=conformal_diagnostics(0.,proper_time,y0,v0,x0,xx0,h,H,kappa,sstar)
    rows=[first];history=[[first[k] for k in CONFORMAL_CENTER_COLUMNS]]
    acceleration=xx0+conformal_forcing(y0,v0,x0,sstar,kappa)
    prev=y0.copy();cur=(y0+h*v0+.5*h*h*acceleration)[:,1:-1]
    tau=proper_time+h*(1+math.exp(float(cur[1,len(cur[0])//2])))/2
    prior_a=float(cur[1,len(cur[0])//2]);check=sample_every;step=1
    parity=np.array([1,1,1,-1])[:,None]
    events=[];neck_minimum=True;prior_Rproper=float(v0[0,n//2])
    status='running';message=None;max_iterations=0;profile=None
    try:
        while cur.shape[1]>=25:
            T=step*h
            if T>endpoint+32*np.finfo(float).eps*max(1.,endpoint):
                status='requested_endpoint';break
            E,W,S=cur[:,2:],cur[:,:-2],prev[:,2:-2];C=cur[:,1:-1]
            dx=(E-W)/(2*h);dxx=(E-2*C+W)/(h*h)
            guess=E+W-S+h*h*conformal_forcing(C,(C-S)/h,dx,sstar,kappa)
            for it in range(16):
                updated=E+W-S+h*h*conformal_forcing((guess+E+W+S)/4,(guess-S)/(2*h),dx,sstar,kappa)
                error=float(np.max(abs(updated-guess)/(1+abs(updated))))
                guess=updated
                if error<3e-13:break
            else:raise RuntimeError('implicit diamond iteration did not converge')
            max_iterations=max(max_iterations,it+1)
            # The supplied one-mode state has this exact reflection symmetry.
            guess=.5*(guess+parity*guess[:,::-1]);dt=(guess-S)/(2*h)
            center=len(C[0])//2;em=math.exp(-float(C[1,center]))
            Rproper=em*float(dt[0,center]);Rll=em*em*float(dxx[0,center])
            history.append([T,tau,float(C[0,center]),float(C[2,center]),em*float(dt[2,center]),
                            2*Rproper/(float(C[0,center])*H),Rproper])
            profile=dict(T=T,center_proper_time=tau,x=h*np.arange(-center,center+1),
                         y=C.copy(),yt=dt.copy(),yx=dx.copy(),yxx=dxx.copy())
            if neck_minimum and Rll<=0:
                events.append(dict(kind='spatial_neck_type_change_on_conformal_slices',conformal_time=T,proper_time=tau))
                neck_minimum=False
            if prior_Rproper<0 and Rproper>=0:
                events.append(dict(kind='center_proper_time_bounce',conformal_time=T,proper_time=tau))
                status='temporal_bounce'
            prior_Rproper=Rproper
            final=T+h>endpoint+32*np.finfo(float).eps*max(1.,endpoint)
            if T>=check or final or status!='running':
                row=conformal_diagnostics(T,tau,C,dt,dx,dxx,h,H,kappa,sstar)
                rows.append(row);check+=sample_every
                if row['central_half_Hamiltonian_RMS_over_H2']>.1:
                    status='numerical_validity_stop';message='central half-patch Hamiltonian RMS/H²>.1; not a physical endpoint'
                if row['min_radius']<.15:
                    status='radius_safety_stop';message='sampled areal radius<.15; no continuation past this guard'
            if status!='running':break
            prev,cur=cur,guess;step+=1
            current_a=float(cur[1,len(cur[0])//2])
            tau+=h*(math.exp(prior_a)+math.exp(current_a))/2;prior_a=current_a
            if time.perf_counter()-started>budget_seconds:
                status='compute_budget';message='bounded conformal evolution budget exceeded';break
        else:status='causal_patch_end'
    except (ValueError,RuntimeError,FloatingPointError,OverflowError) as exc:
        status='numerical_failure';message=str(exc)
    return dict(points=n,step=h,initial_half_width=float(initial['half_width']),
                requested_endpoint=endpoint,status=status,message=message,rows=rows,
                center_history=np.asarray(history),center_columns=CONFORMAL_CENTER_COLUMNS,
                end_profile=profile,events=events,maximum_implicit_iterations=max_iterations,
                elapsed_seconds=time.perf_counter()-started,
                boundary='Shrinking characteristic Cauchy diamond; no supplied boundary data.',
                scope='Finite numerical continuation with empirical errors, not a continuum enclosure, physical singularity, or regular-child proof.')


class ConformalValve:
    """Local continuation of an in-memory reflection-symmetric Cauchy state.

    No filesystem, cache, source selection, or new incoming state belongs
    to this class. Provenance and regeneration are owned by the producer.
    """
    def __init__(self,system,state,*,initial_proper_time,u_max=.9,mapping_order=8193):
        if not isinstance(system,SphericalValve):raise TypeError('SphericalValve system required')
        state=np.asarray(state,dtype=float).copy();system.coefficients(state)
        parity=np.array([1,1,1,1,1,-1,1,-1])[:,None]
        if np.max(abs(state-parity*state[:,::-1]))>1e-10*(1+np.max(abs(state))):
            raise ValueError('the declared reflection-symmetric Cauchy state is required')
        u_max=float(u_max);initial_proper_time=float(initial_proper_time)
        if not math.isfinite(u_max) or not 0<u_max<1 or not math.isfinite(initial_proper_time):
            raise ValueError('finite proper time and u_max strictly between0 and1 required')
        if isinstance(mapping_order,bool) or not isinstance(mapping_order,int) or mapping_order<65 or mapping_order%2!=1:
            raise ValueError('odd proper-radial mapping order >=65 required')
        self.system=system;self.state=state;self.state.setflags(write=False)
        self.initial_proper_time=initial_proper_time;self.u_max=u_max
        uu=np.linspace(0.,u_max,mapping_order);z=system.pieces(state,uu)
        proper=cumulative_simpson(z['q']/np.sqrt(1-uu*uu),x=uu,initial=0.)
        self.half_width=float(proper[-1]);self._inverse=PchipInterpolator(proper,uu)

    def initial_grid(self,points):
        points=_conformal_points(points);X=np.linspace(-self.half_width,self.half_width,points)
        u=-np.sign(X)*self._inverse(abs(X));u[points//2]=0.
        z=self.system.pieces(self.state,u);q=z['q'];R=z['radius'];v=1-u*u
        y=np.array([R,np.zeros(points),z['y'][4],z['y'][5]])
        yt=np.array([-z['kb']*R,-z['ka'],z['y'][6],z['y'][7]])
        yx=np.array([z['rr']/q,np.zeros(points),-np.sqrt(v)/q*z['yu'][4],-np.sqrt(v)/q*z['yu'][5]])
        yxx=np.array([-z['Rr']*R/2,np.zeros(points),z['lapS']-2*yx[0]*yx[2]/R,z['lapC']-2*yx[0]*yx[3]/R])
        return dict(x=X,u=u,y=y,yt=yt,yx=yx,yxx=yxx,step=float(X[1]-X[0]),half_width=self.half_width)

    def evolve(self,points,endpoint,*,sample_every=.05,budget_seconds=60.):
        result=_conformal_evolve(self.initial_grid(points),H=self.system.H,kappa=self.system.kappa,
                                 sstar=self.system.sstar,proper_time=self.initial_proper_time,
                                 endpoint=endpoint,sample_every=sample_every,budget_seconds=budget_seconds)
        result['initial_u_max']=self.u_max
        return result


def conformal_de_sitter_control(points,*,endpoint=.5,system=None):
    """Analytic same-action closed de Sitter control, with no incoming data."""
    points=_conformal_points(points);system=system or SphericalValve(24,epsilon=0.)
    H=system.H;x=np.linspace(-1.,1.,points)
    y=np.zeros((4,points));y[0]=np.cos(H*x)/H
    yt=np.zeros_like(y);yx=np.zeros_like(y);yx[0]=-np.sin(H*x)
    yxx=np.zeros_like(y);yxx[0]=-H*np.cos(H*x)
    initial=dict(x=x,y=y,yt=yt,yx=yx,yxx=yxx,step=float(x[1]-x[0]),half_width=1.)
    result=_conformal_evolve(initial,H=H,kappa=system.kappa,sstar=system.sstar,
                             proper_time=0.,endpoint=endpoint)
    profile=result['end_profile']
    if profile is None:raise RuntimeError('de Sitter control did not produce a finite level')
    T=profile['T'];exactR=np.cos(H*profile['x'])/(H*np.cos(H*T));exacta=-math.log(math.cos(H*T))
    return dict(points=points,conformal_time=T,status=result['status'],
                max_radius_error=float(np.max(abs(profile['y'][0]-exactR))),
                max_metric_log_error=float(np.max(abs(profile['y'][1]-exacta))),
                Hamiltonian_RMS_over_H2=result['rows'][-1]['Hamiltonian_RMS_over_H2'],
                momentum_RMS_over_H2=result['rows'][-1]['momentum_RMS_over_H2'],
                exact='a=-log cos(H T), R=cos(H X)/(H cos(H T)), S=chi=0',
                scope='Field convergence control; small constraint residuals are reported including nonmonotone finite-difference noise.')


def align_conformal_center(result,proper_time):
    """Align a center observation; quantify full-vs-decimated PCHIP change.

    The change is an empirical interpolation indicator, not an error bound.
    No extrapolation outside the actually stored proper-time history is allowed.
    """
    trace=np.asarray(result['center_history']);columns=tuple(result['center_columns'])
    clock=trace[:,columns.index('center_proper_time')];proper_time=float(proper_time)
    if not math.isfinite(proper_time) or not clock[0]<=proper_time<=clock[-1]:
        raise ValueError('proper time lies outside the actual center history')
    if len(clock)<4:raise ValueError('at least four actual center samples required')
    index=np.unique(np.r_[np.arange(0,len(clock),2),len(clock)-1])
    values={};indicators={}
    for i,key in enumerate(columns):
        full=float(PchipInterpolator(clock,trace[:,i])(proper_time))
        coarse=float(PchipInterpolator(clock[index],trace[index,i])(proper_time))
        values[key]=full;indicators[key]=abs(full-coarse)
    return dict(values=values,interpolation_indicator=indicators,
                meaning='Full-vs-every-other actual-step PCHIP difference; empirical, not certified.')
