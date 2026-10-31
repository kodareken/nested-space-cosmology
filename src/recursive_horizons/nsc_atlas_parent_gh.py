"""Nonlinear covariant-GH evolution of the same nonspherical compact parent.

The prescribed source and reduced two-field action are unchanged. A regular
S3 frame and its exact preserved proper-flip symmetries reduce computation;
there is no spherical averaging or suppression of physical modes.
"""
from dataclasses import dataclass
import math
import time

import numpy as np
from numba import set_num_threads
from scipy.integrate import DOP853
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from .nsc_atlas_finite_parent_data import solve_initial_data,source_profile
from .nsc_atlas_gh_constraints import adm_constraints
from .nsc_atlas_gh_pointwise import rhs_nodes
from .nsc_atlas_gh_observer import observer_readout
from .nsc_atlas_gh_matter_work import matter_work
from .nsc_atlas_passive_dop853 import PassiveDOP853Quadrature
from .nsc_atlas_parent_backreaction import anchor_point
from .nsc_atlas_parent_circulation import Parent
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_s3_frame import build_frame_basis,point_values

PAIRS=tuple((a,b) for a in range(4) for b in range(a,4))
FRAME_CHAR=np.array([[-1,1],[-1,-1],[1,-1]])
FRAME_MATRICES=np.array([
    [[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]],
    [[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]],
    [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]],dtype=float)


def epsilon3(i,j,k):
    if len({i,j,k})<3:return 0
    return 1 if (i,j,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1


def reference_connection(r):
    Gamma=np.zeros((4,4,4))
    for a in range(3):
        for i in range(3):
            for b in range(3):Gamma[a+1,i+1,b+1]=r*epsilon3(i,b,a)
    return Gamma


def difference_connection(g,v,dg,Gamma):
    D=np.empty((4,4,4));D[0]=v
    for i in range(3):
        D[i+1]=dg[i]-Gamma[:,i+1,:].T@g-g@Gamma[:,i+1,:]
    lower=np.empty((4,4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):lower[a,b,c]=(D[b,a,c]+D[c,a,b]-D[a,b,c])/2
    return np.einsum('ad,dbc->abc',np.linalg.inv(g),lower)


@dataclass
class Jets:
    g: np.ndarray
    v: np.ndarray
    dg: np.ndarray
    d2g: np.ndarray
    dv: np.ndarray
    q: np.ndarray
    p: np.ndarray
    dq: np.ndarray
    d2q: np.ndarray
    dp: np.ndarray


class ParentGH:
    def __init__(self,degree,epsilon,*,quadrature_order=None,gauge='lapse_controlled',backend='hopf',kick_amplitude=0.,kick_source=1):
        if degree<2 or degree>32 or not 0<=epsilon<=1:
            raise ValueError('2<=degree<=32 and 0<=epsilon<=1 required')
        self.parent=Parent.retained();self.epsilon=epsilon
        self.kick_amplitude=float(kick_amplitude);self.kick_source=kick_source;self.initial_kick_metadata=None
        if not math.isfinite(self.kick_amplitude) or abs(self.kick_amplitude)>1 or isinstance(kick_source,bool) or kick_source not in (1,2):
            raise ValueError('finite |kick_amplitude|<=1 and kick_source1 or2 required')
        if gauge not in ('harmonic','vacuum_reference','lapse_controlled'):raise ValueError('unknown prescribed GH gauge')
        self.gauge=gauge
        model=VolumeParent();self.minimum=next(v for v in model.vacuum_roots() if v['kind']=='minimum')
        self.sstar=self.minimum['s'];self.K=model.K;self.r=self.parent.H
        if backend not in ('hopf','dense'):raise ValueError('unknown spatial transform backend')
        self.backend=backend
        if backend=='hopf':
            from .nsc_atlas_s3_hopf import build_hopf_frame_basis
            order=(math.ceil(.75*degree)+2,2*math.ceil((3*degree+1)/2)) if quadrature_order is None else quadrature_order
            self.basis=build_hopf_frame_basis(degree,order)
        else:
            n=math.ceil(1.5*degree)+3 if quadrature_order is None else quadrature_order
            self.basis=build_frame_basis(degree,n)
        self.N=len(self.basis.points)
        self.Gamma=reference_connection(self.r);self.calls=0
        spacetime=np.r_[[[1,1]],FRAME_CHAR]
        metric_chars=np.array([spacetime[a]*spacetime[b] for a,b in PAIRS])
        self.chars=np.r_[metric_chars,metric_chars,[[1,1],[-1,-1],[1,1],[-1,-1]]]
        self.sectors={}
        for char in ((1,1),(1,-1),(-1,1),(-1,-1)):
            idx=np.flatnonzero(np.all(self.basis.characters==char,axis=1))
            Y=None if backend=='hopf' else np.vstack([self.basis.values[:,idx],np.zeros((1,len(idx)))])
            self.sectors[char]=(idx,Y)
        self.indices=[self.sectors[tuple(char)][0] for char in self.chars]
        self.offsets=np.r_[0,np.cumsum([len(i) for i in self.indices])]
        self.dof=int(self.offsets[-1]);self.started=time.perf_counter();self.budget=math.inf
        self.invalid_stages=0
        self.passive_work=None;self.work_valid=True

    def pack(self,C):return np.concatenate([C[i,idx] for i,idx in enumerate(self.indices)])

    def unpack(self,y):
        C=np.zeros((24,len(self.basis.labels)))
        for i,idx in enumerate(self.indices):C[i,idx]=y[self.offsets[i]:self.offsets[i+1]]
        return C

    def update_point(self,X):
        values=point_values(X,self.basis.labels)
        self.point_row=values
        if self.backend=='dense':
            for idx,Y in self.sectors.values():Y[-1]=values[idx]

    def evaluate(self,C,chars):
        if self.backend=='hopf':return np.vstack([self.basis.evaluate(C.T),self.point_row@C.T])
        output=np.empty((self.N+1,len(C)))
        for char,(idx,Y) in self.sectors.items():
            fields=np.flatnonzero(np.all(chars==char,axis=1))
            if len(fields):output[:,fields]=Y@C[np.ix_(fields,idx)].T
        return output

    def project(self,values,chars):
        if self.backend=='hopf':
            C=self.basis.project(values[:self.N]).T
            for i,char in enumerate(chars):C[i,~np.all(self.basis.characters==char,axis=1)]=0.
            return C
        C=np.zeros((values.shape[1],len(self.basis.labels)))
        for char,(idx,Y) in self.sectors.items():
            fields=np.flatnonzero(np.all(chars==char,axis=1))
            if len(fields):C[np.ix_(fields,idx)]=(Y[:-1].T@(self.basis.weights[:,None]*values[:self.N,fields])).T
        return C

    def initial(self):
        self.X0=anchor_point(self.parent)
        self.update_point(self.X0)
        if self.kick_amplitude:
            from .nsc_atlas_gravitating_kick import build_kicked_gh_state
            return build_kicked_gh_state(self)
        C=np.zeros((24,len(self.basis.labels)))
        C[0,0]=-math.sqrt(2*math.pi**2)
        if not self.epsilon:
            for k,(a,b) in enumerate(PAIRS):
                if a==b and a:C[k,0]=math.sqrt(2*math.pi**2)
            return np.r_[self.pack(C),self.X0,[1.,0.,0.,0.],0.,0.]
        if self.backend=='hopf':
            from .nsc_atlas_s3_hopf import build_hopf_frame_basis
            initial_degree=max(8,min(self.basis.degree,20))
            data=solve_initial_data(initial_degree,epsilon=self.epsilon)
            D=4*initial_degree+self.basis.degree
            initial_basis=build_hopf_frame_basis(self.basis.degree,(math.ceil(D/4)+2,2*math.ceil((D+1)/2)))
            psi=np.empty(len(initial_basis.points))
            for start in range(0,len(psi),8192):psi[start:start+8192]=data.values(initial_basis.points[start:start+8192])
            Fcoef=initial_basis.project(psi**4)
            Fcoef[~np.all(self.basis.characters==[1,1],axis=1)]=0.
            for k,(a,b) in enumerate(PAIRS):
                if a==b and a:C[k]=Fcoef
            f,_=source_profile(initial_basis.points,self.parent)
            C[21]=self.epsilon*initial_basis.project(f)
            C[21,~np.all(self.basis.characters==[-1,-1],axis=1)]=0.
            F=initial_basis.evaluate(Fcoef)
            if min(F)<=0:raise ValueError('projected initial spatial metric lost positivity')
            for i,Di in enumerate(self.basis.generators):
                velocity=-.5*initial_basis.evaluate(self.r*(Di@Fcoef))/F
                k=10+PAIRS.index((0,i+1));C[k]=initial_basis.project(velocity)
                C[k,~np.all(self.basis.characters==self.chars[k],axis=1)]=0.
            return np.r_[self.pack(C),self.X0,[1.,0.,0.,0.],0.,0.]
        if self.epsilon:
            initial=solve_initial_data(max(8,self.basis.degree),epsilon=self.epsilon)
            psi=initial.values(self.basis.points)
        else:psi=np.ones(self.N)
        f,_=source_profile(self.basis.points,self.parent)
        vals=np.zeros((self.N,24))
        for k,(a,b) in enumerate(PAIRS):
            if a==b and a:vals[:,k]=psi**4
        vals[:,21]=self.epsilon*f
        projected=self.project(vals,self.chars);projected[0,0]=-math.sqrt(2*math.pi**2)
        C[:]=projected
        Fcoef=C[PAIRS.index((1,1))]
        F=self.evaluate(Fcoef[None,:],np.array([[1,1]]))[:,0]
        for i,D in enumerate(self.basis.generators):
            grad=self.evaluate((self.r*(D@Fcoef))[None,:],FRAME_CHAR[i:i+1])[:,0]
            k=PAIRS.index((0,i+1));values=np.zeros((self.N+1,1))
            values[:,0]=-.5*grad/F
            C[k+10]=self.project(values,self.chars[k+10:k+11])[0]
        obs=np.r_[self.X0,[1.,0.,0.,0.],0.,0.]
        return np.r_[self.pack(C),obs]

    def jets(self,y):
        C=self.unpack(y);self.update_point(y[self.dof:self.dof+4])
        V=self.evaluate(C,self.chars)
        first=[]
        for i,D in enumerate(self.basis.generators):
            first.append(self.evaluate(self.r*(D@C.T).T,self.chars*FRAME_CHAR[i]))
        selected=np.r_[np.arange(10),[20,21]];Cs=C[selected];chars=self.chars[selected]
        second=np.empty((3,3,self.N+1,12))
        for i,Di in enumerate(self.basis.generators):
            for j in range(i,3):
                Dj=self.basis.generators[j]
                val=self.evaluate(self.r**2*(Di@(Dj@Cs.T)).T,chars*FRAME_CHAR[i]*FRAME_CHAR[j])
                second[i,j]=val
                if i!=j:
                    k=3-i-j
                    second[j,i]=val-2*self.r*epsilon3(i,j,k)*first[k][:,selected]
        N=self.N+1;g=np.empty((N,4,4));v=np.empty_like(g)
        dg=np.empty((N,3,4,4));dv=np.empty_like(dg);d2g=np.empty((N,3,3,4,4))
        for k,(a,b) in enumerate(PAIRS):
            g[:,a,b]=g[:,b,a]=V[:,k];v[:,a,b]=v[:,b,a]=V[:,k+10]
            for i in range(3):
                dg[:,i,a,b]=dg[:,i,b,a]=first[i][:,k]
                dv[:,i,a,b]=dv[:,i,b,a]=first[i][:,k+10]
                for j in range(3):d2g[:,i,j,a,b]=d2g[:,i,j,b,a]=second[i,j,:,k]
        dq=np.stack([f[:,20:22] for f in first],axis=1)
        dp=np.stack([f[:,22:24] for f in first],axis=1)
        d2q=np.empty((N,3,3,2))
        for i in range(3):
            for j in range(3):d2q[:,i,j]=second[i,j,:,10:12]
        return Jets(g,v,dg,d2g,dv,V[:,20:22],V[:,22:24],dq,d2q,dp)

    def potential(self,q):
        s=self.sstar+q[:,0]/math.sqrt(self.K)
        return 7*np.exp(-12*s)-7.5*np.exp(-14*s)+25/32*np.exp(-24*s)+200*np.exp(-12*s)*q[:,1]**2

    def observer_angle(self,y):
        C=self.unpack(y);row=point_values(y[self.dof:self.dof+4],self.basis.labels)
        V=C@row;g=np.empty((4,4));dq=np.empty((3,2))
        for k,(a,b) in enumerate(PAIRS):g[a,b]=g[b,a]=V[k]
        for i,D in enumerate(self.basis.generators):dq[i]=self.r*row@(D@C[20:22].T)
        reading=observer_readout(g,V[20:22],V[22:24],dq,y[self.dof+4:self.dof+8],
            potential=self.potential(V[20:22][None,:])[0],omega_ref=self.parent.reference_frequency)
        return float(reading['phase_wrapped']) if reading['defined'] else None

    def gauge_source(self,t):
        if self.gauge=='harmonic':return 0.,0.
        return 3*self.r*math.tanh(self.r*t),3*self.r**2/math.cosh(self.r*t)**2

    def observer_rhs(self,y,J):
        obs=y[self.dof:];X=obs[:4]/np.linalg.norm(obs[:4]);u=obs[4:8]
        g=J.g[-1];Cup=difference_connection(g,J.v[-1],J.dg[-1],self.Gamma)
        dX=self.r*sum(u[i+1]/u[0]*(FRAME_MATRICES[i]@X) for i in range(3))
        du=-np.einsum('abc,b,c->a',Cup,u,u)/u[0]
        # Phase is read from the resolved field at accepted steps. An
        # auxiliary angular variable must not control the PDE time step.
        return np.r_[dX,du,1/u[0],0.]

    def rhs(self,t,y):
        self.calls+=1
        if self.calls%50==0 and time.perf_counter()-self.started>self.budget:
            raise TimeoutError('nonlinear GH compute budget exceeded')
        if not np.all(np.isfinite(y)):
            self.invalid_stages+=1
            if self.passive_work is not None:self.passive_work.record(0.)
            return np.full_like(y,np.nan)
        J=self.jets(y);H0,Hdot=self.gauge_source(t)
        accg,accq,C,valid=rhs_nodes(J.g,J.v,J.dg,J.d2g,J.dv,J.q,J.p,J.dq,J.d2q,J.dp,
            r=self.r,Hsource0=H0,HsourceDot0=Hdot,
            sstar=self.sstar,kappa=self.parent.kappa,gamma0=self.r,
            lapse_feedback=self.r if self.gauge=='lapse_controlled' else 0.)
        if not np.all(valid):
            # A rejected intermediate RK stage is not an accepted physical
            # state. Let the adaptive integrator reduce its step first.
            self.invalid_stages+=1
            if self.passive_work is not None:self.passive_work.record(0.)
            return np.full_like(y,np.nan)
        vals=np.empty((self.N+1,24))
        for k,(a,b) in enumerate(PAIRS):vals[:,k]=J.v[:,a,b];vals[:,k+10]=accg[:,a,b]
        vals[:,20:22]=J.p;vals[:,22:24]=accq
        derivative=np.r_[self.pack(self.project(vals,self.chars)),self.observer_rhs(y,J)]
        if self.passive_work is not None:
            try:
                work=matter_work(J.g,J.v,J.dg,J.q,J.p,J.dq,r=self.r,potential=self.potential(J.q))
                integrand=work['volume_weight']*(work['work_integrand']
                    +work['lapse']*self.parent.vacuum*work['extrinsic_trace'])
                value=float(self.basis.weights@integrand[:self.N])/(self.parent.vacuum*2*math.pi**2)
                if not math.isfinite(value):raise ValueError('nonfinite passive work')
            except ValueError:
                self.work_valid=False;value=0.
            self.passive_work.record(value)
        return derivative

    def diagnostic(self,t,y):
        J=self.jets(y);W=self.potential(J.q)
        D=adm_constraints(J.g,J.v,J.dg,J.d2g,J.dv,J.q,J.p,J.dq,r=self.r,
                          kappa=self.parent.kappa,potential=W)
        H0,Hdot=self.gauge_source(t)
        _,accq,C,valid=rhs_nodes(J.g,J.v,J.dg,J.d2g,J.dv,J.q,J.p,J.dq,J.d2q,J.dp,
            r=self.r,Hsource0=H0,HsourceDot0=Hdot,
            sstar=self.sstar,kappa=self.parent.kappa,gamma0=self.r,
            lapse_feedback=self.r if self.gauge=='lapse_controlled' else 0.)
        spatial=J.g[:self.N,1:,1:];measure=np.sqrt(np.linalg.det(spatial))*self.basis.weights
        V=float(sum(measure));rho=D['rho'][:self.N];hi=np.linalg.inv(spatial)
        nC=(C[:self.N,0]-np.sum(D['shift'][:self.N]*C[:self.N,1:],axis=1))/D['lapse'][:self.N]
        normC=nC*nC+np.einsum('ni,nij,nj->n',C[:self.N,1:],hi,C[:self.N,1:])
        M=D['M'][:self.N];normM=np.einsum('ni,nij,nj->n',M,hi,M)
        u=y[self.dof+4:self.dof+8]
        reading=observer_readout(J.g[-1],J.q[-1],J.p[-1],J.dq[-1],u,
            potential=W[-1],omega_ref=self.parent.reference_frequency)
        coeff=self.unpack(y);ell=np.array([label[0] for label in self.basis.labels])
        high=ell>=max(1,self.basis.degree-1)
        power=coeff[21]**2+coeff[23]**2/self.parent.mass_squared
        volume_power=coeff[20]**2+coeff[22]**2/(self.minimum['stiffness']/self.K)
        aq_coeff=self.project(accq,self.chars[20:22]);resolved=self.evaluate(aq_coeff,self.chars[20:22])
        rhs_fractions=[]
        for k in range(2):
            num=float(self.basis.weights@((accq[:self.N,k]-resolved[:self.N,k])**2))
            den=float(self.basis.weights@(accq[:self.N,k]**2))
            rhs_fractions.append(math.sqrt(num/den) if den>1e-28 else 0.)
        return dict(coordinate_t=float(t),observer_proper_duration=float(y[self.dof+8]),
            observer_sampled_unwrapped_phase=float(y[self.dof+9]),
            observer_phase_wrapped=float(reading['phase_wrapped']) if reading['defined'] else None,
            observer_phase_winding_certified=False,
            observer_unit_norm_error=float(reading['input_norm_error']),
            observer_chi=float(J.q[-1,1]),
            observer_chi_proper_derivative=float(reading['chi_proper_derivative']),
            observer_volume_mode=float(self.sstar+J.q[-1,0]/math.sqrt(self.K)),
            geodesic_observer_total_density_over_initial_vacuum=float(reading['density_total']/self.parent.vacuum),
            slice_normal_density_at_observer_position_over_initial_vacuum=float(D['rho'][-1]/self.parent.vacuum),
            scalar_last_two_shell_power_fraction=float(sum(power[high])/sum(power)) if sum(power)>1e-28 else 0.,
            volume_last_two_shell_power_fraction=float(sum(volume_power[high])/sum(volume_power)) if sum(volume_power)>1e-28 else 0.,
            scalar_acceleration_unresolved_fraction=rhs_fractions[1],
            volume_acceleration_unresolved_fraction=rhs_fractions[0],
            slice_proper_volume_over_initial_round_reference=V/(2*math.pi**2),
            excitation_normal_energy_over_initial_round_vacuum_energy=float(measure@(rho-self.parent.vacuum))/(self.parent.vacuum*2*math.pi**2),
            slice_mean_normal_total_density_over_initial_vacuum=float(measure@rho/V/self.parent.vacuum),
            slice_max_sampled_normal_total_density_over_initial_vacuum=float(max(rho)/self.parent.vacuum),
            maximum_to_mean_slice_normal_density_ratio=float(max(rho)/(measure@rho/V)),
            volume_mode_min=float(min(self.sstar+J.q[:self.N,0]/math.sqrt(self.K))),
            volume_mode_max=float(max(self.sstar+J.q[:self.N,0]/math.sqrt(self.K))),
            ADM_H_RMS_over_H_squared=math.sqrt(float(measure@(D['H'][:self.N]**2)/V))/self.r**2,
            ADM_M_RMS_over_H_squared=math.sqrt(max(0.,float(measure@normM/V)))/self.r**2,
            GH_C_RMS_over_H=math.sqrt(max(0.,float(measure@normC/V)))/self.r,
            minimum_lapse=float(min(D['lapse'])),maximum_lapse=float(max(D['lapse'])),metric_valid=bool(np.all(valid)))


def evolve_parent(*,degree=6,epsilon=.025,proper_endpoint=None,threads=8,budget_seconds=600,
                  rtol=3e-8,atol=3e-10,max_step=.05,progress=False,gauge='lapse_controlled',state_sink=None,backend='hopf',quadrature_order=None,measure_work=True,dense_sink=None,kick_amplitude=0.,kick_source=1):
    set_num_threads(threads)
    from scipy.fft import set_workers
    with threadpool_limits(limits=threads),set_workers(threads):
        system=ParentGH(degree,epsilon,gauge=gauge,backend=backend,quadrature_order=quadrature_order,
            kick_amplitude=kick_amplitude,kick_source=kick_source);y0=system.initial()
        endpoint=4.364851930251994 if proper_endpoint is None else proper_endpoint
        if endpoint<=0:raise ValueError('positive proper endpoint required')
        system.started=time.perf_counter();system.budget=budget_seconds
        snapshots=[system.diagnostic(0.,y0)];targets=list(np.linspace(0,endpoint,5)[1:])
        initial_energy=snapshots[0]['excitation_normal_energy_over_initial_round_vacuum_energy']
        if measure_work:
            system.passive_work=PassiveDOP853Quadrature()
            snapshots[0]['integrated_excitation_work_normalized']=0.
            snapshots[0]['excitation_work_balance_residual']=0.
        if state_sink is not None:state_sink(0.,y0)
        solver=DOP853(system.rhs,0.,y0,t_bound=8/system.r,rtol=rtol,atol=atol,max_step=max_step)
        if measure_work:system.passive_work.initialize(solver)
        steps=0;norm_drift=0.
        failure=None
        angle_prev=system.observer_angle(y0);phase=0.;clock_valid=True
        try:
            while solver.status=='running' and targets:
                message=solver.step();steps+=1
                if solver.status=='failed':
                    failure=message or 'GH time integrator failed';break
                if measure_work:system.passive_work.accepted_step(solver)
                x=solver.y[system.dof:system.dof+4];norm_drift=max(norm_drift,abs(np.linalg.norm(x)-1))
                angle=system.observer_angle(solver.y)
                before_phase=phase;before_angle=angle_prev
                if angle is not None and angle_prev is not None:
                    increment=math.atan2(math.sin(angle-angle_prev),math.cos(angle-angle_prev))
                    if abs(increment)>.9*math.pi:clock_valid=False
                    phase+=increment
                elif epsilon:clock_valid=False
                angle_prev=angle
                dense=None
                if dense_sink is not None:
                    dense=solver.dense_output()
                    dense_sink(solver.t_old,solver.t,dense)
                if solver.y[system.dof+8]>=targets[0]:
                    if dense is None:dense=solver.dense_output()
                    while targets and solver.y[system.dof+8]>=targets[0]:
                        target=targets.pop(0)
                        t=brentq(lambda v:dense(v)[system.dof+8]-target,solver.t_old,solver.t,xtol=2e-12)
                        state=dense(t);a=system.observer_angle(state);partial=0.
                        if a is not None and before_angle is not None:
                            partial=math.atan2(math.sin(a-before_angle),math.cos(a-before_angle))
                        state[-1]=before_phase+partial
                        point=system.diagnostic(t,state);point['scalar_clock_sampling_check_passed']=clock_valid
                        if measure_work:
                            integral=system.passive_work.dense_value(solver,t)
                            point['integrated_excitation_work_normalized']=integral
                            point['excitation_work_balance_residual']=point['excitation_normal_energy_over_initial_round_vacuum_energy']-initial_energy-integral
                        if not clock_valid:point['observer_sampled_unwrapped_phase']=None
                        snapshots.append(point)
                        if state_sink is not None:state_sink(t,state)
                        if progress:print(f'GH L={degree} epsilon={epsilon:g}: proper {target:.6g}/{endpoint:.6g}, H-residual={snapshots[-1]["ADM_H_RMS_over_H_squared"]:.3g}',flush=True)
        except (TimeoutError,FloatingPointError) as e:
            failure=str(e)
        if targets:
            if failure is None:failure='proper-time endpoint not reached inside coordinate domain'
            if solver.t>snapshots[-1]['coordinate_t']:
                last=solver.y.copy();last[-1]=phase
                point=system.diagnostic(solver.t,last);point['scalar_clock_sampling_check_passed']=clock_valid
                if measure_work:
                    integral=system.passive_work.value
                    point['integrated_excitation_work_normalized']=integral
                    point['excitation_work_balance_residual']=point['excitation_normal_energy_over_initial_round_vacuum_energy']-initial_energy-integral
                if not clock_valid:point['observer_sampled_unwrapped_phase']=None
                snapshots.append(point)
                if state_sink is not None:state_sink(solver.t,last)
        return dict(schema='nsc-atlas-parent-gh-v3' if kick_amplitude else 'nsc-atlas-parent-gh-v2',degree=degree,epsilon=epsilon,
            spatial_transform=backend,spatial_quadrature_points=system.N,
            status='incomplete' if targets else 'endpoint_reached',stop_reason=failure,
            requested_observer_proper_endpoint=endpoint,invalid_RK_stages=system.invalid_stages,
            source='Inherited two-mode field plus constrained canonical momentum kick' if kick_amplitude else 'Unchanged inherited nonspherical two-mode field, finite nonlinear constrained initial geometry',
            initial_kick=system.initial_kick_metadata,
            formulation='Second-order covariant GH, static round reference, canonical volume and real scalar',
            prescribed_GH_gauge=gauge,
            lapse_feedback_coefficient=system.r if gauge=='lapse_controlled' else 0.,
            gauge_scope='Metric-dependent coordinate source; no added physical stress or physical-mode suppression.',
            model_inputs=dict(vacuum=system.parent.vacuum,H_reference=system.r,
                scalar_mass_squared=system.parent.mass_squared,volume_kinetic=system.K,sstar=system.sstar,
                reference_proper_endpoint_origin='Previously specified 300-phase-pass reference event',
                **(dict(kick_amplitude=float(kick_amplitude),kick_source=int(kick_source)) if kick_amplitude else {})),
            gauge_constraint_damping=system.r,physical_mode_suppression=False,
            symmetry='Exact proper paired-flip characters; no spherical average',
            scalar_phase_scope='Along the initially normal geodesic, at fixed proper duration; no finite-amplitude phase-event shift is inferred by linearization.',
            snapshots=snapshots,steps=steps,rhs_calls=system.calls,
            observable_quaternion_max_norm_drift=norm_drift,
            elapsed_seconds=time.perf_counter()-system.started,threads=threads,
            temporal_tolerances=dict(rtol=rtol,atol=atol,max_step=max_step),
            passive_matter_work=dict(enabled=measure_work,valid=system.work_valid,
                scope='Vacuum-subtracted normal excitation energy and geometric work; passive DOP853 stage quadrature never controls the PDE step size.'),
            SI_scale_selected=False,child_formation_proved=False)
