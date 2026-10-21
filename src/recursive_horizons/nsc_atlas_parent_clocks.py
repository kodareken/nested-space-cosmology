"""Operational radial clock transfer for the existing protected parent packet.

Receiver proper age and emitter initial proper separation define the events.
The two amplitudes belong to the existing density/phase preparation family;
they are inferred inputs, not a newly selected primordial spectrum.
"""
from dataclasses import dataclass,field
from functools import cached_property
import math

import numpy as np
from scipy.integrate import quad_vec,quad,solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_parent_perturbations import ParentPacket


@dataclass
class ParentClocks:
    packet: ParentPacket = field(default_factory=ParentPacket)
    cutoff: int = 128
    observation_efold: float = .8
    rtol: float = 2e-10
    atol: float = 2e-12

    def __post_init__(self):
        if not 0<self.observation_efold<=.9:
            raise ValueError('receiver must leave a finite margin inside the one-e-fold solution')
        self.evolution,_,self.weights=self.packet.evolve(self.cutoff,rtol=self.rtol,atol=self.atol)
        self.initial=self.packet.model.initial
        self.modes=self.evolution.modes
        self.to=brentq(lambda t:self.evolution.solution.sol(t)[2]-self.observation_efold,
            0,self.evolution.solution.t[-1],xtol=1e-13)

    def background(self,t):
        if not 0<=t<=self.evolution.solution.t[-1]+1e-12:
            raise ValueError('clock evaluation outside evolved parent')
        chi,Pi,N,H,_,eta,_=self.evolution.solution.sol(t)[:7]
        return dict(chi=chi,Pi=Pi,N=N,H=H,eta=eta,a=self.initial['a']*math.exp(N))

    def harmonics(self,theta):
        m=self.modes+1.;lam=self.modes*(self.modes+2)
        if abs(theta)<1e-5:
            Q=1-lam*theta**2/6+lam*(3*lam-4)*theta**4/360
            Qr=-lam*theta/3+lam*(3*lam-4)*theta**3/90
        else:
            sint=math.sin(theta);cost=math.cos(theta)
            Q=np.sin(m*theta)/(m*sint)
            Qr=np.cos(m*theta)/sint-np.sin(m*theta)*cost/(m*sint*sint)
        return Q,Qr

    def potential(self,t,theta):
        y=self.evolution.solution.sol(t)
        chi,Pi,_,H=y[:4]
        X,_,Phi=y[7:].reshape(3,len(self.modes),2)
        Phid=-H*Phi+self.initial['kappa']*Pi*X/2
        Q,Qr=self.harmonics(theta)
        return dict(value=(self.weights*Q)@Phi,
            radial=(self.weights*Qr)@Phi,time=(self.weights*Q)@Phid)

    @cached_property
    def observer_time_shift(self):
        return -quad_vec(lambda t:self.potential(t,0.)['value'],0,self.to,epsabs=2e-12,epsrel=2e-11)[0]

    def response(self,separation):
        bo=self.background(self.to);s=float(separation)
        if not 0<s<bo['eta']:
            raise ValueError('emitter light must be emitted after the retained initial slice')
        te=brentq(lambda t:bo['eta']-self.background(t)['eta']-s,0,self.to,xtol=1e-13)
        be=self.background(te)
        displacement0=quad_vec(lambda theta:self.potential(0,theta)['value'],0,s,epsabs=1e-13)[0]
        def source_rhs(t,y):
            b=self.background(t);v=y[:2]
            return np.r_[-b['H']*v-self.potential(t,s)['radial']/b['a'],v/b['a']]
        source=solve_ivp(source_rhs,(0,self.to),np.r_[0.,0.,displacement0],method='DOP853',
            rtol=2e-11,atol=2e-13,dense_output=True,max_step=.025)
        if not source.success:raise RuntimeError(source.message)
        ve,dx=source.sol(te).reshape(2,2)
        def integrands(t):
            b=self.background(t);theta=bo['eta']-b['eta'];pot=self.potential(t,theta)
            return np.r_[pot['value']/b['a'],pot['time']]
        integrals=quad_vec(integrands,te,self.to,epsabs=2e-12,epsrel=2e-11)[0]
        dto=self.observer_time_shift
        dte=be['a']*(dto/bo['a']+2*integrals[:2]-dx)
        endpoint=self.potential(te,s)['value']-self.potential(self.to,0.)['value']
        parts=dict(emission_event=be['H']*dte,receiver_event=-bo['H']*dto,
            endpoint_potential=endpoint,emitter_Doppler=-ve,integrated_time_potential=2*integrals[2:])
        total=sum(parts.values())
        return dict(initial_angle=s,initial_proper_separation=self.initial['a']*s,
            receiver_proper_age=self.to,receiver_background_efold=self.observation_efold,
            background_emission_time=te,background_emission_efold=float(be['N']),
            background_frequency_ratio=math.exp(be['N']-bo['N']),
            initial_coordinate_displacement=displacement0.tolist(),
            emitter_displacement=dx.tolist(),emitter_velocity=ve.tolist(),
            receiver_time_shift=dto.tolist(),emission_time_shift=dte.tolist(),
            log_frequency_response=total.tolist(),parts={k:v.tolist() for k,v in parts.items()})

    def finite_metric_clock(self,separation,amplitudes):
        """Independent nonlinear ray check on a prescribed perturbed metric.

        The metric is a finite deformation of a linear solution, not an exact
        nonlinear Einstein-scalar solution. Its first variation validates
        the observer formula, including proper-distance/age event selection.
        """
        amplitudes=np.asarray(amplitudes,dtype=float)
        if amplitudes.shape!=(2,) or max(abs(amplitudes))>.002:
            raise ValueError('two finite metric variation amplitudes within .002 required')
        s=float(separation)
        if not 0<s<self.background(self.to)['eta']:
            raise ValueError('source outside post-initial light cone')
        def metric(t,theta):
            b=self.background(t);pot=self.potential(t,theta)
            F=float(pot['value']@amplitudes);Ft=float(pot['time']@amplitudes);Fr=float(pot['radial']@amplitudes)
            A=1+2*F;B=b['a']**2*(1-2*F)
            if A<=0 or B<=0:raise RuntimeError('metric variation lost Lorentzian signature')
            return b,A,B,F,Ft,Fr
        def proper_age(t):
            return quad(lambda x:math.sqrt(metric(x,0.)[1]),0,t,epsabs=2e-13,epsrel=2e-13)[0]
        to=brentq(lambda t:proper_age(t)-self.to,.98*self.to,1.02*self.to,xtol=2e-13)
        def distance(theta):
            return quad(lambda x:math.sqrt(metric(0,x)[2]),0,theta,epsabs=2e-13,epsrel=2e-13)[0]
        theta0=brentq(lambda theta:distance(theta)-self.initial['a']*s,.98*s,1.02*s,xtol=2e-13)
        def source_rhs(t,y):
            theta,v=y;b,A,B,F,Ft,Fr=metric(t,theta)
            return (math.sqrt(A/B)*v,-(1-v*v)*(Fr/math.sqrt(A*B)+v*(b['H']-Ft/(1-2*F))))
        source=solve_ivp(source_rhs,(0,to),(theta0,0.),method='DOP853',dense_output=True,
            rtol=2e-12,atol=2e-14,max_step=.02)
        if not source.success:raise RuntimeError(source.message)
        def null_rhs(t,y):
            b,A,B,F,Ft,_=metric(t,y[0])
            # log(-p_t), derived from the exact null Hamilton equations.
            return (-math.sqrt(A/B),-b['H']+2*Ft/(1-4*F*F))
        def emission(t,y):return y[0]-source.sol(t)[0]
        emission.terminal=True
        photon=solve_ivp(null_rhs,(to,0),(0.,0.),method='DOP853',events=emission,
            rtol=2e-12,atol=2e-14,max_step=.02)
        if not photon.success or not len(photon.t_events[0]):raise RuntimeError('exact radial emission event not found')
        te=float(photon.t_events[0][0]);theta,logpt=photon.y_events[0][0]
        velocity=float(source.sol(te)[1]);Ao=metric(to,0.)[1];Ae=metric(te,theta)[1]
        value=-float(logpt)+.5*math.log(Ae/Ao)-math.atanh(velocity)
        return dict(log_frequency_ratio=value,receiver_coordinate_time=to,
            emission_coordinate_time=te,source_initial_angle=theta0,
            source_emission_angle=float(theta),source_velocity=velocity,
            source_rhs_calls=source.nfev,null_rhs_calls=photon.nfev)


def infer_two_amplitudes(rows,measurements):
    matrix=np.asarray([row['log_frequency_response'] for row in rows[:2]])
    if matrix.shape!=(2,2) or np.linalg.cond(matrix)>1e5:
        raise ValueError('declared two-clock calibration is not resolved')
    return np.linalg.solve(matrix,np.asarray(measurements,dtype=float))
