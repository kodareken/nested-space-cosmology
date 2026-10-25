"""Finite fixed-clock optics from the actual formed scalar envelope.

The weak averaged metric is explicit. Field currents provide time
derivatives; no differentiated saved waveform or new lens profile is used.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.fft import dst,dct
from scipy.integrate import solve_ivp,quad
from scipy.interpolate import CubicSpline,PPoly
from scipy.optimize import brentq

from .nsc_atlas_formation_window import FormationBackground
from .nsc_atlas_seed_collapse import boundary_potential
from .nsc_atlas_hierarchy_evolution import kinetic_step


def weighted_primitive(spline,power):
    """Exact PPoly antiderivative of r^power times a cubic spline."""
    asc=spline.c[::-1];left=spline.x[:-1];coeff=np.zeros((4+power,len(left)))
    for j in range(power+1):
        coeff[j:j+4]+=asc*(math.comb(power,j)*left**(power-j))
    return PPoly(coeff[::-1],spline.x).antiderivative()


def spectral_derivative(u,outer):
    N=len(u)+1;k=np.pi*np.arange(1,N)/outer
    coefficients=np.r_[0.,k*dst(u,type=1,norm='ortho'),0.]
    return dct(coefficients,type=1)/math.sqrt(2*N)


class PotentialSlice:
    def __init__(self,r,u,A,length,background):
        dr=r[1]-r[0];outer=r[-1]+dr
        derivative=spectral_derivative(u,outer)
        f=1+u/r;fr=(derivative[1:-1]-u/r)/r
        center=1+derivative[0]
        x=np.r_[0.,r,outer]
        density=np.r_[abs(center)**2-1,abs(f)**2-1,0.]
        current=np.r_[0.,np.imag(np.conj(f)*fr),np.imag(derivative[-1])/outer]
        rho=CubicSpline(x,density,bc_type=((1,0.),'not-a-knot'))
        j=CubicSpline(x,current)
        M=weighted_primitive(rho,2);T=weighted_primitive(rho,1);J=j.antiderivative()
        self.pref=background.excitation*length*length/A
        self.a=1/(background.frequency*length*length*A*A*float(background.H(1/A)))
        self.M0=float(M(0));self.T0=float(T(0));self.Jend=float(J(outer))
        self.Mend=float(M(outer)-M(0));self.Tend=float(T(outer)-T(0))
        self.outer=outer;self.A=A
        # Only retain the local query region's polynomials. Constants retain
        # the full outer density/current integrals and the same boundary clock.
        count=int(np.searchsorted(x,.35,side='right'))
        def trim(p):return PPoly(p.c[:,:count].copy(),p.x[:count+1].copy(),extrapolate=False)
        self.mass=trim(M);self.tail=trim(T);self.current=trim(j);self.current_primitive=trim(J)

    def evaluate(self,r):
        r=float(r)
        if not 0<=r<=self.mass.x[-1]:raise ValueError('ray left the stored radial query domain')
        M=float(self.mass(r)-self.M0);T=float(self.tail(r)-self.T0)
        V=-.5*((M/r if r else 0.)+self.Tend-T-self.Mend/self.outer)
        phi=self.pref*V
        phir=self.pref*.5*M/(r*r) if r else 0.
        phiN=-phi+self.pref*.5*self.a*(self.Jend-float(self.current_primitive(r)))
        phiNr=-phir-self.pref*.5*self.a*float(self.current(r))
        return np.array([phi,phir,phiN,phiNr])


def advance_source(r,u,N0,N1,length,background,max_step):
    """Same symmetric field flows; reversible for a signed log-scale step."""
    steps=max(1,math.ceil(abs(N1-N0)/max_step));h=(N1-N0)/steps
    dr=r[1]-r[0];outer=r[-1]+dr;k=np.pi*np.arange(1,len(r)+1)/outer
    u=u.copy()
    for i in range(steps):
        N=N0+i*h;A=np.exp(N+h/2+np.array([-1.,1.])*h/(2*math.sqrt(3)));H=background.H(1/A)
        Q=h*np.mean(1/(background.frequency*length**2*A*A*H))
        P=h*np.mean(background.frequency*background.excitation*length**2/(A*H))
        V,_=boundary_potential(u,r,dr,outer)
        u=(u+r)*np.exp(-.5j*P*V)-r;u=kinetic_step(u,k,Q)
        V,_=boundary_potential(u,r,dr,outer)
        u=(u+r)*np.exp(-.5j*P*V)-r
    return u


class SourceSlab:
    def __init__(self,r,u,*,length,endpoint_scale=1.2,spacing=1e-5,field_step=2e-6):
        self.bg=FormationBackground.inherited();self.length=length
        self.center=math.log(endpoint_scale);self.start=self.center-.0002;self.end=self.center+.0003
        # Include the authenticated endpoint exactly; integrate both ways from it.
        negative=np.linspace(self.center,self.start,math.ceil((self.center-self.start)/spacing)+1)
        positive=np.linspace(self.center,self.end,math.ceil((self.end-self.center)/spacing)+1)
        records={self.center:PotentialSlice(r,u,endpoint_scale,length,self.bg)}
        for times in (negative,positive):
            state=u.copy()
            for a,b in zip(times[:-1],times[1:]):
                state=advance_source(r,state,a,b,length,self.bg,field_step)
                records[float(b)]=PotentialSlice(r,state,math.exp(b),length,self.bg)
        self.times=np.array(sorted(records));self.slices=[records[float(n)] for n in self.times]
        back=advance_source(r,u,self.center,self.start,length,self.bg,field_step)
        echo=advance_source(r,back,self.start,self.center,length,self.bg,field_step)
        self.echo_relative_error=float(np.linalg.norm(echo-u)/np.linalg.norm(u))
        self.c=self.bg.curvature*length*length

    def evaluate(self,N,r,*,homogeneous=False,frozen=False):
        if homogeneous:return np.zeros(3)
        if frozen:
            value=self.slices[int(np.argmin(abs(self.times-self.center)))].evaluate(r)
            return np.array([value[0],value[1],0.])
        if N<self.start-1e-12 or N>self.end+1e-12:raise ValueError('ray left short source interval')
        j=min(len(self.times)-2,max(0,int(np.searchsorted(self.times,N))-1))
        left,right=self.times[j:j+2];h=right-left;x=(N-left)/h
        l=self.slices[j].evaluate(r);rr=self.slices[j+1].evaluate(r)
        H00=2*x**3-3*x*x+1;H10=x**3-2*x*x+x;H01=-2*x**3+3*x*x;H11=x**3-x*x
        phi=H00*l[0]+h*H10*l[2]+H01*rr[0]+h*H11*rr[2]
        phir=H00*l[1]+h*H10*l[3]+H01*rr[1]+h*H11*rr[3]
        phiN=(6*x*x-6*x)/h*l[0]+(3*x*x-4*x+1)*l[2]+(-6*x*x+6*x)/h*rr[0]+(3*x*x-2*x)*rr[2]
        return np.array([phi,phir,phiN])


@dataclass(frozen=True)
class Observer:
    radius: float=.25
    separation: float=math.pi-.08
    ray_rtol: float=2e-12
    ray_atol: float=2e-14


def clock_age(slab,N,observer=Observer(),*,homogeneous=False,frozen=False):
    def integrand(n):
        phi=slab.evaluate(n,observer.radius,homogeneous=homogeneous,frozen=frozen)[0]
        return math.sqrt(1+2*phi)/float(slab.bg.H(math.exp(-n)))
    return quad(integrand,slab.start,N,epsabs=2e-17,epsrel=2e-12)[0]


def emission_event(slab,age,observer,*,homogeneous=False,frozen=False):
    return brentq(lambda n:clock_age(slab,n,observer,homogeneous=homogeneous,frozen=frozen)-age,
                  slab.start,slab.end,xtol=5e-17,rtol=1e-15)


def ray(slab,Ne,angle,observer=Observer(),*,homogeneous=False,frozen=False):
    rs=observer.radius;B0=1-slab.c*rs*rs
    j=rs*math.sin(angle);p0=-math.cos(angle)/math.sqrt(B0)
    phie=slab.evaluate(Ne,rs,homogeneous=homogeneous,frozen=frozen)[0]
    def rhs(q,y):
        r,az,p,N,I,age=y
        phi,phir,phiN=slab.evaluate(N,r,homogeneous=homogeneous,frozen=frozen)
        B=1-slab.c*r*r;Br=-2*slab.c*r
        E=math.sqrt(B*p*p+j*j/(r*r));F=math.sqrt((1+2*phi)/(1-2*phi))
        Fr=F*2*phir/(1-4*phi*phi)
        Nq=math.exp(N)*slab.length*float(slab.bg.H(math.exp(-N)))
        screenphi=slab.evaluate(N,rs,homogeneous=homogeneous,frozen=frozen)[0]
        return [F*B*p/E,F*j/(E*r*r),-Fr*E-F*(Br*p*p-2*j*j/r**3)/(2*E),
                Nq,2*phiN*Nq/(1-4*phi*phi),math.sqrt(1+2*screenphi)*math.exp(N)*slab.length]
    def arrival(q,y):return y[0]-rs
    arrival.direction=1;arrival.terminal=True
    sol=solve_ivp(rhs,(0,2*rs*1.02),(rs,0.,p0,Ne,0.,0.),method='DOP853',
        rtol=observer.ray_rtol,atol=observer.ray_atol,events=arrival,max_step=.005)
    if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('receiver sphere not reached')
    r,az,p,N,I,flight=sol.y_events[0][0];phio=slab.evaluate(N,rs,homogeneous=homogeneous,frozen=frozen)[0]
    E=math.sqrt(B0*p*p+j*j/(rs*rs))
    # Cancellation-safe frequency contrast removes the homogeneous expansion.
    Z=.5*(math.log1p(2*phie)-math.log1p(2*phio))+I
    endpoint_Z=math.log(E)+.5*(math.log1p(-2*phie)-math.log1p(-2*phio))
    angle_o=math.atan2(j/rs,math.sqrt(B0)*p)
    return dict(emission_N=Ne,arrival_N=float(N),arrival_angle=float(az),launch_angle=angle,
        observed_direction=float(angle_o),proper_flight_time=float(flight),
        log_frequency=float(Ne-N+Z),expansion_removed_log_frequency=float(Z),
        integrated_time_potential=float(I),endpoint_frequency_identity_gap=float(endpoint_Z-Z),
        null_conformal_length=float(sol.t_events[0][0]),ray_function_evaluations=sol.nfev)


def fixed_observer(slab,age,observer=Observer(),*,homogeneous=False,frozen=False):
    Ne=emission_event(slab,age,observer,homogeneous=homogeneous,frozen=frozen)
    rs=observer.radius;theta=observer.separation
    reference=math.atan2(math.cos(theta/2),math.sqrt(1-slab.c*rs*rs)*math.sin(theta/2))
    calls={}
    def endpoint(angle):
        value=ray(slab,Ne,float(angle),observer,homogeneous=homogeneous,frozen=frozen)
        calls[float(angle)]=value
        return value['arrival_angle']-theta
    if homogeneous:
        angle=reference;endpoint(angle)
    else:angle=brentq(endpoint,reference-.0001,reference+.0001,xtol=2e-14)
    result=calls.get(float(angle)) or ray(slab,Ne,float(angle),observer,homogeneous=homogeneous,frozen=frozen)
    result['fixed_receiver_residual']=result['arrival_angle']-theta
    result['source_proper_age']=age
    return result


def comparison(slab,age=None,observer=Observer(),*,frozen=False):
    if age is None:age=clock_age(slab,slab.center,observer,homogeneous=True)
    base=fixed_observer(slab,age,observer,homogeneous=True)
    signal=fixed_observer(slab,age,observer,frozen=frozen)
    delay=signal['proper_flight_time']-base['proper_flight_time']
    dlog=signal['log_frequency']-base['log_frequency']
    phio=slab.evaluate(signal['arrival_N'],observer.radius,frozen=frozen)[0]
    beta=signal['observed_direction'];beta0=base['observed_direction']
    log_radius_change=signal['arrival_N']-base['arrival_N']+.5*math.log1p(-2*phio)
    angular=math.exp(base['arrival_N'])*slab.length*observer.radius*(
        math.expm1(log_radius_change)*math.sin(beta)+2*math.cos((beta+beta0)/2)*math.sin((beta-beta0)/2))
    return dict(signal=signal,homogeneous=base,proper_delay_excess=delay,
        image_direction_change=signal['observed_direction']-base['observed_direction'],
        log_frequency_change=dlog,
        endpoint_delay_derivative_prediction=math.exp(-base['log_frequency'])*math.expm1(-dlog),
        endpoint_angular_derivative_prediction=angular,
        expansion_removed_log_frequency=signal['expansion_removed_log_frequency'])
