"""Timelike calibration and a data-only finite-endpoint lens prediction.

Prediction consumes measured velocity changes and kinematics, not the
source potential or the withheld null result. All data here are synthetic.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

from .nsc_atlas_formed_observer import clock_age

TARGET=math.pi-.12


class RadialFall:
    def __init__(self,slab,radius=.25,*,start=None,end=None,velocity=0.,homogeneous=False):
        self.slab=slab;self.ell=math.exp(slab.center)*slab.length
        self.homogeneous=homogeneous
        self.start=slab.start if start is None else start
        self.end=slab.end if end is None else end
        phi=slab.evaluate(self.start,radius,homogeneous=homogeneous)[0]
        B=1-slab.c*radius**2;s=math.exp(self.start)*slab.length*math.sqrt(1-2*phi)
        p0=s/self.ell*velocity/math.sqrt(B*(1-velocity**2))
        self.solution=solve_ivp(self.rhs,(self.start,self.end),(radius,p0,0.),method='DOP853',
            rtol=2e-12,atol=(2e-14,2e-16,2e-16),dense_output=True,max_step=5e-5)
        if not self.solution.success:raise RuntimeError(self.solution.message)

    def rhs(self,N,y):
        r,p,tau=y;phi,phir,_=self.slab.evaluate(N,r,homogeneous=self.homogeneous)
        H=float(self.slab.bg.H(math.exp(-N)));alpha=math.sqrt(1+2*phi)
        s2=(math.exp(N)*self.slab.length)**2*(1-2*phi)
        B=1-self.slab.c*r*r;K=self.ell**2/s2;G=B*K
        Gr=K*(-2*self.slab.c*r+2*B*phir/(1-2*phi))
        E=math.sqrt(1+G*p*p)
        return [alpha*G*p/(H*self.ell*E),
            -(phir*E/alpha+alpha*Gr*p*p/(2*E))/(H*self.ell),alpha/(H*E)]

    def at(self,N):
        r,p,tau=self.solution.sol(N)
        phi=self.slab.evaluate(N,r,homogeneous=self.homogeneous)[0]
        s=math.exp(N)*self.slab.length*math.sqrt(1-2*phi)
        w=math.sqrt(1-self.slab.c*r*r)*self.ell*p/s
        return dict(radius=float(r),momentum=float(p),proper_age=float(tau),velocity=float(w/math.sqrt(1+w*w)))

    def event(self,age):
        return brentq(lambda N:self.at(N)['proper_age']-age,self.start,self.end,xtol=5e-17,rtol=1e-15)


def curvature_distance(re,ro,theta,c):
    Be=math.sqrt(1-c*re*re);Bo=math.sqrt(1-c*ro*ro)
    chord2=(re-ro)**2+2*re*ro*(1-math.cos(theta))+c*(re*re-ro*ro)**2/(Be+Bo)**2
    if c==0:return math.sqrt(chord2)
    return 2*math.asin(math.sqrt(c*chord2)/2)/math.sqrt(c)


def geometric_angle(re,ro,theta,c):
    return math.atan2(re*math.sin(theta),math.sqrt(1-c*re*re)*ro-math.sqrt(1-c*ro*ro)*re*math.cos(theta))


def aberrate(angle,velocity):
    return math.atan2(math.sin(angle),
        (math.cos(angle)-velocity)/math.sqrt(1-velocity*velocity))


def impulse(slab,N,r,half_width,velocity=0.):
    plus=RadialFall(slab,r,start=N,end=N+half_width,velocity=velocity).at(N+half_width)
    minus=RadialFall(slab,r,start=N,end=N-half_width,velocity=velocity).at(N-half_width)
    duration=plus['proper_age']-minus['proper_age']
    accel=-(plus['velocity']-minus['velocity'])/duration
    return dict(N=N,radius=r,initial_velocity=velocity,
        velocity_before=minus['velocity'],velocity_after=plus['velocity'],proper_interval=duration,
        inferred_relative_acceleration=accel)


def precheck(slab,radii=(.0073,.015,.05,.25)):
    snapshot=slab.slices[int(np.argmin(abs(slab.times-slab.center)))];A=math.exp(slab.center)
    H=float(slab.bg.H(1/A));out=[]
    for r in radii:
        phi,force,_,forceN=snapshot.evaluate(r)
        M=float(snapshot.mass(r)-snapshot.M0);Mr=float(snapshot.mass.derivative()(r))
        force_r=snapshot.pref*.5*(Mr/r**2-2*M/r**3)
        log_rate=H*(forceN-force-r*force_r)/force
        omega=math.sqrt(abs(force)/(A*A*slab.length**2*r))
        out.append(dict(radius=r,force_log_rate_at_fixed_physical_radius=log_rate,
            dynamical_frequency=omega,force_change_per_dynamical_time=abs(log_rate)/omega,
            force_change_per_nominal_orbit=2*math.pi*abs(log_rate)/omega))
    return out


def calibrate(slab,*,theta=TARGET,nodes_per_panel=24,half_width=5e-5,release_error=0.):
    # Source/receiver motion is calibration information. No photon equation
    # or target optical result is called here.
    observers=RadialFall(slab)
    age=clock_age(slab,slab.center,homogeneous=True)
    Ne=observers.event(age);source=observers.at(Ne);re=source['radius']
    def conformal_to(No):
        return quad(lambda n:1/(math.exp(n)*slab.length*float(slab.bg.H(math.exp(-n)))),Ne,No,epsabs=2e-13)[0]
    No=brentq(lambda n:conformal_to(n)-curvature_distance(re,observers.at(n)['radius'],theta,slab.c),
              Ne,slab.end,xtol=5e-17,rtol=1e-15)
    receiver=observers.at(No);ro=receiver['radius']
    d=math.sqrt((re-ro)**2+2*re*ro*(1-math.cos(theta)))
    impact=re*ro*math.sin(theta)/d
    ze=-math.sqrt(re*re-impact*impact);zo=math.sqrt(ro*ro-impact*impact)
    n_of_q=solve_ivp(lambda q,y:[math.exp(y[0])*slab.length*float(slab.bg.H(math.exp(-y[0])))],
        (0,zo-ze),(Ne,),method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True)
    z0,w0=roots_legendre(nodes_per_panel);data=[]
    for left,right in ((ze,0.),(0.,zo)):
        for z,weight in zip(left+(z0+1)*(right-left)/2,w0*(right-left)/2):
            r=math.sqrt(impact*impact+z*z);N=float(n_of_q.sol(z-ze)[0])
            point=impulse(slab,N,r,half_width)
            factor=math.exp(N)*slab.length/math.sqrt(1-slab.c*r*r)
            point.update(chord_coordinate=float(z),quadrature_weight=float(weight),
                force_conversion_factor=factor,born_weight=2*(z-ze)/d*impact/r)
            if release_error:
                alternatives=[impulse(slab,N,r,half_width,v)['inferred_relative_acceleration'] for v in (-release_error,release_error)]
                point['release_acceleration_envelope']=max(abs(a-point['inferred_relative_acceleration']) for a in alternatives)
            else:point['release_acceleration_envelope']=0.
            data.append(point)
    return dict(receiver_azimuth=theta,source_proper_age=age,emission_N=Ne,
        reference_arrival_N=No,source=source,receiver=receiver,curvature=slab.c,length=slab.length,
        background_angle=geometric_angle(.25,.25,theta,slab.c),
        nodes_per_panel=nodes_per_panel,half_width=half_width,release_error=release_error,probes=data)


def predict_from_motion(data):
    """Only calibration measurements and geometry enter this prediction."""
    lens=sum(p['quadrature_weight']*p['born_weight']*p['force_conversion_factor']*
             p['inferred_relative_acceleration'] for p in data['probes'])
    state_error=sum(abs(p['quadrature_weight']*p['born_weight']*p['force_conversion_factor'])*
                    p['release_acceleration_envelope'] for p in data['probes'])
    theta=data['receiver_azimuth'];re=data['source']['radius'];ro=data['receiver']['radius']
    geometric=geometric_angle(re,ro,theta,data['curvature']);v=data['receiver']['velocity']
    image=aberrate(geometric+lens,v)
    return dict(predicted_image_direction=image,predicted_image_shift=image-data['background_angle'],
        normal_frame_lens_contribution=lens,geometric_direction=geometric,
        kinematic_direction_change=aberrate(geometric,v)-data['background_angle'],
        release_state_envelope=state_error)


def moving_ray(slab,worldline,Ne,angle,*,theta=TARGET,rtol=2e-12,atol=2e-14):
    """Withheld optical calculation: full null boundary problem, no Born data."""
    emitter=worldline.at(Ne);re=emitter['radius'];j=re*math.sin(angle)
    p0=-math.cos(angle)/math.sqrt(1-slab.c*re*re)
    phie=slab.evaluate(Ne,re)[0]
    def rhs(q,y):
        r,az,p,N,I=y;phi,phir,phiN=slab.evaluate(N,r)
        B=1-slab.c*r*r;Br=-2*slab.c*r;E=math.sqrt(B*p*p+j*j/r**2)
        F=math.sqrt((1+2*phi)/(1-2*phi));Fr=F*2*phir/(1-4*phi*phi)
        Nq=math.exp(N)*slab.length*float(slab.bg.H(math.exp(-N)))
        return [F*B*p/E,F*j/(E*r*r),-Fr*E-F*(Br*p*p-2*j*j/r**3)/(2*E),Nq,
                2*phiN*Nq/(1-4*phi*phi)]
    def meet(q,y):return y[0]-worldline.at(y[3])['radius']
    meet.direction=1;meet.terminal=True
    sol=solve_ivp(rhs,(0,.51),(re,0.,p0,Ne,0.),method='DOP853',rtol=rtol,atol=atol,
        events=meet,max_step=.005)
    if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('moving receiver not intercepted')
    r,az,p,N,I=sol.y_events[0][0];observer=worldline.at(N)
    phio=slab.evaluate(N,r)[0];B=1-slab.c*r*r
    normal_angle=math.atan2(j/r,math.sqrt(B)*p)
    v_e=emitter['velocity'];v_o=observer['velocity'];mu_e=-math.cos(angle);mu_o=math.cos(normal_angle)
    Z=.5*(math.log1p(2*phie)-math.log1p(2*phio))+I
    doppler=-.5*math.log1p(-v_o*v_o)+math.log1p(-v_o*mu_o)+.5*math.log1p(-v_e*v_e)-math.log1p(-v_e*mu_e)
    return dict(arrival_angle=float(az),arrival_N=float(N),normal_observed_direction=normal_angle,
        observed_direction=aberrate(normal_angle,v_o),launch_angle=angle,
        source_velocity=v_e,receiver_velocity=v_o,
        proper_flight_time=observer['proper_age']-emitter['proper_age'],
        log_frequency=float(Ne-N+Z+doppler),log_endpoint_doppler=doppler,
        expansion_removed_normal_frequency=float(Z),function_evaluations=sol.nfev)


def evaluate_withheld(slab,*,theta=TARGET,rtol=2e-12,atol=2e-14,age=None):
    fall=RadialFall(slab)
    if age is None:age=clock_age(slab,slab.center,homogeneous=True)
    Ne=fall.event(age)
    guess=(math.pi-theta)/2;calls={}
    def residual(angle):
        row=moving_ray(slab,fall,Ne,float(angle),theta=theta,rtol=rtol,atol=atol)
        calls[float(angle)]=row;return row['arrival_angle']-theta
    angle=brentq(residual,guess-.0001,guess+.0001,xtol=2e-14)
    row=calls.get(float(angle)) or moving_ray(slab,fall,Ne,angle,theta=theta,rtol=rtol,atol=atol)
    row['image_shift']=row['observed_direction']-geometric_angle(.25,.25,theta,slab.c)
    row['receiver_residual']=row['arrival_angle']-theta
    row['source_proper_age']=age
    return row
