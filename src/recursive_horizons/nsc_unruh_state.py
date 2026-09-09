"""Phase-resolved fermionic horizon covariance for the parent-state interface.

The affine-horizon prescription is imported. Its application to the full
NSC domain and the finite-cutoff CTP functional remains a separate condition.
"""
from dataclasses import dataclass
from functools import lru_cache
from math import exp, log, pi, sqrt

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit, loggamma
from scipy.integrate import quad
import sympy as sp
from numpy.polynomial.legendre import leggauss

from .nsc_gauge_source import Exterior


def horizon_covariance(frequency,surface_gravity,reflection,incoming_occupation=0.):
    """Occupied covariance for k_z=-frequency in a real horizon spin frame.

    Incoming parent modes are empty. The exterior up-mode reflection uses
    H_ext=-i sigma3 d_x+V sigma1. Interior H=m sigma1+(k_z/b)sigma3.
    The cross entry is i*sqrt(f*(1-f))*conj(R), not merely sqrt(f*(1-f))*R*.
    """
    if frequency<0 or surface_gravity<=0:
        raise ValueError("nonnegative frequency and positive surface gravity required")
    if abs(reflection)>1+1e-10:raise ValueError("current-normalized reflection required")
    if not 0<=incoming_occupation<=1:raise ValueError("fermionic incoming occupation must be in [0,1]")
    f=float(expit(-2*pi*frequency/surface_gravity));s=sqrt(f*(1-f))
    return np.array([[1-f,1j*s*np.conj(reflection)],
                     [-1j*s*reflection,f*abs(reflection)**2+(1-abs(reflection)**2)*incoming_occupation]],dtype=complex)


def rindler_reflection(frequency,surface_gravity,mass_coefficient):
    """Known decaying Dirac/Morse solution for V=M exp(kappa*x)."""
    if frequency<0 or surface_gravity<=0 or mass_coefficient<=0:
        raise ValueError("positive Rindler scales and nonnegative frequency required")
    nu=frequency/surface_gravity
    phase=2*loggamma(.5+1j*nu).imag-2*nu*log(mass_coefficient/(2*surface_gravity))
    return -1j*np.exp(1j*phase)


def horizon_frame(frequency,surface_gravity,angular,proper_distance,tortoise,timelike=False,terms=10):
    """Uniform Dirac Frobenius frame at a small horizon distance.

    Exterior columns carry e^(+/-i omega x) and preserve sigma3 current.
    Interior columns use k_z=-omega, carry the same phases, and are unitary.
    Modified-Bessel versus Bessel series distinguishes the two domains.
    """
    if frequency<0 or surface_gravity<=0 or angular<0 or proper_distance<0:
        raise ValueError("declared nonnegative channel and positive horizon scale required")
    nu=frequency/surface_gravity;z=angular*proper_distance
    sign=-1 if timelike else 1
    def series(parameter):
        value=1.+0j;term=1.+0j
        for j in range(1,terms):
            term*=sign*z*z/(4*j*(parameter+j-1));value+=term
        return value
    phase=np.exp(1j*frequency*tortoise)
    a=phase*series(.5+1j*nu)
    b=(1j if not timelike else -1j)*z/(1+2j*nu)*phase*series(1.5+1j*nu)
    return np.array([[a,-b.conjugate() if timelike else b.conjugate()],
                     [b,a.conjugate()]],dtype=complex)


def rindler_ratio_control(frequency,surface_gravity,mass_coefficient,minimum=1e-10,maximum=35.,tolerance=1e-11):
    """Independent first-order ratio solution of the same decaying Rindler mode."""
    nu=frequency/surface_gravity
    if maximum<=nu:raise ValueError("start in the decaying region")
    initial=nu/maximum-1j*sqrt(1-(nu/maximum)**2)
    def rhs(y,state):
        z=exp(y);r=state[0]
        return [1j*z*(1+r*r)-2j*nu*r]
    run=solve_ivp(rhs,(log(maximum),log(minimum)),[initial],method="DOP853",rtol=tolerance,atol=tolerance*.02,max_step=.1)
    if not run.success:raise RuntimeError(run.message)
    phase=2*nu*log(minimum*surface_gravity/mass_coefficient)
    numerical=run.y[0,-1]*np.exp(1j*phase)
    return {"reflection":complex(numerical),"analytic":complex(rindler_reflection(frequency,surface_gravity,mass_coefficient)),
            "error":float(abs(numerical-rindler_reflection(frequency,surface_gravity,mass_coefficient))),
            "current_defect":float(abs(abs(numerical)**2-1))}


def outgoing_ratio_series(frequency,angular,radius,order=8):
    """Outgoing large-parent-radius expansion of f_minus/f_plus.

    Coefficients solve A r'=i V(1+r²)-2i omega r with the recorded asymptotic A.
    This fixes a phase-compatible boundary condition without a new metric fit.
    """
    if frequency<=0 or angular<0 or radius<=0 or order<2:
        raise ValueError("positive frequency/radius and a resolved nonnegative channel required")
    A=np.zeros(order+1);A[0]=1
    for n in range((order+1)//2):
        power=2*n+1
        if power<=order:A[power]=-6*(-1)**n/((2*n+1)*(2*n+3))
    root=np.zeros(order+1);root[0]=1
    for n in range(1,order+1):root[n]=(A[n]-sum(root[j]*root[n-j] for j in range(1,n)))/2
    sphere=np.zeros(order+1);sphere[0]=1;coefficient=1.
    for n in range(1,order//2+1):
        coefficient*=(-.5-(n-1))/n;sphere[2*n]=coefficient
    V=np.zeros(order+1);V[1:]=angular*np.convolve(root,sphere)[:order]
    ratio=np.zeros(order+1,dtype=complex)
    for n in range(1,order+1):
        derivative=-sum(A[n-j-1]*j*ratio[j] for j in range(1,n))
        square=np.convolve(ratio,ratio)[:order+1];square[0]+=1
        forcing=-1j*np.convolve(V,square)[n]
        ratio[n]=-(derivative+forcing)/(2j*frequency)
    return complex(sum(ratio[n]/radius**n for n in range(1,order+1)))


@dataclass(frozen=True)
class ParentDirac:
    horizon_rho: float
    surface_gravity: float

    def _near_coefficients(self):
        h=self.horizon_rho
        second=-6*np.arctan(1/h)+6*h/(1+h*h)
        third=12/(1+h*h)**2
        fourth=-48*h/(1+h*h)**3
        return 2*self.surface_gravity,second,third,fourth

    def A_from_offset(self,offset):
        """Avoid losing the tiny horizon offset when forming rho_h+offset."""
        if offset<=0:raise ValueError("exterior offset required")
        if offset<1e-4:
            first,second,third,fourth=self._near_coefficients()
            return first*offset+second*offset**2/2+third*offset**3/6+fourth*offset**4/24
        return Exterior(self.horizon_rho,self.surface_gravity).A(self.horizon_rho+offset)

    def near_tortoise(self,offset,origin=0.):
        """x=log(offset)/(2*kappa)+regular part with regular part zero at horizon."""
        first,second,third,_=self._near_coefficients()
        constant=-second/(2*first*first)
        slope=second*second/(4*first**3)-third/(6*first*first)
        if offset==0:raise ValueError("nonzero signed horizon offset required")
        return log(abs(offset))/first+constant*offset+slope*offset**2/2+origin

    def reflection(self,frequency,angular,outer_factor=60.,inner_offset=1e-10,tolerance=1e-10,order=8,origin=0.):
        """Actual first-order exterior reflection, preserving phase and current."""
        if frequency<=0 or angular<0:raise ValueError("positive frequency and nonnegative angular channel required")
        end=max(60.,outer_factor*max(1.,angular)/frequency)
        initial=outgoing_ratio_series(frequency,angular,end,order)
        initial_current=1-abs(initial)**2
        if initial_current<=0:raise ArithmeticError("outer asymptotic channel is unresolved")
        def rhs(y,state):
            offset=exp(y);rho=self.horizon_rho+offset;A=self.A_from_offset(offset)
            V=angular*sqrt(A)/sqrt(1+rho*rho);factor=offset/A;r=state[0]
            return [factor*(1j*V*(1+r*r)-2j*frequency*r),factor*V*r.imag]
        run=solve_ivp(rhs,(log(end-self.horizon_rho),log(inner_offset)),[initial,0j],
                      method="DOP853",rtol=tolerance,atol=tolerance*.02,max_step=.15)
        if not run.success:raise RuntimeError(run.message)
        raw_ratio=run.y[0,-1];log_amplitude=float(run.y[1,-1].real)
        rh=sqrt(1+self.horizon_rho**2)
        distance=sqrt(2*inner_offset/self.surface_gravity)/rh
        frame=horizon_frame(frequency,self.surface_gravity,angular,distance,self.near_tortoise(inner_offset,origin))
        coefficients=np.linalg.solve(frame,np.array([1.,raw_ratio]))
        raw=coefficients[1]/coefficients[0]
        transmission=float(initial_current*exp(-2*log_amplitude)/abs(coefficients[0])**2)
        if not 0<=transmission<=1+1e-8:raise ArithmeticError("unresolved transmission normalization")
        magnitude=sqrt(max(0.,1-transmission))
        reflection=raw/abs(raw)*magnitude if abs(raw)>0 else 0j
        return {"frequency":frequency,"angular":angular,"reflection":complex(reflection),
                "raw_reflection":complex(raw),"transmission":transmission,
                "current_defect":float(abs(abs(raw)**2+transmission-1)),
                "log_amplitude":log_amplitude,"outer_radius":end,"inner_offset":inner_offset,
                "tortoise_origin":origin,"function_evaluations":int(run.nfev)}

    @property
    def horizon_q(self):
        return pi-np.arctan(1/self.horizon_rho)

    def interior_W(self,delta_q):
        """The same conformal child, evaluated stably at its past horizon."""
        from .nsc_angular_stress import profile_jets
        if delta_q<=0:raise ValueError("positive interior horizon distance required")
        if delta_q<1e-4:
            q=self.horizon_q
            second=-6*np.sin(2*q)-2*np.cos(2*q)
            third=-12*np.cos(2*q)+4*np.sin(2*q)
            fourth=24*np.sin(2*q)+8*np.cos(2*q)
            return 2*self.surface_gravity*delta_q+second*delta_q**2/2-third*delta_q**3/6+fourth*delta_q**4/24
        return float(profile_jets(self.horizon_q-delta_q)[0])

    def initial_covariance(self,frequency,angular,reflection,delta_q=1e-10,origin=0.,incoming_occupation=0.):
        rh2=1+self.horizon_rho**2
        signed_offset=-rh2*delta_q+self.horizon_rho*rh2*delta_q**2
        x=self.near_tortoise(signed_offset,origin)
        distance=sqrt(2*delta_q/self.surface_gravity)
        F=horizon_frame(frequency,self.surface_gravity,angular,distance,x,timelike=True)
        C=horizon_covariance(frequency,self.surface_gravity,reflection,incoming_occupation)
        return F@C@F.conj().T

    def interior_covariances(self,frequency,angular,reflection,points=(.1,.5,pi/2),
                             delta_q=1e-10,origin=0.,incoming_occupation=0.,tolerance=2e-11):
        """Transport the matched covariance, not a new instantaneous vacuum."""
        initial=self.initial_covariance(frequency,angular,reflection,delta_q,origin,incoming_occupation)
        targets=sorted(points,reverse=True)
        times=np.array([log(self.horizon_q-q) for q in targets])
        if times[0]<=log(delta_q):raise ValueError("points must lie to the future of the starting collar")
        sigma1=np.array([[0.,1.],[1.,0.]]);sigma3=np.diag([1.,-1.])
        def rhs(y,state):
            delta=exp(y);W=self.interior_W(delta)
            H=delta*(angular/sqrt(W)*sigma1-frequency/W*sigma3)
            U=state.reshape(2,2)
            return (-1j*H@U).reshape(4)
        run=solve_ivp(rhs,(log(delta_q),float(times[-1])),np.eye(2,dtype=complex).reshape(4),
                      t_eval=times,method="DOP853",rtol=tolerance,atol=tolerance*.02,max_step=.1)
        if not run.success:raise RuntimeError(run.message)
        rows=[]
        for q,state in zip(targets,run.y.T):
            U=state.reshape(2,2);C=U@initial@U.conj().T
            rows.append({"q":float(q),"covariance":C,"trace":float(np.trace(C).real),
                         "eigenvalues":np.linalg.eigvalsh(C),
                         "unitarity_defect":float(np.linalg.norm(U.conj().T@U-np.eye(2)))})
        return rows


def frequency_grid(edges=(0.,.25,1.,4.,16.,64.),points_per_interval=24):
    nodes,weights=leggauss(points_per_interval)
    frequencies=[];quadrature=[]
    for low,high in zip(edges,edges[1:]):
        frequencies.extend(low+(nodes+1)*(high-low)/2)
        quadrature.extend(weights*(high-low)/2)
    return np.array(frequencies),np.array(quadrature)


def reflection_grid(background,frequencies,angular_max,phase_cutoff=4.,tolerance=1e-11,existing=None):
    """Retain the phase wherever its thermal covariance can contribute."""
    reflection=np.zeros((angular_max,len(frequencies)),complex)
    transmission=np.ones_like(reflection.real)
    worst=0.;calls=0;start=0
    if existing is not None:
        old_frequency=np.asarray(existing["frequencies"])
        shared=min(len(old_frequency),len(frequencies))
        assert np.array_equal(old_frequency[:shared],frequencies[:shared])
        assert existing["horizon_rho"]==background.horizon_rho and existing["surface_gravity"]==background.surface_gravity
        assert existing["phase_cutoff"]==phase_cutoff and existing["tolerance"]==tolerance
        start=min(angular_max,len(existing["reflection"]))
        reflection[:start,:shared]=np.asarray(existing["reflection"])[:start,:shared]
        transmission[:start,:shared]=np.asarray(existing["transmission"])[:start,:shared]
        if shared<len(frequencies):assert np.all(frequencies[shared:]>phase_cutoff)
        worst=existing["maximum_current_defect"]
    for j in range(start,angular_max):
        for i,w in enumerate(frequencies):
            if w>phase_cutoff:continue
            data=background.reflection(float(w),j+1,tolerance=tolerance)
            reflection[j,i]=data["reflection"];transmission[j,i]=data["transmission"]
            worst=max(worst,data["current_defect"]);calls+=1
    return {"reflection":reflection,"transmission":transmission,"maximum_current_defect":worst,
            "phase_cutoff":phase_cutoff,"new_solves":calls,"frequencies":frequencies,
            "horizon_rho":background.horizon_rho,"surface_gravity":background.surface_gravity,"tolerance":tolerance}


def parent_angular_source(background,angular_max=8,points_per_interval=24,
                          frequency_edges=(0.,.25,1.,4.,16.,64.),points=(.1,.2,.5,1.,pi/2),
                          steps=12000,delta_q=1e-10,reflection_data=None,phase_cutoff=4.,origin=0.):
    """Canonical Unruh boundary candidate and thermal incoming control.

    This computes a prescribed characteristic-state realization. It does
    not assert the global Hadamard extension hypotheses or a finite-Lambda
    CTP completion. No auxiliary reference source is recalculated.
    """
    from .nsc_angular_stress import (_unitary_step,profile_jets,adiabatic_coefficients,
        curvature_squared_tensor,static_cylinder_density,physical_source)
    frequency,weights=frequency_grid(frequency_edges,points_per_interval)
    data=reflection_grid(background,frequency,angular_max,phase_cutoff,existing=reflection_data)
    R=np.asarray(data["reflection"])*np.exp(2j*frequency[None,:]*origin)
    T=np.asarray(data["transmission"])
    f=expit(-2*pi*frequency/background.surface_gravity)
    trace_unruh=1-f[None,:]*T
    masses=np.arange(1,angular_max+1,dtype=float)[:,None]
    rh2=1+background.horizon_rho**2
    offset=-rh2*delta_q+background.horizon_rho*rh2*delta_q**2
    x=background.near_tortoise(offset,origin)
    distance=sqrt(2*delta_q/background.surface_gravity)
    plus=np.empty_like(R);minus=np.empty_like(R);in_plus=np.empty_like(R);in_minus=np.empty_like(R)
    for j in range(angular_max):
        for i,w in enumerate(frequency):
            F=horizon_frame(float(w),background.surface_gravity,j+1,distance,x,timelike=True)
            v=F@np.array([sqrt(1-f[i]),-1j*sqrt(f[i])*R[j,i]])
            incoming=F@np.array([0.,sqrt(f[i]*T[j,i])])
            plus[j,i],minus[j,i]=v
            in_plus[j,i],in_minus[j,i]=incoming
    # Fixed rotation maps sigma1 to sigma2, reusing the checked unitary step.
    plus*=np.exp(-1j*pi/4);minus*=np.exp(1j*pi/4)
    in_plus*=np.exp(-1j*pi/4);in_minus*=np.exp(1j*pi/4)
    begin=log(delta_q);ordered=sorted(points,reverse=True)
    targets=[log(background.horizon_q-q) for q in ordered]
    step_bound=(targets[-1]-begin)/steps
    a=(3-2*sqrt(3))/12;b=(3+2*sqrt(3))/12;c1=.5-sqrt(3)/6;c2=.5+sqrt(3)/6
    power=float(sum(4*masses[:,0]/pi*(T@(weights*frequency*f))))
    records=[]
    def generator(y):
        delta=exp(y);W=background.interior_W(delta)
        return masses*delta/sqrt(W),-frequency[None,:]*delta/W
    for q,end in zip(ordered,targets):
        count=max(1,int(np.ceil((end-begin)/step_bound)));step=(end-begin)/count
        for i in range(count):
            y=begin+i*step;y1,z1=generator(y+c1*step);y2,z2=generator(y+c2*step)
            for hy,hz in ((b*y1+a*y2,b*z1+a*z2),(a*y1+b*y2,a*z1+b*z2)):
                plus,minus=_unitary_step(plus,minus,hy,hz,step)
                in_plus,in_minus=_unitary_step(in_plus,in_minus,hy,hz,step)
        canonical_plus=plus*np.exp(1j*pi/4);canonical_minus=minus*np.exp(-1j*pi/4)
        incoming_plus=in_plus*np.exp(1j*pi/4);incoming_minus=in_minus*np.exp(-1j*pi/4)
        jets=profile_jets(q);scale=sqrt(jets[0]);p=-frequency[None,:]/scale
        omega=np.sqrt(masses*masses+p*p)
        cosine=masses/np.sqrt(2*omega*(omega-p));sine=np.sqrt((omega-p)/(2*omega))
        beta=cosine*canonical_plus+sine*canonical_minus
        alpha=-sine*canonical_plus+cosine*canonical_minus
        beta_in=cosine*incoming_plus+sine*incoming_minus
        alpha_in=-sine*incoming_plus+cosine*incoming_minus
        E2,P2,E4,P4=adiabatic_coefficients(p/masses,1.,jets)
        factor=masses[:,0]/(pi*pi*scale)
        henergy,hz=curvature_squared_tensor(jets)
        rho_static=static_cylinder_density();harmonic=float(__import__('mpmath').euler)/2
        rho_local=harmonic*henergy/(480*pi*pi);pz_local=-harmonic*hz/(480*pi*pi)
        W,W1,W2,W3,W4=jets;R2=-W2;box=-W*W4-W1*W3
        anomaly=(-(R2-2)**2/60-11*R2/90-box/30)/(16*pi*pi)
        states={}
        for name,occupation,tangent,trace in (
            ("Unruh_boundary_candidate",abs(beta)**2,-2*np.real(beta.conj()*alpha),trace_unruh),
            ("thermal_incoming_control",abs(beta)**2+abs(beta_in)**2,
             -2*np.real(beta.conj()*alpha+beta_in.conj()*alpha_in),np.ones_like(trace_unruh))):
            excess=2*omega*occupation+omega*(1-trace)
            pressure=(p*p/omega**2)*excess+(p*masses/omega)*tangent
            er=factor*((excess-E2/masses-E4/masses**3)@weights)
            pr=factor*((pressure-P2/masses-P4/masses**3)@weights)
            row={"q":q,"rho_coordinate":-1/np.tan(q),"sphere_radius":1/np.sin(q),"W_jets":jets.tolist(),
                 "bar_rho":rho_static+rho_local+float(sum(er)),
                 "bar_p_parallel":-rho_static+pz_local+float(sum(pr)),
                 "bar_trace":anomaly,"angular_energy_remainders":er.tolist(),"angular_pressure_remainders":pr.tolist()}
            row["bar_p_sphere"]=(row["bar_rho"]-row["bar_p_parallel"]-anomaly)/2
            tensor=physical_source(row)
            tensor["T_hatT_hatz"]=-power*np.sin(q)**4/(4*pi*W) if name=="Unruh_boundary_candidate" else 0.
            states[name]={"barred":row,"physical":tensor}
        norm_defect=float(np.max(abs(abs(plus)**2+abs(minus)**2-trace_unruh)))
        records.append({"q":q,"states":states,"mode_norm_defect":norm_defect})
        begin=end
    return {"parameters":{"angular_max":angular_max,"points_per_interval":points_per_interval,
             "frequency_edges":list(frequency_edges),"steps":steps,"delta_q":delta_q,"phase_cutoff":phase_cutoff},
            "parent_power":power,"maximum_scattering_current_defect":data["maximum_current_defect"],
            "rows":list(reversed(records)),"reflection_data":data}


@lru_cache(maxsize=1)
def _sixth_order_functions():
    """The existing normalized Dirac adiabatic recursion, continued to order6."""
    p,m=sp.symbols("p m",real=True)
    hs=sp.symbols("H0:6",real=True)
    omega=sp.sqrt(m*m+p*p)
    derivative=lambda v:-hs[0]*p*sp.diff(v,p)+sum(hs[i+1]*sp.diff(v,hs[i]) for i in range(5))
    theta=m*p*hs[0]/omega**2
    A=[sp.Integer(1)];B=[sp.Integer(0)];C=[sp.Integer(0)]
    for n in range(1,7):
        A.append(sp.factor(-sum(A[i]*A[n-i]+B[i]*B[n-i]+C[i]*C[n-i] for i in range(1,n))/2))
        B.append(sp.factor(derivative(C[n-1])/(2*omega)))
        C.append(sp.factor(-(derivative(B[n-1])+theta*A[n-1])/(2*omega)))
    energy=sp.factor(-omega*A[6]);pressure=sp.factor(-p*(A[6]*p-B[6]*m)/omega)
    w=sp.symbols("w0:7",real=True)
    dq=lambda v:-sp.sqrt(w[0])*sum(sp.diff(v,w[i])*w[i+1] for i in range(6))
    rates=[-w[1]/(2*sp.sqrt(w[0]))]
    for _ in range(5):rates.append(sp.factor(dq(rates[-1])))
    return (sp.lambdify((p,m,*hs),(energy,pressure),"numpy",cse=True),
            sp.lambdify(w,rates,"numpy",cse=True))


def sixth_order_angular_coefficient(q):
    """Leading kappa^-3 coefficient of the already-subtracted 4D angular sum.

    This is used only after checking the actual high angular channels against
    it. It is a numerical tail approximation, not a new physical counterterm.
    """
    from .nsc_angular_stress import profile_jets
    jets=list(profile_jets(q))
    jets += [48*np.cos(2*q)-16*np.sin(2*q),-96*np.sin(2*q)-32*np.cos(2*q)]
    f,rates=_sixth_order_functions();h=rates(*jets)
    values=[quad(lambda p:f(p,1.,*h)[i]/(pi*pi),0,np.inf,epsabs=1e-9,epsrel=1e-9)[0] for i in (0,1)]
    return {"q":q,"angular_energy_coefficient":float(values[0]),"angular_pressure_coefficient":float(values[1])}


def thermal_phase_omission_bound(angular_max,surface_gravity,frequency_cutoff,minimum_b):
    """Bound finite-grid source error from omitting exponentially small R coherences.

    The theoretical boundary state retains R at all frequencies. This bound
    is for the stated finite angular calculation and finite observation points.
    """
    alpha=pi/surface_gravity
    def integrals(a):
        e=exp(-a*frequency_cutoff)
        return e/a,e*(frequency_cutoff/a+1/a**2)
    a0,a1=integrals(alpha);b0,b1=integrals(2*alpha)
    first=angular_max*(angular_max+1)/2
    second=angular_max*(angular_max+1)*(2*angular_max+1)/6
    density=(second*(2*a0+b0)+first*(2*a1+b1)/minimum_b)/(pi*pi*minimum_b)
    return {"diagonal_source_bound":density,"power_bound":4*first*b1/pi}
