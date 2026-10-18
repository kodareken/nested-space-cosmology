"""Complete cylinder covariance on supplied one-sided retarded maps.

This is integration of the existing analytic root kernel, not another
field solver or state selector. Polynomial response error and angular
completion are external. Keep the noncoincident split until singular
terms have been combined. The map is in the chiral eta frame.
"""
from dataclasses import dataclass,field
import math

import numpy as np
from numpy.polynomial.legendre import leggauss,legvander,legder
from scipy.integrate import quad_vec
from scipy.special import gammaln,lqmn


@dataclass(frozen=True)
class NullResponse:
    """L f = sum_a D_a f_a(L_a)+D'_a f'_a(L_a)+integral G_a(s) f_a(s) ds.

    Each G is represented by its Legendre coefficients on its OWN closed
    leg interval, not by a zero extension or a finite CAR projector.
    """
    lengths: tuple[float,float]
    green: tuple[np.ndarray,np.ndarray]
    direct: np.ndarray
    direct_derivative: np.ndarray=field(default_factory=lambda:np.zeros((2,2),complex))

    def __post_init__(self):
        lengths=np.asarray(self.lengths,dtype=float)
        direct=np.asarray(self.direct,dtype=complex)
        derivative=np.asarray(self.direct_derivative,dtype=complex)
        green=tuple(np.asarray(c,dtype=complex) for c in self.green)
        if (lengths.shape!=(2,) or not np.isfinite(lengths).all() or np.any(lengths<=0)
                or direct.shape!=(2,2) or not np.isfinite(direct).all()
                or derivative.shape!=(2,2) or not np.isfinite(derivative).all()
                or len(green)!=2 or any(c.ndim!=2 or c.shape[1]!=2 or not len(c)
                    or not np.isfinite(c).all() for c in green)):
            raise ValueError('positive leg lengths, two finite (order,2) Green arrays and 2x2 direct map required')
        object.__setattr__(self,'lengths',tuple(lengths))
        object.__setattr__(self,'green',green)
        object.__setattr__(self,'direct',direct)
        object.__setattr__(self,'direct_derivative',derivative)

    def value(self,leg,s,*,derivative=0):
        z=2*np.asarray(s)/self.lengths[leg]-1
        coefficients=self.green[leg]
        if derivative==1:
            coefficients=legder(coefficients,axis=0)*2/self.lengths[leg]
        elif derivative!=0:
            raise ValueError('only the first response derivative is supplied')
        return np.einsum('...n,ni->...i',legvander(z,len(coefficients)-1),coefficients)

    def hilbert(self,leg,s,*,derivative=0):
        """PV integral G(t)/(s-t) dt, with the actual one-sided endpoints."""
        z=2*float(s)/self.lengths[leg]-1
        if abs(z)==1:
            raise ValueError('endpoint Hilbert value needs a combined finite-part limit')
        # Legendre Q includes the real principal value on (-1,1).
        if derivative not in (0,1):
            raise ValueError('only the first Hilbert derivative is supplied')
        Q=lqmn(0,len(self.green[leg])-1,z)[derivative][0]
        return (2 if derivative==0 else 4/self.lengths[leg])*(Q@self.green[leg])


def _moments(response,leg,count):
    """Polynomial and logarithmic moments without a quadrature of log(0)."""
    coefficients=response.green[leg];n=len(coefficients)
    ordinary=np.zeros((count,n));logarithmic=np.zeros_like(ordinary)
    harmonic=np.r_[0.,np.cumsum(1/np.arange(1,count+n+1))]
    for j in range(count):
        for m in range(n):
            if j>=m:
                value=math.exp(2*gammaln(j+1)-gammaln(j-m+1)-gammaln(j+m+2))
                ordinary[j,m]=value
                logarithmic[j,m]=value*(2*harmonic[j]-harmonic[j-m]-harmonic[j+m+1])
            else:
                logarithmic[j,m]=(-1.)**(m-j-1)*math.exp(
                    2*gammaln(j+1)+gammaln(m-j)-gammaln(j+m+2))
    L=response.lengths[leg]
    endpoint_derivative=response.direct_derivative[:,leg]/L
    return (response.direct[:,leg]+np.arange(count)[:,None]*endpoint_derivative+L*(ordinary@coefficients),
            endpoint_derivative+L*(logarithmic@coefficients))


def _direct_green(left,right):
    result=np.zeros((2,2),complex)
    for a in range(2):
        s=left.lengths[a]
        value=right.value(a,s) if s<right.lengths[a] else np.zeros(2)
        conjugate_row=(value/2-1j*right.hilbert(a,s)/(2*math.pi)).conj()
        result+=np.outer(left.direct[:,a],conjugate_row)
        if np.any(left.direct_derivative[:,a]):
            value=right.value(a,s,derivative=1) if s<right.lengths[a] else np.zeros(2)
            row=(value/2-1j*right.hilbert(a,s,derivative=1)/(2*math.pi)).conj()
            result+=np.outer(left.direct_derivative[:,a],row)
    return result


def pair_cylinder_responses(left,right,angular_mass,*,series_terms=24,epsabs=1e-10):
    """Full-root, noncoincident 2x2 covariance on the supplied polynomial maps.

    Returns a Bessel-series norm bound and a numerical PV quadrature
    estimate separately. Neither covers response interpolation, floating
    arithmetic or the angular sum. Derivative maps retain their endpoint
    f' evaluation explicitly. No momentum
    cutoff occurs. The current stable series route requires |mu|sqrt(LR)<=4.
    """
    mu=float(angular_mass)
    if (not isinstance(left,NullResponse) or not isinstance(right,NullResponse)
            or not np.isfinite(mu) or mu==0 or isinstance(series_terms,bool)
            or not isinstance(series_terms,int) or series_terms<8
            or not np.isfinite(epsabs) or epsabs<=0):
        raise ValueError('two null responses, nonzero real mass and finite positive integration controls required')
    if any(a==b for a,b in zip(left.lengths,right.lengths)):
        raise ValueError('noncoincident endpoints required before a finite-part limit')
    result=_direct_green(left,right)+_direct_green(right,left).conj().T
    pv_error=0.
    for a in range(2):
        L,R=left.lengths[a],right.lengths[a]
        D,F=left.direct[:,a],left.direct_derivative[:,a]
        E,G=right.direct[:,a],right.direct_derivative[:,a]
        result+=1j/(2*math.pi)*(np.outer(D,E.conj())/(L-R)
            +(np.outer(D,G.conj())-np.outer(F,E.conj()))/(L-R)**2
            -2*np.outer(F,G.conj())/(L-R)**3)
        nodes,weights=leggauss(max(len(left.green[a]),len(right.green[a])))
        s=(nodes+1)*min(L,R)/2
        result+=np.einsum('s,si,sj->ij',weights*min(L,R)/4,
            left.value(a,s),right.value(a,s).conj())
        integral,error=quad_vec(lambda t:np.outer(left.value(a,t),right.hilbert(a,t).conj()),
            0,L,points=[R] if R<L else None,epsabs=epsabs,epsrel=epsabs,limit=300)
        result+=1j*integral/(2*math.pi);pv_error+=float(error)/(2*math.pi)
    tail=0.
    for a,b in ((0,1),(1,0)):
        L,R=left.lengths[a],right.lengths[b]
        z=abs(mu)*math.sqrt(L*R)
        if z>4:
            raise ValueError('Bessel series requires |mu|sqrt(LR)<=4; no unstable high-sector result returned')
        M,J=_moments(left,a,series_terms);N,K=_moments(right,b,series_terms)
        c=np.euler_gamma+math.log(z/2);coefficient=1.;harmonic=0.
        for j in range(series_terms):
            result-=mu/(2*math.pi)*coefficient*((harmonic-c)*np.outer(M[j],N[j].conj())
                -(np.outer(J[j],N[j].conj())+np.outer(M[j],K[j].conj()))/2)
            harmonic+=1/(j+1);coefficient*=z*z/(4*(j+1)**2)
        gL=L*sum(np.linalg.norm(v) for v in left.green[a])
        gR=R*sum(np.linalg.norm(v) for v in right.green[b])
        dL=np.linalg.norm(left.direct_derivative[:,a])/L
        dR=np.linalg.norm(right.direct_derivative[:,b])/R
        mL=np.linalg.norm(left.direct[:,a])+gL+series_terms*dL
        mR=np.linalg.norm(right.direct[:,b])+gR+series_terms*dR
        logL=gL+dL;logR=gR+dR
        ratio=z*z/(4*(series_terms+1)**2)
        h=harmonic+abs(c);cross=mL*dR+dL*mR;pair=dL*dR
        polynomial=np.array([h*mL*mR+(logL*mR+mL*logR)/2,
            h*cross+mL*mR/(series_terms+1)+(logL*dR+dL*logR)/2,
            h*pair+cross/(series_terms+1),pair/(series_terms+1)])
        powers=np.array([1/(1-ratio),ratio/(1-ratio)**2,
            ratio*(1+ratio)/(1-ratio)**3,ratio*(1+4*ratio+ratio*ratio)/(1-ratio)**4])
        tail+=abs(mu)/(2*math.pi)*coefficient*float(polynomial@powers)
    return dict(kernel=result,bessel_series_tail_bound=tail,pv_quadrature_estimate=pv_error,
        scope='complete_root_covariance_on_supplied_polynomial_responses')
