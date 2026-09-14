"""Global PG Cauchy covariance of the already specified transparent LLL state.

Three characteristic sectors belong to the same massless Dirac field:
outgoing exterior, outgoing interior, and incoming whole line. Their original
horizon/incoming covariance is transported by a unitary flow-coordinate map.
Massive retained channels are not assigned this preparation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp, quad
from scipy.special import expit

from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import MODE_TO_CURRENT
from .nsc_unruh_state import ParentDirac


def gauss(left, right, points):
    x,w=leggauss(points)
    return left+(x+1)*(right-left)/2, w*(right-left)/2


def radial_packet(rho, kind):
    """Original parent/child probes, plus two orthogonal bulk probe shapes."""
    rho=np.asarray(rho)
    if kind=='parent': t=rho; dt=1.; norm=np.sqrt(30.)
    elif kind in ('child','child_bulk'): t=-rho; dt=-1.; norm=np.sqrt(30. if kind=='child' else 210.)
    elif kind=='exterior_bulk': t=rho-3.; dt=1.; norm=np.sqrt(30.)
    else: raise ValueError('unknown fixed spatial probe')
    inside=(t>0)&(t<1)
    value=norm*t*(1-t); derivative=norm*(1-2*t)*dt
    if kind=='child_bulk':
        derivative=norm*((1-2*t)*(2*t-1)+2*t*(1-t))*dt
        value*=2*t-1
    return np.where(inside,value,0.),np.where(inside,derivative,0.)


class LLLFlowAtlas:
    """Same-horizon affine matching; coordinate origins are common translations.

    xplus = log(|rho-rho_h|/rho_h)/kappa + regular(rho).
    One smooth regular primitive is shared by both horizon sides. No relative
    phase or Cauchy-history duration is selected independently.
    """

    def __init__(self, horizon_rho, surface_gravity):
        self.h=float(horizon_rho);self.k=float(surface_gravity)
        self.background=ParentDirac(self.h,self.k)
        first,second,third,fourth=self.background._near_coefficients()
        # Series of a_plus=1-sqrt(1-A) through delta^4, from the owned A series.
        A=np.array([0.,first,second/2,third/6,fourth/24])
        a=.5*A
        for power,coefficient in ((2,1/8),(3,1/16),(4,5/128)):
            poly=np.array([1.])
            for _ in range(power): poly=np.convolve(poly,A)[:5]
            a+=coefficient*poly
        self.regular_series=np.array([
            -a[2]/self.k**2,
            a[2]**2/self.k**3-a[3]/self.k**2,
            -a[2]**3/self.k**4+2*a[2]*a[3]/self.k**3-a[4]/self.k**2,
        ])
        def rhs(rho,_):
            beta=float(geometry(rho)[0]);d=rho-self.h
            if abs(d)<1e-3:
                regular=float(np.polynomial.polynomial.polyval(d,self.regular_series))
            else: regular=1/(1-beta)-1/(self.k*d)
            return [regular,1/(-1-beta)]
        self.rhs=rhs
        self.segments={}
        for side,end in (('negative',-1.),('positive',4.)):
            sol=solve_ivp(rhs,(0.,end),[0.,0.],method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True)
            if not sol.success:raise ArithmeticError(sol.message)
            self.segments[side]=sol

    def velocity(self,sign,rho): return sign-geometry(rho)[0]

    def coordinate(self,sign,rho):
        x=np.asarray(rho,float)
        if sign not in (-1,1) or not np.isfinite(x).all() or (sign==1 and np.any(x==self.h)):
            raise ValueError('finite flow coordinate away from its horizon endpoint required')
        primitive=np.empty_like(x)
        for mask,side in (((x<0)&(x>=-1),'negative'),((x>=0)&(x<=4),'positive')):
            if np.any(mask): primitive[mask]=self.segments[side].sol(x[mask])[0 if sign==1 else 1]
        for index in np.argwhere((x < -1.) | (x > 4.)):
            key=tuple(index);value=float(x[key]);boundary=-1. if value<0 else 4.
            side='negative' if value<0 else 'positive';component=0 if sign==1 else 1
            primitive[key]=self.segments[side].sol(boundary)[component]+quad(
                lambda r:self.rhs(r,None)[component],boundary,value,epsabs=2e-12,epsrel=2e-12)[0]
        if sign==1:primitive=primitive+np.log(abs(x-self.h)/self.h)/self.k
        return primitive


@dataclass(frozen=True)
class PGLLLPreparation:
    horizon_rho: float
    surface_gravity: float
    inheritance_ratio: float

    def __post_init__(self):
        if min(self.horizon_rho,self.surface_gravity,self.inheritance_ratio)<=0:
            raise ValueError('locked positive horizon and inheritance data required')

    @property
    def beta_out(self):return 2*np.pi/self.surface_gravity

    @property
    def beta_in(self):return self.beta_out/self.inheritance_ratio

    def spectral_covariance(self,frequency):
        """Original CH on ext+/int+ and inherited incoming occupation, all real E."""
        frequency=float(frequency)
        f=float(expit(-self.beta_out*frequency));n=float(expit(-self.beta_in*frequency))
        x=.5*self.beta_out*frequency
        coherence=np.exp(-abs(x))/(1+np.exp(-2*abs(x)))
        return np.array([[f,-1j*coherence,0.],[1j*coherence,1-f,0.],[0.,0.,n]],complex)

    def seed_restriction(self,frequency):
        # Inner outgoing and throughgoing incoming; no copy of stored Cj.
        return self.spectral_covariance(frequency)[np.ix_([1,2],[1,2])]

    def off_diagonal_kernel(self,rho,rho_prime,*,atlas=None):
        """Complete equal-PG-time kernel at distinct points, excluding the delta.

        The contact contribution is I delta(rho-rho')/2 and is included by
        project(). This method refuses coincident points instead of returning
        a made-up finite local vacuum value.
        """
        if rho==rho_prime:raise ValueError('coincident covariance is a distribution; smear the complete kernel')
        atlas=atlas or LLLFlowAtlas(self.horizon_rho,self.surface_gravity)
        r=MODE_TO_CURRENT
        Pp=np.outer(r[:,0],r[:,0].conj());Pm=np.outer(r[:,1],r[:,1].conj())
        xp=atlas.coordinate(1,rho)-atlas.coordinate(1,rho_prime)
        xm=atlas.coordinate(-1,rho)-atlas.coordinate(-1,rho_prime)
        ap=np.sqrt(abs(atlas.velocity(1,rho)*atlas.velocity(1,rho_prime)))
        am=np.sqrt(abs(atlas.velocity(-1,rho)*atlas.velocity(-1,rho_prime)))
        exterior=rho>self.horizon_rho;other_exterior=rho_prime>self.horizon_rho
        if exterior==other_exterior:
            plus=(-1j if exterior else 1j)/(2*self.beta_out*ap*np.sinh(np.pi*xp/self.beta_out))
        else:
            plus=(-1j if exterior else 1j)/(2*self.beta_out*ap*np.cosh(np.pi*xp/self.beta_out))
        minus=-1j/(2*self.beta_in*am*np.sinh(np.pi*xm/self.beta_in))
        return plus*Pp+minus*Pm

    def project(self,points=96):
        """Compress the full continuum covariance onto seven orthonormal probes.

        Ordering: original parent current-spin2, original child current-spin2,
        child_bulk current-spin2, exterior_bulk sigma2+ (one component).
        Vacuum delta/PV part is retained; there is no finite-frequency CAR fill.
        """
        atlas=LLLFlowAtlas(self.horizon_rho,self.surface_gravity)
        l,wl=gauss(-1.,0.,points);p,wp=gauss(0.,1.,points)
        x=np.r_[l,p];w=np.r_[wl,wp]
        kinds=('parent','child','child_bulk')
        values=np.array([radial_packet(x,k)[0] for k in kinds])
        derivatives=np.array([radial_packet(x,k)[1] for k in kinds])
        gram=np.einsum('n,in,jn->ij',w,values,values)
        sectors=[]
        for sign,beta,occupation_sign in ((1,self.beta_out,1.),(-1,self.beta_in,-1.)):
            flow=atlas.coordinate(sign,x);speed=atlas.velocity(sign,x)
            dx=flow[:,None]-flow[None,:]
            denominator=np.sinh(np.pi*dx/beta)*np.sqrt(abs(speed[:,None]*speed[None,:]))
            sector=.5*gram.astype(complex)
            for i in range(3):
                for j in range(i+1,3):
                    numerator=values[i,:,None]*values[j,None,:]-values[j,:,None]*values[i,None,:]
                    np.fill_diagonal(denominator,1.)
                    quotient=numerator/denominator
                    diagonal=beta/np.pi*np.sign(speed)*(derivatives[i]*values[j]-values[i]*derivatives[j])
                    np.fill_diagonal(quotient,diagonal)
                    part=occupation_sign*1j/(4*beta)*np.einsum('n,m,nm',w,w,quotient)
                    sector[i,j]+=part;sector[j,i]+=part.conjugate()
            sectors.append(sector)
        plus,minus=sectors
        e,we=gauss(3.,4.,points);fe=radial_packet(e,'exterior_bulk')[0]
        ex=atlas.coordinate(1,e);ix=atlas.coordinate(1,x)
        cross=-1j/(2*self.beta_out)/np.cosh(np.pi*(ex[:,None]-ix[None,:])/self.beta_out)
        cross/=np.sqrt(abs(atlas.velocity(1,e)[:,None]*atlas.velocity(1,x)[None,:]))
        ext_inner=np.einsum('e,n,e,in,en->i',we,w,fe,values,cross)
        c=np.zeros((7,7),complex);g=np.zeros((7,7),complex)
        r=MODE_TO_CURRENT;Pp=np.outer(r[:,0],r[:,0].conj());Pm=np.outer(r[:,1],r[:,1].conj())
        for i in range(3):
            for j in range(3):
                c[2*i:2*i+2,2*j:2*j+2]=plus[i,j]*Pp+minus[i,j]*Pm
                g[2*i:2*i+2,2*j:2*j+2]=gram[i,j]*np.eye(2)
            c[6,2*i:2*i+2]=ext_inner[i]*r[:,0].conj()
            c[2*i:2*i+2,6]=c[6,2*i:2*i+2].conj()
        g[6,6]=np.dot(we,fe*fe);c[6,6]=.5*g[6,6]
        return {'covariance':c,'CAR_gram':g,'lesser_equal_time':1j*c,
                'greater_equal_time':-1j*(g-c),'Keldysh_equal_time':-1j*(g-2*c),
                'atlas':atlas,'probe_kinds':kinds+('exterior_bulk',),
                'complement_is_exhausted_by_three_extra_probes':False}

    def packet_retarded_resolvent(self,z,points=64):
        """Independent characteristic integral compared with locked LLL R(z)."""
        z=complex(z)
        if z.imag<=0:raise ValueError('upper-half-plane control required')
        atlas=LLLFlowAtlas(self.horizon_rho,self.surface_gravity)
        intervals=((0.,1.,'parent'),(-1.,0.,'child'))
        result=np.zeros((4,4),complex)
        for sign,col in ((1,0),(-1,1)):
            spin=np.outer(MODE_TO_CURRENT[:,col],MODE_TO_CURRENT[:,col].conj())
            for i,(il,iu,ik) in enumerate(intervals):
                rr,ww=gauss(il,iu,points);f=radial_packet(rr,ik)[0]
                for j,(jl,ju,jk) in enumerate(intervals):
                    value=0j
                    for rho,weight,fr in zip(rr,ww,f):
                        start=max(rho,jl)
                        if start>=ju:continue
                        ss,vv=gauss(start,ju,points);fs=radial_packet(ss,jk)[0]
                        delay=atlas.coordinate(sign,rho)-atlas.coordinate(sign,ss)
                        kernel=1j*np.exp(1j*z*delay)/np.sqrt(abs(atlas.velocity(sign,rho)*atlas.velocity(sign,ss)))
                        value+=weight*fr*np.dot(vv,fs*kernel)
                    result[2*i:2*i+2,2*j:2*j+2]+=value*spin
        return result


def require_full_retained_preparation(channel_index):
    if channel_index!=0:
        raise ValueError('LLL flow preparation does not reconstruct the 32 massive angular/compact channels')
