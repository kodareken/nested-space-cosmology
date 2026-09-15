"""Oscillatory packet integrals of massive high-energy mode expansions.

Mass and gradient coefficients are computed from the same PG Dirac operator.
Only the characteristic coordinate is shared with the massless geometry
atlas. No LLL covariance or massless mode preparation is assigned here.
"""
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import make_interp_spline

from .nsc_lorentzian import geometry
from .nsc_pg_high_energy import riccati_coefficients
from .nsc_pg_lll_preparation import LLLFlowAtlas,radial_packet
from .nsc_transmitting_dirac_domain import MODE_TO_CURRENT


class SplineFourier:
    """Fourier integral of the interpolating spline, including its jumps.

    Repeated integration by parts is exact for the spline. The highest
    derivative jumps at internal knots must be retained. Interpolation error
    remains a separate numerical check.
    """
    def __init__(self,x,values,degree=7):
        x=np.asarray(x,float);values=np.asarray(values,complex)
        if np.any(np.diff(x)<=0):raise ValueError('ordered Fourier coordinate required')
        self.a=float(x[0]);self.b=float(x[-1]);self.degree=degree
        spline=make_interp_spline(x,values,k=degree)
        self.left=np.array([spline.derivative(j)(self.a) for j in range(degree+1)])
        self.right=np.array([spline.derivative(j)(self.b) for j in range(degree+1)])
        knots=np.unique(spline.t);knots=knots[(knots>=self.a)&(knots<=self.b)]
        mids=(knots[:-1]+knots[1:])/2
        top=spline.derivative(degree)(mids)
        self.jumps=np.diff(top,axis=0);self.knots=knots[1:-1]
        self.value_shape=values.shape[1:]

    def evaluate(self,energy):
        if energy<=0:raise ValueError('positive oscillation frequency required')
        result=np.zeros(self.value_shape,complex);denominator=1j*energy
        for j in range(self.degree+1):
            result+=(-1)**j*(self.right[j]*np.exp(1j*energy*self.b)-self.left[j]*np.exp(1j*energy*self.a))/denominator**(j+1)
        if len(self.knots):
            phase=np.exp(1j*energy*self.knots)
            result+=(-1)**(self.degree+1)*np.tensordot(phase,self.jumps,axes=(0,0))/denominator**(self.degree+1)
        return result


@dataclass
class TailBranch:
    branch: int
    reference: float
    reference_beta: float
    reference_series: np.ndarray
    packets: dict

    def normalization(self,energy):
        S=self.reference_series@(1/(2*energy))**np.arange(1,len(self.reference_series)+1)
        b=self.reference_beta
        if self.branch==1:
            q=1j*(1-b)*S;current=((1-b)+(-1-b)*abs(q)**2)/abs(1-b)
        else:
            q=-1j*(1+b)*S.conjugate();current=((1-b)*abs(q)**2+(-1-b))/abs(-1-b)
        return 1/np.sqrt(abs(current))


class MassiveTailPackets:
    """Vacuum-tail branch projections, with explicit asymptotic order/grid.

    The three columns are exterior outgoing, interior partner, and global
    incoming. Cross-source thermal and reflection corrections are not silently
    declared zero: their control belongs to the tail acceptance record.
    """
    def __init__(self,preparation,mass,angular,*,order=16,points_per_unit=128,minimum_energy=160.):
        if mass<0 or order<4 or points_per_unit<16:raise ValueError('declared channel and resolved tail approximation required')
        self.preparation=preparation;self.mass=mass;self.angular=angular
        self.order=order;self.points_per_unit=points_per_unit;self.minimum_energy=minimum_energy
        self.atlas=LLLFlowAtlas(preparation.horizon_rho,preparation.surface_gravity)
        self.branches=[]
        self.branches.append(self._branch(1,np.linspace(-1.,1.,2*points_per_unit+1),0.,('parent','child','child_bulk')))
        self.branches.append(self._branch(1,np.linspace(3.,4.,points_per_unit+1),3.,('exterior_bulk',)))
        self.branches.append(self._branch(-1,np.linspace(-1.,4.,5*points_per_unit+1),0.,('parent','child','child_bulk','exterior_bulk')))

    def _branch(self,branch,rho,reference,kinds):
        N=self.order;coeff,_=riccati_coefficients(rho,self.mass,self.angular,N)
        beta=geometry(rho)[0];U=self.angular/np.sqrt(1+rho*rho)+1j*self.mass
        if branch==1:
            phase_coeff=-U.conj()[:,None]*coeff;speed=1-beta
            ratio_coeff=1j*speed[:,None]*coeff
        else:
            phase_coeff=U[:,None]*coeff.conj();speed=-1-beta
            ratio_coeff=-1j*(1+beta)[:,None]*coeff.conj()
        primitive=make_interp_spline(rho,phase_coeff,k=7).antiderivative()
        theta=primitive(rho)-primitive(reference)
        exp_coeff=np.zeros((len(rho),N+1),complex);exp_coeff[:,0]=1.
        for n in range(1,N+1):
            exp_coeff[:,n]=sum(1j*k*theta[:,k-1]*exp_coeff[:,n-k] for k in range(1,n+1))/n
        mixed=np.zeros_like(exp_coeff)
        for n in range(1,N+1):mixed[:,n]=sum(ratio_coeff[:,j-1]*exp_coeff[:,n-j] for j in range(1,n+1))
        character=np.stack((exp_coeff,mixed) if branch==1 else (mixed,exp_coeff),axis=-1)
        current=np.einsum('ij,nkj->nki',MODE_TO_CURRENT,character)/np.sqrt(abs(speed))[:,None,None]
        flow=self.atlas.coordinate(branch,rho);packets={}
        domains={'parent':(0.,1.),'child':(-1.,0.),'child_bulk':(-1.,0.),'exterior_bulk':(3.,4.)}
        for kind in kinds:
            a,b=domains[kind];mask=(rho>=a)&(rho<=b)
            x=flow[mask];value=current[mask]*abs(speed[mask])[:,None,None]*radial_packet(rho[mask],kind)[0][:,None,None]
            value[[0,-1]]=0. # exact compact packet endpoints
            if x[0]>x[-1]:x=x[::-1];value=value[::-1]
            packets[kind]=SplineFourier(x,value)
        refcoeff,_=riccati_coefficients([reference],self.mass,self.angular,N)
        return TailBranch(branch,reference,float(geometry(reference)[0]),refcoeff[0],packets)

    def project(self,energy):
        if energy<self.minimum_energy:raise ValueError('energy is below the declared asymptotic tail range')
        powers=(1/(2*energy))**np.arange(self.order+1)
        output=np.zeros((8,3),complex)
        positions={'parent':0,'child':2,'child_bulk':4,'exterior_bulk':6}
        for branch,column in zip(self.branches,(1,0,2)):
            factor=branch.normalization(energy)
            for kind,transform in branch.packets.items():
                value=np.tensordot(powers,transform.evaluate(energy),axes=(0,0))*factor
                start=positions[kind];output[start:start+2,column]=value
        return output
