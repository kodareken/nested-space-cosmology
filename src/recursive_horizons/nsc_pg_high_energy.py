"""High-energy approximation of the same massive PG mode branches.

Coefficients are obtained by applying the already owned Dirac Riccati
recursion to the NSC profiles. The A'/2 term and mass/angular phase are kept.
This is a numerical tail representation, not a new vacuum subtraction or an
assignment of the massless LLL state to a massive sector.
"""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import make_interp_spline

from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import MODE_TO_CURRENT


def _mul(a,b):
    result=np.zeros_like(a,dtype=complex)
    for n in range(a.shape[1]):
        result[:,n]=sum(a[:,j]*b[:,n-j] for j in range(n+1))
    return result


def _derivative(a):
    result=np.zeros_like(a)
    result[:,:-1]=a[:,1:]*np.arange(1,a.shape[1])
    return result


def riccati_coefficients(rho,mass,angular,order=8):
    """Return a_n and a_n' for S_+=sum a_n/(2E)^n, on real geometry.

    a_1=U_+, a_{n+1}=i(A d_rho+A'/2)a_n
       +A U_- sum_{j+k=n}a_j a_k; U_+=lambda/r+i m.
    The incoming coefficients are conjugates, not conjugated frequencies.
    Taylor jets evaluate spatial derivatives of the owned metric/profile.
    """
    x=np.atleast_1d(np.asarray(rho,float))
    if not np.isfinite(x).all() or mass<0 or order<2:
        raise ValueError('finite real geometry and a resolved existing massive channel required')
    size=order+2
    r2=np.zeros((len(x),size),complex);r2[:,0]=1+x*x;r2[:,1]=2*x;r2[:,2]=1.
    inverse=np.zeros_like(r2);inverse[:,0]=1/r2[:,0]
    for n in range(1,size):inverse[:,n]=-sum(r2[:,j]*inverse[:,n-j] for j in range(1,n+1))/r2[:,0]
    sphere=np.zeros_like(r2);sphere[:,0]=np.sqrt(inverse[:,0])
    for n in range(1,size):sphere[:,n]=(inverse[:,n]-sum(sphere[:,j]*sphere[:,n-j] for j in range(1,n)))/(2*sphere[:,0])
    angle=np.zeros_like(r2);angle[:,0]=np.arctan(x)-np.pi/2
    for n in range(1,size):angle[:,n]=inverse[:,n-1]/n
    A=3*_mul(r2,angle);A[:,0]+=1+3*x;A[:,1]+=3.
    Uplus=angular*sphere;Uplus[:,0]+=1j*mass;Uminus=Uplus.conj()
    coefficients=[Uplus];dA=_derivative(A)
    for n in range(1,order):
        last=coefficients[-1]
        product=sum((_mul(coefficients[j-1],coefficients[n-j-1]) for j in range(1,n)),start=np.zeros_like(A))
        coefficients.append(1j*(_mul(A,_derivative(last))+.5*_mul(dA,last))+_mul(_mul(A,Uminus),product))
    values=np.stack([a[:,0] for a in coefficients],axis=1)
    derivatives=np.stack([a[:,1] for a in coefficients],axis=1)
    return values,derivatives


def local_series(rho,energy,mass,angular,order=8):
    if energy<=0:raise ValueError('positive real tail energy required')
    values,derivatives=riccati_coefficients(rho,mass,angular,order)
    powers=(1/(2*energy))**np.arange(1,order+1)
    return values@powers,derivatives@powers


def resummed_local_series(rho,energy,mass,angular,order=8):
    """Gradient expansion about the exact local dispersion root.

    This retains all powers of A*(m^2+lambda^2/r^2)/E^2 in p. It changes
    numerical mode approximation only, not the fourth-order stress reference.
    """
    x=np.atleast_1d(np.asarray(rho,float));size=order+3
    r2=np.zeros((len(x),size),complex);r2[:,0]=1+x*x;r2[:,1]=2*x;r2[:,2]=1
    def inverse(a):
        b=np.zeros_like(a);b[:,0]=1/a[:,0]
        for n in range(1,size):b[:,n]=-sum(a[:,j]*b[:,n-j] for j in range(1,n+1))/a[:,0]
        return b
    def root(a):
        b=np.zeros_like(a);b[:,0]=np.sqrt(a[:,0])
        for n in range(1,size):b[:,n]=(a[:,n]-sum(b[:,j]*b[:,n-j] for j in range(1,n)))/(2*b[:,0])
        return b
    invr2=inverse(r2);sphere=root(invr2)
    angle=np.zeros_like(r2);angle[:,0]=np.arctan(x)-np.pi/2
    for n in range(1,size):angle[:,n]=invr2[:,n-1]/n
    A=3*_mul(r2,angle);A[:,0]+=1+3*x;A[:,1]+=3
    U=angular*sphere;U[:,0]+=1j*mass;Uc=U.conj()
    p2=-_mul(A,_mul(U,Uc));p2[:,0]+=energy*energy
    if np.any(p2[:,0].real<=0):raise ValueError('real propagating local dispersion required')
    p=root(p2);ep=p.copy();ep[:,0]+=energy
    orders=[_mul(U,inverse(ep))];inv2p=.5*inverse(p);dA=_derivative(A)
    for n in range(1,order+1):
        derivative=1j*(_mul(A,_derivative(orders[-1]))+.5*_mul(dA,orders[-1]))
        product=sum((_mul(orders[j],orders[n-j]) for j in range(1,n)),start=np.zeros_like(A))
        orders.append(_mul(inv2p,derivative+_mul(_mul(A,Uc),product)))
    total=sum(orders)
    return total[:,0],total[:,1]


def branch_riccati_residual(rho,energy,mass,angular,branch,order=8,*,resummed=False):
    x=np.atleast_1d(np.asarray(rho,float));S,dS=(resummed_local_series if resummed else local_series)(x,energy,mass,angular,order)
    beta,bp,_=geometry(x);vp=1-beta;vm=-1-beta
    Uplus=angular/np.sqrt(1+x*x)+1j*mass;Uminus=Uplus.conj()
    if branch==1:
        if np.any(abs(vp)<1e-10):raise ValueError('outgoing branch is evaluated on one side of its horizon')
        q=1j*vp*S;dq=1j*(-bp*S+vp*dS)
        rate=1j*energy/vp+bp/(2*vp)-1j*Uminus*S
        y=np.stack((np.ones_like(q),q),axis=1)
        dy=np.stack((rate,dq+q*rate),axis=1)
    elif branch==-1:
        S=S.conj();dS=dS.conj()
        q=-1j*(1+beta)*S;dq=-1j*(bp*S+(1+beta)*dS)
        rate=1j*energy/vm+bp/(2*vm)+1j*Uplus*S
        y=np.stack((q,np.ones_like(q)),axis=1)
        dy=np.stack((dq+q*rate,rate),axis=1)
    else:raise ValueError('characteristic branch is +1 or -1')
    My=np.stack((-1j*Uminus*y[:,1],1j*Uplus*y[:,0]),axis=1)
    residual=-1j*(np.stack((vp,vm),axis=1)*dy-.5*bp[:,None]*y)+My-energy*y
    return np.linalg.norm(residual,axis=1)/np.linalg.norm(y,axis=1)


@dataclass
class BranchEnvelope:
    rho: np.ndarray
    envelope: np.ndarray
    branch: int
    energy: float
    reference_rho: float
    current_residual: float
    equation_residual: float


def branch_envelope(rho,energy,mass,angular,branch,*,order=8,reference=0.,resummed=False):
    """Envelope multiplying exp(i E X_branch); the slow phase is retained.

    The current fixes normalization at `reference`. Its arbitrary common
    phase cancels from one-branch covariance products. Cross-source thermal
    terms must not be approximated by this routine without a separate bound.
    """
    x=np.asarray(rho,float)
    if x.ndim!=1 or len(x)<10 or np.any(np.diff(x)<=0) or not x[0]<=reference<=x[-1]:
        raise ValueError('ordered real interpolation grid containing its reference required')
    evaluator=resummed_local_series if resummed else local_series
    S,_=evaluator(x,energy,mass,angular,order);beta,bp,_=geometry(x)
    Uplus=angular/np.sqrt(1+x*x)+1j*mass;Uminus=Uplus.conj()
    if branch==1:
        speed=1-beta
        if np.any(speed*speed[0]<=0):raise ValueError('outgoing envelope must not cross its horizon')
        q=1j*speed*S;vector=np.stack((np.ones_like(q),q),axis=1);correction=-Uminus*S
    elif branch==-1:
        speed=-1-beta;S=S.conj();q=-1j*(1+beta)*S
        vector=np.stack((q,np.ones_like(q)),axis=1);correction=Uplus*S
    else:raise ValueError('characteristic branch is +1 or -1')
    primitive=make_interp_spline(x,correction,k=7).antiderivative()
    phase=primitive(x)-primitive(reference)
    envelope=np.exp(1j*phase)[:,None]*vector/np.sqrt(abs(speed))[:,None]
    Sref,_=evaluator([reference],energy,mass,angular,order);bref=float(geometry(reference)[0])
    if branch==1:
        qref=1j*(1-bref)*Sref[0];current_ref=((1-bref)+(-1-bref)*abs(qref)**2)/abs(1-bref)
    else:
        qref=-1j*(1+bref)*Sref[0].conjugate();current_ref=((1-bref)*abs(qref)**2+(-1-bref))/abs(-1-bref)
    envelope/=np.sqrt(abs(current_ref))
    current=(1-beta)*abs(envelope[:,0])**2+(-1-beta)*abs(envelope[:,1])**2
    output=envelope@MODE_TO_CURRENT.T
    return BranchEnvelope(x,output,branch,energy,reference,float(np.max(abs(current-np.sign(speed)))),
                          float(np.max(branch_riccati_residual(x,energy,mass,angular,branch,order,resummed=resummed))))
