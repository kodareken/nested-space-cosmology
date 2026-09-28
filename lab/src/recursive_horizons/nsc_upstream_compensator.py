"""Common raw-KS metric directions and their short-slab operator response.

Fundamental identity matrices are operator initial data, not quantum states.
The three real matching coordinates multiply one common geometric bump.
"""
from fractions import Fraction as Q
import math
import numpy as np
from scipy.integrate import solve_ivp

from .nsc_compatible_history_geometry import CompatibleRadiusDirection
from .nsc_retarded_radial_response import ks_generator
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_lorentzian import geometry
from .nsc_transmitting_resolvent import profile
from .nsc_transmitting_dirac_domain import I2,S1,S2,S3

PAULI=np.array([S1,S2,S3]);DIRECTIONS=('original_radius','raw_N','raw_a','raw_r')
MAX_STEP=1/2000


def bump(rho):
    # Exact rational center51/50 and halfwidth1/200; avoid redefining the
    # mathematical support by the rounded sum of two binary parameters.
    return profile(200*float(rho)-204)[0]


def vertices(rho,energies,mass,angular):
    E=np.asarray(energies);beta=float(geometry(rho)[0]);a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
    if E.ndim!=1 or np.iscomplexobj(E) or not np.isfinite(E).all():raise ValueError('real source/output energies required')
    mean=(E[:,None]+E[None,:])/2
    VN=-mass*S1+angular/r*S2-mean[:,:,None,None]/a*S3
    Va=mean[:,:,None,None]/a**2*S3
    Vr=np.broadcast_to(-angular/r**2*S2,VN.shape)
    return np.array([VN,Va,Vr])


def vertex_identities():
    """Exact raw-H derivatives plus the symmetric Fourier momentum identity."""
    import sympy as sp
    N,a,r,m,ell,Eo,Ei,k=sp.symbols('N a r m ell Eo Ei k',real=True,nonzero=True)
    s1=sp.Matrix([[0,1],[1,0]]);s2=sp.Matrix([[0,-sp.I],[sp.I,0]]);s3=sp.diag(1,-1)
    H=N*(-m*s1+ell/r*s2+k/a*s3);mid=-(Eo+Ei)/2
    expected=(-m*s1+ell/r*s2+mid/a*s3,-mid/a**2*s3,-ell/r**2*s2)
    errors=[sp.diff(H,v).subs({N:1,k:mid})-e for v,e in zip((N,a,r),expected)]
    residual=[str(sp.simplify(v)) for matrix in errors for v in matrix]
    # The two derivatives in -i/2{b(z),partial_z} give (Ei+Eo)/2.
    weyl=sp.simplify(Ei+(Eo-Ei)/2-(Ei+Eo)/2)
    return {'raw_vertex_residuals':residual,'Weyl_midpoint_residual':str(weyl)}


def operator_response(upstream,energies,mass,angular,what,*,rtol,atol):
    E=np.asarray(energies);W=np.asarray(what,complex)
    if E.shape!=(4,) or W.shape!=(4,4) or upstream<=1.025:
        raise ValueError('four energy labels and upstream node beyond support required')
    U0=np.broadcast_to(I2,(4,2,2)).copy();B0=np.zeros((4,4,4,2,2),complex)
    normal=CompatibleRadiusDirection(lambda z,n:1. if n==0 else 0.,lambda z,n:0.,.007,.03)
    T1=chart_coordinates(1.)[0]
    def rhs(rho,state):
        U=state[:16].reshape(4,2,2);B=state[16:].reshape(4,4,4,2,2)
        G=ks_generator(rho,E,mass,angular)
        beta=float(geometry(rho)[0]);a=np.sqrt(beta*beta-1)
        V=vertices(rho,E,mass,angular)
        original=normal.value(chart_coordinates(float(rho))[0]-T1,0.)*V[2]
        directions=np.concatenate((original[None],bump(rho)*V))
        dU=np.einsum('iab,ibc->iac',G,U)
        dB=np.einsum('oab,doibc->doiac',G,B)+(1j/a)*W[None,:,:,None,None]*np.einsum('doiab,ibc->doiac',directions,U)
        return np.concatenate((dU.ravel(),dB.ravel()))
    solution=solve_ivp(rhs,(float(upstream),1.),np.concatenate((U0.ravel(),B0.ravel())),
                       method='DOP853',rtol=rtol,atol=atol,max_step=MAX_STEP)
    if not solution.success:raise ArithmeticError(solution.message)
    last=solution.y[:,-1];U=last[:16].reshape(4,2,2);B=last[16:].reshape(4,4,4,2,2)
    K=np.einsum('doiab,icb->doiac',B,U.conj())
    return U,B,K,{'rtol':rtol,'atol':atol,'maximum_step':MAX_STEP,'accepted_steps':len(solution.t)-1,'function_evaluations':solution.nfev}


def fit_selected(K,index=3,tolerance=3e-11):
    diagonal=np.asarray(K)[:,index,index]
    anti=float(np.max(abs(diagonal+diagonal.swapaxes(-1,-2).conj())))
    trace=float(np.max(abs(np.trace(diagonal,axis1=-2,axis2=-1))))
    coefficients=.5j*np.einsum('jab,dba->dj',PAULI,diagonal)
    imaginary=float(np.max(abs(coefficients.imag)))
    if max(anti,trace,imaginary)>=tolerance:raise ArithmeticError('selected operators are not traceless anti-Hermitian at tolerance')
    matrix=coefficients[1:].real.T;target=-coefficients[0].real
    chosen=np.linalg.solve(matrix,target)
    return chosen,{'matrix':matrix.tolist(),'target':target.tolist(),'singular_values':np.linalg.svd(matrix,compute_uv=False).tolist(),
                  'condition_number':float(np.linalg.cond(matrix)),'selected_antiHermitian':anti,'selected_trace':trace,
                  'Pauli_coefficient_imaginary':imaginary,'linear_system_residual':float(np.max(abs(matrix@chosen-target)))}


def combine(K,coefficients):return np.einsum('d,doiab->oiab',np.r_[1.,coefficients],K)


def covariance_response(K,covariance):
    C=np.asarray(covariance)
    return np.array([[K[o,i]@C[i]+C[o]@K[i,o].conj().T for i in range(4)] for o in range(4)])


def small_amplitude_range(coefficients):
    c=list(map(lambda v:abs(Q(float(v))),coefficients))
    limit=min(Q(1,100)/(1+c[0]),Q(1,125)/(1+c[1]),Q(1,2)/(1+c[2]+Q(3,100)))
    # |dN|<.01, |da|<.008 imply N>.99,a>.99*a0; a0>.8.
    numerator=Q(99,100)**2*(1+Q(4,5)**2)-Q(101,100)**2
    assert numerator>0 and limit*c[0]<Q(1,100) and limit*c[1]<Q(1,125)
    assert limit*(c[2]+Q(3,100))<Q(1,2)
    return {'strict_abs_eta_upper':{'numerator':limit.numerator,'denominator':limit.denominator},
            'conservative_display_upper':math.nextafter(float(limit),-math.inf),
            'q_PG_squared_lower_on_compensator':{'numerator':numerator.numerator,'denominator':numerator.denominator},
            'radius_lower':.5,'physical_amplitude_selected':False,
            'bounds':'b,w<=1; rawN>.99; rawa>.99*a0; outside compensator N,a unchanged; r>.5; F=q_PG²>0'}
