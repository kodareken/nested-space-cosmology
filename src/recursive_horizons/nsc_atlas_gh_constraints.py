"""Independent ADM diagnostics from raw regular-frame spacetime jets.

The reference-orthonormal spatial frame obeys [E_i,E_j]=2*r*eps_ijk E_k.
No GH evolution equation or spherical reduction is used here.  Raw second
jets are ORDERED E_i E_j g, not symmetrized coordinate Hessians.  Scalar
momenta supplied to this interface are coordinate-time derivatives.
"""
import math

import numpy as np
from numba import njit,prange

EPS=np.zeros((3,3,3))
for _i,_j,_k in ((0,1,2),(1,2,0),(2,0,1)):
    EPS[_i,_j,_k]=1.
    EPS[_j,_i,_k]=-1.


def _array(value,name,shape):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape!=shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+' must be a finite real array of shape '+str(shape))
    return np.ascontiguousarray(raw,dtype=float)


@njit(cache=True,nogil=True,parallel=True)
def _diagnostics(g,gi,hi,v,dg,d2g,dv,p,dq,W,r,kappa):
    count=len(g)
    lapse=np.empty(count);shift=np.empty((count,3));Kout=np.empty((count,3,3))
    scalarR=np.empty(count);rho=np.empty(count);momentum_density=np.empty((count,3))
    Hout=np.empty(count);Mout=np.empty((count,3))
    bracket=2*r*EPS
    for point in prange(count):
        h=g[point,1:,1:]
        dh=dg[point,:,1:,1:]
        inv=hi[point]
        N=1/math.sqrt(-gi[point,0,0]);lapse[point]=N
        beta=g[point,0,1:]
        for i in range(3):
            shift[point,i]=0.
            for j in range(3):shift[point,i]+=inv[i,j]*beta[j]

        # Koszul formula in the physical metric, including frame brackets.
        lower=np.empty((3,3,3));dlower=np.empty((3,3,3,3))
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    value=dh[i,j,k]+dh[j,i,k]-dh[k,i,j]
                    for ell in range(3):
                        value+=(-h[i,ell]*bracket[j,k,ell]
                                +h[j,ell]*bracket[k,i,ell]+h[k,ell]*bracket[i,j,ell])
                    lower[k,i,j]=value/2
                    for m in range(3):
                        dvalue=(d2g[point,m,i,j+1,k+1]+d2g[point,m,j,i+1,k+1]
                                -d2g[point,m,k,i+1,j+1])
                        for ell in range(3):
                            dvalue+=(-dh[m,i,ell]*bracket[j,k,ell]
                                     +dh[m,j,ell]*bracket[k,i,ell]+dh[m,k,ell]*bracket[i,j,ell])
                        dlower[m,k,i,j]=dvalue/2
        dinv=np.empty((3,3,3))
        for m in range(3):
            for i in range(3):
                for j in range(3):
                    value=0.
                    for a in range(3):
                        for b in range(3):value-=inv[i,a]*dh[m,a,b]*inv[b,j]
                    dinv[m,i,j]=value
        Gamma=np.empty((3,3,3));dGamma=np.empty((3,3,3,3))
        for ell in range(3):
            for i in range(3):
                for j in range(3):
                    Gamma[ell,i,j]=0.
                    for k in range(3):Gamma[ell,i,j]+=inv[ell,k]*lower[k,i,j]
                    for m in range(3):
                        dGamma[m,ell,i,j]=0.
                        for k in range(3):
                            dGamma[m,ell,i,j]+=(dinv[m,ell,k]*lower[k,i,j]+inv[ell,k]*dlower[m,k,i,j])
        Ricci=np.zeros((3,3))
        for k in range(3):
            for j in range(3):
                for i in range(3):
                    Ricci[k,j]+=dGamma[i,i,j,k]-dGamma[j,i,i,k]
                    for m in range(3):
                        Ricci[k,j]+=(Gamma[i,i,m]*Gamma[m,j,k]-Gamma[i,j,m]*Gamma[m,i,k]
                                      -bracket[i,j,m]*Gamma[i,m,k])
        R=0.
        for i in range(3):
            for j in range(3):R+=inv[i,j]*Ricci[i,j]
        scalarR[point]=R

        dN=np.empty(3)
        for m in range(3):
            dinverse00=0.
            for a in range(4):
                for b in range(4):dinverse00-=gi[point,0,a]*dg[point,m,a,b]*gi[point,b,0]
            dN[m]=N**3*dinverse00/2
        A=np.empty((3,3));K=np.empty((3,3));dK=np.empty((3,3,3))
        for i in range(3):
            for j in range(3):
                A[i,j]=v[point,i+1,j+1]-dg[point,i,0,j+1]-dg[point,j,0,i+1]
                for k in range(3):A[i,j]+=(Gamma[k,i,j]+Gamma[k,j,i])*beta[k]
                K[i,j]=-A[i,j]/(2*N)
                for m in range(3):
                    dA=(dv[point,m,i+1,j+1]-d2g[point,m,i,0,j+1]-d2g[point,m,j,0,i+1])
                    for k in range(3):
                        dA+=((dGamma[m,k,i,j]+dGamma[m,k,j,i])*beta[k]
                              +(Gamma[k,i,j]+Gamma[k,j,i])*dg[point,m,0,k+1])
                    dK[m,i,j]=-dA/(2*N)+A[i,j]*dN[m]/(2*N*N)
        Kout[point]=K
        Ku=np.empty((3,3));dKu=np.empty((3,3,3))
        for j in range(3):
            for i in range(3):
                Ku[j,i]=0.
                for k in range(3):Ku[j,i]+=inv[j,k]*K[k,i]
                for m in range(3):
                    dKu[m,j,i]=0.
                    for k in range(3):dKu[m,j,i]+=dinv[m,j,k]*K[k,i]+inv[j,k]*dK[m,k,i]
        trace=Ku[0,0]+Ku[1,1]+Ku[2,2];K2=0.
        for i in range(3):
            for j in range(3):K2+=Ku[i,j]*Ku[j,i]
        normal=np.empty(2);energy=W[point]
        for a in range(2):
            normal[a]=p[point,a]
            for i in range(3):normal[a]-=shift[point,i]*dq[point,i,a]
            normal[a]/=N
            energy+=normal[a]**2/2
            for i in range(3):
                for j in range(3):energy+=inv[i,j]*dq[point,i,a]*dq[point,j,a]/2
        rho[point]=energy
        for i in range(3):
            momentum_density[point,i]=0.
            for a in range(2):momentum_density[point,i]-=normal[a]*dq[point,i,a]
        Hout[point]=R+trace*trace-K2-2*kappa*energy
        for i in range(3):
            value=-kappa*momentum_density[point,i]
            for j in range(3):
                value+=dKu[j,j,i]-dKu[i,j,j]
                for ell in range(3):value+=Gamma[j,j,ell]*Ku[ell,i]-Gamma[ell,j,i]*Ku[j,ell]
            Mout[point,i]=value
    return lapse,shift,Kout,scalarR,rho,momentum_density,Hout,Mout


def adm_constraints(g,v,dg,d2g,dv,q,p,dq,*,r,kappa=1.,potential):
    """Return pointwise ADM geometry, source and constraint residual arrays.

    Trailing shapes: g,v=(4,4), dg,dv=(3,4,4), d2g=(3,3,4,4),
    q,p=(2,), dq=(3,2). All share g's optional leading batch dimensions.
    v=partial_t g, dg_i=E_i g, d2g_ij=E_i E_j g, dv_i=E_i v,
    p=partial_t(S,chi), where S=sqrt(Ks)*(s-s_star) is canonical.
    ``potential`` is the FULL W(S,chi), including its vacuum value. No
    additional cosmological-constant contribution is inserted here.
    """
    raw=np.asarray(g)
    if raw.ndim<2 or raw.shape[-2:]!=(4,4):raise ValueError('g needs trailing shape (4,4)')
    batch=raw.shape[:-2]
    g=_array(g,'g',batch+(4,4));v=_array(v,'v',batch+(4,4))
    dg=_array(dg,'dg',batch+(3,4,4));dv=_array(dv,'dv',batch+(3,4,4))
    d2g=_array(d2g,'ordered d2g',batch+(3,3,4,4))
    q=_array(q,'canonical fields q',batch+(2,));p=_array(p,'coordinate field rates p',batch+(2,))
    dq=_array(dq,'frame field derivatives dq',batch+(3,2))
    for a,name in ((g,'g'),(v,'v'),(dg,'dg'),(d2g,'d2g'),(dv,'dv')):
        tolerance=2e-12*max(1.,float(np.max(abs(a))))
        if np.max(abs(a-np.swapaxes(a,-1,-2)))>tolerance:
            raise ValueError(name+' must be symmetric in its metric indices')
    if isinstance(r,(bool,np.bool_)) or isinstance(kappa,(bool,np.bool_)):
        raise ValueError('finite real frame scale and positive kappa required')
    r,kappa=float(r),float(kappa)
    if not math.isfinite(r) or r<0 or not math.isfinite(kappa) or kappa<=0:
        raise ValueError('require finite r>=0 and kappa>0')
    W=np.asarray(potential)
    if W.dtype.kind not in 'iuf':raise ValueError('potential must be finite and real')
    try:W=np.ascontiguousarray(np.broadcast_to(W,batch),dtype=float)
    except ValueError as e:raise ValueError('potential must broadcast over the point batch') from e
    if not np.all(np.isfinite(W)):raise ValueError('potential must be finite')
    count=int(np.prod(batch)) if batch else 1
    flatg=g.reshape(count,4,4)
    h=flatg[:,1:,1:]
    if np.min(np.linalg.eigvalsh(h))<=0:raise ValueError('spatial metric must be positive definite')
    try:gi=np.linalg.inv(flatg);hi=np.linalg.inv(h)
    except np.linalg.LinAlgError as e:raise ValueError('metric inverse failed') from e
    if not np.all(np.isfinite(gi)) or np.any(gi[:,0,0]>=0):
        raise ValueError('time slices must have a real positive lapse')
    result=_diagnostics(flatg,gi,hi,v.reshape(count,4,4),dg.reshape(count,3,4,4),
        d2g.reshape(count,3,3,4,4),dv.reshape(count,3,4,4),p.reshape(count,2),
        dq.reshape(count,3,2),W.reshape(count),r,kappa)
    if any(not np.all(np.isfinite(value)) for value in result):
        raise ValueError('derived ADM geometry or sources are not finite')
    names=('lapse','shift','Kij','R3','rho','j','H','M')
    shapes=((),(3,),(3,3),(),(),(3,),(),(3,))
    return {name:value.reshape(batch+tail) for name,value,tail in zip(names,result,shapes)}
