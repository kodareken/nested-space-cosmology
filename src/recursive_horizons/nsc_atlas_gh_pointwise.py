"""Parallel pointwise covariant GH and two-scalar acceleration kernel.

Implements Eq. (24) of Lindblom--Szilagyi--Taylor, arXiv:1312.0701v2,
with the trace-reversed canonical matter source. Components use a static
orthonormal S3 reference frame: Gamma_bar^a_(i b)=r epsilon_(i b a).
Raw d2[i,j] means e_i(e_j(component)), not a coordinate Hessian.

This kernel supplies the local evolution equation. Derivative accuracy,
Einstein constraints and time integration belong to the caller. Invalid
spatial metrics, lapse or arithmetic return NaNs and valid=False; no
physical fields are clipped or repaired. Constraint damping vanishes when
the GH constraint is zero and is not a physical dissipation model.
"""
import math

import numpy as np
from numba import njit, prange

_SQRT84=math.sqrt(84.)
_LOG_MAX=math.log(np.finfo(np.float64).max)


@njit(cache=True,nogil=True,inline='always')
def _epsilon3(i,j,k):
    if i==j or j==k or i==k:return 0.
    if (i==0 and j==1 and k==2) or (i==1 and j==2 and k==0) or (i==2 and j==0 and k==1):return 1.
    return -1.


@njit(cache=True,nogil=True)
def _reference_connection(r):
    connection=np.zeros((3,4,4))
    for i in range(3):
        for a in range(3):
            for b in range(3):
                connection[i,a+1,b+1]=r*_epsilon3(i,b,a)
    return connection


@njit(cache=True,nogil=True)
def _slot_action(T,connection):
    out=np.zeros((4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                out[a,b]+=connection[c,a]*T[c,b]+connection[c,b]*T[a,c]
    return out


@njit(cache=True,nogil=True)
def _adm_inverse(g):
    """Cofactor inverse and physical normal; no small-lapse clamp."""
    inv=np.empty((4,4));normal=np.empty(4)
    if not np.isfinite(g).all():return inv,normal,False
    for a in range(4):
        for b in range(a):
            if g[a,b]!=g[b,a]:return inv,normal,False
    A,B,C=g[1,1],g[2,2],g[3,3]
    D,E,F=g[1,2],g[1,3],g[2,3]
    minor=A*B-D*D
    det=A*(B*C-F*F)-D*(D*C-E*F)+E*(D*F-B*E)
    if A<=0 or minor<=0 or det<=0 or not math.isfinite(det):return inv,normal,False
    hi=np.empty((3,3))
    hi[0,0]=(B*C-F*F)/det;hi[1,1]=(A*C-E*E)/det;hi[2,2]=(A*B-D*D)/det
    hi[0,1]=(E*F-D*C)/det;hi[1,0]=hi[0,1]
    hi[0,2]=(D*F-B*E)/det;hi[2,0]=hi[0,2]
    hi[1,2]=(D*E-A*F)/det;hi[2,1]=hi[1,2]
    shift=np.zeros(3)
    for i in range(3):
        for j in range(3):shift[i]+=hi[i,j]*g[0,j+1]
    lapse2=-g[0,0]
    for i in range(3):lapse2+=g[0,i+1]*shift[i]
    if lapse2<=0 or not math.isfinite(lapse2):return inv,normal,False
    lapse=math.sqrt(lapse2)
    inv[0,0]=-1/lapse2;normal[0]=1/lapse
    for i in range(3):
        inv[0,i+1]=shift[i]/lapse2;inv[i+1,0]=inv[0,i+1]
        normal[i+1]=-shift[i]/lapse
        for j in range(3):inv[i+1,j+1]=hi[i,j]-shift[i]*shift[j]/lapse2
    if not np.isfinite(inv).all():return inv,normal,False
    return inv,normal,True


@njit(cache=True,nogil=True)
def _potential(q,sstar):
    volume=sstar+q[0]/_SQRT84
    if volume<0 and -24*volume>_LOG_MAX:return np.nan,np.nan,np.nan
    e12=math.exp(-12*volume);e14=math.exp(-14*volume);e24=math.exp(-24*volume)
    mass2=400*e12;chi=q[1]
    V=7*e12-7.5*e14+(25/32)*e24
    Vs=-84*e12+105*e14-(75/4)*e24
    W=V+mass2*chi*chi/2
    WS=(Vs-6*mass2*chi*chi)/_SQRT84
    Wchi=mass2*chi
    return W,WS,Wchi


@njit(cache=True,nogil=True)
def _node(g,v,dg,d2g,dv,q,p,dq,d2q,dp,connection,r,H0,Hdot,sstar,kappa,gamma0):
    accg=np.full((4,4),np.nan);accq=np.full(2,np.nan);Ccov=np.full(4,np.nan)
    inv,normal,valid=_adm_inverse(g)
    if not valid:return accg,accq,Ccov,False
    if (not np.isfinite(v).all() or not np.isfinite(dg).all() or not np.isfinite(d2g).all()
        or not np.isfinite(dv).all() or not np.isfinite(q).all() or not np.isfinite(p).all()
        or not np.isfinite(dq).all() or not np.isfinite(d2q).all() or not np.isfinite(dp).all()):
        return accg,accq,Ccov,False
    W,WS,Wchi=_potential(q,sstar)
    if not math.isfinite(W) or not math.isfinite(WS) or not math.isfinite(Wchi):
        return accg,accq,Ccov,False

    Dg=np.empty((4,4,4));Dg[0]=v
    for i in range(3):Dg[i+1]=dg[i]-_slot_action(g,connection[i])
    DD=np.zeros((4,4,4,4))
    for i in range(3):
        mixed=dv[i]-_slot_action(v,connection[i])
        DD[0,i+1]=mixed;DD[i+1,0]=mixed
        for j in range(3):
            value=d2g[i,j]-_slot_action(dg[i],connection[j])-_slot_action(Dg[j+1],connection[i])
            for k in range(3):value-=connection[i,k+1,j+1]*Dg[k+1]
            DD[i+1,j+1]=value

    lower=np.empty((4,4,4));upper=np.zeros((4,4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):lower[a,b,c]=(Dg[b,a,c]+Dg[c,a,b]-Dg[a,b,c])/2
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):upper[a,b,c]+=inv[a,d]*lower[d,b,c]
    delta=np.zeros(4)
    for a in range(4):
        for b in range(4):
            for c in range(4):delta[a]+=inv[b,c]*upper[a,b,c]
    for a in range(4):
        Ccov[a]=H0 if a==0 else 0.
        for b in range(4):Ccov[a]+=g[a,b]*delta[b]
    normalC=0.
    for a in range(4):normalC+=normal[a]*Ccov[a]
    normal_lower=np.zeros(4);normal_lower[0]=-1/normal[0]

    Dq=np.empty((4,2));Dq[0]=p
    for i in range(3):Dq[i+1]=dq[i]
    trace_inverse_spatial=inv[1,1]+inv[2,2]+inv[3,3]
    for a in range(4):
        for b in range(a,4):
            rhs=2*H0*upper[0,a,b]
            if a==0 and b==0:rhs-=2*Hdot
            # -2 g^cd Rbar^e_cd(a g_b)e with the stated S3 convention.
            for e in range(1,4):
                if a>0:
                    contracted=inv[a,e]
                    if a==e:contracted-=trace_inverse_spatial
                    rhs-=r*r*contracted*g[b,e]
                if b>0:
                    contracted=inv[b,e]
                    if b==e:contracted-=trace_inverse_spatial
                    rhs-=r*r*contracted*g[a,e]
            quadratic=0.
            for c in range(4):
                for d in range(4):
                    for e in range(4):
                        for f in range(4):
                            quadratic+=inv[c,d]*inv[e,f]*(Dg[e,c,a]*Dg[f,b,d]-lower[a,c,e]*lower[b,d,f])
            rhs+=2*quadratic
            matter=W*g[a,b]
            for field in range(2):matter+=Dq[a,field]*Dq[b,field]
            rhs-=2*kappa*matter
            rhs+=gamma0*(normal_lower[b]*Ccov[a]+normal_lower[a]*Ccov[b]-g[a,b]*normalC)
            known=0.
            for c in range(4):
                for d in range(4):
                    if c!=0 or d!=0:known+=inv[c,d]*DD[c,d,a,b]
            result=(rhs-known)/inv[0,0]
            accg[a,b]=result;accg[b,a]=result

    potential_gradient=np.array([WS,Wchi])
    for field in range(2):
        rhs=potential_gradient[field]
        for i in range(3):
            rhs-=2*inv[0,i+1]*dp[i,field]
            for j in range(3):
                covariant=d2q[i,j,field]
                for k in range(3):covariant-=connection[i,k+1,j+1]*dq[k,field]
                rhs-=inv[i+1,j+1]*covariant
        for a in range(4):rhs+=delta[a]*Dq[a,field]
        accq[field]=rhs/inv[0,0]
    if not np.isfinite(accg).all() or not np.isfinite(accq).all() or not np.isfinite(Ccov).all():
        accg[:]=np.nan;accq[:]=np.nan;Ccov[:]=np.nan
        return accg,accq,Ccov,False
    return accg,accq,Ccov,True


@njit(cache=True,nogil=True,parallel=True)
def _rhs_batch(g,v,dg,d2g,dv,q,p,dq,d2q,dp,r,H0,Hdot,sstar,kappa,gamma0):
    N=g.shape[0]
    accg=np.empty((N,4,4));accq=np.empty((N,2));Ccov=np.empty((N,4));valid=np.empty(N,dtype=np.bool_)
    connection=_reference_connection(r)
    for node in prange(N):
        ag,aq,C,ok=_node(g[node],v[node],dg[node],d2g[node],dv[node],q[node],p[node],dq[node],d2q[node],dp[node],
                         connection,r,H0,Hdot,sstar,kappa,gamma0)
        accg[node]=ag;accq[node]=aq;Ccov[node]=C;valid[node]=ok
    return accg,accq,Ccov,valid


def rhs_nodes(g,v,dg,d2g,dv,q,p,dq,d2q,dp,*,r,Hsource0,HsourceDot0,sstar,kappa=1.,gamma0=0.):
    """Return (accg,accq,Ccov,valid) for regular-frame raw component jets.

    Shapes: g,v (N,4,4), dg,dv (N,3,4,4), d2g (N,3,3,4,4),
    q,p (N,2), dq,dp (N,3,2), d2q (N,3,3,2). Metric component slots
    must be symmetric. The first scalar is canonical S=sqrt(84)(s-sstar).
    H_a=(Hsource0(t),0,0,0); HsourceDot0 is its coordinate-time derivative.
    This wrapper checks array structure; per-node validity is arithmetic
    and metric admissibility, not an Einstein-constraint certificate.
    Numba parallelism respects the caller's configured thread count.
    """
    arrays=[np.ascontiguousarray(x,dtype=np.float64) for x in (g,v,dg,d2g,dv,q,p,dq,d2q,dp)]
    if arrays[0].ndim!=3 or arrays[0].shape[1:]!=(4,4):raise ValueError('g must have shape (N,4,4)')
    N=arrays[0].shape[0]
    tails=[(4,4),(4,4),(3,4,4),(3,3,4,4),(3,4,4),(2,),(2,),(3,2),(3,3,2),(3,2)]
    for value,tail in zip(arrays,tails):
        if value.shape!=(N,)+tail:raise ValueError('incompatible raw jet shape')
    parameters=np.array([r,Hsource0,HsourceDot0,sstar,kappa,gamma0],dtype=np.float64)
    if not np.isfinite(parameters).all():raise ValueError('finite scalar parameters required')
    if r<0 or kappa<0:raise ValueError('nonnegative inverse radius and gravitational coupling required')
    return _rhs_batch(*arrays,float(r),float(Hsource0),float(HsourceDot0),float(sstar),float(kappa),float(gamma0))
