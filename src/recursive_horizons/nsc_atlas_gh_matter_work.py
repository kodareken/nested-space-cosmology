"""First-jet normal-energy balance for the full canonical two-scalar stress.

The reference S3 frame is divergence free and obeys
[E_i,E_j]=2*r*epsilon_ijk E_k.  A static reference metric is used; no GH
acceleration or field evolution is evaluated by this diagnostic.
"""
import math

import numpy as np
from numba import njit,prange

from .nsc_atlas_gh_pointwise import _adm_inverse,_reference_connection,_slot_action


def _array(value,name,shape):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape!=shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+' needs finite real values with shape '+str(shape))
    return np.ascontiguousarray(raw,dtype=float)


@njit(cache=True,nogil=True,parallel=True)
def _work_batch(g,v,dg,p,dq,W,r):
    count=len(g);connection=_reference_connection(r)
    rho=np.full(count,np.nan);work=np.full(count,np.nan);measure=np.full(count,np.nan)
    flux=np.full((count,3),np.nan)
    lapse=np.full(count,np.nan);extrinsic_trace=np.full(count,np.nan)
    for point in prange(count):
        gi,normal,valid=_adm_inverse(g[point])
        if not valid:continue
        N=1/normal[0];beta=-normal[1:]*N
        hi=np.empty((3,3))
        for i in range(3):
            for j in range(3):hi[i,j]=gi[i+1,j+1]+normal[i+1]*normal[j+1]
        D=np.empty((4,4,4));D[0]=v[point]
        for i in range(3):D[i+1]=dg[point,i]-_slot_action(g[point],connection[i])
        K=np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                for a in range(4):
                    K[i,j]-=N*gi[0,a]*(D[i+1,a,j+1]+D[j+1,a,i+1]-D[a,i+1,j+1])/2
        dN=np.zeros(3)
        for i in range(3):
            for a in range(4):
                for b in range(4):dN[i]-=N**3*gi[0,a]*dg[point,i,a,b]*gi[b,0]/2
        pi=np.empty(2);pi2=0.;gradient2=0.
        for A in range(2):
            value=p[point,A]
            for i in range(3):value-=beta[i]*dq[point,i,A]
            pi[A]=value/N;pi2+=pi[A]*pi[A]
            for i in range(3):
                for j in range(3):gradient2+=hi[i,j]*dq[point,i,A]*dq[point,j,A]
        density=(pi2+gradient2)/2+W[point]
        pressure_isotropic=(pi2-gradient2)/2-W[point]
        jcov=np.zeros(3);jup=np.zeros(3);Kupper=np.zeros((3,3));traceK=0.
        for i in range(3):
            for A in range(2):jcov[i]-=pi[A]*dq[point,i,A]
            for j in range(3):
                traceK+=hi[i,j]*K[i,j]
                for k in range(3):
                    for ell in range(3):Kupper[i,j]+=hi[i,k]*hi[j,ell]*K[k,ell]
        stress_work=traceK*pressure_isotropic
        for i in range(3):
            for j in range(3):
                jup[i]+=hi[i,j]*jcov[j]
                for A in range(2):stress_work+=Kupper[i,j]*dq[point,i,A]*dq[point,j,A]
        value=N*stress_work
        for i in range(3):value-=jup[i]*dN[i]
        h=g[point,1:,1:]
        det=(h[0,0]*(h[1,1]*h[2,2]-h[1,2]*h[1,2])
             -h[0,1]*(h[0,1]*h[2,2]-h[0,2]*h[1,2])
             +h[0,2]*(h[0,1]*h[1,2]-h[0,2]*h[1,1]))
        weight=math.sqrt(det)
        rho[point]=density;work[point]=value;measure[point]=weight
        lapse[point]=N;extrinsic_trace[point]=traceK
        for i in range(3):flux[point,i]=weight*(N*jup[i]-beta[i]*density)
    return rho,work,measure,flux,lapse,extrinsic_trace


def matter_work(g,v,dg,q,p,dq,*,r,potential):
    """Return normal density, work and flux from raw spacetime first jets.

    Trailing shapes: g,v=(4,4), dg=(3,4,4), q,p=(2,), dq=(3,2), with common
    optional leading batch axes.  Fields are canonical (S,chi), p=partial_t q
    and dq_i=E_i q; r is already the physical inverse reference radius.
    The full supplied potential W includes its vacuum contribution once.

    ``work_integrand`` is N*K_ij*S^ij-j^i*E_iN BEFORE multiplication by
    ``volume_weight``=sqrt(det h).  ``energy_flux`` is already densitized:
    sqrt(det h)*(N*j^i-beta^i*rho).  With conserved full two-scalar stress,

      partial_t(sqrt(h)*rho)+E_i energy_flux^i
          =sqrt(h)*work_integrand.

    Consequently the compact S3 integral has no boundary-flux or separate
    shift term.  The reference volume measure is supplied by the caller;
    for unit-sphere quadrature it includes the common r^-3 factor.  These
    are slice-normal energies and geometric work, not a conserved geodesic
    energy or a local cosmological-constant measurement.  No subtraction
    of the vacuum or independently coupled field sector is made.

    ``lapse`` and ``extrinsic_trace`` permit an explicit constant-vacuum
    subtraction: for fixed V_star, excitation density is rho-V_star and
    its compact integrated work is work_integrand+N*V_star*tr(K).
    This subtracts the vacuum's own work -N*V_star*tr(K), without changing
    sources, the metric evolution or the full stress balance above.
    """
    raw=np.asarray(g)
    if raw.ndim<2 or raw.shape[-2:]!=(4,4):raise ValueError('g needs trailing shape (4,4)')
    batch=raw.shape[:-2]
    if raw.size==0:raise ValueError('at least one point is required')
    g=_array(g,'g',batch+(4,4));v=_array(v,'v',batch+(4,4))
    dg=_array(dg,'raw frame metric derivatives dg',batch+(3,4,4))
    _array(q,'canonical fields q',batch+(2,))
    p=_array(p,'coordinate rates p',batch+(2,));dq=_array(dq,'raw frame field derivatives dq',batch+(3,2))
    for value,name in ((g,'g'),(v,'v'),(dg,'dg')):
        tolerance=2e-12*max(1.,float(np.max(abs(value))))
        if np.max(abs(value-np.swapaxes(value,-1,-2)))>tolerance:
            raise ValueError(name+' must be symmetric in metric indices')
    if isinstance(r,(bool,np.bool_)):raise ValueError('r must be finite and nonnegative')
    r=float(r)
    if not math.isfinite(r) or r<0:raise ValueError('r must be finite and nonnegative')
    W=np.asarray(potential)
    if W.dtype.kind not in 'iuf' or not np.all(np.isfinite(W)):
        raise ValueError('potential must be finite and real')
    try:W=np.ascontiguousarray(np.broadcast_to(W,batch),dtype=float)
    except ValueError as error:raise ValueError('potential must broadcast to the batch shape') from error
    count=int(np.prod(batch)) if batch else 1
    result=_work_batch(g.reshape(count,4,4),v.reshape(count,4,4),dg.reshape(count,3,4,4),
                       p.reshape(count,2),dq.reshape(count,3,2),W.reshape(count),r)
    if any(not np.all(np.isfinite(value)) for value in result):
        raise ValueError('metric/lapse invalid or derived matter work is not finite')
    return {name:value.reshape(batch+tail) for name,value,tail in
            zip(('rho','work_integrand','volume_weight','energy_flux','lapse','extrinsic_trace'),
                result,((),(),(),(3,),(),()))}
