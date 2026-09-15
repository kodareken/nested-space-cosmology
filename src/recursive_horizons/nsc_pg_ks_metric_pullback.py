"""Raw KS metric directions in the owned fixed-reference PG chart.

This is a metric-coordinate pullback, not a Cauchy-state isometry.  The
reference chart is fixed; history endpoints and a moving surface are not
selected by the map.
"""
import numpy as np

from .nsc_lorentzian import geometry
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_transmitting_resolvent import profile,static_metric_kernel_from_fields
from .nsc_pg_ctp_mode_jets import _columns
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain


def reference_chart(rho):
    beta,bp,_=map(float,geometry(rho))
    if beta<=1:raise ValueError('the owned KS chart is the strictly trapped region')
    a=np.sqrt(beta*beta-1);ap=beta*bp/a
    # (dT,dz) = coordinate_jacobian (d tau,d rho).
    J=np.array([[0.,-1/a],[1.,beta/(a*a)]])
    return beta,bp,a,ap,J


def ks_to_pg(rho,NK,betaK,aK,rK):
    b,_,a,_,_=reference_chart(rho)
    if not np.isfinite([NK,betaK,aK,rK]).all() or min(NK,aK,rK)<=0:
        raise ValueError('finite positive raw KS metric fields required')
    s=b/(a*a)-betaK/a
    q2=aK*aK*s*s-NK*NK/(a*a)
    if q2<=0 or s<=0:raise ValueError('deformation leaves the connected PG Cauchy chart')
    q=np.sqrt(q2);N=NK*aK/(a*q);beta=aK*aK*s/q2
    return np.array([N,beta,q,rK])


def pg_to_ks(rho,N,beta,q,r):
    b,_,a,_,_=reference_chart(rho)
    if min(N,q,r)<=0 or q*beta<=N:raise ValueError('connected trapped PG metric required')
    axial=np.sqrt(q*q*beta*beta-N*N)
    return np.array([a*N*q/axial,b/a-a*q*q*beta/(axial*axial),axial,r])


def pg_log_jacobian(rho,NK,betaK,aK,rK):
    """d(log N_PG,beta_PG,log q_PG,log r)/d(N_KS,beta_KS,a_KS,r)."""
    b,_,a,_,_=reference_chart(rho);N,beta,q,_=ks_to_pg(rho,NK,betaK,aK,rK)
    s=b/(a*a)-betaK/a;F=q*q
    logq=np.array([-NK/(a*a),-aK*aK*s/a,aK*s*s,0.])/F
    logN=np.array([1/NK,0.,1/aK,0.])-logq
    dbeta=np.array([0.,-aK*aK/a,2*aK*s,0.])/F-2*beta*logq
    return np.array([logN,dbeta,logq,[0.,0.,0.,1/rK]])


def deformed_pg_fields(rho,amplitudes):
    """Same compact direction expressed in raw KS fields, including rho jets."""
    amplitude=np.asarray(amplitudes,float)
    if amplitude.shape!=(4,) or not np.isfinite(amplitude).all():raise ValueError('four raw KS amplitudes required')
    b,bp,a,ap,_=reference_chart(rho);shape,shapep=profile(rho);radius=np.sqrt(1+rho*rho)
    NK,bK,aK,rK=np.array([1.,0.,a,radius])+amplitude*shape
    NKp,bKp,aKp,_=np.array([0.,0.,ap,rho/radius])+amplitude*shapep
    N,beta,q,r=ks_to_pg(rho,NK,bK,aK,rK)
    s=b/(a*a)-bK/a
    sp=bp/(a*a)-2*b*ap/(a*a*a)-bKp/a+bK*ap/(a*a)
    F=q*q
    Fp=2*aK*aKp*s*s+2*aK*aK*s*sp-2*NK*NKp/(a*a)+2*NK*NK*ap/(a*a*a)
    qp=Fp/(2*q)
    betap=(2*aK*aKp*s+aK*aK*sp)/F-beta*Fp/F
    Np=N*(NKp/NK+aKp/aK-ap/a-qp/q)
    return N,beta,q,r,Np,betap,qp


def pulled_mode_vertices(fields,mass,angular,*,difference_step=2e-5):
    """Integrate the varying chart Jacobian, not a constant neck matrix."""
    n=fields['fields'].shape[1]*3
    weak=np.zeros((4,n,n),complex);direct=weak.copy();owner=CommonTimeBulkSplit()
    for rho,w,u,up in zip(fields['rho'],fields['weights'],fields['fields'],fields['derivatives'],strict=True):
        U,Up=_columns(u),_columns(up);_,_,axial,_,_=reference_chart(rho)
        J=pg_log_jacobian(rho,1.,0.,axial,np.sqrt(1+rho*rho))
        kernel=static_metric_kernel_from_fields(rho,mass,angular,U,Up,U,Up)
        weak-=w*profile(rho)[0]*np.einsum('AB,Aij->Bij',J,kernel)
        for B in range(4):
            a=np.zeros(4);a[B]=difference_step;values=[]
            for sign in (1,-1):
                N,beta,q,r,Np,bp,qp=deformed_pg_fields(rho,sign*a)
                values.append(np.column_stack([owner.apply_local_expression(v,vp,N=N,q_PG=q,beta=beta,radius=r,
                    N_prime=Np,q_prime=qp,beta_prime=bp,compact_mass=mass,angular_eigenvalue=angular)
                    for v,vp in zip(U.T,Up.T)]))
            direct[B]+=w*U.conj().T@(values[0]-values[1])/(2*difference_step)
    return weak,direct


def coordinate_residuals(rho):
    """Same metric in both charts, plus independently differenced field map."""
    b,_,a,_,J=reference_chart(rho);r=np.sqrt(1+rho*rho)
    # A nearby metric is a coordinate-identity control, not a source solution.
    K=np.array([1.002,.001,a*1.003,r*1.001]);N,beta,q,rp=ks_to_pg(rho,*K)
    NK,bK,aK,rK=K
    gK=np.array([[NK*NK-aK*aK*bK*bK,-aK*aK*bK],[-aK*aK*bK,-aK*aK]])
    gP=np.array([[N*N-q*q*beta*beta,-q*q*beta],[-q*q*beta,-q*q]])
    analytic=pg_log_jacobian(rho,*K);numerical=np.empty((4,4));h=1e-5
    def log_fields(x):
        N,beta,q,r=ks_to_pg(rho,*x);return np.array([np.log(N),beta,np.log(q),np.log(r)])
    for B in range(4):
        e=np.zeros(4);e[B]=h
        numerical[:,B]=(-log_fields(K+2*e)+8*log_fields(K+e)-8*log_fields(K-e)+log_fields(K-2*e))/(12*h)
    return {'metric_tensor_reconstruction':float(np.max(abs(gP-J.T@gK@J))),
            'inverse_field_map':float(np.max(abs(pg_to_ks(rho,N,beta,q,rp)-K))),
            'jacobian_vs_finite_difference':float(np.max(abs(analytic-numerical)))}


def fixed_surface_normal(rho,NK,betaK,aK,rK):
    """Normal and its raw KS derivative on the fixed surface, not an isometry."""
    *_,J=reference_chart(rho);inverse=np.linalg.inv(J)
    normal=inverse@np.array([1/NK,-betaK/NK])
    derivative=inverse@np.array([[-1/NK**2,0.,0.,0.],[betaK/NK**2,-1/NK,0.,0.]])
    N,beta,q,r=ks_to_pg(rho,NK,betaK,aK,rK)
    direct=TransmittingDiracSeamDomain(N,q,beta,r).normal_geometry()
    return {'normal':normal,'raw_KS_derivative':derivative,
            'normal_geometry_residual':float(np.max(abs(normal-direct['future_normal']))),
            'normal_unit_residual':direct['normal_unit_residual'],
            'normal_tangent_residual':direct['normal_tangent_residual']}
