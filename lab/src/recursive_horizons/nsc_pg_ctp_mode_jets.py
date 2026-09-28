"""Physical source-mode insertions for the transmitting PG CTP derivative.

The mode-space kernel acts on the whole spectral field with dE/(2pi).
Finite source-fiber matrices are quadrature controls, not closed packet states
or selected metric histories. Raw KS endpoint coordinates remain a separate
pullback of this PG history derivative.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag

from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import I2,S2,S3
from .nsc_transmitting_resolvent import profile,static_metric_kernel_from_fields
from .nsc_transmitting_ctp_variation import ctp_first_variation


PG_FIELDS=('log_N','beta','log_q_PG','log_r')


def continue_mode_columns(energies,mass,angular,at_zero,*,points=64):
    """Short spatial continuation of authenticated columns; no new preparation."""
    E=np.asarray(energies,float);initial=np.asarray(at_zero,complex)
    if E.ndim!=1 or initial.shape!=(len(E),2,3) or not np.isfinite(initial).all() or not np.isfinite(E).all():
        raise ValueError('physical signed-energy mode columns required')
    x,w=leggauss(points);nodes=np.r_[(x-1)/2,(x+1)/2];weights=np.r_[w/2,w/2]
    values=[];derivatives=[];seam=0.;closures=[]
    def generator(rho):
        beta,bp,_=geometry(rho)
        potential=CommonTimeBulkSplit.local_potential(1.,np.sqrt(1+rho*rho),mass,angular)
        iv=np.linalg.inv(S2-beta*I2)
        return np.einsum('ij,njk->nik',iv,1j*(E[:,None,None]*I2-potential)+.5*bp*I2)
    def rhs(rho,state):return np.einsum('nij,njk->nik',generator(rho),state.reshape(initial.shape)).ravel()
    for left,right in ((0.,-1.),(0.,1.)):
        run=solve_ivp(rhs,(left,right),initial.ravel(),method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True)
        if not run.success:raise ArithmeticError(run.message)
        closures.append(run.sol)
        seam=max(seam,float(np.max(abs(run.sol(0.)-initial.ravel()))))
    for rho in nodes:
        value=closures[int(rho>0)](rho).reshape(initial.shape)
        values.append(value);derivatives.append(np.einsum('nij,njk->nik',generator(rho),value))
    return {'rho':nodes,'weights':weights,'fields':np.asarray(values),
            'derivatives':np.asarray(derivatives),'seam_residual':seam}


def _columns(values):
    # (frequency,spin,source) -> (spin,frequency*source), no packet projection.
    return values.transpose(1,0,2).reshape(2,-1)


def apply_metric_expression(rho,mass,angular,columns,derivatives,parameters):
    """Apply the existing Z1 local Hamiltonian to fixed physical field columns."""
    s,sp=profile(rho);beta,bp,_=geometry(rho)
    N=np.exp(parameters[0]*s);q=np.exp(parameters[2]*s)
    radius=np.sqrt(1+rho*rho)*np.exp(parameters[3]*s)
    owner=CommonTimeBulkSplit()
    return np.column_stack([owner.apply_local_expression(u,up,N=N,q_PG=q,
       beta=float(beta)+parameters[1]*s,radius=radius,
       N_prime=N*parameters[0]*sp,q_prime=q*parameters[2]*sp,
       beta_prime=float(bp)+parameters[1]*sp,compact_mass=mass,angular_eigenvalue=angular)
       for u,up in zip(columns.T,derivatives.T)])


def mode_metric_vertices(fields,mass,angular,*,difference_step=2e-5):
    """M_A(Eo,Ei)=<Phi_Eo|delta H_A|Phi_Ei>, before energy integration.

    The owned weak-resolvent kernel is minus this Hamiltonian insertion.
    Its compact profile makes exterior endpoint terms vanish; the same seam
    trace is used from both sides. The direct control differentiates the
    existing Z1 Hamiltonian, holding the canonical mode columns fixed.
    """
    n=fields['fields'].shape[1]*3;weak=np.zeros((4,n,n),complex);direct=weak.copy()
    for rho,w,u,up in zip(fields['rho'],fields['weights'],fields['fields'],fields['derivatives'],strict=True):
        U,Up=_columns(u),_columns(up)
        weak-=w*profile(rho)[0]*static_metric_kernel_from_fields(rho,mass,angular,U,Up,U,Up)
        for A in range(4):
            plus=np.zeros(4);plus[A]=difference_step
            hp=apply_metric_expression(rho,mass,angular,U,Up,plus)
            hm=apply_metric_expression(rho,mass,angular,U,Up,-plus)
            direct[A]+=w*U.conj().T@(hp-hm)/(2*difference_step)
    return weak,direct


@dataclass(frozen=True)
class PGSourceModeHistoryDerivative:
    """Evaluated kernel of the full-field first-order Dyson insertion.

    In the source representation H_0 is multiplication by E. The derivative
    kernel is -i exp[-i Eo(tf-tau)] M_A(Eo,Ei) exp[-i Ei(tau-ti)]. Times are
    caller coordinates, not a selected physical duration or solution history.
    """
    energies: np.ndarray
    vertices: np.ndarray
    source_covariances: np.ndarray
    source_projectors: np.ndarray

    def __post_init__(self):
        E=np.asarray(self.energies,float);n=3*len(E)
        if (E.ndim!=1 or not np.isfinite(E).all() or self.vertices.shape!=(4,n,n)
            or self.source_covariances.shape!=(len(E),3,3) or self.source_projectors.shape!=(len(E),3,3)):
            raise ValueError('matching physical source energies, fibers and four-field vertices required')
        if not all(np.isfinite(a).all() for a in (self.vertices,self.source_covariances,self.source_projectors)):
            raise ValueError('all mode-space inputs must be evaluated')
        P=block_diag(*self.source_projectors)
        if np.linalg.norm(P@P-P)>3e-11:raise ValueError('canonical open-fiber projector required')
        if max(np.linalg.norm(P@v@P-v) for v in self.vertices)>3e-9:
            raise ValueError('a closed infinity channel was inserted into the physical vertex')

    def derivative_kernel(self,initial_time,final_time,insertion_time):
        times=np.asarray([initial_time,final_time,insertion_time],float)
        if not np.isfinite(times).all() or not initial_time<=insertion_time<=final_time:
            raise ValueError('ordered history coordinates required; this API selects none')
        energy=np.repeat(self.energies,3)
        return (-1j*np.exp(-1j*energy*(final_time-insertion_time))[None,:,None]*self.vertices
                *np.exp(-1j*energy*(insertion_time-initial_time))[None,None,:])

    def source_fiber_ctp_control(self,frequency_weights):
        """Finite Nyström control of the owned CTP differential, not full Gamma.

        The weight factors belong to the sampled operator, not C(E). Closed
        source columns are removed. No finite energy sample supplies a full
        source stress or licenses a finite closed-packet evolution.
        """
        w=np.asarray(frequency_weights,float)
        if w.shape!=self.energies.shape or not np.isfinite(w).all() or np.any(w<=0):
            raise ValueError('positive declared spectral quadrature weights required')
        P=block_diag(*self.source_projectors);active=np.flatnonzero(np.diag(P)>.5)
        C=block_diag(*self.source_covariances)[np.ix_(active,active)]
        if min(np.linalg.eigvalsh(C).min(),np.linalg.eigvalsh(np.eye(len(C))-C).min()) < -3e-11:
            raise ValueError('physical source covariance violates CAR')
        rootw=np.repeat(np.sqrt(w/(2*np.pi)),3)
        result=[];I=np.eye(len(active))
        for vertex in self.vertices:
            M=(rootw[:,None]*vertex*rootw[None,:])[np.ix_(active,active)]
            value=ctp_first_variation(C,I,I,-.5j*M,.5j*M,tolerance=3e-9)
            result.append({**value,'source_trace_residual':float(abs(value['derivative']+np.trace(C@M))),
                           'matrix_hermiticity':float(np.linalg.norm(M-M.conj().T)),
                           'open_source_rank':len(active)})
        return result

    def as_endpoint_branch_jets(self):
        raise ValueError('physical PG mode-history kernel still needs its raw KS endpoint and Cauchy-slice pullback')
