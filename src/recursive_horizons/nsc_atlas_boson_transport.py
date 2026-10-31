"""Finite canonical transport for inherited bosonic matter covariances.

The canonical ordering is z=(q_1,...,q_n,p_1,...,p_n), with [z_a,z_b]
=i*hbar*J_ab.  The Hamiltonian is H=.5*z.T*G(t)*z.  No initial quantum
state is supplied by this module: a covariance must come from the caller.

CF4 is the fourth-order, two-exponential commutator-free Magnus formula.
Each exponential has a Hamiltonian generator when the supplied G is
symmetric.  Roundoff symplectic defects are measured, never projected away.
The covariance is the symmetrized central second moment, not a fermionic
occupation matrix.  In particular, it has no upper bound C<=I.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import cholesky, eigvalsh, expm


def canonical_form(modes):
    """Return J in q-then-p ordering for a positive number of modes."""
    if isinstance(modes,(bool,np.bool_)) or not isinstance(modes,(int,np.integer)) or modes<1:
        raise ValueError('a positive integer number of canonical modes is required')
    eye=np.eye(modes)
    return np.block([[np.zeros_like(eye),eye],[-eye,np.zeros_like(eye)]])


def _matrix(value,name,*,symmetric=False):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.ndim!=2 or raw.shape[0]!=raw.shape[1] or raw.shape[0]<2 or raw.shape[0]%2:
        raise ValueError(name+' requires a real square matrix of positive even size')
    result=np.asarray(raw,dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(name+' must be finite')
    if symmetric:
        scale=max(1.,float(np.max(abs(result))))
        if np.max(abs(result-result.T))>5e-12*scale:
            raise ValueError(name+' must be symmetric')
    return result


def hamiltonian_generator(hamiltonian):
    """Return A=JG without changing the supplied symmetric Hamiltonian."""
    G=_matrix(hamiltonian,'Hamiltonian G',symmetric=True)
    return canonical_form(len(G)//2)@G


def symplectic_diagnostics(matrix):
    """Report the actual canonical-form defect; no group projection is used."""
    M=_matrix(matrix,'transport M');J=canonical_form(len(M)//2)
    defect=M.T@J@M-J
    sign,logdet=np.linalg.slogdet(M)
    return dict(max_abs=float(np.max(abs(defect))),frobenius=float(np.linalg.norm(defect)),
                determinant_sign=float(sign),log_abs_determinant=float(logdet))


def cf4_step(hamiltonian,t,h,matrix=None):
    """Advance M'=JG(t)M by one fourth-order commutator-free step.

    The right exponential weights the early node most strongly and acts
    first. Reversing these exponentials changes the Magnus commutator sign.
    Negative h is allowed; h=0 and nonfinite intervals are rejected.
    """
    t=float(t);h=float(h)
    if not math.isfinite(t) or not math.isfinite(h) or h==0:
        raise ValueError('a finite nonzero integration step is required')
    root=math.sqrt(3.)
    early=hamiltonian_generator(hamiltonian(t+h*(.5-root/6)))
    late=hamiltonian_generator(hamiltonian(t+h*(.5+root/6)))
    if early.shape!=late.shape:raise ValueError('Hamiltonian dimension changed during a step')
    M=np.eye(len(early)) if matrix is None else _matrix(matrix,'transport M')
    if M.shape!=early.shape:raise ValueError('Hamiltonian and transport dimensions differ')
    a1=(3-2*root)/12;a2=(3+2*root)/12
    result=expm(h*(a1*early+a2*late))@(expm(h*(a2*early+a1*late))@M)
    if not np.all(np.isfinite(result)):raise FloatingPointError('nonfinite canonical transport')
    return result


@dataclass(frozen=True)
class TransportResult:
    """Maps from the initial event to requested checkpoints and the endpoint."""
    times: np.ndarray
    matrices: np.ndarray
    steps: int
    hamiltonian_evaluations: int
    max_checkpoint_symplectic_defect: float

    @property
    def matrix(self):
        return self.matrices[-1]


def propagate(hamiltonian,breakpoints,*,max_step,checkpoints=(),initial=None):
    """Transport across a monotone piecewise-smooth background history.

    Every supplied breakpoint and checkpoint is an exact step boundary.
    Only requested checkpoints plus the two endpoints are stored. Forward
    and backward histories are supported. An initial map may be supplied
    for composition; otherwise the initial map is identity. The callback
    is evaluated only inside the corresponding background interval.
    """
    edges=np.asarray(breakpoints,dtype=float)
    requested=np.asarray(checkpoints,dtype=float)
    max_step=float(max_step)
    if edges.ndim!=1 or len(edges)<2 or not np.all(np.isfinite(edges)):
        raise ValueError('at least two finite background breakpoints are required')
    direction=1. if edges[-1]>edges[0] else -1.
    if np.any(direction*np.diff(edges)<=0):raise ValueError('background breakpoints must be strictly monotone')
    if not math.isfinite(max_step) or max_step<=0:raise ValueError('max_step must be finite and positive')
    if requested.ndim!=1 or not np.all(np.isfinite(requested)):
        raise ValueError('checkpoints must be a finite one-dimensional sequence')
    lower,upper=sorted((edges[0],edges[-1]))
    if np.any(requested<lower) or np.any(requested>upper):raise ValueError('checkpoint outside the background history')
    stored=set(map(float,requested));stored.update((float(edges[0]),float(edges[-1])))
    partition=sorted(set(map(float,edges))|stored,reverse=direction<0)
    calls=0
    def counted(t):
        nonlocal calls
        calls+=1
        return hamiltonian(t)
    if initial is None:
        G=_matrix(counted(float(edges[0])),'Hamiltonian G',symmetric=True)
        M=np.eye(len(G))
    else:M=_matrix(initial,'initial transport M').copy()
    times=[float(edges[0])];maps=[M.copy()];steps=0
    for left,right in zip(partition,partition[1:]):
        count=max(1,math.ceil(abs(right-left)/max_step))
        h=(right-left)/count
        for k in range(count):
            M=cf4_step(counted,left+k*h,h,M);steps+=1
        if right in stored:times.append(right);maps.append(M.copy())
    return TransportResult(np.asarray(times),np.asarray(maps),steps,calls,
        max(symplectic_diagnostics(M)['max_abs'] for M in maps))


def covariance_diagnostics(covariance,*,hbar=1.):
    """Check positivity and Gamma+i*hbar*J/2, with no covariance repair."""
    G=_matrix(covariance,'bosonic covariance Gamma',symmetric=True)
    hbar=float(hbar)
    if not math.isfinite(hbar) or hbar<=0:raise ValueError('hbar must be finite and positive')
    J=canonical_form(len(G)//2)
    least=float(eigvalsh(G,check_finite=False)[0])
    uncertainty=float(eigvalsh(G+.5j*hbar*J,check_finite=False)[0])
    scale=max(hbar,1.,float(np.max(abs(G))))
    tolerance=5e-11*scale
    minimum_nu=float(symplectic_eigenvalues(G)[0]) if least>0 else None
    return dict(min_covariance_eigenvalue=least,min_CCR_eigenvalue=uncertainty,
                min_symplectic_eigenvalue=minimum_nu,CCR_tolerance=tolerance,
                symmetry_max_abs=float(np.max(abs(G-G.T))),
                physical=least>0 and uncertainty>=-tolerance
                    and minimum_nu>=hbar*(.5-5e-11))


def validate_covariance(covariance,*,hbar=1.):
    """Return the supplied covariance if its finite CCR test passes."""
    G=_matrix(covariance,'bosonic covariance Gamma',symmetric=True)
    if not covariance_diagnostics(G,hbar=hbar)['physical']:
        raise ValueError('covariance violates positive bosonic CCR uncertainty')
    return G


def covariance_inheritance(matrix,covariance,*,hbar=1.,validate=True):
    """Inherit Gamma_out=M Gamma_parent M.T, preserving cross correlations.

    Validation checks input and output physicality; it does not pick a
    covariance or replace either covariance with a nearby valid matrix.
    """
    M=_matrix(matrix,'transport M')
    G=validate_covariance(covariance,hbar=hbar) if validate else _matrix(covariance,'bosonic covariance Gamma',symmetric=True)
    if M.shape!=G.shape:raise ValueError('transport and covariance dimensions differ')
    inherited=M@G@M.T
    if validate:validate_covariance(inherited,hbar=hbar)
    return inherited


def symplectic_eigenvalues(covariance):
    """Return raw positive symplectic eigenvalues of a positive covariance."""
    G=_matrix(covariance,'bosonic covariance Gamma',symmetric=True)
    try:L=cholesky(G,lower=True,check_finite=False)
    except np.linalg.LinAlgError as error:raise ValueError('covariance must be positive definite') from error
    spectrum=eigvalsh(1j*(L.T@canonical_form(len(G)//2)@L),check_finite=False)
    result=spectrum[len(G)//2:]
    if not np.all(np.isfinite(result)) or np.any(result<=0):
        raise ValueError('finite positive symplectic eigenvalues were not resolved')
    return result


def gaussian_entropy(covariance,*,hbar=1.):
    """Von Neumann entropy in nats for a supplied bosonic Gaussian state.

    The analytic zero-occupation limit is used within roundoff of nu=hbar/2.
    Raw symplectic eigenvalues and the original covariance remain unchanged.
    This is not an entropy production or equilibrium-selection prescription.
    """
    G=validate_covariance(covariance,hbar=hbar)
    occupation=symplectic_eigenvalues(G)/float(hbar)-.5
    tolerance=5e-11
    if np.min(occupation)<-tolerance:raise ValueError('symplectic spectrum violates quantum uncertainty')
    positive=occupation>0
    n=occupation[positive]
    return float(np.sum(np.log1p(n)+n*np.log1p(1/n)))


def reduced_covariance(covariance,modes):
    """Retain selected canonical pairs, in q-then-p ordering."""
    G=_matrix(covariance,'bosonic covariance Gamma',symmetric=True);n=len(G)//2
    raw=np.asarray(modes)
    if raw.ndim!=1 or raw.dtype.kind not in 'iu' or len(raw)==0:
        raise ValueError('a nonempty integer mode selection is required')
    if len(set(map(int,raw)))!=len(raw) or np.any(raw<0) or np.any(raw>=n):
        raise ValueError('mode selection must be distinct and within the covariance')
    indices=np.r_[raw,raw+n]
    return G[np.ix_(indices,indices)]


def retarded_qp(matrix):
    """Return the q response to a unit initial canonical-momentum kick.

    For a forward propagator this equals (i/hbar)[q(t),q(t0)] and is
    independent of the inherited state. A time-ordering step function is
    not inserted because the caller owns the two events.
    """
    M=_matrix(matrix,'transport M');n=len(M)//2
    return M[:n,n:]


def relative_transport(later,earlier):
    """Return M(t,s)=M(t,t0)*M(s,t0)^-1 by an actual linear solve.

    The symplectic transpose identity is not substituted for the inverse:
    doing so could hide an accumulated numerical group defect.
    """
    late=_matrix(later,'later transport M');early=_matrix(earlier,'earlier transport M')
    if late.shape!=early.shape:raise ValueError('the two transport dimensions differ')
    try:return np.linalg.solve(early.T,late.T).T
    except np.linalg.LinAlgError as error:raise ValueError('earlier transport is singular') from error


def retarded_between(later,earlier):
    """Return the canonical momentum-kick response between two saved events."""
    return retarded_qp(relative_transport(later,earlier))
