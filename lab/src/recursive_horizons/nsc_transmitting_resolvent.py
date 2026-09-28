"""Causal packet cross-resolvent of the inherited transmitted PG Dirac operator.

Solve the whole-line resolvent restricted to a strictly trapped probe cell.
Zero inflow is fixed by causal support, not a reflecting finite-box domain.
Project only the solved continuous field. No instantaneous B is supplied.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_lorentzian import geometry
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_transmitting_dirac_domain import I2, S2, S3


FIELDS = ("log_N", "beta", "log_q_PG", "log_r")


def profile(rho):
    """One smooth compact metric-direction probe, normalized to s(0)=1."""
    if abs(rho) >= 1.: return 0., 0.
    d = 1-rho*rho
    s = np.exp(1-1/d)
    return float(s), float(-2*rho*s/d**2)


def probes(rho):
    """Two unit-L2 packet profiles, with disjoint parent/child supports."""
    fp = np.sqrt(30.)*rho*(1-rho) if 0 < rho < 1 else 0.
    fc = np.sqrt(30.)*(-rho)*(1+rho) if -1 < rho < 0 else 0.
    return np.hstack((fp*I2, fc*I2))


def quadrature(points=96):
    x, w = leggauss(points)
    return np.r_[0.5*(x-1), 0.5*(x+1)], np.r_[w/2, w/2]


def static_metric_kernel_from_fields(rho, compact_mass, angular_eigenvalue, u, up, l, lp):
    """Recover the static four-field kernel from archived quadrature fields."""
    radius=np.sqrt(1+rho*rho)
    m=CommonTimeBulkSplit.local_potential(1.,radius,compact_mass,angular_eigenvalue)
    zero=np.zeros((2,2),complex)
    vfields=(S2,-I2,-S2,zero)
    mfields=(m,zero,zero,angular_eigenvalue/radius*S3)
    return np.array([0.5j*(l.conj().T@v@up-lp.conj().T@v@u)-l.conj().T@dm@u for v,dm in zip(vfields,mfields)])


@dataclass
class PacketResolvent:
    owner: object
    z: complex
    parameters: np.ndarray
    segments: dict
    response: np.ndarray
    directional_jets: np.ndarray
    seam_residual: float
    function_evaluations: int

    def field(self, rho):
        return self.segments["parent" if rho >= 0 else "child"].sol(rho)[:8].reshape(2,4)

    def field_derivative(self, rho):
        a, source, _, _ = self.owner.coefficients(rho, self.z, self.parameters)
        return a@self.field(rho)+source

    def as_endpoint_branch_jets(self):
        raise ValueError("a static resolvent jet is not dU on a CTP history; spectral/time reconstruction and the physical state/history are required")


@dataclass(frozen=True)
class TransmittingCrossResolvent:
    compact_mass: float
    angular_eigenvalue: float
    rtol: float = 2e-12
    atol: float = 2e-14

    def __post_init__(self):
        if not np.isfinite([self.compact_mass,self.angular_eigenvalue,self.rtol,self.atol]).all() or min(self.compact_mass,self.angular_eigenvalue) < 0 or min(self.rtol,self.atol) <= 0:
            raise ValueError("finite retained labels and positive solver tolerances required")

    def metric_parts(self, rho, parameters):
        s, sp = profile(rho)
        beta0, bp0, _ = geometry(rho)
        N = np.exp(parameters[0]*s);q = np.exp(parameters[2]*s)
        beta = float(beta0)+parameters[1]*s
        radius = np.sqrt(1+rho*rho)*np.exp(parameters[3]*s)
        a = N/q;ap = a*(parameters[0]-parameters[2])*sp
        v = CommonTimeBulkSplit.principal(N,q,beta)
        vp = ap*S2-(float(bp0)+parameters[1]*sp)*I2
        m = CommonTimeBulkSplit.local_potential(N,radius,self.compact_mass,self.angular_eigenvalue)
        zero = np.zeros((2,2),complex)
        dv = np.array([s*a*S2, -s*I2, -s*a*S2, zero])
        dvp = np.array([(sp*a+s*ap)*S2, -sp*I2, -(sp*a+s*ap)*S2, zero])
        dm = np.array([s*m, zero, zero, N*self.angular_eigenvalue*s/radius*S3])
        return v, vp, m, dv, dvp, dm

    def coefficients(self, rho, z, parameters):
        v, vp, m, dv, dvp, dm = self.metric_parts(rho, parameters)
        iv = np.linalg.inv(v)
        a = iv@(1j*(z*I2-m)-0.5*vp)
        source = 1j*iv@probes(rho)
        da = np.array([-iv@x@a+iv@(-1j*y-0.5*xp) for x,xp,y in zip(dv,dvp,dm)])
        ds = np.array([-iv@x@source for x in dv])
        return a, source, da, ds

    def solve(self, z, *, parameters=None, with_jets=True):
        z = complex(z)
        if not np.isfinite([z.real,z.imag]).all() or z.imag == 0:
            raise ValueError("nonreal resolvent energy required")
        parameters = np.zeros(4) if parameters is None else np.asarray(parameters,float)
        if parameters.shape != (4,) or not np.isfinite(parameters).all():raise ValueError("four finite compact metric-probe amplitudes required")
        # This uniform bound preserves both negative characteristic speeds
        # throughout [-1,1], including all derivative-check perturbations.
        if float(geometry(1.)[0])-abs(parameters[1]) <= np.exp(abs(parameters[0]-parameters[2])):
            raise ValueError("probe deformation leaves the strictly trapped cell")
        count = 5 if with_jets else 1
        orientation = -1. if z.imag > 0 else 1.
        order = ((1.,0.,"parent"),(0.,-1.,"child")) if z.imag > 0 else ((-1.,0.,"child"),(0.,1.,"parent"))
        state = np.zeros(count*24,complex)  # fields: count*2x4, packet integrals: count*4x4
        def rhs(rho,y):
            u = y[:count*8].reshape(count,2,4)
            a,source,da,ds = self.coefficients(rho,z,parameters)
            du = np.einsum('ij,kjl->kil',a,u)
            du[0] += source
            if with_jets:
                du[1:] += np.einsum('kij,jl->kil',da,u[0])+ds
            integrals = orientation*np.einsum('ij,kjl->kil',probes(rho).conj().T,u)
            return np.r_[du.ravel(),integrals.ravel()]
        segments = {}; seam = 0.;nfev = 0
        for start,end,side in order:
            initial = state.copy()
            sol = solve_ivp(rhs,(start,end),initial,method='DOP853',rtol=self.rtol,atol=self.atol,dense_output=True)
            if not sol.success:raise ArithmeticError(sol.message)
            if start == 0.:seam = float(np.max(abs(sol.y[:,0]-state)))
            segments[side] = sol;state=sol.y[:,-1];nfev += sol.nfev
        packet = state[count*8:].reshape(count,4,4)
        return PacketResolvent(self,z,parameters,segments,packet[0],packet[1:] if with_jets else np.empty((0,4,4)),seam,nfev)

    def static_metric_kernel(self, upper, lower, rho):
        """dR_packet = sum_A integral s_A(rho) J_A(rho) d rho.

        Compactly supported log-N/log-q/log-r and additive-shift variations;
        the profile derivative is integrated by parts. Normalized canonical
        probe functions and the coordinate cut are held fixed.
        """
        if abs(lower.z-upper.z.conjugate()) > 1e-14 or np.any(upper.parameters) or np.any(lower.parameters):
            raise ValueError("unperturbed conjugate resolvent solutions required")
        u, up = upper.field(rho), upper.field_derivative(rho)
        l, lp = lower.field(rho), lower.field_derivative(rho)
        return static_metric_kernel_from_fields(rho,self.compact_mass,self.angular_eigenvalue,u,up,l,lp)


def inverse_response(result):
    kernel=np.linalg.solve(result.response,np.eye(len(result.response)))
    jets=np.array([-kernel@dr@kernel for dr in result.directional_jets])
    return kernel,jets
