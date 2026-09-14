"""CTP retarded memory and frequency-transfer vertices of the owned resolvent.

The packet inverse is a nonlocal quadratic kernel, never a full branch unitary.
No state covariance is inferred from its retarded/advanced components.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_transmitting_resolvent import (
    TransmittingCrossResolvent, probes, static_metric_kernel_from_fields,
)


@dataclass(frozen=True)
class CTPRetardedMemory:
    """Convention conversion from (H-z)^-1 to -i theta(t-s) U(t,s).

    The full 4x4 packet response is retained; its cross block is not promoted
    to an independent Hermitian link. Keldysh/state data remain a separate
    two-time object, including the preparation on the original Cauchy slice.
    """

    response: np.ndarray
    advanced: np.ndarray

    def __post_init__(self):
        r, a = (np.array(x, dtype=complex, copy=True) for x in (self.response, self.advanced))
        if r.shape != (4, 4) or a.shape != r.shape or not np.isfinite([r, a]).all():
            raise ValueError("full finite parent/child packet response pair required")
        if np.linalg.norm(a-r.conj().T) > 3e-11*max(1., np.linalg.norm(r)):
            raise ValueError("adjoint advanced response required")
        r.setflags(write=False); a.setflags(write=False)
        object.__setattr__(self, 'response', r)
        object.__setattr__(self, 'advanced', a)

    @property
    def G_retarded(self):
        return -self.response

    @property
    def G_advanced(self):
        return -self.advanced

    @property
    def D_retarded(self):
        return -np.linalg.solve(self.response, np.eye(4))

    @property
    def D_advanced(self):
        return -np.linalg.solve(self.advanced, np.eye(4))

    def relative_vertex(self, input_memory, delta_response):
        """dD_R(zo,zi)=K(zo) dR(zo,zi) K(zi); no constant-B reduction."""
        d = np.asarray(delta_response, dtype=complex)
        if d.shape[-2:] != (4, 4) or not np.isfinite(d).all():
            raise ValueError("finite packet response derivatives required")
        return np.array([self.D_retarded @ x @ input_memory.D_retarded
                         for x in d.reshape(-1, 4, 4)]).reshape(d.shape)

    def as_endpoint_branch_jets(self):
        raise ValueError("retarded memory vertices lack full PG state/preparation, time-history reconstruction and KS endpoint pullback")


class ArchivedResolventField:
    """Interpolate authenticated Gauss-Legendre fields on each side separately.

    Barycentric weights use the known Gauss rule, without a new radial solve.
    The differential equation, not a differentiated interpolation polynomial,
    supplies the derivative when evaluating a metric vertex.
    """

    def __init__(self, owner, z, nodes, weights, values):
        self.owner, self.z = owner, complex(z)
        self.nodes = np.asarray(nodes)
        self.values = np.asarray(values)
        n = len(self.nodes)//2
        if len(self.nodes) != 2*n or self.values.shape != (2*n, 2, 4):
            raise ValueError("two-sided archived 2x4 resolvent fields required")
        self.n = n
        self.barycentric = []
        for x, w in ((self.nodes[:n], weights[:n]), (self.nodes[n:], weights[n:])):
            canonical = 2*x + (1 if x[-1] < 0 else -1)
            self.barycentric.append((-1.)**np.arange(n)*np.sqrt((1-canonical**2)*w))

    def side_value(self, rho, side):
        sl = slice(side*self.n, (side+1)*self.n)
        distances = rho-self.nodes[sl]
        exact = np.flatnonzero(distances == 0.)
        if len(exact): return self.values[sl][exact[0]]
        fractions = self.barycentric[side]/distances
        return np.einsum('n,nij->ij', fractions/fractions.sum(), self.values[sl])

    def field(self, rho):
        if not -1. <= rho <= 1.: raise ValueError("archived field is restricted to the trapped cell")
        return self.side_value(rho, int(rho >= 0.))

    def field_derivative(self, rho):
        a, source, _, _ = self.owner.coefficients(rho, self.z, np.zeros(4))
        return a@self.field(rho)+source

    @property
    def seam_residual(self):
        return float(np.linalg.norm(self.side_value(0., 0)-self.side_value(0., 1)))


def frequency_transfer_weak(owner, outgoing_adjoint, incoming, nodes, weights):
    """-J†R(zo) deltaH R(zi)J for four smooth compact spatial vertices.

    Different real energies implement exp[-i(zo-zi)t]. Their imaginary parts
    must coincide for that harmonic perturbation. The adjoint outgoing field
    is R(zo*)J. Contact terms at fixed canonical endpoints do not vary here.
    """
    if abs(outgoing_adjoint.z.imag+incoming.z.imag) > 1e-13:
        raise ValueError("equal damping for the two real-frequency ports required")
    kernel = []
    for rho in nodes:
        u, up = incoming.field(rho), incoming.field_derivative(rho)
        l, lp = outgoing_adjoint.field(rho), outgoing_adjoint.field_derivative(rho)
        kernel.append(static_metric_kernel_from_fields(
            rho, owner.compact_mass, owner.angular_eigenvalue, u, up, l, lp))
    from .nsc_transmitting_resolvent import profile
    shape = np.array([profile(x)[0] for x in nodes])
    return np.einsum('n,n,naij->aij', weights, shape, np.array(kernel))


def frequency_transfer_driven(owner, z_out, incoming):
    """Solve (H-zo) dU=-deltaH R(zi)J on the same transmitting domain.

    Only a radial response is integrated, not a metric/quantum-state trajectory.
    Inflow is causal; no condition is imposed at the outflow endpoint.
    """
    z_out = complex(z_out)
    if z_out.imag == 0 or abs(z_out.imag-incoming.z.imag) > 1e-13:
        raise ValueError("same nonzero imaginary part required")
    order = ((1., 0.), (0., -1.)) if z_out.imag > 0 else ((-1., 0.), (0., 1.))
    orientation = -1. if z_out.imag > 0 else 1.
    state = np.zeros(4*8+4*16, dtype=complex)

    def rhs(rho, values):
        y = values[:32].reshape(4, 2, 4)
        a, _, _, _ = owner.coefficients(rho, z_out, np.zeros(4))
        v, _, _, dv, dvp, dm = owner.metric_parts(rho, np.zeros(4))
        u, up = incoming.field(rho), incoming.field_derivative(rho)
        delta_h_u = np.array([-1j*(d@up+0.5*dp@u)+m@u
                              for d, dp, m in zip(dv, dvp, dm)])
        dy = np.einsum('ij,ajk->aik', a, y)-1j*np.linalg.solve(v, delta_h_u)
        integral = orientation*np.einsum('ij,ajk->aik', probes(rho).conj().T, y)
        return np.r_[dy.ravel(), integral.ravel()]

    seam = 0.
    for start, stop in order:
        sol = solve_ivp(rhs, (start, stop), state, method='DOP853', rtol=owner.rtol, atol=owner.atol)
        if not sol.success: raise ArithmeticError(sol.message)
        if start == 0.: seam = float(np.max(abs(sol.y[:, 0]-state)))
        state = sol.y[:, -1]
    return state[32:].reshape(4, 4, 4), seam


def ctp_completion_inventory():
    """Precise C1-C4 availability; absent physical data are not numerical zeros."""
    return {
        'retarded_and_advanced_memory': 'evaluated in the inherited packet basis',
        'harmonic_metric_transfer': 'evaluated as two-frequency response, not a chosen history',
        'full_PG_initial_covariance_and_packet_complement_correlations': None,
        'lesser_greater_and_Keldysh': None,
        'unitary_Vplus_Vminus_on_physical_history': None,
        'KS_endpoint_coordinate_pullback': None,
        'EndpointBranchJets_transmitting_fields': None,
        'Gamma_rest_metric_normal_derivative': None,
        'physical_two_sided_Weyl_mismatch': None,
        'C1': 'retarded-memory embedding PASS; full state-dependent CTP embedding OPEN',
        'C2': 'OPEN', 'C3': 'OPEN', 'C4': 'OPEN',
        'metric_timestep_started': False, 'stationary_history_selected': False,
    }
