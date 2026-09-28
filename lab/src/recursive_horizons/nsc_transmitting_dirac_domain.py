"""Retained-channel transparent Dirac trace domain at the stored spacelike seam.

The two traces restrict one field to rho=0 from its two sides. They are not
independent canonical rooms. Geometry, the recorded channel basis and the
known child/PG coframe determine the trace map; no B matrix is accepted.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import sqrt, atanh, cosh, sinh

import numpy as np


I2 = np.eye(2, dtype=complex)
S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
S3 = np.diag([1., -1.]).astype(complex)
MODE_TO_CURRENT = (I2+1j*S1)/sqrt(2)


class DomainRepresentationError(ValueError):
    pass


@dataclass(frozen=True)
class TransmittingDiracSeamDomain:
    """Same-surface finite spectral domain of the existing rho=0 transmission.

PG coframe: theta0=N dtau, theta1=q(d rho+beta dtau).
The rho=0 hypersurface is spacelike when q*beta>N. Its future normal points
toward decreasing rho. The metric is smooth across the cut in one spin/gauge
frame. This defines a fixed-coordinate seam, not a moving-throat selector.
"""

    lapse: float
    radial_scale: float
    shift: float
    radius: float
    surface_id: str = "stored-Bronnikov-rho0-spacelike-seed"

    def __post_init__(self):
        if not np.isfinite([self.lapse, self.radial_scale, self.shift, self.radius]).all():
            raise ValueError("finite seam geometry required")
        if min(self.lapse, self.radial_scale, self.radius) <= 0 or self.radial_scale*self.shift <= self.lapse:
            raise DomainRepresentationError("the declared seed seam requires positive metric factors and q*beta>N")

    @property
    def induced_axial_scale(self):
        return sqrt((self.radial_scale*self.shift)**2-self.lapse**2)

    @property
    def coframe_rapidity(self):
        return atanh(self.lapse/(self.radial_scale*self.shift))

    def normal_geometry(self):
        N, q, beta = self.lapse, self.radial_scale, self.shift
        a = self.induced_axial_scale
        metric = np.array([[N*N-q*q*beta*beta, -q*q*beta], [-q*q*beta, -q*q]])
        # Future normal at the spacelike cut, in (tau,rho) coordinates.
        normal = np.array([q*beta/(N*a), -a/(N*q)])
        tangent = np.array([1., 0.])
        current_metric = (q*beta*I2-N*S2)/a
        return {
            "orbit_metric": metric, "future_normal": normal,
            "parent_outward_normal": normal, "child_outward_normal": -normal,
            "normal_spin_metric": current_metric,
            "normal_unit_residual": float(abs(normal@metric@normal-1)),
            "normal_tangent_residual": float(abs(normal@metric@tangent)),
            "induced_metric_residual": float(abs(tangent@metric@tangent+a*a)),
        }

    def bind_seed(self, mode_state, *, tolerance=3e-11):
        if not isinstance(mode_state, dict) or "cauchy_surfaces" not in mode_state or "channels" not in mode_state:
            raise DomainRepresentationError("a declared channel/Cauchy domain is required; given LinkHistory.values are insufficient")
        seed = mode_state["cauchy_surfaces"]["seed"]
        if abs(self.radius-seed["r"]) > tolerance or abs(self.induced_axial_scale-seed["a_parallel"]) > tolerance:
            raise DomainRepresentationError("this stored covariance belongs to Sigma0; another surface requires its Cauchy transport")

    def trace_map(self, weights):
        """Canonical coefficients -> PG spinor trace coefficients on this seam.

W is scalar frequency quadrature only. The stored parent-neck frequency axis
already carries its inherited state occupations; Omega is not applied again.
The boost is a coframe change on THIS surface, never U_L0 between surfaces.
"""
        w = np.asarray(weights, dtype=float)
        if w.ndim != 1 or len(w) == 0 or not np.isfinite(w).all() or min(w) <= 0:
            raise ValueError("positive finite stored quadrature weights required")
        a, xi = self.induced_axial_scale, self.coframe_rapidity
        inverse_spin_half = cosh(xi/2)*I2+sinh(xi/2)*S2
        spin_half = cosh(xi/2)*I2-sinh(xi/2)*S2
        scale = self.radius*np.sqrt(a*w)
        trace = (inverse_spin_half@MODE_TO_CURRENT)[None, :, :]/scale[:, None, None]
        inverse = scale[:, None, None]*(MODE_TO_CURRENT.conj().T@spin_half)[None, :, :]
        gram = (self.radius**2*a*w)[:, None, None]*self.normal_geometry()["normal_spin_metric"]
        return trace, inverse, gram

    def evaluate_channel(self, weights, covariance):
        trace, inverse, gram = self.trace_map(weights)
        c = np.asarray(covariance, dtype=complex)
        if c.shape != trace.shape or not np.isfinite(c).all():
            raise ValueError("one finite stored 2x2 covariance per frequency node required")
        adjoint = lambda a: a.swapaxes(-1, -2).conj()
        trace_covariance = trace@c@adjoint(trace)
        recovered = inverse@trace_covariance@adjoint(inverse)
        normalized = adjoint(trace)@gram@trace
        maximum_green = 0.0
        for t, g in zip(trace, gram):
            # Oriented parent future boundary, child past boundary.
            graph = np.vstack((t, t))
            form = np.block([[g, np.zeros((2, 2))], [np.zeros((2, 2)), -g]])
            maximum_green = max(maximum_green, float(np.max(np.abs(graph.conj().T@form@graph))))
        eig = np.linalg.eigvalsh(recovered)
        return {
            "trace_map": trace, "trace_inverse": inverse, "trace_metric": gram,
            "coordinate_trace_covariance": trace_covariance,
            "residuals": {
                "normalization": float(np.max(np.abs(normalized-I2))),
                "inverse_map": float(np.max(np.abs(inverse@trace-I2))),
                "same_surface_green_form": maximum_green,
                "canonical_covariance_recovery": float(np.max(np.abs(recovered-c))),
                "recovered_CAR_lower_violation": max(0., -float(eig.min())),
                "recovered_CAR_upper_violation": max(0., float(eig.max())-1.),
            },
            "independent_coefficients": 2*len(c),
            "doubled_trace_dimension": 4*len(c),
            "matching_constraint_rank": 2*len(c),
            "allowed_graph_dimension": 2*len(c),
        }

    def instantaneous_link_history(self):
        raise DomainRepresentationError(
            "spacelike same-field trace sewing does not define an instantaneous parent/child Hamiltonian block B"
        )


def require_transmitting_domain(value):
    if not isinstance(value, TransmittingDiracSeamDomain):
        raise DomainRepresentationError("given B matrices are not a differentiable transmitting Dirac domain")
    return value
