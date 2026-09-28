"""Independent spatial bulk subspaces on the inherited common PG time slice.

The L2 split is orthogonal. The transmitting Dirac operator domain is not a
product of independent half-line domains. This distinction is retained by the
API: no given B, trace duplication, or delta coefficient is accepted as a link.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np

from .nsc_transmitting_dirac_domain import I2, S1, S2, S3


class BulkOperatorDomainError(ValueError):
    pass


@dataclass(frozen=True)
class CommonTimeBulkSplit:
    """H_tau = L2(rho>0,C2) direct_sum L2(rho<0,C2), per retained channel.

    rho=0 has zero L2 measure and is handled by one-sided traces in D(H).
    Quadrature projectors are only finite rank/CAR witnesses for this continuum
    definition; they do not convert the seed frequency covariance into spatial
    Cauchy data or introduce a finite-box Hamiltonian boundary condition.
    """

    seam_rho: float = 0.0

    def __post_init__(self):
        if self.seam_rho != 0.0:
            raise BulkOperatorDomainError("this definition uses the inherited fixed rho=0 seam")

    def projectors(self, rho):
        x = np.asarray(rho, dtype=float)
        if x.ndim != 1 or not np.isfinite(x).all() or np.any(x == self.seam_rho):
            raise BulkOperatorDomainError("distinct off-seam bulk samples required; the seam is trace data")
        if len(np.unique(x)) != len(x) or not (np.any(x > 0) and np.any(x < 0)):
            raise BulkOperatorDomainError("independent spatial samples on both half-lines required")
        parent = np.repeat(x > self.seam_rho, 2)
        child = ~parent
        identity = np.eye(2*len(x))
        rp, rc = identity[parent], identity[child]
        reorder = np.vstack((rp, rc))
        return {"P_parent": rp.T@rp, "P_child": rc.T@rc,
                "R_parent": rp, "R_child": rc, "canonical_reorder": reorder}

    @staticmethod
    def spatial_half_density_map(weights, radial_scale, radius):
        w, q, r = (np.asarray(v, dtype=float) for v in (weights, radial_scale, radius))
        if w.ndim != 1 or q.shape != w.shape or r.shape != w.shape or not np.isfinite(w+q+r).all() or min(np.min(w), np.min(q), np.min(r)) <= 0:
            raise ValueError("positive finite spatial weights, q and r required")
        gram = np.repeat(w*q*r*r, 2)
        return 1/np.sqrt(gram), gram

    @staticmethod
    def principal(N, q_PG, beta):
        if not all(isfinite(x) for x in (N, q_PG, beta)) or min(N, q_PG) <= 0:
            raise ValueError("finite positive lapse and radial scale required")
        return (N/q_PG)*S2-beta*I2

    @staticmethod
    def local_potential(N, radius, compact_mass, angular_eigenvalue):
        """Owned reduced multiplication term in T's sigma2-current convention."""
        if not np.isfinite([N, radius, compact_mass, angular_eigenvalue]).all() or min(N, radius) <= 0:
            raise ValueError("finite retained mode and positive metric required")
        return -N*compact_mass*S1-N*angular_eigenvalue/radius*S3

    def apply_local_expression(self, psi, dpsi, *, N, q_PG, beta, radius,
                               N_prime, q_prime, beta_prime,
                               compact_mass, angular_eigenvalue):
        """Apply the imported symmetric Dirac expression to a smooth local jet."""
        p, dp = np.asarray(psi, complex), np.asarray(dpsi, complex)
        if p.shape != (2,) or dp.shape != (2,) or not np.isfinite(p+dp).all():
            raise ValueError("a finite two-component smooth field and derivative required")
        if not np.isfinite([N_prime, q_prime, beta_prime]).all():
            raise ValueError("finite metric derivatives required")
        v = self.principal(N, q_PG, beta)
        dv = (N_prime/q_PG-N*q_prime/q_PG**2)*S2-beta_prime*I2
        return -1j*(v@dp+0.5*dv@p)+self.local_potential(N, radius, compact_mass, angular_eigenvalue)@p

    def projection_domain_defect(self, *, N, q_PG, beta):
        v = self.principal(N, q_PG, beta)
        singular = np.linalg.svd(v, compute_uv=False)
        return {"delta_coefficient": 1j*v,
                "smallest_singular_value": float(singular.min()),
                "largest_singular_value": float(singular.max()),
                "identity": "H(P_child psi)-P_child Hpsi = i v(0) psi(0) delta(rho)",
                "P_child_preserves_full_transmitting_domain": False,
                "ordinary_link_matrix": None}

    def ordinary_hamiltonian_link(self):
        raise BulkOperatorDomainError(
            "P_child is a valid L2 projector but does not preserve D(H) at a noncharacteristic transmitting seam; choose a domain-respecting weak/resolvent realization before forming a finite B"
        )

    def attach_seed_trace_covariance(self, _covariance):
        raise BulkOperatorDomainError(
            "the Sigma0 frequency covariance is not a covariance on a PG tau spatial slice; a physical Cauchy-state map is required"
        )


def require_common_time_bulk_split(value):
    if not isinstance(value, CommonTimeBulkSplit):
        raise BulkOperatorDomainError("bulk support projectors are required; neither a trace graph nor a given B defines them")
    return value
