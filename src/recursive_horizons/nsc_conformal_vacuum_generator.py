"""Vacuum-volume comparison direction for the owned spherical action.

This does not add a term to the evolved theory. In its (+---) convention
S = integral sqrt(|g|) (-A R_L - 2 A Lambda + ...), A=1/(16 pi G).
N=rL, q=rQ give H_Lambda = 8 pi A Lambda integral r^4 LQ dx.
The existing r^2 term in feedback_V is sphere curvature, not this generator.
"""
from __future__ import annotations

import numpy as np

from .nsc_spherical_coupling import chart_failure, PositiveChartExit
from .nsc_spherical_galerkin_coupling import (
    active_fine_system, prolong_state, pull_geometry,
)


def ungauged_volume_jet(r, L, Q, A):
    """Unit-Lambda Hamiltonian density and derivatives before choosing L=Q."""
    r, L, Q = np.broadcast_arrays(np.asarray(r, dtype=float),
                                  np.asarray(L, dtype=float), np.asarray(Q, dtype=float))
    if (not np.isfinite(A) or A <= 0 or
            any(not np.isfinite(v).all() or np.min(v) <= 0 for v in (r, L, Q))):
        raise ValueError("positive finite volume factors and Einstein coefficient required")
    factor = 8 * np.pi * A
    return {
        "density": factor * r**4 * L * Q,
        "dL": factor * Q * r**4,
        "dQ": factor * L * r**4,
        "dr": 4 * factor * L * Q * r**3,
    }


def lambda_generator(grid, state):
    """Energy and canonical log-geometry vector of one unit of Lambda.

    Adjoint pullback precedes the nodal canonical change. In particular,
    multiplying by r on the quadrature grid then pulling back would give
    a different finite Hamiltonian. No matter or coordinate rate changes.
    """
    if grid.gauge != "conformal":
        raise ValueError("this generator uses the owned conformal L=Q gauge")
    fine = prolong_state(grid, state)
    system = active_fine_system(grid, fine)
    failure = chart_failure(system, fine)
    if failure:
        raise PositiveChartExit(failure, np.nan, state.copy())
    jet = ungauged_volume_jet(fine.r, system.length_density, fine.Q, system.A)
    zeros = np.zeros(grid.ng)
    # The total Q derivative includes dL/dQ=1 after varying the constraint.
    p_Q_dot = -pull_geometry(grid, jet["dQ"] + jet["dL"])
    p_r_dot = -pull_geometry(grid, jet["dr"])
    return {
        "energy": float(grid.dx_q * np.sum(jet["density"])),
        "canonical_rates": {"a": zeros.copy(), "b": zeros.copy(), "chi": zeros.copy(),
                            "pi_a": state.r * p_r_dot, "pi_b": state.Q * p_Q_dot,
                            "p_chi": zeros.copy()},
        "lapse_constraint_density": jet["dL"],
        "shift_constraint_density": np.zeros(grid.nq),
        "matter_rate_is_zero": True,
    }
