"""First variation of the existing finite Gaussian transmitting CTP overlap.

The matrix differential is imported; the NSC application binds it to the
retained Cauchy channels and their existing endpoint covector. A differential
with respect to a supplied unitary is not a law selecting that unitary.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np

from .nsc_influence import _covariance
from .nsc_transmitting_boundary_history import EndpointVariation, FIELDS


def _matrix(value, shape, label):
    matrix = np.asarray(value, dtype=complex)
    if matrix.shape != shape or not np.isfinite(matrix).all():
        raise ValueError(f"finite {label} with shape {shape} required")
    return matrix


def ctp_first_variation(covariance, plus, minus, dplus, dminus, *, tolerance=3e-11):
    """d[-i log det Q], Q=I-C+C Uminus^dagger Uplus, with fixed C0.

Both arguments dplus/dminus must be tangents at their respective unitaries.
The returned number is a finite bare Gaussian action differential. Reference,
local, boundary and coordinate-embedding terms are separate owned inputs.
"""
    if not isfinite(tolerance) or tolerance <= 0:
        raise ValueError("positive finite tolerance required")
    c = _covariance(covariance)
    shape = c.shape
    up, um, dp, dm = (_matrix(v, shape, name) for v, name in zip(
        (plus, minus, dplus, dminus), ("Uplus", "Uminus", "dUplus", "dUminus")))
    identity = np.eye(len(c), dtype=complex)
    unitary_residual = max(float(np.max(np.abs(u.conj().T@u-identity))) for u in (up, um))
    tangent_residual = max(float(np.max(np.abs(du.conj().T@u+u.conj().T@du)))
                           for u, du in ((up, dp), (um, dm)))
    if unitary_residual > tolerance or tangent_residual > tolerance:
        raise ValueError("unitary branches and tangent-compatible variations required")
    q = identity-c+c@um.conj().T@up
    dq = c@(dm.conj().T@up+um.conj().T@dp)
    solved = np.linalg.solve(q, dq)
    linear_residual = float(np.max(np.abs(q@solved-dq)))
    derivative = complex(-1j*np.trace(solved))
    return {
        "derivative": derivative,
        "unitarity_residual": unitary_residual,
        "tangent_residual": tangent_residual,
        "linear_solve_residual": linear_residual,
    }


@dataclass(frozen=True)
class EndpointBranchJets:
    """Supplied physical chain rule, on the existing raw KS endpoint basis.

Shape: (four fields, initial/final endpoint, channel row, channel column).
Each derivative is with respect to the difference-history coordinate g_Delta.
Its origin in U[g,B,embedding] must be supplied by the transmitting action.
"""

    basis_id: str
    plus: np.ndarray
    minus: np.ndarray
    physical_derivative_owner: str


def endpoint_ctp_pullback(covariance, plus, minus, endpoint_basis, jets, *, tolerance=3e-11):
    """Contract one channel's actual branch jets into EndpointVariation.

This supplies the Gaussian part only. It cannot certify a two-sided metric
match or replace reference/local/interface allocations with zero.
"""
    if jets is None:
        return {
            "status": "OPEN", "gaussian_endpoint_covector": None,
            "maximum_residual": None, "tolerance": tolerance,
            "reason": "dUplus/dg_Delta and dUminus/dg_Delta on the transmitting embedding are not supplied",
        }
    if jets.basis_id != endpoint_basis.basis_id or not jets.physical_derivative_owner:
        raise ValueError("matching KS endpoint basis and derivative provenance required")
    c = _covariance(covariance)
    expected = (len(FIELDS), 2, *c.shape)
    jp, jm = (_matrix(v, expected, name) for v, name in ((jets.plus, "plus branch jets"), (jets.minus, "minus branch jets")))
    if np.max(np.abs(np.asarray(plus)-np.asarray(minus))) > tolerance:
        raise ValueError("endpoint metric source requires the physical equal-history limit")
    values = np.zeros((len(FIELDS), 2), dtype=float)
    maximum = 0.0
    for field in range(len(FIELDS)):
        for side in range(2):
            result = ctp_first_variation(c, plus, minus, jp[field, side], jm[field, side], tolerance=tolerance)
            value = result["derivative"]
            if abs(value.imag) > tolerance:
                raise ArithmeticError("physical equal-history action derivative is not real")
            values[field, side] = value.real
            maximum = max(maximum, result["linear_solve_residual"], result["tangent_residual"], abs(value.imag))
    covector = EndpointVariation(endpoint_basis.basis_id, tuple(tuple(row) for row in values))
    return {
        "status": "CONDITIONAL EVALUATION PASS",
        "gaussian_endpoint_covector": covector.to_dict(),
        "maximum_residual": maximum, "tolerance": tolerance,
        "physical_derivative_owner": jets.physical_derivative_owner,
        "remaining_same_action_derivative": None,
        "two_sided_physical_match": None,
    }


def completion_gate(previous, selection, metric_pullback):
    """Re-enter the extended dependency gate without changing old certificates."""
    old = previous["composition"]
    if not old["old_homogeneous_nonexistence"]:
        raise ValueError("homogeneous certificate dependency changed")
    return {
        "decision": "OPEN",
        "F": selection["status"], "G": metric_pullback["status"],
        "missing_residuals": [
            "physical link/embedding variations that constrain the admissible CTP unitary directions",
            "dU/dg_Delta and the same-action non-Gaussian boundary derivative on the KS endpoint/jet basis",
            "the full history stationarity residual after those derivatives are supplied",
        ],
        "computed_differential": "finite Gaussian first variation, conditional on branch tangents",
        "full_metric_variation": {field: None for field in FIELDS},
        "homogeneous_nonexistence_preserved": True,
        "old_generator_rerun": False, "existence_claimed": False,
        "extended_nonexistence_claimed": False, "optimizer_started": False,
        "finite_stress": None, "nulls": None, "updated_constraints": None,
        "metric_timestep_started": False, "coupled_evolution_reopened": False,
    }
