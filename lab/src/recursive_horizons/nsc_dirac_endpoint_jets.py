"""Endpoint derivatives of the already owned midpoint KS Dirac propagator.

The Fréchet exponential differential is imported from SciPy. This application
supplies only the known frequency-diagonal bulk contribution on a supplied
history. Transmitting link and hypersurface derivatives remain separate,
unevaluated physical inputs; these bulk jets cannot authorize evolution.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm, expm_frechet

from .nsc_transmitting_ctp_variation import EndpointBranchJets


I = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=complex)


def bulk_hamiltonian_and_vertices(hx0, hy0, momentum, N, beta, a, r, seed_r):
    """The exact raw-field vertices of nsc_landau_cauchy_isometry's H.

Field order is N,beta,q_ADM(=a),r. No derivative of an unspecified
transmitting B[g,embedding] or of the Cauchy normal is inserted here.
"""
    if not np.isfinite([hx0, hy0, momentum, N, beta, a, r, seed_r]).all() or min(N, a, r, seed_r) <= 0:
        raise ValueError("finite mode data and positive metric fields required")
    base = hx0*PAULI[0] + hy0*seed_r/r*PAULI[1] + momentum/a*PAULI[2]
    h = N*base-beta*momentum*I
    vertices = np.array([base, -momentum*I,
                         -N*momentum/a**2*PAULI[2],
                         -N*hy0*seed_r/r**2*PAULI[1]])
    return h, vertices


def bulk_endpoint_jets(arrays, history, *, seed_a, seed_r, basis_id):
    """Differentiate the fixed midpoint product at its two endpoint nodes.

Each endpoint affects half of the first/last midpoint value. The derivative
of a step exponential is propagated in the same matrix order as U. There is
no change of duration, no metric stepping and no history selector here.
"""
    history.validate()
    count = len(arrays["hx_seed"])
    steps = len(history.time)-1
    total = np.empty((count, 2, 2), dtype=complex)
    derivatives = np.empty((count, 4, 2, 2, 2), dtype=complex)
    linear_error = tangent_error = unitarity_error = 0.0
    for mode in range(count):
        u = I.copy()
        du = np.zeros((4, 2, 2, 2), dtype=complex)
        for node, dt in enumerate(np.diff(history.time)):
            N, beta, a, r = (float((getattr(history, name)[node]+getattr(history, name)[node+1])/2)
                             for name in ("lapse", "shift", "a_parallel", "radius"))
            h, vertices = bulk_hamiltonian_and_vertices(
                arrays["hx_seed"][mode], arrays["hy_seed"][mode],
                arrays["hz_seed"][mode]*seed_a, N, beta, a, r, seed_r)
            z = -1j*dt*h
            step = expm(z)
            next_du = np.einsum("ij,fejk->feik", step, du)
            for side in range(2):
                if node != (0 if side == 0 else steps-1):
                    continue
                for field in range(4):
                    # Existing midpoint convention: dH_mid/dg_end = vertex/2.
                    dstep = expm_frechet(z, -0.5j*dt*vertices[field], compute_expm=False)
                    next_du[field, side] += dstep@u
                    linear_error = max(linear_error, float(np.max(np.abs(dstep.conj().T@step+step.conj().T@dstep))))
            u = step@u
            du = next_du
        total[mode], derivatives[mode] = u, du
        unitarity_error = max(unitarity_error, float(np.max(np.abs(u.conj().T@u-I))))
        for row in du.reshape(-1, 2, 2):
            tangent_error = max(tangent_error, float(np.max(np.abs(row.conj().T@u+u.conj().T@row))))
    jets = [EndpointBranchJets(
        basis_id, 0.5*du, -0.5*du,
        "known KS bulk on frozen diagnostic history; transmitting link and embedding NOT supplied",
    ) for du in derivatives]
    return {
        "unitary": total, "dU_dendpoint": derivatives, "branch_jets": jets,
        "step_tangent_residual": linear_error,
        "final_tangent_residual": tangent_error,
        "unitarity_residual": unitarity_error,
        "physical_jets_complete": False,
        "transmitting_link_derivative": None,
        "transmitting_embedding_derivative": None,
    }
