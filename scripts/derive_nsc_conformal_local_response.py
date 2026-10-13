#!/usr/bin/env python3
"""Finite causal memory primitives for the retained conformal response chain."""
from __future__ import annotations

import argparse
import json

import numpy as np


def finite_memory_kernel(link, outside_hamiltonian, times):
    """Return B exp(-i H_out t) B† without evolving any geometry."""
    link = np.asarray(link, dtype=complex)
    outside = np.asarray(outside_hamiltonian, dtype=complex)
    if (outside.ndim != 2 or outside.shape[0] != outside.shape[1]
            or link.ndim != 2 or link.shape[1] != outside.shape[0]):
        raise ValueError("square exterior Hamiltonian and compatible link required")
    if not np.allclose(outside, outside.conj().T, atol=1e-12, rtol=0):
        raise ValueError("Hermitian exterior Hamiltonian required")
    energies, vectors = np.linalg.eigh(outside)
    rotated = link @ vectors
    return np.einsum("ik,tk,jk->tij", rotated,
                     np.exp(-1j*np.asarray(times)[:, None]*energies), rotated.conj())


def initial_exterior_drive(link, outside_hamiltonian, exterior_initial, times):
    """Retain the exterior initial-state drive alongside the memory kernel."""
    outside = np.asarray(outside_hamiltonian, dtype=complex)
    initial = np.asarray(exterior_initial, dtype=complex)
    if (outside.ndim != 2 or outside.shape[0] != outside.shape[1]
            or initial.shape[0] != outside.shape[0]
            or not np.allclose(outside, outside.conj().T, atol=1e-12, rtol=0)):
        raise ValueError("Hermitian exterior and compatible initial state required")
    energies, vectors = np.linalg.eigh(outside)
    coefficients = vectors.conj().T @ initial
    phases = np.exp(-1j*np.asarray(times)[:, None]*energies)
    return (np.asarray(link) @ vectors @ (phases*coefficients).T).T


def memory_chain():
    """Keep the original measured domains and omitted-channel interpretation."""
    return {
        "pilot_geometry_time": 0.05,
        "successor_window": [0.05, 0.085],
        "readout": "stored-geometry occupation response; geometry is not reintegrated",
        "causal_equation": "i psi_A'=H_AA psi_A-i integral B exp[-i H_BB(t-s)] B† psi_A(s) ds+B exp[-i H_BB t] psi_B(0)",
        "state_average_requires": ["exterior preparation", "initial cross correlations"],
        "controls": ["outside drive omitted", "memory omitted", "initial cross omitted"],
        "comparison": "conditional reduction versus stored autonomous response",
        "occupation_error_target": "one percent of the measured occupation change",
        "physical_local_gate": "OPEN",
        "data": "https://doi.org/10.5281/zenodo.23090015",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--describe", action="store_true", help="show the finite chain and controls")
    parser.parse_args()
    print(json.dumps(memory_chain(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
