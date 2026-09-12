"""Blockwise Dirac evolution for declared homogeneous KS history providers.

This is the conditional U[history] owner.  It proves unitarity and CAR
push-forward for any supplied finite history.  It does not select a physical
history between non-isometric endpoint geometries.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np


@dataclass(frozen=True)
class KSCauchyHistory:
    time: np.ndarray
    a_parallel: np.ndarray
    radius: np.ndarray
    lapse: np.ndarray
    shift: np.ndarray
    label: str

    def validate(self):
        arrays = tuple(np.asarray(getattr(self, name), dtype=float) for name in (
            "time", "a_parallel", "radius", "lapse", "shift"
        ))
        count = len(arrays[0])
        if count < 3 or any(value.shape != (count,) for value in arrays):
            raise ValueError("one equal-length scalar history per metric field required")
        if not all(np.isfinite(value).all() for value in arrays):
            raise ValueError("finite KS history required")
        if np.any(np.diff(arrays[0]) <= 0):
            raise ValueError("strictly increasing history time required")
        if np.min(arrays[1]) <= 0 or np.min(arrays[2]) <= 0 or np.min(arrays[3]) <= 0:
            raise ValueError("positive scale factors and lapse required")
        return count


def endpoint_control_history(
    *, duration, points, seed_radius, selected_radius, a0,
    H_parallel, landau_velocity, landau_gamma, label,
):
    """A smooth endpoint-identical diagnostic history, not a field solution."""
    if not all(isfinite(x) for x in (
        duration, seed_radius, selected_radius, a0, H_parallel,
        landau_velocity, landau_gamma,
    )) or duration <= 0 or points < 3:
        raise ValueError("finite positive control-history data required")
    time = np.linspace(0.0, duration, points)
    s = time/duration
    smooth = 3*s*s-2*s*s*s
    # f(0)=f(1)=0 and df/dt=1 at both endpoints.  This retains the
    # recorded a and H_parallel endpoint jet for every duration.
    hermite = s-3*s*s+2*s*s*s
    a = a0*(1.0+H_parallel*duration*hermite)
    radius = seed_radius+(selected_radius-seed_radius)*smooth
    lapse = 1.0+(1.0/landau_gamma-1.0)*smooth
    endpoint_shift = -landau_velocity/a0
    shift = endpoint_shift*smooth
    history = KSCauchyHistory(time, a, radius, lapse, shift, label)
    history.validate()
    return history


def propagate_blockwise(arrays, history: KSCauchyHistory, *, seed_a, seed_r):
    """Midpoint time-ordered exponential of the declared Hermitian H_ADM."""
    history.validate()
    hx0 = np.asarray(arrays["hx_seed"], dtype=float)
    hy0 = np.asarray(arrays["hy_seed"], dtype=float)
    hz0 = np.asarray(arrays["hz_seed"], dtype=float)
    if not (hx0.shape == hy0.shape == hz0.shape):
        raise ValueError("mode Hamiltonian arrays differ")
    count = len(hx0)
    unitary = np.broadcast_to(np.eye(2, dtype=complex), (count, 2, 2)).copy()
    momentum = hz0*seed_a
    identity = np.eye(2, dtype=complex)
    for index, dt in enumerate(np.diff(history.time)):
        a = (history.a_parallel[index]+history.a_parallel[index+1])/2
        radius = (history.radius[index]+history.radius[index+1])/2
        lapse = (history.lapse[index]+history.lapse[index+1])/2
        shift = (history.shift[index]+history.shift[index+1])/2
        hx = lapse*hx0
        hy = lapse*hy0*seed_r/radius
        hz = lapse*momentum/a
        h_identity = -shift*momentum
        norm = np.sqrt(hx*hx+hy*hy+hz*hz)
        angle = dt*norm
        cosine = np.cos(angle)
        sine_over = np.sinc(angle/np.pi)*dt
        phase = np.exp(-1j*dt*h_identity)
        step = np.empty((count, 2, 2), complex)
        step[:, 0, 0] = phase*(cosine-1j*sine_over*hz)
        step[:, 1, 1] = phase*(cosine+1j*sine_over*hz)
        step[:, 0, 1] = phase*(-1j*sine_over*hx-sine_over*hy)
        step[:, 1, 0] = phase*(-1j*sine_over*hx+sine_over*hy)
        unitary = np.einsum("fij,fjk->fik", step, unitary)
    covariance = np.einsum(
        "fij,fjk,flk->fil", unitary,
        arrays["covariance_seed"], unitary.conj(),
    )
    defect = unitary.swapaxes(1, 2).conj()@unitary-identity
    return {
        "unitary": unitary,
        "covariance": covariance,
        "maximum_unitarity_residual": float(np.max(abs(defect))),
        "minimum_covariance_eigenvalue": float(np.linalg.eigvalsh(covariance).min()),
        "maximum_covariance_eigenvalue": float(np.linalg.eigvalsh(covariance).max()),
    }
