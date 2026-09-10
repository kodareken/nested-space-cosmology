"""Bind a scale-inherited transparent LLL state to the actual NSC neck.

This applies the already imported conformal stress and PG projection.  It is
an interface calculation: no gravity equation, scattering problem or angular
sum is rerun.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import pi

import numpy as np

from .nsc_horizon_source import conformal_stress, pg_source


@dataclass(frozen=True)
class RecursiveLLLSource:
    magnetic_flux: int
    omega: float
    surface_gravity: float
    outgoing_state_constant: float
    incoming_state_constant: float
    parent_killing_power: float
    pg_source: dict[str, float]
    required_null_sign: str
    source_null_signs_match: bool

    def to_dict(self):
        return asdict(self)


def inherited_state_constants(magnetic_flux: int, surface_gravity: float, omega: float):
    """State constants when H_child=Omega H_parent in one energy frame."""
    if (isinstance(magnetic_flux, bool) or not isinstance(magnetic_flux, int)
            or magnetic_flux == 0):
        raise ValueError("nonzero integer magnetic flux required")
    if not np.isfinite(surface_gravity) or surface_gravity <= 0:
        raise ValueError("positive parent surface gravity required")
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError("positive inherited scale ratio required")
    outgoing = abs(magnetic_flux) * surface_gravity**2 / (48.0 * pi)
    incoming = omega**2 * outgoing
    return outgoing, incoming


def bind_recursive_lll_to_neck(
    *,
    magnetic_flux: int,
    omega: float,
    surface_gravity: float,
    metric_A: float,
    metric_A_prime: float,
    metric_A_second: float,
    beta: float,
    radius: float,
) -> RecursiveLLLSource:
    """Evaluate the recursively scaled LLL state on one nonhorizon neck."""
    outgoing, incoming = inherited_state_constants(
        magnetic_flux, surface_gravity, omega
    )
    t_uu, t_uv, t_vv = conformal_stress(
        metric_A,
        metric_A_prime,
        metric_A_second,
        outgoing,
        incoming,
        central_charge=abs(magnetic_flux),
    )
    projected = pg_source(
        metric_A, beta, radius, t_uu, t_uv, t_vv
    )
    projected = {key: float(value) for key, value in projected.items()}
    signs_match = (
        projected["fourD_null_plus"] < 0
        and projected["fourD_null_minus"] < 0
    )
    return RecursiveLLLSource(
        magnetic_flux=magnetic_flux,
        omega=omega,
        surface_gravity=surface_gravity,
        outgoing_state_constant=outgoing,
        incoming_state_constant=incoming,
        parent_killing_power=projected["parent_Killing_power"],
        pg_source=projected,
        required_null_sign="both negative on the stored black-universe neck",
        source_null_signs_match=signs_match,
    )

