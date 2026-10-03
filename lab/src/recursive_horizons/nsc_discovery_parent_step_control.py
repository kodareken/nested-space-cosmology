"""Rank-general numerical admission for the unchanged leading-action RHS.

Source columns may be nonorthonormal.  Their actual amplitude weights stay
in the source/Jv owner; this module changes neither state nor dynamics.
"""
from __future__ import annotations

import math

import numpy as np

from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_leading_step_control as legacy

NUMERICAL_MODE = "rank_general_physically_scaled_local_reaction"
physical_scales = legacy.physical_scales
local_system = legacy.local_system


def source_rank(grid, fine, system):
    """Validate the live source width without assuming normalized columns."""
    phi0, phi1 = np.asarray(fine.phi0), np.asarray(fine.phi1)
    if phi0.ndim != 2 or phi0.shape != phi1.shape or phi0.shape[0] != grid.nq:
        raise ValueError("fine spinor arrays must have equal nq-by-rank shapes")
    rank = phi0.shape[1]
    weights = np.asarray(system.occupations)
    if weights.shape != (rank,) or not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("live source requires one finite nonnegative weight per column")
    if not np.isfinite(phi0).all() or not np.isfinite(phi1).all():
        raise ValueError("source columns must be finite")
    if not math.isfinite(float(grid.dx_q)) or grid.dx_q <= 0 or system.dx != grid.dx_q:
        raise ValueError("fine source spacing must equal positive grid dx_q")
    return rank


def scaled_local_reaction(grid, fine, system):
    """Return all row sums of S^-1 J0 S, with 4+4R real components.

    Field order is phi0.real, phi0.imag, phi1.real, phi1.imag.  Zero spatial
    operators isolate reaction terms, retaining both reciprocal source and
    field links.  R=0 has only the four geometric components.
    """
    rank = source_rank(grid, fine, system)
    scales, frequency, kg, kf = physical_scales(grid, fine, system)
    rows = np.zeros((grid.nq, 4 + 4 * rank))
    sqrt_dx = np.sqrt(grid.dx_q)
    local = local_system(system)
    for component in range(rows.shape[1]):
        tangent = leading.State(*(np.zeros_like(getattr(fine, name)) for name in leading.FIELDS))
        if component < 4:
            setattr(tangent, leading.REAL_FIELDS[component], scales[:, component].copy())
        else:
            group, column = divmod(component - 4, rank)
            field = tangent.phi0 if group < 2 else tangent.phi1
            field[:, column] = sqrt_dx * (1 if group % 2 == 0 else 1j)
        image = leading.fine_jvp(local, fine, tangent)
        geometric = np.column_stack([getattr(image, name) for name in leading.REAL_FIELDS]) / scales
        fields = np.column_stack((image.phi0.real, image.phi0.imag,
                                  image.phi1.real, image.phi1.imag)) / sqrt_dx
        rows += abs(np.column_stack((geometric, fields)))
    if not np.isfinite(rows).all():
        raise ValueError("local scaled reaction must be finite")
    return rows, scales, frequency, kg, kf


def scaled_step_admission(pair, state, cap, *, control_mode="coupled"):
    """Pure metadata record with ``dt``; rates use the actual projected state.

    The optional frozen control uses zero geometric scale velocities.  This
    is an admission indicator, not a global stability or error certificate.
    """
    cap = float(cap)
    if not math.isfinite(cap) or cap <= 0:
        raise ValueError("step cap must be finite and positive")
    rate, bundle = leading.rates(pair, state, return_bundle=True, control_mode=control_mode)
    fine, system, grid = bundle["fine_state"], bundle["fine_system"], pair.grid
    rank = source_rank(grid, fine, system)
    rows, _scales, frequency, kg, kf = scaled_local_reaction(grid, fine, system)
    qt = leading.galerkin.prolong_geometry(grid, pair.geometry_map @ rate.Q)
    rt = leading.galerkin.prolong_geometry(grid, pair.geometry_map @ rate.r)
    log_rates = np.column_stack((qt / fine.Q, rt / fine.r,
                                2 * rt / fine.r - qt / fine.Q, rt / fine.r))
    frozen_norm = float(np.max(rows))
    moving_norm = float(np.max(abs(log_rates)))
    bound = frozen_norm + moving_norm
    gradient = 2 * max(float(np.max(abs(system.derivative @ fine.Q / fine.Q))),
                       float(np.max(abs(system.derivative @ fine.r / fine.r))))
    kquad = leading.episode.quadrature_principal_wavenumber(grid)
    omega = max(kquad, kf) + bound + gradient
    if not math.isfinite(omega) or omega < 0:
        raise ValueError("combined step frequency must be finite and nonnegative")
    dt = min(cap, leading.episode.RK4_HALF_STABILITY / max(omega, 1.))
    return {
        "dt": dt, "step_cap": cap, "combined_omega": omega,
        "numerical_mode": NUMERICAL_MODE, "source_rank": rank, "column_rank": rank,
        "local_real_variables": 4 + 4 * rank,
        "frozen_scaled_reaction_norm": frozen_norm,
        "absolute_moving_scale_rate": moving_norm,
        "reaction_majorant": bound, "gradient_majorant": gradient,
        "actual_geometry_band_k": kg, "field_k": kf, "quadrature_k": kquad,
        "scale_frequency": frequency,
        "scales": "Q,r,2|a|r^2 omega*/Q,2|a|r omega*; psi=Phi/sqrt(dx_q)",
        "moving_scale_sign_cancellation_used": False, "physical_growth_removed": False,
        "principal_coordinate_speed": 1., "field_cfl_alone": False,
        "stability_certificate": False, "control_mode": bundle["control_mode"],
        "raw_coordinate_norm_used": False,
        "restriction_owner": "nsc_discovery_parent_step_control.scaled_step_admission",
        "scope": "all source-rank local weighted components plus moving-scale triangle, "
                 "principal and gradient terms; not a global stability certificate",
    }


def step_restriction(pair, state, step_cap, *, control_mode="coupled"):
    """Compatibility API returning ``(dt, record)`` without runtime patching."""
    record = scaled_step_admission(pair, state, step_cap, control_mode=control_mode)
    return record["dt"], record


def stable_timestep(pair, state, step_cap, *, control_mode="coupled"):
    """Float-only convenience API; ``scaled_step_admission`` keeps the record."""
    return scaled_step_admission(pair, state, step_cap, control_mode=control_mode)["dt"]


# The parent episode's tuple-returning name can delegate to this shared owner.
scaled_step_restriction = step_restriction
