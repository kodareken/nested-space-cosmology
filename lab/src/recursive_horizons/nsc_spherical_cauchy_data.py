"""Initial Cauchy data for a nonzero spherical shift current.

The ansatz is the one the existing constraints actually allow: chi = 0,
p_r = 0, p_chi = 0, and Q equal to the positive gauge ratio b0/a0. Then

    D + j = -Q d_x p_Q + j,
    C + rho is independent of p_Q,

with j = force_beta / dx from the existing column source. p_Q is the
mean-zero antiderivative, in the geometry subspace, of the mean-free part of
j/Q. The radius is the existing projected Newton solve for rho = force_L / dx.

This writes the initial point only. It does not evolve, reset a later
state, or change the action.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

from .nsc_spherical_coupling import (
    PERIOD,
    TOL_INITIAL_CONSTRAINT,
    CauchyState,
    hamilton_constraint,
    shift_constraint,
)
from .nsc_spherical_feedback_action import partial_F
from .nsc_spherical_galerkin_coupling import (
    blank_state,
    build_grid,
    constraint_diagnostics,
    geometric_rates,
    load_physical_columns,
    normal_velocities,
    prolong_state,
    pull_geometry,
    solve_initial_radius,
    source_from_columns,
)

SCHEMA = "NSC-SPHERICAL-CAUCHY-DATA-v1"
_LAB_ROOT = Path(__file__).resolve().parents[2]
RECORD_PATH = _LAB_ROOT / "results" / "development" / "nsc-spherical-cauchy-data-v1.json"
_V5_PATH = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-refinement-v5.json"
_ANTIDERIVATIVE_RESIDUAL_LIMIT = 1e-8


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _max_abs(values):
    array = np.asarray(values)
    if array.size == 0:
        return 0.0
    return float(np.max(np.abs(array)))


def gauge_ratio(grid):
    """Spatially constant Q used by the existing radius residual."""
    calibration = grid.fine.calibration
    q0 = float(calibration["b0"] / calibration["a0"])
    if not np.isfinite(q0) or q0 <= 0.0:
        raise ValueError("gauge ratio b0/a0 is not a positive finite constant")
    return q0


def regional_pair_occupations(occupations):
    values = [float(v) for v in np.asarray(occupations, dtype=float).reshape(-1)]
    pairs = []
    equal = len(values) == 6
    for start in (0, 2, 4):
        if start + 1 >= len(values):
            equal = False
            break
        pairs.append(values[start] == values[start + 1])
    return {
        "values": values,
        "equal_within_regional_pairs": bool(equal and all(pairs)),
    }


def _geometry_derivative(grid):
    """Coarse derivative (ng/nq) A_g.T D A_g. Constants are its kernel."""
    return grid.weight * grid.A_g.T @ grid.derivative_on_geometry


def shift_momentum(grid, current):
    """Mean-zero geometry-band p_Q with d_x p_Q matching the mean-free part of j/Q.

    The current array is not modified. A nonzero mean stays in the shift
    residual because a periodic derivative has mean zero.
    """
    current = np.asarray(current, dtype=float)
    if current.shape != (grid.nq,):
        raise ValueError("current must be the quadrature-node array j = force_beta / dx")
    if not np.isfinite(current).all():
        raise ValueError("current is not finite")
    q0 = gauge_ratio(grid)
    desired = pull_geometry(grid, current / q0)
    desired_mean = float(np.mean(desired))
    # Drop only the derivative target’s mean. A periodic p_Q cannot match it.
    rhs = desired - desired_mean
    derivative = _geometry_derivative(grid)
    augmented = np.vstack((derivative, np.ones(grid.ng)))
    target = np.concatenate((rhs, np.zeros(1)))
    solved, residuals, rank, _singular = np.linalg.lstsq(augmented, target, rcond=None)
    momentum = solved - float(np.mean(solved))
    derivative_gap = _max_abs(derivative @ momentum - rhs)
    accepted = bool(derivative_gap <= _ANTIDERIVATIVE_RESIDUAL_LIMIT and rank >= grid.ng)
    return momentum, {
        "Q": q0,
        "current_mean": float(np.mean(current)),
        "current_integral": float(np.mean(current) * grid.length),
        "current_max": _max_abs(current),
        "projected_derivative_mean": desired_mean,
        "source_mean_subtracted": False,
        "p_Q_mean_declared": 0.0,
        "p_Q_mean": float(np.mean(momentum)),
        "p_Q_max": _max_abs(momentum),
        "derivative_gap": derivative_gap,
        "derivative_rank": int(rank),
        "lstsq_residual": float(residuals[0]) if len(residuals) else 0.0,
        "antiderivative_accepted": accepted,
    }


def _fine_ansatz(grid, p_coarse, radius=None):
    q0 = gauge_ratio(grid)
    points = grid.nq
    if radius is None:
        radius = np.ones(points)
    return CauchyState(
        Q=np.full(points, q0),
        r=np.asarray(radius, dtype=float).copy(),
        chi=np.zeros(points),
        p_Q=grid.A_g @ np.asarray(p_coarse, dtype=float),
        p_r=np.zeros(points),
        p_chi=np.zeros(points),
        phi0=np.zeros((points, 1), dtype=np.complex128),
        phi1=np.zeros((points, 1), dtype=np.complex128),
    )


def shift_residual(grid, p_coarse, current):
    """Owner shift constraint plus a supplied current, on the ansatz."""
    return shift_constraint(grid.fine, _fine_ansatz(grid, p_coarse)) + np.asarray(current, dtype=float)


def _require_prepared_state(grid, state):
    occupations = int(np.asarray(grid.fine.occupations).size)
    if state.phi0.shape != state.phi1.shape:
        raise ValueError("phi0 and phi1 shapes differ")
    if state.phi0.shape != (grid.nf, occupations):
        raise ValueError(
            f"columns must have shape {(grid.nf, occupations)}; "
            "sample prepared columns on the fermion grid"
        )
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi"):
        if np.shape(getattr(state, name)) != (grid.ng,):
            raise ValueError(f"{name} must be a coarse geometry array of length {grid.ng}")


def _chi_rate_report(grid, state):
    fine = prolong_state(grid, state)
    _q_dot, r_dot, chi_dot, _p_q_dot, _p_r_dot, _p_chi_dot = geometric_rates(grid.fine, fine)
    _force_r, f_chi = partial_F(fine.r, grid.fine.A, grid.fine.C_W)
    f_chi = float(f_chi)
    if f_chi == 0.0:
        raise ValueError("F_chi vanished; the auxiliary chi rate is outside this chart")
    predicted = grid.fine.length_density * fine.p_Q / (2.0 * f_chi)
    radius_derivative = grid.derivative @ fine.r
    normal_gap = r_dot - grid.fine.shift * radius_derivative
    velocities = normal_velocities(grid, state)
    return {
        "F_chi": f_chi,
        "fine_max": _max_abs(chi_dot),
        "coarse_max": _max_abs(pull_geometry(grid, chi_dot)),
        "gap_against_geometric_rates": _max_abs(chi_dot - predicted),
        "formula": "L p_Q / (2 F_chi) at chi = p_r = p_chi = 0",
        "chart_radius_rate_gap": _max_abs(normal_gap),
        "velocities": velocities,
    }


def initial_cauchy_data(grid, state=None, *, phi0=None, phi1=None):
    """Constraint-reduced initial data for prepared columns on ``grid``.

    Pass a coarse ``CauchyState`` or the fermion-grid columns. Geometry is
    replaced by the ansatz. Columns, occupations, and multiplicity stay on
    the bindings already stored in the state and the grid.
    """
    if state is not None and (phi0 is not None or phi1 is not None):
        raise ValueError("pass a prepared state or columns, not both")
    if state is None:
        if phi0 is None or phi1 is None:
            raise ValueError("pass a prepared state or both columns")
        state = blank_state(grid, phi0, phi1)
    _require_prepared_state(grid, state)
    q0 = gauge_ratio(grid)
    caller_bytes = tuple(
        np.ascontiguousarray(array).tobytes()
        for array in (
            state.Q, state.r, state.chi, state.p_Q, state.p_r, state.p_chi, state.phi0, state.phi1,
        )
    )
    column_bytes = None
    if phi0 is not None:
        column_bytes = (
            np.ascontiguousarray(phi0).tobytes(),
            np.ascontiguousarray(phi1).tobytes(),
        )
    prepared = state.copy()
    prepared.Q = np.full(grid.ng, q0)
    prepared.chi = np.zeros(grid.ng)
    prepared.p_r = np.zeros(grid.ng)
    prepared.p_chi = np.zeros(grid.ng)
    prepared.p_Q = np.zeros(grid.ng)
    fine = prolong_state(grid, prepared)
    source = source_from_columns(grid.fine, fine)
    if source["multiplicity_applied_once"] is not True:
        raise ValueError("source multiplicity was not applied once")
    if int(grid.fine.multiplicity) != 4 * int(grid.fine.kappa):
        raise ValueError("4 kappa multiplicity drifted")
    current = source["force_beta"] / grid.dx_q
    current_snapshot = np.array(current, dtype=float, copy=True)
    momentum, momentum_info = shift_momentum(grid, current_snapshot)
    if _max_abs(current_snapshot - current) != 0.0:
        raise RuntimeError("shift momentum altered the current")
    solved, radius_info = solve_initial_radius(grid, prepared)
    if radius_info.get("converged"):
        data = solved
    else:
        data = prepared
    data.p_Q = momentum
    data.chi = np.zeros(grid.ng)
    data.p_r = np.zeros(grid.ng)
    data.p_chi = np.zeros(grid.ng)
    data.Q = np.full(grid.ng, q0)
    diagnostics = constraint_diagnostics(grid, data)
    fine_data = prolong_state(grid, data)
    owned_momentum = shift_constraint(grid.fine, fine_data) + current_snapshot
    projected_momentum = pull_geometry(grid, owned_momentum)
    beyond_mean = projected_momentum - float(np.mean(projected_momentum))
    zero_momentum = data.copy()
    zero_momentum.p_Q = np.zeros(grid.ng)
    hamilton_gap = _max_abs(
        hamilton_constraint(grid.fine, fine_data)
        - hamilton_constraint(grid.fine, prolong_state(grid, zero_momentum))
    )
    rates = None
    if bool(np.min(fine_data.r) > 0.0 and np.min(fine_data.Q) > 0.0):
        rates = _chi_rate_report(grid, data)
    current_mean = float(np.mean(current_snapshot))
    mean_gap = abs(float(np.mean(projected_momentum)) - current_mean)
    restrictions = []
    if abs(current_mean) > TOL_INITIAL_CONSTRAINT:
        restrictions.append(
            "With chi = 0, p_r = 0, p_chi = 0 and spatially constant Q, "
            "the shift residual has mean equal to the current mean. "
            "No periodic p_Q removes that mean, and the source was left unchanged."
        )
    if diagnostics["held_out_momentum_max"] > TOL_INITIAL_CONSTRAINT:
        restrictions.append(
            "p_Q lies in the geometry subspace. Current outside that subspace "
            "remains in the full shift residual and is not filtered out."
        )
    if not momentum_info["antiderivative_accepted"]:
        restrictions.append(
            "The geometry derivative did not accept a periodic antiderivative "
            "of the mean-free projected current."
        )
    if not radius_info.get("converged"):
        restrictions.append(
            "The existing radius Newton did not converge, so the Hamilton constraint "
            "is not solved on this slice."
        )
    caller_unchanged = caller_bytes == tuple(
        np.ascontiguousarray(array).tobytes()
        for array in (
            state.Q, state.r, state.chi, state.p_Q, state.p_r, state.p_chi, state.phi0, state.phi1,
        )
    )
    if column_bytes is not None:
        caller_unchanged = caller_unchanged and column_bytes == (
            np.ascontiguousarray(phi0).tobytes(),
            np.ascontiguousarray(phi1).tobytes(),
        )
    report = {
        "ansatz": {
            "chi": 0.0,
            "p_r": 0.0,
            "p_chi": 0.0,
            "Q": q0,
            "Q_spatially_constant": True,
            "p_Q_mean_declared": 0.0,
        },
        "multiplicity": int(grid.fine.multiplicity),
        "kappa": int(grid.fine.kappa),
        "multiplicity_applied_once": True,
        "force_beta_plus_multiplicity_Pmom_max": _max_abs(
            source["force_beta"] + grid.fine.multiplicity * source["Pmom"]
        ),
        "occupations": regional_pair_occupations(grid.fine.occupations),
        "columns_unchanged": bool(
            _max_abs(data.phi0 - state.phi0) == 0.0 and _max_abs(data.phi1 - state.phi1) == 0.0
        ),
        "caller_state_unchanged": bool(caller_unchanged),
        "r_min_quadrature": float(np.min(fine_data.r)),
        "r_max_quadrature": float(np.max(fine_data.r)),
        "source_mean_subtracted": False,
        "evolution_performed": False,
        "current_mean": current_mean,
        "current_integral": float(current_mean * grid.length),
        "current_max": _max_abs(current_snapshot),
        "mean_compatible": bool(abs(current_mean) <= TOL_INITIAL_CONSTRAINT),
        "projected_momentum_mean_gap": mean_gap,
        "projected_shift_beyond_mean_max": _max_abs(beyond_mean),
        "momentum": momentum_info,
        "constraints": diagnostics,
        "hamilton_gap_against_zero_p_Q": hamilton_gap,
        "radius": {
            "converged": bool(radius_info.get("converged")),
            "blocker": radius_info.get("blocker"),
            "projected_residual_max": radius_info.get("residual_max"),
            "rho_independent_of_r_max": radius_info.get("rho_independent_of_r_max"),
            "r_min_coarse": None if not radius_info.get("converged") else float(np.min(data.r)),
            "r_max_coarse": None if not radius_info.get("converged") else float(np.max(data.r)),
            "imposed_radius": False,
            "filtered_v1_radius": False,
        },
        "positive_r": bool(np.min(fine_data.r) > 0.0),
        "positive_Q": bool(np.min(fine_data.Q) > 0.0),
        "chi_rate": rates,
        "restrictions": restrictions,
        "solves_projected_shift_except_mean": bool(
            momentum_info["antiderivative_accepted"]
            and _max_abs(beyond_mean) <= TOL_INITIAL_CONSTRAINT
        ),
        "solves_full_shift": bool(diagnostics["full_momentum_max"] <= TOL_INITIAL_CONSTRAINT),
        "solves_projected_hamilton": bool(
            radius_info.get("converged")
            and diagnostics["projected_hamilton_max"] <= TOL_INITIAL_CONSTRAINT
        ),
        "solves_full_hamilton": bool(diagnostics["full_hamilton_max"] <= TOL_INITIAL_CONSTRAINT),
    }
    return data, report


def finite_window_columns(grid):
    """Completed v1 lobe profiles, as fermion-grid columns.

    Quadrature samples use the embedding's half-density and spin frame.
    ``U_f†`` then keeps the antiperiodic band the Galerkin state can hold.
    """
    from .nsc_finite_window_embedding import (
        COLUMN_PHASE,
        K_CARRIER,
        MINUS_LEFT,
        MINUS_RIGHT,
        MINUS_RIGHT_PHASE,
        PLUS_LEFT,
        PLUS_RIGHT,
        PLUS_RIGHT_PHASE,
        _carrier,
        _lobe,
    )

    sqrt_dx = np.sqrt(grid.dx_q)
    phi0 = np.zeros((grid.nq, 6), dtype=np.complex128)
    phi1 = np.zeros((grid.nq, 6), dtype=np.complex128)
    for region in range(3):
        plus_env, plus_d = _lobe(grid.xi_q, region, PLUS_LEFT, PLUS_RIGHT, PLUS_RIGHT_PHASE)
        minus_env, minus_d = _lobe(grid.xi_q, region, MINUS_LEFT, MINUS_RIGHT, MINUS_RIGHT_PHASE)
        plus, _plus_d = _carrier(grid.xi_q, region, plus_env, plus_d, K_CARRIER, 0.0)
        minus, _minus_d = _carrier(
            grid.xi_q, region, minus_env, minus_d, -K_CARRIER, COLUMN_PHASE,
        )
        phi0[:, 2 * region] = sqrt_dx * plus / np.sqrt(2.0)
        phi1[:, 2 * region] = sqrt_dx * 1j * plus / np.sqrt(2.0)
        phi0[:, 2 * region + 1] = sqrt_dx * minus / np.sqrt(2.0)
        phi1[:, 2 * region + 1] = sqrt_dx * (-1j) * minus / np.sqrt(2.0)
    coarse0 = grid.U_f.conj().T @ phi0
    coarse1 = grid.U_f.conj().T @ phi1
    retained = float(np.linalg.norm(grid.U_f @ coarse0) / np.linalg.norm(phi0))
    return coarse0, coarse1, {
        "owner": "nsc_finite_window_embedding discrete lobe sample",
        "half_density": "sqrt(dx) on quadrature nodes, then U_f†",
        "band_retained_phi0": retained,
        "plus_right_phase": float(PLUS_RIGHT_PHASE),
        "minus_right_phase": float(MINUS_RIGHT_PHASE),
        "column_phase": float(COLUMN_PHASE),
        "carrier_k": float(K_CARRIER),
    }


def classical_embedding_momentum_integral(quad=256):
    """Occupation-weighted classical ∫ ψ†(-i∂)ψ on the v1 piecewise derivative.

    This is the continuum integral of the lobe profiles. It is not the
    spectral current used by ``source_from_columns``.
    """
    from .nsc_conformal_adm_source import sector_multiplicity
    from .nsc_finite_window_embedding import continuum_modes
    from .nsc_spherical_coupling import KAPPA, OCCUPATIONS

    packet = continuum_modes(int(quad))
    weights = packet["wt"]
    total = 0.0
    columns = []
    for mode, occupation in zip(packet["modes"], OCCUPATIONS):
        density = (-1j * np.conj(mode["f"]) * mode["df"]).real
        integral = float(np.sum(weights * density))
        total += float(occupation) * integral
        columns.append({
            "spin": mode["spin"],
            "region": int(mode["region"]),
            "occupation": float(occupation),
            "mass": float(np.sum(weights * np.abs(mode["f"]) ** 2)),
            "momentum_integral": integral,
        })
    multiplicity = int(sector_multiplicity(KAPPA))
    return {
        "quadrature_panels": int(quad),
        "multiplicity": multiplicity,
        "momentum_integral": float(total),
        "mean_of_minus_multiplicity_density": float(-multiplicity * total / PERIOD),
        "columns": columns,
        "pairwise_occupations": regional_pair_occupations(OCCUPATIONS),
    }


def _closed_form_controls(grid):
    q0 = gauge_ratio(grid)
    wavenumber = 2.0 * np.pi / grid.length
    sine = np.sin(wavenumber * grid.xi_q)
    mixed = 0.3 + sine
    mixed_snapshot = mixed.copy()
    momentum, info = shift_momentum(grid, mixed)
    closed = pull_geometry(grid, -np.cos(wavenumber * grid.xi_q) / (wavenumber * q0))
    residual = shift_residual(grid, momentum, mixed)
    shifted = momentum + 0.2
    residual_with_constant = shift_residual(grid, shifted, mixed)
    base = _fine_ansatz(grid, momentum)
    moved = _fine_ansatz(grid, shifted)
    hamilton_gap = _max_abs(
        hamilton_constraint(grid.fine, base) - hamilton_constraint(grid.fine, moved)
    )
    moved.p_chi = np.full(grid.nq, 0.01)
    base.p_chi = np.full(grid.nq, 0.01)
    conditional_gap = _max_abs(
        hamilton_constraint(grid.fine, base) - hamilton_constraint(grid.fine, moved)
    )
    _force_r, f_chi = partial_F(base.r, grid.fine.A, grid.fine.C_W)
    rate_base = geometric_rates(grid.fine, _fine_ansatz(grid, momentum))[2]
    rate_moved = geometric_rates(grid.fine, _fine_ansatz(grid, shifted))[2]
    rate_gap = _max_abs((rate_moved - rate_base) - grid.fine.length_density * 0.2 / (2.0 * float(f_chi)))
    pure, pure_info = shift_momentum(grid, np.ones(grid.nq))
    pure_residual = shift_residual(grid, pure, np.ones(grid.nq))
    return {
        "grid": {"nf": grid.nf, "ng": grid.ng, "nq": grid.nq},
        "mixed_current_unchanged": bool(_max_abs(mixed - mixed_snapshot) == 0.0),
        "source_mean_subtracted": info["source_mean_subtracted"],
        "closed_form_p_Q_gap": _max_abs(momentum - closed),
        "residual_max": _max_abs(residual),
        "residual_mean": float(np.mean(residual)),
        "residual_minus_mean_max": _max_abs(residual - np.mean(residual)),
        "constant_p_Q_shift_gap": _max_abs(residual_with_constant - residual),
        "constant_p_Q_hamilton_gap_at_zero_p_chi": hamilton_gap,
        "constant_p_Q_hamilton_gap_at_nonzero_p_chi": conditional_gap,
        "constant_p_Q_chi_rate_gap": rate_gap,
        "pure_mean_p_Q_max": pure_info["p_Q_max"],
        "pure_mean_residual_gap": _max_abs(pure_residual - 1.0),
        "antiderivative_accepted": info["antiderivative_accepted"],
    }


def _positive_carrier_columns(grid, mode=2.5):
    wavenumber = 2.0 * np.pi * float(mode) / grid.length
    phase = np.exp(1j * wavenumber * grid.xi_f) / np.sqrt(grid.nf)
    count = int(np.asarray(grid.fine.occupations).size)
    phi0 = np.zeros((grid.nf, count), dtype=np.complex128)
    phi1 = np.zeros((grid.nf, count), dtype=np.complex128)
    phi0[:, 0] = phase / np.sqrt(2.0)
    phi1[:, 0] = 1j * phase / np.sqrt(2.0)
    expected = (
        -float(grid.fine.multiplicity)
        * float(grid.fine.occupations[0])
        * wavenumber
        / grid.length
    )
    return phi0, phi1, {
        "mode": float(mode),
        "wavenumber": wavenumber,
        "expected_constant_current": expected,
    }


def _compact_initial(report, data):
    rates = report["chi_rate"] or {}
    velocities = rates.get("velocities") or {}
    constraints = report["constraints"]
    return {
        "current_mean": report["current_mean"],
        "current_integral": report["current_integral"],
        "current_max": report["current_max"],
        "mean_compatible": report["mean_compatible"],
        "source_mean_subtracted": report["source_mean_subtracted"],
        "p_Q_max": report["momentum"]["p_Q_max"],
        "p_Q_mean": report["momentum"]["p_Q_mean"],
        "projected_hamilton_max": constraints["projected_hamilton_max"],
        "full_hamilton_max": constraints["full_hamilton_max"],
        "held_out_hamilton_max": constraints["held_out_hamilton_max"],
        "projected_momentum_max": constraints["projected_momentum_max"],
        "full_momentum_max": constraints["full_momentum_max"],
        "held_out_momentum_max": constraints["held_out_momentum_max"],
        "projected_shift_beyond_mean_max": report["projected_shift_beyond_mean_max"],
        "projected_momentum_mean_gap": report["projected_momentum_mean_gap"],
        "hamilton_gap_against_zero_p_Q": report["hamilton_gap_against_zero_p_Q"],
        "solves_projected_shift_except_mean": report["solves_projected_shift_except_mean"],
        "solves_full_shift": report["solves_full_shift"],
        "solves_projected_hamilton": report["solves_projected_hamilton"],
        "solves_full_hamilton": report["solves_full_hamilton"],
        "positive_r": report["positive_r"],
        "positive_Q": report["positive_Q"],
        "r_min_quadrature": report["r_min_quadrature"],
        "r_max_quadrature": report["r_max_quadrature"],
        "Q_min_quadrature": constraints["Q_min_quadrature"],
        "hamilton_unresolved_fraction": constraints["hamilton_modes"]["unresolved_fraction"],
        "momentum_unresolved_fraction": constraints["momentum_modes"]["unresolved_fraction"],
        "caller_state_unchanged": report["caller_state_unchanged"],
        "radius_converged": report["radius"]["converged"],
        "radius_projected_residual_max": report["radius"]["projected_residual_max"],
        "rho_independent_of_r_max": report["radius"]["rho_independent_of_r_max"],
        "chart_proper_max": velocities.get("chart_proper_max"),
        "lifted_proper_max": velocities.get("lifted_proper_max"),
        "chi_rate_fine_max": rates.get("fine_max"),
        "chi_rate_gap": rates.get("gap_against_geometric_rates"),
        "chi_max": _max_abs(data.chi),
        "p_r_max": _max_abs(data.p_r),
        "p_chi_max": _max_abs(data.p_chi),
        "columns_unchanged": report["columns_unchanged"],
        "occupations_equal_within_regional_pairs": report["occupations"]["equal_within_regional_pairs"],
        "force_beta_plus_multiplicity_Pmom_max": report["force_beta_plus_multiplicity_Pmom_max"],
        "restrictions": report["restrictions"],
    }


def _without_grid(report):
    copied = dict(report)
    copied.pop("_grid", None)
    return copied


def build_record(fermions=64):
    """Exact shift controls plus rank6 and completed-embedding initial slices."""
    control_grid = build_grid(32)
    controls = _closed_form_controls(control_grid)
    carrier_grid = build_grid(32)
    phi0, phi1, carrier_meta = _positive_carrier_columns(carrier_grid)
    carrier_data, carrier_report = initial_cauchy_data(carrier_grid, phi0=phi0, phi1=phi1)
    carrier = _compact_initial(carrier_report, carrier_data)
    carrier["expected_constant_current"] = carrier_meta["expected_constant_current"]
    carrier["current_minus_expected"] = abs(
        carrier["current_mean"] - carrier_meta["expected_constant_current"]
    )

    grid = build_grid(int(fermions))
    rank_phi0, rank_phi1, _preparation, problems = load_physical_columns(grid.nf)
    rank_data, rank_report = initial_cauchy_data(grid, phi0=rank_phi0, phi1=rank_phi1)
    rank6 = _compact_initial(rank_report, rank_data)
    rank6["preparation_problems"] = problems

    embed_phi0, embed_phi1, embed_sampling = finite_window_columns(grid)
    embed_data, embed_report = initial_cauchy_data(grid, phi0=embed_phi0, phi1=embed_phi1)
    embedding = _compact_initial(embed_report, embed_data)
    embedding["sampling"] = embed_sampling

    finer = build_grid(int(fermions), quadrature=8 * int(fermions))
    fine_phi0, fine_phi1, fine_sampling = finite_window_columns(finer)
    fine_prepared = blank_state(finer, fine_phi0, fine_phi1)
    fine_prepared.Q = np.full(finer.ng, gauge_ratio(finer))
    fine_source = source_from_columns(finer.fine, prolong_state(finer, fine_prepared))
    fine_current = fine_source["force_beta"] / finer.dx_q
    fine_momentum, fine_info = shift_momentum(finer, fine_current)
    fine_residual = shift_residual(finer, fine_momentum, fine_current)
    fine_projected = pull_geometry(finer, fine_residual)
    embedding_finer_quadrature = {
        "nf": finer.nf,
        "nq": finer.nq,
        "radius_solved": False,
        "reason": "momentum diagnostic only; the radius owner is unchanged and was not rerun",
        "sampling": fine_sampling,
        "current_mean": fine_info["current_mean"],
        "current_max": fine_info["current_max"],
        "source_mean_subtracted": fine_info["source_mean_subtracted"],
        "p_Q_max": fine_info["p_Q_max"],
        "full_momentum_max": _max_abs(fine_residual),
        "projected_momentum_max": _max_abs(fine_projected),
        "projected_shift_beyond_mean_max": _max_abs(fine_projected - np.mean(fine_projected)),
        "antiderivative_accepted": fine_info["antiderivative_accepted"],
    }
    saved_v5 = json.loads(_V5_PATH.read_text())
    record = {
        "schema": SCHEMA,
        "status": "INITIAL_DATA_MEAN_AND_HELDOUT_RECORDED",
        "role": "initial momentum and radius data; not an evolution step or a source edit",
        "owners": {
            "shift_and_hamilton": "nsc_spherical_coupling.shift_constraint, hamilton_constraint",
            "source": "nsc_spherical_coupling.source_from_columns",
            "radius": "nsc_spherical_galerkin_coupling.solve_initial_radius",
            "chi_rate": "nsc_spherical_coupling.geometric_rates",
            "embedding_profiles": "nsc_finite_window_embedding",
        },
        "ansatz": {
            "chi": 0.0,
            "p_r": 0.0,
            "p_chi": 0.0,
            "Q": "b0/a0, spatially constant",
            "p_Q": "mean-zero geometry-band antiderivative of j/Q",
            "p_Q_mean_declared": 0.0,
            "mean_declared_zero_does_not_enter_C_or_D": True,
            "mean_declared_zero_fixes_chi_rate_zero_mode": True,
        },
        "tolerance_initial_constraint": TOL_INITIAL_CONSTRAINT,
        "classical_embedding_integral": classical_embedding_momentum_integral(256),
        "closed_form_controls": controls,
        "positive_carrier": carrier,
        "rank6_packet": rank6,
        "finite_window_embedding": embedding,
        "finite_window_embedding_nq_8nf_momentum": embedding_finer_quadrature,
        "saved_v5_rank6_witness": {
            "path": "lab/results/development/nsc-spherical-coupling-refinement-v5.json",
            "schema": saved_v5["schema"],
            "nq1024_current_max": saved_v5["nq1024"]["diagnostics"]["current_max"],
            "nq1024_full_momentum_max": saved_v5["nq1024"]["full_D_momentum_max"],
            "nq1024_chart_proper_velocity_max": saved_v5["nq1024"]["chart_proper_velocity_max"],
            "used_as": "saved near-zero current of the same rank6 packet; state not reloaded",
        },
        "field_equations_changed": False,
        "action_changed": False,
        "coefficients_changed": False,
    }
    payload = _jsonable(record)
    RECORD_PATH.parent.mkdir(parents=True, exist_ok=True)
    RECORD_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    payload["source_hashes"] = {
        "src/recursive_horizons/nsc_spherical_cauchy_data.py": _sha256(Path(__file__)),
        "tests/test_nsc_spherical_cauchy_data.py": _sha256(
            _LAB_ROOT / "tests" / "test_nsc_spherical_cauchy_data.py"
        ),
        "src/recursive_horizons/nsc_spherical_coupling.py": _sha256(
            _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_coupling.py"
        ),
        "src/recursive_horizons/nsc_spherical_galerkin_coupling.py": _sha256(
            _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_galerkin_coupling.py"
        ),
        "src/recursive_horizons/nsc_spherical_feedback_action.py": _sha256(
            _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_feedback_action.py"
        ),
        "src/recursive_horizons/nsc_finite_window_embedding.py": _sha256(
            _LAB_ROOT / "src" / "recursive_horizons" / "nsc_finite_window_embedding.py"
        ),
    }
    RECORD_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def _jsonable(value):
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


if __name__ == "__main__":
    build_record()
