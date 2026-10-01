"""Weak initial geometry residual from the spherical action.

The production evolution stays in ``nsc_spherical_galerkin_coupling``. This
module does not step the state, retune a tolerance, or replace the owner
residual. Under constant ``Q``, ``chi = p_r = p_chi = 0`` and ``r = y^2``,

    G = Q^2 r_mag^2 + Q rho / (8 pi A),
    Ry = -4 y'' + Q^2 y - G / y^3,
    C_action = -(8 pi A / Q) y^3 Ry,

with ``rho = force_L / dx``. ``C_action`` equals the owner lapse residual
``hamilton_constraint + rho`` when the spectral product rule holds. The
defect ``D^2(r^2) - 3 (Dr)^2 - 4 y^3 D^2 y`` is kept when it does not.

For ``G > 0`` and ``y > 0``, with ``G`` independent of ``y``,

    J(y) = ∫ [2 (y')^2 + Q^2 y^2 / 2 + G / (2 y^2)] dx

is strictly convex. Its Hessian is

    H(v, v) = 4 ||v'||^2 + ∫ (Q^2 + 3 G / y^4) v^2.

The dual residual of ``R(v) = ∫ v Ry`` controls the energy-norm distance to
a critical point only on the space over which the supremum is taken, and
only while ``Q^2 + 3 G / z^4`` stays positive on the segment. A geometry-band
weak number omits held-out modes. The quadrature grid omits the continuum
tail. No tail majorant is constructed here.

``assemble_record`` only measures. ``verify_saved`` reads a sealed record and
does not write it. ``write_record`` stores a fresh assembly only at an
explicit new path.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import time
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

from .nsc_spherical_cauchy_data import shift_momentum, shift_residual
from .nsc_spherical_coupling import (
    TOL_INITIAL_CONSTRAINT,
    _radius_residual,
    hamilton_constraint,
    magnetic_radius_square,
)
from .nsc_spherical_galerkin_coupling import (
    blank_state,
    build_grid,
    prolong_geometry,
    prolong_state,
    pull_geometry,
    source_from_columns,
    unresolved_fraction,
)

SCHEMA = "NSC-SPHERICAL-CAUCHY-WEAK-v1"
STATUS = "PARTIAL_TERMS_NO_TOTAL_BOUND"
CPU_BUDGET_S = 60.0
_LAB_ROOT = Path(__file__).resolve().parents[2]
_REPO_ROOT = _LAB_ROOT.parent
RECORD_PATH = _LAB_ROOT / "results" / "development" / "nsc-spherical-cauchy-weak-v1.json"
V5_JSON = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-refinement-v5.json"
V5_NPZ = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-refinement-v5.npz"
VOLATILE_EXECUTION_FIELDS = frozenset({
    "checkpoint_head",
    "cpu_seconds",
    "wall_seconds",
    "cpu_budget_exceeded",
})
SCIENTIFIC_SOURCE_PATHS = (
    "src/recursive_horizons/nsc_spherical_coupling.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_spherical_feedback_action.py",
    "src/recursive_horizons/nsc_spherical_cauchy_data.py",
)
GENERATOR_SOURCE_PATHS = (
    "src/recursive_horizons/nsc_spherical_cauchy_weak.py",
    "tests/test_nsc_spherical_cauchy_weak.py",
    "scripts/derive_nsc_spherical_cauchy_weak.py",
)
REPLAY_ABSOLUTE_TOLERANCE = 0.0
REPLAY_RELATIVE_TOLERANCE = 1e-12

STATEMENT = {
    "functional": "J(y)=integral[2 (y')^2 + Q^2 y^2/2 + G/(2 y^2)] dx",
    "euler_lagrange": "Ry=-4 y'' + Q^2 y - G/y^3",
    "owner_link": "C_tot=-(8 pi A/Q) y^3 Ry when r=y^2 and the spectral product rule holds",
    "algebra_without_product_rule": (
        "C_tot=(8 pi A/Q)[D^2(r^2)-3 (Dr)^2 - Q^2 r^2 + G] "
        "with G=Q^2 r_mag^2 + Q rho/(8 pi A)"
    ),
    "hessian": "H(v,v)=4 ||v'||^2 + integral (Q^2 + 3 G/y^4) v^2",
    "energy_norm": "||v||_E^2 = ||v'||^2 + ||v||^2",
    "dual": "||R||_*=sup |integral v Ry dx| / ||v||_E",
    "error_inequality": (
        "If y_* is critical in the tested space, G is independent of y, "
        "and c=Q^2+3G/z^4 stays at least c_seg>0 on the segment, then "
        "||y-y_*||_E <= ||R_y||_* / min(4, c_seg)."
    ),
    "conservative_mu": (
        "Declared G>0 and rho_independent is True give "
        "H>=4||v'||^2+Q^2||v||^2, so mu=min(4,Q^2) does not depend on y_*. "
        "Unknown rho independence does not establish this. "
        "A total bound still needs the full residual, continuum positivity of G, "
        "and a rounding enclosure."
    ),
    "polynomial_owner": (
        "P=2*r*r''-(r')^2-Q^2*r^2+G and C=(8*pi*A/Q)*P. "
        "Degree at most nf-1 from the antiperiodic column band and the odd radius band."
    ),
    "not_sufficient": (
        "A represented weak number does not control held-out modes, the "
        "continuum tail, the product-rule defect, quadrature, roundoff, or "
        "the shift current. This helper stores those partial terms and does "
        "not emit a total enclosure."
    ),
}


def _sha256(path):
    file_path = Path(path)
    if not file_path.is_file():
        return None
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _max_abs(values):
    array = np.asarray(values)
    if array.size == 0:
        return 0.0
    return float(np.max(np.abs(array)))


def _nodal(system, values):
    return float(system.dx * np.sum(np.asarray(values, dtype=float)))


def _head():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=_REPO_ROOT,
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def gauge_constants(system):
    """Positive constant Q and the audited magnetic radius used by the owner."""
    calibration = system.calibration
    q_value = float(calibration["b0"] / calibration["a0"])
    if not math.isfinite(q_value) or q_value <= 0.0:
        raise ValueError("gauge ratio b0/a0 is not a positive finite constant")
    area = float(system.A)
    if not math.isfinite(area) or area == 0.0:
        raise ValueError("area coefficient A is not a finite nonzero constant")
    rmag2 = float(magnetic_radius_square(system.coefficients))
    eight_pi_a = 8.0 * math.pi * area
    return {
        "Q": q_value,
        "A": area,
        "rmag2": rmag2,
        "eight_pi_A": eight_pi_a,
    }


def bracket_G(rho, constants):
    """G = Q^2 r_mag^2 + Q rho / (8 pi A), pointwise on the supplied rho."""
    rho = np.asarray(rho, dtype=float)
    return constants["Q"] ** 2 * constants["rmag2"] + constants["Q"] * rho / constants["eight_pi_A"]


def algebraic_total(system, radius, rho):
    """Owner residual rewritten with D(F), without the y product rule."""
    constants = gauge_constants(system)
    radius = np.asarray(radius, dtype=float)
    rho = np.asarray(rho, dtype=float)
    if radius.shape != rho.shape or radius.shape != (system.points,):
        raise ValueError("radius and rho must be quadrature-node arrays")
    derivative = system.derivative
    radius_derivative = derivative @ radius
    radius_square_second = derivative @ (derivative @ (radius * radius))
    G = bracket_G(rho, constants)
    q_value = constants["Q"]
    total = (constants["eight_pi_A"] / q_value) * (
        radius_square_second - 3.0 * radius_derivative ** 2 - q_value ** 2 * radius ** 2 + G
    )
    return total, G


def owner_total(system, radius, rho):
    """Existing ``_radius_residual``: hamilton_constraint + rho on the ansatz."""
    return np.asarray(_radius_residual(system, np.asarray(radius, dtype=float), rho), dtype=float)


def action_ry(system, y, G, constants=None):
    """Ry = -4 y'' + Q^2 y - G/y^3 on the spectral derivative."""
    if constants is None:
        constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    G = np.asarray(G, dtype=float)
    if np.min(y) <= 0.0:
        raise ValueError("action residual requires y > 0")
    second = system.derivative @ (system.derivative @ y)
    return -4.0 * second + constants["Q"] ** 2 * y - G / y ** 3


def action_total(y, Ry, constants):
    """C_action = -(8 pi A / Q) y^3 Ry."""
    y = np.asarray(y, dtype=float)
    return -(constants["eight_pi_A"] / constants["Q"]) * y ** 3 * np.asarray(Ry, dtype=float)


def product_rule_defect(system, y):
    """D^2(r^2) - 3 (Dr)^2 - 4 y^3 D^2 y, with r = y^2. Zero when the chain rule holds."""
    y = np.asarray(y, dtype=float)
    derivative = system.derivative
    radius = y * y
    radius_derivative = derivative @ radius
    radius_square_second = derivative @ (derivative @ (radius * radius))
    y_second = derivative @ (derivative @ y)
    return radius_square_second - 3.0 * radius_derivative ** 2 - 4.0 * y ** 3 * y_second


def action_functional(system, y, G, constants=None):
    """Nodal integral of 2 (y')^2 + Q^2 y^2/2 + G/(2 y^2)."""
    if constants is None:
        constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    if np.min(y) <= 0.0:
        raise ValueError("action requires y > 0")
    slope = system.derivative @ y
    G = np.asarray(G, dtype=float)
    return _nodal(system, 2.0 * slope ** 2 + 0.5 * constants["Q"] ** 2 * y ** 2 + G / (2.0 * y ** 2))


def weak_pairing(system, y, G, test, constants=None):
    """First variation R(test) = ∫ [4 y' test' + Q^2 y test - G test / y^3]."""
    if constants is None:
        constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    test = np.asarray(test, dtype=float)
    G = np.asarray(G, dtype=float)
    derivative = system.derivative
    return _nodal(
        system,
        4.0 * (derivative @ y) * (derivative @ test)
        + constants["Q"] ** 2 * y * test
        - G * test / y ** 3,
    )


def hessian_quadratic(system, y, G, test, constants=None):
    """H(test, test) = 4 ||test'||^2 + ∫ (Q^2 + 3 G/y^4) test^2."""
    if constants is None:
        constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    if np.min(np.abs(y)) == 0.0:
        raise ValueError("hessian requires y ≠ 0")
    slope = system.derivative @ np.asarray(test, dtype=float)
    coefficient = constants["Q"] ** 2 + 3.0 * np.asarray(G, dtype=float) / y ** 4
    return _nodal(system, 4.0 * slope ** 2 + coefficient * np.asarray(test, dtype=float) ** 2)


def energy_norm(system, test):
    slope = system.derivative @ np.asarray(test, dtype=float)
    return math.sqrt(_nodal(system, slope ** 2 + np.asarray(test, dtype=float) ** 2))


def coefficient_c(y, G, constants):
    return constants["Q"] ** 2 + 3.0 * np.asarray(G, dtype=float) / np.asarray(y, dtype=float) ** 4


def _mu_fields(constants, *, positive_y, positive_G, rho_independent):
    """Uniform mu is conditional. It is not a total error bound."""
    declared = rho_independent is True
    available = bool(positive_y and positive_G and declared)
    mu = min(4.0, constants["Q"] ** 2)
    return {
        "uniform_mu": float(mu) if available else None,
        "uniform_mu_requires_y_star": False,
        "uniform_mu_declared": available,
        "total_bound": None,
        "total_bound_gaps": [
            "full residual dual",
            "continuum positivity of G",
            "rounding enclosure",
        ],
        "continuum_positivity_enclosed": False,
        "rounding_enclosed": False,
    }


def coercivity_witness(system, y, G, *, rho_independent=None):
    """Pointwise sign of G and of c = Q^2 + 3 G/y^4.

    ``rho_independent is True`` is required. An omitted flag is not evidence
    that G is independent of y. Declared positive G then gives
    ``mu = min(4, Q^2)`` without ``y_*``. The sharper segment constant is
    optional, and this witness still emits no total bound.
    """
    constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    G = np.asarray(G, dtype=float)
    positive_y = bool(np.min(y) > 0.0)
    positive_G = bool(np.min(G) > 0.0)
    declared = rho_independent is True
    reason = (
        "Declared G>0 and rho_independent is True give uniform mu=min(4,Q^2) "
        "without y_*. Unknown rho independence does not establish that hypothesis. "
        "No total bound is emitted: the full residual, continuum positivity of G, "
        "and a rounding enclosure are still open."
    )
    if not positive_y:
        report = {
            "positive_y": False,
            "positive_G": positive_G,
            "c_min": None,
            "c_max": None,
            "frozen_alpha": None,
            "frozen_alpha_is_not_the_segment_constant": True,
            "strict_convexity_assumption": False,
            "error_inequality_applicable": False,
            "error_inequality_reason": reason,
            "rho_independent_of_y": rho_independent,
            "blocker": "y is not positive",
        }
        report.update(_mu_fields(
            constants, positive_y=False, positive_G=positive_G, rho_independent=rho_independent,
        ))
        return report
    values = coefficient_c(y, G, constants)
    c_min = float(np.min(values))
    c_max = float(np.max(values))
    assumption = bool(positive_G and c_min > 0.0 and declared)
    alpha = min(4.0, c_min) if c_min > 0.0 else None
    blockers = []
    if not positive_G:
        blockers.append("G is not positive")
    if not declared:
        blockers.append(
            "rho independence of y was not established"
            if rho_independent is None else "rho depends on y"
        )
    if positive_G and declared and c_min <= 0.0:
        blockers.append("nodal c_min is not positive")
    report = {
        "positive_y": True,
        "positive_G": positive_G,
        "G_min": float(np.min(G)),
        "G_max": float(np.max(G)),
        "c_min": c_min,
        "c_max": c_max,
        "frozen_alpha": None if alpha is None else float(alpha),
        "frozen_alpha_is_not_the_segment_constant": True,
        "strict_convexity_assumption": assumption,
        "error_inequality_applicable": False,
        "error_inequality_reason": reason,
        "rho_independent_of_y": rho_independent,
        "blocker": None if assumption else "; ".join(blockers),
    }
    report.update(_mu_fields(
        constants, positive_y=True, positive_G=positive_G, rho_independent=rho_independent,
    ))
    return report


def dual_norm_E_parts(system, density, mode_limit):
    """Fourier energy dual. Nyquist symbol matches ``periodic_derivative`` (k = 0)."""
    density = np.asarray(density, dtype=float)
    count = int(system.points)
    if density.shape != (count,):
        raise ValueError("dual residual expects one quadrature array")
    modes = np.fft.fftfreq(count) * count
    wavenumber = 2.0 * math.pi * modes / float(system.length)
    nyquist = np.abs(modes) == (count / 2.0)
    wavenumber = np.array(wavenumber, dtype=float, copy=True)
    wavenumber[nyquist] = 0.0
    transformed = np.fft.fft(density)
    weight = (np.abs(transformed) ** 2) / (wavenumber ** 2 + 1.0)
    resolved = np.abs(modes) <= float(mode_limit) + 1e-9
    resolved[nyquist] = False
    factor = float(system.dx) / count
    resolved_sq = factor * float(np.sum(weight[resolved]))
    unresolved_sq = factor * float(np.sum(weight[~resolved]))

    def _sqrt(value):
        if value < 0.0 and value > -1e-18:
            value = 0.0
        return float(math.sqrt(value))

    return {
        "represented": _sqrt(resolved_sq),
        "unresolved": _sqrt(unresolved_sq),
        "full": _sqrt(resolved_sq + unresolved_sq),
    }


def represented_dual_matrix(grid, density):
    """Energy dual over prolongations of the odd geometry nodes."""
    derivative_on_columns = grid.derivative @ grid.A_g
    gram = float(grid.dx_q) * (derivative_on_columns.T @ derivative_on_columns)
    gram = gram + float(grid.dx_g) * np.eye(grid.ng)
    gram = 0.5 * (gram + gram.T)
    load = float(grid.dx_g) * pull_geometry(grid, np.asarray(density, dtype=float))
    try:
        solved = np.linalg.solve(gram, load)
    except np.linalg.LinAlgError as error:
        return {
            "dual_E": None,
            "pull_max": _max_abs(load / grid.dx_g),
            "solver": str(error),
        }
    value = float(load @ solved)
    if value < 0.0 and value > -1e-18:
        value = 0.0
    if value < 0.0:
        return {
            "dual_E": None,
            "pull_max": _max_abs(load / grid.dx_g),
            "solver": "energy Gram lost positivity",
            "gram_quadratic": value,
        }
    return {
        "dual_E": float(math.sqrt(value)),
        "pull_max": _max_abs(pull_geometry(grid, density)),
        "solver": None,
    }


def strong_split(grid, values):
    projected = pull_geometry(grid, values)
    held = np.asarray(values, dtype=float) - prolong_geometry(grid, projected)
    modes = unresolved_fraction(grid, values)
    transformed = np.fft.fft(np.asarray(values, dtype=float)) / values.size
    return {
        "full_max": _max_abs(values),
        "represented_pull_max": _max_abs(projected),
        "unresolved_max": _max_abs(held),
        "unresolved_power": modes["unresolved_power"],
        "unresolved_fraction": modes["unresolved_fraction"],
        "unresolved_max_mode": modes["unresolved_max_mode"],
        "highest_mode_amplitude": float(abs(transformed[values.size // 2])),
        "tail_bound": None,
    }


def observe_discrete_inequality(system, y, y_star, G, c_seg):
    """Check the energy-norm inequality on one supplied discrete segment.

    ``c_seg`` is an input. This function does not estimate it from samples
    and does not promote the result to a continuum certificate.
    """
    constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    y_star = np.asarray(y_star, dtype=float)
    G = np.asarray(G, dtype=float)
    if not math.isfinite(c_seg) or c_seg <= 0.0 or np.min(y) <= 0.0 or np.min(y_star) <= 0.0:
        return {
            "holds": False,
            "upper_bound": None,
            "reason": "c_seg > 0 and positive y, y_* are required",
            "is_continuum_certificate": False,
            "total_certified": False,
        }
    critical = action_ry(system, y_star, G, constants)
    residual = action_ry(system, y, G, constants)
    dual = dual_norm_E_parts(system, residual, mode_limit=system.points)["full"]
    alpha = min(4.0, float(c_seg))
    error_norm = energy_norm(system, y - y_star)
    upper = dual / alpha
    critical_gap = _max_abs(critical)
    holds = bool(error_norm <= upper + 1e-9 and critical_gap <= 1e-8)
    return {
        "holds": holds,
        "error_norm_E": error_norm,
        "dual_E": dual,
        "alpha": alpha,
        "c_seg_supplied": float(c_seg),
        "upper_bound": upper,
        "critical_gap": critical_gap,
        "critical_gap_gate": 1e-8,
        "critical_gap_gate_is_stationarity_of_this_example": True,
        "is_continuum_certificate": False,
        "total_certified": False,
        "reason": None if holds else "inequality or stationarity of the supplied critical point failed",
    }


def compare_owner_and_action(system, y, rho):
    """Sign check: algebraic form versus owner, and y-form versus that defect."""
    constants = gauge_constants(system)
    y = np.asarray(y, dtype=float)
    radius = y * y
    rho = np.asarray(rho, dtype=float)
    owner = owner_total(system, radius, rho)
    algebraic, G = algebraic_total(system, radius, rho)
    Ry = action_ry(system, y, G, constants)
    claimed = action_total(y, Ry, constants)
    defect = product_rule_defect(system, y)
    scaled_defect = (constants["eight_pi_A"] / constants["Q"]) * defect
    return {
        "algebraic_gap_max": _max_abs(owner - algebraic),
        "action_gap_max": _max_abs(owner - claimed),
        "scaled_product_rule_defect_max": _max_abs(scaled_defect),
        "defect_explains_action_gap_max": _max_abs((owner - claimed) - scaled_defect),
        "G_min": float(np.min(G)),
        "G_max": float(np.max(G)),
        "owner_max": _max_abs(owner),
        "action_max": _max_abs(claimed),
    }


def constructor_current(grid, state):
    """Mean and held-out shift residual from the existing Cauchy constructor.

    The radius Newton inside ``initial_cauchy_data`` is not called. The source
    array passed to ``shift_momentum`` is copied and checked unchanged.
    """
    fine = prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    current = np.array(source["force_beta"] / grid.dx_q, dtype=float, copy=True)
    snapshot = current.copy()
    momentum, info = shift_momentum(grid, current)
    if _max_abs(current - snapshot) != 0.0:
        raise RuntimeError("shift momentum altered the current")
    residual = shift_residual(grid, momentum, current)
    split = strong_split(grid, residual)
    projected = pull_geometry(grid, residual)
    projected_mean = float(np.mean(projected))
    fine_momentum = prolong_state(grid, _with_momentum(state, momentum))
    zero_momentum = prolong_state(grid, _with_momentum(state, np.zeros(grid.ng)))
    hamilton_gap = _max_abs(
        hamilton_constraint(grid.fine, fine_momentum) - hamilton_constraint(grid.fine, zero_momentum)
    )
    current_mean = float(info["current_mean"])
    return {
        "radius_newton_rerun": False,
        "source_mean_subtracted": False,
        "current_unchanged": True,
        "current_mean": current_mean,
        "current_integral": float(info["current_integral"]),
        "current_max": float(info["current_max"]),
        "mean_compatible_at_constructor_tolerance": bool(abs(current_mean) <= TOL_INITIAL_CONSTRAINT),
        "constructor_tolerance": TOL_INITIAL_CONSTRAINT,
        "tolerance_is_not_a_proof": True,
        "p_Q_max": float(info["p_Q_max"]),
        "p_Q_mean": float(info["p_Q_mean"]),
        "antiderivative_accepted": bool(info["antiderivative_accepted"]),
        "full_momentum_max": split["full_max"],
        "represented_momentum_max": split["represented_pull_max"],
        "held_out_momentum_max": split["unresolved_max"],
        "projected_momentum_mean": projected_mean,
        "projected_momentum_mean_gap": abs(projected_mean - current_mean),
        "projected_shift_beyond_mean_max": _max_abs(projected - projected_mean),
        "hamilton_gap_against_zero_p_Q": hamilton_gap,
        "multiplicity_applied_once": source["multiplicity_applied_once"] is True,
    }


def _with_momentum(state, momentum):
    updated = state.copy()
    updated.p_Q = np.asarray(momentum, dtype=float)
    updated.chi = np.zeros_like(updated.chi)
    updated.p_r = np.zeros_like(updated.p_r)
    updated.p_chi = np.zeros_like(updated.p_chi)
    return updated


def _ansatz_report(state, constants):
    return {
        "Q_minus_gauge_max": _max_abs(state.Q - constants["Q"]),
        "chi_max": _max_abs(state.chi),
        "p_r_max": _max_abs(state.p_r),
        "p_chi_max": _max_abs(state.p_chi),
        "p_Q_max": _max_abs(state.p_Q),
        "holds": bool(
            _max_abs(state.Q - constants["Q"]) <= 1e-12
            and _max_abs(state.chi) <= 1e-12
            and _max_abs(state.p_r) <= 1e-12
            and _max_abs(state.p_chi) <= 1e-12
        ),
    }


def assess_radius(grid, radius, rho, *, rho_independent=None):
    """Separate represented, unresolved, quadrature-ready, and roundoff pieces.

    ``radius`` is the coarse geometry array. ``rho`` is the quadrature array
    used as the owner source. No total bound is returned.
    """
    constants = gauge_constants(grid.fine)
    radius = np.asarray(radius, dtype=float)
    rho = np.asarray(rho, dtype=float)
    if radius.shape != (grid.ng,) or rho.shape != (grid.nq,):
        raise ValueError("radius is coarse and rho is on the quadrature nodes")
    fine_radius = prolong_geometry(grid, radius)
    if np.min(fine_radius) <= 0.0:
        raise ValueError("prolonged radius left the positive chart")
    y = np.sqrt(fine_radius)
    owner = owner_total(grid.fine, fine_radius, rho)
    algebraic, G = algebraic_total(grid.fine, fine_radius, rho)
    second_radius_square = grid.fine.derivative @ (grid.fine.derivative @ (fine_radius * fine_radius))
    scale = constants["eight_pi_A"] / constants["Q"]
    largest_term = _max_abs(second_radius_square)
    bracket_max = _max_abs(algebraic) / abs(scale)
    cancellation = largest_term / max(bracket_max, 1e-300)
    unit_roundoff = float(np.finfo(float).eps) * cancellation * abs(scale)
    Ry = action_ry(grid.fine, y, G, constants)
    claimed = action_total(y, Ry, constants)
    defect = product_rule_defect(grid.fine, y)
    scaled_defect = (constants["eight_pi_A"] / constants["Q"]) * defect
    owner_split = strong_split(grid, owner)
    action_split = strong_split(grid, Ry)
    fourier = dual_norm_E_parts(grid.fine, Ry, grid.ng // 2)
    matrix = represented_dual_matrix(grid, Ry)
    witness = coercivity_witness(grid.fine, y, G, rho_independent=rho_independent)
    matrix_dual = matrix["dual_E"]
    fourier_gap = None if matrix_dual is None else abs(matrix_dual - fourier["represented"])
    return {
        "total_certified": False,
        "represented_weak_certifies_total": False,
        "geometry_error_upper_bound": None,
        "tail_bound": None,
        "tail_bound_reason": (
            "no decay hypothesis was proved for this residual; "
            "highest-mode amplitude and held-out power are grid samples only"
        ),
        "Q": constants["Q"],
        "A": constants["A"],
        "rmag2": constants["rmag2"],
        "y_min": float(np.min(y)),
        "y_max": float(np.max(y)),
        "r_min": float(np.min(fine_radius)),
        "r_max": float(np.max(fine_radius)),
        "coercivity": witness,
        "roundoff": {
            "algebraic_identity_gap_max": _max_abs(owner - algebraic),
            "defect_explains_action_gap_max": _max_abs((owner - claimed) - scaled_defect),
            "largest_D2_r2_term": largest_term,
            "cancellation_condition": cancellation,
            "unit_roundoff_times_condition": unit_roundoff,
            "gap_is_not_a_new_acceptance_line": True,
        },
        "product_rule": {
            "scaled_defect_max": _max_abs(scaled_defect),
            "action_minus_owner_max": _max_abs(claimed - owner),
        },
        "represented": {
            "owner_pull_max": owner_split["represented_pull_max"],
            "action_pull_max": action_split["represented_pull_max"],
            "weak_dual_E_matrix": matrix_dual,
            "weak_dual_E_fourier": fourier["represented"],
            "matrix_fourier_gap": fourier_gap,
            "matrix_solver": matrix["solver"],
        },
        "unresolved": {
            "owner_max": owner_split["unresolved_max"],
            "action_ry_max": action_split["unresolved_max"],
            "weak_dual_E": fourier["unresolved"],
            "power": owner_split["unresolved_power"],
            "fraction": owner_split["unresolved_fraction"],
            "max_mode": owner_split["unresolved_max_mode"],
            "highest_mode_amplitude": owner_split["highest_mode_amplitude"],
        },
        "full_strong": {
            "owner_max": owner_split["full_max"],
            "owner_l2": math.sqrt(_nodal(grid.fine, owner ** 2)),
            "action_max": _max_abs(claimed),
            "ry_max": action_split["full_max"],
            "weak_dual_E": fourier["full"],
        },
        "partial_terms": {
            "represented_weak_dual_E": matrix_dual,
            "represented_owner_pull_max": owner_split["represented_pull_max"],
            "unresolved_owner_max": owner_split["unresolved_max"],
            "unresolved_weak_dual_E": fourier["unresolved"],
            "full_strong_owner_max": owner_split["full_max"],
            "full_weak_dual_E": fourier["full"],
            "product_rule_scaled_defect_max": _max_abs(scaled_defect),
            "algebraic_roundoff_gap_max": _max_abs(owner - algebraic),
            "highest_mode_amplitude": owner_split["highest_mode_amplitude"],
            "tail_bound": None,
        },
    }


def _source_rho(grid, state):
    fine = prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    rho = np.asarray(source["force_L"] / grid.dx_q, dtype=float)
    varied = state.copy()
    varied.r = state.r * 1.7
    other = source_from_columns(grid.fine, prolong_state(grid, varied))
    variation = _max_abs(other["force_L"] - source["force_L"])
    return rho, variation, source


def quadrature_strong(grid, state):
    """Owner strong residual of one coarse state on this grid's source."""
    rho, variation, _source = _source_rho(grid, state)
    owner = owner_total(grid.fine, prolong_geometry(grid, state.r), rho)
    split = strong_split(grid, owner)
    return {
        "nq": int(grid.nq),
        "full_strong_owner_max": split["full_max"],
        "represented_owner_pull_max": split["represented_pull_max"],
        "unresolved_owner_max": split["unresolved_max"],
        "rho_independent_of_r_max": variation,
        "tail_bound": None,
    }


def _load_seed(payload, prefix):
    def take(name):
        key = prefix + name
        if key not in payload:
            raise KeyError(key)
        return payload[key]

    radius = np.asarray(take("radius"), dtype=float)
    rho = np.asarray(take("rho"), dtype=float)
    phi0 = np.asarray(take("columns_phi0"))
    phi1 = np.asarray(take("columns_phi1"))
    quadrature_radius = np.asarray(take("quadrature_radius"), dtype=float)
    geometry_Q = np.asarray(take("geometry_Q"), dtype=float)
    return {
        "nf": int(take("nf")),
        "ng": int(take("ng")),
        "nq": int(take("nq")),
        "radius": radius,
        "rho": rho,
        "phi0": phi0,
        "phi1": phi1,
        "quadrature_radius": quadrature_radius,
        "geometry_Q": geometry_Q,
        "chi_max": _max_abs(take("geometry_chi")),
        "p_r_max": _max_abs(take("geometry_p_r")),
        "p_chi_max": _max_abs(take("geometry_p_chi")),
        "p_Q_max": _max_abs(take("geometry_p_Q")),
        "stored_projected_residual": float(take("projected_residual")),
        "stored_full_operator": float(take("full_operator")),
    }


def fft_periodic_derivative(values, length, order):
    """Fourier derivative. Nyquist symbol is 0, matching ``periodic_derivative``."""
    values = np.asarray(values, dtype=float)
    count = values.size
    modes = np.fft.fftfreq(count) * count
    wavenumber = 2.0 * np.pi * modes / float(length)
    wavenumber = np.array(wavenumber, copy=True)
    wavenumber[np.abs(modes) == (count / 2.0)] = 0.0
    factor = (1j * wavenumber) ** int(order)
    return np.fft.ifft(factor * np.fft.fft(values)).real


def resample_trigonometric(values, count_out):
    """Zero-pad a periodic trigonometric polynomial. The Nyquist sample is split."""
    values = np.asarray(values, dtype=float)
    count = values.size
    if count_out == count:
        return values.copy()
    if count_out < count or count_out % 2 or count % 2:
        raise ValueError("resample expects a finer even periodic grid")
    spectrum = np.fft.fft(values)
    output = np.zeros(count_out, dtype=complex)
    half = count // 2
    output[:half] = spectrum[:half]
    output[half] = 0.5 * spectrum[half]
    output[-half] = 0.5 * spectrum[half]
    output[-half + 1 :] = spectrum[half + 1 :]
    return np.fft.ifft(output * (count_out / count)).real


def alias_degree_guard(nf, ng, nq):
    """Degree cap from the actual column band and the odd radius band.

    Antiperiodic modes stop at ``nf/2 - 1/2``. Products, hence ``G``, have
    integer degree at most ``nf - 1``. The radius degree is at most ``ng//2``,
    so ``P`` has degree at most ``nf - 1``. The samples are alias-free when
    the quadrature Nyquist is strictly above that degree.
    """
    degree_r = int(ng) // 2
    degree_source = int(nf) - 1
    degree_p = max(2 * degree_r, degree_source)
    nyquist = int(nq) // 2
    return {
        "degree_r_band": degree_r,
        "degree_source_band": degree_source,
        "degree_P": int(degree_p),
        "nyquist": nyquist,
        "alias_free": bool(nyquist > degree_p),
    }


def _mode_summary(values, band_limit, predicted_degree):
    values = np.asarray(values, dtype=float)
    count = values.size
    coefficients = np.fft.fft(values) / count
    modes = np.rint(np.fft.fftfreq(count) * count).astype(int)
    amplitude = np.abs(coefficients)
    power = amplitude ** 2

    def highest(threshold):
        occupied = amplitude > threshold
        if not np.any(occupied):
            return 0
        return int(np.max(np.abs(modes[occupied])))

    represented = (np.abs(modes) <= int(band_limit)) & (np.abs(modes) != (count // 2))
    above = np.abs(modes) > int(predicted_degree)
    return {
        "highest_above_1e-8": highest(1e-8),
        "highest_above_1e-12": highest(1e-12),
        "represented_max_coeff": float(np.max(amplitude[represented])) if np.any(represented) else 0.0,
        "heldout_max_coeff": float(np.max(amplitude[~represented])) if np.any(~represented) else 0.0,
        "power_above_predicted_degree": float(np.sum(power[above])),
        "nyquist_coeff": float(amplitude[count // 2]),
    }


def polynomial_owner_residual(radius, rho, constants, length):
    """``P = 2 r r'' - (r')^2 - Q^2 r^2 + G`` and ``C = (8 pi A / Q) P``.

    Derivatives are Fourier derivatives of ``r``. This does not replace ``rho``.
    """
    radius = np.asarray(radius, dtype=float)
    slope = fft_periodic_derivative(radius, length, 1)
    second = fft_periodic_derivative(radius, length, 2)
    G = bracket_G(rho, constants)
    numerator = (
        2.0 * radius * second
        - slope ** 2
        - constants["Q"] ** 2 * radius ** 2
        + G
    )
    total = (constants["eight_pi_A"] / constants["Q"]) * numerator
    return total, G, numerator


def polynomial_case(grid, fine_radius, rho, constants):
    """Dense owner residual against the Fourier product form. Scalars only."""
    rho = np.asarray(rho, dtype=float)
    dense = owner_total(grid.fine, fine_radius, rho)
    fft_total, G, numerator = polynomial_owner_residual(
        fine_radius, rho, constants, grid.length,
    )
    guard = alias_degree_guard(grid.nf, grid.ng, grid.nq)
    band = grid.ng // 2
    return {
        "scalars": {
            "nf": int(grid.nf),
            "ng": int(grid.ng),
            "nq": int(grid.nq),
            "dense_max": _max_abs(dense),
            "fft_product_max": _max_abs(fft_total),
            "dense_minus_fft_max": _max_abs(dense - fft_total),
            "G_min": float(np.min(G)),
            "G_max": float(np.max(G)),
            "alias_guard": guard,
            "degree_r": _mode_summary(fine_radius, band, guard["degree_r_band"]),
            "degree_G": _mode_summary(G, band, guard["degree_source_band"]),
            "degree_numerator": _mode_summary(numerator, band, guard["degree_P"]),
            "degree_dense_C": _mode_summary(dense, band, guard["degree_P"]),
            "source_rho_modified": False,
            "current_modified": False,
            "rounding_enclosure": None,
            "conditioning_indicator": True,
        },
        "fft_total": fft_total,
        "dense": dense,
    }


def assess_saved_seed(seed, *, doubled_budget_s, started):
    """Replay one NPZ seed with the current source. No Newton and no evolution."""
    grid = build_grid(seed["nf"], quadrature=seed["nq"])
    if (grid.nf, grid.ng, grid.nq) != (seed["nf"], seed["ng"], seed["nq"]):
        raise RuntimeError(f"grid {(grid.nf, grid.ng, grid.nq)} does not match the seed")
    state = blank_state(grid, seed["phi0"], seed["phi1"])
    state.r = seed["radius"].copy()
    state.Q = seed["geometry_Q"].copy()
    constants = gauge_constants(grid.fine)
    column_bytes = (
        np.ascontiguousarray(state.phi0).tobytes(),
        np.ascontiguousarray(state.phi1).tobytes(),
    )
    source_rho, variation, _source = _source_rho(grid, state)
    frozen = assess_radius(grid, state.r, seed["rho"], rho_independent=variation <= 1e-8)
    fresh = assess_radius(grid, state.r, source_rho, rho_independent=variation <= 1e-8)
    current = constructor_current(grid, state)
    prolonged = prolong_geometry(grid, state.r)
    baseline_polynomial = polynomial_case(grid, prolonged, seed["rho"], constants)
    quadrature = {
        "computed": False,
        "reason": None,
        "stored_nq": None,
        "doubled_nq": None,
        "difference_of_full_strong_maxima": None,
        "tail_bound": None,
    }
    polynomial = {
        "baseline_stored_rho": baseline_polynomial["scalars"],
        "doubled_resampled_rho": None,
        "doubled_recomputed_source": None,
        "resampled_minus_recomputed_rho_max": None,
        "baseline_polynomial_on_fine_max": None,
        "fine_fft_minus_resampled_baseline_max": None,
        "rounding_enclosure": None,
    }
    spent = time.process_time() - started
    if spent > doubled_budget_s:
        quadrature["reason"] = (
            f"cpu {spent} already past the doubled-grid allowance {doubled_budget_s}"
        )
    else:
        doubled = build_grid(seed["nf"], quadrature=2 * seed["nq"])
        doubled_state = blank_state(doubled, seed["phi0"], seed["phi1"])
        doubled_state.r = seed["radius"].copy()
        doubled_state.Q = seed["geometry_Q"].copy()
        stored_view = {
            "nq": int(grid.nq),
            "full_strong_owner_max": fresh["full_strong"]["owner_max"],
            "represented_owner_pull_max": fresh["represented"]["owner_pull_max"],
            "unresolved_owner_max": fresh["unresolved"]["owner_max"],
            "rho_independent_of_r_max": variation,
            "tail_bound": None,
        }
        doubled_rho, doubled_variation, _doubled_source = _source_rho(doubled, doubled_state)
        doubled_radius = prolong_geometry(doubled, doubled_state.r)
        doubled_owner = owner_total(doubled.fine, doubled_radius, doubled_rho)
        doubled_split = strong_split(doubled, doubled_owner)
        doubled_view = {
            "nq": int(doubled.nq),
            "full_strong_owner_max": doubled_split["full_max"],
            "represented_owner_pull_max": doubled_split["represented_pull_max"],
            "unresolved_owner_max": doubled_split["unresolved_max"],
            "rho_independent_of_r_max": doubled_variation,
            "tail_bound": None,
        }
        recomputed = polynomial_case(doubled, doubled_radius, doubled_rho, gauge_constants(doubled.fine))
        resampled_rho = resample_trigonometric(source_rho, doubled.nq)
        resampled = polynomial_case(
            doubled, doubled_radius, resampled_rho, gauge_constants(doubled.fine),
        )
        resampled_baseline = resample_trigonometric(baseline_polynomial["fft_total"], doubled.nq)
        quadrature = {
            "computed": True,
            "reason": "same coarse seed on nq and 2 nq; source recomputed, not prolonged",
            "stored_nq": stored_view,
            "doubled_nq": doubled_view,
            "difference_of_full_strong_maxima": float(
                doubled_view["full_strong_owner_max"] - stored_view["full_strong_owner_max"]
            ),
            "tail_bound": None,
            "difference_is_not_a_tail_bound": True,
        }
        polynomial = {
            "baseline_stored_rho": baseline_polynomial["scalars"],
            "doubled_resampled_rho": resampled["scalars"],
            "doubled_recomputed_source": recomputed["scalars"],
            "resampled_minus_recomputed_rho_max": _max_abs(resampled_rho - doubled_rho),
            "baseline_polynomial_on_fine_max": _max_abs(resampled_baseline),
            "fine_fft_minus_resampled_baseline_max": _max_abs(
                resampled["fft_total"] - resampled_baseline
            ),
            "rounding_enclosure": None,
            "source_rho_modified": False,
            "current_modified": False,
        }
    unchanged = column_bytes == (
        np.ascontiguousarray(state.phi0).tobytes(),
        np.ascontiguousarray(state.phi1).tobytes(),
    )
    return {
        "nf": seed["nf"],
        "ng": seed["ng"],
        "nq": seed["nq"],
        "evolution_performed": False,
        "radius_newton_rerun": False,
        "source_columns_unchanged": bool(unchanged),
        "ansatz": _ansatz_report(state, constants),
        "stored_quadrature_radius_gap": _max_abs(prolonged - seed["quadrature_radius"]),
        "stored_npz_projected_residual": seed["stored_projected_residual"],
        "stored_npz_full_operator": seed["stored_full_operator"],
        "frozen_rho_gap_max": _max_abs(source_rho - seed["rho"]),
        "rho_independent_of_r_max": variation,
        "rho_independent_assumption": bool(variation <= 1e-8),
        "frozen_rho_assessment": frozen,
        "source_rho_assessment": fresh,
        "constructor_current": current,
        "quadrature": quadrature,
        "polynomial": polynomial,
        "total_certified": False,
        "geometry_error_upper_bound": None,
    }


def _constant_current_control(grid):
    current = np.full(grid.nq, 0.3)
    snapshot = current.copy()
    momentum, info = shift_momentum(grid, current)
    residual = shift_residual(grid, momentum, current)
    return {
        "current_unchanged": bool(_max_abs(current - snapshot) == 0.0),
        "source_mean_subtracted": info["source_mean_subtracted"],
        "current_mean": float(info["current_mean"]),
        "residual_mean": float(np.mean(residual)),
        "residual_minus_mean_max": _max_abs(residual - np.mean(residual)),
        "p_Q_max": float(info["p_Q_max"]),
        "mean_remains": True,
    }


def _high_mode_omission(grid):
    constants = gauge_constants(grid.fine)
    mode = int(grid.ng // 2 + 5)
    if mode >= grid.nq // 2:
        raise ValueError("quadrature grid cannot hold the omitted mode without aliasing")
    coordinate = grid.xi_q
    amplitude = 0.02
    y_star = np.full(grid.nq, 2.0)
    G = np.full(grid.nq, constants["Q"] ** 2 * y_star ** 4)
    omitted = amplitude * np.cos(2.0 * math.pi * mode * coordinate / grid.length)
    y = y_star + omitted
    Ry = action_ry(grid.fine, y, G, constants)
    fourier = dual_norm_E_parts(grid.fine, Ry, grid.ng // 2)
    matrix = represented_dual_matrix(grid, Ry)
    split = strong_split(grid, Ry)
    return {
        "mode": mode,
        "amplitude": amplitude,
        "represented_weak_dual_E": matrix["dual_E"],
        "represented_fourier_dual_E": fourier["represented"],
        "unresolved_weak_dual_E": fourier["unresolved"],
        "unresolved_strong_max": split["unresolved_max"],
        "full_strong_max": split["full_max"],
        "omitted_energy_norm_E": energy_norm(grid.fine, omitted),
        "represented_pull_max": split["represented_pull_max"],
        "total_certified": False,
        "represented_weak_certifies_total": False,
        "geometry_error_upper_bound": None,
        "tail_bound": None,
        "represented_piece_is_nonlinear_image": True,
    }


def manufactured_controls():
    """Small-grid checks. These do not touch the v5 seed."""
    grid = build_grid(32, quadrature=128)
    system = grid.fine
    constants = gauge_constants(system)
    coordinate = system.xi
    length = grid.length
    y = 1.7 + 0.08 * np.cos(2.0 * math.pi * coordinate / length)
    y = y + 0.02 * np.cos(4.0 * math.pi * coordinate / length)
    rho = 0.15 + 0.05 * np.cos(2.0 * math.pi * coordinate / length)
    smooth = compare_owner_and_action(system, y, rho)
    rough_mode = 40
    rough = 1.4 + 0.2 * np.cos(2.0 * math.pi * rough_mode * coordinate / length)
    rough_gap = compare_owner_and_action(system, rough, rho)
    G = bracket_G(rho, constants)
    test = 0.03 * np.cos(2.0 * math.pi * coordinate / length)
    derivative_gap = weak_pairing(system, y, G, test, constants) - _nodal(
        system, action_ry(system, y, G, constants) * test
    )
    left = action_functional(system, y + 1e-5 * test, G, constants)
    center = action_functional(system, y, G, constants)
    right = action_functional(system, y - 1e-5 * test, G, constants)
    hessian_gap = (left - 2.0 * center + right) / 1e-10 - hessian_quadratic(system, y, G, test, constants)
    positive = coercivity_witness(system, y, G, rho_independent=True)
    hessian_value = hessian_quadratic(system, y, G, test, constants)
    y_star = np.full(system.points, 2.0)
    G_star = np.full(system.points, constants["Q"] ** 2 * 16.0)
    y_near = y_star + 1e-3 * np.cos(2.0 * math.pi * coordinate / length)
    y_hi = float(np.max(y_near))
    c_seg = constants["Q"] ** 2 + 3.0 * float(G_star[0]) / y_hi ** 4
    inequality = observe_discrete_inequality(system, y_near, y_star, G_star, c_seg)
    negative_G = np.full(system.points, -1.0)
    negative_y = np.full(system.points, 1.0)
    negative = coercivity_witness(system, negative_y, negative_G, rho_independent=True)
    negative_hessian = hessian_quadratic(system, negative_y, negative_G, np.ones(system.points), constants)
    omission = _high_mode_omission(grid)
    mean_control = _constant_current_control(grid)
    return {
        "grid": {"nf": grid.nf, "ng": grid.ng, "nq": grid.nq},
        "smooth_equivalence": smooth,
        "rough_mode": rough_mode,
        "rough_equivalence": rough_gap,
        "weak_versus_strong_pairing_gap": derivative_gap,
        "hessian_second_difference_gap": hessian_gap,
        "positive_G": positive,
        "positive_hessian_on_mode": hessian_value,
        "discrete_segment_inequality": inequality,
        "negative_G": negative,
        "negative_hessian_on_constants": negative_hessian,
        "omitted_unresolved_mode": omission,
        "constructor_constant_current": mean_control,
        "total_certified": False,
    }


def _saved_reference():
    if not V5_JSON.is_file():
        return None
    payload = json.loads(V5_JSON.read_text())
    nf512 = payload.get("nf512") or {}
    nq512 = nf512.get("nq2048") or {}
    return {
        "schema": payload.get("schema"),
        "nf256_nq1024_full_C": payload["nq1024"]["full_C_hamilton_max"],
        "nf256_nq1024_frozen_full_operator": payload["nq1024"]["frozen_rho_full_operator_max"],
        "nf256_nq1024_projected_C": payload["nq1024"]["projected_C_hamilton_max"],
        "nf256_nq2048_full_C": payload["nq2048"]["full_C_hamilton_max"],
        "nf256_current_max": payload["nq1024"]["diagnostics"]["current_max"],
        "nf512_nq2048_full_C": nq512.get("full_C_hamilton_max"),
        "nf512_nq2048_projected_C": nq512.get("projected_C_hamilton_max"),
        "nf512_current_max": (nq512.get("diagnostics") or {}).get("current_max"),
    }


def _polynomial_control(nf256, nf512):
    """Four saved-grid comparisons. Not a rounding enclosure."""
    cases = []
    for seed in (nf256, nf512):
        block = seed.get("polynomial") or {}
        baseline = block.get("baseline_stored_rho")
        doubled = block.get("doubled_resampled_rho")
        if baseline is not None:
            cases.append({"label": f"nf{seed['nf']}_nq{baseline['nq']}", "rho": "stored", **baseline})
        if doubled is not None:
            cases.append({
                "label": f"nf{seed['nf']}_nq{doubled['nq']}",
                "rho": "fft_resample_of_baseline_source_rho",
                **doubled,
            })
    dense_base = ((nf512.get("polynomial") or {}).get("baseline_stored_rho") or {}).get("dense_max")
    dense_double = ((nf512.get("polynomial") or {}).get("doubled_recomputed_source") or {}).get("dense_max")
    product_base = ((nf512.get("polynomial") or {}).get("baseline_stored_rho") or {}).get("fft_product_max")
    product_double = (
        (nf512.get("polynomial") or {}).get("doubled_resampled_rho") or {}
    ).get("fft_product_max")
    return {
        "temporary_script_required": False,
        "comparisons": cases,
        "comparison_count": len(cases),
        "source_rho_modified": False,
        "current_modified": False,
        "rounding_enclosure": None,
        "conditioning_indicator": {
            "nf512_dense_nq2048": dense_base,
            "nf512_dense_nq4096": dense_double,
            "nf512_product_nq2048": product_base,
            "nf512_product_nq4096": product_double,
            "statement": (
                "The dense owner maximum moves from about 9.98e-6 at nq=2048 "
                "to about 4.54e-5 at nq=4096, while the Fourier product residual "
                "stays near 5.07e-6 to 5.40e-6. This is a conditioning indicator, "
                "not a rounding enclosure."
            ),
        },
    }


def _replay_gap(measured, saved):
    if measured is None or saved is None:
        return None
    return float(measured) - float(saved)


class InputBindingError(ValueError):
    """Declared scientific input does not match the bytes available now."""

    def __init__(self, message, limits):
        super().__init__(message)
        self.limits = list(limits)


def _source_hashes():
    paths = GENERATOR_SOURCE_PATHS + SCIENTIFIC_SOURCE_PATHS
    return {relative: _sha256(_LAB_ROOT / relative) for relative in paths}


def assemble_record():
    """Measure the saved v5 seeds. Does not write."""
    started = time.process_time()
    wall = time.perf_counter()
    npz_hash_before = _sha256(V5_NPZ)
    json_hash_before = _sha256(V5_JSON)
    with np.load(V5_NPZ) as payload:
        seeds = {
            "nf256": _load_seed(payload, ""),
            "nf512": _load_seed(payload, "nf512_"),
        }
    manufactured = manufactured_controls()
    # Leave room for the nf512 doubled grid inside the same 60 s CPU budget.
    nf256 = assess_saved_seed(seeds["nf256"], doubled_budget_s=CPU_BUDGET_S * 0.45, started=started)
    nf512 = assess_saved_seed(seeds["nf512"], doubled_budget_s=CPU_BUDGET_S * 0.85, started=started)
    reference = _saved_reference()
    cpu = time.process_time() - started
    record = {
        "schema": SCHEMA,
        "status": STATUS,
        "role": "weak initial residual of the saved spherical seed; not an evolution step",
        "evolution_performed": False,
        "production_evolution_owner": "nsc_spherical_galerkin_coupling",
        "source_changed": False,
        "coefficients_changed": False,
        "action_changed": False,
        "inheritance_recomputed": False,
        "tolerance_used_as_proof": False,
        "total_certified": False,
        "represented_weak_certifies_total": False,
        "geometry_error_upper_bound": None,
        "statement": STATEMENT,
        "checkpoint_expected": "9a9090a",
        "checkpoint_head": _head(),
        "cpu_budget_seconds": CPU_BUDGET_S,
        "cpu_seconds": cpu,
        "wall_seconds": time.perf_counter() - wall,
        "cpu_budget_exceeded": bool(cpu > CPU_BUDGET_S),
        "v5_npz_sha256_before": npz_hash_before,
        "v5_json_sha256_before": json_hash_before,
        "saved_v5_reference": reference,
        "manufactured": manufactured,
        "nf256": nf256,
        "nf512": nf512,
        "polynomial_control": _polynomial_control(nf256, nf512),
    }
    if reference is not None:
        record["replay_gaps"] = {
            "nf256_frozen_minus_saved_frozen_operator": _replay_gap(
                nf256["frozen_rho_assessment"]["full_strong"]["owner_max"],
                reference["nf256_nq1024_frozen_full_operator"],
            ),
            "nf256_source_minus_saved_full_C": _replay_gap(
                nf256["source_rho_assessment"]["full_strong"]["owner_max"],
                reference["nf256_nq1024_full_C"],
            ),
            "nf256_doubled_minus_saved_nq2048": _replay_gap(
                None if not nf256["quadrature"]["computed"] else nf256["quadrature"]["doubled_nq"]["full_strong_owner_max"],
                reference["nf256_nq2048_full_C"],
            ),
            "nf512_source_minus_saved_full_C": _replay_gap(
                nf512["source_rho_assessment"]["full_strong"]["owner_max"],
                reference["nf512_nq2048_full_C"],
            ),
            "nf256_constructor_current_minus_saved_current_max": _replay_gap(
                nf256["constructor_current"]["current_max"],
                reference["nf256_current_max"],
            ),
        }
    payload = _jsonable(record)
    payload["source_hashes"] = _source_hashes()
    payload["v5_npz_sha256_after"] = _sha256(V5_NPZ)
    payload["v5_json_sha256_after"] = _sha256(V5_JSON)
    payload["v5_bytes_unchanged"] = bool(
        payload["v5_npz_sha256_before"] == payload["v5_npz_sha256_after"]
        and payload["v5_json_sha256_before"] == payload["v5_json_sha256_after"]
    )
    return payload


def build_record():
    """Assemble a fresh record. Does not write."""
    return assemble_record()


def write_record(payload, path):
    """Write a fresh assembly to an explicit new path.

    The sealed development record is never the destination. Verification
    does not call this function.
    """
    destination = Path(path)
    if destination.resolve() == RECORD_PATH.resolve():
        raise FileExistsError(
            f"{destination} is the sealed weak record; write an explicit new path"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(_jsonable(payload), indent=2) + "\n")
    return destination


def _declared_hash(value):
    if isinstance(value, str) and value:
        return value
    return None


def _identity_entry(path, declared, *, role, current):
    """Describe one declared hash without inventing a missing digest."""
    declared_hash = _declared_hash(declared)
    available = current is not None
    matches = bool(available and declared_hash is not None and current == declared_hash)
    return {
        "path": path,
        "role": role,
        "available": available,
        "historical_bytes_available": matches,
        "declared_hash": declared_hash,
        "current_hash": current,
        "matches_declared": matches,
        "fabricated_hash": False,
    }


def declared_input_context(saved):
    """Historical inputs declared by a sealed record.

    ``checkpoint_head`` is the execution HEAD stored with that record. It is
    not a scientific input and it is not compared with the live HEAD.
    """
    hashes = saved.get("source_hashes")
    if not isinstance(hashes, dict):
        hashes = {}
    return {
        "schema": saved.get("schema"),
        "v5_json_sha256": saved.get("v5_json_sha256_before"),
        "v5_npz_sha256": saved.get("v5_npz_sha256_before"),
        "v5_json_sha256_after": saved.get("v5_json_sha256_after"),
        "v5_npz_sha256_after": saved.get("v5_npz_sha256_after"),
        "v5_bytes_unchanged": saved.get("v5_bytes_unchanged"),
        "source_hashes": hashes,
        "historical_checkpoint_head": saved.get("checkpoint_head"),
        "checkpoint_expected": saved.get("checkpoint_expected"),
    }


def binding_report(saved, *, source_root=None):
    """Check declared v5 bytes and source hashes. Does not write or remeasure."""
    context = declared_input_context(saved)
    root = _LAB_ROOT if source_root is None else Path(source_root)
    entries = []
    if context["schema"] != SCHEMA:
        entries.append({
            "path": "schema",
            "role": "scientific_setting",
            "available": True,
            "historical_bytes_available": False,
            "declared_hash": None,
            "current_hash": None,
            "matches_declared": False,
            "fabricated_hash": False,
            "reason": f"schema is {context['schema']!r}",
        })
    blob_specs = (
        ("v5_json", V5_JSON, context["v5_json_sha256"]),
        ("v5_npz", V5_NPZ, context["v5_npz_sha256"]),
    )
    for name, file_path, declared in blob_specs:
        current = _sha256(file_path)
        entries.append(_identity_entry(name, declared, role="scientific_input", current=current))
    if context["v5_bytes_unchanged"] is not True:
        entries.append({
            "path": "v5_bytes_unchanged",
            "role": "scientific_setting",
            "available": True,
            "historical_bytes_available": False,
            "declared_hash": None,
            "current_hash": None,
            "matches_declared": False,
            "fabricated_hash": False,
            "reason": "declared v5 bytes were not unchanged",
        })
    after_pairs = (
        ("v5_json_sha256_after", context["v5_json_sha256_after"], context["v5_json_sha256"]),
        ("v5_npz_sha256_after", context["v5_npz_sha256_after"], context["v5_npz_sha256"]),
    )
    for name, after, before in after_pairs:
        if _declared_hash(after) != _declared_hash(before):
            entries.append({
                "path": name,
                "role": "scientific_setting",
                "available": True,
                "historical_bytes_available": False,
                "declared_hash": _declared_hash(before),
                "current_hash": _declared_hash(after),
                "matches_declared": False,
                "fabricated_hash": False,
                "reason": "declared before and after hashes differ",
            })
    hashes = context["source_hashes"]
    known = set(SCIENTIFIC_SOURCE_PATHS + GENERATOR_SOURCE_PATHS)
    ordered = list(SCIENTIFIC_SOURCE_PATHS) + list(GENERATOR_SOURCE_PATHS)
    ordered.extend(sorted(set(hashes) - known))
    for relative in ordered:
        role = "generator" if relative in GENERATOR_SOURCE_PATHS else "scientific_source"
        file_path = root / relative
        current = _sha256(file_path)
        entries.append(_identity_entry(
            relative, hashes.get(relative), role=role, current=current,
        ))
    blocking = [
        item for item in entries
        if item["role"] != "generator" and item["matches_declared"] is not True
    ]
    limits = [item for item in entries if item["matches_declared"] is not True]
    message = "declared input mismatch: " + ", ".join(item["path"] for item in blocking[:12])
    return {
        "ok": not blocking,
        "message": message if blocking else "declared inputs match",
        "head_equality_required": False,
        "historical_checkpoint_head": context["historical_checkpoint_head"],
        "entries": entries,
        "blocking": blocking,
        "limits": limits,
        "fabricated_hash": False,
        "scientific_sources_match": not any(
            item["role"] == "scientific_source" and item["matches_declared"] is not True
            for item in entries
        ),
    }


def _numbers_differ(saved, fresh):
    left = float(saved)
    right = float(fresh)
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) > REPLAY_ABSOLUTE_TOLERANCE + REPLAY_RELATIVE_TOLERANCE * scale


def _compare_scientific(saved, fresh, path, mismatches, root=False):
    if isinstance(fresh, dict):
        if not isinstance(saved, dict):
            mismatches.append(path or "root")
            return
        saved_keys = set(saved)
        fresh_keys = set(fresh)
        if root:
            saved_keys -= VOLATILE_EXECUTION_FIELDS
            fresh_keys -= VOLATILE_EXECUTION_FIELDS
            saved_keys.discard("source_hashes")
            fresh_keys.discard("source_hashes")
        for key in sorted(saved_keys - fresh_keys):
            mismatches.append(f"{path}.{key}" if path else key)
        for key in sorted(fresh_keys - saved_keys):
            mismatches.append(f"{path}.{key}" if path else key)
        for key in sorted(saved_keys & fresh_keys):
            child = f"{path}.{key}" if path else key
            _compare_scientific(saved[key], fresh[key], child, mismatches)
        return
    if isinstance(fresh, list):
        if not isinstance(saved, list) or len(saved) != len(fresh):
            mismatches.append(path or "root")
            return
        for index, (left, right) in enumerate(zip(saved, fresh)):
            _compare_scientific(left, right, f"{path}[{index}]", mismatches)
        return
    if isinstance(fresh, bool) or isinstance(saved, bool):
        if saved is not fresh:
            mismatches.append(path or "root")
        return
    if isinstance(fresh, (int, float)) and isinstance(saved, (int, float)):
        if _numbers_differ(saved, fresh):
            mismatches.append(path or "root")
        return
    if fresh is None or saved is None or isinstance(fresh, str):
        if saved != fresh:
            mismatches.append(path or "root")
        return
    mismatches.append(path or "root")


def scientific_mismatches(saved, fresh):
    """Paths whose scientific values differ.

    Live checkpoint, CPU, and wall time are omitted. Source identity is
    ``binding_report``, not this comparison.
    """
    mismatches = []
    _compare_scientific(saved, fresh, "", mismatches, root=True)
    return mismatches


def verify_saved(path=RECORD_PATH):
    """Remeasure scientific fields against a sealed record. Does not write."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    saved_bytes = path.read_bytes()
    v5_json_bytes = V5_JSON.read_bytes()
    v5_npz_bytes = V5_NPZ.read_bytes()
    saved = json.loads(saved_bytes)
    report = binding_report(saved)
    if not report["ok"]:
        raise InputBindingError(report["message"], report["blocking"])
    fresh = assemble_record()
    mismatches = scientific_mismatches(saved, fresh)
    record_unchanged = path.read_bytes() == saved_bytes
    v5_unchanged = V5_JSON.read_bytes() == v5_json_bytes and V5_NPZ.read_bytes() == v5_npz_bytes
    if not record_unchanged or not v5_unchanged:
        raise RuntimeError("verification rewrote a sealed input")
    if mismatches:
        raise AssertionError("scientific field mismatch: " + ", ".join(mismatches[:12]))
    if float(fresh["cpu_seconds"]) > CPU_BUDGET_S:
        raise AssertionError("CPU limit exceeded during remeasurement")
    generator_limits = [item for item in report["limits"] if item["role"] == "generator"]
    return {
        "status": saved.get("status"),
        "ok": True,
        "wrote": False,
        "head_equality_required": False,
        "historical_checkpoint_head": saved.get("checkpoint_head"),
        "current_checkpoint_head": fresh.get("checkpoint_head"),
        "cpu_seconds": fresh["cpu_seconds"],
        "v5_bytes_unchanged": True,
        "record_bytes_unchanged": True,
        "limits": generator_limits,
        "fabricated_hash": False,
        "measured": fresh,
        "sealed": saved,
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not math.isfinite(number):
            return {"nonfinite": True}
        return number
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None or isinstance(value, str):
        return value
    return str(value)
