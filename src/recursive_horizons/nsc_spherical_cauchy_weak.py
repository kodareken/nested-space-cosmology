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
tail. ``assemble_record`` builds no tail majorant and emits no total bound.

``initial_lapse_component_bound`` encloses one conditional initial-lapse
component for the saved ``n_f=512`` state. It does not step the state and
does not write the sealed weak record or the v5 bytes.

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
from mpmath import iv as _iv

from .provenance import resolve_pinned_source_bytes
from .nsc_conformal_adm_source import sector_multiplicity
from .nsc_spherical_cauchy_data import shift_momentum, shift_residual
from .nsc_spherical_coupling import (
    CALIBRATION,
    KAPPA,
    OCCUPATIONS,
    PERIOD,
    TOL_INITIAL_CONSTRAINT,
    profile_coordinate,
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
ERROR_RECORD_PATH = _LAB_ROOT / "results" / "development" / "nsc-spherical-cauchy-error-v1.json"
ERROR_SCHEMA = "NSC-SPHERICAL-CAUCHY-INITIAL-LAPSE-v1"
LAPSE_COMPONENT_STATUS = "FINITE_CONDITIONAL_INITIAL_LAPSE_BOUND"
LAPSE_BOUND_CPU_S = 60.0
LAPSE_BOUND_DPS = 25
LAPSE_FINE_COUNT = 65536
EPISODE_JSON = _LAB_ROOT / "results" / "development" / "nsc-spherical-feedback-episode-v1.json"
V5_JSON = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-refinement-v5.json"
V5_NPZ = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-refinement-v5.npz"
LAPSE_BOUND_SOURCE_PATHS = (
    "src/recursive_horizons/nsc_spherical_coupling.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_conformal_adm_source.py",
    "src/recursive_horizons/nsc_covariant_operator.py",
    "src/recursive_horizons/nsc_spherical_feedback_action.py",
    "src/recursive_horizons/nsc_spherical_cauchy_weak.py",
)
# The helper is the historical generator of the sealed lapse record. A later
# edit is a readonly limit, not a reason to rehash that record.
LAPSE_GENERATOR_PATHS = frozenset({
    "src/recursive_horizons/nsc_spherical_cauchy_weak.py",
})
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


def binding_report(saved, *, source_root=None, source_ref=None):
    """Check v5 bytes and formula hashes in a named current or historical domain.

    Current requests remain strict. A historical ``source_ref`` authenticates
    exact formula bytes at that full commit without rebinding the record.
    Generator differences remain limits. Does not write or remeasure.
    """
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
        entry = _identity_entry(relative, hashes.get(relative), role=role, current=current)
        if source_ref is not None and role == "scientific_source":
            try:
                raw = resolve_pinned_source_bytes(
                    root.parent, root.name + "/" + relative, hashes.get(relative),
                    commit=source_ref)
            except (ValueError, RuntimeError):
                entry["matches_declared"] = False
                entry["historical_bytes_available"] = False
            else:
                entry["available"] = True
                entry["historical_bytes_available"] = True
                entry["historical_hash"] = hashlib.sha256(raw).hexdigest()
                entry["matches_declared"] = True
                entry["source_origin"] = source_ref
        elif role == "scientific_source":
            entry["source_origin"] = "working-tree"
        entries.append(entry)
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
        "source_ref": source_ref,
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


def verify_saved(path=RECORD_PATH, *, source_ref=None):
    """Remeasure scientific fields against a sealed record. Does not write."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    saved_bytes = path.read_bytes()
    v5_json_bytes = V5_JSON.read_bytes()
    v5_npz_bytes = V5_NPZ.read_bytes()
    saved = json.loads(saved_bytes)
    report = binding_report(saved, source_ref=source_ref)
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
        "source_ref": source_ref,
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


def polynomial_degree_caps(fermions, geometry):
    """Degree caps traced from the declared bands, not from a sampled tail.

    Antiperiodic fermion modes are the half-integers from ``-(nf/2-1/2)`` to
    ``nf/2-1/2``. A product of two such modes, including a conjugate, is a
    periodic mode of integer degree at most ``nf-1``. The odd geometry band
    has degree ``ng//2``, so ``r r''`` and ``(r')^2`` have degree at most
    ``ng-1`` when ``ng`` is odd. With ``ng = nf-1`` the source degree is the
    larger one, and ``P`` has degree at most ``nf-1``.
    """
    fermions = int(fermions)
    geometry = int(geometry)
    if fermions < 10 or fermions % 2 or geometry != fermions - 1:
        raise ValueError("degree caps expect even nf and odd geometry ng = nf-1")
    degree_radius = geometry // 2
    degree_source = fermions - 1
    degree_radius_products = 2 * degree_radius
    return {
        "degree_r": degree_radius,
        "degree_source": degree_source,
        "degree_radius_products": degree_radius_products,
        "degree_P": max(degree_radius_products, degree_source),
        "highest_antiperiodic_mode": fermions / 2 - 0.5,
    }


def _numpy_ap_coefficients(samples):
    samples = np.asarray(samples)
    count = samples.shape[0]
    index = np.arange(count)
    demodulated = np.exp(-1j * math.pi * index / count).reshape(
        (-1,) + (1,) * (samples.ndim - 1)
    )
    coefficients = np.fft.fft(samples * demodulated, axis=0) / count
    integers = np.fft.fftfreq(count) * count
    return integers, coefficients


def _numpy_synthesize_ap(integers, coefficients, nodes, length):
    phase = np.exp(
        2j * np.pi * (np.asarray(integers) + 0.5)[:, None] * np.asarray(nodes)[None, :] / length
    )
    return phase.T @ np.asarray(coefficients)


def bandlimited_matter_density(phi0, phi1, length, count, gauge):
    """Continuum ``rho = force_L / dx`` of the antiperiodic interpolant.

    ``phi`` are half-density samples on the fermion nodes. ``P = -i d/dx``
    multiplies mode ``m`` by ``2π m / length``. This is the symbol of
    ``CovariantStaticMetric.momentum_matrix`` at ``eta = 1/2``. ``force_L``
    does not read ``r``, ``chi``, the momenta, the lapse profile, or the shift.
    """
    phi0 = np.asarray(phi0)
    phi1 = np.asarray(phi1)
    if phi0.ndim != 2 or phi1.shape != phi0.shape:
        raise ValueError("phi0 and phi1 must share the shape (nf, columns)")
    if phi0.shape[1] != len(OCCUPATIONS):
        raise ValueError("column count does not match the declared occupations")
    fermions = phi0.shape[0]
    count = int(count)
    if count < fermions or count % 2:
        raise ValueError("evaluation grid must be even and at least nf")
    length = float(length)
    sqrt_dx = math.sqrt(length / fermions)
    integers, coeff0 = _numpy_ap_coefficients(phi0)
    _, coeff1 = _numpy_ap_coefficients(phi1)
    symbol = 2.0 * math.pi * (integers + 0.5) / length
    nodes = np.arange(count) * (length / count)
    psi0 = _numpy_synthesize_ap(integers, coeff0 / sqrt_dx, nodes, length)
    psi1 = _numpy_synthesize_ap(integers, coeff1 / sqrt_dx, nodes, length)
    momentum0 = _numpy_synthesize_ap(integers, coeff0 / sqrt_dx * symbol[:, None], nodes, length)
    momentum1 = _numpy_synthesize_ap(integers, coeff1 / sqrt_dx * symbol[:, None], nodes, length)
    occupations = np.asarray(OCCUPATIONS, dtype=float)
    kinetic = np.conjugate(psi0) * (-1j * momentum1) + np.conjugate(psi1) * (1j * momentum0)
    axial = np.conjugate(psi0) * psi1 + np.conjugate(psi1) * psi0
    k_density = np.sum(kinetic * occupations[None, :], axis=1).real
    s1_density = np.sum(axial * occupations[None, :], axis=1).real
    multiplicity = float(sector_multiplicity(KAPPA))
    return multiplicity * (k_density / float(gauge) + float(KAPPA) * s1_density)


def _float_hi(value):
    """Outward float64 image of an interval's right endpoint."""
    number = float(value.b)
    return math.nextafter(number, math.inf)


def _float_lo(value):
    """Outward float64 image of an interval's left endpoint."""
    number = float(value.a)
    return math.nextafter(number, -math.inf)


def _as_mpc(value):
    if isinstance(value, complex) or np.iscomplexobj(value):
        return _iv.mpc(_iv.mpf(float(np.real(value))), _iv.mpf(float(np.imag(value))))
    return _iv.mpc(_iv.mpf(float(value)))


def _conjugate(value):
    """Interval conjugate. ``iv.mpc.conjugate`` is not usable in mpmath 1.3."""
    return _iv.mpc(value.real, -value.imag)


def _bitreverse(values):
    count = len(values)
    reversed_values = list(values)
    index = 0
    for position in range(1, count):
        bit = count >> 1
        while index & bit:
            index ^= bit
            bit >>= 1
        index ^= bit
        if position < index:
            reversed_values[position], reversed_values[index] = (
                reversed_values[index],
                reversed_values[position],
            )
    return reversed_values


def _interval_fft(values, inverse=False):
    """Radix-2 DFT with mpmath outward rounding. Length must be a power of two."""
    count = len(values)
    if count < 2 or count & (count - 1):
        raise ValueError("interval FFT length must be a power of two")
    spectrum = _bitreverse(values)
    length = 2
    sign = 1 if inverse else -1
    while length <= count:
        half = length // 2
        step = _iv.exp(sign * 2j * _iv.pi / length)
        for start in range(0, count, length):
            twiddle = _iv.mpc(1)
            for offset in range(half):
                left = spectrum[start + offset]
                right = spectrum[start + offset + half] * twiddle
                spectrum[start + offset] = left + right
                spectrum[start + offset + half] = left - right
                twiddle *= step
        length *= 2
    if inverse:
        scale = _iv.mpf(1) / count
        spectrum = [item * scale for item in spectrum]
    return spectrum


def _chirp(index, count, sign):
    return _iv.exp(sign * 1j * _iv.pi * index * index / count)


def _interval_dft(values):
    """Unnormalized forward DFT of any length, by a power-of-two Bluestein convolution."""
    count = len(values)
    if count < 1:
        raise ValueError("empty DFT")
    if count & (count - 1) == 0:
        return _interval_fft(values)
    pad = 1
    while pad < 2 * count - 1:
        pad *= 2
    series = [values[index] * _chirp(index, count, -1) for index in range(count)]
    series.extend(_iv.mpc(0) for _ in range(pad - count))
    kernel = [_iv.mpc(0) for _ in range(pad)]
    kernel[0] = _chirp(0, count, 1)
    for index in range(1, count):
        value = _chirp(index, count, 1)
        kernel[index] = value
        kernel[pad - index] = value
    transformed = _interval_fft(series)
    kernel_transform = _interval_fft(kernel)
    convolved = _interval_fft(
        [left * right for left, right in zip(transformed, kernel_transform)],
        inverse=True,
    )
    return [_chirp(index, count, -1) * convolved[index] for index in range(count)]


def _closest_to_zero(interval):
    if float(interval.a) <= 0.0 <= float(interval.b):
        return _iv.mpf(0)
    if float(interval.b) < 0.0:
        return interval.b
    return interval.a


def _abs2_bounds(value):
    """Lower and upper bounds of ``|z|^2`` on a complex interval rectangle."""
    upper = _iv.mpf(0)
    for real_end in (value.real.a, value.real.b):
        for imag_end in (value.imag.a, value.imag.b):
            corner = _iv.mpf(real_end) ** 2 + _iv.mpf(imag_end) ** 2
            if float(corner.b) > float(upper.b):
                upper = corner
    lower = _iv.mpf(_closest_to_zero(value.real)) ** 2 + _iv.mpf(_closest_to_zero(value.imag)) ** 2
    return lower, upper


def _place_modes(coefficients, modes, count, scale_mode):
    """Inverse FFT of normalized coefficients. ``scale_mode(mode)`` is exact."""
    spectrum = [_iv.mpc(0) for _ in range(count)]
    for coefficient, mode in zip(coefficients, modes):
        mode = int(mode)
        if abs(mode) > count // 2:
            raise ValueError("mode does not fit the evaluation grid")
        slot = mode if mode >= 0 else count + mode
        spectrum[slot] = coefficient * scale_mode(mode) * count
    return _interval_fft(spectrum, inverse=True)


def _sobolev_factor(length):
    """``||v||_∞ ≤ factor ||v||_E`` on the circle of this length.

    With ``v = Σ c_k e^{2π i k x/L}`` and ``||v||_E^2 = L Σ |c_k|^2 (1+(2π k/L)^2)``,
    Cauchy-Schwarz gives the factor ``sqrt(coth(L/2) / 2)``.
    """
    half = _iv.mpf(length) / 2
    exponential = _iv.exp(2 * half)
    hyperbolic = (exponential + 1) / (exponential - 1)
    return _iv.sqrt(hyperbolic / 2)


def _array_sha256(array):
    values = np.ascontiguousarray(np.asarray(array))
    digest = hashlib.sha256()
    digest.update(str(values.shape).encode())
    digest.update(str(values.dtype).encode())
    digest.update(values.tobytes())
    return digest.hexdigest()


def _lapse_source_hashes():
    return {relative: _sha256(_LAB_ROOT / relative) for relative in LAPSE_BOUND_SOURCE_PATHS}


def _areal_radius_effects():
    if not EPISODE_JSON.is_file():
        return None
    payload = json.loads(EPISODE_JSON.read_text())
    indicators = payload.get("refinement_indicators")
    if not isinstance(indicators, list):
        return None
    found = {}
    for item in indicators:
        name = item.get("name") if isinstance(item, dict) else None
        if name in ("areal_radius_min", "areal_radius_max"):
            found[name] = item
    if set(found) != {"areal_radius_min", "areal_radius_max"}:
        return None
    return {
        "path": "results/development/nsc-spherical-feedback-episode-v1.json",
        "sha256": _sha256(EPISODE_JSON),
        "areal_radius_min_effect": float(found["areal_radius_min"]["space_effect_scale"]),
        "areal_radius_max_effect": float(found["areal_radius_max"]["space_effect_scale"]),
        "bar": "one percent of the saved areal-radius effect scale",
        "bar_fraction": 0.01,
    }


def initial_lapse_binding_report(saved, *, source_ref=None):
    """Reject reuse when the saved nf512 source, settings, or formula files differ.

    The episode hash guards only the effect comparison. A mismatch there does
    not by itself erase a lapse bound whose v5 and formula hashes still match.
    The weak-helper hash is the historical generator. A difference is reported
    as a limit and does not rehash the sealed scientific record.
    Current requests are strict by default. Historical replay requires an
    explicit full ``source_ref`` and authenticates every scientific formula
    source at that commit; numerical inputs still bind the resident bytes.
    """
    if not isinstance(saved, dict) or saved.get("schema") != ERROR_SCHEMA:
        return {
            "ok": False,
            "effect_comparison_ok": False,
            "reason": "schema",
            "generator_limits": [],
            "fabricated_hash": False,
        }
    reasons = []
    if _sha256(V5_NPZ) != saved.get("inputs", {}).get("v5_npz_sha256"):
        reasons.append("v5 npz")
    if _sha256(V5_JSON) != saved.get("inputs", {}).get("v5_json_sha256"):
        reasons.append("v5 json")
    if _sha256(RECORD_PATH) != saved.get("inputs", {}).get("sealed_weak_sha256"):
        reasons.append("sealed weak record")
    current_sources = _lapse_source_hashes()
    declared_sources = saved.get("inputs", {}).get("source_hashes") or {}
    limits = []
    origins = {}
    for path in sorted(set(current_sources) | set(declared_sources)):
        declared_hash = declared_sources.get(path)
        current_hash = current_sources.get(path)
        entry = {
            "path": path,
            "declared_hash": declared_hash,
            "current_hash": current_hash,
            "role": "generator" if path in LAPSE_GENERATOR_PATHS else "scientific_source",
            "fabricated_hash": False,
        }
        if path in LAPSE_GENERATOR_PATHS:
            if declared_hash != current_hash:
                limits.append(entry)
            continue
        if source_ref is not None:
            try:
                resolve_pinned_source_bytes(
                    _LAB_ROOT.parent, "lab/" + path, declared_hash,
                    commit=source_ref)
            except (ValueError, RuntimeError):
                reasons.append(path)
            else:
                origins[path] = source_ref
        elif declared_hash != current_hash:
            reasons.append(path)
        else:
            origins[path] = "working-tree"
    settings = saved.get("settings") or {}
    if settings.get("occupations") != [float(value) for value in OCCUPATIONS]:
        reasons.append("occupations")
    if settings.get("kappa") != int(KAPPA):
        reasons.append("kappa")
    if settings.get("multiplicity") != int(sector_multiplicity(KAPPA)):
        reasons.append("multiplicity")
    effects = saved.get("effect_comparison") or {}
    episode_ok = effects.get("sha256") == _sha256(EPISODE_JSON)
    return {
        "ok": not reasons,
        "effect_comparison_ok": bool(episode_ok and not reasons),
        "reason": None if not reasons else ", ".join(reasons),
        "generator_limits": limits,
        "source_ref": source_ref,
        "source_origins": origins,
        "fabricated_hash": False,
    }


def _invalid_lapse_record(reason, *, started, diagnostics, status="MISSING_PRIMITIVE"):
    cpu = time.process_time() - started
    return {
        "schema": ERROR_SCHEMA,
        "status": status,
        "reason": reason,
        "geometry_error_upper_bound": None,
        "pointwise_radius_error_upper_bound": None,
        "shift_constraint_error_bound": None,
        "evolution_error_bound": None,
        "minimizer_existence_proved": None,
        "total_certified": False,
        "continuum_solution_claimed": False,
        "mean_current_removed": False,
        "cpu_seconds": cpu,
        "cpu_budget_seconds": LAPSE_BOUND_CPU_S,
        "cpu_budget_exceeded": bool(cpu > LAPSE_BOUND_CPU_S),
        "diagnostics": diagnostics,
    }


def initial_lapse_component_bound():
    """Conditional H1 bound for the saved ``n_f=512`` initial lapse geometry.

    The mathematical object is the odd trigonometric interpolant of the stored
    radius and the antiperiodic interpolant of the stored fermion columns.
    ``P = 2 r r'' - (r')^2 - Q^2 r^2 + G`` is enclosed as a trigonometric
    polynomial. ``||Ry||_* ≤ ||P||_L2 / y_min^3`` and ``μ = min(4, Q^2)`` then
    bound the distance to a positive critical point of ``J``, if one exists.
    The shift constraint is not part of ``J`` and stays null.
    """
    started = time.process_time()
    previous_dps = _iv.dps
    _iv.dps = LAPSE_BOUND_DPS
    v5_npz_before = V5_NPZ.read_bytes()
    v5_json_before = V5_JSON.read_bytes()
    sealed_before = RECORD_PATH.read_bytes()
    try:
        return _initial_lapse_component_bound(started, v5_npz_before, v5_json_before, sealed_before)
    finally:
        _iv.dps = previous_dps
        if V5_NPZ.read_bytes() != v5_npz_before or V5_JSON.read_bytes() != v5_json_before:
            raise RuntimeError("initial lapse bound rewrote v5 bytes")
        if RECORD_PATH.read_bytes() != sealed_before:
            raise RuntimeError("initial lapse bound rewrote the sealed weak record")


def _initial_lapse_component_bound(started, v5_npz_before, v5_json_before, sealed_before):
    del v5_npz_before, v5_json_before, sealed_before
    with np.load(V5_NPZ) as payload:
        fermions = int(payload["nf512_nf"])
        geometry = int(payload["nf512_ng"])
        quadrature = int(payload["nf512_nq"])
        length = float(payload["nf512_length"])
        phi0 = np.array(payload["nf512_columns_phi0"], copy=True)
        phi1 = np.array(payload["nf512_columns_phi1"], copy=True)
        radius = np.array(payload["nf512_radius"], dtype=float, copy=True)
        gauge_samples = np.array(payload["nf512_geometry_Q"], dtype=float, copy=True)
        chi_max = _max_abs(payload["nf512_geometry_chi"])
        p_r_max = _max_abs(payload["nf512_geometry_p_r"])
        p_chi_max = _max_abs(payload["nf512_geometry_p_chi"])
        p_q_max = _max_abs(payload["nf512_geometry_p_Q"])
        occupations = np.array(payload["nf512_occupations"], dtype=float, copy=True)
        phi0_hash = _array_sha256(phi0)
        phi1_hash = _array_sha256(phi1)
        radius_hash = _array_sha256(radius)
        gauge_hash = _array_sha256(gauge_samples)
    caps = polynomial_degree_caps(fermions, geometry)
    if quadrature // 2 <= caps["degree_P"] or length != 8.0:
        return _invalid_lapse_record(
            "saved grid is not the alias-free nf512 quadrature",
            started=started,
            diagnostics={"caps": caps, "nq": quadrature, "length": length},
        )
    if not np.array_equal(occupations, np.asarray(OCCUPATIONS, dtype=float)):
        return _invalid_lapse_record(
            "saved occupations are not the declared source occupations",
            started=started,
            diagnostics=None,
        )
    if chi_max != 0.0 or p_r_max != 0.0 or p_chi_max != 0.0 or p_q_max != 0.0:
        return _invalid_lapse_record(
            "saved initial momenta or chi are not exactly zero",
            started=started,
            diagnostics=None,
        )
    if float(np.min(gauge_samples)) != float(np.max(gauge_samples)) or float(gauge_samples[0]) <= 0.0:
        return _invalid_lapse_record(
            "saved Q is not one positive constant",
            started=started,
            diagnostics=None,
        )
    gauge = _iv.mpf(float(gauge_samples[0]))
    constants = gauge_constants(build_grid(32, quadrature=128).fine)
    if float(constants["Q"]) != float(gauge_samples[0]):
        return _invalid_lapse_record(
            "saved Q is not the gauge ratio b0/a0",
            started=started,
            diagnostics=None,
        )
    area = _iv.mpf(constants["A"])
    rmag2 = _iv.mpf(constants["rmag2"])
    eight_pi_a = 8 * _iv.pi * area
    multiplicity = _iv.mpf(sector_multiplicity(KAPPA))
    kappa = _iv.mpf(KAPPA)
    sqrt_dx = _iv.sqrt(_iv.mpf(length) / fermions)
    wavenumber = 2 * _iv.pi / _iv.mpf(length)

    def ap_columns(samples):
        count = samples.shape[0]
        index = np.arange(count)
        demod = np.exp(-1j * math.pi * index / count)
        integers = np.fft.fftfreq(count) * count
        columns = []
        scale = _iv.mpf(1) / count
        for column in range(samples.shape[1]):
            modulated = [_as_mpc(value) for value in samples[:, column] * demod]
            columns.append([item * scale for item in _interval_fft(modulated)])
        return integers, columns

    integers, coeff0 = ap_columns(phi0)
    _, coeff1 = ap_columns(phi1)
    symbols = [wavenumber * (_iv.mpf(int(mode)) + _iv.mpf("0.5")) for mode in integers]
    modulation = [_iv.exp(1j * _iv.pi * _iv.mpf(index) / quadrature) for index in range(quadrature)]

    def synthesize(columns, with_momentum):
        output = []
        for column in columns:
            spectrum = [_iv.mpc(0) for _ in range(quadrature)]
            for slot, mode in enumerate(integers):
                mode = int(mode)
                coefficient = column[slot] / sqrt_dx
                if with_momentum:
                    coefficient *= symbols[slot]
                slot_out = mode if mode >= 0 else quadrature + mode
                spectrum[slot_out] = coefficient * quadrature
            values = _interval_fft(spectrum, inverse=True)
            output.append([value * factor for value, factor in zip(values, modulation)])
        return output

    psi0 = synthesize(coeff0, False)
    psi1 = synthesize(coeff1, False)
    momentum0 = synthesize(coeff0, True)
    momentum1 = synthesize(coeff1, True)
    occupations_iv = [_iv.mpf(float(value)) for value in OCCUPATIONS]
    rho = []
    imaginary_upper = _iv.mpf(0)
    for node in range(quadrature):
        kinetic = _iv.mpc(0)
        axial = _iv.mpc(0)
        for column in range(6):
            weight = occupations_iv[column]
            kinetic += weight * (
                _conjugate(psi0[column][node]) * (-1j * momentum1[column][node])
                + _conjugate(psi1[column][node]) * (1j * momentum0[column][node])
            )
            axial += weight * (
                _conjugate(psi0[column][node]) * psi1[column][node]
                + _conjugate(psi1[column][node]) * psi0[column][node]
            )
        rho.append(multiplicity * (kinetic.real / gauge + kappa * axial.real))
        for part in (kinetic.imag, axial.imag):
            imaginary_upper = max(imaginary_upper, part.b, -part.a)
    capital_g = [gauge ** 2 * rmag2 + gauge * value / eight_pi_a for value in rho]
    radius_coefficients = [
        item / geometry for item in _interval_dft([_as_mpc(value) for value in radius])
    ]
    radius_modes = np.fft.fftfreq(geometry) * geometry

    def derive(power):
        def scale_mode(mode):
            if power == 0:
                return _iv.mpc(1)
            return (1j * wavenumber * _iv.mpf(mode)) ** power

        values = _place_modes(radius_coefficients, radius_modes, quadrature, scale_mode)
        return [item.real for item in values]

    radius_grid = derive(0)
    slope = derive(1)
    second = derive(2)
    q_square = gauge ** 2
    numerator = []
    for node in range(quadrature):
        numerator.append(
            2 * radius_grid[node] * second[node]
            - slope[node] ** 2
            - q_square * radius_grid[node] ** 2
            + capital_g[node]
        )
    numerator_spectrum = _interval_fft([_iv.mpc(item) for item in numerator])
    coefficient_scale = _iv.mpf(1) / quadrature
    numerator_coefficients = [item * coefficient_scale for item in numerator_spectrum]
    g_coefficients = [
        item * coefficient_scale for item in _interval_fft([_iv.mpc(item) for item in capital_g])
    ]
    modes = np.fft.fftfreq(quadrature) * quadrature
    l2_lower = _iv.mpf(0)
    l2_upper = _iv.mpf(0)
    tail_upper = _iv.mpf(0)
    for coefficient, mode in zip(numerator_coefficients, modes):
        lower, upper = _abs2_bounds(coefficient)
        l2_lower += lower
        l2_upper += upper
        if abs(int(mode)) > caps["degree_P"]:
            tail_upper += upper
    length_iv = _iv.mpf(length)
    numerator_l2_lower = _iv.sqrt(length_iv * l2_lower)
    numerator_l2_upper = _iv.sqrt(length_iv * l2_upper)
    tail_l2_upper = _iv.sqrt(length_iv * tail_upper)

    def lipschitz(coefficients, mode_list):
        total = _iv.mpf(0)
        for coefficient, mode in zip(coefficients, mode_list):
            _lower, upper = _abs2_bounds(coefficient)
            total += abs(_iv.mpf(int(mode))) * wavenumber * _iv.sqrt(upper)
        return total

    radius_lip = lipschitz(radius_coefficients, radius_modes)
    g_lip = lipschitz(g_coefficients, modes)

    def continuum_range(coefficients, mode_list, derivative_lip):
        values = _place_modes(
            coefficients,
            mode_list,
            LAPSE_FINE_COUNT,
            lambda _mode: _iv.mpc(1),
        )
        sample_lower = min(item.real.a for item in values)
        sample_upper = max(item.real.b for item in values)
        half_step = length_iv / LAPSE_FINE_COUNT / 2
        return sample_lower - half_step * derivative_lip, sample_upper + half_step * derivative_lip

    radius_lower, radius_upper = continuum_range(radius_coefficients, radius_modes, radius_lip)
    g_lower, g_upper = continuum_range(g_coefficients, modes, g_lip)
    positive_radius = float(radius_lower.a) > 0.0
    positive_g = float(g_lower.a) > 0.0
    rho_float = bandlimited_matter_density(phi0, phi1, length, quadrature, float(gauge_samples[0]))
    rho_gap = max(abs(float(value.mid) - exact) for value, exact in zip(rho, rho_float))
    float_constants = {
        "Q": float(constants["Q"]),
        "A": float(constants["A"]),
        "rmag2": float(constants["rmag2"]),
        "eight_pi_A": 8.0 * math.pi * float(constants["A"]),
    }
    float_g = bracket_G(rho_float, float_constants)
    radius_coeff_float = np.fft.fft(radius) / geometry
    float_wavenumber = 2.0 * math.pi / length
    float_radius = _embed_normalized(radius_coeff_float, radius_modes, quadrature)
    float_slope = _embed_normalized(
        radius_coeff_float * (1j * radius_modes * float_wavenumber),
        radius_modes,
        quadrature,
    )
    float_second = _embed_normalized(
        radius_coeff_float * (-(radius_modes * float_wavenumber) ** 2),
        radius_modes,
        quadrature,
    )
    float_numerator = (
        2.0 * float_radius * float_second
        - float_slope ** 2
        - float(constants["Q"]) ** 2 * float_radius ** 2
        + float_g
    )
    float_l2 = math.sqrt(length * float(np.sum(np.abs(np.fft.fft(float_numerator) / quadrature) ** 2)))
    slope_gap = max(abs(float(value.mid) - exact) for value, exact in zip(slope, float_slope))
    radius_gap = max(abs(float(value.mid) - exact) for value, exact in zip(radius_grid, float_radius))
    effects = _areal_radius_effects()
    diagnostics = {
        "rho_float_gap_max": float(rho_gap),
        "radius_sample_float_gap_max": float(radius_gap),
        "radius_slope_float_gap_max": float(slope_gap),
        "numerator_l2_float_parseval": float_l2,
        "bilinear_imaginary_upper": _float_hi(imaginary_upper),
        "numerator_width_max": max(_float_hi(item.delta) for item in numerator),
        "identity_gate_is_not_a_geometry_tolerance": True,
        "arithmetic": "mpmath.iv outward interval arithmetic",
        "endpoint_rounding": "math.nextafter outward from the mpmath endpoint",
        "fft_or_refinement_difference_is_not_the_certificate": True,
    }
    if rho_gap > 1e-6 or radius_gap > 1e-6 or slope_gap > 1e-6 or float_l2 > _float_hi(numerator_l2_upper):
        return _invalid_lapse_record(
            "interval polynomial does not cover the independent float evaluation",
            started=started,
            diagnostics=diagnostics,
            status="IMPLEMENTATION_IDENTITY_FAILED",
        )
    if not positive_radius or not positive_g:
        record = _invalid_lapse_record(
            "continuum positivity of r and G is the missing primitive",
            started=started,
            diagnostics=diagnostics,
        )
        record["radius_lower"] = None if not positive_radius else _float_lo(radius_lower)
        record["G_lower"] = None if not positive_g else _float_lo(g_lower)
        return record
    y_min = _iv.sqrt(radius_lower)
    y_max = _iv.sqrt(radius_upper)
    mu = min(_iv.mpf(4), gauge ** 2)
    residual_l2 = numerator_l2_upper / (y_min ** 3)
    energy = residual_l2 / mu
    sobolev = _sobolev_factor(length)
    y_infinity = energy * sobolev
    if float(y_infinity.b) >= float(y_min.a):
        pointwise = None
    else:
        pointwise = y_infinity * (2 * y_max + y_infinity)
    numerator_l2_hi = _float_hi(numerator_l2_upper)
    numerator_l2_lo = _float_lo(numerator_l2_lower)
    energy_hi = _float_hi(energy)
    pointwise_hi = None if pointwise is None else _float_hi(pointwise)
    smaller_effect = None
    resolves = None
    effect_report = None
    if effects is not None:
        smaller_effect = min(effects["areal_radius_min_effect"], effects["areal_radius_max_effect"])
        bar = effects["bar_fraction"] * smaller_effect
        resolves = None if pointwise_hi is None else bool(pointwise_hi < bar)
        effect_report = dict(effects)
        effect_report["smaller_effect_scale"] = smaller_effect
        effect_report["one_percent_bar"] = bar
        effect_report["pointwise_radius_over_effect"] = (
            None if pointwise_hi is None else pointwise_hi / smaller_effect
        )
        effect_report["pointwise_radius_over_one_percent"] = (
            None if pointwise_hi is None else pointwise_hi / bar
        )
        effect_report["resolves_areal_radius_effect"] = resolves
        effect_report["comparison_is_initial_chart_only"] = True
    cpu = time.process_time() - started
    within_budget = cpu <= LAPSE_BOUND_CPU_S
    bound = energy_hi if within_budget else None
    return {
        "schema": ERROR_SCHEMA,
        "status": LAPSE_COMPONENT_STATUS if within_budget else "MISSING_PRIMITIVE",
        "reason": None if within_budget else "numerical CPU budget exceeded before a bound was accepted",
        "controlled_uncertainty": (
            "energy distance of y=sqrt(r) to a positive critical point of the "
            "initial lapse functional J, and the induced pointwise areal-radius distance"
        ),
        "geometry_error_upper_bound": bound,
        "pointwise_radius_error_upper_bound": pointwise_hi if within_budget else None,
        "pointwise_y_error_upper_bound": _float_hi(y_infinity) if within_budget else None,
        "shift_constraint_error_bound": None,
        "evolution_error_bound": None,
        "dense_operator_rounding_certificate": None,
        "minimizer_existence_proved": None,
        "minimizer_constructed": None,
        "total_certified": False,
        "continuum_solution_claimed": False,
        "mean_current_removed": False,
        "represented_weak_certifies_total": False,
        "cpu_seconds": cpu,
        "cpu_budget_seconds": LAPSE_BOUND_CPU_S,
        "cpu_budget_exceeded": not within_budget,
        "precision_dps": LAPSE_BOUND_DPS,
        "fine_positivity_nodes": LAPSE_FINE_COUNT,
        "settings": {
            "nf": fermions,
            "ng": geometry,
            "nq": quadrature,
            "length": length,
            "occupations": [float(value) for value in OCCUPATIONS],
            "kappa": int(KAPPA),
            "multiplicity": int(sector_multiplicity(KAPPA)),
            "Q": float(constants["Q"]),
            "A": float(constants["A"]),
            "rmag2": float(constants["rmag2"]),
            "chi_max": chi_max,
            "p_r_max": p_r_max,
            "p_chi_max": p_chi_max,
            "p_Q_max": p_q_max,
        },
        "inputs": {
            "v5_npz_sha256": _sha256(V5_NPZ),
            "v5_json_sha256": _sha256(V5_JSON),
            "sealed_weak_sha256": _sha256(RECORD_PATH),
            "nf512_phi0_sha256": phi0_hash,
            "nf512_phi1_sha256": phi1_hash,
            "nf512_radius_sha256": radius_hash,
            "nf512_geometry_Q_sha256": gauge_hash,
            "source_hashes": _lapse_source_hashes(),
            "mpmath_version": __import__("mpmath").__version__,
        },
        "hypotheses": {
            "Q_constant": True,
            "chi_pr_pchi_zero": True,
            "p_Q_zero_in_the_lapse_functional": True,
            "phi_is_the_stored_antiperiodic_interpolant": True,
            "momentum_symbol": "P=-i d/dx, mode m multiplied by 2*pi*m/length",
            "force_L": "multiplicity * (K/Q + kappa * real S1); rho = force_L/dx",
            "force_L_does_not_read_r": True,
            "G_formula": "Q^2*rmag2 + Q*rho/(8*pi*A)",
            "G_independent_of_y": True,
            "rho_independence_numerical_flag_used": False,
            "P_formula": "2*r*r''-(r')^2-Q^2*r^2+G = -y^3*Ry for y=sqrt(r)>0",
            "degree_cap": caps,
            "alias_free_quadrature": True,
            "continuum_r_positive": positive_radius,
            "continuum_G_positive": positive_g,
            "mu": "min(4, Q^2), from H>=4||v'||^2+Q^2||v||^2 on the positive cone",
            "dual_majorant": "||Ry||_* <= ||Ry||_L2 <= ||P||_L2 / y_min^3",
            "error_inequality": "||y-y_*||_E <= ||Ry||_* / mu if a positive critical point exists",
            "existence_not_proved": True,
            "shift_constraint_not_in_J": True,
            "pi": "interval pi in the certificate; math.pi remains the owner float path",
        },
        "degree": caps,
        "numerator_l2_lower": numerator_l2_lo,
        "numerator_l2_upper": numerator_l2_hi,
        "numerator_l2_directed_gap": numerator_l2_hi - numerator_l2_lo,
        "numerator_tail_above_degree_l2_upper": _float_hi(tail_l2_upper),
        "radius_lower": _float_lo(radius_lower),
        "radius_upper": _float_hi(radius_upper),
        "G_lower": _float_lo(g_lower),
        "G_upper": _float_hi(g_upper),
        "y_min": _float_lo(y_min),
        "y_max": _float_hi(y_max),
        "mu": _float_lo(mu),
        "sobolev_factor_upper": _float_hi(sobolev),
        "residual_l2_upper": _float_hi(residual_l2) if within_budget else None,
        "radius_derivative_l1_upper": _float_hi(radius_lip),
        "G_derivative_l1_upper": _float_hi(g_lip),
        "effect_comparison": effect_report,
        "diagnostics": diagnostics,
        "v5_bytes_unchanged": True,
        "sealed_weak_record_unchanged": True,
        "evolution_performed": False,
        "radius_newton_rerun": False,
        "source_formula_changed": False,
    }


def _embed_normalized(coefficients, modes, count):
    """Samples of ``Σ c_k e^{2π i k x/L}`` from normalized coefficients."""
    spectrum = np.zeros(count, dtype=complex)
    for coefficient, mode in zip(np.asarray(coefficients), np.asarray(modes)):
        mode = int(mode)
        slot = mode if mode >= 0 else count + mode
        spectrum[slot] = coefficient * count
    return np.fft.ifft(spectrum).real


def write_initial_lapse_record(payload, path):
    """Write a new lapse record. Refuses every sealed record, including error-v1."""
    destination = Path(path)
    forbidden = {
        RECORD_PATH.resolve(),
        ERROR_RECORD_PATH.resolve(),
        V5_JSON.resolve(),
        V5_NPZ.resolve(),
    }
    if destination.resolve() in forbidden:
        raise FileExistsError(f"{destination} is sealed and is not rewritten")
    if payload.get("schema") != ERROR_SCHEMA:
        raise ValueError("lapse record schema required")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(_jsonable(payload), indent=2) + "\n")
    return destination


def barrier_increment_bound(length):
    """Cauchy-Schwarz on an arc: |y(x)-y(z)|^2 is at most arc length times ||y'||_L2^2.

    The distance is the shorter arc, hence at most ``length/2``. Equality holds
    for a function whose derivative is constant on that arc.
    """
    length = float(length)
    if not math.isfinite(length) or length <= 0.0:
        raise ValueError("positive period required")
    return {
        "inequality": "|y(x)-y(z)|^2 <= arc(x,z) * ||y'||_L2^2",
        "max_arc": length / 2.0,
        "equality_case": "y' constant on the arc",
    }


def inverse_square_log_divergence(norm_squared, inner, outer):
    """Lower bound ``∫_inner^outer dx /(M x) = log(outer/inner)/M``."""
    if norm_squared <= 0.0 or inner <= 0.0 or outer <= inner:
        raise ValueError("positive derivative norm and 0 < inner < outer required")
    return math.log(outer / inner) / norm_squared


def lapse_existence_theorem(*, length, q_squared, g_lower, g_trigonometric=True):
    """Direct-method existence for the periodic lapse functional.

    No coefficient is recomputed. The argument uses only the stated hypotheses.
    """
    length = float(length)
    q_squared = float(q_squared)
    g_lower = float(g_lower)
    hypotheses_hold = bool(
        math.isfinite(length) and length > 0.0
        and math.isfinite(q_squared) and q_squared > 0.0
        and math.isfinite(g_lower) and g_lower > 0.0
    )
    steps = (
        "Domain A is the open set of H^1(R/LZ) functions with positive continuous minimum.",
        "J >= 2||y'||_L2^2 + (Q^2/2)||y||_L2^2, and the constant y=1 is a finite sublevel.",
        "Cauchy-Schwarz gives |y(x)-y(z)|^2 <= arc(x,z)||y'||_L2^2, so sublevels are equicontinuous.",
        "The Fourier embedding ||y||_inf <= sqrt(coth(L/2)/2) ||y||_E bounds the sublevel in C^0.",
        "Arzela-Ascoli plus weak compactness of the Hilbert ball produce y_n -> y uniformly and y_n weak to y in H^1.",
        "If y(x0)=0 then y(x)^2 <= |x-x0| ||y'||_L2^2, so the integral of 1/y^2 diverges logarithmically.",
        "Fatou's lemma on the nonnegative densities G/(2 y_n^2) sends a minimizing sequence to infinity if the limit touches zero.",
        "The limit therefore stays in A. Hilbert weak lower semicontinuity of ||y'|| and uniform convergence of the potential give a minimizer.",
        "t |-> t^2 and t |-> t^{-2} are strictly convex on (0, infinity); 6 t^{-4} is the second derivative of t^{-2}.",
        "Q^2>0 and G>=delta>0 make J strictly convex and C^1 on A, so the minimizer is the unique positive critical point.",
        "The weak equation gives y''=(Q^2 y - G/y^3)/4. Trigonometric G is C^infinity, and the Euler equation bootstraps to y in C^infinity.",
    )
    return {
        "domain": "H^1(R/LZ), positive continuous representative, periodic, no boundary term",
        "functional": "J=integral[2 (y')^2 + Q^2 y^2/2 + G/(2 y^2)] dx",
        "hypotheses_hold": hypotheses_hold,
        "g_lower": g_lower,
        "q_squared": q_squared,
        "length": length,
        "g_trigonometric": bool(g_trigonometric),
        "minimizer_exists": hypotheses_hold,
        "unique_positive_critical_point": hypotheses_hold,
        "smooth_critical_point": bool(hypotheses_hold and g_trigonometric),
        "evolution_error_bound": None,
        "initial_propagation_limited": True,
        "steps": list(steps),
        "external_lemmas": ("Fatou nonnegative", "Hilbert weak lower semicontinuity", "Arzela-Ascoli"),
    }


def saved_lapse_existence():
    """Apply the theorem to the sealed positivity enclosure. Does not rewrite it."""
    if not ERROR_RECORD_PATH.is_file():
        raise FileNotFoundError(ERROR_RECORD_PATH)
    sealed_bytes = ERROR_RECORD_PATH.read_bytes()
    sealed = json.loads(sealed_bytes)
    settings = sealed.get("settings") or {}
    theorem = lapse_existence_theorem(
        length=float(settings.get("length", 0.0)),
        q_squared=float(settings.get("Q", 0.0)) ** 2,
        g_lower=float(sealed.get("G_lower") or 0.0),
        g_trigonometric=sealed.get("G_lower") is not None,
    )
    if ERROR_RECORD_PATH.read_bytes() != sealed_bytes:
        raise RuntimeError("existence check rewrote the sealed lapse record")
    theorem.update({
        "sealed_record": "results/development/nsc-spherical-cauchy-error-v1.json",
        "sealed_sha256": hashlib.sha256(sealed_bytes).hexdigest(),
        "sealed_minimizer_existence_proved_field": sealed.get("minimizer_existence_proved"),
        "historical_field_left_null": sealed.get("minimizer_existence_proved") is None,
        "recomputed_numerator": False,
        "numerator_enclosure_reused": True,
    })
    return theorem


def _exact_unit_modulus_defect(real_part, imag_part):
    """Exact ``1-(a^2+b^2)`` for the binary values of a stored phase factor."""
    from fractions import Fraction

    real_ratio = Fraction(*float(real_part).as_integer_ratio())
    imag_ratio = Fraction(*float(imag_part).as_integer_ratio())
    return 1 - (real_ratio * real_ratio + imag_ratio * imag_ratio)


def _interval_abs_sum(coefficients):
    total = _iv.mpf(0)
    for coefficient in coefficients:
        _lower, upper = _abs2_bounds(coefficient)
        total += _iv.sqrt(upper)
    return total


def _interval_weighted_square_sum(coefficients, weights):
    total = _iv.mpf(0)
    for coefficient, weight in zip(coefficients, weights):
        lower, upper = _abs2_bounds(coefficient)
        span = _iv.mpf([_float_lo(lower), _float_hi(upper)])
        total += weight * span
    return total


def saved_shift_constraint_bound():
    """Enclose the saved shift current and the fixed-profile momentum correction.

    The stored columns satisfy ``partner = ω conj(plus)`` with equal pair
    occupations. The current is then exactly the factor ``1-|ω|^2`` times the
    plus-column momentum density. An interval containing zero is not reported
    as a zero mean. The source array is not modified and its mean is not deleted.
    """
    started = time.process_time()
    previous = _iv.dps
    _iv.dps = 40
    sealed_before = ERROR_RECORD_PATH.read_bytes()
    weak_before = RECORD_PATH.read_bytes()
    v5_before = V5_NPZ.read_bytes()
    try:
        return _saved_shift_constraint_bound(started, sealed_before)
    finally:
        _iv.dps = previous
        if ERROR_RECORD_PATH.read_bytes() != sealed_before or RECORD_PATH.read_bytes() != weak_before:
            raise RuntimeError("shift enclosure rewrote a sealed record")
        if V5_NPZ.read_bytes() != v5_before:
            raise RuntimeError("shift enclosure rewrote v5 bytes")


def _saved_shift_constraint_bound(started, sealed_before):
    sealed = json.loads(sealed_before)
    with np.load(V5_NPZ) as payload:
        phi0 = np.array(payload["nf512_columns_phi0"], copy=True)
        phi1 = np.array(payload["nf512_columns_phi1"], copy=True)
        radius = np.array(payload["nf512_radius"], dtype=float, copy=True)
        length = float(payload["nf512_length"])
        gauge = float(payload["nf512_geometry_Q"][0])
    phase = float(CALIBRATION["phase"])
    omega = np.exp(1j * np.float64(phase))
    symmetry_defect = 0.0
    for spinor in (phi0, phi1):
        for region in range(3):
            defect = spinor[:, 2 * region + 1] - omega * np.conjugate(spinor[:, 2 * region])
            symmetry_defect = max(symmetry_defect, float(np.max(np.abs(defect))))
    pair_occupations = all(
        float(OCCUPATIONS[2 * region]) == float(OCCUPATIONS[2 * region + 1]) for region in range(3)
    )
    factor = _exact_unit_modulus_defect(omega.real, omega.imag)
    if symmetry_defect != 0.0 or not pair_occupations or factor == 0:
        cpu = time.process_time() - started
        return {
            "status": "SHIFT_IDENTITY_UNAVAILABLE",
            "reason": "saved columns do not have the exact conjugate pair symmetry",
            "symmetry_defect_max": symmetry_defect,
            "modulus_defect": str(factor),
            "mean_contains_zero": None,
            "shift_residual_upper_bound": None,
            "source_modified": False,
            "mean_deleted": False,
            "cpu_seconds": cpu,
            "evolution_error_bound": None,
        }
    factor_iv = _iv.mpf(factor.numerator) / _iv.mpf(factor.denominator)
    sqrt_dx = _iv.sqrt(_iv.mpf(length) / phi0.shape[0])
    wavenumber = 2 * _iv.pi / _iv.mpf(length)
    multiplicity = _iv.mpf(sector_multiplicity(KAPPA))
    integers = np.fft.fftfreq(phi0.shape[0]) * phi0.shape[0]
    symbols = [wavenumber * (_iv.mpf(int(mode)) + _iv.mpf("0.5")) for mode in integers]
    index = np.arange(phi0.shape[0])
    demod = np.exp(-1j * math.pi * index / phi0.shape[0])
    scale = _iv.mpf(1) / phi0.shape[0]
    moment = _iv.mpf(0)
    infinity_factor = _iv.mpf(0)
    for region in range(3):
        occupation = _iv.mpf(float(OCCUPATIONS[2 * region]))
        column_moment = _iv.mpf(0)
        column_product = _iv.mpf(0)
        for spinor in (phi0, phi1):
            samples = spinor[:, 2 * region] * demod
            coefficients = [
                item * scale / sqrt_dx
                for item in _interval_fft([_as_mpc(value) for value in samples])
            ]
            column_moment += _interval_weighted_square_sum(coefficients, symbols)
            value_sum = _interval_abs_sum(coefficients)
            momentum_sum = _interval_abs_sum([
                coefficient * symbol for coefficient, symbol in zip(coefficients, symbols)
            ])
            column_product += value_sum * momentum_sum
        moment += occupation * column_moment
        infinity_factor += occupation * column_product
    mean_j = -multiplicity * factor_iv * moment
    infinity = multiplicity * abs(factor_iv) * infinity_factor
    radius_coefficients = [
        item / radius.size for item in _interval_dft([_as_mpc(value) for value in radius])
    ]
    radius_modes = np.fft.fftfreq(radius.size) * radius.size
    slope_energy = _iv.mpf(0)
    slope_lip = _iv.mpf(0)
    for coefficient, mode in zip(radius_coefficients, radius_modes):
        weight = wavenumber * _iv.mpf(int(mode))
        lower, upper = _abs2_bounds(coefficient)
        span = _iv.mpf([_float_lo(lower), _float_hi(upper)])
        slope_energy += weight * weight * span
        slope_lip += abs(weight) * _iv.sqrt(upper)
    slope_energy *= _iv.mpf(length)
    zee = -24 * _iv.pi * _iv.mpf(float(sealed["settings"]["A"]))
    q_value = _iv.mpf(gauge)
    denominator = abs(4 * zee * q_value)
    lam = -(mean_j * _iv.mpf(length)) / slope_energy
    delta_c = (lam * lam) * (slope_lip * slope_lip) / denominator
    energy_upper = _iv.mpf(float(sealed["geometry_error_upper_bound"]))
    y_min = _iv.mpf(float(sealed["y_min"]))
    y_max = _iv.mpf(float(sealed["y_max"]))
    slope_l2 = _iv.sqrt(slope_energy)
    delta_r_prime = (slope_lip / y_min) * energy_upper + 2 * y_max * energy_upper
    leftover = abs(lam) * slope_l2 * delta_r_prime
    one_percent = float(sealed["effect_comparison"]["one_percent_bar"])
    mean_lower = _float_lo(mean_j)
    mean_upper = _float_hi(mean_j)
    contains_zero = mean_lower <= 0.0 <= mean_upper
    infinity_upper = _float_hi(infinity)
    cpu = time.process_time() - started
    return {
        "status": "SHIFT_MEAN_ENCLOSED",
        "primitive": "periodic shift constraint for chi=p_r=p_chi=0 and constant Q",
        "identity": (
            "partner = omega conj(plus), equal pair occupations, "
            "j = -M (1-|omega|^2) Re(conj(psi_plus) P psi_plus)"
        ),
        "symmetry_defect_max": symmetry_defect,
        "pair_occupations_equal": True,
        "modulus_defect_numerator": int(factor.numerator),
        "modulus_defect_denominator": int(factor.denominator),
        "mean_j_lower": mean_lower,
        "mean_j_upper": mean_upper,
        "mean_contains_zero": contains_zero,
        "mean_exactly_zero": False,
        "mean_deleted": False,
        "source_modified": False,
        "current_infinity_upper": infinity_upper,
        "geometry_band_tail_l2_upper": math.sqrt(length) * infinity_upper,
        "pr_zero_residual_after_periodic_antiderivative": (
            "the constant mean(j); a periodic derivative cannot cancel it"
        ),
        "correction": {
            "form": "p_r = lambda * derivative(saved r), chi = p_chi = 0, Q constant, j unchanged",
            "admissible": True,
            "from_the_existing_action": True,
            "lambda_lower": _float_lo(lam),
            "lambda_upper": _float_hi(lam),
            "delta_C_formula": "p_r^2 / (4 Z Q), Z = feedback_Z(A)",
            "delta_C_infinity_upper": _float_hi(delta_c),
            "continuum_shift_residual_after_correction": 0.0,
            "critical_point_mean_leftover_upper": _float_hi(leftover),
            "source_modified": False,
            "mean_deleted": False,
        },
        "one_percent_bar": one_percent,
        "current_infinity_over_one_percent": infinity_upper / one_percent,
        "resolves_initial_shift_against_areal_bar": bool(infinity_upper < one_percent),
        "evolution_error_bound": None,
        "initial_propagation_limited": True,
        "cpu_seconds": cpu,
        "cpu_budget_seconds": 30.0,
        "recomputed_lapse_numerator": False,
    }


REGENERATION_NPZ = _LAB_ROOT / "results" / "development" / "nsc-regeneration-episode-v1.npz"
REGENERATION_CASES = (
    "nf256_dtmax_0_0005",
    "nf256_dtmax_0_00025",
    "nf512_dtmax_0_00025",
    "nf512_dtmax_0_0005",
)
_LAMBDA_PROBE = 2e-15
_MISSING_OBSERVABLE_PRIMITIVE = (
    "the map from the source-normalized constraint budget and a stored "
    "Hermitian geometry mismatch delta H(t) to the physical observables; "
    "a constraint norm is not that map"
)


def matter_shift_to_g(pr_square, A, Z):
    """Convert ``p_r**2/(4 Z Q)`` in the lapse constraint into an addition to ``G``.

    ``C`` gains ``p_r**2/(4 Z Q)``. The owner identity is
    ``C = (8 pi A/Q)(P_geom + G)``, so the same addition is
    ``G -> G + p_r**2/(32 pi A Z)``. ``Z = -24 pi A`` is negative.
    """
    return pr_square / (32.0 * math.pi * float(A) * float(Z))


def geometric_constraint_transport():
    """Off-shell transport of the geometric constraints under ``H = ∫(L C + β D)``.

    The identity is the Hamilton flow of the owned first-order action, with
    lapse and shift prescribed. It does not include the Dirac force. The
    principal coupling ``(L/Q**2) D_x`` does not cancel in the ``L^2`` energy
    of ``(C, D)``.
    """
    import sympy as sp
    from .nsc_spherical_feedback_action import (
        JetBundle,
        hamiltonian_right_hand_sides,
        lapse_constraint,
        shift_constraint,
    )

    raw = hamiltonian_right_hand_sides()
    jets = JetBundle(("Q", "r", "chi", "L", "beta", "pQ", "pr", "pchi"), max_order=3)
    flat = {}
    for table in jets.table.values():
        for symbol in table.values():
            flat[symbol.name] = symbol

    def on_jets(expr):
        mapping = {
            symbol: flat[symbol.name]
            for symbol in expr.free_symbols if symbol.name in flat
        }
        return expr.xreplace(mapping)

    velocities = {
        field: on_jets(raw[key]) for field, key in (
            ("Q", "Q_t"), ("r", "r_t"), ("chi", "chi_t"),
            ("pQ", "pQ_t"), ("pr", "pr_t"), ("pchi", "pchi_t"),
        )
    }
    A, C_W, C_F, flux = sp.symbols("A C_W C_F Phi", nonzero=True)
    C = lapse_constraint(
        jets("pQ"), jets("pr"), jets("pchi"), jets("Q"), jets("r"), jets("chi"),
        jets("r", 0, 1), jets("r", 0, 2), jets("chi", 0, 1), jets("chi", 0, 2),
        jets("Q", 0, 1), A, C_W, C_F, flux,
    )
    D = shift_constraint(
        jets("pQ", 0, 1), jets("pr"), jets("pchi"), jets("Q"),
        jets("r", 0, 1), jets("chi", 0, 1),
    )

    def flow(expr):
        acc = sp.Integer(0)
        for field, velocity in velocities.items():
            for (nt, nx), symbol in jets.table[field].items():
                if nt != 0 or symbol not in expr.free_symbols:
                    continue
                derivative = velocity
                for _ in range(nx):
                    derivative = jets.total(derivative, "x")
                acc += sp.diff(expr, symbol) * derivative
        return sp.expand(acc)

    C_t = flow(C)
    D_t = flow(D)
    beta = jets("beta")
    betax = jets("beta", 0, 1)
    L = jets("L")
    Q = jets("Q")
    Lx = jets("L", 0, 1)
    Qx = jets("Q", 0, 1)
    guess_C = (
        beta * jets.total(C, "x") + betax * C
        + (L / Q**2) * jets.total(D, "x")
        + (2 * Lx / Q**2 - 2 * L * Qx / Q**3) * D
    )
    guess_D = beta * jets.total(D, "x") + 2 * betax * D + Lx * C
    force, force_x = sp.symbols("F_Q F_Qx")
    extra_C = sp.factor(sp.diff(C, jets("pQ")) * (-force))
    extra_D = sp.factor(sp.diff(D, jets("pQ", 0, 1)) * (-force_x))
    return {
        "geometric_identity": bool(sp.expand(C_t - guess_C) == 0 and sp.expand(D_t - guess_D) == 0),
        "momentum_equations_match_euler": bool(raw["momentum_rhs_matches_euler"]),
        "C_t": "beta*C_x + beta_x*C + (L/Q**2)*D_x + (2*L_x/Q**2 - 2*L*Q_x/Q**3)*D",
        "D_t": "beta*D_x + 2*beta_x*D + L_x*C",
        "includes_dirac_force": False,
        "force_Q_extra_C": str(sp.factor(extra_C)),
        "force_Q_extra_D": str(sp.factor(extra_D)),
        "principal_term_in_C_t": "(L/Q**2)*D_x",
        "D_t_contains_C_x": False,
        "l2_energy_of_C_and_D_closes": False,
        "l2_energy_of_C_and_D_over_Q_closes": False,
        "energy_obstruction": (
            "(L/Q**2)*partial_x D survives in both the L2 product of (C, D) "
            "and the L2 product of (C, D/Q)"
        ),
        "finite_band_transport": (
            "On a trigonometric band of degree N, ||D_x||_2 <= (2*pi*N/L) ||D||_2, "
            "so that L2 energy closes with a coefficient proportional to N. "
            "The saved series stores maxima, not those coefficients. "
            "The resulting estimate would still be a constraint budget, not an observable error."
        ),
        "controls_physical_observables": False,
        "missing_observable_primitive": _MISSING_OBSERVABLE_PRIMITIVE,
        "evolution_error_bound": None,
    }


def _saved_radius_slope_bounds():
    """Interval Parseval bounds for the saved radius derivative."""
    with np.load(V5_NPZ) as payload:
        radius = np.array(payload["nf512_radius"], dtype=float, copy=True)
        length = float(payload["nf512_length"])
    coefficients = [item / radius.size for item in _interval_dft([_as_mpc(value) for value in radius])]
    modes = np.fft.fftfreq(radius.size) * radius.size
    wavenumber = 2 * _iv.pi / _iv.mpf(length)
    energy = _iv.mpf(0)
    lip = _iv.mpf(0)
    for coefficient, mode in zip(coefficients, modes):
        weight = wavenumber * _iv.mpf(int(mode))
        lower, upper = _abs2_bounds(coefficient)
        span = _iv.mpf([_float_lo(lower), _float_hi(upper)])
        energy += weight * weight * span
        lip += abs(weight) * _iv.sqrt(upper)
    energy *= _iv.mpf(length)
    return _float_lo(energy), _float_hi(energy), _float_hi(lip), length


def nearby_exact_initial_state(lambda_probe=_LAMBDA_PROBE):
    """Exact nearby initial data for the unchanged source, by a scalar IVT.

    The saved state still has ``p_r = 0``. This does not reset that record.
    """
    started = time.process_time()
    previous = _iv.dps
    _iv.dps = 40
    sealed_before = ERROR_RECORD_PATH.read_bytes()
    weak_before = RECORD_PATH.read_bytes()
    v5_before = V5_NPZ.read_bytes()
    try:
        shift = saved_shift_constraint_bound()
        sealed = json.loads(sealed_before)
        energy_lo, energy_hi, lip, length = _saved_radius_slope_bounds()
        report = _nearby_from_bounds(
            shift, sealed, energy_lo, energy_hi, lip, length, float(lambda_probe),
        )
        report["cpu_seconds"] = time.process_time() - started
        return report
    finally:
        _iv.dps = previous
        if (
            ERROR_RECORD_PATH.read_bytes() != sealed_before
            or RECORD_PATH.read_bytes() != weak_before
            or V5_NPZ.read_bytes() != v5_before
        ):
            raise RuntimeError("nearby-state proof rewrote a sealed record")


def _nearby_from_bounds(shift, sealed, energy_lo, energy_hi, lip, length, probe):
    # Enclosed constants describe the saved y, not the unknown critical y_0.
    # Transfer their range through the sealed H1/Sobolev distance before
    # evaluating the perturbation residual delta_G / y_0**3.
    g_lower = _iv.mpf(sealed["G_lower"])
    y_min = _iv.mpf(sealed["y_min"])
    y_max = _iv.mpf(sealed["y_max"])
    radius_error = _iv.mpf(sealed["pointwise_radius_error_upper_bound"])
    state_error = _iv.mpf(sealed["geometry_error_upper_bound"])
    gauge = _iv.mpf(sealed["settings"]["Q"])
    area = _iv.mpf(sealed["settings"]["A"])
    length_iv = _iv.mpf(length)
    probe_iv = _iv.mpf(probe)
    if probe <= 0.0 or energy_lo < 0.0 or energy_hi < energy_lo:
        raise ValueError("positive lambda probe and ordered nonnegative slope bounds required")
    sobolev = _sobolev_factor(length_iv)
    pointwise_y_error = sobolev * state_error
    critical_y_min = y_min - pointwise_y_error
    critical_y_max = y_max + pointwise_y_error
    if float(critical_y_min.a) <= 0.0:
        raise ValueError("initial distance does not certify a positive critical-point minimum")
    mu = min(_iv.mpf(4), gauge ** 2)
    denominator = _iv.mpf(768) * _iv.pi ** 2 * area ** 2
    delta_g = (probe_iv * _iv.mpf(lip)) ** 2 / denominator
    g_lambda_lower = g_lower - delta_g
    f_norm = _iv.sqrt(_iv.mpf(energy_hi))
    saved_y_prime_norm = f_norm / (2 * y_min)
    # r_0' - r_s' = 2 y_0 (y_0-y_s)' + 2 (y_0-y_s) y_s'.
    slope_factor = 2 * critical_y_max + 2 * sobolev * saved_y_prime_norm
    d0_lower = _iv.mpf(energy_lo) - f_norm * slope_factor * state_error
    delta_y = _iv.sqrt(length_iv) * delta_g / critical_y_min ** 3 / mu
    perturbed_y_max = critical_y_max + sobolev * delta_y
    critical_y_prime_norm = saved_y_prime_norm + state_error
    perturbed_slope_factor = 2 * perturbed_y_max + 2 * sobolev * critical_y_prime_norm
    theta = f_norm * perturbed_slope_factor * delta_y
    integral_lo = _iv.mpf(shift["mean_j_lower"]) * length_iv
    integral_hi = _iv.mpf(shift["mean_j_upper"]) * length_iv
    margin = d0_lower - theta
    h_plus_lower = integral_lo + probe_iv * margin
    h_minus_upper = integral_hi - probe_iv * margin
    root = bool(
        _float_lo(g_lambda_lower) > 0.0 and _float_lo(margin) > 0.0
        and _float_lo(h_plus_lower) > 0.0 and _float_hi(h_minus_upper) < 0.0
    )
    return {
        "status": "EXACT_NEARBY_INITIAL_STATE" if root else "IVT_HYPOTHESIS_FAILED",
        "not_actual_reset": True,
        "correction_installed_in_records": False,
        "actual_saved_pr": 0.0,
        "source_unchanged": True,
        "records_unchanged": True,
        "mean_deleted": False,
        "fixed_profile": "f = derivative of the saved radius",
        "G_lambda_formula": "G + p_r**2 / (32*pi*A*Z)",
        "Z_sign": "negative",
        "delta_g_upper": _float_hi(delta_g),
        "g_lambda_lower": _float_lo(g_lambda_lower),
        "critical_y_lower": _float_lo(critical_y_min),
        "critical_y_upper": _float_hi(critical_y_max),
        "arithmetic": "outward interval arithmetic through the IVT margins",
        "unique_critical_radius_for_each_probe_lambda": bool(_float_lo(g_lambda_lower) > 0.0),
        "D0_lower": _float_lo(d0_lower),
        "lambda_probe": probe,
        "H_plus_lower": _float_lo(h_plus_lower),
        "H_minus_upper": _float_hi(h_minus_upper),
        "root_in_open_interval": root,
        "periodic_p_Q_after_root": bool(root),
        "constraints_C_and_D_exact": bool(root),
        "radius_infinity_distance_upper": _float_hi(
            radius_error + sobolev * delta_y * (2 * critical_y_max + sobolev * delta_y)
        ),
        "pr_infinity_upper": _float_hi(probe_iv * _iv.mpf(lip)),
        "evolution_error_bound": None,
        "initial_propagation_limited": True,
    }


def source_constraint_normalization(force_l, force_beta, radius, radial_density, spacing):
    """Map coordinate constraint densities to the owned observer stresses.

    ``observer_stresses`` uses ``rho = force_L / (4 pi r^4 Q dx)`` and
    ``j = -force_beta / (4 pi r^4 Q^2 dx)``. The quadrature constraints are
    ``C_hat = C_geom + force_L/dx`` and ``D_hat = D_geom + force_beta/dx``,
    so the physical residuals are ``C_hat / (4 pi r^4 Q)`` and
    ``-D_hat / (4 pi r^4 Q^2)``. Their absolute integrals against the proper
    volume ``4 pi r^3 Q dx`` are ``integral |C_hat|/r dx`` and
    ``integral |D_hat|/(r Q) dx``.
    """
    radius = float(radius)
    radial_density = float(radial_density)
    spacing = float(spacing)
    if radius <= 0.0 or radial_density <= 0.0 or spacing == 0.0:
        raise ValueError("positive radius and Q and nonzero spacing required")
    sphere = 4.0 * math.pi * radius ** 4
    coordinate_l = float(force_l) / spacing
    coordinate_beta = float(force_beta) / spacing
    density = coordinate_l / (sphere * radial_density)
    current = -coordinate_beta / (sphere * radial_density ** 2)
    volume = 4.0 * math.pi * radius ** 3 * radial_density
    return {
        "physical_density": density,
        "physical_current": current,
        "normal_energy_sample": volume * density * spacing,
        "current_volume_sample": volume * abs(current) * spacing,
        "force_l_over_radius": float(force_l) / radius,
    }


def _static_gauge_integral(length, count):
    spacing = float(length) / int(count)
    nodes = (np.arange(int(count)) + 0.5) * spacing
    coordinate = profile_coordinate(nodes, length)
    scale = np.exp(np.log(float(CALIBRATION["Omega"])) * coordinate)
    lapse = float(CALIBRATION["b0"]) * scale
    shift = float(CALIBRATION["beta0"]) * scale
    return float(spacing * np.sum(np.abs(lapse))), float(spacing * np.sum(np.abs(shift)))


def static_constraint_smearing_integrals(length=PERIOD):
    """Upper indicators for ``integral |L|`` and ``integral |beta|``.

    Two midpoint grids bound the quadrature gap. This is the static gauge
    weight, not a constraint residual.
    """
    coarse_l, coarse_b = _static_gauge_integral(length, 4096)
    fine_l, fine_b = _static_gauge_integral(length, 8192)
    return {
        "length": float(length),
        "lapse_abs_integral_upper": max(coarse_l, fine_l) + abs(fine_l - coarse_l),
        "shift_abs_integral_upper": max(coarse_b, fine_b) + abs(fine_b - coarse_b),
        "quadrature_gap_lapse": abs(fine_l - coarse_l),
        "quadrature_gap_shift": abs(fine_b - coarse_b),
        "indicator_not_continuum_certificate": True,
    }


def _window_case_budget(payload, case, length, lapse_integral, shift_integral):
    prefix = case + "_"
    time_values = np.array(payload[prefix + "time"], dtype=float)
    hamilton = np.abs(np.array(payload[prefix + "full_hamilton_max"], dtype=float))
    momentum = np.abs(np.array(payload[prefix + "full_momentum_max"], dtype=float))
    radius = np.array(payload[prefix + "r_min"], dtype=float)
    radial = np.array(payload[prefix + "Q_min"], dtype=float)
    field = np.array(payload[prefix + "field_energy"], dtype=float)
    if np.min(radius) <= 0.0 or np.min(radial) <= 0.0:
        raise ValueError("stored radius and Q must stay positive")
    normal = hamilton * length / radius
    current = momentum * length / (radius * radial)
    smearing = hamilton * lapse_integral + momentum * shift_integral
    exchange = abs(float(field[-1] - field[0]))
    return {
        "case": case,
        "time_start": float(time_values[0]),
        "time_end": float(time_values[-1]),
        "samples": int(time_values.size),
        "full_hamilton_max": float(np.max(hamilton)),
        "held_out_hamilton_max": float(np.max(np.abs(payload[prefix + "held_out_hamilton_max"]))),
        "full_momentum_max": float(np.max(momentum)),
        "held_out_momentum_max": float(np.max(np.abs(payload[prefix + "held_out_momentum_max"]))),
        "ward_max": float(np.max(np.abs(payload[prefix + "ward_max"]))),
        "source_gap_max": float(np.max(np.abs(payload[prefix + "source_gap_max"]))),
        "unitarity_max": float(np.max(np.abs(payload[prefix + "unitarity"]))),
        "field_energy_exchange": exchange,
        "normal_energy_residual_majorant": float(np.max(normal)),
        "current_volume_residual_majorant": float(np.max(current)),
        "coordinate_smearing_majorant": float(np.max(smearing)),
        "smearing_over_field_exchange": None if exchange == 0.0 else float(np.max(smearing) / exchange),
        "majorant": (
            "recorded quadrature maximum times a positive weight from the stored "
            "r_min and Q_min; not a continuum maximum"
        ),
        "is_trajectory_error": False,
        "state_error_bound": None,
    }


def regeneration_window_coverage():
    """Source-normalized residual budgets for the four saved continuations.

    Raw constraint maxima are kept and are not compared with radius or
    proper-velocity changes. The budgets are not a state-error certificate.
    """
    if not REGENERATION_NPZ.is_file():
        return {
            "status": "WINDOW_UNAVAILABLE",
            "state_error_bound": None,
            "constraint_norm_is_state_accuracy": False,
        }
    weights = static_constraint_smearing_integrals()
    with np.load(REGENERATION_NPZ) as payload:
        cases = [
            _window_case_budget(
                payload, case, weights["length"],
                weights["lapse_abs_integral_upper"],
                weights["shift_abs_integral_upper"],
            )
            for case in REGENERATION_CASES
        ]
    return {
        "status": "SOURCE_RESIDUAL_BUDGET",
        "normalization": (
            "delta_rho = C_hat/(4*pi*r^4*Q), normal-energy integral |C_hat|/r dx; "
            "delta_j = -D_hat/(4*pi*r^4*Q^2), current-volume integral |D_hat|/(r*Q) dx; "
            "coordinate budget integral |L*C_hat+beta*D_hat|"
        ),
        "weights": weights,
        "cases": cases,
        "constraint_norm_is_state_accuracy": False,
        "state_error_bound": None,
        "missing_observable_primitive": _MISSING_OBSERVABLE_PRIMITIVE,
        "evolution_error_bound": None,
    }
