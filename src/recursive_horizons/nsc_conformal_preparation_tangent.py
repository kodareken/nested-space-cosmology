"""One paired-preparation initial tangent for the conformal rank-six slice.

The unitary is applied on the fermion nodes, then prolonged. The existing
pre-gauge projected radius Jacobian supplies the one dense factorization.
This is not a root solve, a time step, a retarded kernel, or a Lambda fit.
A small factorization residual is not a certificate of invertibility.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from scipy.linalg.lapack import get_lapack_funcs

from .nsc_spherical_coupling import (
    CALIBRATION,
    OCCUPATIONS,
    TOL_INITIAL_CONSTRAINT,
    TOL_NEWTON,
    shift_constraint,
    source_from_columns,
)
from .nsc_spherical_galerkin_coupling import (
    _check_state_dimensions,
    _constraint_directional,
    _projected_radius_jacobian,
    prolong_columns,
    prolong_geometry,
    prolong_state,
    projected_radius_operator,
    pull_geometry,
)


# Same absolute acceptance as preparation_problems. Not a continuum bound.
PAIR_GATE = 1e-8
# Initial chart values are stored zeros and the calibration ratio.
ZERO_GATE = 1e-12
# Prior-art numerical consistency of the projected linear solve.
RELATIVE_SOLVE_GATE = 1e-8
PAIRS = ((0, 1), (2, 3), (4, 5))


@dataclass
class PreparationTangent:
    """Coarse radius response and the nf spinor tangent. Other geometry variations are zero."""

    delta_r: np.ndarray
    phi0_tangent: np.ndarray
    phi1_tangent: np.ndarray
    report: dict


def preparation_generator(xi, length):
    """f = cos(2*pi*xi/length) on the nodes where the unitary is applied."""
    nodes = np.asarray(xi, dtype=float)
    period = float(length)
    if (nodes.ndim != 1 or not np.isfinite(period) or period <= 0.0
            or not np.isfinite(nodes).all()):
        raise ValueError("preparation generator needs finite nodes and a positive period")
    return np.cos(2.0 * np.pi * nodes / period)


def _generator_values(generator, count):
    values = np.asarray(generator, dtype=float)
    if values.shape != (count,) or not np.isfinite(values).all():
        raise ValueError("generator must be a finite real value on every fermion node")
    return values


def apply_preparation(phi0, phi1, generator, epsilon):
    """Finite U_eps = exp(i*eps*diag(f)*sigma2), the same factor on every column.

    phi0' = cos(eps*f)*phi0 + sin(eps*f)*phi1
    phi1' = cos(eps*f)*phi1 - sin(eps*f)*phi0
    """
    phi0 = np.asarray(phi0)
    phi1 = np.asarray(phi1)
    if phi0.shape != phi1.shape or phi0.ndim != 2:
        raise ValueError("spinor blocks must share the column shape")
    if not np.isfinite(epsilon):
        raise ValueError("preparation epsilon must be finite")
    angle = float(epsilon) * _generator_values(generator, phi0.shape[0])
    cosine = np.cos(angle)[:, None]
    sine = np.sin(angle)[:, None]
    return cosine * phi0 + sine * phi1, cosine * phi1 - sine * phi0


def spinor_tangent(phi0, phi1, generator):
    """Derivative of U_eps at eps = 0: (f*phi1, -f*phi0)."""
    phi0 = np.asarray(phi0)
    phi1 = np.asarray(phi1)
    if phi0.shape != phi1.shape or phi0.ndim != 2:
        raise ValueError("spinor blocks must share the column shape")
    factor = _generator_values(generator, phi0.shape[0])[:, None]
    return factor * phi1, -factor * phi0


def _pair_problems(phi0, phi1, occupations):
    problems = []
    if phi0.shape != phi1.shape or phi0.ndim != 2 or phi0.shape[1] != 6:
        return ["spinors must have shape (nf, 6)"]
    if not np.isfinite(phi0).all() or not np.isfinite(phi1).all():
        return ["spinor columns must be finite"]
    if not np.array_equal(np.asarray(occupations), OCCUPATIONS):
        problems.append("occupations are not the declared equal pairs 0.75, 0.5 and 0.25")
    recorded = np.exp(1j * float(CALIBRATION["phase"]))
    for region, (plus, minus) in enumerate(PAIRS):
        plus0, plus1 = phi0[:, plus], phi1[:, plus]
        minus0, minus1 = phi0[:, minus], phi1[:, minus]
        if float(np.max(np.abs(plus1 - 1j * plus0))) > PAIR_GATE:
            problems.append(f"pair {region} plus column does not have phi1 = i*phi0")
        if float(np.max(np.abs(minus0 - recorded * np.conjugate(plus0)))) > PAIR_GATE:
            problems.append(
                f"pair {region} minus phi0 is not the recorded phase times conjugate(plus)"
            )
        if float(np.max(np.abs(minus1 - recorded * np.conjugate(plus1)))) > PAIR_GATE:
            problems.append(
                f"pair {region} minus phi1 is not the recorded phase times conjugate(plus)"
            )
    return problems


def _geometry_problems(grid, state):
    problems = []
    q0 = float(CALIBRATION["b0"] / CALIBRATION["a0"])
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi"):
        value = np.asarray(getattr(state, name))
        if value.shape != (grid.ng,) or value.dtype.kind != "f" or not np.isfinite(value).all():
            problems.append(f"{name} must be a finite real array of shape (ng,)")
    if problems:
        return problems
    if float(np.max(np.abs(state.Q - q0))) > ZERO_GATE:
        problems.append("Q is not the constant calibration ratio b0/a0")
    if float(np.max(np.abs(state.chi))) > ZERO_GATE:
        problems.append("chi is not zero")
    for name in ("p_Q", "p_r", "p_chi"):
        if float(np.max(np.abs(getattr(state, name)))) > ZERO_GATE:
            problems.append(f"{name} is not zero")
    if float(np.min(state.r)) <= 0.0:
        problems.append("radius is not positive")
    return problems


def _fail(prefix, problems):
    if problems:
        raise ValueError(prefix + ": " + "; ".join(problems))


def _abs_max(values):
    array = np.asarray(values)
    if np.iscomplexobj(array):
        imaginary = float(np.max(np.abs(np.imag(array))))
        if imaginary > PAIR_GATE:
            raise ValueError(f"a real tangent acquired imaginary part {imaginary}")
        array = np.real(array)
    return float(np.max(np.abs(array)))


def _relative(absolute, scale):
    if scale == 0.0:
        return 0.0 if absolute == 0.0 else math.inf
    return absolute / scale


def _polarization(system, fine_state, source, radius_tangent, phi0_tangent, phi1_tangent):
    zero = np.zeros(fine_state.r.shape[0])
    tangent = (zero, radius_tangent, zero, zero, zero, zero, phi0_tangent, phi1_tangent)
    return _constraint_directional(system, fine_state, tangent, source)


def _constraint_density(h_dot, lapse):
    if np.iscomplexobj(h_dot):
        imaginary = float(np.max(np.abs(np.imag(np.asarray(h_dot)))))
        if imaginary > PAIR_GATE:
            raise ValueError(f"Hamilton polarization left the real line by {imaginary}")
        h_dot = np.real(h_dot)
    density = np.asarray(h_dot, dtype=float) / np.asarray(lapse, dtype=float)
    if density.shape != (lapse.shape[0],) or not np.isfinite(density).all():
        raise ValueError("constraint density tangent is not finite on the quadrature grid")
    return density


def _factor_once(matrix, rhs):
    """Factor one square real matrix and solve it. Do not retry or form an inverse."""
    matrix = np.asarray(matrix, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or rhs.shape != (matrix.shape[0],):
        raise ValueError("radius factorization needs one square Jacobian and one matching right-hand side")
    if not np.isfinite(matrix).all() or not np.isfinite(rhs).all():
        raise ValueError("radius factorization received a nonfinite matrix or right-hand side")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", LinAlgWarning)
        lu, piv = lu_factor(matrix)
    singular = [str(item.message) for item in caught if issubclass(item.category, LinAlgWarning)]
    u_min = float(np.min(np.abs(np.diag(lu))))
    if singular or u_min == 0.0:
        detail = singular[0] if singular else "U diagonal contains an exact zero"
        raise ValueError("singular projected radius Jacobian: " + detail)
    gecon, = get_lapack_funcs(("gecon",), (lu,))
    rcond, info = gecon(lu, np.linalg.norm(matrix, 1), norm="1")
    rcond, info = float(rcond), int(info)
    if info != 0 or not math.isfinite(rcond) or rcond == 0.0:
        raise ValueError(
            f"LAPACK rejected the projected radius Jacobian: gecon info {info}, rcond {rcond}"
        )
    # The residual against the original matrix is computed by the caller.
    # rcond and that residual do not certify invertibility.
    solution = lu_solve((lu, piv), rhs)
    return np.asarray(solution, dtype=float), {
        "rcond": rcond,
        "gecon_info": info,
        "u_diagonal_min_abs": u_min,
        "singular_warning": False,
        "certified_invertibility": False,
        "norm": "1",
        "factorizations": 1,
    }


def preparation_tangent(grid, state, *, checkpoint=None):
    """Differentiate the saved constraint level for the declared chiral family.

    The initial defect is measured and preserved, not solved away. This
    derivative does not require the saved constraint level to be zero.
    Premises and source cancellation precede the single factorization.
    An optional budget checkpoint runs between expensive operations.
    """
    if checkpoint is not None and not callable(checkpoint):
        raise ValueError("budget checkpoint must be callable")

    def check():
        if checkpoint is not None:
            checkpoint()

    check()
    if grid.gauge != "conformal":
        raise ValueError("conformal gauge is required for this preparation tangent")
    _check_state_dimensions(grid, state)
    if grid.nq % 8 or abs(float(grid.xi_q[grid.nq // 8]) - 1.0) > ZERO_GATE:
        raise ValueError("quadrature grid does not contain the observer at x=1")
    _fail("preparation premise failed", _geometry_problems(grid, state))
    _fail(
        "preparation premise failed",
        _pair_problems(state.phi0, state.phi1, grid.fine.occupations),
    )
    generator = preparation_generator(grid.xi_f, grid.length)
    phi0_tangent, phi1_tangent = spinor_tangent(state.phi0, state.phi1, generator)
    if phi0_tangent.shape != state.phi0.shape or phi1_tangent.shape != state.phi1.shape:
        raise ValueError("spinor tangent does not match the fermion columns")

    fine_state = prolong_state(grid, state)
    q0 = float(CALIBRATION["b0"] / CALIBRATION["a0"])
    chart_problems = []
    if float(np.max(np.abs(fine_state.Q - q0))) > ZERO_GATE:
        chart_problems.append("prolonged Q left the calibration ratio")
    if float(np.min(fine_state.r)) <= 0.0:
        chart_problems.append("prolonged radius is not positive")
    _fail("preparation premise failed", chart_problems)

    # Pre-gauge system: force_L and the Hamilton density do not use L or beta.
    # This is the operator differentiated by _projected_radius_jacobian.
    source = source_from_columns(grid.fine, fine_state)
    check()
    rho = source["force_L"] / grid.dx_q
    projected_root, full_root = projected_radius_operator(grid, state.r, rho)
    if not np.isfinite(full_root).all() or not np.isfinite(projected_root).all():
        raise ValueError("original constraint level is not finite")
    root_max = float(np.max(np.abs(projected_root)))
    momentum = shift_constraint(grid.fine, fine_state) + source["force_beta"] / grid.dx_q
    source_problems = []
    s1_max = _abs_max(source["S1"])
    pmom_max = _abs_max(source["Pmom"])
    current_max = _abs_max(momentum)
    if s1_max > PAIR_GATE:
        source_problems.append(f"S1 max {s1_max} exceeds {PAIR_GATE}")
    if pmom_max > PAIR_GATE:
        source_problems.append(f"Pmom max {pmom_max} exceeds {PAIR_GATE}")
    if current_max > TOL_INITIAL_CONSTRAINT:
        source_problems.append(f"momentum/current max {current_max} exceeds {TOL_INITIAL_CONSTRAINT}")
    _fail("background source premise failed", source_problems)
    check()

    fine_phi0 = prolong_columns(grid, phi0_tangent)
    fine_phi1 = prolong_columns(grid, phi1_tangent)
    zero_radius = np.zeros(grid.nq)
    original_h, _d, phi_h_dot, phi_d_dot = _polarization(
        grid.fine, fine_state, source, zero_radius, fine_phi0, fine_phi1,
    )
    original_density = _constraint_density(original_h, fine_state.Q)
    original_projected = pull_geometry(grid, original_density)
    original_held = original_density - prolong_geometry(grid, original_projected)
    phi_momentum = _abs_max(phi_d_dot)
    if phi_momentum > TOL_INITIAL_CONSTRAINT:
        raise ValueError(
            f"momentum variation is not zero for this paired family: max {phi_momentum}"
        )
    delta_density = _constraint_density(phi_h_dot, fine_state.Q)
    projected_source = pull_geometry(grid, delta_density)
    if projected_source.shape != (grid.ng,) or not np.isfinite(projected_source).all():
        raise ValueError("projected source tangent is not a finite coarse vector")

    # Shapes are fixed before the only dense factorization.
    check()
    jacobian = _projected_radius_jacobian(grid, state.r)
    check()
    if jacobian.shape != (grid.ng, grid.ng):
        raise ValueError("projected radius Jacobian is not ng by ng")
    delta_r, factor = _factor_once(jacobian, -projected_source)
    check()
    response = jacobian @ delta_r
    residual = response + projected_source
    residual_max = float(np.max(np.abs(residual)))
    source_scale = float(np.max(np.abs(projected_source)))
    geometric_scale = float(np.max(np.abs(response)))
    scale = max(source_scale, geometric_scale)
    relative = _relative(residual_max, scale)
    if relative > RELATIVE_SOLVE_GATE:
        raise ValueError(
            "projected solve residual "
            f"{residual_max} exceeds {RELATIVE_SOLVE_GATE} relative to response scale {scale}"
        )

    radius_tangent = prolong_geometry(grid, delta_r)
    _h, _d, full_h_dot, full_d_dot = _polarization(
        grid.fine, fine_state, source, radius_tangent, fine_phi0, fine_phi1,
    )
    check()
    full_density = _constraint_density(full_h_dot, fine_state.Q)
    projected_full = pull_geometry(grid, full_density)
    held = full_density - prolong_geometry(grid, projected_full)
    projected_max = float(np.max(np.abs(projected_full)))
    disagreement = float(np.max(np.abs(projected_full - residual)))
    full_max = float(np.max(np.abs(full_density)))
    held_max = float(np.max(np.abs(held)))
    full_momentum = _abs_max(full_d_dot)
    if _relative(projected_max, scale) > RELATIVE_SOLVE_GATE:
        raise ValueError(
            "directional projected cancellation "
            f"{projected_max} exceeds {RELATIVE_SOLVE_GATE} relative to response scale {scale}"
        )
    if _relative(disagreement, scale) > RELATIVE_SOLVE_GATE:
        raise ValueError(
            "constraint polarization and factored radius Jacobian disagree by "
            f"{disagreement} on response scale {scale}"
        )
    if full_momentum > TOL_INITIAL_CONSTRAINT:
        raise ValueError(f"full momentum variation is not zero: max {full_momentum}")

    observer = grid.nq // 8
    radius_q = prolong_geometry(grid, state.r)
    delta_q = prolong_geometry(grid, delta_r)
    r_obs = float(radius_q[observer])
    dr_obs = float(delta_q[observer])
    if r_obs <= 0.0:
        raise ValueError("observer radius is not positive")
    fraction = dr_obs / r_obs
    ratio_q = delta_q / radius_q
    ratio_g = delta_r / state.r

    def rms(values):
        return float(np.sqrt(np.mean(np.square(values))))

    gram = state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1
    report = {
        "generator": "cos(2*pi*xi_f/length)",
        "applied_on": "nf_before_prolongation",
        "recorded_phase": float(CALIBRATION["phase"]),
        "pairs": 3,
        "columns": 6,
        "ng": int(grid.ng),
        "nf": int(grid.nf),
        "nq": int(grid.nq),
        "constraint_target": "preserve_original_projected_constraint_level",
        "gram_defect_max_abs": float(np.max(np.abs(gram - np.eye(6)))),
        "background": {
            "projected_hamilton_max_abs": float(np.max(np.abs(original_projected))),
            "full_hamilton_max_abs": float(np.max(np.abs(original_density))),
            "held_out_hamilton_max_abs": float(np.max(np.abs(original_held))),
            "constant_Q_radius_operator_projected_max_abs": root_max,
            "difference_from_constant_Q_operator_max_abs": _abs_max(original_density-full_root),
            "momentum_max_abs": current_max,
            "S1_max_abs": s1_max,
            "Pmom_max_abs": pmom_max,
            "Q": q0,
            "original_newton_threshold": float(TOL_NEWTON),
            "original_newton_threshold_met": bool(root_max < TOL_NEWTON),
            "newton_threshold_role": "historical constant-Q radius operator diagnostic",
            "initial_defect_removed": False,
        },
        "source_tangent": {
            "phi_only_density_max_abs": float(np.max(np.abs(delta_density))),
            "phi_only_momentum_max_abs": phi_momentum,
            "projected_source_max_abs": source_scale,
        },
        "lu": factor,
        "solve": {
            "residual_max_abs": residual_max,
            "source_response_max_abs": source_scale,
            "geometric_response_max_abs": geometric_scale,
            "relative_residual": float(relative),
            "tolerance_relative": RELATIVE_SOLVE_GATE,
            "certified_invertibility": False,
        },
        "directional_check": {
            "projected_max_abs": projected_max,
            "full_max_abs": full_max,
            "held_out_max_abs": held_max,
            "held_out_rms": rms(held),
            "momentum_max_abs": full_momentum,
            "disagreement_with_factored_residual_max_abs": disagreement,
        },
        "radius_response_per_epsilon": {
            "label": "derivative per epsilon; not an observed length",
            "coarse_max_abs": float(np.max(np.abs(delta_r))),
            "coarse_rms": rms(delta_r),
            "quadrature_max_abs": float(np.max(np.abs(delta_q))),
            "quadrature_rms": rms(delta_q),
            "dimensionless_coarse_max_abs": float(np.max(np.abs(ratio_g))),
            "dimensionless_coarse_rms": rms(ratio_g),
            "dimensionless_quadrature_max_abs": float(np.max(np.abs(ratio_q))),
            "dimensionless_quadrature_rms": rms(ratio_q),
        },
        "observer": {
            "x": 1.0,
            "quadrature_index": int(observer),
            "sample": "prolonged initial observer radius",
            "original_radius": r_obs,
            "Q_fixed": q0,
            "delta_r_per_epsilon": dr_obs,
            "delta_r_over_r_per_epsilon": fraction,
            "delta_N_over_N_per_epsilon": fraction,
            "clock_identity": "Q fixed, so delta_N/N = delta_r/r",
            "label": (
                "derivative per epsilon at the original initial observer; "
                "not an observed length, Lambda, or a continuum or formation result"
            ),
            "light_cone_evolution": False,
        },
        "geometry_variations_per_epsilon": {
            "Q": 0.0, "chi": 0.0, "p_Q": 0.0, "p_r": 0.0, "p_chi": 0.0,
        },
        "fixed_assumptions": [
            "same U on every column, applied on nf before prolongation",
            "plus phi1 = i phi0 and minus = recorded phase times conjugate(plus)",
            "equal pair occupations 0.75, 0.5 and 0.25",
            "delta Q, delta chi and all geometric momentum variations are zero",
            "J_r is the existing pre-gauge projected radius Jacobian",
            "full quadrature remainder is reported and is not removed",
            "the saved nonzero projected constraint level is held fixed",
            "initial Q is held fixed, so the observer clock variation is delta_r/r",
        ],
        "incomplete_continuum_constraints": {
            "continuum_constraint_certified": False,
            "time_integral_certified": False,
            "full_quadrature_constraint_removed": False,
            "certified_invertibility": False,
            "retarded_kernel": False,
            "later_trajectory_keeps_this_chirality": False,
            "formation_history": False,
        },
    }
    return PreparationTangent(
        delta_r=np.array(delta_r, copy=True),
        phi0_tangent=np.array(phi0_tangent, copy=True),
        phi1_tangent=np.array(phi1_tangent, copy=True),
        report=report,
    )
