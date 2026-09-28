"""Branch-owned SGB-L constraint continuity and conditional propagation.

Pointwise identities (unmodified scalars, affine source, hat-wave principal
form) remain.  This module additionally re-evaluates, on SGB-L data with
``F=Mpl^2>0`` and the locked hat factor,

* the exact invertible map from unredefined ``(H,M)`` on the MHG shell to
  ``(nabla_t C^t, nabla_t C^r)``;
* both radial hat-cone characteristic branches; and
* the FO1 kinematic reduction subsidiary identity.

Those identities support one *conditional* boundary-free propagation
statement.  They do not construct a solution, initial slice, centre, or
IBVP, and they do not open ``SGBL_branch_owned_and_healthy``.  Neighboring
CON4/HYP1 aggregate certificates are not imported.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from hashlib import sha256
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .exact_linear_algebra import add as poly_add
from .exact_linear_algebra import matrix_det, matrix_inverse
from .exact_linear_algebra import multiply as poly_mul
from .exact_linear_algebra import poly
from .exact_linear_algebra import scale as poly_scale
from .modified_harmonic import (
    _physical_metric,
    _trace_reversal_projector,
    auxiliary_inverse_metric,
    inverse_metric_null_polynomial,
)
from .modified_harmonic_constraints import physical_constraint_projections
from .modified_harmonic_reduction_subsidiary import (
    ReductionDifferentialFieldJet,
    ReductionDifferentialState,
    reduction_subsidiary_identity,
)
from .modified_harmonic_reference import (
    modified_harmonic_extension_residual,
    modified_harmonic_full_residuals,
    modified_harmonic_gauge_constraint,
)
from .reference_connection import flat_spherical_annulus_reference
from .sgb1_ctl1_source import (
    HAT_NORMAL_FACTOR,
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceInputs,
    SOURCE_ACCELERATION_ORDER,
    TILDE_NORMAL_FACTOR,
    sgbl_source_coefficients,
    sgbl_source_residual,
    sgbl_source_state,
)
from .spherical_reduction import BASE_FIELD_ORDER, residuals
from .spherical_symbol import quadratic_companion_certificate


Q = Fraction
N = 4
ZERO_ACCELERATION = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
HOLDOUT_ACCELERATION = (Q(2), Q(-1), Q(1, 2), Q(3), Q(-2), Q(1))
REFERENCE = flat_spherical_annulus_reference(
    radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
)
FORBIDDEN_NEIGHBOR_IMPORTS = (
    "modified_harmonic_constraint_propagation_certificate",
    "constraint_system_certificate",
    "physical_to_normal_gauge_map",
    "gauge_constraint_characteristics",
    "weak_coupling_health_certificate",
    "background_from_spherical_state",
    "FGCQRActionParameters",
    "activated_compatible_state",
)
ACTIVE_SPHERICAL_GAUGE_COMPONENTS = ("C^t", "C^r")
CONDITIONAL_PROPAGATION_ASSUMPTIONS = (
    "sufficiently_smooth_compatible_solution_assumed",
    "C_H_M_and_radial_reduction_zero_on_initial_slice_assumed",
    "complete_MHG_metric_equations_assumed",
    "both_unmodified_scalar_equations_hold_along_the_solution_assumed",
    "linear_normally_hyperbolic_uniqueness_for_the_hat_wave_assumed",
    "FO1_kinematic_equations_D_and_K_hold_as_differentiable_identities_assumed",
    "boundary_free_domain_of_dependence",
)
UNCONDITIONAL_NONCLAIMS = (
    "smooth_solution_exists",
    "compatible_initial_data_supplied",
    "regular_center_supplied",
    "boundary_condition_supplied",
    "IBVP_proved",
    "unconditional_domain_propagation",
    "SGBL_branch_owned_and_healthy",
)
CONDITIONAL_LOGICAL_CHAIN = (
    "C_zero_on_initial_slice_implies_spatial_nabla_C_zero",
    "H_equals_M_equals_zero_and_MHG_metric_equations_imply_normal_nabla_C_zero_by_invertible_map",
    "unmodified_scalars_plus_MHG_give_homogeneous_hat_wave_subsidiary",
    "named_linear_normally_hyperbolic_uniqueness_gives_C_zero_in_the_boundary_free_domain",
    "vanishing_C_kills_the_extension_so_MHG_recovers_unredefined_equations_and_H_M_remain_zero",
    "FO1_kinematic_identity_preserves_the_six_radial_reduction_constraints",
)
NAMED_HAT_WAVE_UNIQUENESS_PREMISE = (
    "On any sufficiently smooth solution of the complete MHG metric equations "
    "and both unmodified scalar equations, if F=Mpl^2>0, the locked-hat metric "
    "is smooth Lorentzian, the initial slice is hat-spacelike, and both C and "
    "its hat-normal derivative vanish on that slice, then standard uniqueness "
    "for the homogeneous normally hyperbolic hat-wave system gives C=0 "
    "throughout the boundary-free domain of dependence. This is a named "
    "linear-PDE premise. It is not an imported CON3 certificate, not an "
    "energy estimate derived in this module, and not an existence theorem."
)
CONDITIONAL_STATEMENT = (
    "If a sufficiently smooth compatible SGB-L solution exists, C=H=M="
    "radial-reduction constraints=0 on a spacelike initial slice, the "
    "complete MHG metric equations and both unmodified scalar equations "
    "hold, then the exact invertible normal map, the homogeneous hat-wave "
    "subsidiary system, linear normally-hyperbolic uniqueness, and the "
    "FO1 kinematic identity preserve those constraints throughout the "
    "boundary-free domain of dependence of the slice."
)


class SGBLContinuityStop(ArithmeticError):
    """Typed continuity stop for a broken identity or chart/action error."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


def _six(name: str, value: object) -> tuple[Fraction, ...]:
    if not isinstance(value, (tuple, list)) or len(value) != 6:
        raise TypeError(f"{name} must be a 6-tuple")
    return tuple(Fraction(entry) for entry in value)


def _vector_sha256(value: Sequence[Fraction]) -> str:
    payload = ",".join(f"{entry.numerator}/{entry.denominator}" for entry in value)
    return sha256(payload.encode()).hexdigest()


def _mhg_full(point: SGBLSourceInputs, accelerations: Sequence[Fraction]) -> dict[str, Any]:
    state = sgbl_source_state(point, accelerations)
    return modified_harmonic_full_residuals(
        state,
        reference=REFERENCE,
        coordinate_radius=point.coordinate_radius,
        tilde_normal_factor=TILDE_NORMAL_FACTOR,
        hat_normal_factor=HAT_NORMAL_FACTOR,
    )


def sgbl_scalar_equations_unmodified(
    point: SGBLSourceInputs,
    accelerations: Sequence[Fraction | int] = ZERO_ACCELERATION,
) -> dict[str, Any]:
    """Check that MHG adds nothing to the SGB-L scalar rows."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    if point.alpha_gb == 0:
        raise ValueError("scalar-continuity check requires nonzero alpha_gb")
    full = _mhg_full(point, _six("accelerations", accelerations))
    unredefined = residuals(sgbl_source_state(point, accelerations))
    scalar_full = full["full_residual_vector"][4:]
    scalar_original = (unredefined["phi"], unredefined["chi"])
    unmodified = (
        full["scalar_equations_unmodified"] is True
        and scalar_full == scalar_original
        and full["extension_residual_vector"][4:] == (Q(0), Q(0))
    )
    if not unmodified:
        raise SGBLContinuityStop(
            "broken_gauge",
            "MHG extension modified an SGB-L scalar equation",
            {"full": scalar_full, "unredefined": scalar_original},
        )
    return {
        "scalar_equations_unmodified": True,
        "phi_residual": scalar_original[0],
        "chi_residual": scalar_original[1],
        "extension_scalar_rows": (Q(0), Q(0)),
        "F_prime": Q(0),
        "F": point.planck_mass ** 2,
    }


def sgbl_source_jacobian_identity(
    point: SGBLSourceInputs,
    *,
    holdout: Sequence[Fraction | int] = HOLDOUT_ACCELERATION,
) -> dict[str, Any]:
    """Re-read ``R(a)=R0+J a`` on SGB-L without borrowing a neighboring result."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    coefficients = sgbl_source_coefficients(point)
    holdout_a = _six("holdout", holdout)
    reconstructed = coefficients.evaluate(holdout_a)
    actual = sgbl_source_residual(point, holdout_a)
    if reconstructed != actual:
        raise SGBLContinuityStop(
            "source_identity_failed",
            "affine source reconstruction does not match the full residual",
            {"reconstructed": reconstructed, "actual": actual},
        )
    return {
        "affine_identity_holds": True,
        "holdout": holdout_a,
        "residual_sha256": _vector_sha256(actual),
        "jacobian_sha256": _vector_sha256(
            tuple(entry for row in coefficients.jacobian for entry in row)
        ),
        "acceleration_order": SOURCE_ACCELERATION_ORDER,
    }


def sgbl_physical_constraints_acceleration_independent(
    point: SGBLSourceInputs,
) -> dict[str, Any]:
    """Check that unredefined H,M do not depend on ADM accelerations."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    zero = physical_constraint_projections(sgbl_source_state(point, ZERO_ACCELERATION))
    holdout = physical_constraint_projections(
        sgbl_source_state(point, HOLDOUT_ACCELERATION)
    )
    if zero["H"] != holdout["H"] or zero["M"] != holdout["M"]:
        raise SGBLContinuityStop(
            "source_identity_failed",
            "physical H,M changed under a pure acceleration shift",
            {"zero": (zero["H"], zero["M"]), "holdout": (holdout["H"], holdout["M"])},
        )
    return {
        "physical_constraints_independent_of_accelerations": True,
        "H": zero["H"],
        "M": zero["M"],
        "source_is_unredefined_residual": zero["source_is_unredefined_residual"],
    }


def sgbl_mhg_hat_wave_identity(
    point: SGBLSourceInputs,
    *,
    mutate_trace_sign: int = -1,
) -> dict[str, Any]:
    """Independent hat-projector contraction at SGB-L ``F=Mpl^2``.

    On the scalar equations, Noether removes the unredefined metric
    divergence.  What remains is ``(F/2) hatq H``.  The contraction is
    performed here; a CON4 certificate is not imported.
    """

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    if mutate_trace_sign not in (-1, 1):
        raise ValueError("trace-sign mutation must be +1 or -1")
    state = sgbl_source_state(point, ZERO_ACCELERATION)
    if state.beta != 0 or state.eta != 0:
        raise SGBLContinuityStop(
            "action_identity_error",
            "SGB-L hat-wave identity requires beta=eta=0",
        )
    _, physical_inverse = _physical_metric(state)
    hat = auxiliary_inverse_metric(physical_inverse, HAT_NORMAL_FACTOR)
    effective_planck = point.planck_mass ** 2
    xi = (poly((0, -1)), poly((1,)), poly((0,)), poly((0,)))

    def _projector(alpha: int, beta: int, mu: int, nu: int) -> Fraction:
        if mutate_trace_sign == -1:
            return _trace_reversal_projector(hat, alpha, beta, mu, nu)
        return (
            Q(int(alpha == mu)) * hat[nu][beta]
            + Q(int(alpha == nu)) * hat[mu][beta]
            + Q(mutate_trace_sign) * Q(int(alpha == beta)) * hat[mu][nu]
        ) / 2

    derived = tuple(
        tuple(
            _poly_accumulate(
                poly_scale(
                    poly_mul(xi[mu], xi[beta]),
                    effective_planck * _projector(alpha, beta, mu, nu),
                )
                for mu in range(N)
                for beta in range(N)
            )
            for alpha in range(N)
        )
        for nu in range(N)
    )
    expected_factor = poly_scale(
        inverse_metric_null_polynomial(hat), effective_planck / 2
    )
    expected = tuple(
        tuple(expected_factor if nu == alpha else poly((0,)) for alpha in range(N))
        for nu in range(N)
    )
    holds = derived == expected
    if mutate_trace_sign != -1:
        if holds:
            raise SGBLContinuityStop(
                "broken_gauge",
                "trace-sign mutation did not break the hat-wave identity",
            )
        raise SGBLContinuityStop(
            "broken_gauge",
            "trace-sign mutation breaks the SGB-L hat-wave identity",
            {"mutate_trace_sign": mutate_trace_sign},
        )
    if not holds:
        raise SGBLContinuityStop(
            "broken_gauge",
            "SGB-L hat-wave projector contraction failed",
        )
    return {
        "hat_wave_identity_holds": holds,
        "F": effective_planck,
        "F_prime": Q(0),
        "mutate_trace_sign": mutate_trace_sign,
        "expected_hat_wave_factor": expected_factor,
        "uses_locked_auxiliary_hat": HAT_NORMAL_FACTOR,
        "independent_of_con4_certificate": True,
    }


def _poly_accumulate(values) -> tuple[Fraction, ...]:
    total = poly((0,))
    for value in values:
        total = poly_add(total, value)
    return total


def _project_hm(
    tensor: Sequence[Sequence[Fraction]], shift: Fraction
) -> tuple[Fraction, Fraction]:
    hamiltonian = (
        tensor[0][0] - 2 * shift * tensor[0][1] + shift**2 * tensor[1][1]
    )
    momentum = tensor[0][1] - shift * tensor[1][1]
    return hamiltonian, momentum


def _zero_matrix4() -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(Q(0) for _ in range(N)) for _ in range(N))


def sgbl_physical_to_normal_gauge_map(
    point: SGBLSourceInputs,
    *,
    mutate_analytic_sign: int = 1,
    mutate_hat_factor: int | Fraction | None = None,
    mutate_computed_scale: int | Fraction | None = None,
) -> dict[str, Any]:
    """Exact invertible ``(H,M) -> (nabla_t C^t, nabla_t C^r)`` map on the MHG shell.

    The two columns of the forward matrix are the unredefined projections of
    ``E=-X`` with unit ``nabla_t C`` inputs.  They are compared entrywise
    with the analytic ADM formula at ``F=Mpl^2`` and the locked hat factor.
    The inverse of that matrix is the requested map from unredefined
    ``(H,M)`` to ``(nabla_t C^t, nabla_t C^r)``.  CON4's aggregate
    certificate is not imported.
    """

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    if mutate_analytic_sign not in (-1, 1):
        raise ValueError("analytic-sign mutation must be +1 or -1")
    state = sgbl_source_state(point, ZERO_ACCELERATION)
    if state.branch != "SGB-L" or state.beta != 0 or state.eta != 0:
        raise SGBLContinuityStop(
            "action_identity_error",
            "SGB-L normal map requires the linear branch with beta=eta=0",
        )
    effective_planck = point.planck_mass ** 2
    if effective_planck <= 0:
        raise SGBLContinuityStop(
            "action_identity_error",
            "SGB-L normal map requires F=Mpl^2>0",
        )
    hat_factor = (
        HAT_NORMAL_FACTOR
        if mutate_hat_factor is None
        else Fraction(mutate_hat_factor)
    )
    if hat_factor <= 1:
        raise SGBLContinuityStop(
            "broken_hat",
            "hat normal factor must exceed one",
            {"hat_normal_factor": hat_factor},
        )
    h_tt = state.h_tt.value
    h_tr = state.h_tr.value
    h_rr = state.h_rr.value
    if h_rr <= 0:
        raise SGBLContinuityStop(
            "broken_map",
            "coordinate slice must have positive radial metric",
        )
    lapse_squared = -h_tt + h_tr**2 / h_rr
    if lapse_squared <= 0:
        raise SGBLContinuityStop(
            "broken_map",
            "coordinate slice must be spacelike",
        )
    shift = h_tr / h_rr
    gauge = modified_harmonic_gauge_constraint(
        state,
        reference=REFERENCE,
        coordinate_radius=point.coordinate_radius,
        tilde_normal_factor=TILDE_NORMAL_FACTOR,
    )
    columns: list[tuple[Fraction, Fraction]] = []
    for component in range(2):
        derivative = [list(row) for row in _zero_matrix4()]
        derivative[0][component] = Q(1)
        basis_gauge = replace(
            gauge,
            covariant_constraint_derivative=tuple(tuple(row) for row in derivative),
        )
        extension = modified_harmonic_extension_residual(
            state, gauge=basis_gauge, hat_normal_factor=hat_factor
        )["covariant"]
        physical_residual = tuple(
            tuple(-entry for entry in row) for row in extension
        )
        columns.append(_project_hm(physical_residual, shift))
    computed = tuple(
        tuple(columns[column][row] for column in range(2)) for row in range(2)
    )
    if mutate_computed_scale is not None:
        scale = Fraction(mutate_computed_scale)
        computed = tuple(tuple(scale * entry for entry in row) for row in computed)
    locked_analytic = (
        (effective_planck * HAT_NORMAL_FACTOR * lapse_squared / 2, Q(0)),
        (
            -effective_planck * HAT_NORMAL_FACTOR * h_tr / 2,
            -effective_planck * HAT_NORMAL_FACTOR * h_rr / 2,
        ),
    )
    analytic = (
        (mutate_analytic_sign * locked_analytic[0][0], Q(0)),
        (
            mutate_analytic_sign * locked_analytic[1][0],
            mutate_analytic_sign * locked_analytic[1][1],
        ),
    )
    if computed != analytic:
        if mutate_hat_factor is not None:
            reason = "broken_hat"
        elif mutate_analytic_sign != 1:
            reason = "broken_sign"
        else:
            reason = "broken_map"
        raise SGBLContinuityStop(
            reason,
            "direct MHG extension and locked-hat analytic physical-to-gauge map differ",
            {"computed": computed, "analytic": analytic},
        )
    determinant = matrix_det(computed)
    analytic_determinant = (
        -(effective_planck**2)
        * HAT_NORMAL_FACTOR**2
        * lapse_squared
        * h_rr
        / 4
    )
    if determinant != analytic_determinant:
        raise SGBLContinuityStop(
            "broken_map",
            "normal-map determinant disagrees with the locked-hat analytic formula",
            {"determinant": determinant, "analytic_determinant": analytic_determinant},
        )
    if determinant >= 0:
        raise SGBLContinuityStop(
            "broken_sign",
            "normal-map determinant is not strictly negative",
            {"determinant": determinant},
        )
    inverse_map = matrix_inverse(computed)
    return {
        "forward_input_order": ("nabla_t_C^t", "nabla_t_C^r"),
        "forward_output_order": (
            "unredefined_H_on_MHG_shell",
            "unredefined_M_on_MHG_shell",
        ),
        "input_order": ("nabla_t_C^t", "nabla_t_C^r"),
        "output_order": ("unredefined_H_on_MHG_shell", "unredefined_M_on_MHG_shell"),
        "computed_matrix": computed,
        "analytic_matrix": analytic,
        "normal_to_physical_matrix": computed,
        "physical_to_normal_derivative_matrix": inverse_map,
        "map_from_physical_HM_to_normal_nabla_C": inverse_map,
        "direct_and_analytic_matrices_equal": True,
        "determinant": determinant,
        "analytic_determinant": analytic_determinant,
        "analytic_determinant_formula": "-F^2 q^2 ell^2 h_rr / 4",
        "determinant_strictly_negative": True,
        "inverse_matrix": inverse_map,
        "F": effective_planck,
        "F_prime": Q(0),
        "hat_normal_factor": HAT_NORMAL_FACTOR,
        "lapse_squared": lapse_squared,
        "h_rr": h_rr,
        "h_tr": h_tr,
        "shift": shift,
        "C_zero_on_slice_required": True,
        "spatial_covariant_derivative_of_C_zero_on_slice_follows": True,
        "H_and_M_zero_iff_normal_nabla_C_zero_on_MHG_shell": True,
        "independent_of_con4_certificate": True,
    }


def sgbl_hat_cone_characteristics(
    point: SGBLSourceInputs,
    *,
    mutate_hat_factor: int | Fraction | None = None,
    mutate_root_signs: bool = False,
) -> dict[str, Any]:
    """Both radial hat-cone branches of the homogeneous gauge subsidiary system."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    state = sgbl_source_state(point, ZERO_ACCELERATION)
    if state.branch != "SGB-L" or state.beta != 0 or state.eta != 0:
        raise SGBLContinuityStop(
            "action_identity_error",
            "SGB-L hat characteristics require the linear branch",
        )
    hat_factor = (
        HAT_NORMAL_FACTOR
        if mutate_hat_factor is None
        else Fraction(mutate_hat_factor)
    )
    if hat_factor <= 1:
        raise SGBLContinuityStop(
            "broken_hat",
            "hat normal factor must exceed one",
            {"hat_normal_factor": hat_factor},
        )
    _, physical_inverse = _physical_metric(state)
    locked_hat = auxiliary_inverse_metric(physical_inverse, HAT_NORMAL_FACTOR)
    locked_polynomial = inverse_metric_null_polynomial(locked_hat)
    hat = auxiliary_inverse_metric(physical_inverse, hat_factor)
    polynomial = inverse_metric_null_polynomial(hat)
    if polynomial != locked_polynomial:
        raise SGBLContinuityStop(
            "broken_hat",
            "hat-null polynomial is not the locked-hat radial cone",
            {"polynomial": polynomial, "locked_polynomial": locked_polynomial},
        )
    companion = quadratic_companion_certificate(polynomial)
    roots = companion["roots"]

    def _sign(root: Mapping[str, Fraction]) -> int:
        if "exact" in root:
            value = root["exact"]
            return -1 if value < 0 else 1 if value > 0 else 0
        lower, upper = root["lower"], root["upper"]
        if upper < 0:
            return -1
        if lower > 0:
            return 1
        return 0

    signs = tuple(_sign(root) for root in roots)
    if mutate_root_signs:
        signs = tuple(1 for _ in signs)
    if any(value == 0 for value in signs) or set(signs) != {-1, 1}:
        raise SGBLContinuityStop(
            "broken_root",
            "hat-cone characteristics do not supply one incoming and one outgoing radial root",
            {"root_signs": signs},
        )
    leading = polynomial[-1]
    constant = polynomial[0]
    product_sign = -1 if (constant / leading) < 0 else 1 if (constant / leading) > 0 else 0
    if product_sign >= 0:
        raise SGBLContinuityStop(
            "broken_root",
            "hat-null polynomial does not have opposite-sign radial roots",
        )
    exact_roots = tuple(root["exact"] for root in roots if "exact" in root)
    return {
        "principal_equation": "(F/2)*hat_g^ab*nabla_a*nabla_b*C^mu plus lower order equals zero",
        "radial_covector_convention": "xi_a=(-c,1)",
        "hat_null_polynomial": polynomial,
        "root_isolations": roots,
        "exact_roots": exact_roots if len(exact_roots) == len(roots) else None,
        "root_signs": signs,
        "has_both_radial_branches": True,
        "discriminant_positive": companion["discriminant"] > 0,
        "discriminant": companion["discriminant"],
        "hat_normal_factor": HAT_NORMAL_FACTOR,
        "F": point.planck_mass ** 2,
        "active_spherical_component_order": ACTIVE_SPHERICAL_GAUGE_COMPONENTS,
        "independent_of_con4_certificate": True,
    }


def _off_shell_reduction_state() -> ReductionDifferentialState:
    """Six-field off-shell jet with commuting-partial identity still exact."""

    probe = ReductionDifferentialFieldJet(
        u_t=1, p=2, q=3, u_r=4, q_t=5, p_r=6, u_tr=7
    )
    return ReductionDifferentialState(
        **{field: probe for field in BASE_FIELD_ORDER}
    )


def _source_point_reduction_state(point: SGBLSourceInputs) -> ReductionDifferentialState:
    state = sgbl_source_state(point, ZERO_ACCELERATION)
    fields = {}
    for name in BASE_FIELD_ORDER:
        jet = getattr(state, name)
        fields[name] = ReductionDifferentialFieldJet(
            u_t=jet.dt,
            p=jet.dt,
            q=jet.dr,
            u_r=jet.dr,
            q_t=jet.dtr,
            p_r=jet.dtr,
            u_tr=jet.dtr,
        )
    return ReductionDifferentialState(**fields)


def sgbl_fo1_reduction_identity(
    point: SGBLSourceInputs,
    *,
    mutate_identity: bool = False,
) -> dict[str, Any]:
    """Bind the FO1 kinematic identity independently on SGB-L jets.

    The off-shell probe has nonzero ``D``, ``C`` and ``K`` but still satisfies
    ``d_t C + d_r D - K = 0`` field-by-field.  The source-point two-jet is the
    commuting embedding of the SGB-L spherical jet.
    """

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    off_shell = reduction_subsidiary_identity(_off_shell_reduction_state())
    on_point = reduction_subsidiary_identity(_source_point_reduction_state(point))
    residuals = off_shell["off_shell_identity_residual_d_t_C_plus_d_r_D_minus_K"]
    if mutate_identity:
        residuals = (residuals[0] + 1,) + residuals[1:]
    field_by_field_zero = all(value == 0 for value in residuals)
    if not field_by_field_zero:
        raise SGBLContinuityStop(
            "broken_assumption",
            "FO1 reduction subsidiary identity failed field-by-field",
            {"residuals": residuals},
        )
    if not off_shell["off_shell_identity_exact"] or not on_point["off_shell_identity_exact"]:
        raise SGBLContinuityStop(
            "broken_assumption",
            "FO1 reduction identity is not exact on the SGB-L embedding",
        )
    return {
        "off_shell_identity_exact": True,
        "source_point_identity_exact": True,
        "field_by_field_zero": True,
        "field_order": BASE_FIELD_ORDER,
        "off_shell_identity_residual": off_shell[
            "off_shell_identity_residual_d_t_C_plus_d_r_D_minus_K"
        ],
        "off_shell_D": off_shell["u_time_definition_residual_D"],
        "off_shell_C": off_shell["radial_reduction_constraint_C"],
        "off_shell_K": off_shell["mixed_partial_compatibility_residual_K"],
        "off_shell_constraints_need_not_vanish": True,
        "source_point_is_commuting_two_jet_embedding": True,
        "subsidiary_characteristic_speeds": off_shell["subsidiary_characteristic_speeds"],
        "incoming_reduction_constraint_fields": (),
        "independent_of_con4_certificate": True,
    }


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLContinuityRecord:
    """Algebraic continuity/compatibility record. Domain propagation stays closed."""

    point: SGBLSourceInputs
    branch: str
    action: Mapping[str, Fraction]
    scalar_equations_unmodified: bool
    source_affine_identity: bool
    physical_constraints_acceleration_independent: bool
    hat_wave_identity: bool
    mhg_extension_metric_only: bool
    constraint_continuity_compatible: bool
    classification: str
    scalar_payload: Mapping[str, Any]
    source_payload: Mapping[str, Any]
    physical_payload: Mapping[str, Any]
    hat_wave_payload: Mapping[str, Any]

    def __post_init__(self) -> None:
        if type(self.point) is not SGBLSourceInputs:
            raise TypeError("point must be SGBLSourceInputs")
        if self.branch != "SGB-L":
            raise ValueError("continuity record is owned by the linear branch")
        if self.classification not in {
            "constraint_continuity_compatible",
            "continuity_identity_failed",
        }:
            raise ValueError("unknown continuity classification")
        expected = (
            self.scalar_equations_unmodified
            and self.source_affine_identity
            and self.physical_constraints_acceleration_independent
            and self.hat_wave_identity
            and self.mhg_extension_metric_only
        )
        if self.constraint_continuity_compatible != expected:
            raise ValueError("compatibility flag does not match the stored identities")
        if self.classification == "constraint_continuity_compatible" and not expected:
            raise ValueError("compatible classification requires all identities")
        if self.classification == "continuity_identity_failed" and expected:
            raise ValueError("failed classification cannot retain all identities")
        object.__setattr__(self, "action", MappingProxyType(dict(self.action)))
        object.__setattr__(self, "scalar_payload", MappingProxyType(dict(self.scalar_payload)))
        object.__setattr__(self, "source_payload", MappingProxyType(dict(self.source_payload)))
        object.__setattr__(self, "physical_payload", MappingProxyType(dict(self.physical_payload)))
        object.__setattr__(self, "hat_wave_payload", MappingProxyType(dict(self.hat_wave_payload)))

    @property
    def constraint_propagation_on_a_domain_proven(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def strongly_hyperbolic(self) -> bool:
        return False

    @property
    def cone_certificate(self) -> str:
        return "unqualified"


def sgbl_constraint_continuity(point: SGBLSourceInputs) -> SGBLContinuityRecord:
    """Assemble the branch-owned continuity identities on one SGB-L point."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    if point.alpha_gb == 0:
        raise ValueError("production continuity requires nonzero alpha_gb")
    state = sgbl_source_state(point, ZERO_ACCELERATION)
    if state.branch != "SGB-L" or state.beta != 0 or state.eta != 0:
        raise SGBLContinuityStop(
            "action_identity_error",
            "continuity requires the linear SGB-L action",
        )
    scalar = sgbl_scalar_equations_unmodified(point)
    source = sgbl_source_jacobian_identity(point)
    physical = sgbl_physical_constraints_acceleration_independent(point)
    hat = sgbl_mhg_hat_wave_identity(point)
    compatible = (
        scalar["scalar_equations_unmodified"]
        and source["affine_identity_holds"]
        and physical["physical_constraints_independent_of_accelerations"]
        and hat["hat_wave_identity_holds"]
    )
    return SGBLContinuityRecord(
        point=point,
        branch="SGB-L",
        action={
            "planck_mass": point.planck_mass,
            "alpha_gb": point.alpha_gb,
            "scalar_mass": point.scalar_mass,
            "quartic_coupling": point.quartic_coupling,
            "beta": Q(0),
            "eta": Q(0),
        },
        scalar_equations_unmodified=scalar["scalar_equations_unmodified"],
        source_affine_identity=source["affine_identity_holds"],
        physical_constraints_acceleration_independent=physical[
            "physical_constraints_independent_of_accelerations"
        ],
        hat_wave_identity=hat["hat_wave_identity_holds"],
        mhg_extension_metric_only=True,
        constraint_continuity_compatible=compatible,
        classification=(
            "constraint_continuity_compatible"
            if compatible
            else "continuity_identity_failed"
        ),
        scalar_payload=scalar,
        source_payload=source,
        physical_payload=physical,
        hat_wave_payload=hat,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLConstraintPropagationContract:
    """Conditional boundary-free SGB-L constraint-propagation contract.

    The stored boolean is the implication, not existence of a solution.
    """

    point: SGBLSourceInputs
    branch: str
    action: Mapping[str, Fraction]
    compatibility: SGBLContinuityRecord
    physical_to_normal_map: Mapping[str, Any]
    hat_cone: Mapping[str, Any]
    fo1_reduction: Mapping[str, Any]
    map_determinant: Fraction
    map_determinant_strictly_negative: bool
    physical_to_normal_derivative_matrix: tuple[tuple[Fraction, ...], ...]
    hat_cone_root_signs: tuple[int, ...]
    assumptions: Mapping[str, bool]
    logical_chain: tuple[str, ...]
    named_uniqueness_premise: str
    conditional_statement: str
    conditional_boundary_free_constraint_propagation_proved: bool
    classification: str

    def __post_init__(self) -> None:
        if type(self.point) is not SGBLSourceInputs:
            raise TypeError("point must be SGBLSourceInputs")
        if type(self.compatibility) is not SGBLContinuityRecord:
            raise TypeError("compatibility must be SGBLContinuityRecord")
        if self.branch != "SGB-L":
            raise ValueError("propagation contract is owned by the linear branch")
        if self.classification not in {
            "conditional_boundary_free_constraint_propagation_proved",
            "continuity_identity_failed",
        }:
            raise ValueError("unknown propagation classification")
        if self.logical_chain != CONDITIONAL_LOGICAL_CHAIN:
            raise ValueError("logical chain is frozen for the SGB-L contract")
        if self.named_uniqueness_premise != NAMED_HAT_WAVE_UNIQUENESS_PREMISE:
            raise ValueError("uniqueness premise must remain the named non-imported statement")
        if self.conditional_statement != CONDITIONAL_STATEMENT:
            raise ValueError("conditional statement is frozen for the SGB-L contract")
        identities = (
            self.compatibility.constraint_continuity_compatible
            and self.map_determinant_strictly_negative
            and self.map_determinant < 0
            and self.physical_to_normal_map.get("direct_and_analytic_matrices_equal")
            is True
            and self.physical_to_normal_map.get("hat_normal_factor") == HAT_NORMAL_FACTOR
            and self.hat_cone.get("has_both_radial_branches") is True
            and self.hat_cone.get("hat_normal_factor") == HAT_NORMAL_FACTOR
            and set(self.hat_cone_root_signs) == {-1, 1}
            and 0 not in self.hat_cone_root_signs
            and self.fo1_reduction.get("off_shell_identity_exact") is True
            and self.fo1_reduction.get("field_by_field_zero") is True
            and all(
                self.assumptions.get(name) is True
                for name in CONDITIONAL_PROPAGATION_ASSUMPTIONS
            )
        )
        if self.conditional_boundary_free_constraint_propagation_proved != identities:
            raise ValueError(
                "conditional proof flag does not match the re-evaluated identities"
            )
        if (
            self.classification == "conditional_boundary_free_constraint_propagation_proved"
            and not identities
        ):
            raise ValueError("proved classification requires every exact identity")
        if (
            self.classification == "continuity_identity_failed"
            and self.conditional_boundary_free_constraint_propagation_proved
        ):
            raise ValueError("failed classification cannot retain a proved flag")
        object.__setattr__(self, "action", MappingProxyType(dict(self.action)))
        object.__setattr__(
            self, "physical_to_normal_map", MappingProxyType(dict(self.physical_to_normal_map))
        )
        object.__setattr__(self, "hat_cone", MappingProxyType(dict(self.hat_cone)))
        object.__setattr__(self, "fo1_reduction", MappingProxyType(dict(self.fo1_reduction)))
        object.__setattr__(self, "assumptions", MappingProxyType(dict(self.assumptions)))

    @property
    def constraint_propagation_on_a_domain_proven(self) -> bool:
        return False

    @property
    def smooth_solution_exists(self) -> bool:
        return False

    @property
    def compatible_initial_data_supplied(self) -> bool:
        return False

    @property
    def regular_center_supplied(self) -> bool:
        return False

    @property
    def boundary_condition_supplied(self) -> bool:
        return False

    @property
    def IBVP_proved(self) -> bool:
        return False

    @property
    def unconditional_domain_propagation(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
        return False

    @property
    def holdout(self) -> bool:
        return False

    @property
    def remaining_nonclaims(self) -> Mapping[str, bool]:
        return MappingProxyType(
            {
                name: False
                for name in (
                    *UNCONDITIONAL_NONCLAIMS,
                    "FRZ1",
                    "PREF1",
                    "holdout",
                    "constraint_preserving_boundary_map_supplied",
                    "constraint_propagation_on_a_domain_proven",
                )
            }
        )


def sgbl_constraint_propagation_contract(
    point: SGBLSourceInputs,
) -> SGBLConstraintPropagationContract:
    """Assemble the strongest honest conditional SGB-L propagation theorem."""

    compatibility = sgbl_constraint_continuity(point)
    normal_map = sgbl_physical_to_normal_gauge_map(point)
    hat_cone = sgbl_hat_cone_characteristics(point)
    fo1 = sgbl_fo1_reduction_identity(point)
    assumptions = {name: True for name in CONDITIONAL_PROPAGATION_ASSUMPTIONS}
    proved = (
        compatibility.constraint_continuity_compatible
        and normal_map["determinant_strictly_negative"]
        and normal_map["direct_and_analytic_matrices_equal"]
        and normal_map["hat_normal_factor"] == HAT_NORMAL_FACTOR
        and hat_cone["has_both_radial_branches"]
        and hat_cone["hat_normal_factor"] == HAT_NORMAL_FACTOR
        and set(hat_cone["root_signs"]) == {-1, 1}
        and fo1["off_shell_identity_exact"]
        and fo1["field_by_field_zero"]
        and all(assumptions[name] is True for name in CONDITIONAL_PROPAGATION_ASSUMPTIONS)
    )
    return SGBLConstraintPropagationContract(
        point=point,
        branch="SGB-L",
        action=dict(compatibility.action),
        compatibility=compatibility,
        physical_to_normal_map=normal_map,
        hat_cone=hat_cone,
        fo1_reduction=fo1,
        map_determinant=normal_map["determinant"],
        map_determinant_strictly_negative=normal_map["determinant_strictly_negative"],
        physical_to_normal_derivative_matrix=normal_map[
            "physical_to_normal_derivative_matrix"
        ],
        hat_cone_root_signs=hat_cone["root_signs"],
        assumptions=assumptions,
        logical_chain=CONDITIONAL_LOGICAL_CHAIN,
        named_uniqueness_premise=NAMED_HAT_WAVE_UNIQUENESS_PREMISE,
        conditional_statement=CONDITIONAL_STATEMENT,
        conditional_boundary_free_constraint_propagation_proved=proved,
        classification=(
            "conditional_boundary_free_constraint_propagation_proved"
            if proved
            else "continuity_identity_failed"
        ),
    )


__all__ = [
    "CONDITIONAL_LOGICAL_CHAIN",
    "CONDITIONAL_PROPAGATION_ASSUMPTIONS",
    "CONDITIONAL_STATEMENT",
    "FORBIDDEN_NEIGHBOR_IMPORTS",
    "HOLDOUT_ACCELERATION",
    "NAMED_HAT_WAVE_UNIQUENESS_PREMISE",
    "SGBLConstraintPropagationContract",
    "SGBLContinuityRecord",
    "SGBLContinuityStop",
    "UNCONDITIONAL_NONCLAIMS",
    "sgbl_constraint_continuity",
    "sgbl_constraint_propagation_contract",
    "sgbl_fo1_reduction_identity",
    "sgbl_hat_cone_characteristics",
    "sgbl_mhg_hat_wave_identity",
    "sgbl_physical_constraints_acceleration_independent",
    "sgbl_physical_to_normal_gauge_map",
    "sgbl_scalar_equations_unmodified",
    "sgbl_source_jacobian_identity",
]
