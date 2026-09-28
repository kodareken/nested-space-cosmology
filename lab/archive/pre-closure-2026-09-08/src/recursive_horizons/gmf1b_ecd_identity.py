"""GMF-1B-ECD-INT1: action/gamma interaction identity gate.

This bounded artifact defines the ECD action in a canonical ``(+,-,-,-)``
layer, translates it completely to GMF-1B-ECD-SYM1's ``(-,+,+,+)``
Ventrella--Choptuik gamma convention, and then checks only the sphere-reduced
axial contact interaction, its Wirtinger sources, and its direct Hehl--Datta
gamma-matrix projection onto the retained chiral doublet.

    It is not a time-evolution normalization, global spherical-harmonic
    projection theorem, full ECD stress calculation, constraint solution,
    curvature calculation, bounce, or global spacetime.
"""

from __future__ import annotations

from math import cos, isfinite, sin, sqrt
from numbers import Complex, Real
from typing import Final, Sequence

from .gmf1b_ecd_symmetry import (
    geometric_normalization,
    minimal_ecd_contact_scalar,
    ventrella_chiral_axial_bilinears,
)
from .gmf1b_ecd_symmetry import _matrix_product, _matrix_scale, _matrix_vector, _vc_gamma_matrices


_HD_COEFFICIENT: Final[float] = 3.0 / 8.0
_CONTACT_COEFFICIENT: Final[float] = 3.0 / 16.0


def _real(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _positive(name: str, value: Real) -> float:
    result = _real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _complex(name: str, value: Complex) -> complex:
    if isinstance(value, bool) or not isinstance(value, Complex):
        raise ValueError(f"{name} must be a finite complex number")
    result = complex(value)
    if not isfinite(result.real) or not isfinite(result.imag):
        raise ValueError(f"{name} must be a finite complex number")
    return result


def _finite(name: str, value: float) -> float:
    if not isfinite(value):
        raise ValueError(f"{name} is outside the finite range")
    return value


def vc_reduced_contact_density(
    F: Complex, G: Complex, lapse: Real, radial_metric: Real, areal_radius: Real, kappa: Real
) -> float:
    """Return the sphere-reduced contact density ``L4_red``.

    For the VC pair, ``L4_red=-3*kappa*lapse*|F|^2*|G|^2/(4*pi*r^2*a)``.
    It is an action density in ``dt dr``, not a Hamiltonian or energy density.
    """

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    alpha = _positive("lapse", lapse)
    a_value = _positive("radial_metric", radial_metric)
    radius = _positive("areal_radius", areal_radius)
    coupling = _positive("kappa", kappa)
    result = -3.0 * coupling * alpha * abs(f_value) ** 2 * abs(g_value) ** 2 / (4.0 * 3.141592653589793 * radius * radius * a_value)
    return _finite("reduced contact density", result)


def vc_contact_action_sources(
    F: Complex, G: Complex, lapse: Real, radial_metric: Real, areal_radius: Real, kappa: Real
) -> dict[str, complex | float]:
    """Return Wirtinger sources ``dL4/dF*`` and ``dL4/dG*``."""

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    alpha = _positive("lapse", lapse)
    a_value = _positive("radial_metric", radial_metric)
    radius = _positive("areal_radius", areal_radius)
    coupling = _positive("kappa", kappa)
    common = -3.0 * coupling * alpha / (4.0 * 3.141592653589793 * radius * radius * a_value)
    source_f = common * abs(g_value) ** 2 * f_value
    source_g = common * abs(f_value) ** 2 * g_value
    for name, value in (("source_F", source_f), ("source_G", source_g)):
        if not isfinite(value.real) or not isfinite(value.imag):
            raise ValueError(f"{name} is outside the finite range")
    return {"source_F": source_f, "source_G": source_g, "common": _finite("source common factor", common)}


def vc_signature_bridge(
    F: Complex,
    G: Complex,
    radial_metric: Real,
    areal_radius: Real,
    kappa: Real,
    theta: Real,
    phi: Real,
) -> dict[str, object]:
    """Certify the Cabral--Lobo--Rubiera-Garcia signature translation.

    The canonical action layer uses ``(+---)`` matrices ``Gamma`` with the
    action and field-equation signs displayed by Cabral et al.  VC's imaginary
    ``(-+++)`` matrices are related by ``Gamma^a=-i*gamma_VC^a``.
    With ``bar(psi)=psi^dagger*(-i*gamma_VC^0)``, the physical real axial
    current is ``J^a=A^a=-i*bar(psi)*gamma_VC^a*gamma_star*psi``.  This
    routine reports both that raw imaginary bilinear and the real current,
    together with the complete interaction/action and equation translation.

    It changes only an overall signature/matrix convention; it does not
    change a line element or construct a spacetime metric.
    """

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    a_value = _positive("radial_metric", radial_metric)
    radius = _positive("areal_radius", areal_radius)
    coupling = _positive("kappa", kappa)
    angle_theta = _real("theta", theta)
    angle_phi = _real("phi", phi)
    C = geometric_normalization(radius, a_value)
    gamma = _vc_gamma_matrices()
    gamma_star = _matrix_scale(1j, _matrix_product(gamma[0], gamma[1], gamma[2], gamma[3]))
    canonical_gamma = tuple(_matrix_scale(-1j, matrix) for matrix in gamma)
    canonical_gamma5 = _matrix_scale(
        1j,
        _matrix_product(
            canonical_gamma[0], canonical_gamma[1], canonical_gamma[2], canonical_gamma[3]
        ),
    )
    vc_clifford_residual = _clifford_residual(gamma, (-1.0, 1.0, 1.0, 1.0))
    canonical_clifford_residual = _clifford_residual(
        canonical_gamma, (1.0, -1.0, -1.0, -1.0)
    )
    gamma5_residual = _matrix_difference_residual(canonical_gamma5, gamma_star)
    adjoint_residual = _matrix_difference_residual(
        canonical_gamma[0], _matrix_scale(-1j, gamma[0])
    )
    plus, minus = _vc_pair_spinors(f_value, g_value, C, angle_theta, angle_phi)
    raw = tuple(
        _raw_vc_axial_bilinear(plus, gamma[index], gamma_star)
        + _raw_vc_axial_bilinear(minus, gamma[index], gamma_star)
        for index in range(4)
    )
    real_current = ventrella_chiral_axial_bilinears(f_value, g_value, radius, a_value)
    translated = (real_current["A_hat_0"], 0.0, 0.0, real_current["A_hat_r"])
    raw_to_real_residual = max(abs(-1j * raw[index] - translated[index]) for index in range(4))
    raw_minus_norm = -raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2] + raw[3] * raw[3]
    A_squared_minus = real_current["axial_squared"]
    J_squared_plus = -A_squared_minus
    norm_residual = max(abs(raw_minus_norm - J_squared_plus), abs(J_squared_plus + A_squared_minus))
    # Cabral Eq. 83 is brought to the canonical left-hand form
    # ``i Gamma.D psi - c J^a Gamma_a gamma5 psi = 0``.  Substituting
    # Gamma=-i gamma_VC gives exactly the VC left-hand term below.
    canonical_covariant_current_gamma = tuple(
        tuple(
            translated[0] * canonical_gamma[0][row][column]
            - translated[3] * canonical_gamma[3][row][column]
            for column in range(4)
        )
        for row in range(4)
    )
    canonical_lhs_contact = _matrix_scale(
        -_HD_COEFFICIENT * coupling,
        _matrix_product(canonical_covariant_current_gamma, canonical_gamma5),
    )
    vc_lhs_contact = _matrix_scale(
        -1j * _HD_COEFFICIENT * coupling,
        _matrix_product(_covariant_axial_gamma_operator(gamma, translated[0], translated[3]), gamma_star),
    )
    equation_residual = max(
        abs(canonical_lhs_contact[row][column] - vc_lhs_contact[row][column])
        for row in range(4)
        for column in range(4)
    )
    contact_plus = -_CONTACT_COEFFICIENT * coupling * J_squared_plus
    contact_minus = _CONTACT_COEFFICIENT * coupling * A_squared_minus
    contact_residual = abs(contact_plus - contact_minus)
    scale = max(1.0, abs(A_squared_minus), abs(contact_plus))
    maximum_bridge_residual = max(
        vc_clifford_residual,
        canonical_clifford_residual,
        gamma5_residual,
        adjoint_residual,
        raw_to_real_residual,
        norm_residual,
        contact_residual,
        equation_residual,
    )
    return {
        "source": (
            "Cabral-Lobo-Rubiera-Garcia 2019 arXiv:1902.02222 "
            "Eqs. 8-11, 18, 20, 23, 83"
        ),
        "canonical_signature": "+---",
        "vc_signature": "-+++",
        "gamma_conversion": "Gamma^a=-i*gamma_VC^a",
        "canonical_action_defined_before_translation": True,
        "vc_clifford_residual": _finite("VC Clifford residual", vc_clifford_residual),
        "canonical_clifford_residual": _finite(
            "canonical Clifford residual", canonical_clifford_residual
        ),
        "gamma5_translation_residual": _finite("gamma5 translation residual", gamma5_residual),
        "adjoint_translation_residual": _finite("adjoint translation residual", adjoint_residual),
        "real_current_definition": "J^a=A^a=-i*bar(psi)*gamma_VC^a*gamma_star*psi",
        "raw_vc_axial_bilinear": raw,
        "translated_real_current": translated,
        "raw_to_real_current_residual": _finite("raw to real current residual", raw_to_real_residual),
        "raw_bilinear_minus_signature_quadratic_form_not_a_physical_norm": raw_minus_norm,
        "J_squared_plus": _finite("J squared plus", J_squared_plus),
        "A_squared_minus": _finite("A squared minus", A_squared_minus),
        "norm_translation_residual": _finite("norm translation residual", norm_residual),
        "contact_lagrangian_plus": _finite("contact plus", contact_plus),
        "contact_lagrangian_minus": _finite("contact minus", contact_minus),
        "contact_translation_residual": _finite("contact translation residual", contact_residual),
        "canonical_lhs_to_vc_equation_residual": _finite("equation translation residual", equation_residual),
        "full_effective_contact_coefficient_plus": "-3*kappa/16",
        "translated_contact_coefficient_minus": "+3*kappa/16",
        "hehl_datta_coefficient": "+3*kappa/8",
        "maximum_signature_bridge_residual": _finite(
            "maximum signature bridge residual", maximum_bridge_residual
        ),
        "signature_bridge_closed": maximum_bridge_residual <= 1.0e-12 * scale,
        "line_element_physically_changed": False,
    }


def vc_hehl_datta_projection(
    F: Complex,
    G: Complex,
    radial_metric: Real,
    areal_radius: Real,
    kappa: Real,
    theta: Real,
    phi: Real,
) -> dict[str, complex | float | bool]:
    """Project the direct minimal-ECD cubic term onto both VC retained modes.

    The returned coefficients multiply the *opposite-chirality output* basis.
    They are covariant Dirac-equation terms, not time derivatives of ``F,G``.
    """

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    a_value = _positive("radial_metric", radial_metric)
    radius = _positive("areal_radius", areal_radius)
    coupling = _positive("kappa", kappa)
    angle_theta = _real("theta", theta)
    angle_phi = _real("phi", phi)
    sin_half = sin(angle_theta / 2.0)
    cos_half = cos(angle_theta / 2.0)
    if abs(sin_half) < 1.0e-12 or abs(cos_half) < 1.0e-12:
        raise ValueError("theta must be nondegenerate for retained-mode projection")
    C = geometric_normalization(radius, a_value)
    bridge = vc_signature_bridge(f_value, g_value, a_value, radius, coupling, angle_theta, angle_phi)
    axial = ventrella_chiral_axial_bilinears(f_value, g_value, radius, a_value)
    gamma = _vc_gamma_matrices()
    gamma_star = _matrix_scale(1j, _matrix_product(gamma[0], gamma[1], gamma[2], gamma[3]))
    # VC's imaginary matrices obey ``Gamma^a=-i*gamma_VC^a``, where
    # ``Gamma`` is the standard Clifford field entering the ECD equation.
    # The SYM1 axial components are likewise the real current
    # ``A^a=-i*bar(psi)*gamma_VC^a*gamma_star*psi``.  The explicit ``-i``
    # below is therefore a representation conversion, not a fitted phase.
    operator = _matrix_scale(
        -1j * _HD_COEFFICIENT * coupling,
        _matrix_product(_covariant_axial_gamma_operator(gamma, axial["A_hat_0"], axial["A_hat_r"]), gamma_star),
    )
    plus, minus = _vc_pair_spinors(f_value, g_value, C, angle_theta, angle_phi)
    direct_plus = _matrix_vector(operator, plus)
    direct_minus = _matrix_vector(operator, minus)
    plus_basis = _vc_opposite_chiral_basis(C, angle_theta, angle_phi, plus_sign=True)
    minus_basis = _vc_opposite_chiral_basis(C, angle_theta, angle_phi, plus_sign=False)
    plus_f, plus_g, plus_residual = _project_two_component(direct_plus, plus_basis)
    minus_f, minus_g, minus_residual = _project_two_component(direct_minus, minus_basis)
    expected_f = -1.5 * coupling * C * C * abs(g_value) ** 2 * f_value
    expected_g = -1.5 * coupling * C * C * abs(f_value) ** 2 * g_value
    coefficient_residual = max(
        abs(plus_f - expected_f), abs(plus_g - expected_g), abs(minus_f - expected_f), abs(minus_g - expected_g)
    )
    scale = max(1.0, abs(expected_f), abs(expected_g))
    return {
        "normalization": C,
        "A_hat_0": axial["A_hat_0"],
        "A_hat_r": axial["A_hat_r"],
        "axial_squared": axial["axial_squared"],
        "signature_bridge_closed": bridge["signature_bridge_closed"],
        "equation_translation_residual": bridge["canonical_lhs_to_vc_equation_residual"],
        "plus_retained_F": plus_f,
        "plus_retained_G": plus_g,
        "minus_retained_F": minus_f,
        "minus_retained_G": minus_g,
        "expected_retained_F": expected_f,
        "expected_retained_G": expected_g,
        "local_retained_basis_reconstruction_residual": _finite(
            "local retained-basis reconstruction residual", max(plus_residual, minus_residual)
        ),
        "maximum_coefficient_residual": _finite("coefficient residual", coefficient_residual),
        "local_retained_basis_reconstructed": max(plus_residual, minus_residual) <= 1.0e-12 * scale,
        "coefficient_sign_phase_closed": coefficient_residual <= 1.0e-12 * scale,
    }


def vc_action_projection_identity(
    F: Complex,
    G: Complex,
    lapse: Real,
    radial_metric: Real,
    areal_radius: Real,
    kappa: Real,
    angles: Sequence[tuple[Real, Real]],
) -> dict[str, float | bool | str]:
    """Compare direct covariant cubic coefficients with reduced-action sources.

    The weak-form identity is ``dL/dF*=2*lapse*N_F`` and similarly for G,
    where ``N`` is the retained covariant Hehl--Datta coefficient.  This is not
    a claim about the normalization of a time-evolution equation.
    """

    if not isinstance(angles, Sequence) or isinstance(angles, (str, bytes)) or not angles:
        raise ValueError("angles must be a nonempty sequence of (theta, phi) pairs")
    f_value = _complex("F", F)
    g_value = _complex("G", G)
    alpha = _positive("lapse", lapse)
    sources = vc_contact_action_sources(f_value, g_value, alpha, radial_metric, areal_radius, kappa)
    maximum_identity = 0.0
    maximum_omitted = 0.0
    maximum_coefficient = 0.0
    for item in angles:
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)) or len(item) != 2:
            raise ValueError("each angle entry must be a (theta, phi) pair")
        projection = vc_hehl_datta_projection(f_value, g_value, radial_metric, areal_radius, kappa, item[0], item[1])
        pairs = (
            (projection["plus_retained_F"], sources["source_F"]),
            (projection["plus_retained_G"], sources["source_G"]),
            (projection["minus_retained_F"], sources["source_F"]),
            (projection["minus_retained_G"], sources["source_G"]),
        )
        for coefficient, target in pairs:
            derived = 2.0 * alpha * coefficient
            maximum_identity = max(maximum_identity, abs(derived - target))
        maximum_omitted = max(
            maximum_omitted,
            projection["local_retained_basis_reconstruction_residual"],
        )
        maximum_coefficient = max(maximum_coefficient, projection["maximum_coefficient_residual"])
    scale = max(1.0, abs(sources["source_F"]), abs(sources["source_G"]))
    return {
        "maximum_weak_identity_residual": _finite("weak identity residual", maximum_identity),
        "maximum_tested_angle_reconstruction_residual": _finite(
            "tested-angle reconstruction residual", maximum_omitted
        ),
        "maximum_coefficient_residual": _finite("identity coefficient residual", maximum_coefficient),
        "weak_form_identity_closed": maximum_identity <= 1.0e-12 * scale,
        "all_tested_angle_reconstructions_closed": maximum_omitted <= 1.0e-12 * scale
        and maximum_coefficient <= 1.0e-12 * scale,
        "scope": (
            "Weak action/projection identity and pointwise retained-basis reconstruction "
            "at declared tested angles only; not a global spherical-harmonic closure theorem "
            "or time-evolution normalization."
        ),
    }


def contact_tetrad_stress(
    kappa: Real,
    axial_squared: Real,
    metric_covariant: Sequence[Sequence[Real]],
) -> tuple[tuple[float, ...], ...]:
    """Return the torsion-eliminated contact-sector tetrad stress.

    In the VC ``(-+++)`` layer the reduced interaction is
    ``L_4=(3*kappa/16)*A^I*A_I``.  The axial current carries Lorentz indices,
    so ``A^I*A_I`` is tetrad independent when the spinor components are held
    fixed.  Varying the reduced contact action therefore gives
    ``T_mn=L_4*g_mn``.  This is the complete variation of the algebraic contact
    sector, but it is not the complete ECD stress: the torsion-free Dirac
    kinetic and mass terms remain separate.
    """

    coupling = _positive("kappa", kappa)
    invariant = _real("axial_squared", axial_squared)
    metric = _metric4(metric_covariant)
    coefficient = minimal_ecd_contact_scalar(coupling, invariant)
    return tuple(tuple(_finite("contact stress component", coefficient * value) for value in row) for row in metric)


def metric_null_contact_contraction(
    kappa: Real,
    axial_squared: Real,
    metric_covariant: Sequence[Sequence[Real]],
    null_vector: Sequence[Real],
    tolerance: Real = 1.0e-12,
) -> dict[str, float | bool | str]:
    """Validate a metric-null vector and contract the physical contact stress."""

    metric = _metric4(metric_covariant)
    vector = _vector4("null_vector", null_vector)
    tolerance_value = _positive("tolerance", tolerance)
    norm = sum(metric[row][column] * vector[row] * vector[column] for row in range(4) for column in range(4))
    scale = max(1.0, max(abs(value) for row in metric for value in row) * sum(value * value for value in vector))
    if abs(norm) > tolerance_value * scale:
        raise ValueError("null_vector is not metric-null within tolerance")
    stress = contact_tetrad_stress(kappa, axial_squared, metric)
    contraction = sum(stress[row][column] * vector[row] * vector[column] for row in range(4) for column in range(4))
    return {
        "metric_norm": _finite("metric norm", norm),
        "contact_null_contraction": _finite("contact null contraction", contraction),
        "validated_metric_null": True,
        "contact_term_alone_derives_null_defocusing": False,
        "contact_sector_tetrad_variation_complete": True,
        "full_ecd_stress_derived": False,
    }


def required_nonclaims() -> dict[str, bool]:
    """Return explicit GMF-1B-ECD-INT1 scope boundaries."""

    return {
        "time_evolution_normalization_derived": False,
        "full_ecd_stress_verified": False,
        "full_ecd_effective_stress_verified": False,
        "complete_effective_stress_tensor_derived": False,
        "full_ecd_tetrad_dependence_varied": False,
        "torsion_free_dirac_stress_derived": False,
        "full_ecd_source_so3_equivariance_verified": False,
        "einstein_equations_solved": False,
        "constraint_solution_constructed": False,
        "null_constraint_solution_constructed": False,
        "constraint_compatible_vacuum_region_verified": False,
        "regular_centre_constructed": False,
        "regular_centre_series_verified": False,
        "curvature_invariants_verified": False,
        "metric_null_defocusing_derived": False,
        "bounce_constructed": False,
        "horizon_regular_evolution_constructed": False,
        "finite_mass_exterior_constructed": False,
        "finite_mass_exterior_or_match_proven": False,
        "radial_stability_proven": False,
        "global_solution_constructed": False,
        "child_topology_proven": False,
        "external_Q_derived": False,
        "late_dark_energy_derived": False,
        "variable_c_derived": False,
        "quantum_ensemble_equivalence_proven": False,
    }


def _covariant_axial_gamma_operator(gamma: tuple[tuple[tuple[complex, ...], ...], ...], A0: float, Ar: float) -> tuple[tuple[complex, ...], ...]:
    return tuple(tuple(-A0 * gamma[0][row][column] + Ar * gamma[3][row][column] for column in range(4)) for row in range(4))


def _raw_vc_axial_bilinear(
    spinor: tuple[complex, ...],
    gamma: tuple[tuple[complex, ...], ...],
    gamma_star: tuple[tuple[complex, ...], ...],
) -> complex:
    """Return the unconverted imaginary VC bilinear before the ``-i`` map."""

    hermitizer = _matrix_scale(-1j, _vc_gamma_matrices()[0])
    adjoint = _matrix_vector(hermitizer, spinor)
    transformed = _matrix_vector(_matrix_product(gamma, gamma_star), spinor)
    return sum(adjoint[index].conjugate() * transformed[index] for index in range(4))


def _clifford_residual(
    gamma: tuple[tuple[tuple[complex, ...], ...], ...],
    diagonal_signature: tuple[float, float, float, float],
) -> float:
    """Return the maximum residual in ``{gamma^a,gamma^b}=2 eta^ab``."""

    residual = 0.0
    for first in range(4):
        for second in range(4):
            forward = _matrix_product(gamma[first], gamma[second])
            reverse = _matrix_product(gamma[second], gamma[first])
            target = 2.0 * diagonal_signature[first] if first == second else 0.0
            for row in range(4):
                for column in range(4):
                    expected = target if row == column else 0.0
                    residual = max(
                        residual,
                        abs(forward[row][column] + reverse[row][column] - expected),
                    )
    return _finite("Clifford residual", residual)


def _matrix_difference_residual(
    first: tuple[tuple[complex, ...], ...],
    second: tuple[tuple[complex, ...], ...],
) -> float:
    """Return the maximum componentwise difference of two 4 by 4 matrices."""

    return _finite(
        "matrix difference residual",
        max(abs(first[row][column] - second[row][column]) for row in range(4) for column in range(4)),
    )


def _vc_pair_spinors(F: complex, G: complex, C: float, theta: float, phi: float) -> tuple[tuple[complex, ...], tuple[complex, ...]]:
    phase_plus = complex(cos(phi / 2.0), sin(phi / 2.0))
    phase_minus = complex(cos(phi / 2.0), -sin(phi / 2.0))
    sin_half = sin(theta / 2.0)
    cos_half = cos(theta / 2.0)
    plus = tuple(C * phase_plus * item for item in (F * sin_half, G * cos_half, F * sin_half, G * cos_half))
    minus = tuple(C * phase_minus * item for item in (F * cos_half, -G * sin_half, F * cos_half, -G * sin_half))
    return plus, minus


def _vc_opposite_chiral_basis(C: float, theta: float, phi: float, plus_sign: bool) -> tuple[tuple[complex, ...], tuple[complex, ...]]:
    phase = complex(cos(phi / 2.0), sin(phi / 2.0) if plus_sign else -sin(phi / 2.0))
    sin_half = sin(theta / 2.0)
    cos_half = cos(theta / 2.0)
    if plus_sign:
        basis_f = (sin_half, 0j, -sin_half, 0j)
        basis_g = (0j, cos_half, 0j, -cos_half)
    else:
        basis_f = (cos_half, 0j, -cos_half, 0j)
        basis_g = (0j, -sin_half, 0j, sin_half)
    return tuple(tuple(C * phase * value for value in basis) for basis in (basis_f, basis_g))


def _project_two_component(vector: tuple[complex, ...], basis: tuple[tuple[complex, ...], tuple[complex, ...]]) -> tuple[complex, complex, float]:
    coefficients: list[complex] = []
    reconstruction = [0j, 0j, 0j, 0j]
    for element in basis:
        norm = sum(value.conjugate() * value for value in element).real
        if norm <= 0.0:
            raise ValueError("retained basis has zero norm")
        coefficient = sum(element[index].conjugate() * vector[index] for index in range(4)) / norm
        coefficients.append(coefficient)
        for index in range(4):
            reconstruction[index] += coefficient * element[index]
    residual = max(abs(vector[index] - reconstruction[index]) for index in range(4))
    return coefficients[0], coefficients[1], _finite("projection residual", residual)


def _metric4(metric: Sequence[Sequence[Real]]) -> tuple[tuple[float, ...], ...]:
    if not isinstance(metric, Sequence) or isinstance(metric, (str, bytes)) or len(metric) != 4:
        raise ValueError("metric_covariant must be a 4 by 4 finite real matrix")
    rows = tuple(_vector4(f"metric_covariant[{index}]", row) for index, row in enumerate(metric))
    for row in range(4):
        for column in range(4):
            if abs(rows[row][column] - rows[column][row]) > 1.0e-12:
                raise ValueError("metric_covariant must be symmetric")
    return rows


def _vector4(name: str, vector: Sequence[Real]) -> tuple[float, ...]:
    if not isinstance(vector, Sequence) or isinstance(vector, (str, bytes)) or len(vector) != 4:
        raise ValueError(f"{name} must have four finite real components")
    return tuple(_real(f"{name}[{index}]", value) for index, value in enumerate(vector))
