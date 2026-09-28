"""GMF-1B-ECD-SYM1: minimal ECD spherical axial-current preflight.

This module evaluates the Ventrella--Choptuik fixture in its published
``(-,+,+,+)`` imaginary gamma representation.  Its real axial current is the
complete convention translation of a canonical ``(+,-,-,-)`` minimal
Einstein--Cartan--Dirac action, as certified by GMF-1B-ECD-INT1.  It checks two
narrow axial-current controls: the equal-profile Finster-type
classical pair has a cancelling total axial current, whereas the massless
Ventrella--Choptuik chiral doublet has an SO(3)-compatible nonzero axial
current when both radial amplitudes are populated.

It is not a proof of full ECD source equivariance, a stress-tensor calculation,
a Hehl--Datta angular-closure
proof, a curvature calculation, a constraint solver, or an evolution/global
black-hole-to-child solution.  The fields here are classical c-number spinor
fields; this does not establish equivalence to a quantum ensemble.
"""

from __future__ import annotations

from math import cos, exp, isfinite, pi, sin, sqrt
from numbers import Complex, Real
from typing import Final


_THREE_KAPPA_OVER_SIXTEEN: Final[float] = 3.0 / 16.0


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


def geometric_normalization(areal_radius: Real, radial_metric: Real) -> float:
    """Return ``C=1/(2*sqrt(pi)*r*sqrt(a))`` for the VC chiral pair.

    ``r`` is an areal radius and ``a`` is the positive polar-areal radial
    metric coefficient.  The formula is not a horizon-regular gauge choice;
    it is retained only to certify the published polar-areal spinor ansatz.
    """

    radius = _positive("areal_radius", areal_radius)
    metric = _positive("radial_metric", radial_metric)
    result = _finite("geometric normalization", 1.0 / (2.0 * sqrt(pi) * radius * sqrt(metric)))
    if result <= 0.0:
        raise ValueError("geometric normalization is outside the finite nonzero range")
    return result


def finster_equal_profile_axial_control(
    upper_amplitude: Complex,
    lower_amplitude: Complex,
    normalization: Real = 1.0,
    theta: Real = 0.0,
    phi: Real = 0.0,
) -> dict[str, float | bool | str]:
    """Directly evaluate axial cancellation for the Finster Pauli doublet.

    The primary ansatz is ``Psi_a=(u e_a, v sigma^r e_a)^T`` for ``a=1,2``
    (its common stationary phase cancels from bilinears).  This routine uses
    the source's ``(+,-,-,-)`` Dirac-Pauli matrices to evaluate the pair
    rather than returning zero by construction.  The result follows from the
    complete Pauli-doublet identity ``sum_a e_a e_a^dagger=I`` and
    ``tr(sigma_i)=0``.  It is not a theorem that every spherical Dirac
    multiplet has a vanishing axial current.
    """

    upper = _complex("upper_amplitude", upper_amplitude)
    lower = _complex("lower_amplitude", lower_amplitude)
    scale = _positive("normalization", normalization)
    angle_theta = _real("theta", theta)
    angle_phi = _real("phi", phi)
    gamma = _finster_gamma_matrices()
    gamma5 = _matrix_scale(1j, _matrix_product(gamma[0], gamma[1], gamma[2], gamma[3]))
    sigma_r = _pauli_radial(angle_theta, angle_phi)
    e1 = (1 + 0j, 0j)
    e2 = (0j, 1 + 0j)

    def spinor(basis: tuple[complex, complex]) -> tuple[complex, ...]:
        lower_two = (
            lower * (sigma_r[0][0] * basis[0] + sigma_r[0][1] * basis[1]),
            lower * (sigma_r[1][0] * basis[0] + sigma_r[1][1] * basis[1]),
        )
        return (scale * upper * basis[0], scale * upper * basis[1], scale * lower_two[0], scale * lower_two[1])

    spinors = (spinor(e1), spinor(e2))
    individual_axial = tuple(
        tuple(_finster_axial_bilinear(field, gamma[index], gamma5) for index in range(4))
        for field in spinors
    )
    axial = tuple(sum(current[index] for current in individual_axial) for index in range(4))
    residual = max(abs(value) for component in axial for value in (component.real, component.imag))
    first_current_size = max(abs(component) for component in individual_axial[0])
    first_current_imaginary_residual = max(abs(component.imag) for component in individual_axial[0])
    populated = abs(upper) > 0.0 or abs(lower) > 0.0
    return {
        "construction": "equal_profile_classical_finster_control",
        "gamma5_definition": "i*gamma0*gamma1*gamma2*gamma3",
        "evaluation_basis": "Cartesian_orthonormal_tetrad",
        "normalization": scale,
        "profiles_populated": populated,
        "A_hat_0": _finite("Finster A_hat_0", axial[0].real),
        "A_hat_x": _finite("Finster A_hat_x", axial[1].real),
        "A_hat_y": _finite("Finster A_hat_y", axial[2].real),
        "A_hat_z": _finite("Finster A_hat_z", axial[3].real),
        "axial_squared": 0.0,
        "axial_current_nonzero": False,
        "algebraic_torsion_source_nonzero": False,
        "axial_contact_invariant_nonzero": False,
        "single_spinor_1_axial_current_nonzero": first_current_size > 1.0e-14,
        "single_spinor_1_max_abs_axial_component": _finite(
            "Finster single-spinor axial size", first_current_size
        ),
        "single_spinor_1_max_imaginary_residual": _finite(
            "Finster single-spinor imaginary residual", first_current_imaginary_residual
        ),
        "maximum_direct_axial_residual": _finite("Finster axial residual", residual),
        "pauli_doublet_cancellation_certified_to_binary64_roundoff": residual <= 1.0e-12 * max(1.0, scale * scale * (abs(upper) ** 2 + abs(lower) ** 2)),
        "scope": "Specific equal-profile classical-pair cancellation control only.",
    }


def ventrella_chiral_axial_bilinears(
    F: Complex,
    G: Complex,
    areal_radius: Real,
    radial_metric: Real,
) -> dict[str, float | bool | str]:
    """Evaluate the SO(3)-compatible axial current of the VC chiral doublet.

    In the VC ``(-,+,+,+)`` orthonormal-frame convention, the translated real
    current is ``A^a=-i*bar(psi)*gamma_VC^a*gamma_star*psi``.  The published
    massless left-chiral pair gives ``A0=2*C^2*(|F|^2+|G|^2)`` and
    ``Ar=2*C^2*(|F|^2-|G|^2)``, with no angular components.  Hence
    ``A_I A^I=-16*C^4*|F|^2*|G|^2``.  This establishes SO(3), not O(3),
    compatibility of the axial-current sector; it does not compute the full
    tetrad variation of the ECD contact stress.
    """

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    C = geometric_normalization(areal_radius, radial_metric)
    f_squared = abs(f_value) ** 2
    g_squared = abs(g_value) ** 2
    C_squared = C * C
    A0 = _finite("A_hat_0", 2.0 * C_squared * (f_squared + g_squared))
    Ar = _finite("A_hat_r", 2.0 * C_squared * (f_squared - g_squared))
    axial_squared = _finite("axial_squared", -16.0 * C_squared * C_squared * f_squared * g_squared)
    axial_current_nonzero = f_squared > 0.0 or g_squared > 0.0
    contact_invariant_nonzero = f_squared > 0.0 and g_squared > 0.0
    return {
        "construction": "ventrella_choptuik_massless_left_chiral_doublet",
        "canonical_action_signature": "+---",
        "numerical_vc_signature": "-+++",
        "real_current_definition": "A^a=-i*bar(psi)*gamma_VC^a*gamma_star*psi",
        "normalization": C,
        "F_abs_squared": f_squared,
        "G_abs_squared": g_squared,
        "A_hat_0": A0,
        "A_hat_r": Ar,
        "A_hat_theta": 0.0,
        "A_hat_phi": 0.0,
        "axial_squared": axial_squared,
        "both_radial_amplitudes_populated": contact_invariant_nonzero,
        "axial_current_nonzero": axial_current_nonzero,
        "algebraic_torsion_source_nonzero": axial_current_nonzero,
        "axial_contact_invariant_nonzero": contact_invariant_nonzero,
        "contact_action_density_nonzero": contact_invariant_nonzero,
        "so3_axial_current_compatible": True,
        "o3_parity_compatibility_proven": False,
        "classical_c_number_field_interpretation": True,
        "full_ecd_contact_stress_verified": False,
    }


def minimal_ecd_contact_scalar(kappa: Real, axial_squared: Real) -> float:
    """Return the VC-signature interaction ``(3*kappa/16)*A_I*A^I``.

    This is the translation of the canonical ``(+---)`` interaction
    ``-(3*kappa/16)*J_I*J^I`` with ``J_+^2=-A_-^2``; see INT1.
    """

    coupling = _positive("kappa", kappa)
    invariant = _real("axial_squared", axial_squared)
    return _finite("minimal ECD contact scalar", _THREE_KAPPA_OVER_SIXTEEN * coupling * invariant)


def double_null_metric_null_cone(
    inverse_conformal_factor: Real,
    xi_u: Real,
    xi_v: Real,
) -> dict[str, float | bool | str]:
    """Certify the minimal ECD matter principal cone in a double-null base.

    For ``ds^2=-2*exp(-f) du dv+R^2 dOmega^2``, pass
    ``inverse_conformal_factor=exp(f)>0``.  The radial characteristic form is
    ``g^AB xi_A xi_B=-2*exp(f)*xi_u*xi_v``.  Algebraic Cartan torsion and the
    Hehl--Datta contact term are lower order, so they do not alter this cone.
    """

    factor = _positive("inverse_conformal_factor", inverse_conformal_factor)
    u_value = _real("xi_u", xi_u)
    v_value = _real("xi_v", xi_v)
    quadratic = _finite("double-null metric quadratic", -2.0 * factor * u_value * v_value)
    determinant_factor = _finite("Dirac principal determinant factor", quadratic * quadratic)
    return {
        "inverse_conformal_factor": factor,
        "metric_null_quadratic": quadratic,
        "dirac_principal_determinant_factor": determinant_factor,
        "metric_null_characteristic": quadratic == 0.0,
        "algebraic_torsion_changes_metric_cone": False,
        "scope": "Matter principal-symbol statement only; no gauge/Einstein constraint hyperbolicity proof.",
    }


def normalized_annular_bump(radius: Real, inner_radius: Real, outer_radius: Real) -> float:
    """Return a normalized ``C-infinity`` bump with support ``(inner, outer)``.

    It is exactly zero at and outside the specified radii and equals one at
    their midpoint.  It is a profile generator, not a solved ECD initial-data
    constraint construction.
    """

    r_value = _real("radius", radius)
    inner = _nonnegative("inner_radius", inner_radius)
    outer = _positive("outer_radius", outer_radius)
    if outer <= inner:
        raise ValueError("outer_radius must be greater than inner_radius")
    if r_value <= inner or r_value >= outer:
        return 0.0
    x = (r_value - inner) / (outer - inner)
    exponent = -1.0 / (x * (1.0 - x)) + 4.0
    return _finite("normalized annular bump", exp(exponent))


def chiral_annular_packet(
    radius: Real,
    inner_radius: Real,
    outer_radius: Real,
    F_amplitude: Complex,
    G_amplitude: Complex,
) -> dict[str, float | bool | str]:
    """Construct a JSON-safe common-profile smooth annular VC packet point."""

    bump = normalized_annular_bump(radius, inner_radius, outer_radius)
    F0 = _complex("F_amplitude", F_amplitude)
    G0 = _complex("G_amplitude", G_amplitude)
    F_value = bump * F0
    G_value = bump * G0
    return {
        "radius": _real("radius", radius),
        "bump": bump,
        "F_real": F_value.real,
        "F_imag": F_value.imag,
        "G_real": G_value.real,
        "G_imag": G_value.imag,
        "F_abs_squared": abs(F_value) ** 2,
        "G_abs_squared": abs(G_value) ** 2,
        "inside_open_support": bump > 0.0,
        "profile_contract": "C-infinity compact annular profile; constraints and evolution are not solved.",
    }


def audit_chiral_annular_packet(
    radius: Real,
    inner_radius: Real,
    outer_radius: Real,
    F_amplitude: Complex,
    G_amplitude: Complex,
    radial_metric: Real,
    kappa: Real,
) -> dict[str, object]:
    """Audit one annular packet point for taper and axial-contact semantics."""

    if _nonnegative("inner_radius", inner_radius) <= 0.0:
        raise ValueError("audit requires a strictly positive annular inner_radius")
    packet = chiral_annular_packet(radius, inner_radius, outer_radius, F_amplitude, G_amplitude)
    F_value = complex(packet["F_real"], packet["F_imag"])
    G_value = complex(packet["G_real"], packet["G_imag"])
    axial = ventrella_chiral_axial_bilinears(F_value, G_value, radius, radial_metric)
    contact = minimal_ecd_contact_scalar(kappa, axial["axial_squared"])
    outside = not packet["inside_open_support"]
    return {
        "packet": packet,
        "axial": axial,
        "minimal_ecd_contact_scalar": contact,
        "outside_support": outside,
        "axial_current_vanishes_at_this_outside_support_point": (
            outside and axial["axial_current_nonzero"] is False
        ),
        "smooth_taper_contract": "True only for a profile used as a C-infinity function, not a proof of constrained spacetime matching.",
        "curvature_invariants_verified": False,
        "full_ecd_stress_verified": False,
        "hehl_datta_angular_closure_verified": False,
        "evolution_solution_constructed": False,
    }


def required_nonclaims() -> dict[str, bool]:
    """Return explicit GMF-1B-ECD-SYM1 scope boundaries."""

    return {
        "full_ecd_contact_stress_verified": False,
        "full_ecd_source_so3_equivariance_verified": False,
        "hehl_datta_angular_closure_verified": False,
        "curvature_invariants_verified": False,
        "constraint_solution_constructed": False,
        "constraint_compatible_vacuum_region_verified": False,
        "horizon_regular_evolution_constructed": False,
        "global_solution_constructed": False,
        "child_topology_proven": False,
        "external_Q_derived": False,
        "late_dark_energy_derived": False,
        "variable_c_derived": False,
        "quantum_ensemble_equivalence_proven": False,
    }


def _nonnegative(name: str, value: Real) -> float:
    result = _real(name, value)
    if result < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return result


def ventrella_gamma_matrix_certificate(
    F: Complex,
    G: Complex,
    areal_radius: Real,
    radial_metric: Real,
    theta: Real,
    phi: Real,
) -> dict[str, float | bool | str]:
    """Direct complex-4x4 certificate of the VC axial-current formulas.

    The matrices are the imaginary ``(-,+,+,+)`` representation used in the
    source pair.  The calculation certifies bilinears only; it does not
    calculate a full ECD tetrad variation or prove angular closure.
    """

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    C = geometric_normalization(areal_radius, radial_metric)
    angle_theta = _real("theta", theta)
    angle_phi = _real("phi", phi)
    gamma = _vc_gamma_matrices()
    gamma_star = _matrix_scale(1j, _matrix_product(gamma[0], gamma[1], gamma[2], gamma[3]))
    phase_plus = complex(cos(angle_phi / 2.0), sin(angle_phi / 2.0))
    phase_minus = complex(cos(angle_phi / 2.0), -sin(angle_phi / 2.0))
    sin_half = sin(angle_theta / 2.0)
    cos_half = cos(angle_theta / 2.0)
    plus = tuple(C * phase_plus * item for item in (f_value * sin_half, g_value * cos_half, f_value * sin_half, g_value * cos_half))
    minus = tuple(C * phase_minus * item for item in (f_value * cos_half, -g_value * sin_half, f_value * cos_half, -g_value * sin_half))
    axial = tuple(
        _axial_bilinear(plus, gamma[index], gamma_star) + _axial_bilinear(minus, gamma[index], gamma_star)
        for index in range(4)
    )
    expected = ventrella_chiral_axial_bilinears(f_value, g_value, areal_radius, radial_metric)
    residuals = (
        axial[0].real - expected["A_hat_0"],
        axial[1].real,
        axial[2].real,
        axial[3].real - expected["A_hat_r"],
        axial[0].imag,
        axial[1].imag,
        axial[2].imag,
        axial[3].imag,
    )
    maximum_residual = max(abs(value) for value in residuals)
    return {
        "representation": "VC_imaginary_gamma_minus_plus_plus_plus",
        "gamma_star_eigenvalue_for_pair": 1.0,
        "A_hat_0_direct": _finite("direct A_hat_0", axial[0].real),
        "A_hat_r_direct": _finite("direct A_hat_r", axial[3].real),
        "A_hat_theta_direct": _finite("direct A_hat_theta", axial[2].real),
        "A_hat_phi_direct": _finite("direct A_hat_phi", axial[1].real),
        "maximum_formula_residual": _finite("gamma certificate residual", maximum_residual),
        "formula_certified_to_binary64_roundoff": maximum_residual <= 1.0e-12 * max(1.0, abs(expected["A_hat_0"])),
        "scope": "Axial-bilinear certificate only; no full ECD stress or angular-closure claim.",
    }


def _vc_gamma_matrices() -> tuple[tuple[tuple[complex, ...], ...], ...]:
    i = 1j
    zero = 0j
    sigma1 = ((zero, 1), (1, zero))
    sigma2 = ((zero, -i), (i, zero))
    sigma3 = ((1, zero), (zero, -1))
    gamma0 = ((i, zero, zero, zero), (zero, i, zero, zero), (zero, zero, -i, zero), (zero, zero, zero, -i))

    def spatial(sigma: tuple[tuple[complex, complex], tuple[complex, complex]]) -> tuple[tuple[complex, ...], ...]:
        return (
            (zero, zero, i * sigma[0][0], i * sigma[0][1]),
            (zero, zero, i * sigma[1][0], i * sigma[1][1]),
            (-i * sigma[0][0], -i * sigma[0][1], zero, zero),
            (-i * sigma[1][0], -i * sigma[1][1], zero, zero),
        )

    # VC uses frame order (0,1,2,3), with radial gamma equal to gamma^3.
    return (gamma0, spatial(sigma1), spatial(sigma2), spatial(sigma3))


def _finster_gamma_matrices() -> tuple[tuple[tuple[complex, ...], ...], ...]:
    """Return the ``(+,-,-,-)`` Dirac-Pauli matrices used by Finster et al."""

    zero = 0j
    sigma1 = ((zero, 1), (1, zero))
    sigma2 = ((zero, -1j), (1j, zero))
    sigma3 = ((1, zero), (zero, -1))
    gamma0 = ((1, zero, zero, zero), (zero, 1, zero, zero), (zero, zero, -1, zero), (zero, zero, zero, -1))

    def spatial(sigma: tuple[tuple[complex, complex], tuple[complex, complex]]) -> tuple[tuple[complex, ...], ...]:
        return (
            (zero, zero, sigma[0][0], sigma[0][1]),
            (zero, zero, sigma[1][0], sigma[1][1]),
            (-sigma[0][0], -sigma[0][1], zero, zero),
            (-sigma[1][0], -sigma[1][1], zero, zero),
        )

    return (gamma0, spatial(sigma1), spatial(sigma2), spatial(sigma3))


def _pauli_radial(theta: float, phi: float) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    """Return ``sigma dot r_hat`` in a fixed Cartesian Pauli basis."""

    nx = sin(theta) * cos(phi)
    ny = sin(theta) * sin(phi)
    nz = cos(theta)
    return ((nz, complex(nx, -ny)), (complex(nx, ny), -nz))


def _matrix_product(*matrices: tuple[tuple[complex, ...], ...]) -> tuple[tuple[complex, ...], ...]:
    result = matrices[0]
    for matrix in matrices[1:]:
        result = tuple(
            tuple(sum(result[row][index] * matrix[index][column] for index in range(4)) for column in range(4))
            for row in range(4)
        )
    return result


def _matrix_scale(scale: complex, matrix: tuple[tuple[complex, ...], ...]) -> tuple[tuple[complex, ...], ...]:
    return tuple(tuple(scale * value for value in row) for row in matrix)


def _matrix_vector(matrix: tuple[tuple[complex, ...], ...], vector: tuple[complex, ...]) -> tuple[complex, ...]:
    return tuple(sum(matrix[row][column] * vector[column] for column in range(4)) for row in range(4))


def _axial_bilinear(
    spinor: tuple[complex, ...],
    gamma: tuple[tuple[complex, ...], ...],
    gamma_star: tuple[tuple[complex, ...], ...],
) -> complex:
    # In the source representation, bar(psi)=psi^dagger*(-i*gamma^0).
    # Its gamma matrices are imaginary, so the real physical bilinear carries
    # the corresponding -i conversion from this anti-Hermitian convention.
    adjoint_metric = _matrix_scale(-1j, _vc_gamma_matrices()[0])
    transformed = _matrix_vector(_matrix_product(gamma, gamma_star), spinor)
    adjoint = _matrix_vector(adjoint_metric, spinor)
    return -1j * sum(adjoint[index].conjugate() * transformed[index] for index in range(4))


def _finster_axial_bilinear(
    spinor: tuple[complex, ...],
    gamma: tuple[tuple[complex, ...], ...],
    gamma5: tuple[tuple[complex, ...], ...],
) -> complex:
    """Return ``bar(psi) gamma^I gamma5 psi`` in the Finster convention."""

    transformed = _matrix_vector(_matrix_product(gamma, gamma5), spinor)
    # Here gamma^0 is Hermitian and bar(psi)=psi^dagger gamma^0.
    adjoint = _matrix_vector(_finster_gamma_matrices()[0], spinor)
    return sum(adjoint[index].conjugate() * transformed[index] for index in range(4))
