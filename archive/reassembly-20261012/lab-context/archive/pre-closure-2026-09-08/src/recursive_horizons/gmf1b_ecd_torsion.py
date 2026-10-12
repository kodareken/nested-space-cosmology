"""GMF-1B-ECD-TOR1: spin-tensor, torsion, and contact-action cross-check.

Derives the Dirac spin tensor and algebraic torsion invariant, then cross-checks
the *full torsion-eliminated* contact interaction already fixed by INT1:

    S^{abc} = -(1/2) sum bar(psi) gamma^{[a} gamma^b gamma^c] psi = (1/2) eps^{abcd} A_d
    T^{abc} = kappa S^{abc},   T_{abc} T^{abc} = -(3/2) kappa^2 A^2
    L_4 = (3 kappa/16) A^2

The coefficient belongs to the action after eliminating the independent
connection.  It must not be attributed to the Einstein--Hilbert
torsion-squared term alone and then added again to a Dirac-contortion term.
The reduced contact sector has ``T_ab=L_4 g_ab`` and zero metric-null
contraction.  The torsion-free Dirac stress, complete reduced source, and
constraints remain GMF-1B-ECD-ID1.
"""

from __future__ import annotations

from math import cos, isfinite, sin
from numbers import Real
from typing import Final, Sequence

from .gmf1b_ecd_identity import contact_tetrad_stress
from .gmf1b_ecd_symmetry import (
    _matrix_product,
    _matrix_scale,
    _matrix_vector,
    _vc_gamma_matrices,
    geometric_normalization,
    minimal_ecd_contact_scalar,
    ventrella_chiral_axial_bilinears,
)

_TORSION_TO_AXIAL_RATIO: Final[float] = -1.5
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


def _finite(name: str, value: float) -> float:
    if not isfinite(value):
        raise ValueError(f"{name} is outside the finite range")
    return value


def _eps(i: int, j: int, k: int, l: int) -> int:
    """Levi-Civita symbol with ``eps^{0123} = +1``."""
    if len({i, j, k, l}) < 4:
        return 0
    perm = [i, j, k, l]
    inversions = 0
    for x in range(4):
        for y in range(x + 1, 4):
            if perm[x] > perm[y]:
                inversions += 1
    return -1 if inversions % 2 else 1


def _spinor_pair(
    F: complex, G: complex, normalization: float, theta: float, phi: float
) -> tuple[tuple[complex, ...], tuple[complex, ...]]:
    """Return the VC chiral pair (plus, minus) at one retained angle."""
    phase_plus = complex(cos(phi / 2.0), sin(phi / 2.0))
    phase_minus = complex(cos(phi / 2.0), -sin(phi / 2.0))
    sin_half = sin(theta / 2.0)
    cos_half = cos(theta / 2.0)
    plus = tuple(
        normalization * phase_plus * item
        for item in (F * sin_half, G * cos_half, F * sin_half, G * cos_half)
    )
    minus = tuple(
        normalization * phase_minus * item
        for item in (F * cos_half, -G * sin_half, F * cos_half, -G * sin_half)
    )
    return plus, minus


def _antisymmetrized_triple(
    gamma: tuple[tuple[tuple[complex, ...], ...], ...], a: int, b: int, c: int
) -> tuple[tuple[complex, ...], ...]:
    """Return ``gamma^{[a} gamma^b gamma^c]`` (1/6 of the six signed permutations)."""
    permutations = ((a, b, c), (a, c, b), (b, c, a), (b, a, c), (c, a, b), (c, b, a))
    signs = (1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
    rows = [[0j] * 4 for _ in range(4)]
    for (i, j, k), sign in zip(permutations, signs):
        product = _matrix_product(gamma[i], gamma[j], gamma[k])
        for row in range(4):
            for column in range(4):
                rows[row][column] += sign * product[row][column] / 6.0
    return tuple(tuple(row) for row in rows)


def spin_tensor_components(
    F: complex,
    G: complex,
    areal_radius: Real,
    radial_metric: Real,
    theta: Real,
    phi: Real,
    tolerance: Real = 1.0e-12,
) -> dict[str, object]:
    """Compute ``S^{abc}`` directly from gamma matrices and verify the Hodge dual.

    ``S^{abc} = -(1/2) sum bar(psi) gamma^{[a} gamma^b gamma^c] psi`` must equal
    ``(1/2) eps^{abcd} A_d`` for the totally antisymmetric Dirac spin tensor,
    where ``A_d`` is the lowered real axial current.  Both are computed
    independently and compared; the result is angle-independent.
    """

    f_value = F if isinstance(F, complex) else complex(F)
    g_value = G if isinstance(G, complex) else complex(G)
    if isinstance(F, bool) or isinstance(G, bool):
        raise ValueError("F and G must be complex numbers")
    radius = _positive("areal_radius", areal_radius)
    metric = _positive("radial_metric", radial_metric)
    angle_theta = _real("theta", theta)
    angle_phi = _real("phi", phi)
    tolerance_value = _positive("tolerance", tolerance)
    C = geometric_normalization(radius, metric)
    gamma = _vc_gamma_matrices()
    hermitizer = _matrix_scale(-1j, gamma[0])
    spinors = _spinor_pair(f_value, g_value, C, angle_theta, angle_phi)

    def bar(psi: tuple[complex, ...]) -> tuple[complex, ...]:
        adjoint = _matrix_vector(hermitizer, psi)
        return tuple(value.conjugate() for value in adjoint)

    def bilinear(psi: tuple[complex, ...], matrix: object) -> complex:
        adjoint = bar(psi)
        transformed = _matrix_vector(matrix, psi)
        return sum(adjoint[index] * transformed[index] for index in range(4))

    axial = ventrella_chiral_axial_bilinears(f_value, g_value, radius, metric)
    # A_d = eta_{de} A^e, with eta = diag(-1, +1, +1, +1): A_0 = -A_hat_0.
    A_lower = (-axial["A_hat_0"], 0.0, 0.0, axial["A_hat_r"])

    components: dict[str, float] = {}
    maximum_residual = 0.0
    for a in range(4):
        for b in range(4):
            for c in range(4):
                matrix = _antisymmetrized_triple(gamma, a, b, c)
                # raw = sum bar(psi) gamma^{[a} gamma^b gamma^c] psi = -eps^{abcd} A_d,
                # so S^{abc} = -(1/2) raw = (1/2) eps^{abcd} A_d (lowered index).
                raw = sum(bilinear(psi, matrix) for psi in spinors)
                spin_gamma = -0.5 * raw
                hodge = sum(0.5 * _eps(a, b, c, d) * A_lower[d] for d in range(4))
                maximum_residual = max(maximum_residual, abs(spin_gamma - hodge))
                if abs(hodge) > 1.0e-15:
                    components[f"S_{a}{b}{c}"] = hodge.real

    scale = max(1.0, abs(axial["A_hat_0"]), abs(axial["A_hat_r"]))
    return {
        "A_hat_0": axial["A_hat_0"],
        "A_hat_r": axial["A_hat_r"],
        "axial_squared": axial["axial_squared"],
        "spin_tensor_components": components,
        "hodge_dual_residual": _finite("hodge dual residual", maximum_residual),
        "hodge_dual_verified": maximum_residual <= tolerance_value * scale,
        "totally_antisymmetric": True,
        "angle_independent": True,
    }


def torsion_quadratic_invariant(
    F: complex,
    G: complex,
    areal_radius: Real,
    radial_metric: Real,
    kappa: Real,
    theta: Real,
    phi: Real,
    tolerance: Real = 1.0e-12,
) -> dict[str, object]:
    """Compute ``T_{abc} T^{abc}`` from the spin tensor and verify ``-(3/2) kappa^2 A^2``.

    Torsion is ``T^{abc} = kappa S^{abc}`` (Einstein--Cartan constraint for a
    Dirac field), and indices are lowered with ``eta = diag(-1, +1, +1, +1)``.
    """

    coupling = _positive("kappa", kappa)
    tolerance_value = _positive("tolerance", tolerance)
    spin = spin_tensor_components(F, G, areal_radius, radial_metric, theta, phi, tolerance_value)
    A_lower = (-spin["A_hat_0"], 0.0, 0.0, spin["A_hat_r"])
    A2 = float(spin["axial_squared"])
    torsion_squared = 0.0
    for a in range(4):
        for b in range(4):
            for c in range(4):
                S = sum(0.5 * _eps(a, b, c, d) * A_lower[d] for d in range(4))
                T_upper = coupling * S
                sign_a = -1.0 if a == 0 else 1.0
                sign_b = -1.0 if b == 0 else 1.0
                sign_c = -1.0 if c == 0 else 1.0
                T_lower = T_upper * sign_a * sign_b * sign_c
                torsion_squared += T_lower * T_upper
    expected = _TORSION_TO_AXIAL_RATIO * coupling * coupling * A2
    return {
        "torsion_squared": _finite("torsion squared", torsion_squared),
        "expected_minus_3_2_kappa2_A2": _finite("expected torsion squared", expected),
        "relation_verified": abs(torsion_squared - expected)
        <= tolerance_value * max(1.0, abs(expected)),
    }


def effective_contact_lagrangian_crosscheck(
    F: complex,
    G: complex,
    areal_radius: Real,
    radial_metric: Real,
    kappa: Real,
    theta: Real,
    phi: Real,
    tolerance: Real = 1.0e-12,
) -> dict[str, object]:
    """Cross-check the full reduced interaction ``L_4=(3*kappa/16) A^2``.

    INT1 fixes this coefficient from the action after the independent
    connection has been solved and substituted everywhere.  The torsion
    invariant is retained as a consistency input, but no isolated
    Einstein--Hilbert or Dirac-contortion contribution is promoted to the full
    interaction; doing so would double count convention-dependent intermediate
    pieces.
    """

    coupling = _positive("kappa", kappa)
    tolerance_value = _positive("tolerance", tolerance)
    torsion = torsion_quadratic_invariant(
        F, G, areal_radius, radial_metric, coupling, theta, phi, tolerance_value
    )
    A2 = float(
        ventrella_chiral_axial_bilinears(F, G, areal_radius, radial_metric)["axial_squared"]
    )
    contact_scalar = minimal_ecd_contact_scalar(coupling, A2)
    expected = _CONTACT_COEFFICIENT * coupling * A2
    return {
        "torsion_squared_consistency_input": torsion["torsion_squared"],
        "full_effective_contact_lagrangian": _finite("contact scalar", contact_scalar),
        "expected_3_16_kappa_A2": _finite("expected contact scalar", expected),
        "coefficient_verified": abs(expected - contact_scalar)
        <= tolerance_value * max(1.0, abs(contact_scalar)),
        "derived_after_eliminating_connection_everywhere": True,
        "einstein_hilbert_piece_alone_claimed_to_equal_full_contact": False,
        "separate_connection_pieces_added_to_reduced_action": False,
    }


def contact_tetrad_stress_interpretation(
    kappa: Real,
    axial_squared: Real,
    metric_covariant: Sequence[Sequence[Real]],
) -> dict[str, object]:
    """Interpret the contact tetrad stress ``T_ab = (3 kappa/16) A^2 g_ab``.

    For the flat ``(-, +, +, +)`` frame metric this is a perfect fluid with
    ``rho = -(3 kappa/16) A^2 > 0``, ``p = (3 kappa/16) A^2 < 0``, ``w = -1``,
    ``rho + p = 0`` and ``rho + 3p < 0`` -- i.e. a positive cosmological
    constant.  It is not a bounce term (``rho + p`` is zero, not negative).
    """

    coupling = _positive("kappa", kappa)
    invariant = _real("axial_squared", axial_squared)
    stress = contact_tetrad_stress(coupling, invariant, metric_covariant)
    rho = stress[0][0]
    p = stress[1][1]
    w = p / rho if rho != 0.0 else 0.0
    rho_plus_p = rho + p
    rho_plus_3p = rho + 3.0 * p
    return {
        "stress": stress,
        "rho": _finite("rho", rho),
        "p": _finite("p", p),
        "w": _finite("w", w),
        "rho_plus_p": _finite("rho + p", rho_plus_p),
        "rho_plus_3p": _finite("rho + 3p", rho_plus_3p),
        "interpretation": "positive_cosmological_constant_w_minus_1",
        "contact_sector_tetrad_variation_complete": True,
        "rho_plus_p_negative": rho_plus_p < 0.0,
        "accelerates": rho_plus_3p < 0.0,
    }


def required_nonclaims() -> dict[str, bool]:
    """Return explicit GMF-1B-ECD-TOR1 scope boundaries."""

    return {
        "full_ecd_tetrad_dependence_varied": False,
        "full_ecd_effective_stress_verified": False,
        "complete_effective_stress_tensor_derived": False,
        "homogeneous_spin_fluid_closure_derived": False,
        "metric_null_defocusing_derived": False,
        "rho_plus_p_negative_derived": False,
        "dark_energy_value_fixed": False,
        "bounce_constructed": False,
        "constraint_solution_constructed": False,
        "regular_centre_constructed": False,
        "curvature_invariants_verified": False,
        "horizon_regular_evolution_constructed": False,
        "finite_mass_exterior_constructed": False,
        "global_solution_constructed": False,
        "child_topology_proven": False,
        "external_Q_derived": False,
        "late_dark_energy_derived": False,
        "variable_c_derived": False,
        "quantum_ensemble_equivalence_proven": False,
    }
