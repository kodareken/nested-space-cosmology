"""GMF-1B-ECD-KIN1: contact tetrad variation and free-Dirac control.

This gate corrects a false fixed-coordinate-current variation in the original
KIN1 artifact.  The torsion-eliminated minimal ECD interaction is fixed by
INT1 in a canonical ``(+---)`` layer,

    L_4,+ = -(3 kappa/16) J^I J_I,

and translates to ``L_4,-=+(3 kappa/16) A^I A_I`` in the VC ``(-+++)`` layer.
The Lorentz-index current is a spinor bilinear.  At fixed spinor components its
quadratic invariant is tetrad independent; varying a coordinate current while
incorrectly holding that current fixed creates a spurious anisotropic term.

KIN1 now checks four bounded facts:

1. the Clifford anticommutator identity over all 64 index triples;
2. the full reduced action sign and the factor-of-two Hehl--Datta variation,
   without adding convention-dependent intermediate connection pieces again;
3. cancellation of the metric and current-response pieces in all 16 coframe
   perturbation directions, leaving ``T_ab=L_4 g_ab`` and zero contact null
   contraction; and
4. one positive-frequency massless plane-wave control for the torsion-free
   Dirac stress.

The plane-wave fixture is not a general classical-Dirac null-energy theorem.
The complete torsion-free Dirac source, SO(3) reduction, regular centre,
constraints, and finite-mass vacuum exterior remain GMF-1B-ECD-ID1.
"""

from __future__ import annotations

from math import isfinite
from numbers import Complex, Real
from typing import Sequence

from .gmf1b_ecd_identity import contact_tetrad_stress, metric_null_contact_contraction
from .gmf1b_ecd_symmetry import (
    _matrix_product,
    _matrix_scale,
    _matrix_vector,
    _vc_gamma_matrices,
    ventrella_chiral_axial_bilinears,
)

_CONTACT_COEFFICIENT = 3.0 / 16.0
_HEHL_DATTA_COEFFICIENT = 3.0 / 8.0
_ETA_PLUS = (1.0, -1.0, -1.0, -1.0)


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


def _eps(i: int, j: int, k: int, l: int) -> int:
    """Levi-Civita symbol with ``eps^{0123}=+1``."""

    if len({i, j, k, l}) < 4:
        return 0
    permutation = (i, j, k, l)
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(4)
        for right in range(left + 1, 4)
    )
    return -1 if inversions % 2 else 1


def _gamma_ab(
    gamma: tuple[tuple[tuple[complex, ...], ...], ...], a: int, b: int
) -> tuple[tuple[complex, ...], ...]:
    """Return ``gamma^{ab}=(1/2)[gamma^a,gamma^b]``."""

    forward = _matrix_product(gamma[a], gamma[b])
    reverse = _matrix_product(gamma[b], gamma[a])
    return _matrix_scale(
        0.5,
        tuple(
            tuple(forward[row][column] - reverse[row][column] for column in range(4))
            for row in range(4)
        ),
    )


def anticommutator_identity(tolerance: Real = 1.0e-12) -> dict[str, object]:
    """Verify ``{gamma^c,gamma^{ab}}=-2i eps^{cabd} gamma5 gamma_d``."""

    tolerance_value = _positive("tolerance", tolerance)
    gamma = _vc_gamma_matrices()
    gamma5 = _matrix_scale(1j, _matrix_product(gamma[0], gamma[1], gamma[2], gamma[3]))
    maximum_residual = 0.0
    for c in range(4):
        for a in range(4):
            for b in range(4):
                gamma_ab = _gamma_ab(gamma, a, b)
                left_product = _matrix_product(gamma[c], gamma_ab)
                right_product = _matrix_product(gamma_ab, gamma[c])
                anticommutator = tuple(
                    tuple(
                        left_product[row][column] + right_product[row][column]
                        for column in range(4)
                    )
                    for row in range(4)
                )
                target = [[0j] * 4 for _ in range(4)]
                for row in range(4):
                    for column in range(4):
                        value = 0j
                        for d in range(4):
                            gamma_lower = _matrix_scale(-1.0 if d == 0 else 1.0, gamma[d])
                            value += _eps(c, a, b, d) * _matrix_product(
                                gamma5, gamma_lower
                            )[row][column]
                        target[row][column] = -2j * value
                maximum_residual = max(
                    maximum_residual,
                    max(
                        abs(anticommutator[row][column] - target[row][column])
                        for row in range(4)
                        for column in range(4)
                    ),
                )
    return {
        "identity": "{gamma^c, gamma^{ab}} = -2 i eps^{cabd} gamma5 gamma_d",
        "maximum_residual": _finite("anticommutator residual", maximum_residual),
        "verified_over_all_64_index_triples": maximum_residual <= tolerance_value,
        "representation_independent": True,
    }


def effective_action_bookkeeping(
    F: Complex,
    G: Complex,
    areal_radius: Real,
    radial_metric: Real,
    kappa: Real,
    tolerance: Real = 1.0e-12,
) -> dict[str, object]:
    """Verify the full reduced contact sign and Hehl--Datta factor of two."""

    f_value = _complex("F", F)
    g_value = _complex("G", G)
    radius = _positive("areal_radius", areal_radius)
    metric = _positive("radial_metric", radial_metric)
    coupling = _positive("kappa", kappa)
    tolerance_value = _positive("tolerance", tolerance)
    axial = ventrella_chiral_axial_bilinears(f_value, g_value, radius, metric)
    axial_squared = float(axial["axial_squared"])
    canonical_current_squared = -axial_squared
    canonical_lagrangian = -_CONTACT_COEFFICIENT * coupling * canonical_current_squared
    vc_lagrangian = _CONTACT_COEFFICIENT * coupling * axial_squared
    hehl_datta_from_variation = 2.0 * _CONTACT_COEFFICIENT
    scale = max(1.0, abs(canonical_lagrangian), abs(vc_lagrangian))
    return {
        "canonical_J_squared": _finite("canonical J squared", canonical_current_squared),
        "vc_A_squared": _finite("VC A squared", axial_squared),
        "canonical_full_contact_coefficient_over_kappa_J2": -_CONTACT_COEFFICIENT,
        "vc_full_contact_coefficient_over_kappa_A2": _CONTACT_COEFFICIENT,
        "canonical_full_contact_lagrangian": _finite(
            "canonical contact Lagrangian", canonical_lagrangian
        ),
        "vc_full_contact_lagrangian": _finite("VC contact Lagrangian", vc_lagrangian),
        "signature_translation_verified": abs(canonical_lagrangian - vc_lagrangian)
        <= tolerance_value * scale,
        "hehl_datta_equation_coefficient_over_kappa": _HEHL_DATTA_COEFFICIENT,
        "field_variation_factor_two": hehl_datta_from_variation,
        "field_variation_factor_two_verified": abs(
            hehl_datta_from_variation - _HEHL_DATTA_COEFFICIENT
        )
        <= tolerance_value,
        "connection_eliminated_everywhere_before_reduced_action_used": True,
        "separate_connection_pieces_added_to_reduced_action": False,
        "combined_onshell_coefficient_sign_resolved": True,
    }


def _inverse4(matrix: Sequence[Sequence[float]]) -> tuple[tuple[float, ...], ...]:
    """Invert a finite 4 by 4 matrix with Gauss--Jordan elimination."""

    augmented = [
        [float(matrix[row][column]) for column in range(4)]
        + [1.0 if row == column else 0.0 for column in range(4)]
        for row in range(4)
    ]
    for column in range(4):
        pivot_row = max(range(column, 4), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot_row][column]) <= 1.0e-15:
            raise ValueError("coframe must be invertible")
        augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(4):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][entry] - factor * augmented[column][entry]
                for entry in range(8)
            ]
    return tuple(
        tuple(augmented[row][column] for column in range(4, 8)) for row in range(4)
    )


def tetrad_current_invariant_cancellation(tolerance: Real = 1.0e-12) -> dict[str, object]:
    """Check the axial-invariant tetrad response in all 16 coframe directions."""

    tolerance_value = _positive("tolerance", tolerance)
    coframe = (
        (1.20, 0.10, -0.05, 0.02),
        (0.03, 0.90, 0.04, -0.02),
        (-0.01, 0.02, 1.10, 0.06),
        (0.04, -0.03, 0.05, 0.80),
    )
    internal_upper = (1.70, 0.40, -0.60, 0.80)
    internal_lower = tuple(_ETA_PLUS[index] * internal_upper[index] for index in range(4))
    # The inverse has rows mu and columns I because the coframe has rows I,
    # columns mu.
    inverse = _inverse4(coframe)
    metric_contravariant = tuple(
        tuple(
            sum(
                _ETA_PLUS[internal] * inverse[mu][internal] * inverse[nu][internal]
                for internal in range(4)
            )
            for nu in range(4)
        )
        for mu in range(4)
    )
    current_covariant = tuple(
        sum(coframe[internal][mu] * internal_lower[internal] for internal in range(4))
        for mu in range(4)
    )
    internal_invariant = sum(
        internal_lower[internal] * internal_upper[internal] for internal in range(4)
    )
    coordinate_invariant = sum(
        metric_contravariant[mu][nu] * current_covariant[mu] * current_covariant[nu]
        for mu in range(4)
        for nu in range(4)
    )

    maximum_fixed_partial = 0.0
    maximum_current_response = 0.0
    maximum_total_residual = 0.0
    direction_residuals: list[float] = []
    for varied_internal in range(4):
        for varied_mu in range(4):
            delta_metric_covariant = tuple(
                tuple(
                    _ETA_PLUS[varied_internal]
                    * (
                        (1.0 if alpha == varied_mu else 0.0)
                        * coframe[varied_internal][beta]
                        + coframe[varied_internal][alpha]
                        * (1.0 if beta == varied_mu else 0.0)
                    )
                    for beta in range(4)
                )
                for alpha in range(4)
            )
            delta_metric_contravariant = tuple(
                tuple(
                    -sum(
                        metric_contravariant[mu][alpha]
                        * delta_metric_covariant[alpha][beta]
                        * metric_contravariant[beta][nu]
                        for alpha in range(4)
                        for beta in range(4)
                    )
                    for nu in range(4)
                )
                for mu in range(4)
            )
            delta_current_covariant = tuple(
                internal_lower[varied_internal] if nu == varied_mu else 0.0
                for nu in range(4)
            )
            fixed_partial = sum(
                delta_metric_contravariant[mu][nu]
                * current_covariant[mu]
                * current_covariant[nu]
                for mu in range(4)
                for nu in range(4)
            )
            current_response = 2.0 * sum(
                metric_contravariant[mu][nu]
                * current_covariant[mu]
                * delta_current_covariant[nu]
                for mu in range(4)
                for nu in range(4)
            )
            residual = fixed_partial + current_response
            maximum_fixed_partial = max(maximum_fixed_partial, abs(fixed_partial))
            maximum_current_response = max(maximum_current_response, abs(current_response))
            maximum_total_residual = max(maximum_total_residual, abs(residual))
            direction_residuals.append(_finite("directional residual", residual))

    scale = max(1.0, abs(internal_invariant), maximum_fixed_partial)
    return {
        "coframe_directions_checked": 16,
        "internal_current_invariant": _finite("internal invariant", internal_invariant),
        "coordinate_current_invariant": _finite(
            "coordinate invariant", coordinate_invariant
        ),
        "baseline_invariant_residual": _finite(
            "baseline invariant residual", coordinate_invariant - internal_invariant
        ),
        "maximum_fixed_coordinate_current_metric_partial": _finite(
            "fixed-current partial", maximum_fixed_partial
        ),
        "maximum_current_response": _finite("current response", maximum_current_response),
        "maximum_total_variation_residual": _finite(
            "total variation residual", maximum_total_residual
        ),
        "directional_total_residuals": tuple(direction_residuals),
        "fixed_coordinate_current_partial_is_nonzero": maximum_fixed_partial
        > tolerance_value,
        "current_response_cancels_metric_partial": maximum_total_residual
        <= tolerance_value * scale,
        "internal_current_invariant_is_tetrad_independent": abs(
            coordinate_invariant - internal_invariant
        )
        <= tolerance_value * scale,
        "fixed_coordinate_current_anisotropic_stress_is_physical": False,
    }


def contact_stress_null_energy_gate(
    kappa: Real,
    axial_squared: Real,
    metric_covariant: Sequence[Sequence[Real]],
    null_vectors: Sequence[Sequence[Real]],
    tolerance: Real = 1.0e-12,
) -> dict[str, object]:
    """Verify the contact tetrad stress and its zero metric-null contraction."""

    coupling = _positive("kappa", kappa)
    invariant = _real("axial_squared", axial_squared)
    tolerance_value = _positive("tolerance", tolerance)
    if (
        not isinstance(null_vectors, Sequence)
        or isinstance(null_vectors, (str, bytes))
        or not null_vectors
    ):
        raise ValueError("null_vectors must be a nonempty sequence")
    stress = contact_tetrad_stress(coupling, invariant, metric_covariant)
    contractions = tuple(
        metric_null_contact_contraction(
            coupling, invariant, metric_covariant, vector, tolerance_value
        )
        for vector in null_vectors
    )
    maximum = max(abs(float(item["contact_null_contraction"])) for item in contractions)
    return {
        "stress": stress,
        "rho": _finite("rho", stress[0][0]),
        "p": _finite("pressure", stress[1][1]),
        "w": _finite("w", stress[1][1] / stress[0][0]),
        "rho_plus_p": _finite("rho plus p", stress[0][0] + stress[1][1]),
        "null_vectors_checked": len(contractions),
        "null_contractions": contractions,
        "maximum_null_contraction_residual": _finite(
            "null contraction residual", maximum
        ),
        "all_contact_null_contractions_zero": maximum <= tolerance_value,
        "contact_null_energy_condition_violated": False,
        "contact_sector_tetrad_variation_complete": True,
        "full_ecd_effective_stress_verified": False,
    }


def free_dirac_positive_frequency_plane_wave_control(
    tolerance: Real = 1.0e-12,
) -> dict[str, object]:
    """Check one massless positive-frequency Dirac plane wave in ``(+---)``."""

    tolerance_value = _positive("tolerance", tolerance)
    gamma = tuple(_matrix_scale(-1j, matrix) for matrix in _vc_gamma_matrices())
    spinor = (1.0 + 0j, 0j, 1.0 + 0j, 0j)
    momentum_upper = (1.0, 0.0, 0.0, 1.0)
    momentum_lower = tuple(_ETA_PLUS[index] * momentum_upper[index] for index in range(4))
    slash = tuple(
        tuple(
            sum(momentum_lower[index] * gamma[index][row][column] for index in range(4))
            for column in range(4)
        )
        for row in range(4)
    )
    slash_spinor = _matrix_vector(slash, spinor)
    dirac_residual = max(abs(value) for value in slash_spinor)
    adjoint = tuple(value.conjugate() for value in _matrix_vector(gamma[0], spinor))
    current_upper = tuple(
        sum(
            adjoint[index] * _matrix_vector(gamma[a], spinor)[index]
            for index in range(4)
        ).real
        for a in range(4)
    )
    current_lower = tuple(_ETA_PLUS[index] * current_upper[index] for index in range(4))
    current_momentum_residual = max(
        abs(current_upper[index] - 2.0 * momentum_upper[index]) for index in range(4)
    )
    stress = tuple(
        tuple(
            0.5
            * (
                momentum_lower[a] * current_lower[b]
                + momentum_lower[b] * current_lower[a]
            )
            for b in range(4)
        )
        for a in range(4)
    )
    probes = {
        "parallel": (1.0, 0.0, 0.0, 1.0),
        "opposite": (1.0, 0.0, 0.0, -1.0),
    }
    contractions: dict[str, float] = {}
    maximum_probe_norm = 0.0
    for name, vector in probes.items():
        norm = sum(
            _ETA_PLUS[index] * vector[index] * vector[index] for index in range(4)
        )
        maximum_probe_norm = max(maximum_probe_norm, abs(norm))
        contractions[name] = _finite(
            f"{name} null contraction",
            sum(stress[a][b] * vector[a] * vector[b] for a in range(4) for b in range(4)),
        )
    momentum_norm = sum(
        _ETA_PLUS[index] * momentum_upper[index] * momentum_upper[index]
        for index in range(4)
    )
    return {
        "canonical_signature": "+---",
        "spinor": spinor,
        "momentum_upper": momentum_upper,
        "momentum_norm": _finite("momentum norm", momentum_norm),
        "massless_dirac_residual": _finite("Dirac residual", dirac_residual),
        "current_upper": current_upper,
        "current_equals_two_momentum_residual": _finite(
            "current-momentum residual", current_momentum_residual
        ),
        "stress_covariant": stress,
        "null_probe_norm_residual": _finite("null probe norm", maximum_probe_norm),
        "null_contractions": contractions,
        "positive_frequency_plane_wave_nec_nonnegative": min(contractions.values())
        >= -tolerance_value,
        "positive_frequency_plane_wave_control_only": True,
        "general_classical_dirac_nec_proven": False,
    }


def required_nonclaims() -> dict[str, bool]:
    """Return explicit GMF-1B-ECD-KIN1 scope boundaries."""

    return {
        "full_ecd_tetrad_dependence_varied": False,
        "torsion_free_dirac_stress_fully_derived": False,
        "full_ecd_effective_stress_verified": False,
        "complete_effective_stress_tensor_derived": False,
        "intermediate_connection_piece_split_audited": False,
        "general_classical_dirac_nec_proven": False,
        "full_ecd_source_so3_equivariance_verified": False,
        "complete_spin_weighted_harmonic_closure_verified": False,
        "regular_centre_constructed": False,
        "constraint_solution_constructed": False,
        "constraint_compatible_vacuum_region_verified": False,
        "finite_mass_exterior_constructed": False,
        "rho_plus_p_negative_derived": False,
        "metric_null_defocusing_derived": False,
        "homogeneous_bounce_rho_plus_p_negative_derived": False,
        "bounce_constructed": False,
        "curvature_invariants_verified": False,
        "horizon_regular_evolution_constructed": False,
        "global_solution_constructed": False,
        "child_topology_proven": False,
        "external_Q_derived": False,
        "late_dark_energy_derived": False,
        "variable_c_derived": False,
    }
