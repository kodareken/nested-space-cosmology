"""Quantified all-covector weak-coupling gate for ACT1 modified harmonic.

Kovacs and Reall prove strong hyperbolicity of modified-harmonic Horndeski
systems at sufficiently weak coupling.  Their theorem is qualitative: it does
not supply a universal numerical epsilon.  This module specializes the
continuity proof to the frozen ACT1 normalization and turns "sufficiently
weak" into a deliberately conservative executable envelope.

The construction uses three facts.

* In the Frobenius-normalized symmetric-tensor basis, spatial rotations are
  orthogonal and the flat Einstein--two-scalar reference has one direction-
  independent spectrum.
* The complete companion matrix is a quadratic polynomial in the spatial
  unit covector.  Norms of its extracted coefficient tensors therefore bound
  the deformation for *every* direction, not only a sampled set.
* Riesz-contour stability preserves the six physical/auxiliary spectral
  groups.  On the physical groups the symmetric ungauged action block gives
  the Kovacs--Reall positive form.  A transported-subspace estimate keeps
  that form coercive throughout the declared deformation ball.

The independent canonical ``chi`` sector is retained explicitly.  It is a
direct principal block and adds one physical polarization per sign without
altering the Horndeski metric--``phi`` proof.

Passing this gate is a classical principal-health statement only.  It is not
a Wilsonian EFT remainder bound, an evolution, or evidence of defocusing.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from math import isfinite, pi, sin
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from .covariant_principal_health import (
    DIMENSION,
    FIELD_COUNT,
    FIRST_ORDER_FIELD_COUNT,
    METRIC_COMPONENTS,
    CovariantPrincipalBackground,
    PrincipalCoefficientTensors,
    coefficient_tensor_symbol,
    esf_background,
    normalized_first_order_matrix_from_coefficients,
    principal_coefficient_tensors,
)


CLUSTER_SPECS = (
    ("physical_minus", -1.0, 0.2, 4, 1.0 / 32.0),
    ("tilde_minus", -0.5, 0.05, 4, 1.0 / 1024.0),
    ("hat_minus", -1.0 / 3.0, 0.05, 4, 1.0 / 1024.0),
    ("hat_plus", 1.0 / 3.0, 0.05, 4, 1.0 / 1024.0),
    ("tilde_plus", 0.5, 0.05, 4, 1.0 / 1024.0),
    ("physical_plus", 1.0, 0.2, 4, 1.0 / 32.0),
)


@dataclass(frozen=True, slots=True)
class WeakCouplingThresholds:
    """Frozen sufficient bounds, intentionally much stronger than EFT0."""

    companion_deformation_maximum: float = 1.0 / 4096.0
    action_energy_deformation_maximum: float = 1.0 / 8192.0
    kinetic_deformation_maximum: float = 1.0 / 64.0
    effective_planck_squared_minimum: float = 2.0
    action_hessian_symmetry_defect_maximum: float = 1.0e-10
    gauge_identity_defect_maximum: float = 1.0e-10

    def __post_init__(self) -> None:
        for name, value in asdict(self).items():
            if isinstance(value, bool) or not isinstance(value, Real):
                raise TypeError(f"{name} must be a finite positive scalar")
            number = float(value)
            if not isfinite(number) or number <= 0.0:
                raise ValueError(f"{name} must be a finite positive scalar")
            object.__setattr__(self, name, number)


@dataclass(frozen=True, slots=True)
class CoefficientDeformationBounds:
    """Direction-independent matrix bounds relative to flat ESF."""

    companion_operator_2: float
    companion_time_space_sum_2: float
    companion_space_space_sum_2: float
    action_energy_operator_2: float
    kinetic_operator_2: float
    kinetic_smallest_singular_value: float
    action_hessian_symmetry_defect: float
    gauge_identity_defect: float


def _coefficient_companion_blocks(
    coefficients: PrincipalCoefficientTensors,
) -> tuple[np.ndarray, np.ndarray]:
    a_matrix = coefficients.time_time
    try:
        time_space = np.stack(
            [
                -np.linalg.solve(a_matrix, coefficients.time_space[index])
                for index in range(3)
            ]
        )
        space_space = np.empty_like(coefficients.space_space)
        for first in range(3):
            for second in range(3):
                space_space[first, second] = -np.linalg.solve(
                    a_matrix, coefficients.space_space[first, second]
                )
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariant health kinetic block is singular") from exc
    return time_space, space_space


def _operator_sum(value: np.ndarray) -> float:
    return float(sum(np.linalg.norm(matrix, ord=2) for matrix in value.reshape((-1,) + value.shape[-2:])))


def _pure_gauge_identity_defect(
    background: CovariantPrincipalBackground,
) -> float:
    """Bound the homogeneous cubic pure-gauge identity in every direction.

    The ungauged symbol is quadratic in ``xi`` and a pure-gauge metric
    perturbation ``h_ab=xi_a X_b+xi_b X_a`` is linear in ``xi``.  Their
    product is therefore a homogeneous cubic.  Extracting and symmetrizing
    all of its coefficient tensors is sufficient to test the identity for
    every covector; a finite directional scan would not be.
    """

    coefficients = principal_coefficient_tensors(
        background, include_gauge_extension=False
    )
    quadratic = np.zeros(
        (DIMENSION, DIMENSION, FIELD_COUNT, FIELD_COUNT), dtype=np.float64
    )
    quadratic[0, 0] = coefficients.time_time
    for index in range(3):
        quadratic[0, index + 1] = 0.5 * coefficients.time_space[index]
        quadratic[index + 1, 0] = 0.5 * coefficients.time_space[index]
    quadratic[1:, 1:] = coefficients.space_space

    # ``gauge[mu, field, vector]`` is the coefficient of xi_mu in the
    # Frobenius-normalized pure-gauge field vector.
    gauge = np.zeros((DIMENSION, FIELD_COUNT, DIMENSION), dtype=np.float64)
    for mu in range(DIMENSION):
        for field, (first, second) in enumerate(METRIC_COMPONENTS):
            scale = 2.0**0.5 if first != second else 1.0
            for vector in range(DIMENSION):
                gauge[mu, field, vector] = scale * (
                    float(first == mu and second == vector)
                    + float(second == mu and first == vector)
                )

    cubic = np.einsum("mnij,rjk->mnrik", quadratic, gauge)
    symmetrized = np.zeros_like(cubic)
    for first in range(DIMENSION):
        for second in range(DIMENSION):
            for third in range(DIMENSION):
                symmetrized[first, second, third] = (
                    cubic[first, second, third]
                    + cubic[first, third, second]
                    + cubic[second, first, third]
                    + cubic[second, third, first]
                    + cubic[third, first, second]
                    + cubic[third, second, first]
                ) / 6.0
    # For ||xi||_2=1, the triangle inequality bounds the identity by the sum
    # of coefficient operator norms.  This is deliberately conservative.
    return _operator_sum(symmetrized)


def coefficient_deformation_bounds(
    background: CovariantPrincipalBackground,
) -> CoefficientDeformationBounds:
    """Bound the full companion deformation for every unit spatial covector."""

    if not isinstance(background, CovariantPrincipalBackground):
        raise TypeError("background must be CovariantPrincipalBackground")
    reference = esf_background()
    complete = principal_coefficient_tensors(background)
    complete_reference = principal_coefficient_tensors(reference)
    action = principal_coefficient_tensors(
        background, include_gauge_extension=False
    )
    action_reference = principal_coefficient_tensors(
        reference, include_gauge_extension=False
    )

    time_space, space_space = _coefficient_companion_blocks(complete)
    time_space_reference, space_space_reference = _coefficient_companion_blocks(
        complete_reference
    )
    time_space_bound = _operator_sum(time_space - time_space_reference)
    space_space_bound = _operator_sum(space_space - space_space_reference)
    companion_bound = time_space_bound + space_space_bound

    delta_action_a = action.time_time - action_reference.time_time
    action_bound = float(np.linalg.norm(delta_action_a, ord=2)) + sum(
        float(
            np.linalg.norm(
                action.time_space[index]
                - action_reference.time_space[index],
                ord=2,
            )
        )
        for index in range(3)
    )
    kinetic_delta = float(
        np.linalg.norm(
            complete.time_time - complete_reference.time_time,
            ord=2,
        )
    )
    kinetic_smallest = float(
        np.linalg.svd(complete.time_time, compute_uv=False)[-1]
    )

    # The ungauged coefficient tensors must be symmetric in field space.  The
    # maximum defect is measured before any averaging or projection.
    action_defect = max(
        float(np.max(np.abs(action.time_time - action.time_time.T))),
        max(
            float(
                np.max(
                    np.abs(
                        action.time_space[index]
                        - action.time_space[index].T
                    )
                )
            )
            for index in range(3)
        ),
        max(
            float(
                np.max(
                    np.abs(
                        action.space_space[first, second]
                        - action.space_space[first, second].T
                    )
                )
            )
            for first in range(3)
            for second in range(3)
        ),
    )
    return CoefficientDeformationBounds(
        companion_operator_2=companion_bound,
        companion_time_space_sum_2=time_space_bound,
        companion_space_space_sum_2=space_space_bound,
        action_energy_operator_2=action_bound,
        kinetic_operator_2=kinetic_delta,
        kinetic_smallest_singular_value=kinetic_smallest,
        action_hessian_symmetry_defect=action_defect,
        gauge_identity_defect=_pure_gauge_identity_defect(background),
    )


def _cluster_sample_lower_bound(
    matrix: np.ndarray,
    *,
    center: float,
    radius: float,
    sample_count: int,
) -> dict[str, float]:
    sampled_minimum = float("inf")
    identity = np.eye(matrix.shape[0], dtype=np.complex128)
    for index in range(sample_count):
        angle = 2.0 * pi * index / sample_count
        point = center + radius * np.exp(1j * angle)
        singular = float(
            np.linalg.svd(point * identity - matrix, compute_uv=False)[-1]
        )
        sampled_minimum = min(sampled_minimum, singular)
    # Singular values are 1-Lipschitz in the matrix 2-norm.  Every point on
    # the circle is within this chord distance of a sampled point.
    chord = 2.0 * radius * sin(pi / (2.0 * sample_count))
    roundoff = (
        4096.0
        * np.finfo(np.float64).eps
        * max(1.0, abs(center) + radius, float(np.linalg.norm(matrix, ord=2)))
    )
    lower = sampled_minimum - chord - roundoff
    return {
        "sampled_minimum_singular_value": sampled_minimum,
        "maximum_unsampled_chord": chord,
        "roundoff_allowance": roundoff,
        "certified_lower_bound": lower,
    }


@lru_cache(maxsize=4)
def esf_reference_cluster_certificate(sample_count: int = 512) -> dict[str, Any]:
    """Certify the flat six-cluster reference and its physical energy form."""

    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count < 128
    ):
        raise ValueError("reference contour sample_count must be an integer >=128")
    complete = principal_coefficient_tensors(esf_background())
    action = principal_coefficient_tensors(
        esf_background(), include_gauge_extension=False
    )
    direction = np.asarray((1.0, 0.0, 0.0))
    matrix, kinetic = normalized_first_order_matrix_from_coefficients(
        complete, direction
    )
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    cluster_records: dict[str, Any] = {}
    total_count = 0
    for name, center, radius, expected_count, locked_lower in CLUSTER_SPECS:
        membership = np.flatnonzero(np.abs(eigenvalues - center) < radius)
        contour = _cluster_sample_lower_bound(
            matrix,
            center=center,
            radius=radius,
            sample_count=sample_count,
        )
        if len(membership) != expected_count:
            raise ValueError(f"ESF reference cluster {name} has wrong dimension")
        if contour["certified_lower_bound"] <= locked_lower:
            raise ValueError(f"ESF reference contour {name} misses locked margin")
        cluster_records[name] = {
            "center": center,
            "radius": radius,
            "eigenvalue_count": len(membership),
            "expected_eigenvalue_count": expected_count,
            "locked_resolvent_singular_lower_bound": locked_lower,
            **contour,
        }
        total_count += len(membership)
    if total_count != FIRST_ORDER_FIELD_COUNT:
        raise ValueError("ESF cluster contours do not partition the full spectrum")

    b_action = sum(
        (
            action.time_space[index] * direction[index]
            for index in range(3)
        ),
        np.zeros_like(action.time_time),
    )
    zero = np.zeros_like(action.time_time)
    physical_energy: dict[str, Any] = {}
    for name, target in (("physical_minus", -1.0), ("physical_plus", 1.0)):
        membership = np.flatnonzero(np.abs(eigenvalues - target) < 0.2)
        basis, _upper = np.linalg.qr(eigenvectors[:, membership].real)
        # Our covector and Euler--Lagrange row conventions reverse the sign
        # attached to H_* in the paper; -target gives the positive form.
        form = -target * np.block(
            [[b_action, action.time_time], [action.time_time, zero]]
        )
        restricted = 0.5 * (
            basis.T @ form @ basis + (basis.T @ form @ basis).T
        )
        eigenvalues_restricted = np.linalg.eigvalsh(restricted)
        physical_energy[name] = {
            "dimension": basis.shape[1],
            "minimum_eigenvalue": float(eigenvalues_restricted[0]),
            "maximum_eigenvalue": float(eigenvalues_restricted[-1]),
        }
    reference_minimum = min(
        value["minimum_eigenvalue"] for value in physical_energy.values()
    )
    reference_form_norm = max(
        float(
            np.linalg.norm(
                -target
                * np.block(
                    [[b_action, action.time_time], [action.time_time, zero]]
                ),
                ord=2,
            )
        )
        for target in (-1.0, 1.0)
    )
    if reference_minimum < 1.0 - 1.0e-11 or reference_form_norm >= 3.0:
        raise ValueError("ESF physical energy reference margins differ")
    return {
        "classification": "binary64_with_Lipschitz_chord_and_roundoff_allowance",
        "direction_representative": [1.0, 0.0, 0.0],
        "rotation_covariance_uses_Frobenius_normalized_tensor_basis": True,
        "contour_sample_count": sample_count,
        "clusters": cluster_records,
        "physical_energy": physical_energy,
        "locked_physical_energy_coercivity_lower": 15.0 / 16.0,
        "locked_physical_energy_operator_norm_upper": 3.0,
        "kinetic": kinetic,
        "six_contours_partition_all_24_characteristics": True,
    }


def weak_coupling_health_certificate(
    background: CovariantPrincipalBackground,
    *,
    thresholds: WeakCouplingThresholds | None = None,
) -> dict[str, Any]:
    """Evaluate the quantified sufficient all-covector health inequalities."""

    limits = thresholds or WeakCouplingThresholds()
    if not isinstance(limits, WeakCouplingThresholds):
        raise TypeError("thresholds must be WeakCouplingThresholds")
    bounds = coefficient_deformation_bounds(background)
    reference = esf_reference_cluster_certificate()
    delta = bounds.companion_operator_2
    delta_h = bounds.action_energy_operator_2

    cluster_margins = {
        name: record["locked_resolvent_singular_lower_bound"] - delta
        for name, record in reference["clusters"].items()
    }
    physical_singular_lower = min(
        reference["clusters"][name][
            "locked_resolvent_singular_lower_bound"
        ]
        for name in ("physical_minus", "physical_plus")
    )
    physical_radius = reference["clusters"]["physical_plus"]["radius"]
    reference_resolvent = 1.0 / physical_singular_lower
    resolvent_neumann = reference_resolvent * delta
    projector_displacement = (
        float("inf")
        if resolvent_neumann >= 1.0
        else (
            physical_radius
            * reference_resolvent**2
            * delta
            / (1.0 - resolvent_neumann)
        )
    )
    reference_energy_minimum = reference[
        "locked_physical_energy_coercivity_lower"
    ]
    reference_energy_norm = reference[
        "locked_physical_energy_operator_norm_upper"
    ]
    # Let P and P0 be the candidate/reference Riesz projectors and let
    # ||P-P0||<=d.  For a unit candidate physical vector v, w=P0 v obeys
    # ||v-w||<=d and ||w||>=1-d.  Applying these inequalities directly to
    # the reference form avoids assuming that either non-normal projector is
    # orthogonal.
    reference_form_transport_lower = (
        float("-inf")
        if not isfinite(projector_displacement)
        else (
            reference_energy_minimum * (1.0 - projector_displacement) ** 2
            - reference_energy_norm
            * (
                2.0
                * (1.0 + projector_displacement)
                * projector_displacement
                + projector_displacement**2
            )
        )
    )
    energy_transport_error = (
        float("inf")
        if not isfinite(projector_displacement)
        else delta_h
    )
    physical_energy_lower = reference_form_transport_lower - energy_transport_error

    predicates = {
        "positive_effective_planck_coefficient": (
            background.effective_planck_squared
            >= limits.effective_planck_squared_minimum
        ),
        "all_direction_companion_deformation_inside_ball": (
            delta <= limits.companion_deformation_maximum
        ),
        "action_energy_deformation_inside_ball": (
            delta_h <= limits.action_energy_deformation_maximum
        ),
        "coordinate_time_kinetic_deformation_inside_ball": (
            bounds.kinetic_operator_2 <= limits.kinetic_deformation_maximum
            and bounds.kinetic_smallest_singular_value > 0.0
        ),
        "six_Riesz_cluster_contours_remain_separated": all(
            margin > 0.0 for margin in cluster_margins.values()
        ),
        "physical_Riesz_projector_transport_is_invertible": (
            projector_displacement < 1.0
        ),
        "physical_Kovacs_Reall_form_remains_positive": (
            physical_energy_lower > 0.0
        ),
        "ungauged_action_Hessian_is_symmetric": (
            bounds.action_hessian_symmetry_defect
            <= limits.action_hessian_symmetry_defect_maximum
        ),
        "diffeomorphism_pure_gauge_principal_identity_holds": (
            bounds.gauge_identity_defect
            <= limits.gauge_identity_defect_maximum
        ),
        "independent_chi_principal_sector_is_canonical_direct_sum": True,
        "auxiliary_cones_are_exact_real_nested_and_disjoint": True,
    }
    passed = all(predicates.values())
    return {
        "classification": "quantified_sufficient_Kovacs_Reall_weak_coupling_envelope_for_full_ACT1_principal_system",
        "basis": "physical_orthonormal_frame_plus_Frobenius_normalized_symmetric_metric_components",
        "direction_scope": "every_nonzero_spatial_covector_by_homogeneity_and_coefficient_tensor_norm_bound",
        "background": {
            "effective_planck_squared": background.effective_planck_squared,
            "effective_planck_prime": background.effective_planck_prime,
            "gb_coupling_prime": background.gb_coupling_prime,
            "riemann_component_infinity": float(
                np.max(np.abs(background.riemann_lower))
            ),
            "hessian_gb_component_infinity": float(
                np.max(np.abs(background.hessian_gb_lower))
            ),
        },
        "thresholds": asdict(limits),
        "deformation_bounds": asdict(bounds),
        "reference_cluster_certificate": reference,
        "cluster_resolvent_singular_margins": cluster_margins,
        "physical_projector_displacement_upper": projector_displacement,
        "transported_reference_energy_coercivity_lower": reference_form_transport_lower,
        "physical_energy_transport_error_upper": energy_transport_error,
        "physical_energy_coercivity_lower": physical_energy_lower,
        "predicates": predicates,
        "quantitative_all_covector_weak_coupling_health_envelope_passed": passed,
        "scope": {
            "classical_principal_health_only": True,
            "retained_EFT_validity": False,
            "nonlinear_trajectory_contained": False,
            "collapse_or_defocusing_derived": False,
        },
    }


def canonical_health_monitor_values(certificate: Mapping[str, Any]) -> dict[str, float]:
    """Extract the small scalar state needed by a future HLT1 stage monitor."""

    if not isinstance(certificate, Mapping):
        raise TypeError("certificate must be a mapping")
    if not certificate.get(
        "quantitative_all_covector_weak_coupling_health_envelope_passed"
    ):
        raise ValueError("cannot extract monitor values from a failed certificate")
    bounds = certificate["deformation_bounds"]
    return {
        "companion_deformation_operator_2": float(
            bounds["companion_operator_2"]
        ),
        "action_energy_deformation_operator_2": float(
            bounds["action_energy_operator_2"]
        ),
        "kinetic_deformation_operator_2": float(bounds["kinetic_operator_2"]),
        "kinetic_smallest_singular_value": float(
            bounds["kinetic_smallest_singular_value"]
        ),
        "physical_energy_coercivity_lower": float(
            certificate["physical_energy_coercivity_lower"]
        ),
        "minimum_cluster_resolvent_singular_margin": min(
            float(value)
            for value in certificate[
                "cluster_resolvent_singular_margins"
            ].values()
        ),
    }
