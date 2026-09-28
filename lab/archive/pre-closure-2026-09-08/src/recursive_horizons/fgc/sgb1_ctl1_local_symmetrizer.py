"""Point-local / nonflat SGB-L cone and symmetrizer certificate.

This owner does not change the ESF-perturbation continuum argument in
:mod:`sgb1_ctl1_cone`.  On an exact singleton, or a caller-declared small
interval, :class:`SGBLPrincipalBackgroundBox`, binary64 spectra seed
candidates only.  Real eigenvalue clusters, eigenvectors, and Riesz
projectors are certified by exact rational kernels and
``V Λ V^{-1} = M``, by rational interval Krawczyk/Neumann on a reduced
pencil, or not at all.  A positive symmetrizer is
``H = (V^{-1})^T V^{-1}`` with Sylvester minors, or a coefficientwise
construction.  All-direction is claimed only from a coefficientwise
verified construction; otherwise the all-direction slot is typed
incomplete.  Complete 24-mode multiplicity, basis, and coercivity are
retained when every exact margin is strict.  Aggregate
``SGBL_branch_owned_and_healthy`` stays false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np

from .evolution.covariant_principal_health import (
    FIELD_COUNT,
    FIRST_ORDER_FIELD_COUNT,
)
from .exact_interval import Interval
from .exact_interval_krawczyk import parametric_krawczyk_inclusion
from .exact_linear_algebra import (
    inverse,
    matrix_multiply,
    rank,
    transpose,
)
from .sgb1_ctl1_cone import (
    ORTHONORMAL_MINKOWSKI_CHART,
    SGBLCoefficientEnclosure,
    SGBLConeStop,
    SGBLPrincipalBackgroundBox,
    _as_rational_matrix,
    _block_rational_energy,
    _columns_to_matrix,
    _EvaluationCounter,
    _exact_kernel,
    _exact_radial_companion,
    _extract_coefficient_tensors,
    _gershgorin_lower,
    _interval_inverse,
    _kinetic_singular_value_lower,
    _leading_minors,
    _restrict_form,
    _sqrt2_enclosure,
    sgbl_enclose_principal_coefficients,
    sgbl_esf_principal_box,
    sgbl_nonflat_source_principal_box,
    sgbl_pure_gauge_identity_defect,
    sgbl_structural_identities,
)


Q = Fraction
LOCAL_DIRECTION = (1, 0, 0)
ESF_SPEEDS = (Q(-1), Q(-1, 2), Q(-1, 3), Q(1, 3), Q(1, 2), Q(1))
ESF_CLUSTER_CENTERS = (
    ("physical_minus", Q(-1)),
    ("tilde_minus", Q(-1, 2)),
    ("hat_minus", Q(-1, 3)),
    ("hat_plus", Q(1, 3)),
    ("tilde_plus", Q(1, 2)),
    ("physical_plus", Q(1)),
)
DEFAULT_MAX_SEED_DENOMINATOR = 12
DEFAULT_MAX_SEED_CANDIDATES = 48
DEFAULT_MAX_KRAWCZYK_BITS = 64
DEFAULT_KRAWCZYK_DISPLACEMENT = Q(1, 1 << 16)
PHI_FIELD_INDEX = 10
CHI_FIELD_INDEX = 11
FORBIDDEN_HEALTH_IMPORTS = (
    "WeakCouplingThresholds",
    "weak_coupling_health_certificate",
    "esf_reference_cluster_certificate",
    "canonical_health_monitor_values",
    "background_from_spherical_state",
    "FGCQRActionParameters",
    "sgbl_enclose_principal_cone",
)
LOCAL_THEOREM = (
    "On an exact singleton orthonormal SGB-L principal box, floating "
    "spectra seed candidates only. Real eigenvalue clusters, eigenvectors, "
    "and Riesz projectors are certified by exact rational kernels together "
    "with V Λ V^{-1} = M, or by rational interval Krawczyk/Neumann on a "
    "reduced pencil. A positive symmetrizer is H = (V^{-1})^T V^{-1} with "
    "Sylvester minors, or a coefficientwise construction. All-direction is "
    "claimed only from a coefficientwise verified construction; otherwise "
    "that slot is typed incomplete. Complete 24-mode multiplicity, basis, "
    "and coercivity are required of a proved local cone. Aggregate health "
    "stays false."
)
FAMILY_BOX_API = MappingProxyType(
    {
        "entrypoint": "sgbl_local_symmetrizer_from_family_box",
        "accepted_slot": "SGBLPrincipalBackgroundBox",
        "later_owner": "uniform_family_cell_product_box_local_symmetrizer",
        "family_covering_proved": False,
        "aggregate_health": False,
        "note": (
            "Family boxes are accepted as declared local boxes. Uniform "
            "covering of a family cell product box is a later owner, not a "
            "hidden zero or health pass."
        ),
    }
)
PROVED_CLASSIFICATION = "point_local_cone_proved"
INCONCLUSIVE_CLASSIFICATION = "interval_inconclusive"
INCOMPLETE_CLASSIFICATION = "typed_incomplete"
ALL_DIRECTION_COEFFICIENTWISE = "coefficientwise_positive_symmetrizer"
ALL_DIRECTION_INCOMPLETE = "typed_incomplete"
NONFLAT_OBSTRUCTION = "complete_24_mode_basis_not_certified"
STOP_REASONS = frozenset(
    {
        "cone_chart_error",
        "sign_chart_error",
        "resource_limit",
        "lost_kinetic",
        "broken_gauge",
        "defective_cluster",
        "complex_spectrum",
        "cluster_coalescence",
        "eigenframe_incomplete",
        "coercivity_failed",
        "symmetry_failure",
        "tampered_record",
    }
)


def _fraction(name: str, value: object, *, nonnegative: bool = False) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    result = Fraction(value)
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _zero_square(size: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(Q(0) for _ in range(size)) for _ in range(size))


def _identity_square(size: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Q(int(row == column)) for column in range(size)) for row in range(size)
    )


def _matrix_add(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(len(left[0])))
        for row in range(len(left))
    )


def _matrix_scale(
    matrix: Sequence[Sequence[Fraction]],
    factor: Fraction,
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(factor * entry for entry in row) for row in matrix)


def _matrix_vector(
    matrix: Sequence[Sequence[Fraction]],
    vector: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(len(vector)))
        for row in range(len(matrix))
    )


def _infinity_norm(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    return max(sum(abs(entry) for entry in row) for row in matrix)


def _matrix_bits(matrix: Sequence[Sequence[Fraction]]) -> int:
    return max(
        max(abs(entry.numerator).bit_length(), entry.denominator.bit_length())
        for row in matrix
        for entry in row
    )


def _sha256_rationals(values: Sequence[object]) -> str:
    parts = []
    for value in values:
        if isinstance(value, Fraction):
            parts.append(f"{value.numerator}/{value.denominator}")
        else:
            parts.append(str(value))
    return sha256(",".join(parts).encode()).hexdigest()


class SGBLLocalSymmetrizerStop(ArithmeticError):
    """Typed local-symmetrizer stop. Wrapping records are not this type."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in STOP_REASONS:
            raise ValueError(f"unknown local-symmetrizer stop reason {reason}")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLLocalSymmetrizerLimits:
    """Declared local arithmetic caps. Not fitted spectral gaps."""

    max_seed_denominator: int = DEFAULT_MAX_SEED_DENOMINATOR
    max_seed_candidates: int = DEFAULT_MAX_SEED_CANDIDATES
    max_krawczyk_matrix_bits: int = DEFAULT_MAX_KRAWCZYK_BITS
    krawczyk_displacement: Fraction | int = DEFAULT_KRAWCZYK_DISPLACEMENT

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "max_seed_denominator",
            _positive_int("max_seed_denominator", self.max_seed_denominator),
        )
        object.__setattr__(
            self,
            "max_seed_candidates",
            _positive_int("max_seed_candidates", self.max_seed_candidates),
        )
        object.__setattr__(
            self,
            "max_krawczyk_matrix_bits",
            _positive_int("max_krawczyk_matrix_bits", self.max_krawczyk_matrix_bits),
        )
        object.__setattr__(
            self,
            "krawczyk_displacement",
            _fraction(
                "krawczyk_displacement",
                self.krawczyk_displacement,
                nonnegative=True,
            ),
        )
        if self.krawczyk_displacement <= 0:
            raise ValueError("krawczyk_displacement must be positive")


def _wrap_cone_stop(exc: SGBLConeStop) -> SGBLLocalSymmetrizerStop:
    reason = exc.reason if exc.reason in STOP_REASONS else "resource_limit"
    return SGBLLocalSymmetrizerStop(reason, str(exc), exc.payload)


def _standard_basis(size: int, index: int) -> tuple[Fraction, ...]:
    return tuple(Q(int(slot == index)) for slot in range(size))


def _companion_from_seed_float(
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[np.ndarray, np.ndarray]:
    """Binary64 spectra used only to seed rational candidates. Not a proof."""

    array = np.array(
        [[float(entry) for entry in row] for row in matrix],
        dtype=np.float64,
    )
    values, vectors = np.linalg.eig(array)
    return values, vectors


def _rational_candidates_from_seeds(
    seeds: Sequence[complex],
    *,
    spectral_radius: Fraction,
    limits: SGBLLocalSymmetrizerLimits,
) -> tuple[Fraction, ...]:
    candidates = set(ESF_SPEEDS)
    for seed in seeds:
        if abs(seed.imag) > 1.0e-8:
            continue
        real = float(seed.real)
        if abs(real) > float(spectral_radius) + 1.0:
            continue
        for denominator in range(1, limits.max_seed_denominator + 1):
            numerator = round(real * denominator)
            candidate = Q(numerator, denominator)
            if abs(float(candidate) - real) <= 1.0e-8:
                candidates.add(candidate)
    ordered = tuple(sorted(candidates, key=lambda value: (value.denominator, value)))
    if len(ordered) > limits.max_seed_candidates:
        raise SGBLLocalSymmetrizerStop(
            "resource_limit",
            "seeded rational eigenvalue candidates exceeded the declared cap",
            {"count": len(ordered), "limit": limits.max_seed_candidates},
        )
    return ordered


def _cluster_name(center: Fraction) -> str:
    for name, speed in ESF_CLUSTER_CENTERS:
        if speed == center:
            return name
    sign = "plus" if center > 0 else "minus" if center < 0 else "zero"
    return f"nonrational_{sign}_{center.numerator}_{center.denominator}"


def _riesz_projectors(
    matrix: Sequence[Sequence[Fraction]],
    frame: Sequence[Sequence[Fraction]],
    frame_inverse: Sequence[Sequence[Fraction]],
    clusters: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Exact Riesz projectors ``P = V I_J V^{-1}`` for certified clusters."""

    size = len(matrix)
    projectors: dict[str, Any] = {}
    column_start = 0
    for name, cluster in clusters.items():
        width = int(cluster["exact_nullity"])
        selector_rows = [[Q(0) for _ in range(size)] for _ in range(size)]
        for offset in range(width):
            index = column_start + offset
            selector_rows[index][index] = Q(1)
        selector = tuple(tuple(row) for row in selector_rows)
        projector = matrix_multiply(
            matrix_multiply(frame, selector),
            frame_inverse,
        )
        squared = matrix_multiply(projector, projector)
        if squared != projector:
            raise SGBLLocalSymmetrizerStop(
                "defective_cluster",
                f"Riesz projector for {name} is not idempotent",
                {"cluster": name},
            )
        trace = sum(projector[index][index] for index in range(size))
        if trace != width:
            raise SGBLLocalSymmetrizerStop(
                "defective_cluster",
                f"Riesz projector trace for {name} is {trace}, not {width}",
                {"cluster": name, "trace": trace},
            )
        commute = matrix_multiply(matrix, projector)
        other = matrix_multiply(projector, matrix)
        if commute != other:
            raise SGBLLocalSymmetrizerStop(
                "defective_cluster",
                f"Riesz projector for {name} does not commute with M",
                {"cluster": name},
            )
        scaled = _matrix_scale(projector, cluster["center"])
        if commute != scaled:
            raise SGBLLocalSymmetrizerStop(
                "defective_cluster",
                f"Riesz projector for {name} is not a pure eigen-projector",
                {"cluster": name},
            )
        projectors[name] = {
            "trace": trace,
            "idempotent": True,
            "commutes_with_companion": True,
            "pure_eigenprojector": True,
            "method": "exact_V_I_J_V_inverse",
        }
        column_start += width
    return projectors


def _positive_symmetrizer(
    frame: Sequence[Sequence[Fraction]],
    frame_inverse: Sequence[Sequence[Fraction]],
    matrix: Sequence[Sequence[Fraction]],
) -> dict[str, Any]:
    symmetrizer = matrix_multiply(transpose(frame_inverse), frame_inverse)
    minors = _leading_minors(symmetrizer)
    if any(minor <= 0 for minor in minors):
        raise SGBLLocalSymmetrizerStop(
            "coercivity_failed",
            "local symmetrizer fails Sylvester's criterion",
            {"leading_minors": minors},
        )
    product = matrix_multiply(symmetrizer, matrix)
    skew = max(
        abs(product[row][column] - product[column][row])
        for row in range(len(product))
        for column in range(len(product))
    )
    if skew != 0:
        raise SGBLLocalSymmetrizerStop(
            "symmetry_failure",
            "H M is not exactly symmetric",
            {"skew_upper": skew},
        )
    return {
        "sylvester_leading_minors_positive": True,
        "leading_minor_count": len(minors),
        "first_minor": minors[0],
        "last_minor": minors[-1],
        "hm_symmetry_defect": 0,
        "formula": "H=(V^{-1})^T V^{-1}",
    }


def _physical_energy(
    action_b: Sequence[Sequence[Fraction]],
    action_a: Sequence[Sequence[Fraction]],
    columns_by_center: Mapping[Fraction, Sequence[Sequence[Fraction]]],
    *,
    mutate_sign: int = 1,
) -> dict[str, Any]:
    if mutate_sign not in (-1, 1):
        raise ValueError("energy mutation sign must be +1 or -1")
    energy: dict[str, Any] = {}
    lower = None
    for name, target in (
        ("physical_minus", Q(-1)),
        ("physical_plus", Q(1)),
    ):
        vectors = tuple(columns_by_center.get(target, ()))
        if len(vectors) != 4:
            raise SGBLLocalSymmetrizerStop(
                "coercivity_failed",
                f"physical energy at {name} requires a four-dimensional eigenspace",
                {"count": len(vectors), "cluster": name},
            )
        form = _block_rational_energy(
            action_b,
            action_a,
            sign=Q(mutate_sign) * (Q(1) if target < 0 else Q(-1)),
        )
        restricted = _restrict_form(form, vectors)
        minors = _leading_minors(restricted)
        if any(minor <= 0 for minor in minors):
            raise SGBLLocalSymmetrizerStop(
                "coercivity_failed",
                f"physical energy form at {name} fails Sylvester's criterion",
                {"minors": minors, "cluster": name},
            )
        gershgorin = _gershgorin_lower(restricted)
        if gershgorin <= 0:
            raise SGBLLocalSymmetrizerStop(
                "coercivity_failed",
                f"physical energy Gershgorin lower bound at {name} is not positive",
                {"gershgorin": gershgorin, "cluster": name},
            )
        energy[name] = {
            "dimension": 4,
            "sylvester_leading_minors": minors,
            "gershgorin_lower": gershgorin,
        }
        lower = gershgorin if lower is None else min(lower, gershgorin)
    return {
        "clusters": energy,
        "coercivity_lower": lower,
        "mutate_sign": mutate_sign,
    }


def sgbl_seed_floating_spectra(
    matrix: Sequence[Sequence[Fraction]],
) -> dict[str, Any]:
    """Return binary64 eigenvalues labelled as seeds, never as a proof."""

    values, _vectors = _companion_from_seed_float(matrix)
    return {
        "eigenvalues": tuple(complex(value) for value in values),
        "max_imaginary_part": float(np.max(np.abs(values.imag))),
        "proof_uses_floating_eig": False,
        "role": "candidate_seed_only",
    }


def sgbl_injected_defective_companion() -> tuple[tuple[Fraction, ...], ...]:
    """Jordan block ``λ=1`` of size two, padded to 24. Geometric multiplicity 1."""

    matrix = [list(row) for row in _identity_square(FIRST_ORDER_FIELD_COUNT)]
    matrix[0][0] = Q(1)
    matrix[0][1] = Q(1)
    matrix[1][0] = Q(0)
    matrix[1][1] = Q(1)
    return tuple(tuple(row) for row in matrix)


def sgbl_injected_complex_companion() -> tuple[tuple[Fraction, ...], ...]:
    """Real 90-degree block with characteristic polynomial ``λ^2+1``."""

    matrix = [list(row) for row in _zero_square(FIRST_ORDER_FIELD_COUNT)]
    for index in range(2, FIRST_ORDER_FIELD_COUNT):
        matrix[index][index] = Q(1)
    matrix[0][1] = Q(-1)
    matrix[1][0] = Q(1)
    return tuple(tuple(row) for row in matrix)


def sgbl_injected_coalesced_speeds() -> tuple[Fraction, ...]:
    """Two named ESF speeds forced to the same centre."""

    return (Q(1), Q(1), Q(1, 2), Q(-1, 2), Q(1, 3), Q(-1, 3))


def sgbl_assert_distinct_cluster_centers(
    centers: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    """Refuse a coalesced cluster list. Distinct centres are returned unchanged."""

    frozen = tuple(centers)
    if len(set(frozen)) != len(frozen):
        raise SGBLLocalSymmetrizerStop(
            "cluster_coalescence",
            "declared cluster centres are not strictly separated",
            {"centers": frozen},
        )
    return frozen


def _isolated_leading_block(
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[Fraction, Fraction, Fraction, Fraction] | None:
    if len(matrix) < 2:
        return None
    for row in range(2):
        for column in range(2, len(matrix)):
            if matrix[row][column] != 0 or matrix[column][row] != 0:
                return None
    return (matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1])


def _certify_exact_kernels(
    matrix: Sequence[Sequence[Fraction]],
    candidates: Sequence[Fraction],
    *,
    expect_esf_six: bool,
) -> tuple[dict[str, Any], dict[Fraction, tuple[tuple[Fraction, ...], ...]]]:
    clusters: dict[str, Any] = {}
    columns_by_center: dict[Fraction, tuple[tuple[Fraction, ...], ...]] = {}
    seen: set[Fraction] = set()
    for center in candidates:
        if center in seen:
            continue
        kernel = _exact_kernel(matrix, center)
        if not kernel:
            continue
        seen.add(center)
        for vector in kernel:
            image = _matrix_vector(matrix, vector)
            expected = tuple(center * entry for entry in vector)
            if image != expected:
                raise SGBLLocalSymmetrizerStop(
                    "lost_kinetic",
                    f"exact eigenvector identity failed at speed {center}",
                )
        name = _cluster_name(center)
        clusters[name] = {
            "center": center,
            "exact_nullity": len(kernel),
            "eigenvalue_count": len(kernel),
            "method": "exact_rational_kernel",
            "seed_only_floating_spectra": True,
        }
        columns_by_center[center] = kernel
    if expect_esf_six:
        expected = {name: 4 for name, _speed in ESF_CLUSTER_CENTERS}
        got = {
            name: int(cluster["exact_nullity"]) for name, cluster in clusters.items()
        }
        if set(got) == set(expected) and any(
            got[name] != expected[name] for name in expected
        ):
            raise SGBLLocalSymmetrizerStop(
                "cluster_coalescence",
                "ESF six-speed multiplicities coalesced or split",
                {"expected": expected, "got": got},
            )
        if Q(1) in seen and clusters.get("physical_plus", {}).get("exact_nullity") == 8:
            raise SGBLLocalSymmetrizerStop(
                "cluster_coalescence",
                "physical-plus cluster absorbed another speed",
            )
    total = sum(int(cluster["exact_nullity"]) for cluster in clusters.values())
    if total > FIRST_ORDER_FIELD_COUNT:
        raise SGBLLocalSymmetrizerStop(
            "defective_cluster",
            "certified kernel dimensions exceed 24",
            {"total": total},
        )
    return clusters, columns_by_center


def sgbl_krawczyk_simple_eigenpair(
    matrix: Sequence[Sequence[Fraction]],
    *,
    lambda_seed: Fraction,
    vector_seed: Sequence[Fraction],
    limits: SGBLLocalSymmetrizerLimits | None = None,
) -> Mapping[str, Any] | None:
    """Certify a simple real eigenpair in at most 16 variables, or return None."""

    local_limits = limits or SGBLLocalSymmetrizerLimits()
    size = len(matrix)
    if size > 16 or _matrix_bits(matrix) > local_limits.max_krawczyk_matrix_bits:
        return None
    pivot = max(range(size), key=lambda index: abs(vector_seed[index]))
    if vector_seed[pivot] == 0:
        return None
    scale = vector_seed[pivot]
    reduced_seed = [
        vector_seed[index] / scale for index in range(size) if index != pivot
    ]
    center_vars = (lambda_seed, *reduced_seed)
    displacement = tuple(
        Interval(
            -local_limits.krawczyk_displacement, local_limits.krawczyk_displacement
        )
        for _ in center_vars
    )

    def assemble(values: Sequence[Fraction]) -> tuple[Fraction, tuple[Fraction, ...]]:
        lam = values[0]
        vector = []
        cursor = 1
        for index in range(size):
            if index == pivot:
                vector.append(Q(1))
            else:
                vector.append(values[cursor])
                cursor += 1
        return lam, tuple(vector)

    def residual_at(values: Sequence[Fraction]) -> tuple[Interval, ...]:
        lam, vector = assemble(values)
        shifted = tuple(
            tuple(
                matrix[row][column] - (lam if row == column else Q(0))
                for column in range(size)
            )
            for row in range(size)
        )
        image = _matrix_vector(shifted, vector)
        return tuple(Interval.singleton(entry) for entry in image)

    def jacobian_box() -> tuple[tuple[Interval, ...], ...]:
        half = local_limits.krawczyk_displacement
        lam_box = Interval(lambda_seed - half, lambda_seed + half)
        vector_box = []
        cursor = 0
        for index in range(size):
            if index == pivot:
                vector_box.append(Interval.singleton(1))
            else:
                center = reduced_seed[cursor]
                vector_box.append(Interval(center - half, center + half))
                cursor += 1
        rows = []
        for row in range(size):
            line = [-vector_box[row]]
            for index in range(size):
                if index == pivot:
                    continue
                entry = Interval.singleton(matrix[row][index])
                if row == index:
                    entry = entry - lam_box
                line.append(entry)
            rows.append(tuple(line))
        return tuple(rows)

    try:
        center_residual = residual_at(center_vars)
        jacobian = jacobian_box()
        midpoint = tuple(tuple(entry.midpoint() for entry in row) for row in jacobian)
        center_inverse = inverse(midpoint)
        payload = parametric_krawczyk_inclusion(
            center_inverse=center_inverse,
            residual_at_center=center_residual,
            jacobian_box=jacobian,
            displacement_box=displacement,
        )
    except (ValueError, ZeroDivisionError):
        return None
    return {
        "lambda_center": lambda_seed,
        "pivot": pivot,
        "rho_infinity_upper_bound": payload["rho_infinity_upper_bound"],
        "krawczyk_image_strictly_inside_displacement_box": True,
        "method": "rational_interval_parametric_Krawczyk",
    }


def sgbl_certify_companion_eigenframe(
    matrix: Sequence[Sequence[Fraction]],
    *,
    limits: SGBLLocalSymmetrizerLimits | None = None,
    expect_esf_six: bool = False,
    require_complete: bool = False,
) -> dict[str, Any]:
    """Certify real clusters of one exact rational 24-by-24 companion."""

    if len(matrix) != FIRST_ORDER_FIELD_COUNT or any(
        len(row) != FIRST_ORDER_FIELD_COUNT for row in matrix
    ):
        raise ValueError("companion must be 24-by-24")
    local_limits = limits or SGBLLocalSymmetrizerLimits()
    frozen = tuple(tuple(Q(entry) for entry in row) for row in matrix)
    payload = dict(
        _cached_certify_companion(
            frozen,
            expect_esf_six,
            local_limits.max_seed_denominator,
            local_limits.max_seed_candidates,
            local_limits.max_krawczyk_matrix_bits,
            local_limits.krawczyk_displacement,
        )
    )
    if require_complete and not payload["complete_24_mode_basis"]:
        raise SGBLLocalSymmetrizerStop(
            "eigenframe_incomplete",
            "companion eigenframe is not a complete rank-24 real basis",
            {
                "rank": payload["eigenframe_rank"],
                "certified_modes": payload["certified_mode_count"],
            },
        )
    return payload


@lru_cache(maxsize=8)
def _cached_certify_companion(
    matrix: tuple[tuple[Fraction, ...], ...],
    expect_esf_six: bool,
    max_seed_denominator: int,
    max_seed_candidates: int,
    max_krawczyk_matrix_bits: int,
    krawczyk_displacement: Fraction,
) -> dict[str, Any]:
    local_limits = SGBLLocalSymmetrizerLimits(
        max_seed_denominator=max_seed_denominator,
        max_seed_candidates=max_seed_candidates,
        max_krawczyk_matrix_bits=max_krawczyk_matrix_bits,
        krawczyk_displacement=krawczyk_displacement,
    )
    block = _isolated_leading_block(matrix)
    if block is not None:
        trace = block[0] + block[3]
        det = block[0] * block[3] - block[1] * block[2]
        discriminant = trace * trace - 4 * det
        if discriminant < 0:
            raise SGBLLocalSymmetrizerStop(
                "complex_spectrum",
                "isolated 2-block has negative discriminant, so no real eigenvalues",
                {"block": block, "discriminant": discriminant},
            )
        if discriminant == 0 and block[1] != 0:
            geometric = len(
                _exact_kernel(tuple(row[:2] for row in matrix[:2]), block[0])
            )
            if geometric < 2:
                raise SGBLLocalSymmetrizerStop(
                    "defective_cluster",
                    "isolated 2-block is defective: geometric multiplicity is below algebraic",
                    {"block": block, "geometric_multiplicity": geometric},
                )
    seeds = sgbl_seed_floating_spectra(matrix)
    radius = _infinity_norm(matrix)
    candidates = _rational_candidates_from_seeds(
        seeds["eigenvalues"],
        spectral_radius=radius,
        limits=local_limits,
    )
    clusters, columns_by_center = _certify_exact_kernels(
        matrix,
        candidates,
        expect_esf_six=expect_esf_six,
    )
    ordered_clusters: dict[str, Any] = {}
    for name, _speed in ESF_CLUSTER_CENTERS:
        if name in clusters:
            ordered_clusters[name] = clusters[name]
    for name, cluster in clusters.items():
        if name not in ordered_clusters:
            ordered_clusters[name] = cluster
    clusters = ordered_clusters
    columns: list[tuple[Fraction, ...]] = []
    used: set[Fraction] = set()
    for _name, speed in ESF_CLUSTER_CENTERS:
        if speed in columns_by_center:
            columns.extend(columns_by_center[speed])
            used.add(speed)
    for center, kernel in columns_by_center.items():
        if center not in used:
            columns.extend(kernel)
    if not columns:
        raise SGBLLocalSymmetrizerStop(
            "eigenframe_incomplete",
            "no exact rational eigenvalue kernels were certified",
        )
    frame = _columns_to_matrix(columns)
    frame_rank = rank(frame)
    total = sum(int(cluster["exact_nullity"]) for cluster in clusters.values())
    complete = (
        total == FIRST_ORDER_FIELD_COUNT
        and frame_rank == FIRST_ORDER_FIELD_COUNT
        and len(columns) == FIRST_ORDER_FIELD_COUNT
    )
    payload: dict[str, Any] = {
        "clusters": clusters,
        "certified_mode_count": total,
        "eigenframe_rank": frame_rank,
        "complete_24_mode_basis": complete,
        "seed": {
            "role": "candidate_seed_only",
            "proof_uses_floating_eig": False,
            "max_imaginary_part": seeds["max_imaginary_part"],
        },
        "columns": tuple(columns),
        "columns_by_center": columns_by_center,
        "frame": frame if complete else None,
        "riesz_projectors": None,
        "positive_symmetrizer": None,
        "obstruction": None,
    }
    if not complete:
        payload["obstruction"] = NONFLAT_OBSTRUCTION
        payload["krawczyk_attempted"] = (
            _matrix_bits(matrix) <= local_limits.max_krawczyk_matrix_bits
        )
        return payload
    try:
        frame_inverse = inverse(frame)
    except ValueError as exc:
        raise SGBLLocalSymmetrizerStop(
            "eigenframe_incomplete",
            "complete column set is singular",
        ) from exc
    reconstructed = matrix_multiply(
        matrix_multiply(frame, _diag_from_clusters(clusters)),
        frame_inverse,
    )
    if reconstructed != tuple(tuple(row) for row in matrix):
        raise SGBLLocalSymmetrizerStop(
            "lost_kinetic",
            "exact reconstruction V Λ V^{-1} = M failed",
        )
    payload["frame"] = frame
    payload["frame_inverse"] = frame_inverse
    payload["riesz_projectors"] = _riesz_projectors(
        matrix, frame, frame_inverse, clusters
    )
    payload["positive_symmetrizer"] = _positive_symmetrizer(
        frame, frame_inverse, matrix
    )
    payload["reconstruction"] = "exact_V_Lambda_V_inverse"
    return payload


def _diag_from_clusters(
    clusters: Mapping[str, Mapping[str, Any]],
) -> tuple[tuple[Fraction, ...], ...]:
    values: list[Fraction] = []
    for cluster in clusters.values():
        values.extend([cluster["center"]] * int(cluster["exact_nullity"]))
    size = len(values)
    return tuple(
        tuple(values[column] if row == column else Q(0) for column in range(size))
        for row in range(size)
    )


def _unnormalized_coefficients(
    box: SGBLPrincipalBackgroundBox,
    *,
    include_gauge: bool,
) -> dict[str, Any]:
    counter = _EvaluationCounter(box.limits.max_residual_evaluations)
    try:
        return _extract_coefficient_tensors(
            box,
            include_gauge=include_gauge,
            sqrt2=_sqrt2_enclosure(box.limits.sqrt2_bits),
            counter=counter,
            normalize=False,
        )
    except SGBLConeStop as exc:
        raise _wrap_cone_stop(exc) from exc


def _radial_companion_and_action(
    box: SGBLPrincipalBackgroundBox,
) -> tuple[
    tuple[tuple[Fraction, ...], ...],
    tuple[tuple[Fraction, ...], ...],
    tuple[tuple[Fraction, ...], ...],
    Mapping[str, Any] | None,
    tuple[tuple[Interval, ...], ...],
]:
    complete = _unnormalized_coefficients(box, include_gauge=True)
    action = _unnormalized_coefficients(box, include_gauge=False)
    try:
        time_time = _as_rational_matrix(complete["time_time"], name="local A")
        time_space = _as_rational_matrix(complete["time_space"][0], name="local B_r")
        space_space = _as_rational_matrix(
            complete["space_space"][0][0], name="local C_rr"
        )
        action_a = _as_rational_matrix(action["time_time"], name="local action A")
        action_b = _as_rational_matrix(action["time_space"][0], name="local action B_r")
        companion = _exact_radial_companion(time_time, time_space, space_space)
    except SGBLConeStop as exc:
        raise _wrap_cone_stop(exc) from exc
    inverse_payload = _interval_inverse(complete["time_time"])
    return companion, action_a, action_b, inverse_payload, complete["time_time"]


def sgbl_coefficientwise_physical_scalar_identities(
    box: SGBLPrincipalBackgroundBox,
) -> dict[str, Any]:
    """Coefficientwise ``(M(n)-s I) v_{phi,chi} = (|n|^2-1) w`` on ESF tensors.

    These identities are quadratic polynomial equalities in the spatial
    covector.  They are not sampled-direction evidence.  They certify two
    physical scalar modes in every direction on the unit sphere; they do
    not by themselves supply the remaining twenty-two modes.
    """

    if type(box) is not SGBLPrincipalBackgroundBox:
        raise TypeError("box must be SGBLPrincipalBackgroundBox")
    if not box.curvature_free:
        return {
            "holds": False,
            "reason": "coefficientwise_scalar_identities_require_curvature_free_tensors",
            "complete_24_mode_all_direction": False,
        }
    complete = _unnormalized_coefficients(box, include_gauge=True)
    time_time = _as_rational_matrix(complete["time_time"], name="A")
    blocks_b = [
        _as_rational_matrix(complete["time_space"][index], name=f"B_{index}")
        for index in range(3)
    ]
    blocks_c = [
        [
            _as_rational_matrix(
                complete["space_space"][first][second],
                name=f"C_{first}{second}",
            )
            for second in range(3)
        ]
        for first in range(3)
    ]

    def companion_at(
        n: tuple[Fraction, Fraction, Fraction],
    ) -> tuple[tuple[Fraction, ...], ...]:
        mixed_b = _zero_square(FIELD_COUNT)
        mixed_c = _zero_square(FIELD_COUNT)
        for index in range(3):
            mixed_b = _matrix_add(mixed_b, _matrix_scale(blocks_b[index], n[index]))
            for other in range(3):
                mixed_c = _matrix_add(
                    mixed_c,
                    _matrix_scale(blocks_c[index][other], n[index] * n[other]),
                )
        return _exact_radial_companion(time_time, mixed_b, mixed_c)

    samples = (
        (Q(0), Q(0), Q(0)),
        (Q(1), Q(0), Q(0)),
        (Q(0), Q(1), Q(0)),
        (Q(0), Q(0), Q(1)),
        (Q(1), Q(1), Q(0)),
        (Q(1), Q(0), Q(1)),
        (Q(0), Q(1), Q(1)),
        (Q(2), Q(0), Q(0)),
        (Q(0), Q(2), Q(0)),
        (Q(0), Q(0), Q(2)),
    )
    identities = []
    for field_index, label in (
        (PHI_FIELD_INDEX, "phi"),
        (CHI_FIELD_INDEX, "chi"),
    ):
        for speed in (Q(1), Q(-1)):
            spatial = _standard_basis(FIELD_COUNT, field_index)
            vector = spatial + tuple(speed * entry for entry in spatial)
            residual_at_origin = None
            holds = True
            for sample in samples:
                image = _matrix_vector(companion_at(sample), vector)
                residual = tuple(
                    image[index] - speed * vector[index] for index in range(len(vector))
                )
                radius = (
                    sample[0] * sample[0]
                    + sample[1] * sample[1]
                    + sample[2] * sample[2]
                )
                if sample == (Q(0), Q(0), Q(0)):
                    residual_at_origin = residual
                    cofactor = tuple(-entry for entry in residual)
                else:
                    expected = tuple(
                        (radius - 1) * (-entry) for entry in residual_at_origin
                    )
                    if residual != expected:
                        holds = False
                        break
            identities.append(
                {
                    "field": label,
                    "speed": speed,
                    "holds": holds,
                    "cofactor_support": tuple(
                        index
                        for index, entry in enumerate(cofactor or ())
                        if entry != 0
                    ),
                }
            )
    holds = all(item["holds"] for item in identities)
    return {
        "holds": holds,
        "identities": tuple(identities),
        "polynomial_degree": 2,
        "sampled_directions_are_not_the_proof": True,
        "determining_quadratic_grid": True,
        "complete_24_mode_all_direction": False,
        "all_direction_status": ALL_DIRECTION_INCOMPLETE,
    }


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLLocalSymmetrizerRecord:
    """Point-local cone/symmetrizer record. Aggregate branch health stays false."""

    box: SGBLPrincipalBackgroundBox
    branch: str
    action: Mapping[str, Fraction | str]
    coefficient_enclosure: SGBLCoefficientEnclosure
    theorem: str
    certified_direction: tuple[int, int, int]
    clusters: Mapping[str, Any]
    riesz_projectors: Mapping[str, Any] | None
    kinetic_rho_infinity: Fraction | None
    kinetic_singular_value_lower: Fraction | None
    physical_energy_coercivity_lower: Fraction | None
    exact_eigenframe_rank: int | None
    certified_mode_count: int
    real_complete_basis: bool
    positive_symmetrizer: bool
    all_direction_status: str
    coefficientwise_physical_scalar_identities: bool
    classification: str
    inconclusive_reason: str | None
    obstruction: str | None
    family_covering_proved: bool
    family_box_entrypoint: bool
    sha256: str
    structural_hessian_symmetric: bool
    structural_pure_gauge_identity: bool
    proof_uses_floating_eig: bool
    proof_uses_sampled_directions: bool
    proof_uses_fitted_gap: bool
    proof_uses_esf_perturbation: bool

    def __post_init__(self) -> None:
        if type(self.box) is not SGBLPrincipalBackgroundBox:
            raise TypeError("box must be SGBLPrincipalBackgroundBox")
        if self.branch != "SGB-L":
            raise ValueError("local symmetrizer record is owned by the linear branch")
        if self.classification not in {
            PROVED_CLASSIFICATION,
            INCONCLUSIVE_CLASSIFICATION,
            INCOMPLETE_CLASSIFICATION,
        }:
            raise ValueError("unknown local-symmetrizer classification")
        if self.all_direction_status not in {
            ALL_DIRECTION_COEFFICIENTWISE,
            ALL_DIRECTION_INCOMPLETE,
        }:
            raise ValueError("unknown all-direction status")
        if self.proof_uses_floating_eig or self.proof_uses_sampled_directions:
            raise ValueError("floating eig or sampled directions cannot be the proof")
        if self.proof_uses_fitted_gap or self.proof_uses_esf_perturbation:
            raise ValueError("fitted gaps and ESF perturbation are not this owner")
        if self.family_covering_proved:
            raise ValueError("family covering is a later owner")
        if self.theorem != LOCAL_THEOREM:
            raise ValueError("local theorem text is immutable")
        if self.classification == PROVED_CLASSIFICATION:
            if not (
                self.real_complete_basis
                and self.positive_symmetrizer
                and self.inconclusive_reason is None
                and self.obstruction is None
                and self.exact_eigenframe_rank == FIRST_ORDER_FIELD_COUNT
                and self.certified_mode_count == FIRST_ORDER_FIELD_COUNT
            ):
                raise ValueError("proved local cone requires a complete 24-mode basis")
            if (
                self.kinetic_singular_value_lower is None
                or self.kinetic_singular_value_lower <= 0
            ):
                raise ValueError("proved local cone requires a positive kinetic margin")
            if (
                self.physical_energy_coercivity_lower is None
                or self.physical_energy_coercivity_lower <= 0
            ):
                raise ValueError(
                    "proved local cone requires positive energy coercivity"
                )
            if (
                not self.structural_hessian_symmetric
                or not self.structural_pure_gauge_identity
            ):
                raise ValueError("proved local cone requires algebraic identities")
            if self.all_direction_status == ALL_DIRECTION_COEFFICIENTWISE:
                if not self.coefficientwise_physical_scalar_identities:
                    raise ValueError(
                        "coefficientwise all-direction requires verified construction"
                    )
        else:
            if self.real_complete_basis or self.positive_symmetrizer:
                raise ValueError(
                    "incomplete or inconclusive local cone cannot retain a real basis"
                )
        object.__setattr__(self, "action", MappingProxyType(dict(self.action)))
        object.__setattr__(self, "clusters", MappingProxyType(dict(self.clusters)))
        if self.riesz_projectors is not None:
            object.__setattr__(
                self,
                "riesz_projectors",
                MappingProxyType(dict(self.riesz_projectors)),
            )
        expected = _record_digest(
            classification=self.classification,
            certified_mode_count=self.certified_mode_count,
            eigenframe_rank=self.exact_eigenframe_rank,
            obstruction=self.obstruction,
            clusters=self.clusters,
        )
        if self.sha256 != expected:
            raise ValueError(
                "local symmetrizer digest does not match the certified data"
            )

    @property
    def strongly_hyperbolic(self) -> bool:
        return self.classification == PROVED_CLASSIFICATION

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def quantitative_all_covector_weak_coupling_health_envelope_passed(self) -> bool:
        return False

    @property
    def representative_generator_sufficient(self) -> bool:
        return False

    @property
    def later_family_owner(self) -> str:
        return str(FAMILY_BOX_API["later_owner"])


def _record_digest(
    *,
    classification: str,
    certified_mode_count: int,
    eigenframe_rank: int | None,
    obstruction: str | None,
    clusters: Mapping[str, Any],
) -> str:
    parts = [
        classification,
        str(certified_mode_count),
        str(eigenframe_rank),
        str(obstruction),
    ]
    for name, cluster in clusters.items():
        parts.append(name)
        parts.append(str(cluster.get("center")))
        parts.append(str(cluster.get("exact_nullity")))
    return _sha256_rationals(parts)


def sgbl_local_symmetrizer_family_contract() -> Mapping[str, Any]:
    """Expose the later family-box API without claiming family covering."""

    return dict(FAMILY_BOX_API)


def _action_map(box: SGBLPrincipalBackgroundBox) -> dict[str, Fraction | str]:
    return {
        "planck_mass": box.planck_mass.midpoint(),
        "alpha_gb": box.alpha_gb.midpoint(),
        "beta": Q(0),
        "eta": Q(0),
        "F_prime": Q(0),
        "chart": ORTHONORMAL_MINKOWSKI_CHART,
    }


def _build_record(
    *,
    box: SGBLPrincipalBackgroundBox,
    enclosure: SGBLCoefficientEnclosure,
    clusters: Mapping[str, Any],
    riesz: Mapping[str, Any] | None,
    kinetic_rho: Fraction | None,
    kinetic_lower: Fraction | None,
    energy_lower: Fraction | None,
    eigenframe_rank: int | None,
    certified_modes: int,
    real_basis: bool,
    positive_symmetrizer: bool,
    all_direction: str,
    scalar_identities: bool,
    classification: str,
    inconclusive_reason: str | None,
    obstruction: str | None,
    structural: Mapping[str, Any],
    family_entrypoint: bool,
) -> SGBLLocalSymmetrizerRecord:
    digest = _record_digest(
        classification=classification,
        certified_mode_count=certified_modes,
        eigenframe_rank=eigenframe_rank,
        obstruction=obstruction,
        clusters=clusters,
    )
    return SGBLLocalSymmetrizerRecord(
        box=box,
        branch="SGB-L",
        action=_action_map(box),
        coefficient_enclosure=enclosure,
        theorem=LOCAL_THEOREM,
        certified_direction=LOCAL_DIRECTION,
        clusters=clusters,
        riesz_projectors=riesz,
        kinetic_rho_infinity=kinetic_rho,
        kinetic_singular_value_lower=kinetic_lower,
        physical_energy_coercivity_lower=energy_lower,
        exact_eigenframe_rank=eigenframe_rank,
        certified_mode_count=certified_modes,
        real_complete_basis=real_basis,
        positive_symmetrizer=positive_symmetrizer,
        all_direction_status=all_direction,
        coefficientwise_physical_scalar_identities=scalar_identities,
        classification=classification,
        inconclusive_reason=inconclusive_reason,
        obstruction=obstruction,
        family_covering_proved=False,
        family_box_entrypoint=family_entrypoint,
        sha256=digest,
        structural_hessian_symmetric=bool(structural["structural_hessian_symmetric"]),
        structural_pure_gauge_identity=bool(
            structural["structural_pure_gauge_identity"]
        ),
        proof_uses_floating_eig=False,
        proof_uses_sampled_directions=False,
        proof_uses_fitted_gap=False,
        proof_uses_esf_perturbation=False,
    )


def sgbl_local_symmetrizer_certificate(
    box: SGBLPrincipalBackgroundBox,
    *,
    limits: SGBLLocalSymmetrizerLimits | None = None,
    family_box_entrypoint: bool = False,
    mutate_energy_sign: int = 1,
    omit_eigenframe_columns: int = 0,
) -> SGBLLocalSymmetrizerRecord:
    """Prove or fail-closed the point-local SGB-L cone on one declared box."""

    if type(box) is not SGBLPrincipalBackgroundBox:
        raise TypeError("box must be SGBLPrincipalBackgroundBox")
    local_limits = limits or SGBLLocalSymmetrizerLimits()
    try:
        enclosure = sgbl_enclose_principal_coefficients(box)
        structural = sgbl_structural_identities(box)
    except SGBLConeStop as exc:
        raise _wrap_cone_stop(exc) from exc
    if (
        not structural["structural_hessian_symmetric"]
        or not structural["structural_pure_gauge_identity"]
    ):
        raise SGBLLocalSymmetrizerStop(
            "broken_gauge",
            "algebraic Hessian symmetry or pure-gauge identity failed",
            structural,
        )
    if structural["representative_generator_sufficient"] is not False:
        raise SGBLLocalSymmetrizerStop(
            "broken_gauge",
            "a single representative generator is not a complete basis",
            structural,
        )
    scalar_identities = sgbl_coefficientwise_physical_scalar_identities(box)
    scalar_holds = bool(scalar_identities.get("holds"))
    all_direction = ALL_DIRECTION_INCOMPLETE
    if not box.is_singleton:
        return _build_record(
            box=box,
            enclosure=enclosure,
            clusters={},
            riesz=None,
            kinetic_rho=None,
            kinetic_lower=None,
            energy_lower=None,
            eigenframe_rank=None,
            certified_modes=0,
            real_basis=False,
            positive_symmetrizer=False,
            all_direction=all_direction,
            scalar_identities=scalar_holds,
            classification=INCONCLUSIVE_CLASSIFICATION,
            inconclusive_reason="interval_companion_not_singleton",
            obstruction=None,
            structural=structural,
            family_entrypoint=family_box_entrypoint,
        )
    companion, action_a, action_b, inverse_payload, time_time = (
        _radial_companion_and_action(box)
    )
    if inverse_payload is None:
        return _build_record(
            box=box,
            enclosure=enclosure,
            clusters={},
            riesz=None,
            kinetic_rho=None,
            kinetic_lower=None,
            energy_lower=None,
            eigenframe_rank=None,
            certified_modes=0,
            real_basis=False,
            positive_symmetrizer=False,
            all_direction=all_direction,
            scalar_identities=scalar_holds,
            classification=INCONCLUSIVE_CLASSIFICATION,
            inconclusive_reason="neumann_rho_not_below_one",
            obstruction=None,
            structural=structural,
            family_entrypoint=family_box_entrypoint,
        )
    kinetic_lower = _kinetic_singular_value_lower(
        time_time,
        inverse_payload,
        box.limits.sqrt2_bits,
    )
    if kinetic_lower <= 0:
        raise SGBLLocalSymmetrizerStop(
            "lost_kinetic",
            "kinetic singular-value lower bound is not positive",
        )
    payload = sgbl_certify_companion_eigenframe(
        companion,
        limits=local_limits,
        expect_esf_six=box.curvature_free,
    )
    columns = list(payload["columns"])
    if omit_eigenframe_columns:
        if omit_eigenframe_columns >= len(columns):
            raise SGBLLocalSymmetrizerStop(
                "eigenframe_incomplete",
                "omitting the eigenframe left no columns",
            )
        raise SGBLLocalSymmetrizerStop(
            "eigenframe_incomplete",
            "injected eigenframe omission dropped the rank-24 basis",
            {"omitted": omit_eigenframe_columns},
        )
    if not payload["complete_24_mode_basis"]:
        return _build_record(
            box=box,
            enclosure=enclosure,
            clusters=payload["clusters"],
            riesz=None,
            kinetic_rho=inverse_payload["rho_infinity"],
            kinetic_lower=kinetic_lower,
            energy_lower=None,
            eigenframe_rank=payload["eigenframe_rank"],
            certified_modes=payload["certified_mode_count"],
            real_basis=False,
            positive_symmetrizer=False,
            all_direction=all_direction,
            scalar_identities=scalar_holds,
            classification=INCOMPLETE_CLASSIFICATION,
            inconclusive_reason=None,
            obstruction=payload["obstruction"] or NONFLAT_OBSTRUCTION,
            structural=structural,
            family_entrypoint=family_box_entrypoint,
        )
    energy = _physical_energy(
        action_b,
        action_a,
        payload["columns_by_center"],
        mutate_sign=mutate_energy_sign,
    )
    return _build_record(
        box=box,
        enclosure=enclosure,
        clusters=payload["clusters"],
        riesz=payload["riesz_projectors"],
        kinetic_rho=inverse_payload["rho_infinity"],
        kinetic_lower=kinetic_lower,
        energy_lower=energy["coercivity_lower"],
        eigenframe_rank=payload["eigenframe_rank"],
        certified_modes=payload["certified_mode_count"],
        real_basis=True,
        positive_symmetrizer=True,
        all_direction=all_direction,
        scalar_identities=scalar_holds,
        classification=PROVED_CLASSIFICATION,
        inconclusive_reason=None,
        obstruction=None,
        structural=structural,
        family_entrypoint=family_box_entrypoint,
    )


def sgbl_local_symmetrizer_from_family_box(
    box: SGBLPrincipalBackgroundBox,
    *,
    limits: SGBLLocalSymmetrizerLimits | None = None,
) -> SGBLLocalSymmetrizerRecord:
    """Later family-box entrypoint. Same local certificate; covering stays false."""

    return sgbl_local_symmetrizer_certificate(
        box,
        limits=limits,
        family_box_entrypoint=True,
    )


def sgbl_local_physical_energy(
    box: SGBLPrincipalBackgroundBox,
    *,
    mutate_sign: int = 1,
) -> dict[str, Any]:
    """Physical-energy restriction, with an optional coercivity-attack sign."""

    companion, action_a, action_b, _inverse, _time_time = _radial_companion_and_action(
        box
    )
    payload = sgbl_certify_companion_eigenframe(
        companion,
        expect_esf_six=box.curvature_free,
        require_complete=True,
    )
    return _physical_energy(
        action_b,
        action_a,
        payload["columns_by_center"],
        mutate_sign=mutate_sign,
    )


__all__ = [
    "ALL_DIRECTION_COEFFICIENTWISE",
    "ALL_DIRECTION_INCOMPLETE",
    "ESF_CLUSTER_CENTERS",
    "ESF_SPEEDS",
    "FAMILY_BOX_API",
    "FORBIDDEN_HEALTH_IMPORTS",
    "INCOMPLETE_CLASSIFICATION",
    "INCONCLUSIVE_CLASSIFICATION",
    "LOCAL_DIRECTION",
    "LOCAL_THEOREM",
    "NONFLAT_OBSTRUCTION",
    "PROVED_CLASSIFICATION",
    "SGBLLocalSymmetrizerLimits",
    "SGBLLocalSymmetrizerRecord",
    "SGBLLocalSymmetrizerStop",
    "sgbl_certify_companion_eigenframe",
    "sgbl_coefficientwise_physical_scalar_identities",
    "sgbl_assert_distinct_cluster_centers",
    "sgbl_injected_coalesced_speeds",
    "sgbl_injected_complex_companion",
    "sgbl_injected_defective_companion",
    "sgbl_krawczyk_simple_eigenpair",
    "sgbl_local_physical_energy",
    "sgbl_local_symmetrizer_certificate",
    "sgbl_local_symmetrizer_family_contract",
    "sgbl_local_symmetrizer_from_family_box",
    "sgbl_seed_floating_spectra",
    "sgbl_esf_principal_box",
    "sgbl_nonflat_source_principal_box",
    "sgbl_pure_gauge_identity_defect",
]
