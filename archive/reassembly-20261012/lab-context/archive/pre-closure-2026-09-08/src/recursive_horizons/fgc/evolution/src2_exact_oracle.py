"""Exact-rational oracle for one captured SRC2 binary64 ADM source point.

This module does not alter the GR-0 source, its tolerance, or an evolution
trajectory.  It treats a supplied binary64 lower two-jet as an exact dyadic
input, evaluates the same full REF1 equations with :class:`Fraction`, and
solves only the resulting pointwise six-by-six affine acceleration system.

It is deliberately narrower than a continuum result: exactness here means
exact evaluation of the *stored binary64 point data*, not exact continuum
initial data or an FGC-QR mechanism conclusion.
"""

from __future__ import annotations

from fractions import Fraction
from math import isfinite
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from ..action import ActionParameters, ModelID
from ..modified_harmonic_reference import modified_harmonic_full_residuals
from ..reference_connection import flat_spherical_annulus_reference
from ..spherical_reduction import Jet2, _state_from_adm_pg_jets
from .gr0_direct_source import gr0_ref1_residual_batch
from .health_monitor import FIELD_ORDER
from .vectorized_source import batch_ref1_residual


Q = Fraction
FIELD_COUNT = len(FIELD_ORDER)
ADM_EXACT_FIELD_ORDER = (
    "alpha",
    "shift",
    "lambda",
    "areal_radius",
    "phi",
    "chi",
)


def exact_dyadic(value: object) -> Fraction:
    """Return the exact rational value represented by one finite binary64 input."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError("SRC2 exact inputs must be finite real binary64 values")
    binary64 = float(value)
    if not isfinite(binary64):
        raise ValueError("SRC2 exact inputs must be finite")
    return Q.from_float(binary64)


def _six(name: str, value: object) -> tuple[float, ...]:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (FIELD_COUNT,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be one finite six-field binary64 row")
    return tuple(float(item) for item in array)


def _radius(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError("coordinate radius must be a finite binary64 scalar")
    answer = float(value)
    if not isfinite(answer) or answer <= 0.0:
        raise ValueError("coordinate radius must be finite and positive")
    return answer


def _as_fraction_vector(value: Sequence[float]) -> tuple[Fraction, ...]:
    return tuple(exact_dyadic(item) for item in value)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction_vector_text(values: Sequence[Fraction]) -> list[str]:
    return [_fraction_text(value) for value in values]


def _float_hex_vector(values: Sequence[float | Fraction]) -> list[str]:
    return [float(value).hex() for value in values]


def _exact_residual(
    u: Sequence[float],
    p: Sequence[float],
    q: Sequence[float],
    p_r: Sequence[float],
    q_r: Sequence[float],
    radius: float,
    acceleration: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    """Evaluate the complete exact GR-0 REF1 vector for one acceleration."""

    values = tuple(_as_fraction_vector(row) for row in (u, p, q, p_r, q_r))
    if len(acceleration) != FIELD_COUNT:
        raise ValueError("exact acceleration must have six fields")
    second = tuple(
        item if isinstance(item, Fraction) else Q(item) for item in acceleration
    )
    jets = {
        name: Jet2(
            values[0][index],
            values[1][index],
            values[2][index],
            second[index],
            values[3][index],
            values[4][index],
        )
        for index, name in enumerate(ADM_EXACT_FIELD_ORDER)
    }
    exact_radius = exact_dyadic(radius)
    state = _state_from_adm_pg_jets(
        {
            "model_id": "GR-0",
            "action_parameters": {
                "planck_mass": Q(2),
                "scalar_mass": Q(3),
                "quartic_coupling": Q(1, 2),
            },
        },
        jets,
    )
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=exact_radius / 2
    )
    residual = modified_harmonic_full_residuals(
        state,
        reference=reference,
        coordinate_radius=exact_radius,
        tilde_normal_factor=4,
        hat_normal_factor=9,
    )["full_residual_vector"]
    return tuple(residual)


def exact_ref1_residual(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radius: object,
    acceleration: object,
) -> tuple[Fraction, ...]:
    """Public exact evaluator for a binary64 lower jet and acceleration."""

    rows = tuple(_six(name, value) for name, value in (
        ("u", u), ("p", p), ("q", q), ("p_r", p_r), ("q_r", q_r)
    ))
    binary_radius = _radius(radius)
    binary_acceleration = _six("acceleration", acceleration)
    return _exact_residual(
        *rows,
        binary_radius,
        _as_fraction_vector(binary_acceleration),
    )


def exact_affine_system(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radius: object,
) -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, ...], ...]]:
    """Extract ``b,A`` from exact zero/unit acceleration evaluations."""

    rows = tuple(_six(name, value) for name, value in (
        ("u", u), ("p", p), ("q", q), ("p_r", p_r), ("q_r", q_r)
    ))
    binary_radius = _radius(radius)
    zero = (Q(0),) * FIELD_COUNT
    constant = _exact_residual(*rows, binary_radius, zero)
    columns: list[tuple[Fraction, ...]] = []
    for column in range(FIELD_COUNT):
        seed = list(zero)
        seed[column] = Q(1)
        evaluated = _exact_residual(*rows, binary_radius, tuple(seed))
        columns.append(tuple(evaluated[row] - constant[row] for row in range(FIELD_COUNT)))
    matrix = tuple(
        tuple(columns[column][row] for column in range(FIELD_COUNT))
        for row in range(FIELD_COUNT)
    )
    return constant, matrix


def exact_gaussian_solve(
    matrix: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    """Solve one exact six-by-six system, failing closed on a zero pivot."""

    if len(matrix) != FIELD_COUNT or len(rhs) != FIELD_COUNT:
        raise ValueError("exact affine solve requires six rows")
    augmented = [
        [Q(value) for value in matrix[row]] + [Q(rhs[row])]
        for row in range(FIELD_COUNT)
    ]
    if any(len(row) != FIELD_COUNT + 1 for row in augmented):
        raise ValueError("exact affine solve requires a six-by-six matrix")
    for column in range(FIELD_COUNT):
        pivot = next(
            (row for row in range(column, FIELD_COUNT) if augmented[row][column] != 0),
            None,
        )
        if pivot is None:
            raise ValueError("exact affine system is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(FIELD_COUNT):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    left - factor * right
                    for left, right in zip(augmented[row], augmented[column], strict=True)
                ]
    return tuple(augmented[row][-1] for row in range(FIELD_COUNT))


def _float_residuals(
    u: tuple[float, ...],
    p: tuple[float, ...],
    q: tuple[float, ...],
    p_r: tuple[float, ...],
    q_r: tuple[float, ...],
    radius: float,
    acceleration: tuple[float, ...],
) -> dict[str, list[float]]:
    arrays = [np.asarray(row, dtype=np.float64)[None, :] for row in (u, p, q, p_r, q_r)]
    seeded = np.asarray(acceleration, dtype=np.float64)[None, None, :]
    radii = np.asarray((radius,), dtype=np.float64)
    direct = gr0_ref1_residual_batch(
        arrays[0], arrays[1], arrays[2], seeded, arrays[3], arrays[4], radii
    ).full_residual[0, 0]
    action = ActionParameters(
        ModelID.GR_0,
        planck_mass=2.0,
        scalar_mass=3.0,
        quartic_coupling=0.5,
        pulse_width=1.0,
    )
    generic = batch_ref1_residual(
        arrays[0], arrays[1], arrays[2], seeded, arrays[3], arrays[4], radii,
        action=action,
    )[0, 0]
    return {"direct": [float(value) for value in direct], "generic": [float(value) for value in generic]}


def diagnose_exact_oracle(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radius: object,
    baseline_acceleration: object,
) -> dict[str, Any]:
    """Return deterministic serializable evidence for a compact SRC2 point fixture."""

    lower = tuple(_six(name, value) for name, value in (
        ("u", u), ("p", p), ("q", q), ("p_r", p_r), ("q_r", q_r)
    ))
    binary_radius = _radius(radius)
    baseline = _six("baseline_acceleration", baseline_acceleration)
    baseline_exact = _exact_residual(*lower, binary_radius, _as_fraction_vector(baseline))
    constant, matrix = exact_affine_system(*lower, binary_radius)
    root = exact_gaussian_solve(matrix, tuple(-item for item in constant))
    root_exact = _exact_residual(*lower, binary_radius, root)
    rounded = tuple(float(item) for item in root)
    rounded_exact = _exact_residual(*lower, binary_radius, _as_fraction_vector(rounded))
    float_baseline = _float_residuals(*lower, binary_radius, baseline)
    float_rounded = _float_residuals(*lower, binary_radius, rounded)

    def infinity(values: Sequence[Fraction | float]) -> float:
        return float(max((abs(value) for value in values), default=Q(0)))

    return {
        "schema_version": 1,
        "field_order": list(FIELD_ORDER),
        "branch": "GR-0",
        "action_parameters": {"planck_mass": 2, "scalar_mass": 3, "quartic_coupling": "1/2"},
        "tilde_normal_factor": 4,
        "hat_normal_factor": 9,
        "radius_binary64_hex": binary_radius.hex(),
        "lower_jet_binary64_hex": {
            name: _float_hex_vector(row)
            for name, row in zip(("u", "p", "q", "p_r", "q_r"), lower, strict=True)
        },
        "baseline_acceleration_binary64_hex": _float_hex_vector(baseline),
        "baseline_exact_residual": _fraction_vector_text(baseline_exact),
        "baseline_exact_residual_infinity": infinity(baseline_exact),
        "exact_affine_nonsingular": True,
        "exact_root_residual_is_zero": all(value == 0 for value in root_exact),
        "exact_root_fraction": _fraction_vector_text(root),
        "rounded_exact_root_binary64_hex": _float_hex_vector(rounded),
        "rounded_exact_root_residual": _fraction_vector_text(rounded_exact),
        "rounded_exact_root_residual_infinity": infinity(rounded_exact),
        "floating_residuals": {
            "baseline": float_baseline,
            "rounded_exact_root": float_rounded,
        },
        "floating_residual_infinities": {
            "baseline_direct": infinity(float_baseline["direct"]),
            "baseline_generic": infinity(float_baseline["generic"]),
            "rounded_exact_root_direct": infinity(float_rounded["direct"]),
            "rounded_exact_root_generic": infinity(float_rounded["generic"]),
        },
        "continuum_equations_changed": False,
        "raw_residual_threshold_changed": False,
        "mechanism_question_answered": False,
    }


__all__ = [
    "ADM_EXACT_FIELD_ORDER",
    "exact_dyadic",
    "exact_ref1_residual",
    "exact_affine_system",
    "exact_gaussian_solve",
    "diagnose_exact_oracle",
]
