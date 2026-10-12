"""Second derivatives of the fixed flat spherical reference connection.

REF1's two-jet evaluator deliberately owns only ``bar_Gamma`` and its first
derivative.  The later CON2 metric-derived divergence needs one additional
derivative.  Keeping it in this downstream module prevents that later proof
surface from changing the hash-bound REF1 data contract.
"""

from __future__ import annotations

from fractions import Fraction

from .reference_connection import (
    ReferenceConnection,
    flat_spherical_annulus_connection,
)


Q = Fraction
N = 4


def _fraction(value: int | Fraction) -> Fraction:
    if isinstance(value, bool):
        raise ValueError("reference-connection values must be rational")
    return value if isinstance(value, Fraction) else Q(value)


def flat_spherical_annulus_connection_second_derivative(
    reference: ReferenceConnection,
    *,
    coordinate_radius: int | Fraction,
) -> tuple[tuple[tuple[tuple[tuple[Fraction, ...], ...], ...], ...], ...]:
    """Return exact ``partial_a partial_b bar_Gamma^c_de`` at the equator.

    The call to the original REF1 evaluator preserves its reference identity
    and strict-annulus validation.  Index order is
    ``[first derivative][second derivative][upper][lower one][lower two]``.
    """

    radius = _fraction(coordinate_radius)
    flat_spherical_annulus_connection(reference, coordinate_radius=radius)
    result = [
        [
            [
                [[Q(0) for _ in range(N)] for _ in range(N)]
                for _ in range(N)
            ]
            for _ in range(N)
        ]
        for _ in range(N)
    ]
    _, r, theta, phi = range(N)

    # Complete nonzero table at theta=pi/2.  The t/r subtable contains only
    # d_r d_r Gamma^theta_(r theta) and its phi counterpart.  The angular
    # entry follows from d_theta^2[-r sin(theta)^2]=2r at the equator.
    result[r][r][theta][r][theta] = result[r][r][theta][theta][r] = Q(2) / radius**3
    result[r][r][phi][r][phi] = result[r][r][phi][phi][r] = Q(2) / radius**3
    result[theta][theta][r][phi][phi] = Q(2) * radius

    return tuple(
        tuple(
            tuple(
                tuple(tuple(component for component in lower) for lower in upper)
                for upper in second_direction
            )
            for second_direction in first_direction
        )
        for first_direction in result
    )
