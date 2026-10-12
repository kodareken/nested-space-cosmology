"""Exact reference-connection data for the spherical modified-harmonic gauge.

The first full reference-gauge MHG gate works on a positive-radius annulus.  Its fixed
reference geometry is Minkowski spacetime in standard spherical coordinates,

``ds_bar^2 = -dt^2 + dr^2 + r^2 dOmega^2``.

Only the connection and its first coordinate derivative at ``theta=pi/2`` are
needed by the exact local two-jet evaluator.  Keeping this object separate from
the evolving areal radius prevents the reference geometry from silently
becoming a field-dependent gauge source.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


Q = Fraction
N = 4

REFERENCE_COORDINATE_ORDER = ("t", "r", "theta", "phi")
SPHERICAL_FLAT_REFERENCE_ID = "flat_spherical_annulus"


def _q(value: int | Fraction) -> Fraction:
    if isinstance(value, bool):
        raise ValueError("reference-connection values must be rational")
    return value if isinstance(value, Fraction) else Q(value)


def _zero_connection() -> list[list[list[Fraction]]]:
    return [
        [[Q(0) for _ in range(N)] for _ in range(N)]
        for _ in range(N)
    ]


def _zero_connection_derivative() -> list[list[list[list[Fraction]]]]:
    return [
        [
            [[Q(0) for _ in range(N)] for _ in range(N)]
            for _ in range(N)
        ]
        for _ in range(N)
    ]


@dataclass(frozen=True, slots=True)
class ReferenceConnection:
    """Frozen reference-connection choice on a declared radial annulus."""

    reference_id: str
    coordinate_order: tuple[str, ...]
    radial_domain_minimum: Fraction
    angular_evaluation: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "radial_domain_minimum", _q(self.radial_domain_minimum)
        )
        if self.reference_id != SPHERICAL_FLAT_REFERENCE_ID:
            raise ValueError("unknown reference connection")
        if self.coordinate_order != REFERENCE_COORDINATE_ORDER:
            raise ValueError("reference coordinate order must be t,r,theta,phi")
        if self.radial_domain_minimum <= 0:
            raise ValueError("reference radial-domain minimum must be positive")
        if self.angular_evaluation != "theta=pi/2":
            raise ValueError("reference angular evaluation must be theta=pi/2")


@dataclass(frozen=True, slots=True)
class ReferenceConnectionData:
    """Exact ``bar_Gamma`` and ``partial bar_Gamma`` at one annulus point."""

    reference_id: str
    coordinate_radius: Fraction
    christoffel: tuple[tuple[tuple[Fraction, ...], ...], ...]
    christoffel_derivative: tuple[
        tuple[tuple[tuple[Fraction, ...], ...], ...], ...
    ]


def flat_spherical_annulus_reference(
    *, radial_domain_minimum: int | Fraction
) -> ReferenceConnection:
    """Declare the fixed flat spherical reference on ``r>r_min>0``."""

    return ReferenceConnection(
        reference_id=SPHERICAL_FLAT_REFERENCE_ID,
        coordinate_order=REFERENCE_COORDINATE_ORDER,
        radial_domain_minimum=_q(radial_domain_minimum),
        angular_evaluation="theta=pi/2",
    )


def flat_spherical_annulus_connection(
    reference: ReferenceConnection,
    *,
    coordinate_radius: int | Fraction,
) -> ReferenceConnectionData:
    """Evaluate the fixed flat spherical connection and its first derivative.

    Indices are ``Gamma[upper][lower_1][lower_2]`` and
    ``dGamma[derivative][upper][lower_1][lower_2]`` in the frozen coordinate
    order.  The angular derivatives retained below are essential: although
    the intrinsic two-sphere Christoffels vanish at the equator, their
    ``theta`` derivatives do not.
    """

    radius = _q(coordinate_radius)
    if radius <= reference.radial_domain_minimum:
        raise ValueError("coordinate radius must lie strictly inside the annulus")

    gamma = _zero_connection()
    derivative = _zero_connection_derivative()

    t, r, theta, phi = range(N)
    del t  # documents the frozen index order without a magic offset below

    gamma[r][theta][theta] = -radius
    gamma[r][phi][phi] = -radius
    gamma[theta][r][theta] = gamma[theta][theta][r] = 1 / radius
    gamma[phi][r][phi] = gamma[phi][phi][r] = 1 / radius

    derivative[r][r][theta][theta] = -1
    derivative[r][r][phi][phi] = -1
    derivative[r][theta][r][theta] = derivative[r][theta][theta][r] = (
        -1 / radius**2
    )
    derivative[r][phi][r][phi] = derivative[r][phi][phi][r] = (
        -1 / radius**2
    )

    # d_theta(-sin(theta) cos(theta))=1 and d_theta(cot(theta))=-1
    # at theta=pi/2.  The corresponding Christoffels themselves vanish there.
    derivative[theta][theta][phi][phi] = 1
    derivative[theta][phi][theta][phi] = -1
    derivative[theta][phi][phi][theta] = -1

    return ReferenceConnectionData(
        reference_id=reference.reference_id,
        coordinate_radius=radius,
        christoffel=tuple(
            tuple(tuple(component for component in lower) for lower in upper)
            for upper in gamma
        ),
        christoffel_derivative=tuple(
            tuple(
                tuple(tuple(component for component in lower) for lower in upper)
                for upper in direction
            )
            for direction in derivative
        ),
    )
