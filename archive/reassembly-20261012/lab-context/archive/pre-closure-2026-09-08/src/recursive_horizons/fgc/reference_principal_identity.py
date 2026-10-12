"""Formal index proof of the REF1-to-MHG gauge principal identity.

The reference connection enters the modified-harmonic gauge vector without
metric derivatives.  It is therefore lower order in the metric principal
linearization.  This module proves the remaining index contraction as an
identity in a free commutative polynomial algebra.  No spacetime point,
fixture, random substitution, floating point, or interval-overlap argument is
used.

The convention is the one implemented by REF1,

``C^a = - tilde(g)^(rho sigma) (Gamma^a_(rho sigma)-barGamma^a_(rho sigma))``.

At principal order, varying the physical connection and differentiating once
more gives the same projector contraction used by MHG1.  Coefficient
variations of ``F``, the three inverse metrics, and index-lowering metrics are
lower order and consequently remain formal background coefficients here.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from numbers import Integral
from typing import Any, Iterable


Q = Fraction
N = 4
Monomial = tuple[str, ...]


class _Formal:
    """A small canonical sparse polynomial over exact rational coefficients."""

    __slots__ = ("terms",)

    def __init__(self, terms: Iterable[tuple[Monomial, Fraction]] = ()) -> None:
        combined: dict[Monomial, Fraction] = {}
        for monomial, coefficient in terms:
            if isinstance(coefficient, bool) or not isinstance(coefficient, Fraction):
                raise TypeError("formal coefficients must be Fractions")
            key = tuple(sorted(monomial))
            combined[key] = combined.get(key, Q(0)) + coefficient
        self.terms = tuple(
            (monomial, coefficient)
            for monomial, coefficient in sorted(combined.items())
            if coefficient
        )

    @classmethod
    def constant(cls, value: Fraction | Integral) -> "_Formal":
        if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
            raise TypeError("formal constants must be exact rationals")
        coefficient = value if isinstance(value, Fraction) else Q(int(value))
        return cls(()) if coefficient == 0 else cls((((), coefficient),))

    @classmethod
    def atom(cls, name: str) -> "_Formal":
        if not isinstance(name, str) or not name:
            raise ValueError("formal atom name must be nonempty")
        return cls((((name,), Q(1)),))

    @staticmethod
    def coerce(value: object) -> "_Formal" | Any:
        if isinstance(value, _Formal):
            return value
        if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
            return NotImplemented
        return _Formal.constant(value)

    def __add__(self, other: object) -> "_Formal" | Any:
        value = self.coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return _Formal(self.terms + value.terms)

    __radd__ = __add__

    def __neg__(self) -> "_Formal":
        return _Formal((monomial, -coefficient) for monomial, coefficient in self.terms)

    def __sub__(self, other: object) -> "_Formal" | Any:
        value = self.coerce(other)
        return NotImplemented if value is NotImplemented else self + (-value)

    def __rsub__(self, other: object) -> "_Formal" | Any:
        value = self.coerce(other)
        return NotImplemented if value is NotImplemented else value - self

    def __mul__(self, other: object) -> "_Formal" | Any:
        value = self.coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return _Formal(
            (left + right, a * b)
            for left, a in self.terms
            for right, b in value.terms
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "_Formal" | Any:
        if isinstance(other, bool) or not isinstance(other, (Fraction, Integral)):
            return NotImplemented
        denominator = other if isinstance(other, Fraction) else Q(int(other))
        if denominator == 0:
            raise ZeroDivisionError("formal rational scale is zero")
        return _Formal(
            (monomial, coefficient / denominator)
            for monomial, coefficient in self.terms
        )

    def __eq__(self, other: object) -> bool:
        value = self.coerce(other)
        return False if value is NotImplemented else self.terms == value.terms

    def canonical(self) -> str:
        return ";".join(
            f"{coefficient.numerator}/{coefficient.denominator}:"
            + ("*".join(monomial) if monomial else "1")
            for monomial, coefficient in self.terms
        )


ZERO = _Formal.constant(0)
ONE = _Formal.constant(1)


def _symmetric_atoms(prefix: str) -> tuple[tuple[_Formal, ...], ...]:
    return tuple(
        tuple(_Formal.atom(f"{prefix}{min(row, column)}{max(row, column)}") for column in range(N))
        for row in range(N)
    )


def _projector(
    inverse: tuple[tuple[_Formal, ...], ...],
    alpha: int,
    beta: int,
    mu: int,
    nu: int,
) -> _Formal:
    """Formal ``P_alpha^(beta mu nu)`` in the repository convention."""

    return (
        int(alpha == mu) * inverse[nu][beta]
        + int(alpha == nu) * inverse[mu][beta]
        - int(alpha == beta) * inverse[mu][nu]
    ) / 2


def _formal_routes(
    *, connection_difference_sign: Fraction
) -> tuple[tuple[tuple[_Formal, ...], ...], tuple[tuple[_Formal, ...], ...]]:
    """Construct the connection-variation and projector routes independently."""

    if isinstance(connection_difference_sign, bool) or not isinstance(
        connection_difference_sign, Fraction
    ):
        raise TypeError("connection_difference_sign must be a Fraction")

    physical = _symmetric_atoms("gU")
    tilde = _symmetric_atoms("tU")
    hat = _symmetric_atoms("hU")
    perturbation = _symmetric_atoms("k")
    xi = tuple(_Formal.atom(f"xi{index}") for index in range(N))
    effective_planck = _Formal.atom("F")

    # Route A starts from the universal principal variation of the physical
    # Levi-Civita connection.  barGamma and coefficient variations carry no
    # second derivative of the perturbation and therefore do not enter.
    delta_gamma = tuple(
        tuple(
            tuple(
                sum(
                    (
                        physical[alpha][lam]
                        * (
                            xi[rho] * perturbation[lam][sigma]
                            + xi[sigma] * perturbation[lam][rho]
                            - xi[lam] * perturbation[rho][sigma]
                        )
                        / 2
                        for lam in range(N)
                    ),
                    ZERO,
                )
                for sigma in range(N)
            )
            for rho in range(N)
        )
        for alpha in range(N)
    )
    delta_constraint = tuple(
        connection_difference_sign
        * sum(
            (
                tilde[rho][sigma] * delta_gamma[alpha][rho][sigma]
                for rho in range(N)
                for sigma in range(N)
            ),
            ZERO,
        )
        for alpha in range(N)
    )
    connection_route = tuple(
        tuple(
            effective_planck
            * sum(
                (
                    _projector(hat, alpha, beta, mu, nu)
                    * xi[beta]
                    * delta_constraint[alpha]
                    for alpha in range(N)
                    for beta in range(N)
                ),
                ZERO,
            )
            for nu in range(N)
        )
        for mu in range(N)
    )

    # Route B is the MHG1 projector formula, constructed without reusing the
    # connection variation above.
    hat_xi = tuple(
        tuple(
            tuple(
                sum(
                    (xi[gamma] * _projector(hat, alpha, gamma, mu, nu) for gamma in range(N)),
                    ZERO,
                )
                for nu in range(N)
            )
            for mu in range(N)
        )
        for alpha in range(N)
    )
    tilde_xi = tuple(
        tuple(
            tuple(
                sum(
                    (xi[delta] * _projector(tilde, beta, delta, rho, sigma) for delta in range(N)),
                    ZERO,
                )
                for sigma in range(N)
            )
            for rho in range(N)
        )
        for beta in range(N)
    )
    projector_route = tuple(
        tuple(
            -effective_planck
            * sum(
                (
                    physical[alpha][beta]
                    * hat_xi[alpha][mu][nu]
                    * tilde_xi[beta][rho][sigma]
                    * perturbation[rho][sigma]
                    for alpha in range(N)
                    for beta in range(N)
                    for rho in range(N)
                    for sigma in range(N)
                ),
                ZERO,
            )
            for nu in range(N)
        )
        for mu in range(N)
    )
    return connection_route, projector_route


def reference_principal_identity_certificate(
    *, connection_difference_sign: Fraction = Q(-1)
) -> dict[str, Any]:
    """Prove the universal REF1 gauge-principal identity and fail closed.

    ``connection_difference_sign`` exists only as an exact mutation probe.  A
    sign inconsistent with REF1 must make the independently assembled routes
    disagree and is rejected.
    """

    connection, projector = _formal_routes(
        connection_difference_sign=connection_difference_sign
    )
    differences = tuple(
        tuple(connection[mu][nu] - projector[mu][nu] for nu in range(N))
        for mu in range(N)
    )
    exact = all(difference == ZERO for row in differences for difference in row)
    if not exact:
        raise ValueError("formal REF1 connection route differs from the MHG projector route")
    connection_text = "|".join(
        connection[mu][nu].canonical() for mu in range(N) for nu in range(N)
    )
    projector_text = "|".join(
        projector[mu][nu].canonical() for mu in range(N) for nu in range(N)
    )
    if connection_text != projector_text:
        raise AssertionError("equal formal tensors have different canonical encodings")
    return {
        "classification": "universal_formal_REF1_gauge_principal_identity_not_a_point_or_sample_test",
        "dimension": N,
        "gauge_vector_convention": "C^a=-tilde_g^(rho sigma)(Gamma^a_rho_sigma-barGamma^a_rho_sigma)",
        "reference_connection_absent_from_principal_variation": True,
        "coefficient_variations_are_lower_order": True,
        "formal_background_atoms_are_algebraically_independent": True,
        "connection_variation_route_equals_projector_route": True,
        "all_sixteen_contravariant_tensor_differences_zero": True,
        "connection_route_term_counts": tuple(
            tuple(len(connection[mu][nu].terms) for nu in range(N)) for mu in range(N)
        ),
        "canonical_identity_sha256": sha256(connection_text.encode("ascii")).hexdigest(),
        "proof_method": "exact_sparse_commutative_polynomial_index_contraction_over_Q",
        "uses_fixtures_sampling_intervals_or_floating_point": False,
    }
