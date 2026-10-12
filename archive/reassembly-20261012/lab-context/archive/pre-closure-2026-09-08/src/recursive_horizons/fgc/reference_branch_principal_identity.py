"""Structural REF1-to-MHG branch-principal identity.

This is a definition-level and abstract-block theorem.  It does not establish
that a particular background has an invertible kinetic block; it says what the
complete REF1 branch matrix is whenever the stated regularity premises hold.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from typing import Any, Iterable

from . import modified_harmonic_reference as ref1
from . import spherical_reduction as red1
from . import spherical_symbol as sym1
from . import modified_harmonic as mhg1
from .reference_principal_identity import reference_principal_identity_certificate


Q = Fraction


class _Expr:
    """Tiny exact noncommutative expression algebra with ``Kinv*K=1``."""

    __slots__ = ("terms",)

    def __init__(self, terms: Iterable[tuple[tuple[str, ...], Fraction]] = ()) -> None:
        combined: dict[tuple[str, ...], Fraction] = defaultdict(lambda: Q(0))
        for word, coefficient in terms:
            if not isinstance(coefficient, Fraction):
                raise TypeError("abstract block coefficients must be Fractions")
            reduced: list[str] = []
            for token in word:
                if reduced and reduced[-1] == "Kinv" and token == "K":
                    reduced.pop()
                else:
                    reduced.append(token)
            combined[tuple(reduced)] += coefficient
        self.terms = tuple(sorted((word, coefficient) for word, coefficient in combined.items() if coefficient))

    @classmethod
    def atom(cls, token: str) -> "_Expr":
        return cls((((token,), Q(1)),))

    @classmethod
    def scalar(cls, value: int) -> "_Expr":
        return cls((((), Q(value)),))

    def __add__(self, other: "_Expr") -> "_Expr":
        return _Expr(self.terms + other.terms)

    def __neg__(self) -> "_Expr":
        return _Expr((word, -coefficient) for word, coefficient in self.terms)

    def __sub__(self, other: "_Expr") -> "_Expr":
        return self + (-other)

    def __mul__(self, other: "_Expr") -> "_Expr":
        return _Expr((left + right, a * b) for left, a in self.terms for right, b in other.terms)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _Expr) and self.terms == other.terms


Z = _Expr.scalar(0)
I = _Expr.scalar(1)
K = _Expr.atom("K")
KI = _Expr.atom("Kinv")
L = _Expr.atom("L")
M = _Expr.atom("M")


def _block_multiply(left, right):
    return tuple(
        tuple(sum((left[row][inner] * right[inner][column] for inner in range(2)), Z) for column in range(2))
        for row in range(2)
    )


def _definition_route() -> dict[str, bool]:
    """Fail closed if the live code-definition handoffs change."""
    full_names = set(ref1.modified_harmonic_full_residuals.__code__.co_names)
    principal_block_names = set(red1.principal_block.__code__.co_names)
    radial_names = set(sym1.radial_symbol.__code__.co_names)
    standard_names = set(mhg1.standard_first_order_principal_system.__code__.co_names)
    facts = {
        "ref1_full_rows_call_unredefined_residuals": "residuals" in full_names,
        "red1_second_jet_block_differentiates_residuals": "residuals" in principal_block_names,
        "radial_symbol_is_constructed_from_principal_matrix": "principal_matrix" in radial_names,
        "standard_reduction_uses_matrix_inverse": "matrix_inverse" in standard_names,
    }
    if not all(facts.values()):
        raise ValueError("live REF1/RED1/MHG definition route has changed")
    return facts


def abstract_branch_block_identity(*, mutate_coefficient_sign: bool = False, mutate_lower_sign: bool = False) -> dict[str, Any]:
    """Prove the Schur/reduction identity from ``P(c)=M-cL+c^2K``.

    The sole algebraic premise is the displayed left inverse relation
    ``Kinv*K=I``.  No numerical block values or background fixtures enter.
    """
    # P(c)=M-cL+c^2K means the code's coefficient-one block is -L.
    coefficient_one = L if mutate_coefficient_sign else -L
    a_tr = -coefficient_one
    a_tt, a_rr = K, M
    a_time_inverse = ((KI, Z), (Z, I))
    a_radial = ((a_tr, a_rr), ((I if mutate_lower_sign else -I), Z))
    standard_route = _block_multiply(a_time_inverse, a_radial)
    ref1_route = ((KI * L, KI * M), (-I, Z))
    if standard_route != ref1_route:
        raise ValueError("abstract REF1 branch and standard first-order reduction differ")
    return {
        "polynomial_convention": "P(c)=M-cL+c^2K",
        "kinetic_premise": "K is invertible with declared left inverse Kinv",
        "standard_a_tt": "K",
        "standard_a_tr": "L=-coefficient_one(P)",
        "standard_a_rr": "M",
        "ref1_branch_blocks": ("Kinv*L", "Kinv*M", "-I", "0"),
        "standard_reduction_blocks": ("Kinv*L", "Kinv*M", "-I", "0"),
        "abstract_block_identity_exact": True,
        "uses_numeric_block_fixture": False,
    }


def reference_branch_principal_identity_certificate() -> dict[str, Any]:
    """Compose live row definitions, universal gauge identity, and reduction."""
    definitions = _definition_route()
    gauge = reference_principal_identity_certificate()
    block = abstract_branch_block_identity()
    if not gauge["connection_variation_route_equals_projector_route"]:
        raise ValueError("universal REF1 gauge principal identity is unavailable")
    if not block["abstract_block_identity_exact"]:
        raise ValueError("abstract branch block identity is unavailable")
    return {
        "classification": "universal_structural_REF1_to_MHG_branch_principal_identity_conditional_on_regular_background_premises",
        "live_definition_route": definitions,
        "unredefined_rows": {
            "source": "spherical_reduction.residuals",
            "independent_order": ("metric_tt", "metric_tr", "metric_rr", "metric_theta_theta", "scalar_phi", "scalar_chi"),
            "ref1_uses_unredefined_rows_directly": True,
            "red1_principal_and_SYM1_radial_symbol_share_that_source": True,
        },
        "gauge_extension": {
            "universal_formal_certificate": gauge,
            "scalar_rows_unmodified_by_definition": True,
            "complete_second_order_symbol_equals_unredefined_radial_symbol_plus_MHG_extension": True,
        },
        "branch_reduction": block,
        "regularity_premises": {
            "reference_connection_contributes_only_lower_order_terms": True,
            "effective_Einstein_coefficient_F_strictly_positive": True,
            "coordinate_time_kinetic_block_K_invertible": True,
        },
        "theorem_scope": {
            "universal_structural_identity": True,
            "numerical_background_existence_proven": False,
            "open_domain_kinetic_regularness_proven": False,
            "characteristics_or_hyperbolicity_proven": False,
        },
    }
