"""The Matryoshka comparison: nesting counts and the recursive weights.

The cosmohedron theorem (arXiv:2603.03425) proves that the faces of the
polytope underlying the Tr(phi^3) cosmological wavefunction are
Matryoshkas: polygons sequentially wrapped into larger polygons. This
module holds the bounded, enumerative side of the "so similar to ours"
question for Nested-Space Cosmology.

* The cosmohedron facet count is the Hipparchus-Schroeder sequence, and
  the count of maximal Matryoshkas is OEIS A177384, with the
  D-algebraic generating function M = x^2 + M M' (1, 2, 10, 72, ...).
* The NSC regional nesting is a chain of concentric spherical shells:
  region, child, parent, ambient. Wrap any contiguous group of adjacent
  shells into the next doll, in order of containment. The count of
  ordered bracketings of a path is the small-Schroeder sequence
  (1, 1, 3, 11, 45, 197, ...).

So both structures are "nested wrappings of contiguous blocks", but the
polygon Matryoshka carries, per triangulation, the additional bracket
freedom that makes its count factorial (72, ...), while the spherical
chain is the Schroeder skeleton (45, 197, ...). That is the measured
relationship: the NSC sphere is the chain special case of the family the
cosmohedron enumerates, not the polygon case itself. One precise
structural claim, with the diverging part named.

The recursive weights are compared second. The recorded NSC boundary
response is

    Gamma_p(x) = K_p(x) - (1/Omega) b Gamma_c(x/Omega)^{-1} b^dagger,

so every nesting level adds exactly one inverse factor and rescales the
energy by 1/Omega. The cosmohedron weights each bracket subgraph by the
inverse total energy entering it, which is additive in the momenta of
the enclosing polygon. The skeleton (one inverse factor per enclosing
doll) matches; the weight type differs: multiplicative scale inheritance
(1/Omega per level) on the sphere side, additive perimeter energies on
the polygon side. Both statements are recorded, neither is a claim that
the two objects are the same.

Scope. Combinatorics and recursion structure only. It computes no
wavefunction, no amplitude, and no measurement of the actual universe,
and it identifies no cosmohedron equivalence of the NSC interaction.
"""
from __future__ import annotations

import math
from typing import Dict, List, Sequence

API_VERSION = "nsc-matryoshka-v1"

#: Values read from the paper's table in the introduction.
MATRYOSHKA_PREFIX = (1, 2, 10, 72)
#: Small Schroeder numbers, the ordered-bracketing count of a path.
SCHRODER_PREFIX = (1, 1, 3, 11, 45, 197)


def matryoshka_numbers(count: int) -> List[int]:
    """Coefficients of M = x^2 + M M' up to x^(count+1).

    Coefficient balance: m_k = delta_{k,2} + sum over a+b=k+1 of
    b m_a m_b, since M M' contributes b m_a m_b to x^{a+b-1} and every
    admissible index starts at 2.
    """
    if count < 2:
        raise ValueError("the cosmohedron sequence starts at x^2")
    value = [0] * (count + 2)
    value[2] = 1
    for k in range(3, count + 2):
        total = 0
        for a in range(2, k):
            b = k + 1 - a
            if b < 2:
                continue
            total += b * value[a] * value[b]
        value[k] = total
    return value[2:]


def schroder_numbers(count: int) -> List[int]:
    """Small Schroeder numbers with the standard exact recurrence.

    s_0 = s_1 = 1 and s_n = ((6n-3) s_{n-1} - (n-2) s_{n-2}) / (n+1).
    These count complete ordered bracketings of a path of n+1 leaves,
    which is the count of NSC contiguous-shell wrappings of the chain.
    """
    if count < 1:
        raise ValueError("a chain needs at least one shell")
    value = [1, 1]
    for n in range(2, count):
        numerator = (6 * n - 3) * value[n - 1] - (n - 2) * value[n - 2]
        assert numerator % (n + 1) == 0
        value.append(numerator // (n + 1))
    return value[:count]


def growth_ratio(sequence: Sequence[int]) -> List[float]:
    """Consecutive ratios; a divergent ratio indicates factorial growth."""
    return [sequence[i] / sequence[i - 1] for i in range(1, len(sequence))]


def resolvent_levels(omega: float, depth: int) -> List[Dict[str, object]]:
    """Structure of the recorded NSC recursive boundary response.

    Gamma(x) = K(x) - (1/Omega) b Gamma(x/Omega)^{-1} b^dagger. Only the
    structural content is extracted: each level contributes exactly one
    inverse factor and rescales the argument by one more 1/Omega. The
    polygon-side rule contributes, per wrapped sub-polygon, one factor
    that is the inverse of the *sum* of the entering energies.
    """
    if depth < 1:
        return []
    if not math.isfinite(omega) or omega <= 0.0:
        raise ValueError("omega must be a positive finite scale")
    levels = []
    for level in range(1, depth + 1):
        levels.append({
            "level": level,
            "argument_rescale": 1.0 / (omega ** level),
            "own_inverse_factors": 1,
            "coefficient_scale_per_link": 1.0 / omega,
            "weight_type": "multiplicative_scale_inheritance",
        })
    return levels


def comparison() -> Dict[str, object]:
    """The enumerative and structural comparison, with the diverging part named."""
    matryoshka = matryoshka_numbers(7)
    schroder = schroder_numbers(len(SCHRODER_PREFIX))
    return {
        "api_version": API_VERSION,
        "matryoshka_prefix": MATRYOSHKA_PREFIX,
        "matryoshka_computed": matryoshka[: len(MATRYOSHKA_PREFIX)],
        "matryoshka_prefix_matches_paper": tuple(matryoshka[: len(MATRYOSHKA_PREFIX)]) == MATRYOSHKA_PREFIX,
        "schroder_prefix": SCHRODER_PREFIX,
        "schroder_computed": schroder,
        "schroder_prefix_matches_reference": tuple(schroder) == SCHRODER_PREFIX,
        "chain_count_growth_ratios": growth_ratio(schroder),
        "matryoshka_growth_ratios": growth_ratio(matryoshka),
        "structural_claim": (
            "both families are nested wrappings of contiguous blocks; the NSC "
            "spherical chain is the small-Schroeder (path bracketing) skeleton, "
            "and the cosmohedron polygon Matryoshka is that skeleton raised by "
            "the triangulation freedom, which is why its count is factorial"
        ),
        "diverging_part": (
            "polygons versus spheres, and therefore an additive perimeter "
            "energy weight versus the multiplicative scale inheritance of the "
            "recorded NSC response"
        ),
        "resolvent_levels_omega_4": resolvent_levels(4.0, 3),
        "weight_comparison": {
            "nsc_per_level": "one inverse factor, argument rescaled by 1/Omega",
            "cosmohedron_per_bracket": "inverse of the sum of the entering energies",
            "skeleton_matches": True,
            "weight_type_matches": False,
            "no_cosmohedron_equivalence_identified": True,
        },
        "scope": (
            "combinatorics and recursion structure only; no wavefunction, no "
            "amplitude, no measurement of the actual universe"
        ),
    }