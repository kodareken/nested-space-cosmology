# The Matryoshka comparison

Owning note for
[nsc-matryoshka-comparison-v1.json](../results/development/nsc-matryoshka-comparison-v1.json).
The bounded, enumerative side of the "so similar to ours" question about
the cosmohedron.

## What was compared

The cosmohedron theorem
([arXiv:2603.03425](https://arxiv.org/abs/2603.03425)) proves the faces of
the polytope behind the Tr(phi^3) cosmological wavefunction are
Matryoshkas: polygons sequentially wrapped into larger polygons. Its facet
count is the Hipparchus–Schroeder sequence and its maximal Matryoshkas are
counted by OEIS A177384, `1, 2, 10, 72, ...`, with the D-algebraic
generating function `M = x^2 + M M'`.

The NSC regional nesting is a chain of concentric spherical shells:
region inside child inside parent inside ambient. Wrapping any contiguous
group of adjacent shells into the next doll, in containment order, is a
bracketing of a path, counted by the small-Schroeder numbers
`1, 1, 3, 11, 45, 197, ...`.

Both computations reproduced their published prefixes exactly.

## The measured relationship

Both families are nested wrappings of contiguous blocks. The NSC
spherical chain is the small-Schroeder (path bracketing) skeleton; the
cosmohedron polygon Matryoshka is that skeleton raised by the
triangulation freedom of the polygon, which is why its count grows
factorially (`72, ...`) while the chain stays exponential (`45, 197, ...`).
The comparison record stores both growth-ratio lists so the ordering is
visible, not asserted.

The diverging part is named in the same record: polygons versus spheres,
and therefore an additive perimeter energy weight (one inverse of the sum
of the entering momenta per bracket) versus the multiplicative scale
inheritance of the recorded NSC response,

    Gamma_p(x) = K_p(x) - (1/Omega) b Gamma_c(x/Omega)^{-1} b^dagger,

which contributes one inverse factor and one `1/Omega` argument rescale
per nesting level.

## What this is not

It computes no wavefunction, no amplitude, and no measurement of the
actual universe. `no_cosmohedron_equivalence_identified` stays true: the
earlier geometric-vertices result already measured that the NSC
interaction carries no inserted cubic vertex, and this note finds the
enumerative skeleton comparable without claiming the two objects are the
same.

## Owners

Module [`nsc_matryoshka.py`](../src/recursive_horizons/nsc_matryoshka.py),
driver [`derive_nsc_matryoshka.py`](../scripts/derive_nsc_matryoshka.py),
tests [`test_nsc_matryoshka.py`](../tests/test_nsc_matryoshka.py).