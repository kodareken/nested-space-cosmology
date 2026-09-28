# FGC-1-HYP1-CON2-MPROP1: local metric-derived subsidiary identities

MPROP1 packages two exact local differential identities in the frozen REF1
spherical annulus chart.  The first derives the metric-defined modified-
harmonic gauge vector through its second coordinate partials from a physical
third jet, including the second derivative of the fixed flat spherical
reference connection.  It verifies

\[
L_{\hat g}(C[g])=\nabla\!\cdot E_{\rm REF1}+E_\phi\nabla\phi+E_\chi\nabla\chi
\]

in the repository's ACT1 Noether-sign convention.  The second is the complete
one-spatial-dimensional kinematic reduction identity

\[
\partial_t C+\partial_r D-K=0,
\]

for \(D=\partial_tu-p\), \(C=q-\partial_ru\), and
\(K=\partial_tq-\partial_rp\).  The controls are a flat zero-third-jet
reference and the activated COMP1 two-jet with deterministic exact third
derivatives.

The flat spherical reference connection really has nonzero second coordinate
derivatives at the equatorial evaluation point.  They are evaluated and passed
through the metric-derived construction.  For the strictly spherical
\(C^t,C^r\) sector used here, their contraction vanishes identically because
the allowed inverse-metric and connection index patterns do not overlap.  The
zero contribution is therefore a proved symmetry cancellation, not an omitted
reference-connection term.

This is a local differential identity, not a Cauchy uniqueness theorem, a
constraint-compatible hypersurface, physical Einstein-constraint propagation,
a boundary map/Kreiss estimate, an IBVP, evolution, collapse, transition, or
singularity-resolution result.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_hyp1_con2_mprop1.py
```
