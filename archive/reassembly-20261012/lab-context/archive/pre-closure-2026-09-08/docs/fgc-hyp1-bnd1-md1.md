# FGC-1-HYP1-BND1-MD1: frozen radial maximal dissipation

**FGC-1-HYP1-BND1-MD1** composes the exact compact radial UHYP1 frame with
the complete one-dimensional kinematic reduction identity in CON2-MPROP1. It
constructs a homogeneous modal boundary operator on a finite spherical
annulus and proves a uniform frozen-coefficient energy-flux estimate for the
eighteen-variable first-order principal system.

Write the principal variables as

\[
U=(u,W),\qquad W=(p,q),
\]

where the six value variables \(u\) have zero radial principal speed and
UHYP1 supplies twelve propagating modes

\[
A=V\Lambda V^{-1},\qquad H=V^{-T}V^{-1}.
\]

The full principal matrices are

\[
A_{\rm full}=\operatorname{diag}(0_6,A),\qquad
H_{\rm full}=\operatorname{diag}(I_6,H).
\]

In modal coordinates \(y=V^{-1}W\), the boundary flux is exactly

\[
W^THAW=y^T\Lambda y=\sum_j c_j y_j^2.
\]

The outer normal is \(+\partial_r\), while the inner normal is
\(-\partial_r\). Exact interval speed enclosures—not sampled floating-point
speeds—select six incoming propagating modes at each end. The boundary
operator sets precisely those incoming amplitudes to zero and no others:

\[
B_\partial W=\Pi^{\partial}_{\rm in}V^{-1}W=0.
\]

Both selectors have exact rank six. The normalized modal Lopatinski matrix is
the identity, the outward flux is nonnegative at both ends, and the
homogeneous frozen principal energy is nonincreasing. Thus the propagating
block has a uniformly maximally dissipative boundary subspace throughout the
same compact REF1 branch graph certified by UHYP1.

CON2-MPROP1 independently proves, for each of the six fields,

\[
D=\partial_tu-p,\qquad C=q-\partial_ru,\qquad
K=\partial_tq-\partial_rp,
\]

and

\[
\partial_tC+\partial_rD-K=0.
\]

On the differentiable kinematic shell \(D=K=0\), all six \(C\) fields have
zero characteristic speed and no incoming boundary component in one spatial
dimension. BND1-MD1 therefore imposes no independent reduction-constraint
boundary data. This is a compatibility statement for the complete
**kinematic** reduction subsystem only.

## Exact claim boundary

This artifact proves uniform frozen radial **main-system** maximal
dissipation plus compatibility with the zero-speed kinematic reduction
subsystem. It does not derive the nonlinear lower-order source evaluator, a
quasilinear existence theorem, the incoming metric-derived gauge or physical
Einstein-constraint fields, a constraint-preserving ACT1 boundary map,
corner conditions, a regular centre, multidirectional boundary stability, an
EFT-valid evolution, collapse, affine defocusing, a finite transition
surface, or singularity resolution.

In particular, a divergence at a proposed transition cannot be called a
physical wall on this evidence. A later transition gate must replace it by a
finite invariant surface carrying a closed, constraint-compatible handoff.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_hyp1_bnd1_md1.py
```

The construction follows the standard symmetric-hyperbolic incoming-field
and energy-flux strategy used in generalized-harmonic boundary analyses; the
ACT1 constraint boundary map remains new work rather than an inherited
result.
