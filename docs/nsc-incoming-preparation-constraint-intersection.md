# Prepared first-order corrections and the inherited current

On the compatible two-function family's exact subfamily `w(z)=0`,

\[
\mathcal E_N=A_0U+S_N,\qquad
\mathcal E_\beta=dU'+S_\beta,
\]
so

\[
\boxed{\mathcal E_\beta-\frac d{A_0}\partial_z\mathcal E_N
=S_\beta<-\frac{5887634380411363}{265464355000000000}<0.}
\]

Both included constraints therefore cannot vanish on an open axial interval
in this subfamily. This conclusion requires no numerical value for the
finite constant `S_N` and no choice of interval length or endpoint data.

## Exact inputs and normalization

The [compatible surface equations](nsc-incoming-compatible-surface.md) are

\[
\mathcal E_N=A(w)U+B(w)w''+C(w')^2+D(w),\qquad
\mathcal E_\beta=dU'+ew'''+F(w)w'+S_\beta.
\]

Here `D(0)=S_N`, and the unchanged canonical covariance and intrinsic
geometry make both source constants independent of axial position.
The [principal certificate](nsc-incoming-surface-principal.md) encloses

\[
A_0=A(0)\in[0.057399277586276304,0.05739927758627644],
\quad
d\in[-0.029705191751120807,-0.029705191751120727].
\]

In particular `A_0` is nonzero. The displayed elimination uses the exact
coefficients, not a rounded ratio of interval midpoints. The
[momentum balance](nsc-incoming-surface-momentum-balance.md) identifies
`S_beta=P_K<-M`, with the exact rational `M` in the boxed relation. This is
the raw action normalization, not the earlier normalized shift miss.

With `w` identically zero, all its spatial derivatives vanish, and the
polynomial equations reduce exactly to the first pair. Differentiating
the lapse equation and subtracting proves the boxed identity. Equivalently,
lapse cancellation fixes the constant `U=-S_N/A_0`, leaving
`E_beta=S_beta<-M`. Allowing a spatially varying U cannot repair that
contradiction while retaining lapse cancellation at every point.

## Connection to actual preparation

The [physical fixed-transfer theorem](nsc-incoming-fixed-transfer.md)
establishes `delta a_T(z)=delta r_T(z)=0` for first-order history tangents
that preserve the complete incoming C0 in its stated class: unchanged
upstream/source preparation, fixed intrinsic data, smooth compact axial
variations, and the existing unbounded energy continuum.

In the two-function ansatz `delta r=s*w(z)+s^3*U(z)/6`, where
`s=T(rho)-T(1)`, this sets the tangent function `w` to zero. The linearized
constraint changes then obey

\[
\delta\mathcal E_N=A_0\delta U,\qquad
\delta\mathcal E_\beta=d\,\partial_z\delta U,
\quad
\delta\mathcal E_\beta-\frac d{A_0}
\partial_z\delta\mathcal E_N=0.
\]

The baseline residual `(S_N,S_beta)` cannot be canceled by such a
first-order correction: its shift component has the strictly nonzero
combination in the box. This restricts the linearized constraint repair,
not just a chosen positive bump or a finite list of Fourier pairs.

The exact exclusion of the `w=0` subfamily and the first-order preparation
restriction are separate statements. The preparation theorem does not
force `w=0` in an arbitrary finite-amplitude history. The
[local compatible existence theorem](nsc-incoming-local-surface-existence.md)
for freely assigned fixed-C0 surface functions remains mathematically
valid; this calculation shows why its parent preparation is an additional
physical condition.

## Included functional when every permitted normal tangent vanishes

The [higher-order preparation theorem](nsc-incoming-jet-rigidity.md)
extends the same physical argument through normal order four. For its
smooth compact histories with lapse/shift jets fixed as functions along
Sigma, all20 permitted a/r normal/mixed tangent slots vanish. Consequently,

\[
\delta\mathcal E_N=\delta\mathcal E_\beta=0
\]

for the included incoming functional. The
[joint assembly](../src/recursive_horizons/nsc_incoming_joint_constraints.py)
is the sum of the matter insertion, local action gradient and reference
action change. Its matter insertion depends on the same canonical C0 and
fixed instantaneous metric/frame; its derivative is zero. The
[local Euler owner](../src/recursive_horizons/nsc_incoming_local_constraints.py)
uses only the raw metric germ through total order four. The
[spatial reference owner](../src/recursive_horizons/nsc_spatial_reference_symbol.py)
likewise uses that germ through its retained order four. Every allowed
argument of these two terms has zero tangent, so their derivatives vanish.

This implication uses no source-accuracy tolerance and requires no numerical
Jacobian scan. It excludes an affine first-order cancellation of the inherited
nonzero shift residual in that tangent domain. It does not establish
nonlinear injectivity of the preparation map or eliminate changes starting
at higher powers of a finite history parameter. The
[declared full history variation](nsc-declared-action-scope.md) still has
its own transmitting and endpoint questions; no independent remaining
force is introduced into this included constraint assembly.

## Reuse and verification

The source receipts
`nsc-incoming-surface-principal.json`, `nsc-incoming-constraint-gate.json`
and `nsc-incoming-fixed-transfer.json` passed their declared owner/input
hash and nested payload authentication during this integration. No source
quadrature, mode propagation, coefficient derivation or old producer was
rerun. The new step is only the exact substitution and derivative above.

Finite-amplitude preparation and the larger transmitting history problem
remain open. No source occupations, couplings, boundary data, physical
initial tuple, action terms or metric evolution are introduced by this
intersection.
