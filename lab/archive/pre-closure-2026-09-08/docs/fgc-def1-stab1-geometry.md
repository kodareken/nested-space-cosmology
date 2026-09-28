# DEF1-STAB1 geometric sensitivity instrument

This outcome-blind subcomponent of the accepted PLAN.md Wave-0 error map
implements local geometry-to-Q propagation. It is not a trajectory reader,
an admission policy, or a completed DEF1-STAB1/DEF1-PREF1 certificate.
The [DEF1-STAB1 owner](fgc-def1-stab1.md) defines the saved component budget,
margin contracts and remaining pre-holdout qualification boundary.

`src/recursive_horizons/fgc/def1_geometry_error.py` accepts the full two-jets
of ADM lapse, shift, radial metric and areal radius, plus a radial tangent.
All 26 named inputs are explicit. Exact closed intervals must keep lapse,
radial metric, areal radius and the future tangent time component positive.
The existing RED1 warped-product curvature and metric instruments supply
the full four-dimensional Ricci contraction.

```text
theta = 2 k^a partial_a R / R
Rkk = R^(4)_ab k^a k^b
Q = -theta^2/2 - Rkk
```

The shear term is exactly zero for the spherical radial congruence. This is
the smooth Q expression on the input space; its physical Raychaudhuri
interpretation still requires independent nullness and affine-transport
evidence. A box whose nullness interval contains zero is not a nullness
proof. The full Ricci expression is retained for off-null perturbations;
the exact-null shortcut `Rkk=-2 Hess(R)(k,k)/R` would omit
`Ric^(2)(k,k)` there.

For supplied exact nominal inputs x and nonnegative error radii epsilon,
construct the box x +/- epsilon before evaluating Q. Existing exact interval
first tangents enclose every partial derivative throughout that box. The
multivariable mean-value bound is

```text
|Q(true)-Q(nominal)| <= sum_i sup_box |partial Q / partial x_i| * epsilon_i.
```

The returned right-hand side bounds the Q error conditional on the true
inputs lying inside the declared box. The segment from nominal to true
inputs stays in the box; the whole-box Lorentzian/annular checks keep the
expression differentiable there. The implementation does not fit the bound
to Q's value or sign. Correlated components may be conservatively summed.
Each radius uses its corresponding input's coordinate/affine units; protocol
normalization and conversion of dimensionless error norms remain caller-owned.

The inherited generic interval determinant can lose correlations between
ADM-generated metric entries. Consequently a wide shift box can fail this
certificate even though every underlying ADM metric is Lorentzian. Such a
refusal means this enclosure construction is inconclusive; it is not evidence
of physical signature change. No smaller uncertainty box is silently used.

This closes a local calculus dependency, not the input-error problem. The
owners of spatial/temporal, constraints, source, affine, interpolation,
extraction and other errors must supply valid radii and their provenance.
This instrument proves neither a global PDE error nor those input premises.
Its numbers are not self-authenticating records or substitutes for the nine
separate final DEF1 booleans. Pre-holdout qualification of the complete error
map remains a Wave-0 obligation; later trajectory application is separate.

Focused controls compare against the separate direct four-dimensional
curvature route, include Minkowski and contracting de Sitter, retain the
off-null base-Ricci correction, and check exact derivative and injected-error
bounds. Missing inputs, floating/bool aliases, negative radii, a lost metric
domain or arithmetic resource limits refuse. Prospective input/output
rational limits are 4,096/262,144 bits; they are resource limits, not fitted
scientific thresholds. No historical source, config or result is modified.
