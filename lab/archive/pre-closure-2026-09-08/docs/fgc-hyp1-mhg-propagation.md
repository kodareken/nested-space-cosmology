# FGC-1-HYP1-MHG4-PROP1: lower-order reference-gauge operator

**FGC-1-HYP1-MHG4-PROP1** closes one equation-level algebraic prerequisite
downstream of REF1 and IMP1. It constructs the exact lower-order covariant
operator acting on an independently represented, spherically symmetric gauge
vector two-jet,

\[
\boxed{
\mathfrak G^\nu[C]
=
\nabla_\mu\!\left(
F\widehat P_\alpha{}^{\beta\mu\nu}
\nabla_\beta C^\alpha
\right).
}
\]

The result is conditional on the ACT1 Noether identity and the two unmodified
scalar equations. It is not a direct componentwise evaluation of the
divergence of the complete metric-generated REF1 residual: that calculation
would require third metric/scalar jets and second derivatives of the reference
connection, which the current exact two-jet representation does not contain.

PROP1 also does not define a nonlinear first-order radial system or prove its
reduction constraints propagate. That logically separate construction is the
next **FGC-1-HYP1-FO1-RC1** gate.

## Frozen residual convention and Noether handoff

ACT1 uses

\[
S=\int d^4x\sqrt{-g}\left[
\frac12F(\phi)R
-\frac12(\nabla\phi)^2-V(\phi)
-\frac12(\nabla\chi)^2
+f(\phi)\mathcal G
\right],
\]

with metric residual

\[
E_{\mu\nu}
=F G_{\mu\nu}
+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)F
-T^{(\phi)}_{\mu\nu}
-T^{(\chi)}_{\mu\nu}
-T^{(\mathrm{GB})}_{\mu\nu}
\]

and scalar residuals

\[
E_\phi
=\Box\phi-V'(\phi)+\frac12F'(\phi)R+f'(\phi)\mathcal G,
\qquad
E_\chi=\Box\chi.
\]

The variation convention is

\[
\delta S
=\frac12\int\!\sqrt{-g}\,E_{\mu\nu}\,\delta g^{\mu\nu}
+\int\!\sqrt{-g}\left(E_\phi\,\delta\phi+E_\chi\,\delta\chi\right).
\]

Under an infinitesimal diffeomorphism,
`delta g^{mu nu}=-2 nabla^(mu xi^nu)`,
`delta phi=xi^mu nabla_mu phi`, and
`delta chi=xi^mu nabla_mu chi`. Integration by parts therefore fixes the
project's exact Noether sign:

\[
\boxed{
\nabla_\mu E^{\mu\nu}
+E_\phi\nabla^\nu\phi
+E_\chi\nabla^\nu\chi
=0.
}
\]

This sign also follows directly in the minimally coupled scalar control:
`nabla T^(phi)=E_phi nabla phi`, while `E_metric=G-T`.

REF1's complete gauge-fixed metric residual is

\[
\mathcal E^{\mu\nu}
=E^{\mu\nu}
+F\widehat P_\alpha{}^{\beta\mu\nu}\nabla_\beta C^\alpha,
\]

with unchanged scalar equations. Hence

\[
\mathcal E^{\mu\nu}=0
\quad\Longrightarrow\quad
\mathfrak G^\nu[C]
=E_\phi\nabla^\nu\phi+E_\chi\nabla^\nu\chi.
\]

Only after imposing `E_phi=E_chi=0` is the gauge equation homogeneous:

\[
\boxed{\mathfrak G^\nu[C]=0.}
\]

PROP1 serializes this conditional handoff. It does not mark the differential
Noether identity as newly evaluated from third jets; VAR1 supplied its
action-level derivation and exact curvature-contraction checks.

## Tensorial reference gauge

The propagated object is REF1's tensorial vector

\[
C^\alpha
=-\widetilde g^{\rho\sigma}
\left(
\Gamma^\alpha{}_{\rho\sigma}
-\bar\Gamma^\alpha{}_{\rho\sigma}
\right).
\]

The fixed reference connection determines what `C=0` means and constrains
compatible gauge/metric initial data. Once `C` is treated as the vector acted on by
`mathfrak G`, the reference connection is not an independent forcing term in
the homogeneous operator. Its influence is already contained in the
metric-derived value and derivatives of `C`.

PROP1's executable input is instead an explicit independent gauge-vector
two-jet. That is sufficient to certify the operator coefficients without
pretending that a metric two-jet contains the third derivatives needed to
construct `nabla nabla C[g]` directly.

## Exact lower-order expansion

The hat projector is

\[
\widehat P_\alpha{}^{\beta\mu\nu}
=\delta_\alpha^{(\mu}\widehat g^{\nu)\beta}
-\frac12\delta_\alpha^\beta\widehat g^{\mu\nu}.
\]

Writing

\[
X^{\mu\nu}
=F\widehat P_\alpha{}^{\beta\mu\nu}\nabla_\beta C^\alpha,
\]

gives

\[
X^{\mu\nu}
=\frac F2\left(
\widehat g^{\nu\beta}\nabla_\beta C^\mu
+\widehat g^{\mu\beta}\nabla_\beta C^\nu
-\widehat g^{\mu\nu}\nabla_\alpha C^\alpha
\right).
\]

Using the frozen curvature convention

\[
[\nabla_\mu,\nabla_\beta]V^\rho
=R^\rho{}_{\sigma\mu\beta}V^\sigma,
\qquad
R_{\sigma\beta}=R^\rho{}_{\sigma\rho\beta},
\]

the first and third second-derivative terms cancel up to one Ricci
commutator. The complete operator is

\[
\boxed{
\begin{aligned}
\mathfrak G^\nu[C]
={}&
\nabla_\mu\!\left(
F\widehat P_\alpha{}^{\beta\mu\nu}
\right)\nabla_\beta C^\alpha
\\
&+\frac F2\widehat g^{\mu\beta}
\nabla_\mu\nabla_\beta C^\nu
+\frac F2\widehat g^{\nu\beta}R_{\alpha\beta}C^\alpha.
\end{aligned}
}
\]

The Ricci term has the displayed plus sign. A possible-looking term

\[
\widehat g^{\mu\beta}R^\nu{}_{\alpha\mu\beta}C^\alpha
\]

vanishes identically because `hat g^{mu beta}` is symmetric while the Riemann
tensor is antisymmetric in `mu,beta`; it is not retained as a second curvature
source.

Thus the principal part is exactly

\[
\boxed{
\frac F2\widehat g^{\mu\beta}
\nabla_\mu\nabla_\beta C^\nu.
}
\]

The lower-order coefficient tensors are

\[
L^\nu{}_\alpha{}^\beta
=\nabla_\mu\!\left(
F\widehat P_\alpha{}^{\beta\mu\nu}
\right),
\qquad
M^\nu{}_\alpha
=\frac F2\widehat g^{\nu\beta}R_{\alpha\beta}.
\]

`L` includes the complete derivatives of dynamic `F`, the hat auxiliary
inverse metric, and every physical-connection contribution for all projector
indices.

## Independent exact operator routes

PROP1 evaluates the same operator in two ways.

1. The **coordinate-divergence route** first constructs
   `X^{mu nu}=F hat_P nabla C`, differentiates its exact coordinate components,
   and adds the two contravariant-tensor connection terms required by
   `nabla_mu X^{mu nu}`.
2. The **expanded covariant route** constructs `nabla(F hat_P)`,
   `nabla C`, `nabla nabla C`, and the Ricci commutator separately and evaluates
   the boxed wave-plus-lower-order formula above.

Both routes use ordinary exact rational arithmetic and must agree entry by
entry. The coordinate route assembles
`partial_direction(nabla_beta C^alpha)` directly from the raw second partial
jet, `partial Gamma`, the vector value, and its first partial jet; it does not
consume or algebraically reconstruct the expanded route's
`nabla_direction nabla_beta C^alpha`. A mutation regression corrupts only that
expanded covariant-second-derivative intermediate and requires the two routes
to disagree. Neither route imports MHG1's principal matrix as its answer.

The executable also evaluates all four output components on basis gauge-vector
two-jets with independent inputs only in `C^t,C^r`. This is the spherical
`1+1` coefficient form after evaluating the full four-dimensional equatorial
angular covariant derivatives. Angular connection contributions are retained
in the stored value- and first-derivative blocks and are not set to zero merely
because `C^theta=C^phi=0`. The certificate serializes four-output-by-two-input
coefficient arrays in the explicit input order `C^t,C^r`; it serializes MHG1's
full four-by-four principal authority separately. Angular input columns are
absent, not numerical zeros masquerading as evaluated coefficients. Its
physical spherical subsystem publishes every coefficient of

\[
\mathfrak G^A
=A^A{}_{B}{}^{ab}\,\partial_a\partial_b C^B
+B^A{}_{B}{}^a\,\partial_a C^B
+D^A{}_{B}C^B,
\qquad A,B,a,b\in\{t,r\}.
\]

For the symmetric mixed slot `dtr`, the stored coefficient includes both
`(t,r)` and `(r,t)` contributions.

## MHG1 principal regression

Coefficient extraction from the complete PROP1 operator must give

\[
\frac F2\widehat g^{ab}\xi_a\xi_b\,\delta^A{}_B.
\]

For the radial covector `xi_A=(-c,1)`, this is MHG1's complete hat-cone
gauge-error block. PROP1 compares every coefficient—not merely its determinant
or roots—with the existing independently constructed MHG1 principal
contraction.

This is a regression after the full lower-order operator is assembled. MHG1
does not supply PROP1's first- or zeroth-order coefficients.

## Declared exact controls

The certificate uses two background classes.

### Flat inactive FGC-QR root

The IMP1 root retains nonzero FGC-QR action couplings while
`phi=chi=0`, `F=4`, `nabla F=0`, and physical/reference flat spherical
connections agree. A zero gauge-vector two-jet must give

\[
\mathfrak G^\nu[0]=0
\]

through both exact routes. A separate nonzero gauge probe prevents that check
from being vacuous and isolates the flat hat-wave plus spherical-coordinate
connection terms.

### Activated coefficient probe

The existing `FGCQR_activated_generic` background is paired with a frozen
independent nonzero `C^t,C^r` two-jet. It is a coefficient probe, not a
solution and not a claim that the supplied two-jet is the metric-derived REF1
constraint through second order. It must exercise nonzero contributions from

- `nabla F`;
- the differentiated hat projector;
- physical connection and curvature;
- first gauge-vector derivatives; and
- second gauge-vector derivatives.

The certificate verifies that the two operator routes, extracted coefficient
map, linear homogeneity and superposition, and MHG1 principal regression all
remain exact on that background.

## Source boundary

Kovacs and Reall derive the modified-harmonic gauge-error wave mechanism and,
for their one-scalar Horndeski system, use the ungauge-fixed scalar equation and
Noether identity to retain it
([arXiv:2003.08398](https://arxiv.org/abs/2003.08398), especially Eqs. 7--12
and 83--88). East and Ripley display a neighboring nonlinear ESGB
constraint-violation equation with projector-derivative, curvature, and
damping terms
([arXiv:2011.03547](https://arxiv.org/abs/2011.03547), Eqs. 1--4).

Those works fix the structural burden; they do not prove this ACT1 result.
PROP1 derives the dynamic-`F`, independent-`chi`, exact spherical
coefficient map from REF1's own extension. The Noether identity and the
general gauge-wave mechanism are standard; the project-specific value is the
explicit ACT1/REF1 coefficient-level closure and its fail-closed evidence.

## Executable record and reproduction

The frozen [configuration](../configs/fgc/fgc-1-hyp1-mhg-propagation.toml)
binds ACT1, VAR1, RED1, SYM1, MHG1, REF1, and IMP1; the exact
[implementation](../src/recursive_horizons/fgc/modified_harmonic_propagation.py)
contains the independent two-jet operator and certificate; and the strict
[reproducer](../scripts/reproduce_fgc_hyp1_mhg_propagation.py) writes the
canonical [machine record](../results/fgc-1-hyp1-mhg-propagation.json). The
record hashes every predecessor configuration, the owner document, and the
transitive implementation/reproducer surface. Its nested certificate digest,
source ledgers, exact-check ledger, gate-status ledger, and false nonclaims are
independently rechecked by `scripts/check_repo.py`.

Reproduce and test only this gate with:

```bash
python3 scripts/reproduce_fgc_hyp1_mhg_propagation.py
python3 -m unittest -v \
  tests.test_fgc_modified_harmonic_propagation \
  tests.test_fgc_hyp1_mhg_propagation_reproduction
```

Or regenerate and verify the complete FGC chain through PROP1 with:

```bash
make verify-fgc
```

The strict loader rejects unknown keys, noncanonical rationals, source-path
traversal, convention/fixture/factor drift, a changed Noether sign or output
order, and any promoted open gate. Exact zero-`C` activated data also verify
that the reference connection is not inserted as an independent forcing term.

## What this changes

MHG1 had only the gauge-error principal block. REF1 supplied the complete
gauge-fixed residual, and IMP1 supplied one local acceleration branch. PROP1
adds the complete exact lower-order operator that the tensorial gauge vector
must satisfy **conditionally** when the gauge-fixed metric equations and both
unmodified scalar equations hold.

That is not yet a propagation theorem for an initial-boundary-value problem.
To infer that initially zero `C` remains zero requires a retained hyperbolic
domain, uniqueness/energy control, compatible physical constraints, and
constraint-preserving boundary or regular-centre data.

The downstream [FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) artifact now supplies
the exact local eighteen-variable first-order implicit lift, all complete
residual argument derivatives, the flat-root linearized branch map, and the
kinematic radial reduction identity for its literal coordinate rows. It does
not turn PROP1's independent `C^mu` operator into a metric-derived gauge
propagation theorem or provide a nonlinear source evaluator.
The later [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) gate certifies a
unique REF1 acceleration root on a full-dimensional local rational box, but
does not change PROP1's independent-operator or constraint-propagation limits.
The later [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) artifact uses the
metric-defined REF1 gauge vector to close one compatible local root; it still
does not supply the third-jet identity or propagation theorem that PROP1
deliberately leaves open.

## Fail-closed nonclaims

PROP1 does not establish:

- a direct third-jet evaluation of the divergence of the complete
  metric-generated REF1 residual;
- that the independent gauge-vector probe is a solution or a metric-derived
  gauge constraint through second order;
- an unconditional homogeneous equation off the scalar equations;
- a nonlinear first-order radial source system;
- reduction-constraint or physical-initial-constraint propagation;
- compatible initial data or boundary/regular-centre conditions;
- a PROP1-native quantified implicit-branch neighborhood (the downstream
  REF1 branch certificate is QIFT1, not a propagation theorem);
- a uniform eigenbasis, symmetrizer, open hyperbolicity domain, or retained EFT
  box;
- permission to evolve;
- collapse, affine metric-null defocusing, singularity resolution, a child
  spacetime, a dark sector, varying local light speed, a particle spectrum, or
  a theory of everything.

The exact result is narrower: PROP1 publishes and cross-checks the complete
lower-order ACT1/REF1 reference-gauge operator on declared exact local
backgrounds, with its Noether/scalar-equation handoff and every downstream gate
kept explicit.
