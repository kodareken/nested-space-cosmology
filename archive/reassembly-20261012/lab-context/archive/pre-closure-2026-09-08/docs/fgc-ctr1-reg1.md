# FGC-1-CTR1-REG1: exact regular-centre REF1 formulation

## Result

**FGC-1-CTR1-REG1** removes the annulus-only centre limitation from the
declared spherical REF1 formulation at the level required for later initial
data and numerical boundary construction.

The gate does not evaluate the old spherical reference connection at
`r=0`. That coordinate object is undefined there. Instead it:

1. changes to centre-smooth variables `v=r V` and `R=r A`;
2. enforces the parity and elementary-flatness conditions required by a
   smooth Cartesian origin;
3. rederives every affected flat-reference connection difference and REF1
   descendant in exact Laurent arithmetic;
4. defines the finite centre equation vector

   ```text
   (E_tt, E_tr/r, E_rr, E_theta_theta/r^2, E_phi, E_chi);
   ```

5. serializes its exact centre Taylor coefficients;
6. verifies exact Minkowski and smooth nontrivial FGC-QR controls;
7. verifies finite Ricci, Ricci-squared, Kretschmann, and Gauss--Bonnet
   series; and
8. confirms that the unchanged positive-radius REF1 evaluator approaches the
   derived centre limits at least quadratically at three exact first-grid
   radii.

The result is local and conditional. It supplies no constraint-compatible
finite-mass initial-data family, solution-existence theorem, nonlinear
boundary operator, collapse trajectory, retained-EFT authorization, or
physical transition.

## 1. Regular variables and parity

Use the spherical ADM form

\[
ds^2=-\alpha^2dt^2+\lambda^2(dr+v\,dt)^2+R^2d\Omega^2.
\]

At a smooth rotational centre, a radial scalar is even and a radial vector
component is odd. CTR1 therefore writes

\[
v(t,r)=rV(t,r),
\qquad
R(t,r)=rA(t,r),
\]

and requires

\[
\alpha,\lambda,A,V,\phi,\chi
\quad\hbox{to be even in }r.
\]

Their first radial derivatives are consequently odd. In the original base
metric variables,

\[
h_{tt}=-\alpha^2+\lambda^2r^2V^2,
\qquad
h_{tr}=\lambda^2rV,
\qquad
h_{rr}=\lambda^2.
\]

Thus `h_tt` and `h_rr` are even, `h_tr` is odd, and `R` is odd. These are the
component parities induced by a smooth spherically symmetric tensor in a
Cartesian chart; they are not optional numerical reflection rules.

## 2. Elementary flatness

The spatial line element near the centre is

\[
dl^2=\lambda^2dr^2+r^2A^2d\Omega^2.
\]

The proper radial distance and areal circumference agree to leading order
only when

\[
A(t,0)=\lambda(t,0)>0.
\]

Because the equality is a condition at every time, its first and second time
derivatives must also agree:

\[
A_t(t,0)=\lambda_t(t,0),
\qquad
A_{tt}(t,0)=\lambda_{tt}(t,0).
\]

CTR1 enforces all three coefficient equalities before evaluating a centre
series. It does not numerically suppress a resulting pole.

The necessity is visible directly in the warped-product Ricci scalar. If
`A_0` and `lambda_0` differ, its leading singular coefficient is

\[
{}^{(4)}\!R
=\frac{2\left(1-A_0^2/\lambda_0^2\right)}{A_0^2r^2}
+O(1).
\]

For the frozen negative control `A_0=11/10`, `lambda_0=1`, the exact
coefficient is `-42/121`. The profile is rejected as a conical defect before
a positive certificate can be constructed.

Parity is independently necessary. If a nominal scalar has

\[
\psi=\psi_0+\psi_1r+O(r^2),
\]

then the angular part of its wave operator contributes

\[
\Box\psi=\frac{2\psi_1}{\lambda_0^2r}+O(1).
\]

The frozen odd-`phi` mutation has `phi_1=1/17` and therefore the exact pole
coefficient `2/17`. It too is rejected before certification.

## 3. The flat spherical reference is regularized by differences

The original REF1 reference geometry is flat spacetime in spherical
coordinates. Its individual connection coefficients include `1/r` and their
derivatives include `-1/r^2`; the coordinate chart itself has no value at the
origin. CTR1 therefore keeps the same fixed geometric reference on the
punctured neighbourhood but rederives the tensorial connection differences
before taking the centre limit.

For the warped angular block,

\[
\Gamma^A{}_{ij}
=-R\,h^{AB}\partial_BR\,\gamma_{ij},
\qquad
\Gamma^i{}_{Aj}
=\frac{\partial_AR}{R}\,\delta^i{}_j.
\]

Subtracting the flat-spherical reference gives

\[
\Delta\Gamma^A{}_{ij}
=\left(-R\,h^{AB}\partial_BR+r\,\delta^A_r\right)\gamma_{ij},
\]

and the only apparently singular mixed quotient becomes

\[
\frac{R_r}{R}-\frac1r
=\frac{A+rA_r}{rA}-\frac1r
=\frac{A_r}{A}.
\]

Its radial derivative is likewise finite:

\[
\partial_r\left(\frac{R_r}{R}-\frac1r\right)
=\partial_r\left(\frac{A_r}{A}\right)
=\frac{A_{rr}}A-\left(\frac{A_r}{A}\right)^2.
\]

The temporal quotient is

\[
\frac{R_t}{R}=\frac{A_t}{A}.
\]

The radial angular contraction additionally contains the difference between
`1/A^2` and `1/lambda^2`, divided by `r`. Both functions are even and share
their centre value by elementary flatness, so their difference begins at
`r^2` and the contracted expression is odd and finite.

These identities are evaluated entry by entry in the exact connection and
connection-derivative tensors. CTR1 requires every certified negative Laurent
coefficient of `Gamma-bar_Gamma` and its derivative to vanish before the
metric-defined gauge vector or gauge extension is assembled.

The completeness of the singularity analysis follows from the standard
warped-product basis. Apart from the finite two-dimensional base curvature,
every four-dimensional curvature component is generated by

\[
\frac{\nabla_A\nabla_BR}{R}
\quad\hbox{and}\quad
\frac{1-(\nabla R)^2}{R^2}.
\]

For even regular primitives, the first quotient is finite. Elementary
flatness cancels the constant numerator of the second, leaving a finite even
limit. The only additional angular scalar-Hessian quotient is

\[
\frac{\nabla R\mathbin{\cdot}\nabla\psi}{R},
\]

which is finite for an even scalar. Einstein, nonminimal, scalar-stress, and
double-dual Gauss--Bonnet terms are algebraic contractions of this finite
basis. The auxiliary metrics, trace-reversal projectors, and `nabla C` that
form REF1 add no independent spherical denominator once the reference
difference identities above are used. This is the analytic closure behind
the component-level Laurent certificate.

This satisfies the frozen protocol's reference policy by rederiving all
affected REF1 descendants. It does not pretend that the annular reference
object itself can be called at `r=0`.

## 4. The regular centre equation vector

Smooth spherical tensor parity implies

\[
E_{tt}=E_{tt}^{(0)}+O(r^2),
\qquad
E_{tr}=rE_{tr}^{(1)}+O(r^3),
\]

\[
E_{rr}=E_{rr}^{(0)}+O(r^2),
\qquad
E_{\theta\theta}=r^2E_\Omega^{(2)}+O(r^4),
\]

while both scalar equations are even. CTR1 therefore defines

\[
\mathcal E_{m reg}=
\left(
E_{tt},
\frac{E_{tr}}r,
E_{rr},
\frac{E_{\theta\theta}}{r^2},
E_\phi,
E_\chi
\right).
\]

The active gauge components have the corresponding form

\[
C^t=C^t_0+O(r^2),
\qquad
C^r=r\,\mathcal C^r_0+O(r^3),
\]

so the regular gauge pair is `(C^t,C^r/r)`.

Every element of both vectors is even. Their centre equations are the exact
zeroth Taylor coefficients

\[
\mathcal E_{{\rm reg},I}^{(0)}=0.
\]

The machine record serializes, for each of the six equations, its centre
coefficient, its next even coefficient, and the complete certified negative-
power ledger. At a regular centre, radial and angular isotropy also requires

\[
E_{rr}^{(0)}=E_\Omega^{(2)}.
\]

The nontrivial exact control satisfies this identity with the common value

```text
-14836176781/62092800000.
```

This centre Taylor system is a regular residual formulation. It is not yet a
solved centre acceleration system or a Cauchy solution; those require the
later initial-data, run-domain, and numerical gates.

## 5. Exact Laurent construction and certified window

The implementation uses standard-library rational arithmetic only. It
constructs the complete ACT1/VAR1 tensor residual through the sealed RED1
evaluator and separately rebuilds the reference differences, auxiliary
metrics, metric-defined gauge vector, its covariant derivative, and the full
REF1 extension in Laurent series.

The internal series window is

```text
-12 <= power of r <= 20.
```

The public coefficient tables display only

```text
-4 <= power of r <= 4.
```

The wider internal window covers every intermediate angular inverse-metric
factor through the curvature-squared and Gauss--Bonnet contractions. The
controls use finite quartic primitive profiles, so omitted higher primitive
coefficients are exactly zero. No decimal approximation, symbolic-algebra
package, or substitution `r=epsilon` participates in the centre limit.

Every computed negative coefficient from `r^-12` through `r^-1` is checked
for zero. The structural warped-product budget shows that a final physical
pole can be no worse than `r^-4` (from a squared curvature); the larger
internal negative window is therefore also an underflow guard, not an
uninspected region. Any arithmetic product or derivative that would cross the
`r^-12` guard edge raises a typed error instead of silently truncating a more
singular coefficient.

The regular gate requires the certified negative coefficients to vanish for:

- every reference-connection difference;
- every first derivative of that difference;
- all six components of `E_reg`;
- both active regular gauge components;
- the Ricci scalar;
- Ricci contraction `R_ab R^ab`;
- Kretschmann scalar `R_abcd R^abcd`; and
- the Gauss--Bonnet invariant.

## 6. Exact controls

### Minkowski

The first profile fixes

\[
\alpha=\lambda=A=1,
\qquad
V=\phi=\chi=0
\]

with all time derivatives zero, while retaining the frozen FGC-QR couplings.
The exact formal centre residual, regular gauge pair, Ricci scalar,
Ricci-squared, Kretschmann scalar, and Gauss--Bonnet invariant all vanish. The
unchanged annular REF1 evaluator also returns exact zero at `r=1/16`.

### Smooth nontrivial series

The second profile gives every primitive field independent rational
coefficients at powers `r^0`, `r^2`, and `r^4`, together with independent
first and second time-derivative series. Its centre is shifted, time dependent,
scalar active, and curved. It is deliberately not a solution.

All regularity checks pass exactly. In particular,

```text
Kretschmann(0) = 18910513/17280000
Ricci(0)       = 2701/2400
```

so the finite-curvature test cannot pass through an accidentally flat control.
The fact that its residual is nonzero is intentional: centre regularity is
being tested independently of solving the equations.

## 7. First-grid approach to the analytic limit

For an even regularized component `f(r)`,

\[
f(r)-f(0)=O(r^2).
\]

CTR1 evaluates the unchanged positive-radius REF1 code at the exact radii

```text
1/16, 1/32, 1/64
```

using jets generated from the same finite polynomial profile. For every
component of `E_reg` and for both regular gauge components, it computes the
exact rational errors against the independent formal centre coefficient.
Whenever the fine error is nonzero, halving the radius must reduce the error
by at least the frozen ratio

\[
\frac{15}{4}.
\]

All eight components pass both halvings. This is a first-grid-point
consistency control for the derived centre limit. It is explicitly not a PDE
discretization-convergence result and does not stand in for NUM1.

## 8. Machine reproduction

```bash
python3 scripts/reproduce_fgc_ctr1_reg1.py \
  --output results/fgc-1-ctr1-reg1.json
```

The canonical result is hash-bound to ACT1, VAR1, the original annular REF1
formulation, SRC1, CON4, the frozen SF1 protocol, RUN1, the new exact series
implementation, and this derivation document.

The gate

```text
regular_center_formulation_verified = true
```

closes only the regular-centre evidence slot. RUN1 combines that slot with
the still-absent finite-mass compatible initial-data family, so the composite
RUN1 predicate remains false and the FGC-QR holdout remains forbidden.

## Nonclaims

CTR1-REG1 does not derive or authorize:

- a constraint-compatible initial-data family;
- finite ADM or Misner--Sharp mass;
- an outer vacuum buffer;
- a nonlinear centre acceleration solve or local existence theorem;
- a constraint-preserving boundary condition or complete IBVP;
- a classical collapse run;
- retained-Wilsonian-EFT validity;
- affine metric-null defocusing or singularity resolution;
- a child domain or topology change;
- a dark-sector mechanism; or
- a varying locally measured speed of light.
