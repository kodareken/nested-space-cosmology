# FGC-1-ID1-FAM1: constraint-compatible finite-mass data family

**Status:** initial-hypersurface certificate; no time evolution and no
mechanism outcome

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Human scientific responsibility:** Douglas Ek

## Question answered

`FGC-1-ID1-FAM1` asks one premise question:

> Does the unredefined FGC-QR classical action admit more than an isolated
> point of regular-centre, finite-Misner--Sharp-mass spherical initial data
> with an independent ingoing `chi` pulse, an explicit nonzero `phi` seed,
> and solved physical and metric-defined gauge constraints?

The certificate answers **yes within its declared initial-slice scope**. It
does not evolve the data, inspect an FGC-QR holdout outcome, establish a
healthy spacetime run domain, or test null defocusing.

This is the successor to the immutable `FGC-1-ID0-PREF1` diagnosis. ID0 found
that PROTO1 omitted the regulator momentum and that its declared matter
amplitudes could not reach its own initial-compactness floor. PROTO2 repaired
only those input premises before any FGC-QR outcome was inspected. ID1 tests
whether that repair produces executable FGC-QR initial data; it does not
reinterpret the protocol diagnosis as mechanism evidence.

## Frozen slice and fields

The slice uses unit lapse, zero shift, polar--areal radius, and maximal
extrinsic curvature,

\[
  \alpha=1,\qquad v=0,\qquad R=r,
\]

\[
  K^r{}_r=-2k,\qquad
  K^\theta{}_\theta=K^\varphi{}_\varphi=k.
\]

The independent matter pulse and regulator seed are compact smooth bumps,

\[
  r\chi(r)=A_\chi B\!\left(\frac{r-12}{w_\chi}\right),
  \qquad
  \partial_t(r\chi)=\partial_r(r\chi),
\]

\[
  \phi(r)=A_\phi B\!\left(\frac{r-12}{2}\right),
  \qquad \Pi_\phi=0,
\]

where

\[
 B(x)=
 \begin{cases}
 \exp\!\left(1-\frac{1}{1-x^2}\right),& |x|<1,\\
 0,& |x|\geq 1.
 \end{cases}
\]

The action parameters are the frozen ACT1 fixture values

\[
 M_{\rm Pl}=2,\quad \beta=-\frac14,\quad \mu=3,
 \quad g_4=\frac12,\quad \eta=\frac12.
\]

These are dimensionless code-unit fixtures. They are not a physical parameter
fit.

## Why the constraint solve is direct

On this slice, the complete ACT1/VAR1 physical normal projections specialize
to

\[
  \mathcal H=H_0(r,\lambda,k,\phi,\chi)
       +H_\lambda(r,\lambda,k,\phi,\chi)\,\partial_r\lambda=0,
\]

\[
  \mathcal M=M_0(r,\lambda,k,\phi,\chi)
       +M_k(r,\lambda,k,\phi,\chi)\,\partial_r k=0.
\]

The two cross coefficients vanish exactly. Therefore, wherever
\(H_\lambda M_k\neq0\),

\[
  \partial_r\lambda=-\frac{H_0}{H_\lambda},
  \qquad
  \partial_r k=-\frac{M_0}{M_k}.
\]

This discovery changes the implementation architecture materially: a generic
block-tridiagonal Newton solve is unnecessary for the declared maximal
polar--areal family. The repository uses a transparent two-component radial
ODE and monitors the exact diagonal constraint Jacobian. It does not replace
the original equations with a second implementation: two nontrivial rational
fixtures compare the specialized kernel exactly with the complete
`physical_constraint_projections` ACT1/VAR1 evaluator.

At zero regulator seed, the specialized FGC-QR constraint pair reduces
exactly to the GR-0 constraint pair. The nonzero seed is reached along the
frozen sequence

\[
  A_\phi/A_{\phi,0}\in
  \left\{0,\frac14,\frac12,\frac34,1\right\}.
\]

This is a premise continuation, not tuning against a collapse result.

## Gauge compatibility

For the frozen `tilde_normal_factor=4`, the two metric-defined gauge
constraints on the support are solved exactly by

\[
  \partial_t h_{tt}=0,
\]

\[
  \partial_t h_{tr}
  =\frac{2\lambda^3-2\lambda+r\,\partial_r\lambda}
         {4\lambda r}.
\]

The rational controls evaluate the full REF1 gauge expression and obtain
\(C^t=C^r=0\) exactly. ID1 does not claim a new theorem for the normal
derivatives of the gauge constraints. Their conditional propagation remains
owned by `FGC-1-CON4-PHY1` and its stated premises.

## Centre and exterior

The compact profiles have an exact centre vacuum buffer. Consequently the
regular solution there is

\[
  \lambda=1,\qquad k=0,
\]

with the centre regularity and elementary-flatness contract inherited from
`FGC-1-CTR1-REG1`.

Outside the combined compact support the exact vacuum continuation is

\[
  k=\frac{J}{r^3},
  \qquad
  \lambda^{-2}=1-\frac{2M}{r}+\frac{J^2}{r^4}.
\]

Here \(M\) is the constant exterior Misner--Sharp mass and \(J=r^3k\) is the
constant vacuum momentum parameter inherited at the support edge. Positivity
of the exterior metric denominator is checked analytically over the whole
declared exterior, not only at its endpoints.

## Certified family and numerical checks

The qualification box is

\[
  2\leq A_\chi\leq\frac94,
  \qquad
  \frac{15}{8}\leq w_\chi\leq2,
  \qquad
  \frac12\leq\frac{A_\phi}{A_{\phi,0}}\leq2,
  \qquad
  A_{\phi,0}=\frac1{131072}.
\]

Its exact three-axis parameter volume is \(3/8388608>0\). A sampled tensor grid
of 27 members is evaluated with both RK4 and SSPRK3 at the report resolution. The
central member is independently checked on three nested grids using the
common-grid infinity norm of both \(\lambda\) and \(k\). Every reported member
must satisfy frozen strict bounds on

- initial compactness and absence of an initially trapped sphere;
- physical-constraint residuals;
- the constraint-Jacobian determinant and condition number;
- the effective Planck coefficient;
- positivity of the complete exterior vacuum denominator;
- cross-method exterior mass and peak-compactness differences;
- regular-centre, exact inner-vacuum, and exact outer-vacuum conditions.

The finite tensor grid is **not** represented as an interval proof over the
entire rectangular box. The positive-family conclusion is narrower and
standard: the ODE right-hand side is smooth at the central member, all its
denominators have strict margins, and continuous dependence therefore gives
a nonzero open neighborhood of constraint-compatible solutions. The grid
then supplies quantified nearby stress points. A future interval enclosure
could strengthen the box-wide statement but is not needed for the existential
RUN1 predicate.

The upper width endpoint is deliberately `2`, preserving the original PROTO2
minimum centre buffer of `10`. Consequently this certificate does **not** qualify
the later held-out width factor `9/8`. PROTO3 prospectively repairs the
inconsistent round-number floor to the exact positive buffer `39/16 L0`, but
does not retroactively enlarge ID1. The expanded slice still requires its own
static constraint-compatible initial-data result before `PRO3-HLD1` or
robustness execution; it is not silently absorbed into ID1.

## Machine contract

Reproduce the canonical result with

```bash
python3 scripts/reproduce_fgc_id1_fam1.py
```

The result binds the configuration, this derivation document, implementation,
ACT1, VAR1, REF1, CON4, CTR1, ID0, RUN1, and the PROTO3 semantic envelope by
SHA-256. RUN1 consumes the result only when the exact artifact identity,
classification, scope bindings, payload, and gate all validate.

## Exact conclusion and nonclaims

The supported conclusion is:

> The unredefined FGC-QR action admits a regular-centre, finite-mass,
> nonzero-width local family of spherical constraint-compatible initial
> hypersurfaces with independent compact `chi` data and an explicit nonzero
> `phi` seed under the declared slice, parameter, and numerical controls.

The certificate does not establish evolution well-posedness on that family,
multidirectional hyperbolicity, collapse, horizon formation, regulator
activation, affine-null defocusing, singularity resolution, a child domain,
retained-EFT validity, a dark-sector mechanism, or varying locally measured
light speed. There was no time evolution performed. In particular, a classical
singularity remains a boundary of
regular classical continuation rather than a derived microscopic mechanism;
ID1 never reaches or interprets such a boundary.
