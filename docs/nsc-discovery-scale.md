# Initial source-selected scale

This note derives the areal radius selected by the initial lapse constraint
of the same conformal action. It does not evolve the state, does not call
the initial-radius Newton, and does not write a production record. A later
source-consistent solve can compare its radius, clock and proper ratio with
the prediction below. That comparison is not a dynamical selection.

The owned sector is constant \(Q=Q_0=b_0/a_0\), \(\chi=0\) and
\(p_Q=p_r=p_\chi=0\). Outside that sector the extra momentum terms in
`hamilton_constraint` remain, and this radius formula does not apply.

## Verified operator

The owned geometric density at \(\chi=0\) and zero momenta is

\[
C_g=\frac{Z r_x^2}{Q}-QV-2\partial_x(F_x/Q),
\]

with \(F=-4\pi A r^2\), \(Z=-24\pi A\) and
\(V=8\pi A r^2-2\pi C_F\mathrm{flux}^2\). The matter density is the owned
nodal force divided by the sample spacing,

\[
\rho=\frac{M}{\Delta x}\left(\frac{K}{Q_0}+\kappa\operatorname{Re}S_1\right),
\qquad M=4\kappa,
\]

applied once. \(K\), \(S_1\) and \(P_{\mathrm{mom}}\) are the bilinears from
`column_moments`. The shift source is the owned sign

\[
j=-\frac{M P_{\mathrm{mom}}}{\Delta x}.
\]

Its mean is retained. Setting the conformal shift to zero does not solve
\(D_g+j=0\).

Set \(r=y^2\) with \(y>0\). Where the spectral product rule holds,

\[
C_g+\rho=-\frac{8\pi A}{Q_0}y^3 R_y,
\qquad
R_y=-4 y''+Q_0^2 y-\frac{G}{y^3},
\]

\[
G=Q_0^2 r_{\mathrm{mag}}^2+\frac{Q_0\rho}{8\pi A},
\qquad
r_{\mathrm{mag}}^2=\frac{C_F\mathrm{flux}^2}{4A}.
\]

So \(R_y=0\) is

\[
(-4\partial_x^2+Q_0^2)y=G\,y^{-3}.
\]

A unit coefficient on \(\partial_x^2\) does not reproduce
`hamilton_constraint`. On the discrete grid the product rule need not hold.
The form that matches the owner residual without that rule is

\[
C=\frac{8\pi A}{Q_0}\bigl(D^2(r^2)-3(Dr)^2-Q_0^2 r^2+G\bigr).
\]

The difference between this expression and the \(y\)-form is
\((8\pi A/Q_0)\) times \(D^2(r^2)-3(Dr)^2-4 y^3 D^2 y\). The projected
radius equation is the geometry-band pullback of this same residual. It is
not a new constraint.

At \(\chi=0\) the Weyl auxiliary \(C_W\) drops out of \(G\). It still
defines the chart. The prescribed lapse profile and the shift \(\beta\) do
not enter \(\rho\). They do enter \(F_Q\). Ambient dependence that does
enter \(G\) is \(Q_0\), \(A\), \(C_F\), \(\mathrm{flux}\), \(\kappa\),
\(M=4\kappa\), and the column bilinears \(K\) and \(\operatorname{Re}S_1\).
No Hubble parameter and no magnetic-field value is adjusted to match a
target radius. \(r_{\mathrm{mag}}^2\) is the locked ratio above.

## Integrated identity and conditional bounds

Multiply the owned equation by \(y\) and integrate over the period. Skewness
of the periodic derivative gives

\[
\int\bigl(4(y')^2+Q_0^2 y^2\bigr)\,dx=\int\frac{G}{y^2}\,dx.
\]

If \(G>0\) and \(y>0\) solves the equation, the maximum principle yields

\[
\frac{G_{\min}}{Q_0^2}\le y(x)^4\le\frac{G_{\max}}{Q_0^2}.
\]

The lower bound is withheld when \(G\) is not positive. These bounds are
properties of a positive solution. They are not a certificate that a
particular numerical solve has converged.

## Uniform source

If \(G\) is constant the positive solution is constant,

\[
r_{\mathrm{flat}}^2=y^4=\frac{G}{Q_0^2}
=r_{\mathrm{mag}}^2+\frac{\int\rho\,dx}{8\pi A Q_0 L}.
\]

On a uniform collocation grid the same constant radius makes every
derivative term vanish, so the owned residual is exactly

\[
C=\rho-\operatorname{mean}(\rho).
\]

It is zero if and only if \(\rho\) is uniform. For a nonuniform source the
same \(r_{\mathrm{flat}}\), built from the mean alone, still solves the
constant mode and leaves the contrast as the residual. Shape therefore
changes the source contrast and the conditional gap below. It does not
change \(r_{\mathrm{flat}}\) when the integral of \(\rho\) is held fixed.

If a positive solution of the initial lapse constraint exists and \(G>0\),

\[
\bigl|r(x)^2-r_{\mathrm{flat}}^2\bigr|
\le\frac{\|G-\operatorname{mean}(G)\|_\infty}{Q_0^2}.
\]

The relative size of that gap is \(\|G-\operatorname{mean}(G)\|_\infty/\operatorname{mean}(G)\).
The gap is the applicability statement. It is not a compulsory one-percent
test, and a small gap on one preparation does not become a pass rule for
another.

The conformal normal clock is \(d\tau=r Q_0\,dt\). The prediction reports

- radius \(r_{\mathrm{flat}}\),
- clock rate \(r_{\mathrm{flat}} Q_0\),
- proper ratio \(d\tau/(Q_0 dt)=r_{\mathrm{flat}}\).

A later initial solver can test those three numbers. This module does not.

## Coordinate weight

The owned affine pullback \(x=b+a\xi\) gives weight 0 to \(r\) and weight 1
to \(Q\). In that chart \(\rho\) scales by \(a\), the coordinate period
scales by \(1/a\), and \(\int\rho\,dx\) is invariant, so \(Q L\) and
\(r_{\mathrm{flat}}\) are invariant. The areal-radius weight of the
covariant control is zero. Stretching only the period, or assigning \(r\)
the weight of \(Q\), changes the constraint and is not a pullback. A
coordinate change cannot be used to select the physical radius.

The clock rate \(r Q\) has weight 1 because \(Q\) does. The proper ratio
\(r\) keeps weight 0.

## Held-out preparation

The withheld family is an override passed through the existing
`columns_override` gate. It is not the saved `original` layout and not the
saved `separated` layout. Those names stay previously saved preparations
even if a caller attaches the new carrier frequency to them.

Each member uses the owned bump envelopes, then an explicit real Löwdin of
those envelopes only. The spinor frame is not repaired afterward, and the
original child columns are not subtracted. The carrier is

\[
e^{i\nu(x-n\ell)},\qquad \nu=1.15,
\]

which is not the saved carrier wavenumber \(k=1\). The plus column is the
chirality \((1,i)/\sqrt{2}\). The minus column is \(e^{i\sigma\alpha}\) times
its conjugate, with \(\sigma=\pm 1\) and \(\alpha=0.37\), not the saved
minus-column phase. Occupations stay the owned CAR weights in \((0,1]\).
A column phase preserves the Gaussian and the Gram, so the two signs are
predicted to share source invariants. For this chirality \(\operatorname{Re}S_1\)
is zero before collocation roundoff. The frequency is not a column phase.
A control at the saved wavenumber, built by the same column code and
measured on the same quadrature grid, is an override frequency control. It
is not labeled held out.

The prediction reads \(\int\rho\), the contrast of \(G\), and the retained
current from `source_from_columns`. It does not copy the probe radius
\(r=1\). Source positivity, the lapse sector, and the shift residual are
reported as three domains. The shift residual is the current, including
its mean.

## Reproduction

From the repository root, with the pinned validation interpreter:

```sh
.venv/validation/bin/python scripts/lab.py -m pytest \
  tests/test_nsc_discovery_scale.py -q
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_scale.py
```

`--solve` and `--output` are refused. No result file is written.
