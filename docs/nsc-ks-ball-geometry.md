# Ball-arithmetic jets of the owned KS chart and flat cutoffs

This owner encloses Taylor coefficients of the existing KS background,
time shift and axial plateau. It does not change the metric, source or
scales, and it does not certify a field residual or the local gate.

On the restricted slab `rho in [1, 33/32]`,

\[
r=\sqrt{1+\rho^2},\qquad
a=\sqrt{3\bigl((1+\rho^2)(\pi/2-\arctan\rho)-\rho\bigr)-1},
\]

with `s=T(rho)-T(1)=-int_1^rho d rho/a`. The plateau is the same
C-infinity window already used by `LocalAxialFunction`,

\[
\eta(u)=\operatorname{expit}\bigl(1/u-1/(1-u)\bigr),\qquad
u=(\mathrm{distance}-\mathrm{inner})/(\mathrm{outer}-\mathrm{inner}).
\]

`arb_series` coefficients are enclosures of `f^{(n)}/n!`. Numerical
samples are containment controls, not proofs of those bounds.

## python-flint 0.9.0 APIs

Interval jets use `arb_series([rho_ball, 1])` at `ctx.cap=order+1`.
Algebraic operations then apply automatic differentiation with interval
coefficients, so each coefficient encloses the corresponding Taylor
term on the whole ball. This is valid enclosure, not a claim that the
interval form is sharp.

The constant term of `s` is not the termwise integral of that series.
`(-inv_a).integral()` supplies only the `n>=1` coefficients. The value
is `acb.integral(inv_a_integrand, 1, mid(rho))` from python-flint
0.9.0. The integrand must forward the integrator's `analytic` flag to
`acb.sqrt(analytic=...)`; omitting it is the documented incorrect
pattern for branch cuts. The analytic callback also restricts its region to
the open right half-plane, where atan has no branch cut. If `rho` is an interval, add the proved
variation `|T'|<=5/4` times the radius. That Lipschitz bound is
`1/a<=5/4` from `a>=4/5`, proved with `pi>3.14` and
`atan(1/65)<1/65` at the right endpoint `rho=33/32=1+1/32`, where
`atan(33/32)=pi/4+atan(1/65)`.

Binary cutoffs and Chebyshev `mapparms()` enter as `arb(float(...))`,
the exact binary64 dyadic, never as a decimal string at working
precision.

## Cauchy flat-strip bound

`eta` is not analytic at the flat endpoints. A cell that meets `u=0`
or `u=1` is never passed through `1/u`. For `u in (0,eps]`,
`eps=1/256`, set `f(z)=exp(-1/z+1/(1-z))` and `1-eta=f/(1+f)`. On the
disk `|z-u|<=u/2`,

\[
\operatorname{Re}(1/z)\ge 2/(3u),\qquad
|1/(1-z)|\le 1/(1-3u/2),
\]

so `|f|<=t(u)=exp(1/(1-3u/2)-2/(3u))`. If `t<1`, then
`|1-eta|<=t/(1-t)` on the disk and Cauchy's estimate gives

\[
\frac{|(1-\eta)^{(n)}(u)|}{n!}
\le\Bigl(\frac{2}{u}\Bigr)^n\frac{t(u)}{1-t(u)}.
\]

For `n<=16`, `eps<=2/(3n)`, so `n\log(2/u)-2/(3u)` is increasing on
`(0,eps]` and the maximum is at `u=eps`:

\[
\frac{|(1-\eta)^{(n)}|}{n!}
\le\Bigl(\frac{2}{\varepsilon}\Bigr)^n\frac{t(\varepsilon)}{1-t(\varepsilon)}.
\]

The inner function extends by `eta=1` (all derivatives of the
complement vanish). The right endpoint uses `eta(u)=1-eta(1-u)`. Those
derivative balls are composed with the nonconstant part of the
distance series, with the constant coefficient of that series replaced
by exact `0` rather than cancelled in interval arithmetic. A cell that
cannot be proved wholly inner, wholly outer, interior to `(0,1)`, or
inside one `eps`-strip raises `SubdivisionNeeded`.

## Unresolved edges

- A rho ball outside the declared slab, or one for which interval evaluation
  does not prove a positive axial scale. Outward-rounded balls may require
  a carefully enclosed boundary cell rather than decimal-radius parsing.
- Plateau cells that meet a flat endpoint and leave its `eps`-strip.
- z balls that straddle the plateau centre and extend outside the inner
  plateau. Cells wholly inside the plateau use its exact constant value.
- Jet order greater than the implemented limit 16. The Cauchy inequality
  above is checked explicitly; 16 is an implementation limit.
- Using the jet as a polynomial away from the expansion point needs a
  separate remainder; this owner returns the jet only.

```sh
python -m pip install -r requirements-validation.txt
PYTHONPATH=src python -m pytest -q tests/test_nsc_ks_ball_geometry.py
```
