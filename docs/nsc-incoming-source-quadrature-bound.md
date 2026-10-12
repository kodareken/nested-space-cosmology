# Directed quadrature for the incoming order24 source

For each eight-unit cell, the existing source kernel is analytically continued
to a radius-eight disk about its center. With halfwidth `h=4`, a bound `M` on
the circle and a positive n-point Gauss rule, the imported Cauchy estimate and
polynomial exactness give

$$
\left|\int_{c-h}^{c+h}f(E)\,dE-Q_n(f)\right|
\le\frac{4hM(h/R)^{2n}}{1-h/R},\qquad R=8,\quad n=48.
$$

The two factors `2h` bound the integral and the positive quadrature functional
on the Taylor remainder after degree `2n-1`. The source-specific work is the
analytic domain, circle enclosure and directed evaluation of the exact rule.
The standard Gaussian rule is imported from [NIST DLMF 3.5(v)](https://dlmf.nist.gov/3.5#v);
analytic Taylor estimates use [DLMF 1.9(iv)](https://dlmf.nist.gov/1.9#iv).

## Bounded pilot and reuse decision

Reuse the existing `window_kernel`, its generated ad4 expressions, and the
pointwise order24 Riccati coefficient recurrence. The pilot is group6 LOW
`[32,40]` and its existing middle interval `[40,160]`. The missing connection
is rigorous quadrature and arithmetic for the direct LOW vacuum and the
middle order24-minus16 correction. Run one interval evaluation of these two
windows, authenticate it, and stop once its lapse-projected numerical error
is below `3e-11` or a concrete interval obstruction is recorded. No radial
geometry cells, mode solve, scattering solve or physical-defect calculation
is repeated. The original source arrays remain intact.

The three material risks are conjugating complex energy, an unexcluded
normalization pole, and mistaking numerical precision agreement for a
directed enclosure. The implementation handles each explicitly.

## Analytic and arithmetic domains

The continuation uses `S(E)*Sbar_coeff(E)`, where the second polynomial
conjugates coefficients at the **same E**. This agrees with `abs(S)**2` on
the real axis without introducing an antiholomorphic dependence. Both
order16 and order24 normalization denominators are pole-free because an
interval triangle majorant proves `abs(a²*S*Sbar_coeff)<1` on the disk.
The disk rectangle satisfies `Re(E)>abs(Im(E))`, so the energy-square gap
has positive real part. Its square roots continue the positive real branch.
The other source coefficients depend on fixed positive real incoming
geometry. Thirty-two directed angular arcs cover each boundary circle.

The generated ad4 code is reused with an interval numerical backend. Its
integer and half-integer literals are exact; a guard rejects other floating
constants. Unsupported complex interval square-root dispatch is replaced by
the interval half-power operation on the proved branch. This changes no
adiabatic expression or source convention.

Approximate Gauss nodes locate brackets only. Directed Legendre endpoint
signs prove one exact root in each of n disjoint brackets inside `(-1,1)`;
degree n accounts for every root. The derivative formula encloses the exact
positive weights. Source evaluation at those root intervals and subsequent
weighted sums enclose arithmetic error. Exact defining channel labels and
their stored binary versions are enclosed together. Binary interval
endpoints are serialized losslessly in the content-addressed artifact.

The record provides a new certified approximant, its integral enclosure and
separate analytic quadrature and arithmetic/parameter bounds. It does not
certify the old floating quadrature merely because it has similar values.
Replay validates roots, weights, cell coverage and saved interval sums
without reevaluating the source or its point coefficients.

## Scope

The [pilot record](../results/development/nsc-incoming-source-quadrature-bound.json)
binds the existing signed panels and all producer/input hashes. This closes
only these specified numerical integrals. Physical projector/thermal errors
still belong to their separate radial and state certificates. In particular,
the middle certificate covers **the order24-minus16 correction only**;
quadrature of the unchanged order16 baseline remains separate. No complete
source, physical IV, constraint solution or metric evolution is claimed.

```sh
python3 scripts/derive_nsc_incoming_source_quadrature_bound.py --prepare
python3 scripts/derive_nsc_incoming_source_quadrature_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_source_quadrature_bound.py
```
