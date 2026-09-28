# Remaining incoming high-band source updates

## Decision and stopping condition

Reuse the frozen source-window helper, authenticated LOW/middle arrays, the
same local Riccati recurrence and unchanged affine-horizon source law. This
batch completes the explicit numerical source updates for groups1–32 except
12,14,22,32, whose updates already have separate records. Group6 is the first
gate: calculate its LOW `[32,40]` and middle `[40,160]` sources and validate
them before dispatching the remaining27 groups.

The missing quantities are the actual signed source values and additive
corrections matching the independent order24 radial batch. No radial,
physical mode or scattering generator is called. Each group receives one
48/64-node quadrature resolution and50/70-digit control; a failed control is
reported, not followed by an unbounded refinement loop. Use at most two CPU
workers. Content-addressed per-group artifacts and authenticated cache
pointers allow interruptions to resume only missing groups.

Three expected obstacles are accidental overlapping intervals, repeated
coefficient/quadrature work, and conflation of small numerical indicators
with physical error bounds. Authenticate actual interval metadata; cache
local geometry, quadrature rules and coefficients; keep physical mode-error
certification OPEN until the matching independent radial certificates exist.

## Exact interval coverage

All ordinary retained families use a selected LOW `[32,40]` cell followed
by their archived middle interval, whose upper endpoint is160 or320.
**Group13 already starts its middle interval at16**, so its only update is
middle `[16,160]`. No LOW `[32,40]` term is added for group13. Group14's
existing early middle interval is excluded with its completed separate work.
The old per-sign arrays and group factors remain unchanged.

LOW updates are direct order24 vacuum approximations minus the selected old
LOW contribution. The physical thermal source remains a nonzero bounded
remainder. Middle updates are explicit order24-minus16 vacuum corrections;
the old thermal insertion is retained, with the existing conservative
thermal/scattering difference bound. No physical covariance, scale, ad4
reference, seed, or action coefficient is selected or fitted.

## Faster evaluation with the same mathematics

The source owner evaluates the same finite polynomials with Horner's rule.
An exact symbolic identity checks the order16 polynomial, the powers17–24
correction and the full order24 polynomial. Prepared physical cases also
compare Horner evaluation to the frozen source-window kernel. The same
stable Bloch difference and raw metric vertices are used. This is a change
of arithmetic organization, not a new approximation order or state.

Positive quadrature nodes and fixed geometry are cached by arithmetic
precision; channel coefficients are shared between its two windows. Both
actual angular signs are retained, except the existing zero-angular families
which have one sign. Negative-frequency folding and degeneracy remain in
`incoming_group_factor`, included once.

## Evidence and replay

Each group artifact contains all original/source/correction arrays, signed
integrals, precision and quadrature indicators, producer/input hashes and
scope. Cache hits authenticate and replay these artifacts; they never
regenerate their source preparation. The final batch record enumerates
exactly28 groups and55 disjoint windows, excluding the four completed groups.

Numerical source checks do not certify physical mode error or rigorous
quadrature remainder. The physical source-bound field remains OPEN here
until the separate radial certificates are composed. Full source convergence,
joint stationarity and metric evolution are not claimed.

```sh
python3 scripts/derive_nsc_incoming_high_band_batch_source.py --prepare --workers 2
python3 scripts/derive_nsc_incoming_high_band_batch_source.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_high_band_batch_source.py
```
