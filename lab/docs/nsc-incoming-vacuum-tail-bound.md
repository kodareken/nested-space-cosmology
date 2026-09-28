# Incoming vacuum tail: a directed Riccati defect enclosure

The calculation compares the exact affine-horizon vacuum projector with the
existing order-16 Riccati projector on the fixed incoming collar, for all 32
retained non-LLL groups and their existing energy splits, 160 or 320. It changes
no state, geometry, physical scale, finite-offset archive, or metric history.

## Decision and stopping condition

Reuse `nsc_pg_high_energy.riccati_coefficients`, the authenticated retained
channel inventory, and the [paired approximate source tail](nsc-incoming-source-tail.md).
The missing quantity is a physical vacuum-projector error above the existing
energy endpoints. Its size decides whether the approximate tail is adequate
in the stated affine-horizon domain. Enclose the recurrence defect over eight
whole radial cells, integrate each positive energy-power bound analytically,
and sum all retained groups. Stop after one all-group pass and focused
verification; do not run scattering, horizon, field, or old mode generators.
The per-group component budget is `3e-12`, and the all-group component budget
is `32 * 3e-12 = 9.6e-11`. A broad enclosure remains OPEN; it is not rounded
down or replaced by measured order differences.

Three material obstacles remain explicit: interval dependency can make the
bound broad; finite-offset initialization is different from a horizon limit;
and this tail estimate cannot resolve finite-band, subgap, or angular/compact
truncation errors. None is addressed by assigning a new stress or changing a
physical parameter.

## Reused projector estimate and finite polynomial defect

Write `A=-a^2`, `U=lambda/r+i*m`, and

$$
S_{16}=\sum_{n=1}^{16}\frac{c_n}{(2E)^n},\qquad
v=(1,-iaS_{16})^T,\qquad P_{16}=\frac{vv^\dagger}{v^\dagger v}.
$$

The owned recurrence cancels powers 0 through 15 in

$$
F=2ES_{16}-U-i(AS_{16}'+A'S_{16}/2)-A\bar U S_{16}^2
 =\sum_{n=16}^{32}\frac{f_n}{(2E)^n}.
$$

For the same canonical Dirac generator, the projector transport defect has
operator norm `|F|/[a(1+a^2|S16|^2)]`. The standard unitary variation-of-constants
estimate bounds the projector difference by its initial difference plus the
integrated transport defect. The rank-two trace-norm estimate then bounds
each vertex contraction by twice its operator norm times that difference.
The scalar/projector equality and recurrence ownership have independent
algebraic tests. No physical transport equation is solved here.

On `1 <= rho <= rho_h`, `A''<0` and `A'(rho)>=A'(rho_h)>0`. With
`rho=1+(rho_h-1)(1-t^2)`,

$$
\frac{|d\rho|}{a}\le
 2\sqrt{\frac{\rho_h-1}{A'(\rho_h)}}\,dt.
$$

Directed interval arithmetic encloses every `f16,...,f32` over each complete
cell, including a numerical enclosure of the existing profile's horizon.
Positive interval upper sums bound the radial integral. Energy powers are
integrated exactly to infinity; no quadrature samples are treated as maxima.
The existing signed-group factors are applied once. The vacuum current error
is exactly zero because both projectors have trace one, and its vertex is a
multiple of the identity.

## Initial condition and precise scope

Both the declared affine-horizon vacuum and `P16` approach the same rank-one
canonical projector at the horizon. Thus their **limiting-domain** initial
projector difference is zero. This is not a statement that the finite-offset
initialization in archived exact mode rows has zero error. For any separately
specified finite-offset vacuum, the same estimate contains an additional
nonnegative initial-projector term. Its energy-weighted integral has not been
bounded here and is recorded as `null`, never zero. The minimal remaining
connection is a bound on that initial-projector difference, with integrable
high-energy decay, using the actual archived initialization prescription.

The certificate concerns vacuum error above the original 160/320 endpoints.
It does not certify the order-16 approximation on the middle bands, subgap
quadrature, the angular/compact complement, changed normal jets, full source
convergence, incoming constraint stationarity, or extended EXISTENCE.
Thermal corrections retain their independent bound in the paired-tail owner.

## Preparation and replay

`--prepare` performs this bounded interval calculation once and writes an
authenticated, content-addressed coefficient-bound artifact and result.
`--check` authenticates its inputs and replays the aggregate and scope from
that artifact without repeating interval coefficient preparation.

```sh
python3 scripts/derive_nsc_incoming_vacuum_tail_bound.py --prepare
python3 scripts/derive_nsc_incoming_vacuum_tail_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_vacuum_tail_bound.py
```
