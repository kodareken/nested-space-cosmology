# Complete validated group32 source on E24–32

This owner connects the uniform real-energy field enclosure to its actual
incoming source integral. The raw degree24 endpoint columns give a degree48
covariance; a linear metric vertex gives degree49. Exact Chebyshev moments
integrate that polynomial, and a separately certified Gauss48 calculation
subtracts the unchanged full ad4 reference.

## Authorized fixed run and stopping condition

Reuse the [successful preflight](nsc-incoming-energy-panel-preflight.md),
its first completed step for both signs, the certified affine-vacuum
initializer and the existing unitary residual theorem. The missing quantity
is the complete group32 source on the original LOW cell `[24,32]`, with
uniform field, integration, arithmetic and thermal errors. Run degree24 in
energy, degree36 in the logarithmic radial coordinate,699 real steps of at
most1/32, analytic radius1/8 and80-digit directed arithmetic. Do not rerun
the pointE24 calculation or the successful preflight step.

The complete window lapse/shift budget is `3e-12`. At most two sign workers
share one locked geometry cache entry per step. Save checkpoints every16
steps. The preflight's176 CPU-minute planning allowance is divided between
the signs; a cost overrun or already excessive cumulative field bound stops
at a checkpoint with an explicit obstruction. No automatic degree increase,
state adjustment or silent restart is permitted.

The principal risks are losing a cached-step boundary, double counting an
old LOW/high source, and claiming quadrature accuracy from a field endpoint
alone. The record binds exact coverage, uses an explicit disjoint replacement,
and carries all numerical and physical remainders through the observable.

## Reuse and continuous field proof

The first real step has exactly the preflight's size1/32; subsequent steps
use that spacing with one final shorter step. Thus all699 intervals form a
single contiguous domain, and the first step is reused byte-for-byte in its
numerical fields. Source/input hashes and a live sign lock prevent a second
worker from duplicating an active run. Completed receipts are authenticated
before reuse. A restarted process resumes only completed checkpoint rows.

The scalar-phase gauge and energy-degree25 residual accounting remain the
frozen preflight method. Both signs use the same geometry Taylor coefficients.
The accumulated column radius includes initial interpolation/representation,
continuous time/energy residuals, centering and actual endpoint-coordinate
uncertainty. The inherited affine-projector initialization error remains
separate. Pointwise real-energy unitarity gives a uniform covariance bound
over the whole interval. No evolved metric or chosen physical history is
introduced: this is field propagation on the existing fixed geometry.

## Polynomial source, reference and exact current

For real x, `T_i T_j=(T_(i+j)+T_|i-j|)/2`. The source uses the raw polynomial
outer product, with no endpoint normalization. The exact moments
`integral_-1^1 T_n dx=0` for odd n and `2/(1-n²)` for even n integrate its
constant and energy-linear vertices; `dE=4 dx` is applied once.

The reference uses the existing generated Bloch expressions through ad4,
with the same frozen binary mass and signed angular labels. Its Gauss48
roots and directed weights are reused from the certified source-quadrature
owner. On the energy disk of center28 and radius8, the gap square has positive
real part. Thirty-two directed boundary arcs bound the analytic remainder
`4hM(h/R)^96/(1-h/R)`. Reference quadrature, interval arithmetic and final
binary source representation are included explicitly.

The true vacuum projector and the ad4 reference each have trace1, so the
vacuum current is exactly zero. This identity is used explicitly; the raw
polynomial trace/current is retained as a diagnostic with its own error
bound. No polynomial or archived column is normalized. Physical thermal is
not zero: the unchanged horizon/incoming occupation law supplies its existing
positive finite-window bound, included separately in all four moments.

## Actual replacement and complete scope

The old arrays are the authenticated `low_ref/32_1` and `low_ref/32_-1`
32-node cell on `[24,32]`. The new source minus that old source is recorded
as an explicit replacement delta. The previously certified high region starts
at32, so these updates are disjoint. The old thermal insertion leaves with
the old cell; the new vacuum approximation retains a positive bound for the
physical thermal remainder. All prior arrays and receipts remain unchanged.

The [result](../results/development/nsc-incoming-validated-energy-window.json)
binds both field artifacts, the observable artifact, complete source vector
and its total error. Other LOW and subgap intervals remain outside this
certificate. A successful cell result is not a complete source or joint
constraint solution. Source couplings, scales, seeds, occupation kappa,
physical IVs, metric state, GitHub and PDF are unchanged.

```sh
python3 scripts/derive_nsc_incoming_validated_energy_window.py --prepare --workers 2
python3 scripts/derive_nsc_incoming_validated_energy_window.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_validated_energy_window.py
```
