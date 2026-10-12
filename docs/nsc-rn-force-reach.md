# Effective reach of the two RN radial forces

Douglas asked two connected questions. Does the larger of the two forces
dominate a larger region and reach further, with nothing pushing back?
And are strongly compressed fields tiny, strongly localised ones rather
than long-reach ones?

This note answers both on the recorded states, and the first answer is
only partly yes.

## Domain

Areal radius on the recorded RN radial grid, out to `32 r_m`, observer the
areal chart of the recorded coupled arm, state a station payload already
written by `nsc-rn-neutral-family-v1`. Six cases, 177 stations. No time
step taken and nothing re-run.

It is a construction check on recorded finite states. A reach here is a
property of one computed collar on one grid, not a force law.

## Definition

A force's **reach** is the areal radius out to which it stays above
`REACH_FLOOR = 1e-3` of its own peak. Two deliberate choices:

* reach is taken at the **last** crossing, scanning inward, so the
  ringing near the throat cannot make a force look like it died early;
* a **supported prefix** is also recorded, the node count inside which
  the force never returns above the floor, so a later dip does not
  shorten the reported reach.

The floor is the same fraction the [balance
wall](nsc-rn-gradient-bubble.md) uses for its live region, so both
answers refer to the same supported shell.

## Does the larger force reach further?

Mostly, but not reliably.

| Case | reach radial | reach angular | longer | radial dominant on |
|---|---|---|---|---|
| `mu-1.002-eps-0.001` | 5.051 | 3.762 | radial | 96.2% of nodes |
| `mu-1.002-eps-0.010` | 5.051 | 3.746 | radial | 96.2% |
| `mu-1.010-eps-0.001` | 7.849 | 6.255 | radial | 98.6% |
| `mu-1.010-eps-0.010` | 8.700 | 5.759 | radial | 98.5% |
| `mu-1.040-eps-0.001` | 8.316 | 9.565 | **angular** | 96.8% |
| `mu-1.040-eps-0.010` | 8.816 | 6.854 | radial | 95.5% |

Amplitude dominance is unambiguous. The radial force is larger on
`95.5%` to `98.6%` of grid nodes in every case. So the answer to "does
the larger force dominate a larger area" is yes.

Reach is a different question. The larger force out-reaches in five of
six cases, but in `mu-1.040, eps=0.001` the **smaller** angular force
reaches further, `9.565` against `8.316`. And no case holds a single
verdict across all of its stations: the ordering flips between stations
in every case.

So "larger, and reaching further with nothing pushing back" is the right
picture for the weak-compression end of the family and the wrong picture
at the strong-compression end.

## Are strongly compressed fields tiny and localised?

Yes, and that is the clean relation here. Across the six cases the
throat compression slope and the radial force's reach correlate at
`r = 0.969`:

| Case | throat slope | reach radial |
|---|---|---|
| `mu-1.002-eps-0.010` | -152.64 | 5.051 |
| `mu-1.002-eps-0.001` | -121.60 | 5.051 |
| `mu-1.010-eps-0.010` | -18.54 | 8.700 |
| `mu-1.010-eps-0.001` | -13.52 | 7.849 |
| `mu-1.040-eps-0.010` | -8.44 | 8.816 |
| `mu-1.040-eps-0.001` | -7.37 | 8.316 |

Steeper compression shortens the dominant force's reach. The most
compressed case is also the shortest-reach one. That is the "strong but
local" reading, and it holds.

This is the mechanism behind the failure above. A force concentrated at
the throat is a steep force, and a steep force dies fast. So the same
property that makes the radial force dominant over most of the collar is
the property that limits its support beyond the collar.

Six points, one family, no error bar. It is a direction, not a law.

## Why reach and dominance come apart

Two reasons are visible in the recorded profiles and neither is removed
here.

The forces ring. The local logarithmic slope of the radial force swings
through `-0.17`, `-20`, `-26`, `+56`, `+2`, `-1`, `-2`, `-5`, `-16`
across one collar, so the profile is not a clean power law and any
single "slope" is a local sample of an oscillating signal. The
oscillation is characterised separately by
[nsc_rn_force_oscillation.md](nsc-rn_force_oscillation.md).

The crossovers sit far out. Every case has between 25 and 52 order
swaps between the two magnitudes, and they cluster at large radius where
both forces are near their floor. Amplitude dominance over the collar is
therefore not contradicted by the swaps; the two measures are asked at
different radii.

## Owners

The measurement is
[src/recursive_horizons/nsc_rn_force_reach.py](../src/recursive_horizons/nsc_rn_force_reach.py),
driven by
[scripts/derive_nsc_rn_force_reach.py](../scripts/derive_nsc_rn_force_reach.py),
with tests in
[tests/test_nsc_rn_force_reach.py](../tests/test_nsc_rn_force_reach.py).
The record is
[results/development/nsc-rn-force-reach-v1.json](../results/development/nsc-rn-force-reach-v1.json),
creation only. The throat slope it reuses comes from
[nsc_rn_gradient_bubble.py](../src/recursive_horizons/nsc_rn_gradient_bubble.py);
the forces themselves come from
[nsc_rn_source.py](../src/recursive_horizons/nsc_rn_source.py).

## What this does not say

It does not say gravity is a fabric, and it does not establish any force
law. It does not separate the ringing from a step-control artifact; that
is the oscillation note's job and it has not reported yet. It does not
cover `mu` above `1.04`, which was never run, and it inherits the
incomplete family: no case has a `summary.json`, and the two `mu=1.002`
cases stop silently at `t = 8`.