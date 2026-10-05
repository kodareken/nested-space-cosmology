# Does time run one way, nest to nest, and never back?

A one-way flow is a claim about causal structure, not about a profile. This
asks it directly from the recorded states, using the characteristic speeds in
the evolved metric.

The test has three parts, and all three must hold:

| test | statement |
|---|---|
| **sealed** | the outgoing speed is strictly negative inside the innermost zero |
| **stable** | the number of zeros does not change over time |
| **one-way in time** | the located zero does not reverse |

## What the recorded states say

**Every run passes all three.**

| run | zeros | sealed | one-way |
|---|---|---|---|
| recorded `mu=1.002`, `eps=0.001` | 1 | yes | **yes** |
| recorded `mu=1.002`, `eps=0.010` | 1 | yes | **yes** |
| recorded `mu=1.010`, `eps=0.001` | 1 | yes | **yes** |
| recorded `mu=1.010`, `eps=0.010` | 1 | yes | **yes** |
| recorded `mu=1.040`, `eps=0.001` | 1 | yes | **yes** |
| recorded `mu=1.040`, `eps=0.010` | 1 | yes | **yes** |
| inner packet, default edge | 1 | yes | **yes** |
| inner packet, edge `0.7555` | 1 | yes | **yes** |

There is exactly one zero of the outgoing speed in every station of every
run. Everything inside it has negative outgoing speed. The count never
changes, so no sealed region opens and none closes.

In the recorded family the zero sits exactly at `r_plus` and never moves:
`1.0653` for `mu = 1.002`, `1.1518` for `1.010`, `1.3257` for `1.040`. The
ingoing characteristic never vanishes, so the zero stays simple.

## Where the trapped region sits

This was not expected and matters for the nesting picture. The zero is at
`r_plus = 1.3257`. The chart's labelled child interval is
`[child_start, child_end] = [1.3614, 1.6113]`, which is **outside** it. The
region that is causally sealed is the one below `1.3257`, which the chart
leaves unnamed.

Every inner-edge this programme has ever placed, from the horizon midpoint at
`1.04` down to `0.7555`, is inside the trapped region, with the outgoing
speed already negative at the boundary. **The compression measured in this
line is compression of matter that cannot escape.**

## What is not resolved

The horizon does move in the inner-packet runs, from `1.3259` to `1.3302`. That
excursion is `4.2e-03`, which is **`0.12` grid cells** on the `869`-point
grids those runs used. The grid cannot resolve it, and no claim is made about
the horizon's dynamics.

There is a small inward step of `1.6e-06` near the end of those runs, `0.0001`
cells. The first version of this test flagged it as a reversal because the
tolerance was an absolute `1e-9`, which treats sub-resolution variation as a
change of direction. The test is now scaled by the grid spacing, and reports
every motion in cells so the reader can judge resolution directly rather than
take a verdict on it.

## What this cannot answer

**Whether the nesting iterates.** This construction defines one child interval
and one parent. Nest-to-nest-to-nest can be tested here only as a single step,
from parent into child, and that step is one-way. Whether there is a further
nest inside that one is a property of the construction, not of these states,
and no claim about a hierarchy is made.

**Whether this is a global event horizon.** A simple zero of the outgoing
characteristic on a sampled slice is the standard numerical proxy for a
trapping surface. A sequence of such slices is not a proof of a global event
horizon and is not offered as one.

## Scope

Construction checks on recorded finite states. **This does not measure the
actual universe.** A one-way causal ordering that survives every recorded
case in one numerical construction is evidence about that construction.

Eight runs, six of them the unmodified recorded family, one mass-ratio sweep,
two inner-packet depths, two source strengths.

## Provenance

- criterion: `nsc_rn_observables.trapping_horizon`, an outermost simple zero
  of the outgoing speed in the evolved metric, never copied from the
  reference horizon
- reductions: `lab/scripts/derive_nsc_rn_trapping_horizon.py`,
  `lab/scripts/derive_nsc_rn_causal_order.py`
- records: `nsc-rn-trapping-horizon-v1.json`, `nsc-rn-causal-order-v2.json`,
  `nsc-rn-causal-order-v3.json`; `v1` and `v2` are retained, `v2` having used
  the absolute reversal tolerance
- the inner edge is reported separately and is never called a horizon

## What would settle it

The `8000`-point runs now marching at the deep edges resolve the edge region
six times better, which turns the `0.12`-cell horizon excursion into roughly
`0.7` cells. That is the resolution needed to say whether the horizon's motion
is monotone or has structure, which is the difference between a fixed
membrane and a living one.