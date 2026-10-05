# The gain saturates at about x1.58

The previous note said the compression gain was "still rising with depth" and
that the trend "has not turned over even now that peak and relaxation are
both resolved". **That was wrong.** It has turned over, at about `x1.58`.

Two new runs at `0.76` and `0.7555`, and two re-runs of the shallow baselines
to `t = 6.0` so all six edges are compared on the same clock:

| `r_inner` | gap to `r_minus` | **gap in cells** | resolvable | **gain** | increment |
|---|---|---|---|---|---|
| `1.04` | `0.2857` | `8.009` | yes | `x1.1972` | — |
| `0.85` | `0.0957` | `2.666` | yes | `x1.3905` | `+0.1933` |
| `0.80` | `0.0457` | `1.270` | yes | `x1.4328` | `+0.0423` |
| `0.77` | `0.0157` | `0.435` | **no** | `x1.5758` | `+0.1430` |
| `0.76` | `0.0057` | `0.157` | **no** | `x1.5825` | `+0.0067` |
| `0.7555` | `0.0012` | `0.032` | **no** | `x1.5849` | `+0.0024` |

The last two increments are thirty and sixty times smaller than the first.
The gain has flattened.

## What the wrong claim got wrong

I read the trend off depth in `r_m`, where the gaps look meaningfully
different: `0.0157`, `0.0057`, `0.0012`. In **grid cells** they are `0.435`,
`0.157` and `0.032` — all less than half a cell. Those runs are not "closer
to the horizon"; on this grid they are the same place, and the small
increment between them is not a depth effect at all.

So the sweep was **resolution-limited before it was depth-limited**, which is
a different statement from "not converged". The correct reading:

- across the **resolvable** range, `8.0` down to `1.27` cells, the gain rises
  monotonically from `x1.20` to `x1.43`;
- below one cell the edge is indistinguishable from the horizon, and the gain
  sits at `x1.58`.

The edge sweep driver now reports the gap in cells and flags which edges are
resolvable, because depth in `r_m` alone overstates how far the sweep went.

## So the midpoint number is a lower bound by x1.32

The `x1.20` reported at the horizon midpoint is not the effect. Deepening the
domain by `0.28 r_m` raises it to `x1.58`. Every headline number measured at
the midpoint excision is truncated by about a third.

## What does not depend on the edge at all

| `r_inner` | `1.04` | `0.85` | `0.80` | `0.77` | `0.76` | `0.7555` |
|---|---|---|---|---|---|---|
| fraction of gain below the edge | `0.943` | `0.942` | `0.942` | `0.941` | `0.941` | `0.941` |

Flat to three figures across a factor-of-two change in depth. The deep domain
**retains** the mass at every edge: `94%` of the accreted mass ends below it.

The relaxation, by contrast, does not organise with depth at all:
`0.167`, `0.490`, `0.144`, `0.635`, `0.527`, `0.372`. It scatters between
`0.14` and `0.64` with no ordering. That is not a trend and none is claimed;
it is a quantity that needs its own study rather than being read off a depth
axis.

## A guard that fired correctly

The first attempt at this sweep was refused: it mixed `t = 3.0` and `t = 6.0`
runs and said `runs differ in target_over_rm: [3.0, 3.0, 6.0, 6.0, 6.0,
6.0]. This is not an edge sweep.` The two shallow baselines were re-marched
to `t = 6.0` rather than the guard being relaxed, which is why all six rows
above are on the same clock.

## Scope

Construction checks on finite marched collars. **This does not measure the
actual universe.**

Every edge here is still a **boundary**, and three of the six sit less than
half a grid cell from the inner horizon at `0.754`, which makes them
indistinguishable from it on this grid but does not make them horizons. No
reading here is a statement about regularity at a horizon.

Six edge positions, one release radius, one mass ratio, one source strength,
one grid resolution. The `x1.58` is the saturation of this construction on
this grid, not a property of anything else.

## Provenance

- marches: `run_nsc_rn_inner_packet.py`, all six at `--target 6.0`,
  `61` stations each
- reduction: `derive_nsc_rn_edge_sweep.py --runs`, which now reports
  `gap_in_cells` and `edge_resolvable_from_horizon`
- records: `nsc-rn-edge-sweep-v2.json`, superseding `v1`; shell records
  `nsc-rn-inner-packet-shell-exc0760-long-v1.json` and
  `-exc07555-long-v1.json`

## What would settle it

The saturation is established at this resolution and cannot be pushed further
by moving the edge. To resolve a genuinely deeper boundary the grid needs
more points near the edge, which is the constraint that stopped the recorded
family at `3921`; four cells between `0.755` and the horizon at the present
outer radius would need roughly `8000` points.

Whether the gain's approach to `x1.58` is a genuine asymptote or a
resolution artefact is therefore still open, and the answer is a resolution
question, not a depth one.