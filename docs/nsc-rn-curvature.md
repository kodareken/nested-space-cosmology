# Does the compressing mass make the geometry singular?

Mass piling up at a boundary is a mass statement. The geometry question is
separate and sharper: does the curvature of the *evolved* metric blow up
there, or stay finite?

This reads the curvature reconstructed from the evolved lapse and shift and
their jets, at each station of a `3.0 r_m` release inner-packet march. That
matters: the vacuum reference on the grid does not know the shell is there,
so it could not answer this even in principle. The recorded family's `R4`
comes from this same path.

All `31` stations carried a reading, and `R4` and `Ricci2` are finite
throughout.

## What the edge curvature does

| | value | when |
|---|---|---|
| initial | `3.1030e-04` | `t = 0`, before the shell arrives |
| **maximum** | `8.3271e-04` | `t = 0.1`, `x2.68` of initial |
| at the mass peak | `1.5108e-04` | `t = 1.4`, the same station as the mass peak |
| final | `4.5076e-06` | `t = 3.0`, `x0.0145` of initial |

Three features, and the first-versus-last comparison hides all of them:

- a **transient maximum** at `t = 0.1`, the packet forming;
- a **second and lower maximum** at `t = 1.4`, which lands on exactly the
  station where the mass density peaks. The curvature responds to the
  compression, by about half the size of the formation transient;
- a **relaxation** of `x0.0145`, so the edge ends about `69` times lower
  than it started, and about `35` times lower than its compression maximum.

The edge layer also stops being sharply confined: the ratio of the edge
value to its immediate neighbour falls from `12.15` to `5.60`.

**The geometry stays finite while the mass compresses.** Nothing diverges at
the excision.

## Why the obvious metric is the wrong one

The first version of this reduction compared the edge value against the bulk
maximum. On this family that ratio grows from `320` to `122`, which looks
like a curvature catastrophe.

It is not. The bulk maximum *decays* over the run while the edge value
barely moves, so the ratio is measuring the interior going quiet, not the
boundary going singular. The record reports the edge on its own terms, its
ratio to the immediate neighbour, and the bulk ratio only with that caveat
attached.

This is the same mistake the shell reduction made twice already: a statistic
that moves for two reasons at once, read as if it moved for one.

## Scope

A construction check on one finite marched collar. It does not measure the
actual universe.

**Regularity at the inner edge of this collar is not regularity at a
horizon.** The excision is the horizon midpoint of the reference chart,
`0.5 * (r_minus + r_plus) = 1.04`. It is a boundary, not an event horizon.
Regularity at a boundary is a far weaker statement than regularity at a
horizon, and nothing here extrapolates through that edge to the inner
horizon at `r_minus = 0.754`.

So this does not settle the regularity question the family was opened to
address. What it does establish is narrower and still worth having: the
accumulation of mass against the excision, which the recorded family could
not see at all, does not make the reconstructed curvature at that boundary
diverge. The mass compresses, the curvature answers, and both relax.

## Provenance

- march: `lab/scripts/run_nsc_rn_inner_packet.py`, now recording `R4` and
  `Ricci2` reconstructed from the evolved metric
- reduction: `lab/scripts/derive_nsc_rn_curvature.py`
- records: `nsc-rn-curvature-v1.json`, whose verdict was built from a
  first-versus-last test that hid the transient maximum, and
  `nsc-rn-curvature-v2.json`, which is current
- bulk margin: `10` nodes at each end, matching the recorded family's choice

## What would settle it

The excision has to move. Regularity *at a horizon* needs the chart to
extend past it, and the inner edge is currently derived rather than free:
`build_grid` sets it to `0.5 * (r_minus + r_plus)`. Putting the grid's edge
inboard of that, between `r_minus = 0.754` and `1.04`, is a real change to
the chart's boundary condition and to the child/parent cut ordering, and it
must be added alongside `build_grid` rather than inside it so the recorded
family is untouched. That is the construction that would let this question
be asked properly, and it is the next thing to build.