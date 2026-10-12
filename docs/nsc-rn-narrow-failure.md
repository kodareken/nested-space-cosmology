# Narrowing the domain does not work, and the reason is instructive

The saturation at `x1.58` was established on one grid. The note that reported
it said the open question was whether the approach is a genuine asymptote or a
resolution artefact, and that answering it needs more points near the edge.

There is an obvious cheaper route: the far field is nearly empty, so narrowing
the outer radius buys resolution near the inner edge for free at the same point
count. Three runs at `r_out = 6` instead of `32`, at the three edges that are
resolvable on the wide grid, everything else held fixed.

**It does not work, and the reason is the finding.**

## What the narrow domain claims

| edge | cells at `r_out = 32` | cells at `r_out = 6` | gain at `r_out = 32` | gain at `r_out = 6` |
|---|---|---|---|---|
| `1.04` | `8.009` | `49.990` | `x1.1972` | `x2.4011` |
| `0.85` | `2.666` | `16.122` | `x1.3905` | `x2.9261` |
| `0.80` | `1.270` | `7.621` | `x1.4328` | `x3.6167` |

The resolution gain is real: six times finer, and all three edges now have
between seven and fifty cells between the edge and the inner horizon instead
of between one and eight. That was the plan and it worked.

Two to two and a half times larger, and still rising with depth: no
saturation at all. The boundary curvature reaches `5.08e-01` at `t = 3.9`,
which is **x310** of its initial value and is order one rather than a
perturbation. The verdict from the curvature reduction is
`edge_curvature_grows_beyond_any_settled_scale`.

That reads like a stronger version of the same physics. It is not.

## The narrow runs are unstable

| run | finite | reversals below `r = 16` | far-field reversals | lapse | stable |
|---|---|---|---|---|---|
| `r_out = 6`, edge `1.04` | yes | **15,268** | 0 | ok | **no** |
| `r_out = 6`, edge `0.85` | yes | **15,095** | 0 | ok | **no** |
| `r_out = 6`, edge `0.80` | yes | **15,036** | 0 | ok | **no** |

The wide-domain runs have **zero** reversals below `r = 16` across all six
edges. The narrow runs have about `15,000` each, and they are spread evenly
over the whole domain rather than sitting anywhere in particular:

| radius band | reversals |
|---|---|
| `0.80 - 1.44` | `1,586` |
| `1.44 - 2.09` | `1,828` |
| `2.09 - 2.73` | `1,975` |
| `2.73 - 3.38` | `2,179` |
| `3.38 - 4.02` | `2,333` |
| `4.02 - 4.67` | `2,343` |
| `4.67 - 5.31` | `1,914` |
| `5.31 - 5.96` | `878` |

**`28.9%` of all node-steps reverse.** Every field is finite and the lapse
stays positive, so the finiteness checks pass; it is the monotonicity that
fails, and it fails everywhere.

So the large gains and the order-one curvature are **not trustworthy**. They
are what an unstable march produces. The stability audit exists precisely to
catch this, and it did.

## What survives, and what does not

**Withdrawn**: every number in `nsc-rn-edge-sweep-narrow-v1.json`,
`nsc-rn-inner-packet-shell-n6-*-v1.json` and
`nsc-rn-curvature-n6-exc080-v1.json`. They are retained and registered, marked
unstable, because a route that looks promising and is not is worth keeping.

**Unchanged**: the `r_out = 32` saturation at `x1.58`, which stands as a
measurement of the construction this programme actually uses. All six wide
runs are stable, finite, and monotone below `r = 16`.

**Still open**: whether `x1.58` is a genuine asymptote or a resolution
artefact. The route that would have answered it is now closed.

## Why the wide domain is the load-bearing one

The wide domain's far field is not decoration. It is where the recorded
oscillation lives, and removing it removes the region that absorbs the
boundary interaction. The empty far field carries almost no mass, but it
carries the room the march needs to stay well-posed.

That is the same thing this line learned about the recorded family, where the
far field at `r >= 16` was dismissed as ringing on negligible amplitude. It
was negligible in mass and apparently load-bearing in numerics.

## Scope

Construction checks on finite marched collars. **This does not measure the
actual universe.**

Both variants use boundaries, not horizons, and the wide-domain edge at
`0.7555` sits `0.032` cells from the inner horizon at `0.754`.

Three edges, one release radius, one mass ratio, one source strength, one
grid resolution, one outer radius — and the outer radius is now known to
matter, which is itself a limit on how converged anything here is.

## Provenance

- marches: `run_nsc_rn_inner_packet.py --r-out 6.0`, three runs at
  `--target 6.0`, `61` stations each
- reductions: edge sweep, shell, curvature, stability, all on the narrow runs
- a first cut of the narrow sweep hardcoded `r_out = 32` when computing the
  spacing and cell counts, so it reported the wide grid's resolution for the
  narrow runs. `v1` of that record is retained with the wrong cells;
  `nsc-rn-edge-sweep-narrow-v2.json` is current and reads the extent from
  each run.
- records: `nsc-rn-edge-sweep-narrow-v1.json` (wrong cells, retained),
  `nsc-rn-edge-stability-narrow-v1.json`,
  `nsc-rn-curvature-n6-exc080-v1.json`, shell records
  `nsc-rn-inner-packet-shell-n6-*-v1.json`

## What would settle it

Not a narrower domain. To resolve the edge properly the grid needs more
points over the same extent, which is the constraint that stopped the
recorded family at `3921`. Four cells between `0.755` and the horizon at
`r_out = 32` needs roughly `8000` points, and the recorded family's memory
failures were at that scale.

So the question is now squarely a resources question, and it should be
treated as one rather than worked around again.