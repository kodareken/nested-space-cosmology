# The march crosses the inner horizon

The last commit recorded what the fix would cost. It has been done.

## What changed

The outgoing speed vanishes at both horizons, so the characteristic
right-hand side carries a simple pole at each. The code cleared the one at
`r_plus` by evaluating the spinor jet on the horizon rather than dividing by
the vanishing speed. The same construction now runs at `r_minus`.

Three things, all defaulted off:

1. `allow_below_inner_horizon` threaded through `_rhs_characteristic`,
   `_point_chart` and `vacuum_chart`, so the solver's right-hand side stops
   refusing radii below the inner horizon when it is asked to.
2. `_horizon_jet` takes a radius, defaulting to `r_plus`, and builds a
   finite jet at the inner horizon.
3. The inward branch is cut at `r_minus` and restarted from that jet, with
   nodes near each horizon filled by Taylor so no branch starts at its own
   singular point.

## It works

| lower edge | before | after |
|---|---|---|
| `1.0400` | ok | ok |
| `0.7550` | ok | ok |
| `0.7000` | step-size failure | **ok** |
| `0.6000` | step-size failure | **ok** |
| `0.5000` | step-size failure | **ok** |

`0.50` is well below the inner horizon at `0.754343`, and close to the shift's
real limit at `Q^2/(2M) = 0.4808`, which the grid still refuses correctly.

## Nothing moved for anyone who did not ask

The prepared packets for three configurations — `mu = 1.002` at centre
`12`, `mu = 1.040` at centre `12`, and `mu = 1.040` at centre `3` — hash
identically before and after the change. The combined digest is
`91b4e51dc5504dc5934feed86b63bb0f` on both sides.

The curated suite is unchanged at `368 passed, 1 skipped`.

That check is the reason this was done as its own piece of work. The
failure mode of a change to this path is not a loud crash; it is a recorded
family that differs from the last one in the eighth digit and nobody
notices for months. Two default-off flags and a branch split are only
defensible with a bit-identity check attached.

## What this does not establish

**That the third region contains anything.** The march now reaches it. That
is not the same as there being something there to find.

**That a march in it is meaningful.** The region below `r_minus` is static
again in this geometry, which is a claim about the reference solution and not
about what a sourced packet does in it. Nothing here tests that.

**That the nesting iterates.** Still one chart, one child interval.

The obvious next step is to march a packet with its edge below the inner
horizon and look at what the mass function does there, which is the first
thing in this line that would be a direct measurement in the third region
rather than a reading of the reference chart. That has not been done and is
not claimed.

## Scope

A production change to the horizon continuation. It does not measure the
actual universe.

The three threading changes are keyword-only and default to the previous
behaviour. The branch split runs only when the grid actually reaches below
`r_minus`, so the recorded family never takes it.

## Provenance

- `lab/src/recursive_horizons/nsc_rn_source.py`
- `lab/docs/nsc-rn-pole-structure.md` located the pole this clears
- `lab/tests/test_nsc_rn_third_region_reach.py`, `23` passing, including the
  default still refusing and no node being dropped by the split
