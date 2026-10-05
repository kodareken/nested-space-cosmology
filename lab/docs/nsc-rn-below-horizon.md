# The third region is real, and the guard that hides it is a convention

The map found three regions rather than two. The third sits below the inner
horizon, and `vacuum_chart` refuses to sample it:

> `static modes stay strictly above the inner horizon`

That is the region a further nest would have to be in. It is worth examining
rather than accepting.

## The guard reaches further in than the geometry does

The shift is `beta = sqrt(2M/r - Q^2/r^2)`, which is real wherever
`r > Q^2/(2M)`. The chart's guard stops at `r_minus`. Those are not the same
place:

| mu | chart guard | real domain of the shift |
|---|---|---|
| 1.002 | `0.9387` | `0.4990` |
| 1.010 | `0.8682` | `0.4950` |
| 1.040 | `0.7543` | `0.4808` |

The guard sits at roughly **three quarters of the way** between the real
limit and the outer horizon. The horizons are where the shift passes one;
they are not where it stops being real.

## The model evaluates there perfectly well

Calling the model's own functions directly, below `r_minus`:

| r | lapse | shift | outgoing | ingoing |
|---|---|---|---|---|
| `0.7400` | `1.000000` | `0.992301` | `+0.007699` | `−1.992301` |
| `0.6000` | `1.000000` | `0.829993` | `+0.170007` | `−1.829993` |
| `0.5000` | `1.000000` | `0.400000` | `+0.600000` | `−1.400000` |

Finite, real, and the same light constant `outgoing + |ingoing| = 2`. Nothing
breaks. The lapse is one throughout.

## The outgoing speed crosses zero twice

| mu | crossings |
|---|---|
| 1.002 | `0.93872`, `1.06528` |
| 1.010 | `0.86823`, `1.15177` |
| 1.040 | `0.75434`, `1.32566` |

Both crossings are at a horizon. So the trapped region is bounded on both
sides, and **below the inner horizon the outgoing speed is positive again**:
matter there can move outward.

Only the outermost crossing is a sign change inside the chart's own domain,
because the chart stops at `r_minus` and so never sees the second one. That
is why every reduction in this line reported exactly one zero, and it was
correct on its domain — the domain just excluded the other one.

## Nothing was relaxed

This commit changes no production guard. Showing that the third region is
describable says nothing about whether it contains anything, and nothing
about whether the nesting iterates. Whether a march can be run there is a
separate question with a separate answer, and changing `vacuum_chart` on the
strength of a reading would be the wrong order of operations.

It does establish that the obstacle is a choice in one function rather than
a property of the geometry. That is worth knowing before spending anything
on the alternative.

## What this still does not say

The third region is reachable in principle. Nothing here says what is in it,
whether it is static in any useful sense, whether matter that crosses the
inner horizon outward does anything interesting, or whether it is the next
nest. All four are open.

It also does not touch the cosmological questions. This is one charged
geometry with one wave packet and three regions.

## Scope

A construction check on the reference model. **It does not measure the
actual universe.** No station payloads, no time step, no marching, no
production guard changed.

Three values of the mass ratio, `20000` samples per entry.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_below_horizon.py`, which calls the
  model's functions directly rather than through the chart wrapper
- record: `lab/results/development/nsc-rn-below-horizon-v1.json`
- tests: `lab/tests/test_nsc_rn_below_horizon.py`