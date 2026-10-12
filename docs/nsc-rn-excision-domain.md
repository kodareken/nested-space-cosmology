# Where this construction's boundary condition makes sense

The march now reaches below the inner horizon. That turns out not to be the
interesting question. The interesting question is whether a *sourced* march
is admissible there, and the answer is no — for a reason that is not a
convention.

## The fifth gate, and it is a different kind

`admit_chart` requires the excision to be a pure outflow:

```
outgoing(edge) < 0  and  ingoing(edge) < 0
```

Measured across the chart, for `mu = 1.04`:

| edge | outgoing | ingoing | pure outflow? |
|---|---|---|---|
| `1.0400` | `-0.03704` | `-2.03704` | yes |
| `0.9000` | `-0.03757` | `-2.03757` | yes |
| `0.8000` | `-0.01858` | `-2.01858` | yes |
| `0.7550` | `-0.00033` | `-2.00033` | yes |
| `0.7400` | `+0.00770` | `-1.99230` | **no** |
| `0.6000` | `+0.17001` | `-1.82999` | **no** |
| `0.5000` | `+0.60000` | `-1.40000` | **no** |

The outgoing speed changes sign at `r_minus = 0.754343`. Above it the
excision absorbs; below it the excision **emits**. The condition holds
strictly between the two horizons.

## This is the answer, not a blockage

The first four gates were choices, and opening them was right. This one is
different. Below the inner horizon the boundary's causal role reverses:
an excision that can emit is not the same boundary condition as one that
cannot, and inventing the second from the first is not a matter of flipping
a guard.

So the natural domain of this construction is not a choice either. **It maps
exactly onto the trapped region** — which is where every measurement in this
line was taken, from the `x2.05` compression to the three-sided shell to the
`2.5%` floor.

## What this settles and what it does not

**Settles:** where the sourced march is admissible, and why. It is admissible
on the trapped region and nowhere else, for a reason internal to the
geometry rather than internal to the code.

**Does not settle:** whether the third region contains anything. That remains
open, and it is now open for the right reason: not because a guard was in the
way, but because there is no boundary condition here for it, and inventing
one would be a new construction rather than an extension of this one.

Whether such a construction is wanted is a question about the programme, not
about the code.

## Scope

A measurement on the reference chart. **It does not measure the actual
universe.**

One mass ratio, the edge swept across the whole real domain of the shift.
No station payloads, no marching.

## Provenance

- measured on the reference chart through `nsc_rn_source._point_chart` with
  the below-horizon flag
- `lab/scripts/run_nsc_rn_inner_packet.py` gained `--below-horizon`, which
  threads the flag to both `build_grid` and `prepare_radial_packet`, and
  refuses an edge below `r_minus` unless that flag is given
