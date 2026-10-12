# What the compression is bounded by

The saturation of the compression gain has been described as *deepening*:
further in, more compression, and it stops near the inner horizon. That is
true, and it says where the saturation happens without saying why.

Read the gain against the two distances that define the region the
compression can actually occupy — the distance to the horizon the matter is
sealed behind, and the distance to the inner horizon, past which the
excision stops being an excision.

| `r_inner` | to horizon | to `r_minus` | % of trapped width | gain |
|---|---|---|---|---|
| `1.0400` | `0.2857` | `0.2857` | **50.0%** | `1.2571` |
| `0.8500` | `0.4757` | `0.0957` | `83.3%` | `1.6039` |
| `0.8000` | `0.5257` | `0.0457` | `92.0%` | `1.8391` |
| `0.7700` | `0.5557` | `0.0157` | `97.3%` | **`2.0481`** |
| `0.7600` | `0.5657` | `0.0057` | `99.0%` | `2.0295` |
| `0.7555` | `0.5702` | `0.0012` | `99.8%` | `1.9952` |

The gain rises monotonically as the boundary recedes from the horizon, peaks
at `97.3%` of the trapped width, and turns over at `99.0%` and `99.8%`.

## The saturation and the boundary condition end in the same place

The deepest boundary tried sits `0.0012` from `r_minus`. At that radius the
excision boundary condition reverses: the outgoing speed changes sign, the
excision stops absorbing and starts emitting, and the pure-outflow condition
that every run in this line relies on stops holding.

So the compression is bounded by **the room left inside the trapped region**,
and the run that came closest to the deepest edge is the run that came closest
to losing the boundary condition altogether.

That is a different statement from resolution saturation, and it is the
better one: it says the limit is structural, not numerical. Pushing the
boundary deeper cannot raise the gain, because past this point there is
nowhere for the matter to be compressed *into* while the construction still
admits an excision.

## What it rules out

It rules out the reading that `x2.05` is merely what a `869`-point grid
happens to produce. The refinement check moved the gain from `1.58` to `2.05`,
and that was real. But the ceiling itself is set by geometry, and it sits
where the boundary condition does.

## What it does not say

It does not say the compression is maximal in the universe. It says it is
maximal **for an excision-based construction on this chart**, at about `97%`
of the trapped width. A construction with a different boundary condition, or
one that does not excise at all, is not bounded here.

## Scope

A re-reading of records already taken. **It adds no new run and does not
measure the actual universe.** Six edges, one mass ratio, one source
strength.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_compression_bounding.py`
- record: `lab/results/development/nsc-rn-compression-bounding-v2.json`; v1
  copied the gains by hand and is retained, with the risk noted below
- gains read from the fine-resolution shell records
  `nsc-rn-inner-packet-shell-n8000*-v1.json`, horizons from
  `nsc-rn-geometry-map-v6.json`, the sign change from
  `nsc-rn-excision-domain.md`

The first version of this reduction copied the six gains by hand and then
compared them against the coarse edge sweep. Those differ by up to thirty
percent, so the comparison would have been meaningless while still looking
like a clean number. The values are now read from the fine records, and a
test checks that each one exceeds the coarse value rather than matching it.
