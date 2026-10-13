# A map of the RN radial family

This line has been chasing a map rather than a discovery: take one space,
describe its structure completely, and see how the structure moves as the
one free parameter changes.

Every quantity here is classical. The geometry is a charged spherical
solution in areal coordinates. The quantities are its two horizons, the two
nesting cuts the chart defines, the shift, and the causal membrane located
from the metric. Nothing quantum is involved and nothing exotic is assumed.

Eleven values of the mass ratio, `mu` from `1.0005` to `1.100`.

| mu | r_minus | excision | r_plus | causal membrane | child start | child end |
|---|---|---|---|---|---|---|
| 1.0005 | 0.9689 | 1.0005 | 1.0321 | 1.0321 | 1.0361 | 1.0638 |
| 1.0020 | 0.9387 | 1.0020 | 1.0653 | 1.0653 | 1.0732 | 1.1286 |
| 1.0100 | 0.8682 | 1.0100 | 1.1518 | 1.1518 | 1.1695 | 1.2935 |
| 1.0400 | 0.7543 | 1.0400 | 1.3257 | 1.3257 | 1.3614 | 1.6113 |
| 1.1000 | 0.6417 | 1.1000 | 1.5583 | 1.5583 | 1.6155 | 2.0165 |

## The ordering holds everywhere

```
r_minus  <  excision  <  r_plus  =  causal membrane  <  child start  <  child end
```

Every feature keeps its place across the whole map. Nothing crosses anything
else. **The nesting has no transition in it** over this range, which is the
first thing a map is for.

The causal membrane sits on `r_plus` at every value, agreeing to within
`0.0095` of a grid cell. So the one real interface in this geometry is the
outer horizon, and it is there for classical reasons: that is where the
outgoing characteristic speed changes sign.

## Three regions, not two

The shift equals one at **both** horizons, so the outgoing speed vanishes at
both and the chart divides in three. For `mu = 1.04`:

| region | interval | outgoing speed | trapped |
|---|---|---|---|
| parent | `[1.3257, 32]` | positive | no |
| **trapped** | `[0.7543, 1.3257]` | **negative** | **yes** |
| below the inner horizon | `[0, 0.7543]` | positive | no |

Only the outermost zero was a sign change inside the original chart's
guarded domain, because that interface stopped at `r_minus`. The later
default-off horizon continuation below distinguishes that guard from the
reference geometry's real PG domain.

**Every inner edge this programme has ever used lies inside the trapped
region** — `1.04`, `0.85`, `0.80`, `0.77`, `0.76`, `0.7555`, all of them
between `0.7543` and `1.3257`. So the whole compression story is the
dynamics of matter that is already behind the membrane: the `x2.05`, the
three-sided shell, the `94%` retention. That is not an accident of where the
grid was placed; it is where the trapped region is.

The third region is a reference-geometric candidate, not an established
next nest. The original guard refused samples below `r_minus`; the later
continuation can evaluate them when explicitly enabled. A sourced evolution
there still needs its own boundary construction.

## Real PG domain and the inner-horizon continuation

The reference shift `beta = sqrt(2M/r-Q^2/r^2)` remains real for
`r > Q^2/(2M)`, below the inner horizon. At `mu=1.04`, the latter is
`0.754343`, while the real-shift limit is `0.4808`. Both horizons have
outgoing characteristic speed zero; the speed is positive again below
the inner horizon within this real domain. This is a description of one
charged reference solution, without a third-region occupation claim.

The inward characteristic equation has a simple pole where the outgoing
speed vanishes. The retained implementation uses the same finite horizon-jet
construction at `r_minus` as at `r_plus`, splits and restarts the branch,
and fills nearby nodes by Taylor jets. `allow_below_inner_horizon` is
keyword-only and default off; the real-shift guard remains. The recorded
checks reached lower edges `0.70`, `0.60`, `0.50` after the former step-size
failures. Three unchanged packet preparations (`mu=1.002` at centre `12`,
`mu=1.040` at centres `12` and `3`) retained combined digest
`91b4e51dc5504dc5934feed86b63bb0f`. The historical curated result was
`368 passed, 1 skipped`; that count is provenance, not a current suite run.
Crossing the reference horizon proves neither a valid sourced march beyond
pure-outflow excision nor iterated nesting.

## The one structural fact the map exposes

**The shift has a maximum at `r = Q²/M`, and the default excision always
sits above it.**

| | |
|---|---|
| shift maximum | `r = 1/mu` |
| default excision | `r = mu` |
| ratio | `0.90` to `0.998`, always below one |

For every entry on the map the shift's turning point lies *inside* the
default grid's inner edge. Every run this programme has ever done has begun
above that turn and has therefore never sampled the region where the shift
is largest.

The deep-edge runs added later are the only ones that cross it, and they are
where the three-sided shell appears. That is not a coincidence worth
ignoring.

## What the map cannot show

**Whether the nesting iterates.** This is one chart with one child interval
and one parent. There is no second nest for anything to flow between. A map
of one space is not a map of nested spaces, and saying otherwise would be
the easy error.

**Whether the ordering persists outside this family.** It holds for eleven
values of one parameter in one construction. That is what was asked.

## Scope

A map of one family. **It does not measure the actual universe.** What it
can say is whether the ordering of the features is stable as the parameter
moves, and it is.

Reference chart only: no station payloads, no time step, no marching.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_geometry_map.py`
- records: `nsc-rn-geometry-map-v4.json` is current. `v1` used an arbitrary
  `1e-6` tolerance on the ordering test and reported the ordering failing for
  nine of eleven entries; `v2` and `v3` carry the located zero's deficit
  before the tolerance was scaled to the grid. All are retained.
- tests: `lab/tests/test_nsc_rn_geometry_map.py`, `16` passing
