# The release-radius sweep: the compression holds, the relaxation does not

The single inner-packet run showed a density peak at the inner boundary
growing by `x1.197` after the shell arrived. One run cannot tell a
construction from an accident. This is the same packet, the same width of
`2 r_m`, the same spectrum, the same chart, the same `mu = 1.04` and
`eps = 0.001`, the same `869` points, the same inner edge at `1.04`. Only
the release radius moves.

## What holds

| release | stations | falls | arrives | formation `r` | rate | `r^2` | **gain** | **rise** |
|---|---|---|---|---|---|---|---|---|
| `2.0` | 21 | no | `0.0` | pinned | — | — | `x1.1543` | `0.70` |
| `3.0` | 21 | 6 | `0.6` | `2.8234` | `-1.8446` | `0.9968` | `x1.1972` | `0.80` |
| `4.5` | 36 | 16 | `1.6` | `4.4285` | `-1.6565` | `0.9985` | `x1.1875` | `0.80` |

**The compression holds.** The gain is `1.154` to `1.197` across release
radii spanning more than a factor of two, a relative spread of `3.6%`. The
density peak at the inner boundary grows by about a fifth every time.

**The rise time holds.** `0.70` to `0.80` time units after arrival. The peak
comes at the same lag regardless of how far the shell fell to get there.

Together these are the first two quantities in this programme that survive a
deliberate change of a construction parameter.

## What does not hold

**The relaxation after the peak does not settle.** The final amplitude as a
fraction of the peak is `0.648`, `0.921` and `0.738`. Those track how long
each march continued past its peak, not the release radius. These runs do
not settle it and **no relaxed value is claimed**. The earlier note's
"slowly relaxes, sitting at 92% of its peak" was a statement about one run's
last station and read as a general result. It is not.

**The width compression is not a fixed factor.** It grows with release
radius: `123` to `100` cells for `3.0`, `136` to `99` for `4.5`. It is a
property of how far the shell fell, not a constant. The shell that falls
further has further to compress.

**The infall rate has only two points.** `-1.8446` at release `3.0` and
`-1.6565` at `4.5`: the deeper shell falls faster. Two points is not a law
and no exponent is fitted.

## The short release is degenerate, and the guard was wrong

At release `2.0` the run has **zero falling stations**. The packet is `2 r_m`
wide, so centred at `2.0` its left edge is at `1.0`, already below the
inner edge at `1.04`. The shell is born pinned to the boundary; there is no
fall to measure.

My guard checked the *centre* against `1.5 * r_inner` and let `2.0` through.
The binding condition is the packet's **left edge**, not its centre. The
march script now refuses a packet whose left edge is at or below the inner
edge, and says why.

The `2.0` run is kept in the record and flagged `starts_pinned_at_inner_edge`.
It still measures the accumulation, and its gain of `1.154` is the smallest
of the three, consistent with the others. But it measures nothing about a
shell moving, and the record does not let a reader mistake one for the other.

## Resolution limits, unchanged

The peak radius stays quantised to the grid, so per-interval speeds remain
unusable and the fitted rate is a mean over the phase. Widths are half-max
widths on a fifteen-cell median. Nothing here resolves what the compressed
mass becomes; it only measures how much it compresses against a boundary
while it is inside the domain.

## Scope

A construction check on finite marched collars. It does not measure the
actual universe. Holding a result across a change of one construction
parameter shows it is not an accident. It does not make it physical.

Three runs, one mass ratio, one source strength, one packet width. The
sweep varies the release radius and nothing else.

## Provenance

- marches: `lab/scripts/run_nsc_rn_inner_packet.py`, runs
  `nsc-rn-inner-packet-c2`, `nsc-rn-inner-packet-v1`, `nsc-rn-inner-packet-c4p5`
- reduction: `lab/scripts/derive_nsc_rn_packet_sweep.py`
- record: `lab/results/development/nsc-rn-packet-sweep-v2.json`
- accumulation logic shared with the single-run reduction through
  `nsc_rn_mass_shell.accumulation_after_arrival`
- tests: `lab/tests/test_nsc_rn_packet_sweep.py`, `20` passing

## What would settle the relaxation

The relaxation is unresolved because each march ended at an arbitrary
station rather than at a fixed lag after the peak. Marching every run to,
say, twice the rise time past its own peak would put them on a common clock
and settle whether the `0.65` to `0.92` spread is a transient or a genuine
dependence. That is a small change to the march target and is the obvious
next step.