# Resolution control, and the relaxation on a common clock

Two open questions from the release-radius sweep, and a third that appeared
while closing them.

## The compression is not a resolution artifact

The growing peak sits at the grid's left edge, and the left edge is exactly
where the one-sided SBP closure lives. A growth of about a fifth measured
while pinned there could be the closure accumulating rather than the
geometry. It had to be checked before it could be claimed.

Same release radius of `3.0`, same packet, same `mu = 1.04` and
`eps = 0.001`, marched at two resolutions:

| points | gain | arrival | rise |
|---|---|---|---|
| `869` | `x1.1972` | `0.6` | `0.80` |
| `1740` | `x1.2285` | `0.6` | `0.90` |

**Relative difference `2.55%`.** The arrival time is identical to the
station. The gain moves by a couple of percent when the grid is doubled.

That is this line's **first same-parameter, different-resolution replay**.
Earlier notes repeatedly recorded that no such replay existed, which had left
every near-throat quantity in this programme without a convergence control.

Two resolutions bracket nothing. If the gain had moved by more than a few
percent the single-resolution result would be unusable; that it does not
means the result is merely not obviously an artifact. No convergence order
is claimed and no Richardson extrapolation is attempted. Agreement at two
resolutions is evidence against a resolution artifact, not evidence of
physics.

## The relaxation, on a common clock

The first sweep ended each run at an arbitrary station, so its after-peak
fractions of `0.648`, `0.921` and `0.738` were measuring how long each march
happened to continue. Two runs were re-marched to about twice the rise time
past their own peaks:

| release | peak | marched to | lag past peak | after-peak fraction |
|---|---|---|---|---|
| `3.0` | `1.4` | `3.0` | `1.6` | `0.506` |
| `4.5` | `2.4` | `4.1` | `1.7` | `0.458` |

**Spread `0.047`, against `0.273` before.** Putting the runs on a common
clock collapses most of the disagreement. So the first sweep's spread really
was an artifact of the stopping rule, and that much is now settled.

What is **not** settled is the relaxation itself. Agreement here bounds the
relaxed amplitude at the lag actually reached; it does not establish that the
amplitude has come to rest there. The runs ended at roughly twice the rise
time, and nothing here says where it goes after that. No single relaxed
value is claimed.

The compression gain and rise time are unchanged by all of this: `1.1972`
and `1.1875`, spread `0.0082`, rise `0.800` at both.

## The glob mixed two resolutions, and said so

The first attempt at the settled sweep used `--pattern "nsc-rn-inner-packet-c*"`,
which matched `c2`, `c4p5` **and `c3-fine`**, the `1740`-point resolution
check. It produced a record claiming to be a release-radius sweep at one
resolution whose gains were `1.1543`, `1.2285`, `1.1875` — the middle one a
different grid entirely.

The driver checked `mu`, `eps`, packet width and `r_inner` for equality and
did not check `points`. It now takes an explicit list of run names and
refuses any set that differs in resolution:

    runs differ in points: [869, 869, 1740, 869]. This is not a
    release-radius sweep; refusing to reduce it as one.

Record `nsc-rn-packet-sweep-v3.json` is retained and is self-evidently
invalid: its own `held_fixed.points` reads `[869, 1740]`.

The prose in the sweep record was also hardcoded from the first run set, so
the settled record inherited a statement about the relaxation that was false
for it. That prose is now computed from the data, and a sweep that finds
nothing failing says so rather than repeating an old claim.

## Where the line stands

Held across a deliberate change of the release radius, spanning more than a
factor of two, at a fixed resolution:

- the density peak at the inner boundary **compresses by about a fifth**,
  gain `1.15` to `1.20`;
- it reaches its peak **about `0.8` time units** after the shell arrives,
  however far the shell fell.

Held under a doubling of the grid:

- the same gain, to `2.55%`.

Not held:

- the width compression, which grows with release radius;
- the relaxed amplitude, which is bounded at the lag reached, not settled.

## Scope

Construction checks on finite marched collars. This does not measure the
actual universe. Three release radii at one resolution, two resolutions at
one radius, one mass ratio, one source strength, one packet width. The peak
is always at the inner boundary and nothing here observes what the compressed
mass becomes below it.

## Provenance

- marches: `lab/scripts/run_nsc_rn_inner_packet.py`
- resolution reduction: `lab/scripts/derive_nsc_rn_resolution_check.py`
- sweep reduction: `lab/scripts/derive_nsc_rn_packet_sweep.py`, now taking
  an explicit run list
- records: `nsc-rn-resolution-check-v1.json`,
  `nsc-rn-packet-sweep-v4.json`, `nsc-rn-packet-sweep-settle-v3.json`,
  and the invalid `nsc-rn-packet-sweep-v3.json`
- tests: `lab/tests/test_nsc_rn_resolution_and_settle.py`

## What would settle it

A third resolution, which would let a Richardson extrapolation say whether
the gain is converging or merely stable so far. And marching past roughly
four rise times after the peak, which is where the question of whether the
amplitude rests at all becomes decidable.