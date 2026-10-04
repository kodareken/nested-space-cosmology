# Putting the packet inside the collar: the accumulation, resolved

The recorded neutral family releases its source packet at `12 r_m`, from
outside the collar it is about to cross. The mass-shell reduction showed
what that costs: the shell falls as a coherent body and reaches the grid's
inner edge at `t = 7.0` in every case, and by the end of every run between
**61% and 99% of the accreted mass lies below that edge**. The accumulation,
which is the part of the picture that matters, happens outside the domain.

This note is the first look at that accumulation.

## The construction

Same packet, same width of `2 r_m`, same spectrum, same chart, same `mu`
and `eps`. Only the release radius changes, from `12 r_m` to `3 r_m`. The
shell then forms inside the collar and falls a much shorter distance, so
its arrival and the profile it leaves behind are resolved rather than
truncated.

It is a different construction from the recorded family and its numbers are
not comparable with that family's except through quantities that depend only
on the construction.

`mu = 1.04`, `eps = 0.001`, `869` points, `r` from `1.04` to `32`, marched to
`t/r_m = 2.0` in `21` stations.

## What the shell does

| | measured |
|---|---|
| forms at | `r = 2.8234` |
| falls at | `-1.8446` per unit `t/r_m` |
| fit quality | `r^2 = 0.99684` |
| reaches the inner edge | `t = 0.6` |
| half-max width, start | `4.387`, or `123` cells |
| half-max width, arrival | `3.567`, or `100` cells |

The shell **compresses as it falls**, narrowing by `19%` over the fall. The
shrink is `23` grid cells against a `15`-cell median window, so it is the
shell narrowing and not the estimator drifting.

## What the accumulation does

This is the part the recorded family could not reach. After the shell reaches
the inner boundary at `t = 0.6`, the density peak there keeps growing:

| | |
|---|---|
| amplitude at arrival | `2.468e-04` |
| peak amplitude | `2.955e-04` at `t = 1.4` |
| gain | **x`1.197`** |
| amplitude at the final station | `0.921` of the peak |

So the mass does not simply pile up against the boundary and stop. It
**compresses by about a fifth** over the `0.8` time units after arrival, then
slowly relaxes, sitting at `92%` of its peak at `t = 2.0`.

That is the first resolved statement about the accumulation geometry that
this programme has produced. In the recorded family the same mass leaves the
domain and is never seen again.

## Resolution limits, stated because they bind

**The peak radius is quantised to the grid.** The shell advances `5.000`,
`4.000`, `6.000`, `5.000`, `6.000` cells per station, exactly. The apparent
alternation in local speed through the falling phase is therefore the grid
and not dynamics. The mean rate over the phase is the reliable quantity; the
per-interval speeds are not evidence of anything, and no test here uses them.

**The fitted rate is a mean, not a law.** Extrapolating it past arrival gives
`1.717` where the shell was last seen at `1.896`. The final segment is
steeper than the mean. With only six falling stations that is not enough to
characterise, and no claim is made about it.

**The width is measured on a smoothed profile.** It is a half-max width on a
fifteen-cell median, so it is a width of the envelope, not of a sharp edge.

**One case, one parameter pair.** The release radius was varied once. The
reading that the compression is real rests on one run, and would be stronger
with the packet released at two or three radii.

## What was corrected on the way

**The coherence test was wrong and it truncated the falling phase.** It
compared the shell's peak against the density at the inner edge and required
a ratio of three. A shell approaching the boundary has a genuinely sharper
peak than the same shell in flight, so the test rejected exactly the
stations where the shell was still moving. The recorded family's falling
phase was cut at `t = 5.0` when the shell was in fact coherent through
`t = 6.5`. Coherence is now decided by **position** against the inner edge
plus amplitude continuity, which is the right discriminator: once a shell is
pinned, a position test catches it and an amplitude test does not.

With that fixed the recorded family's fitted infall rate moves from `-1.444`
to `-1.502` per unit, over `84` falling stations instead of `66`, with
`r^2` from `0.9982` to `0.9985` and arrival still `t = 7.0` in all six cases.
The spread across the family stays at `0.7%`, so the conclusion that the
infall rate is nearly parameter-independent is unchanged.

**The first reduction of this construction scaled its clock by the mass
ratio.** `t/r_m` is a clock and `r_m` is the magnetic radius, which is not
`mu`. The station names carry the clock at the recorded `0.1` resolution, so
the reduction now reads its clock from the name. This is why arrival is
`t = 0.6` and not `0.577`.

Records `v1` and `v2` of both reductions are retained and sealed and contain
those errors.

## Scope

A construction check on one finite marched collar. It does not measure the
actual universe. A shell that compresses by a fifth before relaxing against a
boundary in one numerical construction is not evidence that any real object
does so.

The packet is still released *above* the inner edge. Nothing here observes
what the mass does after it crosses, which remains the same limit the
recorded family had.

## Provenance

- march: `lab/scripts/run_nsc_rn_inner_packet.py`
- reduction: `lab/scripts/derive_nsc_rn_inner_packet.py`
- march record: `lab/results/development/nsc-rn-inner-packet-v1/`, 21 stations
- reduction record: `lab/results/development/nsc-rn-inner-packet-shell-v2.json`
- shell estimator: `lab/src/recursive_horizons/nsc_rn_mass_shell.py`
- tests: `lab/tests/test_nsc_rn_inner_packet.py` and
  `lab/tests/test_nsc_rn_mass_shell.py`

## What would settle it

Two things. First, vary the release radius. If the compression factor and the
arrival behaviour hold across `2 r_m`, `3 r_m` and `4.5 r_m`, they are
properties of the construction rather than of one run.

Second, and larger: start the march *below* `1.04` so the accumulation region
itself is inside the domain. That needs an inner boundary condition at a
radius between `r_minus = 0.754` and `r_in = 1.04`, which is a narrow window
and a genuinely harder construction. It is the only way to see what the
compressed mass becomes rather than where it stops.