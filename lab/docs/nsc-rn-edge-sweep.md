# The edge sweep: the compression has still not converged

Moving the grid's inner edge from the horizon midpoint to `0.85` changed two
numbers in the same direction and implied the earlier ones were truncation
artifacts. This asks whether the trend continues as the edge approaches the
inner horizon, or turns over because the chart stops supporting the region.

It continues. It has not stopped.

| `r_inner` | gap to `r_minus` | arrival | **gain** | peak | `R4` edge at `t=0` | max / initial |
|---|---|---|---|---|---|---|
| `1.04` | `0.2857` | `0.6` | `x1.1972` | `t=1.4` | `3.103e-04` | `2.68` |
| `0.85` | `0.0957` | `0.5` | `x1.3905` | `t=1.6` | `5.644e-04` | `3.39` |
| `0.80` | `0.0457` | `0.5` | `x1.4328` | `t=1.7` | `1.052e-03` | `3.12` |
| `0.77` | `0.0157` | `0.4` | **`x1.5758`** | `t=1.7` | `3.697e-03` | `1.58` |

**The compression gain is monotone increasing across the whole sweep.** It
is still growing at the deepest edge tried, which sits `0.0157` above the
inner horizon at `0.754`.

So the `x1.20` reported at the midpoint is a **lower bound by a factor of
about `1.3`**, and the effect has not converged with respect to domain depth.
Every headline number in this line that was measured at the midpoint
excision is still truncated.

## The relaxation is not the same story

The sweep's verdict reads `effect_does_not_continue_monotonically`, and that
is about the relaxation only. Its final-to-initial edge curvature ratio goes
`0.0145`, `0.6922`, `1.4703`, `1.2457`.

Ratios above one mean the edge curvature at the final station is *larger*
than at the first. At `0.80` and `0.77` the system has not relaxed at all by
`t = 3.0`; it is still rising.

That is not a failure of the trend, it is the trend: a deeper domain holds
more of the accumulation, and more of it is still arriving at the end of the
run. **These runs are too short to measure a relaxation at all.** The
`x0.0145` reported at the midpoint was not a relaxation being observed; it
was a domain that had emptied.

The previous note's claim that relaxation weakens as the edge deepens is
therefore withdrawn. It is not a trend, it is a threshold crossing: below
some depth the run outlasts the dynamics.

## The deeper edge starts in more curved geometry

The initial edge curvature grows by a factor of about `12` across the sweep,
`3.10e-04` to `3.70e-03`, before any matter arrives. The deeper edge is
simply sitting in a more strongly curved part of the reference chart, which
is what approaching a horizon does.

So the gain growing across the sweep is not simply "more mass, more
compression". Part of the extra compression is measured against a stiffer
starting geometry. This reduction does not separate the two, and no claim is
made that it does. Separating them would need a matched comparison at fixed
initial curvature, which this is not.

Curvature stays bounded everywhere. The maximum over initial is `2.68`,
`3.39`, `3.12`, `1.58`; nothing diverges at any edge tried.

## What this does not establish

Every edge in this sweep is still a **boundary**. The deepest, `0.77`, is
`0.0157` above the inner horizon at `0.754`, and a boundary `0.016` above a
horizon is still a boundary. Regularity at a boundary remains a much weaker
statement than regularity at a horizon, and no reading here is a statement
about regularity at one.

What it does establish is that the effect is **not converged**. The honest
next step is not another edge but a longer run at a fixed deep edge, to see
whether the compression ratio keeps climbing or turns over once the domain
finally contains the whole accumulation.

## Scope

Construction checks on finite marched collars. This does not measure the
actual universe. Four edge positions, one release radius, one mass ratio, one
source strength, one grid resolution. A shell compressing by `x1.58` against
a numerical boundary `0.016` from a horizon in one construction is not
evidence that any real object does so.

## Provenance

- march: `lab/scripts/run_nsc_rn_inner_packet.py --excision`
- reduction: `lab/scripts/derive_nsc_rn_edge_sweep.py`, runs named explicitly
- runs: `nsc-rn-inner-packet-curv` (default edge), `-exc085`, `-exc080`,
  `-exc077`
- record: `lab/results/development/nsc-rn-edge-sweep-v1.json`

## What would settle it

Two things, and they are different sizes.

A longer run at `0.77`, to `t` well past the accumulation, to see whether
the gain turns over or keeps climbing.

And a matched comparison that separates the two things the sweep currently
mixes: compression against a *softer* starting geometry at the same depth,
and compression against a *stiffer* one. Without that the gain's growth is
real and its cause is not identified.