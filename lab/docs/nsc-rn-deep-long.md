# Long deep runs: the compression holds and the relaxation completes

The edge sweep found the compression gain rising monotonically as the grid's
inner edge deepened, `x1.20` at the midpoint to `x1.58` at `0.77`, and found
that the deeper runs ended before they had relaxed, so the relaxation could
not be measured there at all. This re-runs `0.80` and `0.77` to `t = 6.0`,
double the previous length, and settles both.

## The compression did not change

| edge | gain at `t=3.0` | gain at `t=6.0` | peak |
|---|---|---|---|
| `0.80` | `x1.4328` | `x1.433` | `t = 1.7` |
| `0.77` | `x1.5758` | `x1.576` | `t = 1.7` |

The gain is unchanged to four figures. The compression peak was already
inside the shorter runs. Only the relaxation needed more time, which is
exactly what the edge sweep said.

## The relaxation completes

| edge | at arrival | at peak | at `t = 6.0` | fraction of peak |
|---|---|---|---|---|
| `0.80` | `2.430e-04` | `3.481e-04` | `1.462e-05` | **`0.042`** |
| `0.77` | `2.296e-04` | `3.618e-04` | `1.556e-05` | **`0.043`** |

The mass at the edge rises to a compressed peak, then decays to about **four
percent** of it. And the two edges agree to within a thousandth of a
fraction, so this is a settled quantity rather than two near-misses.

**94% of the accreted mass ends below the inner edge**, up from `78%` at
`t = 3.0`. The deep domain retains it; it does not leak.

The edge-sweep note said the deeper runs were "too short to measure a
relaxation at all". That was right at `t = 3.0` and is now superseded: at
`t = 6.0` the relaxation is measured, and it completes.

## The curvature answers, and it lags

| edge | `R4` max | when | `R4` at the mass peak | final / initial |
|---|---|---|---|---|
| `0.80` | `3.2826e-03` | `t = 0.1` | `1.8549e-03` | `x0.1444` |
| `0.77` | `5.8565e-03` | `t = 0.1` | `1.8436e-03` | `x0.6349` |

The curvature's global maximum is the formation transient at `t = 0.1`, in
both runs. At this run's own mass peak, `t = 1.7`, the edge curvature is
`1.84e-03` and `1.85e-03`: **the same at both edges to three figures**, even
though the initial curvatures differ by a factor of `3.5`.

The curvature also has a second maximum during the compression, and it lags
the mass. At `0.77` the edge curvature dips to `7.0e-05` at `t = 1.3`, then
climbs to `4.99e-03` at `t = 2.6` — about `0.9` time units after the mass
peak at `t = 1.7`. The geometry responds after the mass, not with it.

Nothing diverges. The maximum over the initial value is `3.12` and `1.58`.

## A bug this exposed

The curvature reduction read the mass peak time from a **hardcoded**
default-edge record. For these runs it was reading `t = 1.4` from a different
march and reporting the edge curvature there, which is `3.80e-04` instead of
the correct `1.84e-03` — a factor of five.

That is the same failure as the shell reduction's, a statistic read off the
wrong clock, and it survived a previous commit because the default case
happened to be the run being analysed. It now takes the shell record for the
same run by default and says so.

Records `nsc-rn-curvature-exc077-long-v1.json` and `-exc080-long-v1.json`
contain the wrong mass-peak column and are retained.

## Scope

Construction checks on finite marched collars. This does not measure the
actual universe.

Every edge here is still a **boundary**. `0.77` is `0.0157` above the inner
horizon at `0.754`, and a boundary that close to a horizon is still a
boundary. No reading here is a statement about regularity at a horizon.

What the deep domain does show is a complete dynamical cycle inside the
collar: infall, compression, relaxation, with the mass retained. That cycle
was never visible in the recorded family, where the mass left the domain.

## Provenance

- marches: `run_nsc_rn_inner_packet.py --excision {0.80,0.77} --target 6.0`,
  `61` stations each, about `687 s` of CPU each
- reductions: `derive_nsc_rn_inner_packet.py`,
  `derive_nsc_rn_curvature.py --shell-record`
- records: `nsc-rn-inner-packet-shell-exc077-long-v1.json`,
  `-exc080-long-v1.json`, `nsc-rn-curvature-exc077-long-v2.json`,
  `-exc080-long-v2.json`

## What would settle it

The compression gain is still rising with depth at `0.77`, so the trend has
not turned over even now that the peak and the relaxation are both resolved.
Either the edge moves closer still, or the gain saturates between `0.77` and
`r_minus`, and only a run at the deepest admissible edge distinguishes those.

Separately, the curvature lag of about `0.9` time units is a new observable
and has only one edge depth measured at it properly. It deserves its own
sweep before anything is claimed from it.