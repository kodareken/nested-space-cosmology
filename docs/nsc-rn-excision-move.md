# Moving the excision: the compression was truncated, and so was the relaxation

The recorded family and every inner-packet run so far started their grid at
the **horizon midpoint**, `excision = 0.5 * (r_minus + r_plus) = 1.04`. Mass
piled up against that edge and was then cut off. The compression and the
relaxation measured there were both bounded by where the domain ended.

`build_grid` now takes the excision as an optional argument. The default is
the midpoint exactly, so a recorded run that does not pass it is
bit-identical to one that predates the change; the curated suite is
unchanged at `368 passed, 1 skipped`.

This note is the first march with the edge moved inboard, to `0.85`, so the
accumulation is inside the domain instead of against its wall.

## The chart admits it

Worth checking before building anything, because the obvious objection is
that the PG chart goes imaginary below the outer horizon. It does not.
`beta^2 = 2M/r - Q^2/r^2` vanishes at `r = Q^2/(2M) = 0.4808`, comfortably
below `r_minus = 0.754`, and `N^2 = f + beta^2 = 1` identically across the
whole interval. The shift and the characteristic speeds stay smooth.

The grid already started *below* `r_plus = 1.326`, so the new region is an
extension of territory the chart already covered, not a new chart.

## The compression is larger

Same packet, same release radius of `3.0`, same `mu`, `eps`, width and
`869` points. Only the inner edge moved, by `0.19 r_m`.

| edge | gain | peak | after peak | gain below edge |
|---|---|---|---|---|
| `1.04` (midpoint) | `x1.1972` | `t = 1.4` | `0.921` | `0.583` |
| **`0.85`** | **`x1.390`** | `t = 1.6` | `0.597` | **`0.795`** |

**The compression grows from `x1.20` to `x1.39`** when the domain is deep
enough to contain it. And `79.5%` of the accreted mass ends up below the
edge, against `58.3%` before.

The earlier `x1.20` was therefore a **truncation**, not the effect. It was
measured at a boundary that cut across the middle of the accumulation, so
it was a lower bound on a construction that is still growing.

This is the same failure the recorded family had in a more severe form, one
level down: there, `61` to `99%` of the mass left the domain; here, `58%`
left it; only with the edge at `0.85` does most of it stay in.

## The relaxation was truncated too

| edge | max `|R4|` at edge | at mass peak | final | final / initial |
|---|---|---|---|---|
| `1.04` | `8.3271e-04` | `1.5108e-04` | `4.5076e-06` | `x0.0145` |
| **`0.85`** | `1.9138e-03` | `9.3333e-04` | `3.9070e-04` | **`x0.6922`** |

The curvature at the deeper edge is larger at every point, which is what a
larger accumulated mass should do. But the striking difference is at the
end: the default-edge run relaxed to `1.4%` of where it started, and this
one relaxes only to `69%`.

So the dramatic relaxation previously reported was **also** a truncation
effect. Mass that crossed the old edge took its curvature with it, and the
remaining profile went quiet. With the mass held in the domain, the
geometry has barely relaxed by `t = 3.0`.

Curvature remains bounded in both: the deeper edge peaks at `x3.39` of its
initial value during formation, and nothing diverges.

## What this does not establish

The edge is still a **boundary**, at `0.85`, and it is still outside the
inner horizon at `r_minus = 0.754`. Regularity at a boundary remains a much
weaker statement than regularity at a horizon. This is a smaller step than
it looks: `0.85` is not the horizon.

It is, however, the first step in this direction that the chart actually
admits, and it shows the previous numbers were domain-limited rather than
converged. Two release radii, two edge positions, one mass ratio, one source
strength.

The obvious next edge is as close to `r_minus = 0.754` as the closure will
admit. Whether the trend continues — compression growing and relaxation
weakening as the edge deepens — is the question, and it is answerable by
running `0.80` and `0.77`.

## Scope

Construction checks on finite marched collars. This does not measure the
actual universe. A shell compressing by a factor of `1.39` against a
numerical boundary in one construction is not evidence that any real object
does so.

## Provenance

- grid: `lab/src/recursive_horizons/nsc_rn_pg.py`, `build_grid` gains an
  optional `excision`; the default is the midpoint and is unchanged
- march: `lab/scripts/run_nsc_rn_inner_packet.py --excision 0.85`
- reductions: `derive_nsc_rn_inner_packet.py`, `derive_nsc_rn_curvature.py`,
  both now taking explicit run names
- records: `nsc-rn-inner-packet-shell-exc085-v1.json`,
  `nsc-rn-curvature-exc085-v1.json`