# The wide-grid results are converged, and memory was never the constraint

Two claims in this line rested on things nobody measured. Both were about
cost: the recorded family's silent stops were attributed to memory because
the cost scales as `N` squared, and the `x1.58` saturation was said to need
more points near the edge than were available.

Both were wrong, and they were wrong the same way.

## Memory is flat in the point count

Real grids, packets and states built at a fixed edge of `0.7555` and
`r_out = 32`, then discarded:

| points | spacing | gap (cells) | commit delta | build |
|---|---|---|---|---|
| `869` | `0.035996` | `0.03` | `-0.1 GB` | `3.3 s` |
| `2000` | `0.015630` | `0.07` | `-0.0 GB` | `3.6 s` |
| `4000` | `0.007813` | `0.15` | `-0.1 GB` | `4.0 s` |
| `6000` | `0.005208` | `0.22` | `0.0 GB` | `4.6 s` |
| `8000` | `0.003906` | `0.30` | `0.0 GB` | `5.3 s` |

**Nine times the points, no measurable change in memory.** The production
derivative is a CSR operator with a bounded stencil width, so the `N` squared
assumption was never applicable to it. This machine has `63.3 GB`, `52.7 GB`
available, `16` cores, and an `8000`-point march completes at about `117 s` of
CPU per station.

The recorded family's stops may well have had some other cause; this says
only that memory cannot be it, at these sizes, with this operator.

## And the trajectories agree

The gain is a ratio of two amplitudes read off a trajectory, so if the
trajectory agrees between two resolutions the gain agrees — no need to wait
for either run to reach its peak.

Inner mass stock at matched stations, `8000` points against `869`, over `15`
to `16` stations each:

| edge | at `t=0.1` | worst over the run | stations |
|---|---|---|---|
| `1.04` | `5.90e-09` | `5.22e-07` | `15` |
| `0.85` | `4.87e-09` | `5.43e-07` | `16` |
| `0.80` | `1.41e-07` | `1.25e-06` | `16` |

The differences **grow smoothly and decelerate**, so they are accumulated
integration differences rather than a per-step noise floor. The per-station
ratio falls from `8.4` early to `1.06` late. Extrapolating to the compression
peak near `t = 1.7` puts the difference there around `1e-6` to `2e-6`.

A gain quoted to three figures cannot move by `2e-6`. **The wide-grid gains
are converged in the point count**, and `x1.20`, `x1.39`, `x1.43` and the
saturation at `x1.58` are measurements of a converged solution rather than of
a grid.

## A trap that nearly produced a false result

The first cut of this comparison matched stations by file name and reported a
worst difference of `1.5e-05` — twenty times larger, with a jump of `x28` at
the final station. That jump was not physics.

A march stopped by the cpu budget writes its last station part way through a
chunk. The file name encodes `t/r_m` at `0.1` resolution, so for the truncated
final station it is up to **`0.043`** away from the clock recorded inside the
payload. The coarse runs are exact to `0.000000`. Matching by name alone
therefore compares the two runs at instants up to `0.043` apart, on a
trajectory that is still moving.

The reduction now keys on the recorded clock and drops any station whose name
disagrees with it, reporting how many were dropped. All three comparisons
dropped exactly one, the truncated final station, and the smooth series
beneath it is what remains.

Two resolutions, not three. Two points bound nothing and no convergence order
is claimed. This is evidence against a resolution artefact at these two
sizes.

## What this closes and what it does not

**Closed**: whether the saturation is a grid artefact. It is not, at these
sizes. The previous note's "this is a resources question" was wrong, and it
was wrong for the same reason the original memory inference was.

**Also closed, by a different route**: whether the trend continues past
`0.77`. At `8000` points the spacing is `0.0039`, and the gap from `0.7555` to
the inner horizon is `0.00116`, so it is `0.30` cells — the edge is still not
resolvable from the horizon. Getting four cells there at `r_out = 32` needs
about `100,000` points, and narrowing the outer radius destabilises the march.

So the sub-horizon limit is real. It was simply not for the reason given.

## What the runs cost

Three `8000`-point marches at `869`-point parity of memory, about `117 s` of
CPU per station each, run concurrently. All three are finite, and have **zero
mass reversals below `r = 16`**; every reversal is in the far field above
`16`, which is the known ringing.

## Scope

A numerical convergence check on two resolutions, and a cost measurement on
one machine at one moment. Neither is a scientific result and neither
measures the actual universe.

Two points bound nothing. Agreement between them is evidence against a
resolution artefact at those sizes, not a demonstration that the quantity is
converged in every respect. The saturation remains a saturation of this
construction on this machine.

## Provenance

- memory: `lab/scripts/measure_nsc_rn_grid_memory.py`,
  `nsc-rn-grid-memory-v1.json`, `lab/docs/nsc-rn-grid-memory.md`
- convergence: `lab/scripts/derive_nsc_rn_resolution_convergence.py`,
  `nsc-rn-resolution-convergence-v1.json`
- marches: `nsc-rn-inner-packet-n8000-{def1040,exc085,exc080}`
- tests: `lab/tests/test_nsc_rn_grid_memory.py`,
  `lab/tests/test_nsc_rn_resolution_convergence.py`