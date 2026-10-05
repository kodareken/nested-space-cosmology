# Memory was never the constraint

The recorded family's `mu = 1.002` cases stop silently at `t = 8`. The
evidence gathered was: the largest grids stopped earliest, all final states
are finite, and one continuation recorded `205 MB` at `869` points. From that
the notes inferred memory pressure, on the grounds that the cost scales as
`N` squared, and kept it as a standing inference.

**It was never tested. On this machine it is false.**

## What was measured

Real grids, packets and states at a fixed inner edge of `0.7555` and outer
radius `32`, built and discarded. Nothing marched.

| points | spacing | gap (cells) | commit delta | commit total | build |
|---|---|---|---|---|---|
| `869` | `0.035996` | `0.03` | `-0.1 GB` | `11.9 GB` | `3.3 s` |
| `2000` | `0.015630` | `0.07` | `-0.0 GB` | `11.8 GB` | `3.6 s` |
| `4000` | `0.007813` | `0.15` | `-0.1 GB` | `11.7 GB` | `4.0 s` |
| `6000` | `0.005208` | `0.22` | `0.0 GB` | `11.8 GB` | `4.6 s` |
| `8000` | `0.003906` | `0.30` | `0.0 GB` | `11.8 GB` | `5.3 s` |

**Nine times the points, no measurable change in memory.** The production
derivative is a CSR operator with a bounded stencil width, so its nonzero
count grows linearly. The `N` squared assumption was never applicable to it.

The machine: `63.3 GB` total, `52.7 GB` available, `16` logical cores.

## It actually runs

An `8000`-point march completes and produces stations at about `104 s` of CPU
each. `30` stations is under an hour, and three of them run concurrently with
memory flat throughout.

So the `8000`-point route that the last note called "a resources question" is
available. The conclusion that it was closed was wrong, and it was wrong for
the same reason the original inference was: nobody measured.

## But it does not fix the thing it was needed for

The reason to want more points was to resolve the gap between the inner edge
and the inner horizon. At `0.7555` that gap is `0.00116` in `r_m`, and at
`8000` points the spacing is `0.0039`, so the gap is still `0.30` cells.

Getting four cells there at `r_out = 32` needs about `100,000` points. The
march cost scales with the node count, so that is roughly `40` hours per run
for one station's worth of a second.

Narrowing the outer radius would shrink it, and that has already been tried:
it destabilises the march, with `29%` of node-steps reversing.

So the sub-horizon gap at `0.7555` is not reachable by either route. That part
of the previous note stands, for a different and now measured reason.

## What does follow

More points at the edges that *are* resolvable. At `8000` points the gap at
`0.80` is `11.7` cells, against `1.27` at `869`. That is the experiment worth
running: re-measure the gain at `1.04`, `0.85` and `0.80` at `8000` points and
see whether `x1.20`, `x1.39` and `x1.43` move. If they do not, the
wide-grid resolution was adequate and the `x1.58` saturation is trustworthy.
If they do, the saturation was a resolution artefact after all.

## Scope

A hardware and cost measurement on one machine at one moment. It says nothing
about the recorded family's silent stops, which happened on a different
occasion and may have had a different cause; it says only that memory cannot
be that cause here, at these sizes, with this operator.

It does not measure the actual universe and it is not a scientific result. It
is a measurement of what this machine can do, made so that the next scientific
step is chosen from evidence.

## Provenance

- measurement: `lab/scripts/measure_nsc_rn_grid_memory.py`, using
  `GlobalMemoryStatusEx` because `psutil` is not installed on this host
- record: `lab/results/development/nsc-rn-grid-memory-v1.json`
- tests: `lab/tests/test_nsc_rn_grid_memory.py`