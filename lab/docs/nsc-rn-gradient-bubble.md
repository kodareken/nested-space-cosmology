# Two-force gradient balance on the RN radial parent

This note owns one finite question asked on recorded states: on the areal
radius of the nonperiodic horizon-regular PG parent, do the enclosing and
the enclosed form a balanced pair, or does one force run away and leave
the other as decoration?

It is a construction check. It does not measure the actual universe, and a
bounded balance here is not evidence about cosmology.

## Domain

Areal radius `r` on the recorded radial grid, from the fixed inner mass at
`grid[0]` to the absorbing outer boundary at `32 r_m`. Observer is the
areal chart of the recorded coupled arm. No coordinate change, no time
slice, no second geometry. The state is a station payload already written
by `nsc-rn-neutral-family-v1`; this step takes no time step and re-runs
nothing.

The two forces are read from the same payload and are not fitted or
reassembled. `radial_pressure` is the inherited radial contribution, the
compressive side. `angular_pressure` is the generated angular
contribution, the resisting side. In `nsc_rn_source` they are

```text
physical_pr = (force_N - angular_force) / (4 pi volume)
physical_pa = -force_r / (2 lapse radial_density r 4 pi)
```

so the radial entry is the matter density with the circulating term
removed, and the angular entry is the transverse pressure.

## Balance and gradient

`B(r) = |F_in| / (|F_in| + |F_out|)` is scale free. `B = 1/2` is two equal
forces and `B -> 1` means the inside acts alone. The gradient is the
logarithmic slope `d ln|F| / d ln r`, which is `-1` for inverse-radius
decay and steeper than that for a force growing without bound inward.

The integrated forces are `4 pi int |F| r^2 dr` on the recorded grid.

## Provenance, and one mistake that was caught

The station payload carries two curvature entries with different origins.
`R4` is reconstructed from the evolved PG jets and belongs to the same
state as the pressures. `Ricci2` is the analytic source-free RN reference
`4 magnetic_r2^2 / r^8`, matching that closed form to about `1.7e-2`
relative on this grid.

An earlier version of the curvature test divided `Ricci2` by the force
density. That manufactured a ratio of about `9.4e+05` and a false
`curvature_dominates` dead end. The reference curvature is not a force and
is never compared with one. The module now rejects any curvature key other
than `R4` by name, and a test pins that rejection. The same-provenance
ratio, `max|R4|` over the peak force density, is about `12.5`.

## Declared limits

The limits are conventions of this check and are recorded in the record
body, not fitted to the outcome.

| Limit | Value | Meaning |
|---|---|---|
| `throat_slope_limit` | `8.0` | steeper inner slope is a collapse signature |
| `single_force_balance` | `0.90` | balance never below this means the outer force is absent |
| `exterior_decoupled` | `0.999` | reported as `r_decouple`, never a dead end |
| `curvature_force_ratio` | `1.0e+03` | same-provenance curvature over force density |
| `throat_fraction` | `0.10` | fraction of the collar the throat test judges |

The throat limit sits close to the data: `mu=1.040, eps=0.001` passes at
`-7.37` and `mu=1.040, eps=0.010` fails at `-8.44`. The slopes are in the
record so the reader can apply a different convention without re-running.

## Measured outcome

Six cases, 177 stations, reduced from what is actually on disk.

| Case | Reached | Balance at throat | Throat slope | `F_out/F_in` | `r_decouple` | Verdict |
|---|---|---|---|---|---|---|
| `mu-1.002-eps-0.001` | 8.00 | 0.247 | -121.60 | 0.819 | 4.66 | `throat_runaway` |
| `mu-1.002-eps-0.010` | 8.00 | 0.292 | -152.64 | 0.820 | 4.64 | `throat_runaway` |
| `mu-1.010-eps-0.001` | 15.50 | 0.084 | -13.52 | 0.226 | 3.72 | `throat_runaway` |
| `mu-1.010-eps-0.010` | 16.00 | 0.557 | -18.54 | 0.224 | 3.31 | `throat_runaway` |
| `mu-1.040-eps-0.001` | 18.50 | 0.915 | -7.37 | 0.156 | 4.04 | `two_force_balance` |
| `mu-1.040-eps-0.010` | 19.50 | 0.607 | -8.44 | 0.197 | 3.11 | `throat_runaway` |

Two relations fall out of the table.

The outside always acts. The integrated outer force is `15.6%` to `82.0%`
of the inner force, and no case falls below a tenth. The non-local content
is present in every recorded state, which is the part of the picture that
a local-only computation would miss.

The two forces trade rather than settle. Where the exterior compression
weakens from `mu=1.040` toward `mu=1.002`, the required throat gradient
steepens from about `-7` to about `-153`, and the outer force grows from
`16%` to `82%` of the total. Balance is a trade along the family, not a
fixed point. Five of six cases then trip `throat_runaway`, and only
`mu=1.040, eps=0.001` reaches a balance verdict. A sustained bubble is not
established across this family; it appears at most in a narrow window near
`mu = 1.04`.

## Coverage, and why it is partial

The status page previously recorded this family as still running. It was
not. On disk, no case has a `summary.json`. Four cases stopped on the CPU
projection at `t = 15.98`, `16.20`, `18.83` and `19.99` with the recorded
blocker `march_pair stopped before the chunk exceeded the CPU
projection`. The two `mu=1.002` cases stopped at `t = 8.00` with no
failure file at all, which is an unrecorded stop and not yet explained.
The record reports `family_complete: false` and refuses to describe a
family that does not exist.

## Named blocker: the family cannot currently be re-marched

`nsc_rn_pg._fornberg` builds finite-difference weights by the standard
Fornberg recursion and then verifies them against the polynomial moments
with a `1e-8` relative tolerance. On this runtime that verification fails
at collection time for every width of 9 or more. The worst relative moment
error is `2.7e-5` at width 11, third derivative. This is round-off
accumulating in the recursion, not a wrong formula, but `2.7e-5` is too
large to accept for a third derivative, so the guard is right to raise.

This is pre-existing and unrelated to this note: `tests/test_nsc_rn_pg.py`
fails collection with or without the new module present. It is also the
reason the CPU-blocked cases cannot simply be finished. Restoring the
family needs a stabilized weight construction, which belongs to the
`nsc_rn_pg` owner and is not attempted here.

## Owners

The reduction is
[src/recursive_horizons/nsc_rn_gradient_bubble.py](../src/recursive_horizons/nsc_rn_gradient_bubble.py),
driven by
[scripts/derive_nsc_rn_gradient_bubble.py](../scripts/derive_nsc_rn_gradient_bubble.py),
with tests in
[tests/test_nsc_rn_gradient_bubble.py](../tests/test_nsc_rn_gradient_bubble.py).
The record is
[results/development/nsc-rn-gradient-bubble-v1.json](../results/development/nsc-rn-gradient-bubble-v1.json),
creation only, `89438` bytes, sha256 prefix `5e7220601b8dd0b2`, produced
against `12d56e23`. The input fields come from
[nsc_rn_source.py](../src/recursive_horizons/nsc_rn_source.py) and the
curvature provenance from
[nsc_rn_reference.py](../src/recursive_horizons/nsc_rn_reference.py).

The traction concept in [nsc-discovery-parent-traction.md](nsc-discovery-parent-traction.md)
is a lead about the same question in a different domain, the discovery
collar `|x| <= 1`. That boundary gradient charge is not this areal-radius
force balance, and the two are not identified here.

## What this does not say

It does not say the universe is a gradient bubble. It measures a trade
between two recorded forces on one finite grid, and it finds that the
balance is narrow rather than general. It does not finish the family,
explain the two silent stops at `t = 8`, repair the Fornberg recursion, or
certify a continuum limit.