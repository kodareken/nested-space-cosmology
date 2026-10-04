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

## The 2005 criterion, and the hard negative result

Douglas's framing pointed at "the 2005 study of black hole". That is
Hayward, [gr-qc/0506126](https://arxiv.org/abs/gr-qc/0506126), "Formation
and evaporation of non-singular black holes", published 28 June 2005. It
was already in the arXiv index as an archived citation; it is now read.

The static region there is "supported by finite density and pressures ...
vanishing rapidly at large radius and behaving as a cosmological constant
at small radius", and the black hole is "a compact space-time region of
trapped surfaces, with inner and outer boundaries". That is the two-force
picture directly: a finite interior stress, an exterior that dies away,
and a bounded region between them.

It also supplies the regularity criterion this check was missing. A force
with logarithmic slope $g$ behaves as $r^g$, so it stays finite as
$r\to0$ only when $g\geq-1$. Anything steeper diverges inward and the
interior is not regular in the sense of that paper. The limit is sourced,
not chosen, which matters because the previously used limit of $8$ was
mine and the data straddled it.

On that criterion **no case in the family has a regular interior**:

| Case | Throat slope | Regular at $g\geq-1$ |
|---|---|---|
| `mu-1.002-eps-0.001` | -121.60 | no |
| `mu-1.002-eps-0.010` | -152.64 | no |
| `mu-1.010-eps-0.001` | -13.52 | no |
| `mu-1.010-eps-0.010` | -18.54 | no |
| `mu-1.040-eps-0.001` | -7.37 | no |
| `mu-1.040-eps-0.010` | -8.44 | no |

The best case is still more than seven times steeper than the criterion
allows. This is the decisive result of the step and it is negative: the
recorded RN radial states are not the regular interior that the 2005
construction describes. `interior_not_regular` is therefore a hard
dead-end flag, and four of the six cases carry it at their final station.

This also explains the numerical blocker. A throat slope of $-153$ is
exactly the regime in which the Fornberg recursion loses its polynomial
moments, so the failing moment check and the failing regularity are the
same fact seen twice.

## Declared limits

The limits are conventions of this check and are recorded in the record
body, not fitted to the outcome. None of them decides the verdict.

| Limit | Value | Meaning |
|---|---|---|
| `throat_slope_limit` | `8.0` | steeper inner slope is reported, never a verdict |
| `single_force_balance` | `0.90` | balance never below this means the outer force is absent |
| `exterior_decoupled` | `0.999` | reported as `r_decouple`, never a dead end |
| `regular_slope_limit` | `-1.0` | sourced from gr-qc/0506126; the hard regularity test |
| `regular_slope_tolerance` | `1e-9` | so an exactly regular power law survives round-off |
| `curvature_force_ratio` | `1.0e+03` | same-provenance curvature over force density |
| `throat_fraction` | `0.10` | fraction of the collar the throat test judges |

The throat limit sits close to the data: `mu=1.040, eps=0.001` reports at
`-7.37` and `mu=1.040, eps=0.010` at `-8.44`. This is exactly why it was
demoted. The slopes stay in the record so a reader can apply a different
convention without re-running.

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
a local-only computation would miss. No hard dead-end flag fires anywhere
in the family.

The balance wall exists but does not settle. A wall is present at
`49.7%` of all stations and appears in all six cases at some station, so
the two forces genuinely become equal somewhere in the collar. But the
wall oscillates rather than holding a fixed radius, and its persistence
falls as the exterior compression strengthens:

| Case | Wall present | $r_{\mathrm{wall}}$ range |
|---|---|---|
| `mu-1.002-eps-0.001` | 65% | 1.073 to 2.054 |
| `mu-1.002-eps-0.010` | 65% | 1.073 to 2.054 |
| `mu-1.010-eps-0.001` | 69% | 1.045 to 2.055 |
| `mu-1.010-eps-0.010` | 61% | 1.038 to 2.055 |
| `mu-1.040-eps-0.001` | 32% | 1.290 to 2.039 |
| `mu-1.040-eps-0.010` | 30% | 1.290 to 2.039 |

A stronger exterior dissolves the wall. Persistence drops from about
`65%` to about `31%`, and the inner excursion of the wall is pushed away
from the throat, from `1.04` to `1.29`. Only the two `mu=1.002` cases
still hold a wall at their final station.

So the honest reading is not that the bubble exists or that it does not.
The inside and outside do reach equality, for roughly half of the
recorded time. What is missing is a settled wall: the balance is
transient, and stronger exterior compression shortens it. A sustained
gradient bubble is not established by this family.

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
[results/development/nsc-rn-gradient-bubble-v3.json](../results/development/nsc-rn-gradient-bubble-v3.json),
creation only. It names
[nsc-rn-gradient-bubble-v2.json](../results/development/nsc-rn-gradient-bubble-v2.json)
in a `supersedes` block with the reason, and lists
[nsc-rn-gradient-bubble-v1.json](../results/development/nsc-rn-gradient-bubble-v1.json)
under `earlier_records`. All three are preserved. The
input fields come from
[nsc_rn_source.py](../src/recursive_horizons/nsc_rn_source.py) and the
curvature provenance from
[nsc_rn_reference.py](../src/recursive_horizons/nsc_rn_reference.py).

The traction concept in [nsc-discovery-parent-traction.md](nsc-discovery-parent-traction.md)
is a lead about the same question in a different domain, the discovery
collar `|x| <= 1`. That boundary gradient charge is not this areal-radius
force balance, and the two are not identified here.

## What this does not say

It does not say the universe is a gradient bubble. It measures a trade
between two recorded forces on one finite grid, and it finds the balance
real but transient, and the interior non-regular on the 2005 criterion.
It does not certify a continuum limit, and it does not separate the wall
oscillation from any step-control artifact. It does not finish the family,
explain the two silent stops at `t = 8`, repair the Fornberg recursion, or
say what happens at `mu` above `1.04`, which was never run. And the
non-regular result is a statement about these recorded states, not a
refutation of the 2005 paper, which describes a different construction.