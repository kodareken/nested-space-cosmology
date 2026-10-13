# The cone map on the RN radial parent

Owning note for [nsc-rn-cone-map-v1.json](../results/development/nsc-rn-cone-map-v1.json).
One finite question in one named domain: on the recorded nonperiodic
horizon-regular PG parent, what does the state's own metric say the light
cone is, and how much has the value of light changed after riding the
space from the outer boundary down to each radius?

## Domain and provenance

Areal radius `r` on the recorded radial grid of
`nsc-rn-neutral-family-v1`. Observer is the chart of the recorded coupled
arm. The state is the payload already on disk; no time step is taken and
no file is rewritten.

The coupled arm's lapse `N` and shift `beta` are not stored in the station
payload, so they are reconstructed with `nsc_rn_pg.reconstruct_geometry`
at weight 1 from the stored `coupled_phi`, `occupations` and
`coupled_inner_mass` — the same route the march used. The guard: the
reconstructed mass must reproduce the stored `coupled_mass`. Across all
177 stations the maximum residual is `0.000e+00`, so every reconstruction
is the state's own metric, not a re-fit.

## The measured cone facts

| Fact | Meaning | Reading |
|---|---|---|
| `local_speed_everywhere = 1` | the null cone is one everywhere | the invariant statement |
| `v_out = -beta + N`, `v_in = -beta - N` | the chart values | they change with radius because the recorded lapse and shift carry the room's stress |
| `ray_frequency_ratio = 1/N(r)` | light's energy after riding the room | the dilution of the carried value |
| `trap_radius` | where the recorded metric's outgoing speed crosses zero | the chart's own one-way interface |
| `redshift_at_excision = 1/N(inner)` | the value of light at the inner edge | the headline invariant |

## What the six recorded cases say

The trapped boundary of the recorded coupled metric sits **outside** the
source-free RN outer horizon in every case, and the offset is ordered:

| Case | stations | redshift@inner (last) | trap − r₊ |
|---|---:|---|---:|
| μ=1.002, ε=0.001 | 17 | 0.999946 | +0.00993 |
| μ=1.002, ε=0.010 | 17 | 0.999471 | +0.07422 |
| μ=1.010, ε=0.001 | 32 | 1.000007 | +0.00774 |
| μ=1.010, ε=0.010 | 33 | 1.000063 | +0.06869 |
| μ=1.040, ε=0.001 | 38 | 1.000004 | +0.00470 |
| μ=1.040, ε=0.010 | 40 | 1.000039 | +0.04569 |

At fixed energy fraction the trap offset **falls** as the mass ratio grows
(+0.074 at μ=1.002 to +0.046 at μ=1.040), while at fixed mass it **grows**
with the energy fraction. The one-way interface is therefore carried by
the stress content, not by the reference horizon alone.

The lapse reading at the inner edge stays within about `6e-4` of 1 in
every final station. In this flat-slice gauge the value of light barely
dilutes across the ride; the barrier is carried by the shift, which is
the recorded `claim-bindings` statement ("the local speed stays one …
the coordinate speed carries the shift") now measured on the coupled
states rather than on the fixed background.

## Horizon-reaching domain (v2)

The same reduction now runs over the inward-march campaign, whose grid stops one
cell above the inner horizon, and the driver takes `--family`/`--record`
arguments so any recorded layout is reduced without editing it. Record:
[nsc-rn-cone-map-v2.json](../results/development/nsc-rn-cone-map-v2.json)
(12 stations, 3 cases, ε=0.001, t = 5, 10, 15, 20).

The mass-reconstruction guard is `0.000e+00` there as well, so the metric is
still the state's own. The measured interface on the deeper domain:

| Case | trap − r₊ (last station) |
|---|---:|
| μ=1.002 | +0.01506 |
| μ=1.010 | +0.00781 |
| μ=1.040 | +0.00472 |

The ordering is preserved and the offsets agree with the family's ε=0.001
values to within a grid cell, so extending the domain toward the horizon did
not move the one-way interface: the offset is a property of the state's stress,
not of where the grid was cut.

## Declared limits

The family never reached a `summary.json`, and the two μ=1.002 cases
stopped at t=8.00 with no failure file; coverage is reported per case and
`family_complete` is false. No extrapolation is performed below the
recorded inner edge, which sits outside the RN inner horizon in every
case. The `no_riding_dilution` flag names the stations where the lapse
reading is at or above 1; it is a flag, not a verdict.

This is a construction check on recorded finite states. It says what the
recorded state's metric says about light. "c is a constant whose value
changes with the space it rides" is a reading of those numbers, and it
stays a reading; this note does not measure the actual universe.

## Owners

Module [`nsc_rn_cone_map.py`](../src/recursive_horizons/nsc_rn_cone_map.py),
driver [`derive_nsc_rn_cone_map.py`](../scripts/derive_nsc_rn_cone_map.py),
tests [`test_nsc_rn_cone_map.py`](../tests/test_nsc_rn_cone_map.py).

This strand also applied the guard correction documented in
[nsc-rn-gradient-bubble.md](nsc-rn-gradient-bubble.md): the Fornberg
polynomial-moment self-check is now a compensated, forward-error-scaled
test that returns the identical weights, so `nsc_rn_pg` imports on this
Windows interpreter. No operator value, calibration hash, or provenance
changed.

## Observer clocks and the finite room diagnostic

On `ds^2=N^2 dt^2-(dr+beta dt)^2-r^2 dOmega^2`, the PG normal is
$u=N^{-1}(\partial_t-\beta\partial_r)$. Its transported clock obeys
$dr/dt=-\beta$, $d\tau/dt=N$, with tangents
$d(\delta r)/dt=-\delta\beta-\beta_r\delta r$ and
$d(\delta\tau)/dt=\delta N+N_r\delta r$. A fixed-radius static clock
instead requires $N^2-\beta^2>0$ and rate $\sqrt{N^2-\beta^2}$.
The old saved normal scalar integrates the lapse at fixed excision radius;
it is not this transported proper time or a child clock. The saved outer
Killing scalar's zero at non-timelike samples is a clamp. The missing
four RK4 metric stages prevent recovering the old coupled observer's proper
clock exactly. The observer adapter refuses a domain exit or non-timelike
static step and keeps the last valid state rather than extrapolating.

The separately recorded finite-resolution diagnostic defines
`omega_zero=N/r`, `omega_cutoff=(N+|beta|)/spacing` and
`Q=omega_zero/omega_cutoff`. The six-case inner `Q` ranges from
`3.9299e-03` to `1.6829e-02`, with inner/outer ratios `18.858–19.909`.
The three deeper-domain final values are `4.1734e-03`, `9.9791e-03`,
`2.2380e-02`, with ratios `21.084`, `22.517`, `25.161`.
These are per-room measurements tied to the stored states and numerical
spacing, invariant under the stated common coordinate rescaling. They do
not determine Planck's constant, absolute zero, or a universal physical cutoff.
The source family remained incomplete, without a `summary.json`.
