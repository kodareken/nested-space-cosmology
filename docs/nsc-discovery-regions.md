# Stage-4 discovery regions

`analyze(pair, state, time)` in
[nsc_discovery_regions.py](../src/recursive_horizons/nsc_discovery_regions.py)
is a postprocessing row for one saved or native nested state. It does not
step the state, call `initial_state`, rebuild `W`, or reset the source
columns. The rate it uses is the existing conformal carrier rate. The
normal-energy slope is the analytic polarization of `F_L/r` owned by
`normal_energy_slope` in
[derive_nsc_nested_parent_child.py](../scripts/derive_nsc_nested_parent_child.py).

```python
from recursive_horizons.nsc_discovery_regions import analyze, read_station_regions

row = analyze(pair, state, time)
stations = read_station_regions(directory, case_id)
```

`control_mode="frozen_geometry"` keeps the current fields and sets the
geometry jets, including `Qdot` and `rdot`, to zero. Field-column rates and
the polarized source kernels `K` and `S` stay. The default is `coupled`.

## Three localizations

These are different questions about the same state.

The six-column total is the occupation-weighted probability of all six
source columns. Its proper width is the same proper-arc packet width the
stage-1 observer reports for that total mass. `q` is not multiplied into
the measure, and the width is not a dispersal time.

Columns `(2, 3)` are the originally child-prepared source pair. Each column
pair has its own probability content, proper width, and retained fraction
inside its preparation support. The support is the recorded
`source_support_closures` when that record has three intervals, and the
separated supports `(0, 1)`, `(1, 3)`, `(3, 4)` otherwise. The ancestry
string is a computational label. The pair is not an independent energy
parcel: the normal shell below is the one `F_L/r` shell.

The gradient boundaries are the connected contours of the total
proper-length probability density at one quarter, one half, and three
quarters of its dominant peak. A contour endpoint is a root of the Fourier
interpolant. The comparison contour keeps Fourier modes with
`|k| ≤ n/4`. The absolute circular distance between those endpoints is the
contour uncertainty. `contour_certified` stays false.

Flags on the native density:

| Flag | Meaning |
|---|---|
| `flat` | Peak contrast below `1e-4`, or the peak is not positive |
| `multiple_dominant_peaks` | More than one local maximum at or above half the dominant peak |
| `merge` | Two or more of those maxima lie in the dominant half-maximum component |
| `split` | The half-maximum superlevel has more than one connected component |

`flat` leaves the other three flags false and does not invent endpoints.
`merge` and `split` can both be true: the dominant half-maximum component can
contain more than one peak while another component remains elsewhere.

## Contour velocity

The peak is the maximum of the current nodal samples. It is not moved off a
node to raise the Fourier interpolant. The contour level is that peak times
one quarter, one half, or three quarters. Differentiating
`ρ(x(t), t) = fraction * peak(t)` gives

```
xdot = (fraction * peak_rate - ρ_t) / ρ_x
```

`peak_rate` is `ρ_t` at the first argmax. Samples within roundoff of that
maximum are the same peak when their rates also agree within roundoff, and
that same nodal rate is used. Incompatible tied rates stay `ambiguous` and
are not divided. A slope shallower than `10^{-3}` times the peak divided by
the node spacing is still refused before the division. `ρ` is the total
probability per proper length. `ρ_t` is the nodal derivative of that density
along the actual carrier rate, including `∂t(rQ)`. `ρ_x` is the owned
periodic derivative, whose Nyquist symbol is zero. Off-node samples and the
moving-ledger integrals use the real Fourier evaluator, with the Nyquist mode
kept as a cosine.

The row records the peak rate, the level rate `fraction * peak_rate`, the
status, and this convention. `continuum_optimized_peak` stays false.

## Sealed v1 speeds

`results/development/nsc-discovery-regions-nf256-v1.json` and
`results/development/nsc-discovery-regions-nf512-v1.json` keep the speeds and
moving half-maximum ledgers produced with `xdot = -ρ_t / ρ_x`. Those bytes
are qualified and are not healed. Fixed windows still use zero boundary
velocity and the same Reynolds identity. A creation-only successor replays
original station chunks with the nodal peak rate:

```sh
python scripts/lab.py scripts/derive_nsc_discovery_regions_successor.py --check
python scripts/lab.py scripts/derive_nsc_discovery_regions_successor.py --read <episode-directory> --case <case-id> --output <new.json>
```

`--check` hashes the current source closure and those two sealed records. It
does not parse them and it does not write. `--output` is exclusive, stays
outside the episode directory, and refuses a file above 64 MiB. Before that
file is created, a case and ordinal that appears in either sealed record must
match the stored chunk hash. A mismatch raises and writes nothing. The stored
v1 speeds are not copied or healed. Root freezes the producer bytes before a
production replay. This note does not run that replay.

## Leading parent observer

`observe(pair, state, time, control_mode="coupled", bundle=None)` in
[nsc_discovery_parent_observer.py](../src/recursive_horizons/nsc_discovery_parent_observer.py)
is the variable-rank adapter for the leading Einstein diagnostic. A passed
bundle must carry that owner's `leading_rate`, source, fine state, and fine
system. The quadrature state flow is `prolong(decode(leading_rate))` for the
geometry and the AP columns. Raw `unprojected_rates` remain a projection
diagnostic. Probability rate, pressure, the analytic normal-energy slope, and
contour speeds use the projected flow. Passed `metric_jets` must match that
same prolonged rate. Nothing in this adapter calls the auxiliary chi carrier.

`control_mode="source_free"` keeps that label. The geometry rate is the
coupled leading rate, so the metric still moves. Occupation weights may be
zero. Zero total occupation reports flat matter curves and does not invent a
packet. A nonfinite column is refused. The episode control-mode normalizer
is not patched.

The holder exposes `.grid`, `.geometry_map`, `.weights`,
`.reference_columns`, `.source_columns`, `.source_metadata`,
`.geometry_metadata`, `.child_interval`, `.parent_interval`, and
`.clock_locations`. `W` may be the identity. Fine occupations are the
holder's actual weights. The centre is `L/2`. Signed distance `s` from that
centre splits the carrier into disjoint measurement cuts:

| Region | Signed distance | Accounted |
|---|---|---|
| child | `|s| <= 0.5` | one interval |
| inner | `0.5 < |s| <= 1.2` | two intervals |
| parent annulus | `1.2 < |s| <= 3` | two intervals |
| ambient | `3 < |s| <= L/2` | wrapped across the branch cut |

The protected collar `|s| <= 1` overlaps child and inner, so it is not a
fifth window and it is not given piston work. The declared parent interval
`|s| <= 3` contains the child; the accounted parent channel is only the
annulus. If the global nodal peak lies in that annulus, the child nodal peak
is still reported on its own and `dominance_forced` stays false. Global
contours keep the global peak.

The returned row is JSON. Its channels reuse the fixed and wrapped region
ledgers: flux, pressure, lapse, and the finite projection defect, each
divided by `dx` once. Metric channels are the passed leading jets: normal
4D tides, `R`, and `W`, with positive proper probability. Proper clocks are
`rQ` at the holder locations and at the cuts. Null rays use coordinate speed
`±1` on the conformal chart. Column ancestry is not an energy parcel.

## Moving normal energy

The shell summand is nodal `F_L/r`. The density is that summand divided by
`dx` once:

```
d/dt ∫_a^b η dx = η(b) ḃ − η(a) ȧ + F(a) − F(b) + pressure + lapse
                  + finite projection defect.
```

`F(a) − F(b)` is the collocation flux `(sample[a] − sample[b]) / dx`.
Pressure and lapse are the existing volumetric integrands. The finite
projection defect is the window integral of

```
slope + ∂x(flux summand) − pressure summand − lapse summand
```

divided by `dx` once, plus the gap between the collocation flux and the
zero-Nyquist derivative. That gap is kept. The derivative matrix is not
altered to remove it.

A fixed window has `ȧ = ḃ = 0`. Its rate equals the window integral of the
analytic slope. A moving or wrapped contour uses the same identity. A
wrapped component is split at the periodic branch point; that artificial
cut is fixed and its two contributions cancel. The measurement cut does
not add piston work `−pV`, and `physical_wall` stays false.

`normal_energy.analytical_continuity_closure_residual_max` is the absolute
difference between this sum and `∫ ∂t η dx + η(b)ḃ − η(a)ȧ`.

## Reader and command

`read_station_regions` accepts `snapshot_kind="station"` chunks through
`nsc_discovery_episode.read_station_states` and `load_checkpoint`. It
reloads the stored `W`, columns, and `Φ`. It does not replace `Φ` by the
prepared source, and it hashes each chunk before and after the row. A
changed hash raises.

```sh
python scripts/lab.py scripts/derive_nsc_discovery_regions.py --check
python scripts/lab.py scripts/derive_nsc_discovery_regions.py --read <episode-directory> --case <case-id>
python scripts/lab.py scripts/derive_nsc_discovery_regions.py --read <episode-directory> --case <case-id> --output <new.json>
```

`--check` records producer hashes, settings, and sealed-input hashes. It
does not parse those npz payloads and it does not write a file. `--output`
uses exclusive creation and must sit outside the episode directory. An
existing output path is refused. The pinned check is

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_regions.py -q
```

from the repository root. The tests are manufactured. They do not read a
production trajectory.

## What renewal classification would still require

This row does not classify renewal. `renewal_asserted`,
`continuum_certified`, and `contour_certified` stay false. A later
classification would need all of the following.

1. Saved station states at the handoff and at the later stations under
   comparison, for the coupled cases and the frozen-geometry controls, at
   each prepared resolution.
2. The quarter, half, and three-quarter contour endpoints and implicit
   velocities on those stations, with the resolution gaps from this row,
   and a pre-registered tolerance for those gaps.
3. A pre-registered distinction between a maintained packet and a
   regenerated boundary. Motion of a measurement cut is not that
   distinction, because the cut carries no wall traction.
4. The stage-1 null characteristic tracks set beside the measured contour
   velocities at the same stations.
5. The inherited local response of the realized source inside the measured
   contour. That response is not a region integral.
6. A continuum or refinement certificate. This row does not estimate a
   state-error bound.

Modal complement content is the probability in the columns other than
`(2, 3)`. Spatial exterior probability is the parent annulus of the total
density. They are stored as different numbers.
`modal_complement_is_spatial_exterior` is false.

## Domain

No new inherited-law evolution is defined here. An imposed contour
trajectory is not regeneration evidence. The fixed coordinate windows
`(1, 3)` and `(0, 4)` remain the prepared windows; they are not the
gradient contours and not the six-column width. The root index was not
edited with this note.
