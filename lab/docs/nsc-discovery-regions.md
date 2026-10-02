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

On a contour, the implicit level motion is

```
xdot = -ρ_t / ρ_x
```

when `|ρ_x|` is at least `10^{-3}` times the dominant peak divided by the
node spacing. `ρ` is the total probability per proper length. `ρ_t` is the
nodal derivative of that density along the actual carrier rate, including
`∂t(rQ)`. `ρ_x` is the owned periodic derivative, whose Nyquist symbol is
zero, evaluated with the same Fourier interpolant. A shallower slope is
reported as `ambiguous` and is not divided.

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
