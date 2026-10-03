# Parent response adapter

Reusable directional adapter for the leading six-field chart. It does not
run a production episode and it does not record a held nonlinear
measurement.

The executable owner is
[nsc_discovery_parent_response.py](../src/recursive_horizons/nsc_discovery_parent_response.py).
Checks live in
[test_nsc_discovery_parent_response.py](../tests/test_nsc_discovery_parent_response.py).
The driver is
[derive_nsc_discovery_parent_response.py](../scripts/derive_nsc_discovery_parent_response.py).

## State and windows

The state is `Q`, `r`, `p_Q`, `p_r`, `phi0`, `phi1`. Canonical momenta are
`π = Δx_g W.T p`. `W` may be the identity. The geometry frame is square, so
the full geometry band stays active, and the fermion columns stay on the
ambient antiperiodic band.

The carrier centre is `L/2` (`4` when `L=8`). `s` is the periodic signed
distance from that centre. The windows are `|s|<=0.5` (child), `|s|<=1`
(protected collar), `|s|<=3` (parent) and `1.2<=|s|<=3` (parent annulus).

`ParentHolder` carries `.grid`, `.geometry_map`, `.weights`,
`.reference_columns`, `.source_columns`, `.source_metadata`,
`.geometry_metadata`, `.child_interval`, `.parent_interval` and
`.clock_locations`. Live occupations are the supplied weights, installed
with `dataclasses.replace` on `grid.fine`. The historical `build_pair`
six-column constructor is not called.

## What is complete

| Function | Role |
|---|---|
| `response_rates` | `leading.jvp` for the six fields, plus the existing force derivative for `δc` |
| `advance_tangent` | One RK4 step of the state and the tangent; centre clock integrated on the four stages |
| `centre_clock_variation` | `δO|τ = δO|t − Ȯ δτ / τ̇` through `matched_tau_correction` |
| `nonlinear_endpoint_clock_bracket` | Same formula from the nonlinear endpoints; it does not read the integrated linear `δτ` |
| `child_regional_content` | Primary readout: probability in the child window |
| `proper_localization` | Child proper length and child, collar and annulus probabilities |
| `fixed_mode_observer` | Thin `V† C V` on the fixed reference columns; not the spatial window |
| `explicit_normalization_derivative` | Analytic Hermitian projection when the frame is orthonormal |
| `preparation_tangent` | Two-step finite difference of the full T0 preparation |
| `hierarchy_rates` | Inherited `grandchild.route_rates` for widths `(2, 2, 4)`; variable rows otherwise |
| `full_vs_reduced` | One full RK4 step against streamed direct and sequential reductions |
| `evolve_prepared_source` | One admitted RK4 step of the fresh prepared source and its centre clock |
| `prepare`, `predict`, `measure`, `check` | Explicit stage contract. A missing seal stays open |

Source columns may be non-orthogonal. The fixed modal observer must be
orthonormal. Those are different objects.

`preparation_tangent` is labelled `finite_difference_full_preparation`.
The probes move the state, momenta, columns and weights. The recorded
components are the Hamilton constraint, the shift constraint, the canonical
momenta and the actual column norms. The radius solve is one component, not
the preparation. `analytic_preparation_jacobian` is false. A singular radius
Jacobian is returned as unavailable and is not replaced by a fabricated
analytic solve. The state derivative is kept.

## Reductions

Retained amplitudes have `R` rows, the first hierarchy block. The source
rank is the weight length and is not fixed at 6. There is one metric.

`full_vs_reduced` evaluates geometry with `leading.rates` on the
reconstructed columns. Column rates are one concatenated `apply_dirac` per
RK4 stage. The dense exterior propagator is not stored. Cross forces are
`force(full) − force(retained) − force(drive) − force(memory)`. Complement
forces are the nodal forces of the exterior columns.

Memory omission is labelled twice. Conditional omission zeros the memory
rate and keeps the exterior drive. Autonomous omission closes the retained
rows alone. Both residuals are against the full Dirac image. Neither is a
regeneration claim. The opening split requires a nonzero prepared cross,
`a diag(c) e_drive†`.

The coupled-memory rank-6 guard is not on this path. A rank other than 6
makes that guard raise; this adapter still builds the holder.

## Contract

`response_rates` and `advance_tangent` do not write. A missing binding, or
an episode path that is not a file, raises `BindingUnavailable` with the
stage left open. That is not a scientific pass.

`prepare` writes an exclusive record. Outside the repository this can be a
manufactured fixture. Under `lab/results/development` it also requires a
frozen producer commit, a physical OWNER1 or OWNER2 binding, and a matching
source closure. Input files and sealed evidence directories are not writable.
The seeded amplitude `0.0013` is the baseline source weight. It is not a
measurement veto.

`predict` locks the tangent forecast and admits the step through
`nsc_discovery_parent_step_control.step_restriction`. It copies a held
direction only when that direction is already sealed on the prepare record.
The seal chooses `population` or `gradient` and a relative change, such as
`+5%`, before any nonlinear measurement.

`measure` stays `open` when that seal is absent. When the seal is present it
evolves the fresh prepared source, and the sealed population or gradient
copy, by one admitted RK4 step at the matched centre clock. Streamed
full-to-reduced forces keep the complement. `scientific_pass` is not minted.
`check` requires prepare, then a locked forecast, then measure, and requires
the input hashes to be unchanged.

## Dependencies

- `nsc_discovery_leading_einstein`: state, `rates`, `rk4_step`, `jvp`, `encode`, `constraint_arrays`
- `nsc_discovery_response`: `matched_tau_correction`, moment, force, image and density directions
- `nsc_discovery_grandchild`: `initial_routes`, `reconstruct_route`, `route_rates` for widths `(2, 2, 4)`
- `nsc_discovery_coupled_memory`: split interpretation only; `_require_rank6` is not called
- `nsc_spherical_coupling`: `apply_dirac`, `source_from_columns`
- `nsc_discovery_backend`: FFT carrier classes, one worker; `make_fft_grid` is not called
- `nsc_discovery_extent`: real interval integral and periodic sample
- `nsc_discovery_parent_step_control`: `step_restriction`; the admission arithmetic is not copied
- OWNER1 `nsc_discovery_parent.prepare_parent`: callback only; default calls do not run it
- OWNER2 `nsc_discovery_parent_episode`: authenticated `load_parent_record(path)`
  returns `(record, arrays)`; `reconstruct_parent_pair(arrays, record)` rebuilds
  the carrier; `state_from_arrays(arrays)` reads canonical `pi_Q` and `pi_r`.
  The callback dictionary keeps `load_parent_preparation` as an alias of this
  same public loader. Callback discovery calls none of these functions.
- OWNER3 observables: not imported; the spatial window is not their fixed-mode observer

## Not in this adapter

An analytic preparation Jacobian, and a scientific pass. The current physical
stage advances one admitted step after freeze, seal and a physical binding.
A longer matched-clock held campaign still requires its episode binding and
bounded continuation; this stage is not evidence that such a campaign is
complete. This note's tests do not run that campaign.

## Reproduction

```sh
python3 scripts/lab.py -m pytest tests/test_nsc_discovery_parent_response.py -q
python3 scripts/lab.py scripts/derive_nsc_discovery_parent_response.py
```

The pytest command uses a manufactured rank-3 chart. The driver preview
prints the contract and writes nothing.
