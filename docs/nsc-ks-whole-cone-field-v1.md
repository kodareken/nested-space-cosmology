# Current-history whole-cone field campaign, v1

This successor is a resumable selected-source driver for representative family
`[14, 1]`. It binds history `0b0e4ced…`, restores the authenticated v4
local-Fourier coefficient payload by hash, and re-evolves the unchanged
retained upstream source selected by
`derive_nsc_ks_source_control_v2.select_source`. Streamed DOP853 segments
are enclosed as history-minus-reference cells through `WholeConeAccumulator`.

The reusable API is `src/recursive_horizons/nsc_ks_current_field_campaign.py`.
The executable owner is `scripts/run_nsc_ks_whole_cone_field_v1.py`.
The accumulator and checkpoint bytes remain those of cone v1. This driver
does not replace the historical v4 one-cell record.

## Production bindings

| Object | Role |
|---|---|
| Declared history | Profile `0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0` |
| Source | Unchanged `RetainedUpstreamArchive` rows for family `[14, 1]` |
| Integrator | DOP853, `rtol=5e-14`, `atol=5e-19`, `max_step=1/8192`, joint step control, zero tangents |
| Profile payload | Authenticated v4 `profile_coefficient_payload`; never recomputed |
| Stream | `stream_ks_trajectory` segments from the same start |
| Enclosure | History-minus-reference before norms; residual integrals are not multiplied by a second cell width |
| Checkpoint | Atomic rewrite after every eight newly enclosed accepted segments |

On resume the driver re-streams from the fixed upstream start and authenticates
every stored dense-segment identity before enclosing the next cell. A changed
source binding, settings digest, payload digest, prefix digest or segment hash
is rejected.

## Selected-source coverage and the remaining family gate

`all_history_cells` is true only when the stream reaches `rho=1` and every
accepted segment has a matching enclosed cell, with contiguous indices
`0 … n-1`. `selected_source_trajectory_complete` additionally requires the
declared solver/profile bindings and `stream_ks_trajectory`. The reused
`select_source` currently selects four energy rows, twelve coherent columns.
It does not cover the entire family. Therefore `family_14_1_scientific_run`
and `full_family_source_coverage` remain false even at temporal completion.
Partial, interrupted, method-failed,
resource-stopped and max-new-cells diagnostic output remains `OPEN` and cannot
be a completed family record, a feasibility verdict, or a physical certificate.

Physical `rho=1` source accuracy stays `None`. Missing Gronwall inputs stay
named nulls; they are not filled with zero. Runtime, memory and remaining-work
forecasts are excluded from the scientific digest.

## Continuous coefficient and norm inputs

Owned inputs reused from existing analytic owners:

- `K0` / off-diagonal integrals from `value_integral_bounds`
- `B_z` from `pure_radius_commutator_integrals`
- radius and axial lowers from `radius_bounds`
- source operator norm from `source_norm_upper`
- signed-family multiplicity from retained copies and degeneracy
- maximum source energy from the selected labels
- exact zero initial `D` from the difference-envelope initialization
- residual integrals from enclosed cells
- `M_integral` and `Mz_integral` from the continuous radius-coupling bounds
- `reference_residual` from anchored reference defects on every streamed cell;
  replay recomputes it across the authenticated checkpoint prefix
- zero `reference_initial` only for the numerical subproblem with the supplied
  binary columns exact; physical preparation error stays separate and unknown

Explicit named nulls, not invented zeros:

- physical `rho=1` source error
- reference residual for an incomplete time prefix (never promoted to a full
  error bound by the campaign)
- continuous sup-norms of `delta F` and the resulting contraction errors

`propagate_difference_error` and `difference_matter_error` run only when every
required slot is owned. The reference method reuses the background Taylor
and Bernstein bounds with exact ball conversion of the anchored polynomial.
It reverses the background time coordinate for decreasing-rho evolution.
The saved cell-122 reference integral is bounded by `2.8866449342721377e-16`;
this is a one-cell numerical bound, not the full field error.

```sh
PYTHONPATH=src .venv/validation/bin/python scripts/derive_nsc_ks_reference_residual_pilot.py --check
```

## Diagnostic and resource stops

`--diagnostic-max-new-cells N` encloses at most `N` new cells, writes a valid
checkpoint and a non-scientific runtime forecast, and cannot emit
`family_14_1_scientific_run`. A free-space floor of 60 GiB and an optional CPU
budget are resource stops with the same rule. Do not launch
`--run-family-14-1` unless the full scientific run is intended.

Reproduce the driver tests with:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_nsc_ks_current_field_campaign.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_nsc_ks_current_field_cone.py tests/test_nsc_ks_local_profile_fourier_bound.py -k "not v4_cell"
```
