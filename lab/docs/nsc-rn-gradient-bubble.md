# RN radial and transverse pressure diagnostics

## Semantic correction, 10 October 2026

The historical “two-force gradient balance” measures radial and transverse
orthonormal **pressures in one recorded state**. It does not identify forces
from separate parent and child regions. Their equal-magnitude crossings are
retained finite diagnostics, not a mechanical-equilibrium condition or a
cosmological boundary. This correction changes the interpretation and live
consumer, without evolving a state or altering a historical result file.

The v1–v5 criterion `g >= -1` and its attribution to Hayward are withdrawn.
[Hayward, gr-qc/0506126](https://arxiv.org/html/gr-qc/0506126v2), equations
3, 7 and 8, supplies a regular-core metric and finite stress components.
It does not supply that generic pressure-slope criterion. For an assumed
nonzero pure power `p ~ r^g`, pointwise finiteness at zero would require
`g >= 0`; finiteness of `integral p dr` instead requires `g > -1` (the
endpoint is logarithmically divergent). Neither is the finite-grid
measurement performed here. No replacement centre-regularity test is imposed.

Reaching the positive RN inner-horizon radius does not reach the centre.
A slope sampled outside or at that horizon cannot establish the behaviour
as `r -> 0`. The v6 consumer therefore never emits `interior_not_regular`
or `throat_extrapolation_diverges`. An inward continuation is not selected
as a prerequisite by this diagnostic.

## Quantities and domain

Inputs are saved station payloads of `nsc-rn-neutral-family-v1`: the areal
grid, `radial_pressure`, `angular_pressure` and same-state reconstructed
`R4`. The grid runs from the recorded inner edge to `32 r_m` and has
strictly positive radius. No second geometry, new evolution or fit is used.

`nsc_rn_source.source_evaluate` defines the pressures as

```text
physical_pr = (force_N - angular_force) / (4 pi volume)
physical_pa = -force_r / (2 lapse radial_density r 4 pi)
volume = radial_density r^2
```

The diagnostic preserves `B = |p_r| / (|p_r| + |p_Omega|)` and the central
log-log differences `d ln|p| / d ln r`. A crossing of `B = 1/2` equates
pressure magnitudes; it does not include stress gradients or the metric
terms needed for a mechanical balance equation. Magnitudes also discard
the signs, so equal signed stresses are not implied.

The old “integrated force” is `4 pi integral |p| r^2 dr`, a
volume-weighted pressure integral, not a surface integral. The new
`radial_pressure_volume_integral` and `transverse_pressure_volume_integral`
aliases give the same values. This chart has flat spatial radial measure;
the expression must not silently be reused as a general proper-volume rule.

`R4` and pressures share state provenance. `Ricci2` in these payloads is an
analytic source-free reference and remains rejected by the consumer.
The old `max|R4| / max|p_r|` diagnostic is retained in the recorded units;
common provenance alone is not a new force or regularity criterion.

## Retained historical measurements

The immutable [v5 record](../results/development/nsc-rn-gradient-bubble-v5.json)
covers six cases and 177 stations. It reports volume-weighted transverse to
radial pressure ratios of `0.156` to `0.820`, local radial log-slopes of
`-7.37` to `-152.64`, and pressure crossings at 49.7 percent of stations.
These numbers are retained, while their former exterior-force, collapse
and regularity interpretations are withdrawn.

| Case | Stations with a crossing | Recorded crossing-radius range |
|---|---|---|
| `mu-1.002-eps-0.001` | 65% | 1.073 to 2.054 |
| `mu-1.002-eps-0.010` | 65% | 1.073 to 2.054 |
| `mu-1.010-eps-0.001` | 69% | 1.045 to 2.055 |
| `mu-1.010-eps-0.010` | 61% | 1.038 to 2.055 |
| `mu-1.040-eps-0.001` | 32% | 1.290 to 2.039 |
| `mu-1.040-eps-0.010` | 30% | 1.290 to 2.039 |

No case wrote a `summary.json`. Four stopped on CPU projection; two
`mu=1.002` cases stopped at `t=8` without a failure record. Their cause
remains unknown. The observed station fractions describe this incomplete
family, not a completed evolution or sustained cosmological structure.

The saved stencil audit agrees with exact-rational weights and reproduces
the sampled slopes. This supports the numerical stencil comparison only;
it does not establish physical central divergence. A gauge-dependent
spinor crossing correlation likewise cannot establish which stress is
“physically driven.” No new numerical audit was run for this correction.

## Compatibility and owners

The [consumer](../src/recursive_horizons/nsc_rn_gradient_bubble.py) now has
API version `nsc-rn-gradient-bubble-v6`. `diagnostic_semantics` travels with
its outputs. Force/wall field names and legacy verdict strings remain
available to existing readers, but describe finite pressure diagnostics.
`single_force` means the declared small-transverse-pressure threshold;
`curvature_dominates` remains a conventional diagnostic threshold. Neither
is a physical exclusion theorem. The steep-slope flag is reporting only.

`slope_meets_regularity_limit`, `regularity_slope_limit` and
`regularity_slope_tolerance` are deprecated and return `None`.
`regularity_decidable_on_this_domain` is deprecated and always `False`.
The old module constants remain importable but are unused. Consumers must
not coerce a missing criterion into a failed physical test. Horizon reach
and an explicitly assumed power-law extrapolation remain separately reported.

The [driver](../scripts/derive_nsc_rn_gradient_bubble.py) targets a new,
creation-only `nsc-rn-gradient-bubble-v6.json` and links v5 as its predecessor.
It was not run for this correction; no v6 family result is claimed. All
v1–v5 records remain intact. The [tests](../tests/test_nsc_rn_gradient_bubble.py)
check unchanged finite pressure numerics, unsupported-criterion handling,
horizon/centre separation and historical-record provenance separately.

The distinct boundary-traction calculation in
[nsc-discovery-parent-traction.md](nsc-discovery-parent-traction.md) is not
identified with this pressure-component comparison. Neither diagnostic
requires becoming the atlas's sole geometric model.
