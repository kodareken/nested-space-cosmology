# FGC-1-HLT10-MON10: PROTO12 runtime preflight and GR-0 authorization

`FGC-1-HLT10-MON10` is the pre-trajectory runtime gate for frozen PROTO12.
Its canonical machine record is `results/fgc-1-hlt10-mon10.json`.

The only positive decisions are:

```text
PROTO12_successor_runtime_implemented = true
PROTO12_fresh_GR0_dynamic_calibration_authorized = true
```

They mean that SRC4 and the pairwise spectral classifier have been bound to
the inherited GR-0 campaign, all frozen inputs and time-zero events have been
rebuilt, and one fresh outcome-neutral calibration may start. They do not say
that a trajectory has run, that a trapped sphere forms, that FGC-QR is healthy,
or that collapse defocuses.

## Immutable evidence boundary

HLT10 reads PROTO12/PRO12-FRZ1, HLT9, CAL8/PREF10, RSP1/PREF12, SRC4/VEC1,
ID2, and the CAL8 plan from immutable commit
`a42780179591012a29eb81a87f59be1f60f9e514`. Every configured blob hash and
every inherited implementation hash must still match the current checkout.

The fresh CAL9 plan is normalized back to immutable CAL8. Exact equality after
removing only the declared PROTO12 source, spectral, ownership, and namespace
delta proves that the action, initial data, amplitudes, grid ladder, methods,
CFL, dissipation, retry budgets, source threshold, constraint thresholds,
event schedule, trapped-sign rule, causal boundary, and physical stops were
not silently changed.

CAL8 and RSP1 remain disclosed diagnostic history. Neither old trajectory is
accepted as a PROTO12 outcome. CAL8 isolated a floating source-evaluator wall
and a separate coarse-pair spectral veto; RSP1 supplied qualification evidence
for the latter. PROTO12 must rebuild and reassess its own states.

## SRC4 runtime binding

The physical state and reference-balanced radial map remain PROTO11:

```text
q(0) = D_h(u-u_ref) + q_ref
q_r  = D_h(q-q_ref)
C_q  = q - [D_h(u-u_ref) + q_ref].
```

The source operator evaluates the same unredefined GR-0 ACT1/VAR1/REF1 affine
systems through SRC4's four-dimensional tensor contractions. Points are
processed in fixed batches of at most 2,048. Batching changes no pointwise
matrix, right-hand side, root branch, tolerance, or equation. The complete raw
residual remains subject to the strict `< 1e-12` gate, with at most 16 declared
refinement iterations and the unchanged kinetic condition limit.

The exact-reference zero-acceleration identity is still selected only after
the complete raw source passes and only under bitwise equality. A one-bit
perturbation takes the ordinary SRC4 path. SRC3's explicit-index evaluator is
retained as an independent differential oracle, not as a fallback selected by
the outcome.

## Pairwise spectral binding

For every field and derivative-weighted tail on both adjacent grid pairs,
direct contraction remains primary:

```text
tail_fine / tail_coarse < 1/4.
```

If a raw ratio fails, PROTO12 may label that pair diagnostically saturated only
when both pair members pass their original absolute budgets, both members'
field, derivative/Nyquist, and top-band-erasure witnesses lie within their own
measured map scales, and the complete three-grid profile and round-trip
residuals contract. Every raw power fraction, ratio, map residual, tail
amplitude, erasure witness, and classification remains serialized.

The round-trip residual is a measured property of the diagnostic map. It is
not a continuum-error bound. Diagnostic saturation is not physical resolution
and is not evidence for the gradient mechanism.

## Rebuilt inputs and attacks

HLT10 reconstructs:

```text
2 amplitudes x 2 methods x (2049, 4097, 8193) = 12 inputs.
```

All twelve `u`, `p`, and `q` arrays must match immutable HLT9 bitwise; all
twelve must pass the strict SRC4 source precheck; and the four amplitude/method
time-zero ladders must pass the complete PROTO12 common event.

The adversarial controls require:

- bitwise stationarity of exact spherical Minkowski on both methods and all
  three grids;
- ordinary source routing after a one-bit perturbation;
- bitwise SRC4/SRC3 agreement on a nontrivial physical sample;
- visible generic reduction and source-derivative defects;
- bitwise agreement between common-event reduction rows and the PROTO11 map;
- equality of every PROTO12 raw ratio with its predecessor raw ratio;
- a demonstrated raw adjacent-pair failure at time zero that remains public
  while the fully guarded route is exercised, with CAL8's later evolved
  coarse-pair failure retained separately as immutable history;
- fail-closed rejection when any saturation guard, threshold, source backend,
  batching contract, branch, claim, or physical input is broadened.

## Decision and nonclaims

HLT10 authorizes only one fresh GR-0 calibration in
`runs/fgc-2-sf1/proto12/calibration`. The holdout namespace remains empty and
closed. In particular:

```text
PROTO12_resolved_holdout_manifest_authorized = false
classical_spherical_diagnostic_authorized = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

HLT10 is therefore an instrument-and-protocol result. It is not evidence for
defocusing, singularity resolution, a child domain, a dark-sector mechanism,
particle ontology, or varying locally measured light speed.

## Reproduction

Before the first trajectory:

```bash
python3 scripts/reproduce_fgc_hlt10_mon10.py
python3 scripts/reproduce_fgc_hlt10_mon10.py --check
python3 scripts/run_fgc_gr0_calibration_v12.py --check-authorization-only
```

The authorization command must leave both PROTO12 output roots absent or
empty. After a calibration exists, canonical verification reuses the stored
pre-launch namespace observation rather than claiming the run never occurred.
