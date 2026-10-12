# FGC-1-RSP1-FRZ1: amplitude-three resolution-spectrum freeze

`FGC-1-RSP1-FRZ1` is the pre-trajectory authorization for one bounded GR-0
resolution study. It does not repeat the full PROTO11 calibration and it does
not open SGB-L or FGC-QR. It asks only whether the evolved amplitude-`3`
`phi/Lambda` derivative tail obeys the already-frozen direct `<1/4`
contraction rule when the ladder is shifted to one genuinely finer grid.

The machine decisions are:

```text
RSP1_runtime_implemented = true
RSP1_execution_authorized = true
amplitude_three_spectral_veto_cleared = false
PROTO12_frozen = false
fresh_GR0_calibration_completed = false
classical_spherical_diagnostic_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

These decisions authorize the experiment, not its result.

## Why this study exists

[CAL8/PREF10](fgc-cal8-pref10.md) found two independent obstructions. The
amplitude-`5/2` run later reached a source-arithmetic floor. SRC2/PREF11 and
[SRC3](fgc-src3.md) localized and repaired that numerical instrument without
changing the equations or the strict `1e-12` source gate. Amplitude `3` never
encountered that source obstruction. It was rejected because the direct
coarse-to-medium derivative-tail ratio for `phi/Lambda` exceeded `1/4`, even
though the then-finest pair passed.

The historical evolved ratios on `(2049, 4097, 8193)` were:

```text
RK4:     (0.326717281113685,   0.22812049822366257)
SSPRK3:  (0.26675704013240803, 0.2142551901915683)
```

RSP1 therefore adds `16385` rather than changing a threshold. Its frozen
ladder is:

```text
4097 -> 8193 -> 16385
```

Both adjacent `phi/Lambda` derivative-tail ratios must be strictly below
`1/4` for both RK4 and SSPRK3 at the single synchronized evolved endpoint
`t=1/16`. No saturation, epsilon floor, guarded reinterpretation, or parameter
tuning is permitted for the target ratios.

## Prospectively frozen inputs

The authorization rebuilds all six combinations of two original numerical
methods and three grids from the same amplitude-`3` matter pulse and regulator
seed. It applies PROTO11's reference-balanced `q` map and SRC3's
reference-covariant source evaluator consistently on every grid. The
pointwise six-by-six affine systems are processed in fixed batches of 2,048 to
bound memory; batching changes no equation and no root.

Every initial source solve passes the unchanged strict residual gate. The
complete physical, gauge, and reduction constraints, medium/fine absolute
spectral budgets, and full-profile contraction remain vetoes. Existing CFL,
retry, boundary, checkpoint, and health rules are inherited unchanged.

## Time-zero evidence is public but is not the endpoint

At `t=0`, both methods pass their constraint, absolute-budget, and
whole-profile premises. The direct target ratios are nevertheless above
`1/4` on both adjacent pairs:

```text
RK4:     (0.8797163045076191, 0.3459398350356097)
SSPRK3:  (0.8846382581520159, 0.3459231990125775)
```

Those failures are serialized rather than hidden: they are not the evolved endpoint.
They do not pre-answer the
frozen question, which concerns the evolved endpoint at `t=1/16`. Treating
the initial ratios as the terminal outcome would make the proposed evolution
logically impossible before it was run.

## Instrument controls

The authorization requires:

- bitwise-stationary spherical Minkowski evolution for both methods;
- bitwise equality between bounded-batch and direct SRC3 pointwise solves on
  a nontrivial 128-point control;
- visibility of a one-bit input mutation;
- typed rejection of an invalid batch, non-centre grid, nonpositive source
  threshold, and promoted physical claim;
- immutable hashes for the plan, all six inputs, implementation, SRC3 result,
  and CAL8 diagnosis;
- an absent or empty fresh raw output namespace without creating it.

The tracked authorization contains no evolved RSP1 outcome.

## Outcome boundary

A positive RSP1 result can clear only CAL8's amplitude-three spectral veto for
successor design. It is not an eligible GR-0 calibration, collapse result,
candidate result, or mechanism result. A negative result preserves the veto
and must identify which frozen pair/method or premise failed. A runtime stop
or invalid/nonconverged execution is neither a positive nor a scientific
negative result.

In every case RSP1 says nothing about regulator activation, trapped-region
defocusing, singularity resolution, a daughter domain, a dark-sector
mechanism, varying locally measured light speed, retained-EFT validity, or a
physical model of nature.

## Reproduction

Before the trajectory is launched:

```bash
python3 scripts/reproduce_fgc_rsp1_frz1.py
python3 scripts/reproduce_fgc_rsp1_frz1.py --check
python3 scripts/run_fgc_rsp1_resolution_study.py --check-authorization-only
```

The authorization commands create no run output. The separately authorized
runner may then populate the ignored namespace exactly once.
