# FGC-1-CAL6-PREF8: proper-spectrum conditioning diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced, pre-trajectory numerical-diagnostic result;
no mechanism result

```text
PROTO9_preflight_conditioning_diagnosis_completed = true
PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence = false
PROTO10_spectral_contract_revision_required = true
PROTO10_fresh_GR0_dynamic_calibration_authorized = false
```

## Question

HLT7 directly constructed the complete `2049 -> 4097 -> 8193` ladder for
both calibration amplitudes and both numerical methods. All twelve raw source
checks passed `1e-12`; all four common events passed the unchanged constraint
admission; and every medium/fine absolute spectrum budget passed. The only
veto was the unchanged requirement that every adjacent top-band power ratio
be below `1/4`.

This certificate asks a narrower question before changing that contract:

> Is the failed `4097 -> 8193` ratio a resolved numerical observable, or is
> it being evaluated below the demonstrated conditioning scale of the
> proper-grid interpolation/FFT map?

No trajectory, SGB-L state, FGC-QR state, trapped-sphere outcome, or
Raychaudhuri outcome is read.

## Immutable input boundary

The diagnosis consumes HLT7 and its run plan directly from immutable commit
`234c5a5cc2c0b6bdc16d22146cdb9a9c6ae5b98d`. It reconstructs all twelve
initial states and requires every projected-state hash to match HLT7. The raw
PROTO9 power fractions, ratios, and four failed admissions must reproduce
without reinterpretation before the conditioning calculation is allowed.

The PROTO9 calibration and holdout namespaces remain absent or empty. The
diagnosis neither creates nor reads a trajectory namespace.

## What “diagnostically saturated” means

Let `T` map native-grid samples to the uniform proper-distance FFT grid and
let `S` map those samples back. HLT7 already records

```text
delta_f = ||S T f - f||_infinity.
```

This is a measured non-idempotence of the diagnostic map. It is explicitly
not asserted to be a bound on the unknown continuum field.

For each field, the new diagnosis also performs an orthogonal FFT operation
on the windowed diagnostic input: it sets the declared top-band coefficients
to zero, transforms back, and records the exact infinity-norm perturbation
needed to erase that band. Call this `epsilon_tail`.

A failed finest-pair ratio is classified as diagnostically saturated only
when all of the following hold:

1. the original coarse-to-medium field and derivative ratios pass `< 1/4`
   directly;
2. the medium and fine grids pass every unchanged absolute spectrum budget;
3. `epsilon_tail <= delta_f` on both the medium and fine grids;
4. the round-trip residual contracts strictly on both refinements, or the
   field and residual are exactly zero;
5. the complete resampled profile difference contracts strictly from the
   first grid pair to the second in both infinity and RMS norm; and
6. every raw fraction, ratio, residual, norm, and erasure witness remains
   public.

There is no epsilon floor, fitted factor, relaxed `1/4`, or changed absolute
ceiling. The direct route remains primary. The saturation route applies only
to the finest pair and adds whole-profile and map-conditioning guards.

## Result

The raw result remains exactly HLT7's result: all four PROTO9 admissions fail.
Across the four cases there are 48 finest-pair field/derivative decisions:
32 pass `< 1/4` directly and 16 fail it.

All 16 failures satisfy the conditioning witness. For the failing fields,
erasing the fine-grid top band requires only about `0.11%` to `0.56%` of the
measured fine-grid round-trip discrepancy. The ratio therefore changes its
pass/fail class inside a much smaller diagnostic-input perturbation than the
map's own observed non-idempotence.

The conclusion is reinforced by independent checks:

- every nonzero complete profile contracts from the first grid pair to the
  second in both infinity and RMS norm;
- every nonzero round-trip residual contracts on both adjacent pairs;
- the worst medium/fine field-power fraction remains more than eight hundred
  thousand times below the unchanged field ceiling; and
- the worst medium/fine derivative-power fraction remains more than
  seventy-four thousand times below the unchanged derivative ceiling.

Under the prospective resolved-or-saturated rule, all four `t=0` cases pass.
That prospective pass is design evidence only. It does not retroactively
authorize PROTO9, authorizes no campaign, and is not a calibration result.

## Smallest justified successor

The evidence supports a separately frozen `FGC-2-SF1-PROTO10` that changes
only the interpretation of the finest-pair nested tail:

```text
coarse -> medium: unchanged direct ratio veto
medium -> fine:   direct ratio pass
                  OR
                  conditioning + whole-profile saturation witness
```

PROTO10 must retain the `2049 -> 4097 -> 8193` ladder, action, physical
inputs, methods, equations, raw `1e-12` source tolerance, CFL, retries,
constraint ownership and thresholds, absolute spectrum budgets, direct
`1/4` ceiling, boundary and health stops, trapped-sign rule, affine-null
observable, outcome language, and all physical nonclaims. A new namespace,
freeze record, and runtime certificate are required before any fresh campaign.

## Scientific boundary

This result diagnoses an ill-conditioned ratio-only test contract. It does
not prove that the initial data are continuum-exact, that future evolved
states will satisfy the successor rule, or that diagnostic saturation equals
physical resolution. It does not test regulator activation, collapse-induced
defocusing, singularity resolution, a child domain, a dark-sector mechanism,
or a varying locally measured speed of light. It supplies no retained-EFT or
physical-model authorization.

## Reproduction

```bash
python3 scripts/reproduce_fgc_cal6_pref8.py --check
```
