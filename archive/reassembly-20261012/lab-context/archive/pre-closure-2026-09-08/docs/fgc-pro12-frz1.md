# FGC-1-PRO12-FRZ1: source-instrument and pairwise-spectrum freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** immutable numerical-premise protocol boundary; no PROTO12 runtime,
fresh calibration, candidate trajectory, or mechanism result

```text
PROTO12_frozen = true
PROTO12_successor_runtime_implemented = false
PROTO12_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

## What this certificate resolves

PROTO11 did two useful things without producing an eligible calibration case.
It proved that both original amplitudes could reach evolved common events, and
it separated two numerical obstructions that had previously been entangled:

1. the complete affine source evaluator crossed the unchanged `1e-12` gate
   through binary64 cancellation even though the exact affine system remained
   nonsingular; and
2. the original three-grid calibration ladder retained a direct
   coarse-to-medium spectral veto.

PREF11, SRC3, and SRC4 resolved the first obstruction without modifying the
equations, branch, reference geometry, or threshold. RSP1 then asked the
second question on the independent `4097 -> 8193 -> 16385` ladder. All six
members reached `t=1/16` with zero source retries, and both numerical methods
passed the prospectively fixed direct `phi/Lambda` derivative-tail target.
The generic all-field raw spectral aggregate nevertheless remained false.

PROTO12 freezes the smallest justified response. It does not declare those
failed raw ratios to be passes. It freezes an explicit classification that
can distinguish a directly contracting tail from a tail already below the
measured sensitivity of the diagnostic map.

## Unchanged source physics

The GR-0 source backend becomes SRC4's tensor-contracted implementation of the
same SRC3/ACT1/VAR1/REF1 equations. The following remain unchanged:

- the six affine acceleration unknowns and physical root branch;
- the exact spherical reference geometry and PROTO11 derivative map;
- the complete unredefined residual evaluated before any exact-reference
  identity branch;
- the raw `1e-12` source ceiling, `1e10` kinetic-condition ceiling, and maximum
  sixteen complete-refinement iterations;
- every coupling, initial profile, amplitude, grid, method, CFL rule,
  dissipation setting, stop, observable, and claim.

SRC4 changes contraction order in binary64. Its independent exact-dyadic,
SRC3 differential, exact-reference, and CAP1 controls establish numerical
instrument equivalence within prospectively frozen bounds. Local timing is
context only and is not a scientific result.

## The pairwise spectral rule

For each field, each derivative-weighted spectrum, and each adjacent pair on
the unchanged PROTO12 calibration ladder

```text
2049 -> 4097 -> 8193,
```

the direct route remains

```text
finer top-band power / coarser top-band power < 1/4.
```

Every numerator, denominator, fraction, and raw ratio remains serialized. A
ratio that fails this strict test may receive only the label
`diagnostically_saturated`, and only when all of these independent predicates
are true:

- both grids in that pair pass their unchanged individual absolute spectral
  budgets;
- both top-band field amplitudes lie below their own measured round-trip map
  scales;
- both derivative-weighted amplitudes lie below the corresponding Nyquist
  multiple of those map scales;
- erasing each top band changes the windowed diagnostic input by no more than
  that grid's own measured map scale;
- the round-trip residual contracts strictly on all three grids, or is exactly
  zero on all three; and
- the complete resampled profile contracts strictly in both infinity and RMS
  norm across the three-grid sequence, or is exactly zero.

This is not an epsilon floor. No small denominator is replaced, no failed
ratio is removed, no asymptote is fitted, and no tolerance is tuned to the
answer.

## Why RSP1 supports this rule

The reproducer restores all six hash-bound RSP1 terminal states and recomputes
the complete spectral budgets, raw ratios, proper-distance profiles,
round-trip residuals, top-band erasure witnesses, and pairwise classifications
for RK4 and SSPRK3.

Every individual absolute budget passes on every RSP1 grid. Every complete
profile and round-trip sequence contracts. Every raw ratio that fails on an
adjacent pair has both pair members below their independently measured map
scales. The top-band erasure is less than `0.0061` of the measured round-trip
non-idempotence in the worst observed case. Both methods therefore pass the
fully guarded pairwise classification while their generic raw all-field
admissions remain explicitly false.

RSP1's high ladder is qualification evidence for the classifier. It does not
replace PROTO12's original calibration ladder, select a physical amplitude,
or answer the collapse mechanism question.

## The epistemic boundary

The round-trip residual measures non-idempotence of the diagnostic map from a
native radial grid to a uniform proper-distance grid and back. It is not a
continuum discretization error bound, a physical uncertainty, or evidence that
the underlying field is resolved in nature. `Diagnostically_saturated` means
only that this specific FFT-based tail ratio cannot support a stricter claim at
the measured map sensitivity.

The classifier can still reject a future state. It fails if a direct ratio is
not below `1/4` and either pair member loses its absolute budget, map-scale
witness, or complete three-grid profile contraction. Constraint, source,
kinetic, characteristic, boundary, scale, trapped-sign, temporal, and every
other inherited stop remain independent and fail closed.

PROTO12 has not constructed a fresh state or advanced a trajectory. It does
not show that GR-0 will calibrate, that a trapped sphere will form, that the
regulator activates, or that the complete affine-null Raychaudhuri expression
becomes positive. It does not reject FGC-QR or validate it. It does not derive
retained-EFT validity, singularity resolution, a child domain, a dark sector,
or a varying locally measured speed of light.

## Next gate

`FGC-1-HLT10-MON10` must implement and adversarially validate the complete
PROTO12 runtime before any fresh calibration can be authorized. It must bind
SRC4, reconstruct all frozen inputs, prove exact-reference and differential
controls, compose the pairwise classifier without altering raw evidence,
inject every failure route, and confirm that both new output namespaces are
fresh. Only a separate all-of authorization may open one GR-0 calibration.

## Reproduction

With the hash-bound RSP1 terminal checkpoint present:

```bash
python3 scripts/reproduce_fgc_pro12_frz1.py --check
```

The canonical result remains inspectable without treating the external raw
checkpoint as a tracked Git artifact. A missing or hash-mismatched checkpoint
prevents recomputation; it does not weaken the stored nonclaim boundary.
