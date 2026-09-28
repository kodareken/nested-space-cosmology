# FGC-1-CAL4-PREF6: PROTO7 common-event ownership diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration admission diagnosis; no
mechanism result

```text
PROTO7_common_event_ownership_contract_obstructed = true
PROTO8_common_event_admission_revision_required = true
fresh_GR0_dynamic_calibration_completed = false
FGCQR_holdout_execution_authorized = false
```

## What the frozen campaign established

HLT5 authorized exactly one fresh GR-0 calibration under PROTO7. The campaign
was launched from immutable commit
`b8f1921985191b9ad488ecfd2f183c19936996b1` and terminated normally after
about 749 seconds. Both amplitudes passed the complete primary/comparator
common event at `t=0`, crossed the source-solver obstruction found by CAL3,
and completed the first nonzero common event at `t=1/16`.

For each amplitude, the finest RK4 member rejected and exactly rolled back two
wholly unaccepted source-only proposals. The unchanged raw source tolerance
remained `1e-12`; no non-source premise failed; the subsequent proposals
committed. PROTO7 therefore did what HLT5 authorized it to do.

Both amplitudes were then rejected by the shared
`common_event_constraint_or_spatial_spectral_admission` gate. The runner
returned its frozen classification, `calibration_failed_no_eligible_GR0_case`,
selected no amplitude, and read no SGB-L or FGC-QR trajectory. The manifest,
eight-line event log, terminal checkpoint, and result are bound here by
SHA-256 before interpretation.

## Why the constraint veto does not own the evidence it used

The numerical projector intentionally fixes four outer rows to their initial
vacuum values. A spatial derivative near that interface necessarily combines
evolved points with fixed points. The old common-event diagnostic nevertheless
asked the auxiliary reduction constraints `q-D_h u` to converge over the
entire grid, including the projector-owned interface.

At `t=1/16`, the physical and gauge constraints peak near the collapsing pulse
and converge with the expected method orders. For amplitude `5/2`, the
ownership-scoped global norms are:

```text
RK4       1.800673e-2  -> 6.001664e-3  -> 7.044068e-4
SSPRK3    2.805935e-2  -> 1.396234e-2  -> 4.722423e-3
```

The finest-pair orders are approximately `3.09` to `3.72` for RK4 and `1.56`
to `1.95` for SSPRK3, all above the unchanged `3/2` minimum. Both the coarse
and fine magnitude guards pass.

The four old nonconvergent reduction components are instead:

```text
reduction_alpha
reduction_shift
reduction_lambda
reduction_R
```

Every full-domain maximum for those components lies beyond the last row owned
by the evolution diagnostic—inside the fixed-row/stencil interface near
`r=128`. CAL4 excludes exactly the four fixed rows plus the derivative
operator's declared stencil reach: three additional rows for RK4 and one for
SSPRK3. It does not use the much smaller `r<=24` measurement region as an
excuse to hide errors elsewhere.

Inside that owned domain, `q` and `D_h u` agree to accumulated binary64
roundoff. The old fixed `4096*epsilon` enclosure covered one diagnostic
evaluation but not the independently accumulated Runge–Kutta histories of
`q` and `u`. CAL4 uses the deterministic upper budget

```text
4096 * binary64_epsilon * (1 + accepted_stage_count)
```

for convergence-order classification only. Raw residuals remain serialized;
they continue to decide the unchanged magnitude guards and are never
subtracted from a physical or observable margin. Under that enclosure all six
reduction components are indistinguishable from accumulated roundoff, while
the physical and gauge orders remain unchanged and finite.

This is an ownership correction, not a relaxed physical constraint. An
interior reduction defect above the declared operation-count envelope still
fails, and projector-owned outer errors remain visible in the full-domain
record.

## Why the spectrum veto also mixed roles

The old nested-spectrum gate required every individual absolute budget on all
three grids and then separately required every adjacent tail ratio to decay.
That assigns two incompatible jobs to the coarsest grid: it must be coarse
enough to witness convergence, yet already meet the final resolved-grid
budget.

For amplitude `5/2`, the only absolute-budget failures are the coarsest-grid
lapse and shift spectra. The medium and fine grids pass every unchanged
absolute field-power, derivative-power, and RMS-scale threshold for all six
fields. Every coarse-to-medium and medium-to-fine tail ratio also remains
below the unchanged `1/4` ceiling. Representative contractions are:

```text
RK4 alpha field tail       9.59e-3 -> 2.12e-5
RK4 shift derivative tail  3.13e-2 -> 1.61e-4
SSPRK3 alpha field tail     1.05e-2 -> 2.81e-5
SSPRK3 shift derivative     3.06e-2 -> 1.54e-4
```

CAL4 therefore assigns the roles explicitly:

```text
coarsest grid
    convergence witness; raw budget remains public

medium and finest grids
    each must satisfy every unchanged absolute spectral budget

both adjacent refinement pairs
    every field and derivative tail must contract below 1/4
```

No cutoff or threshold is moved. A coarse-only miss can pass only when the
entire nested sequence demonstrates contraction and the finest pair is
independently resolved.

Amplitude `3` proves that this is not an outcome-directed bypass. Its
constraint magnitudes narrowly exceed the frozen RK4 and SSPRK3 fine guards,
and its `phi` spectral tail fails the unchanged nested-ratio ceiling. It
remains inadmissible under the prospective rule.

## Deterministic replay and decision

The terminal checkpoint owns amplitude `3`, so CAL4 reconstructs amplitude
`5/2` from the immutable HLT5 inputs and committed PROTO7 implementation. The
replay independently reproduces:

- both public source-only retry records after canonical serialization;
- the complete public `t=1/16` common-event assessment;
- all six final state hashes;
- the old rejection; and
- the corrected RK4 and SSPRK3 admissions.

The replayed states are cached atomically outside Git so a formatting failure
cannot require another evolution. Their hashes and the reduced canonical
evidence are bound into the tracked certificate.

The result is:

> PROTO7's transaction repair succeeded, but its common-event gate combined
> projector-owned boundary data and coarse-grid convergence witnesses with
> final-resolution admission. Amplitude `5/2` clears the corrected ownership
> contract under every unchanged numerical threshold; amplitude `3` does not.

This requires a separately frozen PROTO8 protocol. PROTO8 may alter only the
common-event ownership rules established here. The physical inputs, amplitude
order, source solver and raw tolerance, methods, grids, CFL and retry rules,
constraint magnitude guards, minimum convergence order, spectral thresholds,
trapped-sign test, boundary monitor, and all physical stops remain unchanged.

## Scientific boundary

CAL4 is a test-contract result, not the tangible mechanism result sought by
FGC-2-SF1; the mechanism question remains unanswered. It establishes no
sustained collapse, trapped sphere, regulator
activation, defocusing interval, singularity resolution, child domain,
dark-sector mechanism, or varying locally measured speed of light. It does
not reject FGC-QR or the general gradient route. It authorizes neither SGB-L
nor FGC-QR and does not promote retained-EFT validity.

The honest next step is to freeze PROTO8, implement the two ownership adapters
without changing the trajectory engine, authorize one fresh GR-0 campaign,
and let that campaign pass or fail under the predeclared rules.

## Reproduction

Fast validation of the immutable campaign, replay caches, and stored result:

```bash
python3 scripts/reproduce_fgc_cal4_pref6.py --check
```

Full amplitude-`5/2` reconstruction and comparison:

```bash
python3 scripts/reproduce_fgc_cal4_pref6.py --check --replay
```

Regenerate only the deterministic replay caches:

```bash
python3 scripts/reproduce_fgc_cal4_pref6.py --replay-only
```
