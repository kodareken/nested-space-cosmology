# FGC-1-CAL5-PREF7: PROTO8 resolution-ladder diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration numerical-premise diagnosis;
no mechanism result

```text
PROTO8_campaign_terminated_normally = true
PROTO8_GR0_case_eligible = false
PROTO9_resolution_ladder_revision_required = true
FGCQR_holdout_execution_authorized = false
```

## What the fresh campaign established

HLT6 authorized one fresh GR-0 calibration under PROTO8. The campaign ran
from immutable commit `d667c7465786883b1e5d6afad168b64e7b8921d0`,
terminated normally after about 1,134 seconds, and produced a nine-record
canonical event log. Both amplitudes passed at `t=0`. Amplitude `5/2` also
passed both RK4 and SSPRK3 at `t=1/16`, then both methods stopped at `t=1/8`.
Amplitude `3` stopped in both methods at `t=1/16`.

Each amplitude recorded two rejected RK4-4097 source-only proposals. Every
proposal was wholly unaccepted, its fields and transaction state were
bitwise preserved, and no non-source failure was relabelled as retryable. The
unchanged raw source tolerance remained `1e-12`.

The terminal result is therefore exactly:

```text
classification = calibration_failed_no_eligible_GR0_case
selected_amplitude = null
SGBL_outcome_read = false
FGCQR_outcome_read = false
mechanism_question_answered = false
```

## What PROTO8 repaired—and what it did not

PROTO8's ownership correction is exercised by live evolved data. For
amplitude `5/2`, the complete spatial-spectrum admission passes at `t=1/8`
for both methods. Coarse-grid absolute misses remain public, but the medium
and fine grids satisfy every unchanged absolute budget and both adjacent
tails contract below `1/4`. The PROTO7 spectrum veto was therefore a real
diagnostic-contract defect, and PROTO8 removes it without weakening a
threshold.

The remaining amplitude-`5/2` stop is the physical-constraint calibration
gate. Its owned global residuals are:

```text
RK4       2.963370e-2 -> 9.325011e-3 -> 1.117695e-3
SSPRK3    3.852645e-2 -> 1.864481e-2 -> 6.845100e-3
```

The RK4 finest-pair component order is at least `2.846`; its coarse and fine
norms miss the frozen guards `1/50` and `1/1000` by factors of about `1.48`
and `1.12`. SSPRK3 clears its coarse `1/10` guard but misses its fine `1/200`
guard by about `1.37`; its minimum component order is `1.446`, below the
unchanged `3/2` threshold.

Amplitude `3` is more demanding. At `t=1/16`, its residuals are:

```text
RK4       2.588195e-2 -> 8.634125e-3 -> 1.013440e-3
SSPRK3    4.042782e-2 -> 2.011041e-2 -> 6.799815e-3
```

Both methods have passing finest-pair component orders (`2.513` and `1.555`),
but each misses its fine magnitude guard. In addition, all medium/fine
absolute spectrum budgets pass while the `phi/Lambda` derivative tail fails
the unchanged `1/4` nested-ratio ceiling on the `2049 -> 4097` pair:

```text
RK4       0.326717
SSPRK3    0.266757
```

That failure remains a veto. It is not hidden by calling the fine grids
"resolved."

## Why the next test is a ladder shift, not a threshold change

Every owned constraint norm decreases monotonically with refinement. More
specifically, the existing `2049` states already clear the unchanged
coarsest constraint guards for both methods and both amplitudes. The existing
`4097` states already clear every absolute spectral budget required of a
resolved grid. The smallest experiment that can distinguish ordinary
under-resolution from a formulation obstruction is therefore:

```text
old ladder: 1025 -> 2049 -> 4097
new ladder: 2049 -> 4097 -> 8193
```

The observed finest-pair orders conditionally project the `8193` global
constraint norms below the frozen fine guards. Those projections are recorded
as planning evidence only. They are not observations, do not repair the
SSPRK3 order miss, and do not establish whether the new `4097 -> 8193`
`phi/Lambda` tail contracts. Only new data can answer those questions.

This is why CAL5 requires a separately frozen PROTO9 rather than editing the
finished PROTO8 campaign. PROTO9 may change only the three grid sizes. It must
retain the action, initial physical data, amplitude order, methods, equations,
source solver, raw `1e-12` tolerance, CFL and retry rules, constraint
ownership and guards, spectral budgets and tail ceiling, trapped-sign test,
boundary monitor, health stops, and outcome language. A new namespace is
mandatory.

## Scientific boundary

CAL5/PREF7 is a numerical-premise result. It shows that PROTO8 corrected its
declared ownership defect and that the remaining GR-0 failure is compatible
with an under-resolved convergence ladder. It does not prove under-resolution;
a new `8193` result can still fail. It does not select an
amplitude, complete calibration, observe a trapped sphere, read SGB-L or
FGC-QR, test regulator activation, reject FGC-QR or gradients, validate the
retained EFT, resolve a singularity, derive a child domain or dark sector, or
vary a locally measured speed of light.

A failure under PROTO9 would be more informative than this one. If the new
finest residual does not clear the unchanged guards, the new tail does not
contract, or the order remains inadmissible, the evidence would point away
from mere grid placement and toward the continuum formulation, constraint
propagation, or numerical method itself.

## Reproduction

```bash
python3 scripts/reproduce_fgc_cal5_pref7.py --check
```
