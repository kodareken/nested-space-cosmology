# FGC-1-PRO10-FRZ1: outcome-neutral PROTO10 freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** machine-reproduced protocol and immutable-lineage freeze; no
runtime, trajectory, calibration, or mechanism result

```text
PROTO10_outcome_neutral_protocol_frozen = true
PROTO10_guarded_finest_pair_contract_frozen = true
PROTO10_successor_runtime_implemented = false
FGCQR_holdout_execution_authorized = false
```

## Frozen change

PROTO10 consumes FGC-2-SF1-PROTO9 and FGC-1-CAL6-PREF8 directly from
immutable commit `dd3f5adb22f08aa7950b74cef2595aaca2e8cd6e`. It retains the
complete resolution ladder:

```text
2049 -> 4097 -> 8193
```

It changes one premise-level interpretation. The coarse-to-medium field and
derivative tail ratios remain direct vetoes under the unchanged `< 1/4` rule.
The medium and fine grids must still pass every unchanged absolute spectral
budget. On the finest pair only, the direct `< 1/4` result remains the primary
route; a direct failure may instead be classified as diagnostically saturated
only when every independent conditioning and whole-profile guard passes.

This is not a relaxed ratio. The failed raw ratio remains serialized and
visible. No epsilon floor, fitted factor, new tolerance, or changed absolute
ceiling is introduced.

## Guarded finest-pair route

For every failed finest-pair field or derivative ratio, PROTO10 requires:

1. direct coarse-to-medium field and derivative ratios below `1/4`;
2. unchanged medium and fine absolute spectral-budget passes;
3. an orthogonal top-band erasure perturbation no larger than each grid's
   measured proper-grid round-trip non-idempotence;
4. strict contraction of every nonzero round-trip residual on both adjacent
   refinements, with exact zero handled separately;
5. strict contraction of the complete resampled profile difference from the
   first pair to the second in both infinity and RMS norms; and
6. serialization of all raw fractions, ratios, residuals, profile norms, and
   erasure witnesses.

The round-trip residual is a diagnostic-map conditioning measurement, not a
bound on the unknown continuum solution. The top-band erasure is an orthogonal
FFT diagnostic-input witness, not a physical perturbation of the field.
Passing this route means only that the ratio-only discriminator has saturated
below its demonstrated numerical conditioning scale. It does not prove
physical resolution or convergence of an evolved spacetime.

## Exact inheritance

The following remain byte- or contract-inherited from PROTO9 and its
predecessors:

- the ACT1 action, branches, couplings, profiles, physical initial data, and
  calibration amplitudes;
- the `2049 -> 4097 -> 8193` ladder and all constraint ownership, magnitude,
  and order rules;
- RK4 and SSPRK3, spatial orders, CFL maximum, dissipation, final time, retry
  factor, retry budgets, minimum step, and event interval;
- the source solver, raw `1e-12` gate, accepted/unaccepted proposal ownership,
  non-source veto, and rollback evidence;
- every absolute spectral ceiling and the direct nested-tail `1/4` ceiling;
- causal boundary, branch, kinetic, characteristic, constraint, invariant,
  scale, affine-null, outcome, and robustness rules; and
- every candidate, retained-EFT, physical, global, and cosmological nonclaim.

The only other syntactic change is provenance. Future calibration and holdout
outputs must use fresh `proto10` namespaces, and no PROTO9 state or diagnosis
may be relabelled as a PROTO10 trajectory outcome.

## Immutable evidence and fail-closed boundary

PRO10-FRZ1 validates the exact PROTO9 protocol and canonical CAL6/PREF8
diagnosis from the checkpoint commit. It verifies their declared SHA-256
hashes, the checkpoint ancestry, CAL6's pre-trajectory boundary, all sixteen
preserved raw failures, all sixteen guarded classifications, all thirty-two
direct finest-pair passes, and the four prospective design-only admissions.
Both new output roots are absent when frozen.

Mutation controls reject a changed resolution, altered ratio ceiling, weakened
coarse veto, omitted conditioning guard, promoted round-trip semantics, reused
namespace, or promoted claim.

The freeze performs no evolution and does not implement HLT8. It does not
complete GR-0 calibration, select an amplitude, resolve a holdout manifest,
read SGB-L or FGC-QR, classify a trapped sphere, test regulator activation,
answer the mechanism question, validate the retained EFT, reject gradients,
resolve a singularity, derive a child domain or dark sector, or vary locally
measured `c`.

## Reproduction

```bash
python3 scripts/reproduce_fgc_pro10_frz1.py --check
```
