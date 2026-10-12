# FGC-1-SRC2-CAP1: evolved affine-source replay and arithmetic-capture contract

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** prospectively frozen diagnostic runner and arithmetic comparison;
no SRC2 run result, successor protocol, candidate trajectory, or physical
claim

```text
SRC2_capture_implementation_frozen = true
PROTO11_source_wall_replayed = false
PROTO11_test_instrument_shown_incomplete = false
PROTO12_frozen = false
FGCQR_holdout_execution_authorized = false
```

## Why this gate exists

[CAL8/PREF10](fgc-cal8-pref10.md) localized the amplitude-`5/2`
`RK4-8193` stop to a complete affine-source residual of
`1.0659145473163184e-12` against the unchanged strict `1e-12` ceiling. The
failure survives exact proposal rollback and binary step halving, but the
immutable campaign did not retain the six-by-six pointwise affine matrix or
constant term at the terminal candidate endpoint. CAL8 therefore cannot tell
whether the wall is caused by conditioning, coefficient-extraction
cancellation, binary64 representability, complete-residual evaluation, or a
deeper semidiscrete incompatibility.

SRC2-CAP1 freezes the smallest experiment that can distinguish those cases. It
does not assume that the old solver is wrong, and it does not treat failure of
one arithmetic route as failure of FGC-QR or the wider relational-gradient
programme.

## Immutable replay target

The runner consumes the original PROTO11 manifest, event ledger, checkpoint,
campaign result, HLT9 authorization, CAL8 plan, and CAL8 result by exact hash.
It rebuilds only the implicated GR-0 member:

```text
amplitude = 5/2
method = RK4
point_count = 8193
first evolved event = t=1/16
terminal accepted time = 0.0706352245695038
terminal accepted step/serial = 97/485
source-only rejection records = 164
terminal evaluated stage time = 0.07063522612170793
```

The first evolved state hash, all 164 canonical rejection records, the final
accepted-state hash, and the typed minimum-step exhaustion must reproduce
exactly. A replay that merely reaches a similar residual is rejected.

## Captured object

At the final candidate endpoint, the complete GR-0 REF1 residual is affine in
the six ADM accelerations at each positive-radius node:

\[
   \mathcal R_i(a_i)=b_i+A_i a_i,
   \qquad A_i\in\mathbb R^{6\times6}.
\]

The ignored raw fixture binds the complete evaluated state
`(u,p,q,p_r,q_r,r)`, the forward-extracted `A` and `b`, the baseline
acceleration, the complete baseline residual, and the best diagnostic
candidate by canonical array-content hashes. The fixture also binds the exact
source authority and PROTO11 reference-derivative implementation through the
run manifest.

This is a binary64 representation of the semidiscrete equations. Exact
rational work on its coefficients diagnoses finite-precision linear algebra;
it is not an assertion that the stored coefficients are exact continuum
numbers.

## Frozen arithmetic attacks

All routes must return a binary64 acceleration to the same complete
unredefined REF1 evaluator. Only that evaluator's strict residual
`<1e-12` is admissive.

1. **PROTO11 baseline.** Reproduce the zero-plus-unit coefficient extraction,
   NumPy solve, and existing complete-residual refinement.
2. **Power-of-two equilibration.** Scale rows and columns only by binary powers
   before the solve. Those scalings are exact while they remain in normal
   binary range, but improvement is measured rather than assumed.
3. **Symmetric affine extraction.** Compute each column from
   `[R(+s e_j)-R(-s e_j)]/(2s)` at prospectively frozen power-of-two scales
   `s in {1,16,256}`. This attacks cancellation in coefficient extraction; it
   does not remove physical or coordinate asymmetry.
4. **Exact-binary residual refinement.** At no more than eight nodes whose
   complete residual is at least half the unchanged ceiling, solve the captured
   binary64 correction equations as exact rationals and round the correction
   back to binary64. The complete REF1 evaluator still owns acceptance.

The long-double reconstructed affine residual is diagnostic only. It can show
that a captured linear system was solved accurately, but it can never replace
the complete REF1 residual or its threshold.

## Nontrivial controls

Before the expensive replay, unit controls establish that:

- the captured forward rows reconstruct fresh complete residual evaluations
  to a scale-aware binary64 roundoff bound;
- every stable route preserves a nontrivial amplitude-`5/2` PROTO11 initial
  slice under the unchanged gate;
- selected exact-binary solves do not alter unselected nodes;
- non-power-of-two probes, threshold changes, broadened trajectory scope, and
  promoted claims fail closed; and
- target capture requires bitwise time and residual identity.

Passing these tests establishes instrument competence, not that the evolved
wall is arithmetic.

## Decision boundary

If a frozen arithmetic route clears the terminal wall under the unchanged
complete evaluator and the later certificate verifies all replay, control,
hash, and nonclaim conditions, SRC2 may conclude that the original PROTO11
solve path produced a false numerical rejection. That would justify designing
a separately frozen PROTO12; it would not by itself complete GR-0 calibration
or authorize FGC-QR.

If no route clears the wall, SRC2 must report the measured backward error,
conditioning, coefficient sensitivity, representability, and complete-
residual discrepancy. That is a structural result for this semidiscrete
source representation under the declared arithmetic routes. It is not, by
itself, a general obstruction to FGC-QR, thermodynamic gradients, nested
domains, a particle-as-spectral-pole interpretation, or any physical model.
A broader negative claim requires the converged and explicitly scoped
obstruction burden defined by FGC-2-SF1.

The amplitude-`3` direct coarse `phi/Lambda` derivative-tail veto remains
binding throughout. SRC2 cannot erase it, move the grid, change a threshold,
or inspect SGB-L/FGC-QR outcomes.

## Commands

The cheap prospective contract check creates no output:

```bash
python3 scripts/run_fgc_src2_affine_capture.py --check
```

After the implementation is committed with a clean tracked worktree, the raw
capture is launched explicitly with:

```bash
python3 scripts/run_fgc_src2_affine_capture.py
```

Raw outputs live only beneath the ignored
`runs/fgc-2-sf1/src2/affine-wall-capture/` namespace. They require a separate
hash-bound PREF11 reduction before any public result or successor freeze.
