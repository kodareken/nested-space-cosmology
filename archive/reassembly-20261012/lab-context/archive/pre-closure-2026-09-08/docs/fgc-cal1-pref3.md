# FGC-1-CAL1-PREF3: pre-trajectory semidiscrete composition audit

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced test-contract diagnosis; no trajectory

```text
PROTO4_semidiscrete_common_event_contract_obstructed = true
PROTO5_premise_revision_required = true
PROTO4_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

## Question

ID2 leaves exactly two statically eligible GR-0 amplitudes, (5/2) and (3).
Before opening either trajectory, this artifact performs the composition that
HLT2 and ID2 did not: it maps each complete initial state into each declared
method-of-lines discretization and evaluates all ten PROTO4 constraints on the
three nested grids at the common event (t=0).

This is an input and numerical-premise audit. It reads no fresh calibration,
SGB-L, or FGC-QR evolution output.

## What the composition found

ID2 stores continuum-compatible analytic first derivatives in the auxiliary
field (q). A first-order semidiscrete evolution instead monitors

\[
 C_q=q-D_hu.
\]

Those are not the same finite-grid object. On the original ID2 arrays the
physical and gauge projections are small, but (C_q) contains ordinary
derivative truncation error. For amplitude (5/2), the primary normalized
global values on (1025,2049,4097) points are approximately

\[
  6.19\times10^{-2},\quad
  1.17\times10^{-2},\quad
  1.26\times10^{-3}.
\]

They miss PROTO4's raw coarse and fine guards. This does not mean the
continuum constraints fail; it means the continuum-to-grid map was never
frozen.

Projecting only the auxiliary variable,

\[
 q\leftarrow D_hu,
\]

preserves (u) and (p) bitwise and makes the discrete reduction constraint
zero initially. The remaining physical/gauge residuals then converge. The
fourth-order primary fits the existing magnitude guards. The second-order
comparator has the expected larger truncation error and misses the single
method-independent fine guard even while its finest pair converges.

A separate issue appears at the floating-point floor. Constraints that vanish
analytically can evaluate near (10^{-17}). PROTO4's exact-zero-only special
case assigns these values a finite—and sometimes negative—convergence order.
That turns roundoff into a fake continuum failure.

Therefore the exact PROTO4 common-event rule cannot admit its own declared
semidiscrete experiment at (t=0). No trajectory should be launched under
that contract.

## Prospective repair demonstrated without outcomes

The companion adapter demonstrates, but does not yet authorize, the smallest
repair:

- preserve physical (u,p) and project only (q=D_hu) for the native SBP
  method before the first stage;
- retain raw normalized constraints for all magnitude guards;
- attach a public (4096\epsilon_{64}) normalized floating-operation
  enclosure for zero classification only;
- require monotone refinement across all three grids and order at least
  (3/2) on the finest adjacent pair;
- retain the primary guards (1/50) and (1/1000);
- use prospective lower-order comparator guards (1/10) and (1/200), while
  keeping the same finest-pair order requirement.

Both eligible amplitudes pass that prospective (t=0) composition for both
methods. The looser comparator magnitude guard does not hide its raw values:
all raw component norms and the enclosure remain serialized. The convergence
requirement, independent method, and later cross-method trapped-sign test are
unchanged.

Because PROTO4's revision policy says a premise change requires a new protocol
version, this preview cannot silently repair PROTO4. It justifies PROTO5. The
new version must freeze these rules and a new output namespace before any
fresh trajectory is opened.

## Scientific boundary

This is a false test design caught before it can taint the mechanism test. It
says nothing about whether GR collapse forms a trapped sphere, whether FGC-QR
activates, or whether metric-null defocusing occurs. It neither supports nor
rejects Finite Gradient Closure, black-hole continuation, nested domains, or
any dark-sector interpretation.

The useful result is narrower: continuum constraint compatibility is not by
itself a semidiscrete Cauchy datum, and convergence tests need a declared
floating-point zero enclosure if exact analytic zeros are represented in
binary64. The test is being repaired before its answer is observed.

## Reproduction

```bash
python3 scripts/reproduce_fgc_cal1_pref3.py \
  --config configs/fgc/fgc-1-cal1-pref3.toml \
  --output results/fgc-1-cal1-pref3.json
python3 -m unittest tests.test_fgc_cal1_pref3_reproduction -v
```

The canonical machine record is
[`results/fgc-1-cal1-pref3.json`](../results/fgc-1-cal1-pref3.json).
