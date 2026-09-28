# FGC-1-CAL2-PREF4: first dynamic calibration diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration test-contract diagnosis; no
mechanism result

```text
PROTO5_GR0_source_stop_ownership_contract_obstructed = true
PROTO6_premise_revision_required = true
fresh_GR0_dynamic_calibration_completed = false
FGCQR_holdout_execution_authorized = false
```

## What actually happened

HLT3 authorized one fresh, frozen GR-0 calibration under PROTO5. The campaign
ran from immutable commit `fb7d1dc2499438ac639e6a1686cd3f5561fa93bc`
and terminated normally after about 377 seconds. It did not select an
amplitude.

Both declared amplitudes, `5/2` and `3`, passed the complete primary and
comparator common-event admission at `t=0`. Neither completed the first
nonzero common event at `t=1/16`. Both encountered the same branch-independent
source gate while advancing the finest RK4 member:

```text
amplitude 5/2: residual = 1.0511642392237476e-12
amplitude 3:   residual = 1.0511640449698101e-12
frozen limit: residual = 1.0000000000000000e-12
```

The excess is about 5.12 percent, and the two observed residuals differ by
less than two parts in ten million relative to one another. The campaign
therefore returned its frozen terminal label,
`calibration_failed_no_eligible_GR0_case`. That label is preserved exactly in
the raw record. It is not silently converted into a successful calibration.

## Why this is not a collapse or mechanism result

The terminal checkpoint contains no synchronized state after `t=0`. For the
second amplitude, the two coarser RK4 members had reached `t=1/16`, the finest
RK4 member remained near `t=0.03906`, and the SSPRK3 members had not begun.
The runner correctly refused to construct a mixed-time common event.

Consequently the campaign did not establish any of the following:

- a dynamically trapped sphere;
- an eligible GR-0 amplitude;
- an SGB-L control;
- FGC-QR activation or health;
- a Raychaudhuri sign;
- singularity resolution or a physical transition.

The candidate action was not executed. The result cannot support or reject
FGC-QR, Finite Gradient Closure, or a general gradient mechanism.

## The stop-ownership replay

CAL2-PREF4 freezes the manifest, event log, terminal checkpoint, and terminal
result by SHA-256 before interpreting them. It then restores only the last
accepted `RK4-4097` state of amplitude `3` and replays four *proposals* from
that state. It does not accept any of them or alter the original campaign.

The replay establishes four facts:

1. The source solve at the last accepted state itself remains below the
   original `1e-12` absolute residual limit.
2. Later internal stages of the full trial step miss that raw limit.
3. Repeating the proposal at factors `1/2`, `1/4`, and `1/8` shows that even
   the first half-step retry, `1/2`, passes the unchanged raw source gate at
   every stage.
4. Throughout the replay, the kinetic condition number remains far below the
   frozen `1e10` stop, and the unresolved affine acceleration correction stays
   below `1e-12` relative to the acceleration scale. The largest failed raw
   residual is localized at the first annular node, where spherical-coordinate
   cancellation is most severe.

This distinction matters. An accepted-state source failure cannot be cured by
changing a Runge--Kutta proposal and remains a legitimate terminal premise
failure. An internal-stage failure that disappears when an *unaccepted* trial
step is reduced is a numerical proposal failure. Treating the latter as
evidence that the continuum system or physical mechanism failed assigns a
physical meaning to a state the integrator never accepted.

## Minimal prospective repair

The diagnosis permits one narrow successor and nothing broader:

- preserve the action, initial data, amplitude order, physical parameters,
  raw `1e-12` source residual limit, condition limit, constraints, spectra,
  boundary accounting, and trapped-sign rule;
- evaluate the source gate at the last accepted state before each proposal;
- keep an accepted-state source failure terminal;
- make only an internal, unaccepted-stage source miss retryable using the
  already frozen factor `1/2`, maximum of 32 retries, and minimum step size;
- retain every raw residual and retry count in the public record;
- require a new protocol version and a fresh output namespace.

This is not permission to edit PROTO5 or rerun its namespace. A separately
frozen `FGC-2-SF1-PROTO6` and a machine-tested HLT4 implementation are required
before another calibration may begin.

## Scientific value of the failure

The first dynamic campaign did useful work: it found that the test confused
an adaptive integrator's rejected trial with an accepted continuum state. The
failure closed that ambiguity before the candidate action or a physical
observable could be exposed to it. It therefore strengthens the eventual
falsification test, but it does not answer the mechanism question or yet
advance the physical evidence for or against the gradient proposal.

## Reproduction

With the immutable raw campaign present at its frozen local path:

```bash
python3 scripts/reproduce_fgc_cal2_pref4.py \
  --config configs/fgc/fgc-1-cal2-pref4.toml \
  --output results/fgc-1-cal2-pref4.json
python3 -m unittest tests.test_fgc_cal2_pref4_reproduction -v
```

The tracked result stores the compact diagnosis and hashes every raw campaign
source. Raw checkpoints remain outside Git under the repository's declared
large-run-data boundary.
