# FGC-1-HLT12-MON12 — PROTO14 TDG6 runtime authorization

`FGC-1-HLT12-MON12` is a narrow pre-trajectory authorization for the fresh
`FGC-1-CAL11-RUN1-PLAN` PROTO14 GR-0 calibration. It binds the production TDG6 temporal-admission
runtime to the six immutable `t=23/16` restart payloads. It is not a collapse,
trapped-region, candidate-action, defocusing, EFT, or physical result.

## Question owned by this gate

The question is operational: can the frozen PROTO14 inputs enter the one
successor runner which will require TDG6's one/two/four complete-state test at
every accepted post-restart macro step? The result must prove this without
creating either PROTO14 namespace or advancing any member.

The restart state is discrete calibration input, not evidence that a
pre-restart temporal error bound exists. All six payloads retain their
original state, tracers, 24-event history, monitor, causal debit, accepted
stage count, source-retry count, and CFL-retry count. TDG6 begins at the
restart with an exact-zero 18-channel debit and zero temporal retries; that
new ledger does not reinterpret the inherited history.

## What is bound

HLT12 directly verifies the immutable PRO14 freeze at commit
`11c800cbd6f31d10d2b9578ed150c23a0a19aa75`, its protocol and TDG6 runtime
module, the raw PROTO12/RSP2 checkpoints, and all six complete payload hashes.
It checks that payload restoration leaves the raw checkpoint bytes unchanged.

It additionally records the exact arithmetic fingerprint—CPython `3.14.3`,
NumPy `2.5.1`, Accelerate, Darwin `arm64`—before source evaluation or state
advance. The canonical authorization hashes the TDG6 runtime, the successor
PROTO14 runtime adapter, and the dedicated PROTO14 runner. A changed backend
or implementation cannot silently continue the restart.

The only new true gates are:

- `PROTO14_successor_runtime_implemented`;
- `PROTO14_fresh_GR0_dynamic_calibration_authorized`.

Every candidate, EFT, transition, cosmological, and physical gate remains
false. In particular, this authorization does not say that GR-0 produces a
trapped interval or that FGC-QR has been read.

## Temporal and namespace boundary

The successor runtime must require TDG6 admission for every accepted macro
step, commit only the admitted four-quarter fine path, retain source and CFL
retry ownership, serialize its 18-channel debit and complete rejections in
checkpoints, and keep the old 64-sample history only as a non-veto diagnostic.

At the moment the authorization is generated, both
`runs/fgc-2-sf1/proto14/calibration` and
`runs/fgc-2-sf1/proto14/holdout` must be absent. A clean clone where both are
absent is allowed; partial presence or namespace reuse fails. Once a later
authorized run creates calibration evidence, verification follows the
immutable authorization record rather than re-asking that past absence
question.

Terminal publication is a two-phase recoverable transaction. The atomic
checkpoint first owns the complete canonical prospective result and its hash;
only then may the result file be published atomically. A restart either
publishes that exact checkpoint-owned result, verifies an already matching
publication, or fails closed on any mismatch. A result beside a nonterminal
checkpoint is forbidden. Local authorization generation requires both raw
restart archives and verifies their exact hashes. A portable clean-clone audit
may instead verify only the compact tracked certificate; in that state the raw
restart evidence is explicitly unverified, never reported as locally replayed.

## Failure semantics

- A checkpoint or payload-hash mismatch is provenance failure.
- A backend mismatch fails before source evaluation or state advance.
- Missing or hash-drifted TDG6/runtime/runner code prevents authorization.
- A runtime crash later belongs to the calibration execution and is invalid,
  never a physical negative.
- This gate opens neither SGB-L nor FGC-QR.

The successor runtime and runner are now settled and bound by the canonical
authorization record. Reproduce and verify that no-trajectory result with:

```bash
python3 scripts/reproduce_fgc_hlt12_mon12.py \
  --check --output results/fgc-1-hlt12-mon12.json
python3 -m unittest \
  tests.test_fgc_proto14_runtime \
  tests.test_fgc_hlt12_mon12_reproduction \
  tests.test_fgc_gr0_campaign_runner_v14 -v
```

The result is also consumed by `make verify-fgc-pro14-prelaunch`. The long
calibration trajectory remains a separate recoverable command which is
forbidden until this authorization is sealed in an immutable clean tracked
commit.
