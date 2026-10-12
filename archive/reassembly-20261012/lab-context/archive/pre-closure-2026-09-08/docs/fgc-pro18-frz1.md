# FGC-1-PRO18-FRZ1: production-prelaunch overlay

**Status:** static, compact-only production-prelaunch freeze. This document
does not open a source archive, replay a source event log, create either future
PROTO17 output root, or advance a numerical state.

## What is frozen

PRO18 remains an overlay on `FGC-2-SF1-PROTO17`, not a new physical protocol.
It binds the immutable HLT14 checkpoint `4b27a48facf75b1c3186efed6cdfa752ce160e93`,
the compact HLT14/PROTO17/PROTO15/CAL9/RSP2 artifacts, exactly two shared
source-container names, their two legacy `events.jsonl` names, and the six
restart selectors at `t=23/16`.  Each selector has a fixed source container,
restart and input identity, step index, and transaction serial.

The future launch authority has a distinct identity, `FGC-1-PRO18-AUTH1`.
It must commit an authority blob, build an isolated detached launch image, and
in the same launch invocation recheck actual authority and imports immediately
before atomic genesis. This reduces, but does not eliminate, time-of-check/
time-of-use exposure. The external Git trust root is not cryptographically
authenticated by this repository, and this freeze does not validate a launch.

## Deliberate access boundary

`--write` and explicit `--verify-prelaunch` authenticate compact tracked data
through both `git show` at the sealed commit and live worktree bytes, then may
inspect only whether `runs/fgc-2-sf1/proto17/calibration` and `.../holdout`
exist. Ordinary `--verify` validates the recorded freeze-time absence plus immutable
compact lineage without observing those roots again. The commands never
open, hash, parse, list, or semantically replay either named input NPZ
archive or event log. The absent-root observation is temporal evidence at this
freeze boundary, not a fact to be reinterpreted after a future launch.

## Required future chain

```text
FGC-1-PRO18-PREF26  real raw archive and legacy-history verification, read-only
  -> FGC-1-PRO18-AUTH1  committed external authority and detached launch image
  -> FGC-1-HLT15-GEN1  authority recheck plus atomic generation zero
  -> FGC-1-PRO18-PREF27  actual namespace reuse and foreign-store binding
```

Every stage has to fail closed with its typed stop. No stage may convert an
absent namespace observation into authority, or use a later output manifest as
the authority that permitted it.

## Nonclaims

The freeze earns only five static statements: the PROTO17 identity and compact
bindings are fixed, HLT14's synthetic result is preserved, the production
source/output boundary and PRO18 overlay are frozen, and PREF26 implementation
or design may begin. Every real raw/history/namespace/trajectory/candidate/
physical claim remains false. In particular, this neither authorizes a
pretrajectory operation, fresh GR-0 calibration, SGB-L, FGC-QR, DEF1,
retained-EFT promotion, singularity resolution, child topology, dark-sector
mechanism, nor a variable local speed of light claim.
