# FGC-1-CAL11-PREF22: PROTO14 pre-shadow terminal binder

`FGC-1-CAL11-PREF22` binds the completed terminal transaction of the fresh
PROTO14 amplitude-`3` GR-0 calibration.  It is an outcome-neutral
implementation/nonconvergence result.  It is not a completed GR-0
calibration, a temporal-admission result, a constraint or trapped-surface
result, a candidate result, or a physical conclusion.

## Frozen terminal identity

The raw bundle is bound by the authorization commit
`7534a1662d8078a63c98025ecd464bb068fa012f` and campaign identifier
`038c0128e76f48118b97d985384a1e7239edf7c97349cfc7552ee7d33e289e98`.
The manifest, two-line canonical event log, atomic checkpoint, and external
result have SHA-256 values recorded in the CAL11 configuration.

The checkpoint retains the inherited event-`23` (`t=23/16`) restart boundary.
The first proposed post-restart target, event `24`, aborts while `RK4-2049`
is active with:

```text
ValueError: TDG6 subdivision is not bitwise uniform
```

The abort is classified as `invalid_implementation_or_nonconverged_run`.
The event log explicitly marks its tail as recovery-orphan evidence rather
than a physical classification, and the checkpoint remains the owner of the
accepted state.  No post-restart macro step was accepted for any member; every
TDG6 ledger therefore has zero accepted steps, zero temporal retries, and a
zero eighteen-channel debit.  The external terminal result equals the result
embedded in the terminal checkpoint, and that checkpoint embeds the exact
external event-log bytes.

## Two replay partitions

The portable reconstruction reads the four raw files, checks all hashes,
restores the checkpoint arrays and TDG6 extensions, and proves the event,
checkpoint, and result identities.  It never runs or resumes a campaign.

The source-pinned partition is deliberately smaller: it reads the exact TDG6
runtime blob from the immutable authorization commit, then verifies the
hash-bound `_subdivision_boundaries` guard and its recorded exception.  It does
not depend on the mutable successor runtime path, so a later prospective TDG7
repair cannot rewrite this historical result.  The raw abort does not
serialize the proposed macro width, so a numerical invocation would require
reconstructing and advancing mutable runner state.  CAL11 therefore does not
pretend to have derived the arithmetic root cause.  It establishes only a
deterministic subdivision-validator fault at the pre-shadow boundary.

## Boundary and successor

The historical checkpoint may not resume.  CAL11 does not repair the runtime,
mutate PROTO14, authorize a fresh calibration, or open SGB-L, FGC-QR, DEF1,
EFT, transition, or physical claims.

Only `FGC-1-TDG7-FRZ1`, a prospective no-trajectory subdivision-uniformity
diagnosis design, may begin.  A later repair must be independently qualified
and bound in a new protocol/runtime authorization and output namespace before
any fresh GR-0 calibration is considered.
