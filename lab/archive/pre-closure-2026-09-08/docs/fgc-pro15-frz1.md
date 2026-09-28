# FGC-1-PRO15-FRZ1: campaign-owned retry lifecycle freeze

**Status:** prospective, machine-reproduced protocol freeze. It performs no
numerical evolution, reads no raw historical campaign data, and creates no
run namespace.

This artifact freezes **FGC-2-SF1-PROTO15** and nothing beyond its prospective
lifecycle contract.

```text
PROTO15_frozen = true
PROTO15_campaign_owned_cursor_contract_frozen = true
PROTO15_successor_runtime_implemented = false
PROTO15_fresh_GR0_dynamic_calibration_authorized = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

## Why another freeze was necessary

CAL11 bound PROTO14's terminal result as an **invalid runtime or nonconverged
outcome**, not a GR-0 or FGC-QR physical obstruction. TDG7 then repaired the
shared stage lattice and qualified that repair synthetically. It also exposed a
necessary limitation: a stateless adapter cannot tell a genuinely fresh call
from a caller that discarded an externally recorded retry receipt and presents
otherwise identical inputs again.

PROTO15 resolves only that lifecycle gap. It defines a campaign-owned durable
cursor; it does not change the action, source evaluator, equations, matter
input, reference map, CFL policy, constraints, health/scale/boundary rules,
null observable, calibration classification, or physical claim boundary.

## Frozen cursor contract

The successor runtime must maintain exactly one durable cursor per member in
one of two modes:

```text
FRESH_READY
RETRY_PENDING
```

The cursor binds the protocol and campaign identities, member/method/grid,
committed common-event index,
accepted state hash, rational time **and lowercase binary64 hexadecimal time**,
the complete TDG6-ledger hash, accepted step/transaction serial, event target,
monotone attempt serial and canonical attempt identity, generation, journal
tip, cursor-chain parent/hash, and either a
JSON-null retry payload or the complete exact retry-successor payload. All
durable records use canonical compact sorted UTF-8 JSON. Every journal record
has the envelope `record_kind`, `journal_parent_sha256`, `payload`, and
`record_sha256`; its hash is computed from the canonical record before its
self-hash field. A cursor's `journal_tip_sha256` is the latest durable journal
record immediately preceding that cursor transition. The cursor-transition
record's own hash becomes the checkpoint's new journal tip. The cursor hash is
likewise computed before its self-hash field. The zero hash is the only
journal- or cursor-chain root parent.

The complete ledger hash binds last accepted time in float-hex, accepted-step
count, cumulative and current retry counts, last accepted macro-step retry
count, the complete debit vector in float-hex, and serialized rejections. The
freeze exposes a `cursor_origin_template` for each historical restart only. It
is deliberately **not** a serialized complete cursor or an executed journal
record.

`FRESH_READY` requires both a JSON-null retry payload and TDG6 current retry
count zero. `RETRY_PENDING` requires a positive retry count and the exact
predecessor plan, successor half-plan, durable TDG6 rejection evidence,
complete serialized rejection prefix, member identity, state/time/step/serial,
half-cap, and frozen minimum step. The successor additionally binds the durable
rejection-record hash and predecessor journal tip. Its retry count is exactly
the previous current count plus one; its cap is the exact binary64 half of the
predecessor cap; accepted identity is unchanged; and rejection evidence is the
prefix tail and the payload of the hashed rejection record. A pending cursor
cannot be reset to fresh
at the same accepted boundary. Only an accepted fine commit or an explicitly
recorded common-event rollback may transition it.

## The two failure paths after a TDG6 rejection

The ordering is part of the scientific instrument, rather than an informal
implementation detail.

`TDG6TemporalRetryRequired` and `TDG6TemporalRetryExhausted` remain typed
exceptions, but the campaign owner must catch them before they escape. That
same owner must complete the exact pending or terminal durable closure before
the outcome may propagate to an outer control loop.

For a retryable TDG6 rejection:

```text
append + fsync durable TDG6 rejection
    -> raise typed TDG6 retry-required
    -> derive exact successor
    -> append + fsync RETRY_PENDING cursor transition
    -> adopt that cursor
    -> invoke the retry entrypoint
```

For an exhausted rejection—maximum temporal retries reached or the frozen
minimum macro step reached—the ordering is instead:

```text
append + fsync durable TDG6 rejection
    -> raise typed TDG6 exhaustion
    -> append + fsync terminal invalid/nonconverged record and checkpoint
```

The terminal exhaustion record binds the campaign/member/method/grid/state,
accepted boundary, step/serial, TDG6 ledger hash, durable rejection-record
hash, predecessor cursor, journal tip, exhaustion reason, and exact
`invalid_implementation_or_nonconverged_run` classification. Its terminal lock
forbids both entrypoints until a separately authorized new campaign identity.
Exhaustion creates no `RETRY_PENDING` cursor and is never a
physical, constraint, spatial, or mechanism classification.

If a process dies before any new durable record, recovery restores the exact
last durable cursor—`FRESH_READY` or `RETRY_PENDING`—rather than guessing from
the attempted call. If it dies after a durable rejection but before either successor
cursor or terminal record, recovery must reconstruct the exact typed branch
from the bound record: pending cursor for retry-required, terminal invalid
for exhaustion. It may not recompute or classify the boundary as fresh. An
unmatched rejection, forged successor, or branch mismatch is an invalid stop.
After cursor persistence, recovery restores that exact `RETRY_PENDING` cursor.

Recovery begins from the highest completely valid atomic checkpoint generation
and its exact six-member cursor, ledger, and accepted-state sets. It then scans
only the contiguous canonical journal suffix rooted at that checkpoint tip.
Only an exactly linked retry-required or retry-exhausted suffix may be
reconstructed; every other uncheckpointed, torn, noncontiguous, forged, or
mixed-generation suffix stops invalid. A syntactically present journal record
does not otherwise activate state. If a crash leaves either the old or new
checkpoint directory entry visible, complete hash-chain and generation
validation—not filename recency—selects it.

Every future journal append must write one canonical record, flush, and fsync
before returning. Every future checkpoint must bind the canonical ordered set
of all six member cursors, all six complete-ledger hashes, all six accepted
state hashes, common-event identity, campaign generation, and journal tip; its
write requires file fsync, atomic replace,
and parent-directory fsync. A successful fine step must publish its accepted
state, updated TDG6 ledger, and new `FRESH_READY` cursor as one atomic
checkpoint generation. A partial or mixed six-member set is invalid. A
common-event abort restores the complete verified ancestral six-member
snapshot and appends one campaign-level rollback record. That rollback creates
a new campaign generation and new generation for every member cursor, binds
the aborted/restored cursor, ledger, and state sets, and issues new attempt
identities. It never reuses old cursor bytes or an aborted attempt identity.

## Compact evidence only

The freeze binds sealed commit
`64316d8f80f8532524c306817239f823a73b14be`, the compact tracked PROTO14
freeze, HLT12 authorization, CAL11 invalid-result binder, and TDG7-IMP3
repair. It verifies their bytes both at that immutable commit and in the
current worktree. It also binds the live PROTO15 protocol/configuration and
pure validator by SHA-256.

It derives all six restart records only from the compact PRO14 and HLT12 JSON
artifacts. It does not open the earlier raw run bundle, state arrays, history,
or checkpoints. The six fixed payload identities remain the GR-0 amplitude-3
restart inputs at `t=23/16`; no PROTO14 post-restart or terminal state may be
used.

Before generating its canonical result, the freeze proves that these roots do
not exist and creates neither:

```text
runs/fgc-2-sf1/proto15/calibration
runs/fgc-2-sf1/proto15/holdout
```

That absence is temporal prelaunch evidence. A later runtime binder must bind
this freeze from immutable history, rather than attempting to re-observe the
same absence after execution has started.

## What it does not establish

PROTO15 is not a solver, a run plan, a runtime authorization, a calibration,
or a new trajectory. It does not make GR-0 eligible and does not authorize
SGB-L, FGC-QR, COL1, DEF1, retained-EFT evolution, a physical transition,
singularity resolution, a child domain, dark-sector mechanism, or varying
locally measured light speed.

The successor implementation owner is `FGC-1-HLT13-MON13`. It must implement
and synthetically qualify the durable cursor, both typed post-rejection
branches, recovery, complete six-member checkpoint and common-event rollback,
and fault injection at every journal/cursor/terminal/checkpoint publication
window before it can ask to
create a fresh GR-0 calibration namespace.

## Reproduction

```bash
python3 scripts/reproduce_fgc_pro15_frz1.py --verify
python3 -m unittest tests/test_fgc_pro15_frz1_reproduction.py -v
```

Both commands are freeze-only. They do not launch a campaign.
