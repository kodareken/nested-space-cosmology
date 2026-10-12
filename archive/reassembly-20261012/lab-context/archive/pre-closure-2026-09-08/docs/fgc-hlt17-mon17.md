# FGC-1-HLT17-MON17 — IMP1 durable progression qualification (prospective)

This is the implementation map for the next dependency of PRO20-EV1-FRZ1,
not a completed certificate or permission to run a campaign. The accepted
[PLAN](../PLAN.md) owns the complete SF1 finish line. [IMP1](fgc-tdg11-imp1.md)
qualifies a synthetic temporal-admission kernel; it does not authenticate a
physical source, a historical origin, a durable cursor, or this protocol.
AGENTS.md stays unchanged creative context, not an executable premise.

## Fixed production boundary

The new protocol is `FGC-2-SF1-PROTO19`. Its prospective first-event store is
`runs/fgc-2-sf1/pro20-event1/calibration`. This is not the historical directory
named `proto19`, which belongs to PROTO18/RCV3 and remains sealed.

The only allowed production origin is the synchronized HLT15 `t=23/16`
generation in `runs/fgc-2-sf1/proto17/calibration`; the first common target
is exactly `t=3/2`. RCV3, QA and MSEL states are not origins. IMP1's
`historical_origin_authenticated=False` is not rewritten: HLT17 must provide
and independently verify that external proof before production seeding.

Preserve the existing six members and their method-owned operators:

| Member | Owned rows N | Fallback D01 ceiling | Fallback D12 ceiling |
|---|---:|---:|---:|
| RK4-2049 | 2044 | 16352 | 32704 |
| RK4-4097 | 4092 | 32736 | 65472 |
| RK4-8193 | 8188 | 65504 | 131008 |
| SSPRK3-4097 | 4092 | 32736 | 65472 |
| SSPRK3-8193 | 8188 | 65504 | 131008 |
| SSPRK3-16385 | 16380 | 131040 | 262080 |

Here N=point_count-5 and the ceilings are 8N/16N. HLT17 injects these exact
limits into the named C1R1 evaluator, which keeps the IMP1 enclosure-limit
type as reference-wire. Other numerical thresholds remain unchanged. In
particular all eighteen complete-state channels and p>=3/2 remain required. The IMP1 2,049-point synthetic timing is not a full-ladder
or physical-RHS qualification. Actual source/operator and largest-grid
resource qualification must be measured before the freeze can authorize
production. No grid reduction or diagnostic hybrid is substituted.

## First bounded implementation slice

The first slice is synthetic, state-machine and replay qualification only:

- `evolution/hlt17_imp1_cursor.py`: new strict immutable cursor and canonical
  codec, six-member limits, full TDG7 plans, accepted identities, independent
  source/CFL/temporal counters, and the composite retry state below.
- `evolution/hlt17_imp1_bridge.py`: exact persisted-rejection reconstruction,
  callback/boundary authentication, and a narrowly pinned reduced-plan
  adapter. The in-memory path now prepares with named C1R1; it does not run
  sealed IMP1 `_prepare` and replace that assessment. It does not read a
  real store or authenticate an origin.
- `tests/test_fgc_hlt17_imp1_cursor.py` and
  `tests/test_fgc_hlt17_imp1_bridge.py`: exact/small synthetic success,
  restart, mixed-owner, malformed-input and injected-failure controls.

These modules have new schema identities. They do not coerce an IMP1 ledger
into TDG6, or modify HLT16/PROTO18, IMP1, TDG6 or TDG7. Private C1R1 prepare
helpers and IMP1 ledger/wire helpers are explicitly named and pinned rather
than silently becoming a new historical public API.

## Complete replay identity

Persist the complete TDG7 plan, including requested_cap separately from the
executed macro_width, all five boundaries, quantum, selection budget,
conservative reduction, target_limited and minimum-width field. Width-only
reconstruction is insufficient even when two caps produce equal endpoints.

A pending temporal cursor binds the canonical preceding IMP1 rejection,
the full preceding plan, the predecessor IMP1 ledger and its SHA-256, the
updated ledger and its SHA-256, preparation and assessment hashes, method,
member, target, accepted time, step and transaction serial. It separately
binds the serialized state descriptor hash and the physical u/p/q state hash;
they are not interchangeable. Coordinates, source/operator closure, monitor,
causal state, tracer and guard configuration identities must stay visible.
Checkpoint generations and the journal parent/tip are separate identities.

Replay reconstructs the same rejected proposal using the saved predecessor
ledger and complete plan. At retry depth zero the public C1R1 fresh API
suffices; at deeper depth the adapter uses the explicitly pinned C1R1
preparation instrument. Compare the entire canonical rejection/preparation
receipt, not just its classification. Receipt hashes are C1R1-named; the
rejection event type remains the IMP1 ledger wire. Replay never writes the
old rejection again, spends a new retry, or adopts a shadow. A changed or
failing replay is an invalid premise, not a new source/CFL/temporal attempt.

## Mixed retry contract

The historical HLT16 attempt layer can source/CFL-fail while preparing a
temporal retry, but its exclusive pending-owner cursor cannot persist that
combination. HLT17 must represent it; neither resetting a counter nor making
every mixed retry terminal is an acceptable substitute.

`FRESH_READY` has no active temporal rejection. `RETRY_PENDING` retains the
exact IMP1 half-width successor of the preceding temporal rejection. An
additional source/CFL overlay owns any further strict reduction. Do not
replace `IMP1RetrySuccessor.plan` with that overlay plan: its contract is the
exact temporal half-cap. Keep the updated IMP1 ledger and inherited TDG6
snapshot unchanged through each source/CFL overlay transition.

The overlay records its owner, full actually attempted plan, strictly smaller
next TDG7 plan, rejection evidence and independent counter. Source and CFL
counters may both remain nonzero while the latest overlay identifies the
next attempt. A fresh source/CFL retry uses the same overlay representation
without fabricating temporal evidence. A new temporal rejection after an
overlay appends to IMP1 history and replaces the preceding temporal identity
with the exact actually attempted overlay proposal; source/CFL counters do
not disappear merely because ownership changes.

If a source/CFL reduction cannot produce a strictly smaller admissible
lattice plan, report that owner's exhausted retry with the unchanged IMP1
ledger. Replay failures remain invalid. Successful fine-only adoption resets
the per-macro counters only at the accepted boundary and retains cumulative
counts, historical rejections and exact non-cancelling debit.

## Remaining end-to-end owners

### Joint in-memory checkpoint contract

Before adding filesystem publication, join the reviewed cursor/bridge and
member/codec in `evolution/hlt17_member_checkpoint.py`, with focused tests in
`tests/test_fgc_hlt17_member_checkpoint.py`. This is still a synthetic,
no-I/O integration layer, not a completed MON17 certificate.

Keep five different identities distinct: physical u/p/q state, accepted
descriptor, full kernel cursor, kernel-transition chain, and actual future
store journal. The codec's projected cursor-identity hash is not the hash of
the full cursor bytes. Accepted physical generation counts accepted IMP1
macro steps; cursor/attempt/checkpoint/journal revisions are not that count.
The unsealed cursor now names those prototype hashes
`kernel_transition_parent_sha256` and `kernel_transition_digest`. They are
not a store-journal parent/tip. Actual FUTURE store journal and store
checkpoint revisions stay unbound until the external journal owner exists.

An initial checkpoint first encodes the provided accepted boundary against
its supplied PROTO17 origin identity, then seeds the kernel cursor with the
resulting actual descriptor digest. Live members require an explicit origin
cursor reference; synthetic fixtures may use a labeled placeholder digest
that is not a manufactured real PROTO17 cursor. Matching hash text remains
unauthenticated. For a fine adoption, the new descriptor must name the old
full cursor; the new cursor must name the new descriptor. There is no
circular hash. Finalize the descriptor reference inside one checked
transition; do not leave an old descriptor attached to a new state hash or
repair it by an unchecked dataclass replacement.

Export is a validated immutable generation bundle of descriptor, payload,
and full-cursor bytes. Hashes, revisions, projections, and authorization
flags are derived after codec and cursor reauthentication. Construction
and dataclass replace cannot forge those identities. agree() and export
reauthenticate the accepted encoding through the codec and cross-check the
entire live boundary with the cursor, including callable binding, monitor,
causal, guard, tracer, and event history. A retry overlay may only be the
accepted ledger or an exact current-step rejection suffix of it; source and
CFL totals stay distinct and are not reset to the accepted baseline. Restore
rejects a caller-supplied last accepted plan that differs from the accepted
descriptor.

An isolated provisional copy may be used to construct the prospective
accepted encoding. Recheck the real IMP1 commit against it before adopting
the result. Never call the kernel commit twice on the same live transaction,
scrape an unqualified corrected endpoint, or retain partial monitor/tracer
advancement if encoding, cursor finalization or result construction fails.

At restart, authenticate the accepted descriptor and full cursor bytes,
restore the accepted member, then apply the checked cursor/IMP1 retry overlay.
Source/CFL-only, temporal and mixed retries all retain accepted state/time,
the immutable TDG6 sibling and the last accepted executed plan. Restoring a
retry overlay must not silently erase that plan. The current retry plan is a
different object. Source/CFL totals must agree between member and cursor,
including declared inherited baselines; counters are not silently zeroed.

Joint tests cover initial encode/restore, actual fine commit with a new
descriptor, first acceptance followed by temporal/source/CFL retries and
restart, both original integrators, foreign cursor/descriptor/closure/ledger
swaps, all failed-transition rollback boundaries, and reaching a common
target without granting a new target or invoking a source again. No
durability, source-manifest or historical-origin authentication is inferred
from this in-memory checkpoint.

### In-memory C1R1 admission integration

This in-memory path uses one fixed C1R1 implementation. It is not a
user-extensible or runtime-pluggable backend. IMP1 remains the sealed
mathematical and ledger-wire reference. C1R1 assesses the same binary64
1/2/4 family, preserves all eighteen channels, `p>=3/2`, full
non-cancelling debit, and the original binary64 fine endpoint. It does
not claim that IMP1 executed.

Implementation identity bound into the cursor, member encoding, and
receipt:

| Field | Value |
|---|---|
| Cursor schema | `FGC-1-HLT17-C1R1-cursor-v1` |
| Member codec schema | `FGC-1-HLT17-C1R1-member-state-v1` |
| `implementation_id` | `tdg11_c1r1_integer_exponent_direct_ring_v1` |
| `mathematical_object` | `exact_accumulation_reconstruction_with_debit` |
| `reference_wire_artifact_id` | IMP1 artifact `FGC-1-TDG11-IMP1` |
| `reference_wire_evaluator_id` | `tdg11_imp1_exact_bernstein_with_dual_fallback_v1` |
| `admission_runtime_id` | `hlt17_c1r1_admission_runtime_v1` |

`evolution/hlt17_admission_runtime.py` pins C1R1 `_prepare`, plan,
revalidate, rejection, and receipt helpers. Sealed IMP1 `_prepare` is not
in that pin set. Foreign IMP1 prepared/committed objects, the previous
`FGC-1-HLT17-IMP1-cursor-v1` / `FGC-1-HLT17-member-state-v1` encodings,
and a swapped implementation marker are refused. Mathematical/ledger-wire
compatibility is not type-identity with IMP1 execution.

Source, physical-origin, environment, resource, store, and runner
authority remain external and unqualified. This is synthetic-component
integration, not production activation. `protocol_v19.py` now derives and
persists the five C1R1/reference-wire identity fields in every checkpoint
member state from reauthenticated descriptor/cursor bytes. The store still
does not authenticate the source, origin, environment or production authority.

### Publication and production qualification

This bounded slice is the synthetic `protocol_v19.py` checkpoint/journal
contract and `hlt17_campaign_store.py` publisher. It consumes validated
member bundles of three immutable byte-string components (`descriptor`,
`payload`, `cursor_bytes`) and revalidates them from those bytes; supplied
derived labels are not authority. Tests use synthetic temporary directories
only. This is not a completed MON17 certificate, not a live runner, and not
production write authority.

Public schema APIs live in `evolution/protocol_v19.py`:

- `build_root_manifest` / `validate_root_manifest`
- `build_seed_generation` / `build_successor_generation` /
  `parse_validated_generation` / `authenticate_generation_chain`
- `production_member_identities` (exact six production N and 8N/16N
  ceilings, without executing those sizes)
- session and terminal-lock encoders
- `HLT17ValidatedGeneration`

Public store APIs live in `evolution/hlt17_campaign_store.py`:

- `HLT17CampaignStore.create` / `open`
- `acquire_writer` / `close_writer` / `abandon_writer` (`HLT17WriterCapability`)
- `publish_seed` / `publish_attempt`
- `authenticate` / `authenticate_read_only` (`HLT17AuthenticatedView`)

Atomic layout of one synthetic store, with no mutable latest pointer:

```text
<root>/
  root_manifest.json          # immutable; published once
  writer_exclusion.lock       # process/thread exclusion; not session authority
  session-NNNN-open.json      # durable session ownership
  session-NNNN-close.json     # explicit close; required before a later session
  terminal.lock               # typed terminal; forbids further appends
  generation-NNNNNNNNNNNNNNNN/
    generation_manifest.json  # binds journal, checkpoint, and new blobs
    checkpoint.json           # names journal_record_sha256
    journal.json              # no successor checkpoint hash
    blobs/<sha256>            # newly referenced descriptor/payload/cursor only
```

Identity relations kept distinct:

- checkpoint hash is SHA-256 of `checkpoint.json` bytes
- journal hash is SHA-256 of `journal.json` bytes; the checkpoint points at
  that hash and the journal does not contain a successor checkpoint hash
- member full-cursor hash is SHA-256 of the immutable cursor bytes
- kernel-transition digest remains the cursor-chain identity, not a store
  journal parent or store generation
- store generation and journal sequence count publications; they are not
  accepted generation, cursor generation, attempt serial, or kernel
  checkpoint generation
- accepted fine increments accepted generation and 4 substeps plus 20 RK4
  or 16 SSPRK3 records; source/CFL/temporal retries preserve accepted
  state/time/blobs and the immutable TDG6 sibling while their own
  counters/history advance, keeping inherited source/CFL totals

The complete predecessor/successor/cohort transition is validated before
any directory becomes visible. One member may change per attempt. Seed
uses the truthful RK4-9/SSPRK3-9 cohort at origin `t=23/16` with target
`t=3/2`. First-event completion requires every cohort member fresh and
parked at that target with no pending retry and is a final event boundary.
If the immutable root manifest is visible but writer-exclusion publication
fails, the overall store creation is already post-publication uncertain; the
partial root is preserved and cannot be opened as a store.
Typed terminals preserve the partial history, publish `terminal.lock`, and
forbid further appends. If the generation is already visible and
`terminal.lock` cannot be published, that is post-publication uncertainty:
the generation is preserved, `published=True`, and further mutation is
poisoned. Temporal-exhaustion evidence must be the exact one-rejection
successor of this member's cursor ledger, bound to method, origin,
state/time/step/serial, inherited TDG6, debit, and history prefix; it is
not an executable retry cursor. Derived generation and view fields are
recomputed from immutable documents. Reads reauthenticate root-manifest
and exclusion-lock bytes plus root/held-lock inodes. A free OS lock is not
authority to take over an unclosed session; `abandon_writer` may release
process resources without publishing session close. There is no automatic
stale-writer recovery. Complete visible bytes after uncertainty may be
read as evidence and are not permission to resume an unclosed or uncertain
writer. Reads reject forks, holes, partial generations, unexpected files,
link/path races, and post-publication uncertainty without repairing them.

Remaining duties are outside this slice and stay explicit seams:

- independently verify the sealed HLT15 `t=23/16` origin
- authenticate physical source/operator and configuration closure
- environment, compiler, and largest-grid resource gates
- an independent PRO20 binder that does not trust runner labels
- a live runner or Make target, only after that committed authority

This layer does not authorize writes under
`runs/fgc-2-sf1/pro20-event1/calibration`, does not execute production-size
preparations, and does not claim calibration, eligibility, or physics.

HLT17 qualification still has to exercise physical-source/operator
integration, all six production sizes, mixed retries against a real origin,
and six-member common-event synchronization. PRO20's separate prospective
authority freezes the source closure, environment, limits, exact origin,
new namespace and one-at-a-time execution. Compact verification remains
raw/store/shadow/Git-blind.

No compact result, live target, new raw directory or scientific state
advance is justified by this document or this synthetic publication slice.

## Proof and failure forecast

Require round-trip identity of every persisted field and the full inherited
snapshot; true fresh/retry API separation; temporal -> source -> CFL -> retry
and restart controls; equal-width/different-requested-cap replay controls;
preserved counters/debit through all failed preparations; and late-failure
rollback of accepted state, monitor, causal and tracer histories. Test unknown
keys, bool/int aliases, noncanonical numbers, swapped descriptor/state hashes,
altered plans/closures/generations and attempted state promotion.

The likely failures are lost replay identity (compare whole receipts), an
illegal mixed-owner durable transition (qualify the composite state before
writing), or actual full-ladder resource/physical-source rejection (preserve
its typed evidence, not a smaller grid or a relaxed numerical threshold).
