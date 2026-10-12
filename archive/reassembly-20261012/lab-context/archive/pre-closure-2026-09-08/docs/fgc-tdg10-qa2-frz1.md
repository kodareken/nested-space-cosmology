# FGC-1-TDG10-QA2-FRZ1 — two-width exact complete-C robustness freeze

## Decision

`FGC-1-TDG10-QA2-FRZ1` prospectively freezes one bounded, read-only
two-width robustness test of the already-qualified retry-3 exact complete-C
interval owner on the inherited RK4-2049 SBP4 operator. It does **not**
authorize a campaign-state advance, a production method, or old-member
adoption.

QA1-PREF1 independently bound an all-18-channel pass on retry 3 only. That
result licenses this QA2 method-design freeze. It does **not** license
transplanting the diagnostic endpoint onto another member, adopting an
old-member state, or treating retry-3 success as a general production remedy.

The qualification changes no numerical owner relative to QA1. For each of
retry 4 and retry 5 it reconstructs seven SSPRK3 tableau shadows on the
inherited SBP4 operator from that width's accepted predecessor and assesses
all 18 complete-state channels with the already-qualified tdg10 exact
radius-free complete-C dual-rational implementation. The discriminator
remains

\[
8U_{12}^{2} \le L_{01}^{2},
\]

including equality. Exact ceilings and depth are unchanged. No tolerance or
clipping is permitted.

## Correct retry boundaries

Retry 4 starts from the generation-10 checkpoint
`6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187`,
journal sequence/tip 12
`d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4`,
prior retry count 3, and attempted width `0x1.aaa9612df8000p-12`. Sequence 13
remains the historical retry-4 rejection. Generation 11 is not an executable
retry-4 predecessor.

Retry 5 starts from the generation-11 checkpoint
`11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8`,
journal sequence/tip 14
`0810018134df958b8dab32bd9c0ecd7ebfc4376fbc87869a8bee7dd02eb1e199`,
prior retry count 4, and attempted width `0x1.aaa9612df0000p-13`. Sequence 15
remains the historical retry-5 rejection. Generation 12 is not an executable
retry-5 predecessor.

Generation 9 / retry 3 is outside this freeze. The sealed 115-leaf store
`runs/fgc-2-sf1/tdg8-rcv3/calibration` is authenticated read-only.

Accepted time, physical state, member descriptor, coordinates, cursor mode
`RETRY_PENDING`, pending owner `temporal`, ledger owner `TDG6TemporalLedger`,
tracer owner `NormalFlowTracers`, and transaction owner
`GR0RuntimeStageTransaction` are the committed TI2/LOC2/PREF2 identities.
Per-width `transaction_sha256` values bind the restored monitor, causal, and
ledger payload.

## Authorized diagnostic transaction

The one authorized run:

1. authenticates the sealed store and the two selected RK4-2049 predecessors;
2. restores each finite predecessor without mutating the store;
3. constructs seven SSPRK3 shadows per width on the inherited SBP4 operator;
4. reconstructs and classifies all 18 exact complete-C row streams per width;
5. publishes only canonical `manifest.json` and `terminal.json` under the
   new absent gitignored namespace
   `runs/fgc-2-sf1/tdg10-qa2/retries4-5-ssprk3-sbp4-exact-complete-c`; and
6. proves the historical store is unchanged.

There is no temporal admission commit, no state/cursor/ledger/journal/checkpoint
write, no endpoint transplant, and no diagnostic endpoint serialization. No
retry beyond the two selected widths is authorized.

## Valid terminal classes

- `completed_both_width_all_channel_pass_no_state_advance`
- `completed_one_or_more_width_all_channel_nonpass_no_state_advance`
- `exact_route_or_resource_inconclusive_no_state_advance`
- `shadow_or_scientific_premise_stop_no_state_advance`
- `invalid_provenance_or_implementation`

The terminal may state `width_robustness_passed` only if both widths have all
18 channels pass. No class may mutate the historical store or publish an
accepted campaign descriptor, cursor, journal, checkpoint, temporal ledger,
common event, or diagnostic endpoint.

## Frozen successor rule

If both widths pass, this freeze licenses design of a truthful production
method or admission owner. If any selected width nonpasses, this exact
SSPRK3+SBP4+exact-C remedy is rejected as a general production remedy on the
tested retry neighborhood. No threshold or method change is permitted after
the outcome.

A pass still does not earn a production method, state advance, two-method
agreement, common event, GR-0 calibration, candidate, mechanism, or physics.
QA1 does not license old-member adoption.

## Frozen exclusions

No retry 3, fourth width, new grid, threshold reduction, diagnostic endpoint
serialization, endpoint transplant, old-member adoption, RA1, candidate
branch, campaign-state write, publication, push, release, or external contact
is authorized.

Historical TDG5 through TDG10-QA1, HLT16, PROTO15, SID1, RCV3, and their
hash-bound evidence remain unchanged. `.qdrant-initialized` remains
unrelated.

## Compact verification

Ordinary verification consumes only the tracked compact certificate and is
store-, raw-, and shadow-blind. The compact freeze is unexecuted: zero shadows
measured, no channel outcome, `width_robustness_passed=false`, and no state
advance. Namespace absence is a one-time prelaunch observation frozen into
that certificate.

```bash
make fgc-tdg10-qa2-frz1
make verify-fgc-tdg10-qa2-prelaunch
python scripts/check_repo.py --only-tdg10-qa2-frz1
```

`make run-fgc-tdg10-qa2` is explicit and absent from ordinary verification.
