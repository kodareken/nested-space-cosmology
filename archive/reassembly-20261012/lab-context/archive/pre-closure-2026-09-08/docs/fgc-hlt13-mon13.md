# FGC-1-HLT13-MON13: PROTO15 durable-cursor runtime qualification

**Status:** synthetic implementation and qualification only. This artifact
implements the prospective `FGC-2-SF1-PROTO15` campaign-owned cursor contract
on deterministic synthetic state objects. It reads no historical PROTO14 raw
bundle, state array, history, or checkpoint; it creates neither PROTO15 output
root and advances no numerical trajectory.

## Question owned by this gate

HLT13 asks whether the frozen durable lifecycle can be implemented without
turning an external caller's assertion into authority. The runtime owns one
durable cursor for each of the six fixed members and distinguishes exactly
`FRESH_READY` from `RETRY_PENDING`. Its tests exercise only content-addressed
synthetic objects whose identities are derived from public compact lineage;
they are not restarted physical states.

The implementation is bound to immutable commit
`f3812d0b93e318890e973b4b21f681dbde8f437e`, its compact PRO15 freeze and
protocol, and the sealed TDG7-IMP3 repair. It verifies every bound file both
from that commit and from the live worktree. The historical PROTO14 raw
campaign remains closed and is not an HLT13 input.

## Implemented durable instrument

Every state object is content addressed. Each journal append is one immutable,
canonical compact sorted UTF-8 JSON record with the exact envelope
`record_kind`, `journal_parent_sha256`, `payload`, and `record_sha256`; the
record hash is calculated before its self-hash field. Cursor and complete TDG6
ledger hashes use the same canonical discipline. Time identity contains both
the frozen rational value and lowercase binary64 hexadecimal representation.

The runtime persists attempt and generation high-water marks. Thus a retry,
restart, or campaign rollback cannot reuse an attempt identity. Checkpoint
selection validates complete hash chains and generations rather than file
names or timestamps; competing valid-looking checkpoint packages at one
generation are a fork and stop invalid.

After `TDG6TemporalRetryRequired`, the campaign owner durably records the
rejection, derives the exact successor, durably adopts `RETRY_PENDING`, and
only then permits retry entry. Recovery may reconstruct precisely that
post-rejection branch. It cannot classify an unmatched record as fresh.
Every accepted fine commit carries its complete executed TDG7 plan. A fresh
commit must start at the cursor's exact accepted boundary, share its event
target, and end at the new accepted boundary. A pending retry additionally
requires that this plan equal the already durable successor byte for byte.
Recovery repeats the same checks from the persisted journal/checkpoint graph;
even a rehashed alternative plan with the same endpoint is rejected.
After `TDG6TemporalRetryExhausted`, the campaign durably closes exact
`invalid_implementation_or_nonconverged_run` with a terminal lock. The terminal
checkpoint preserves the unchanged six nonterminal cursor records while
embedding the updated terminal member ledger; it never creates
`RETRY_PENDING`.

An accepted fine commit publishes one complete six-member cursor/ledger/state
checkpoint generation. An uncheckpointed accepted or rollback suffix is
invalid. A common-event rollback restores only a verified ancestral complete
snapshot and creates new campaign, cursor-generation, and attempt identities;
it does not reuse discarded cursor bytes.

PROTO15 names two rollback facts `new_journal_parent_sha256` and
`new_journal_record_sha256`. HLT13 resolves them explicitly from the immutable
journal envelope as `journal_parent_sha256` and `record_sha256`. The latter
cannot also live inside the hash-covered payload without becoming a recursive
self-hash. Recovery checks that the resolved parent is the aborted
checkpoint's journal tip and the resolved record hash is the new checkpoint's
tip; the sealed PROTO15 bytes remain unchanged.

Named deterministic fault hooks cover journal write/fsync, cursor transition,
terminal record, checkpoint file fsync, replacement, parent-directory fsync,
and rollback-record windows. Qualification must demonstrate the stated
recovery or fail-closed outcome at every such publication boundary.

The executed qualification contains **110 deterministic cases**: six nominal
lifecycle controls, fourteen adversarial corruption and self-consistent
semantic forgeries, and ninety fault cuts covering nine publication scenarios
across all ten durability windows. The last four forgeries cover a
nonadvancing accepted time, target overrun, terminal substitution of an older
valid rejection, and altered rollback cursor content. The `after_replace` and
`before_dir_fsync` windows
may expose either an old or new directory entry after a real crash, so all
eighteen such cases are replayed as both filesystem images. This yields 108
recovery images for the ninety fault cuts. Every case passes, and the ordered
case corpus is bound by digest
`0acc88f1783f54f4c2b57fedd629f0a940945d1a01580054c7e35f099e8d8435`.

The hash threat model is deliberately narrower than authentication. HLT13
detects stale hashes, partial writes, forks, and cross-record semantic
inconsistencies relative to its trusted genesis. Because every hash is locally
recomputable, it cannot authenticate wholesale replacement of the entire
campaign store by a different self-consistent genesis and history. Any later
production authorization must bind that genesis to an external immutable
manifest or commit; HLT13 does not claim to provide that authentication.

## Claim and nonclaim boundary

The only HLT13 implementation gates opened are:

```text
PROTO15_successor_runtime_implemented = true
PROTO15_runtime_synthetic_qualification_passed = true
```

HLT13 does **not** authorize a pretrajectory runtime, creation of either
`runs/fgc-2-sf1/proto15/calibration` or `runs/fgc-2-sf1/proto15/holdout`, fresh
GR-0 calibration, an eligible GR-0 case, SGB-L, FGC-QR, DEF1, retained-EFT
evolution, a transition, or any physical result. PROTO14 remains an historical
invalid-runtime/nonconverged outcome, not a physical obstruction.

The canonical result is compact and hash-bound. It is not a run plan, a
production campaign, or evidence about collapse, trapping, defocusing,
singularity resolution, child domains, dark-sector mechanisms, or locally
measured light speed.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt13_mon13.py --verify \
  --output results/fgc-1-hlt13-mon13.json
python3 -m unittest tests/test_fgc_proto15_runtime.py \
  tests/test_fgc_hlt13_mon13_reproduction.py -v
python3 scripts/check_repo.py --only-hlt13-mon13
```

These commands inspect compact tracked lineage and deterministic synthetic
qualification only. They must leave both PROTO15 output roots absent.
