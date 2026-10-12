# FGC-1-PRO17-PREF25: independent PROTO17 construction theorem

**Status:** independent, compact, synthetic construction binder. It authorizes
only implementation and synthetic qualification of a future
`FGC-1-HLT14-MON14`. It does not implement that runtime, create a namespace,
load a production state, resume a campaign, or evaluate a trajectory.

## Question

PREF24 found that PROTO16 removed the successful-event hash cycle but did not
fully specify how to construct generation zero or the event successor. PROTO17
then froze those missing rules in a pure reference constructor. PREF25 asks the
next narrow question:

> Can a separate standard-library implementation, without importing the
> PROTO17 schema or constructor modules, reconstruct the same visible schema,
> generation-zero checkpoint, predecessor-only `COMMON_EVENT_COMMIT`,
> successor checkpoint, and lone-suffix recovery?

The answer is yes for the complete declared synthetic fixture and mutation
domain. This is a theorem about the frozen abstract construction, not about a
physical state or a collapse.

## Immutable authority

The binder reads the PROTO17 protocol, freeze, result, reproducer, schema,
constructor, tests, and document from immutable commit
`f4d6e766aeb6351dc34ea482ddb3bd8eaac8d96d` using `git show`. Every bound blob
must match its literal SHA-256 and the present tracked file. The sealed result
must remain duplicate-safe canonical JSON.

The immutable commit is the binder's input. Current-worktree equality is a
drift check, not a replacement source of authority. Trust in the Git history
is still an explicit external assumption; PREF25 does not authenticate a
future launch authority or signing root.

The one-time PROTO17 namespace-absence observation is read from the sealed
result and is not reobserved against mutable `runs/` paths.

## Independent construction

The theorem module is standard-library-only and its import syntax is audited.
It imports neither `protocol_v17.py` nor
`proto17_pure_construction.py`. It independently fixes and checks:

- the complete `GenesisSpec`, member-descriptor, state-object, physical-array,
  and TDG6-ledger field tuples;
- the canonical six-member order and the PROTO12/RSP2 source split;
- the seven array names, exact shapes, little-endian `float64` C layout,
  unique NPZ storage keys, safe relative paths, and required `metadata_utf8`;
- distinct physical `u,p,q` content hashes, restart-payload hashes, and
  canonical persisted state-object addresses;
- eighteen nonnegative TDG6 debit channels and exact zero-initial ledgers;
- the event-23, `23/16 -> 3/2` generation-zero construction, content-addressed
  state paths, empty journal, and deterministic checkpoint path;
- a six-record hash-contiguous **structural** accepted-fine prefix;
- the sequence-7 predecessor-only common-event receipt;
- the event-24, `3/2 -> 25/16` successor cursors, high-water maps, state and
  ledger preservation, parent checkpoint, journal tip, and checkpoint hash;
- exact reconstruction from one lone durable receipt and rejection of every
  other suffix shape.

The independent constructor reproduces the three sealed PROTO17 synthetic
hashes exactly:

```text
generation-zero checkpoint  01e8214dcbbef025fe7883670043659d3ec66bb47e9b3e9dd4b8089b8ba5f8de
common-event receipt         17f925e8489341392ac5df4fba1cc47b66a27153e64eca4455c3f701f37cbc17
successor checkpoint         c4c31f53a6a97632addafd54d6b3ec08bdcffe4aa4321269f8c8bca957e5380b
```

The accepted-fine prefix is checked only as a canonical, hash-contiguous
record envelope. PREF25 deliberately does not replay the semantics of its
historical payloads; that remains an HLT14 obligation.

## Mutation boundary

The independent theorem rejects the same fourteen named pure-construction
mutation classes recorded by PRO17-FRZ1:

```text
array_storage_key
foreign_protocol
frozen_restart
input_hash
manifest_digest
nonzero_generation_root_journal
physical_state_digest
receipt_cycle_injection
receipt_sequence
state_object_digest
successor_checkpoint_hash
successor_high_water
successor_untouched_method
zero_ledger
```

The compact binder separately attacks duplicate JSON, immutable blob/hash and
commit drift, schema tuple and member-order changes, qualification count/name
changes, false completion of a deferred HLT14 duty, a forbidden import, and
promotion of a runtime or physical claim. These are synthetic structural
controls only.

## What remains open

Four HLT14 requirements remain exactly `required=true`, `completed=false`:

1. recompute raw bundle and physical-array bytes and hashes;
2. validate the external Git authority, tracked blob, and live imports at
   launch;
3. replay and validate the semantics of historical journal payloads;
4. reject namespace reuse and a foreign self-consistent store relative to the
   externally authenticated authority.

PREF25 therefore authorizes only the implementation and synthetic
qualification of HLT14 against the already-sealed construction. It does not
authorize pretrajectory operation, a production adapter or runner, namespace
creation, fresh GR-0 calibration, eligibility, SGB-L, FGC-QR, DEF1,
retained-EFT promotion, a transition, or any physical conclusion.

## Reproduction

The compact proof does not inspect `runs/`:

```bash
python3 scripts/reproduce_fgc_pro17_pref25.py --verify
python3 -m unittest \
  tests/test_fgc_proto17_construction_theorem.py \
  tests/test_fgc_pro17_pref25_reproduction.py -v
python3 scripts/check_repo.py --only-pro17-pref25
```
