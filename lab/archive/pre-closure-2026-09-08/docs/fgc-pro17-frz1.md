# FGC-1-PRO17-FRZ1: exact successor and external generation-zero freeze

**Status:** a premise-only construction contract. This artifact freezes no
runtime, physical input bundle, production `GenesisSpec`, state object,
checkpoint, namespace, trajectory, calibration, or candidate result.

PROTO16 correctly removed the `COMMON_EVENT_COMMIT` receipt hash cycle, while
PREF24 showed that its successor and genesis constructions remained
under-specified. PROTO17 supplies those missing rules prospectively.

## Exact successful common-event edge

The successful event graph is acyclic:

```text
complete predecessor checkpoint
        -> predecessor-only COMMON_EVENT_COMMIT receipt
        -> six derived successor cursors
        -> derived successor checkpoint
```

All six members must be `FRESH_READY`, at the active target, and have zero
current TDG6 retries. The receipt contains predecessor facts only. It has no
own digest, successor cursor, successor set, or checkpoint digest. Its outer
canonical record digest exists before successor construction and becomes every
successor cursor's journal tip. The successor checkpoint advances exactly one
event and one `1/16` target cadence, preserves the state and ledger sets, and
derives new cursors and high-water maps. A recovery can reconstruct exactly one
lone durable common-event receipt; every other uncheckpointed suffix is
invalid.

Receipt `sequence` is also derived: it equals the length of the validated,
hash-contiguous supplied prior journal plus one. It is `1` for the first event
only when the generation-zero journal is empty. This pure gate validates the
record envelope, digest chain, and sequence only; replay of prior journal
payload semantics remains an HLT14 runtime obligation.

## External trusted genesis

A future HLT14 authority must embed a complete visible `GenesisSpec` and bind
it through an **external** tuple of authorization commit, tracked result path,
raw result hash, and GenesisSpec hash. Neither the GenesisSpec nor the derived
generation-zero checkpoint may embed the GenesisSpec hash; the authority result
may not self-reference its own commit, result, or spec hash.

The GenesisSpec must expose six ordered descriptors, complete initial TDG6
ledger mappings, canonical state-object mappings, root cursor inputs, and a
physical-input manifest. Genesis begins from accepted restart boundary
`23/16` with active first target `3/2`. The manifest is a bare visible object whose sibling
`physical_input_manifest_sha256` is computed over its canonical bytes; the
manifest never contains its own digest. Each source bundle names its required
NPZ member keys. Each physical array is bound by source bundle/key,
shape, little-endian `float64` C-layout, and raw-byte identity. The derived
generation-zero store contains only content-addressed state objects, an empty
journal, and its one deterministic checkpoint.

Those raw array identities are declared external inputs at this freeze; the
freeze does **not** open raw NPZ bytes or recompute array, physical-state, or
restart-payload hashes. HLT14's adapter must do that work before any namespace
is materialized.

There are deliberately two distinct hash meanings:

- `physical_state_sha256` is the existing normalized `u,p,q` array-content
  hash. It identifies physical field arrays.
- `state_object_canonical_sha256` is the canonical serialized state-object
  content address. This is the value placed in a PROTO15-style cursor and
  checkpoint as `accepted_state_sha256`.

They must never be silently substituted for each other.

At launch the runner reads the authority result via `git show` at the exact
full commit, verifies its tracked blob and all pinned source bytes, and rejects
import shadowing. The assumed trustworthiness of that Git history is explicit
and external; local hash consistency does not prove it.

The mutable output namespace is not authority. It must be absent at the freeze
boundary, then be created only from the verified external tuple, GenesisSpec,
and physical input. A locally self-consistent foreign store, changed physical
bundle, or reused namespace is rejected rather than adopted.

## Frozen boundary

The freeze binds immutable HLT13, PROTO16, and PREF24 records at commit
`03f1f249a9484ba43513c3ad4292cc993598a65b`. It does not read their raw
checkpoints or create a replacement state. Concrete production descriptors,
physical-bundle hashes, and generation-zero values remain future HLT14/adapter
inputs. Therefore no runtime, pretrajectory, calibration, GR-0 eligibility,
SGB-L, FGC-QR, DEF1, EFT, transition, or physical claim is authorized.

## Executed pure-oracle qualification

The tracked side-effect-free oracle executes one synthetic complete genesis,
one synthetic synchronized common-event successor, and exact lone-suffix
recovery. It also rejects exactly fourteen listed construction mutations:
physical-state digest, state-object digest, manifest digest, array storage key,
input hash, zero-ledger invariant, frozen restart time, foreign protocol,
receipt-cycle injection, receipt sequence, untouched successor method,
successor high-water map, successor checkpoint digest, and nonzero-generation
root-journal inconsistency. These are construction controls, not historical
payload replay or raw-bundle verification.

The remaining HLT14 requirements are frozen but not passed controls in this
record: raw-bundle byte/hash recomputation, external Git authority/blob/live-
import validation, semantic replay of historical journal payloads, and
namespace-reuse/foreign-store rejection. Each is recorded as `required=true`,
`completed=false`; none is relabelled as a present mutation pass.

Reproduce the tracked evidence without touching any run directory:

```bash
python3 scripts/reproduce_fgc_pro17_frz1.py --verify
python3 scripts/reproduce_fgc_pro17_frz1.py --verify-freeze-boundary
python3 -m unittest tests/test_fgc_pro17_frz1_reproduction.py -v
```
