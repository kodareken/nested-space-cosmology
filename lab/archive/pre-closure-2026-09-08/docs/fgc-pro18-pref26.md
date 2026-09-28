# FGC-1-PRO18-PREF26: actual source-evidence binder

**Status:** completed bounded read-only binder. The frozen contract was executed
against the two real legacy source checkpoints and their matching append-only
event logs, producing the canonical `FGC-1-PRO18-PREF26` evidence record. It is
not a launch, a namespace operation, a calibration, or a collapse calculation.

## Purpose

`FGC-1-PRO18-PREF26` follows the sealed `FGC-1-PRO18-FRZ1` checkpoint at
`8a41a70606dbdc0092e45eed215fd108b35fac9b`. It binds the actual legacy
evidence that PRO18 deliberately did not open:

- `PROTO12`: one shared checkpoint and one 49-line legacy event log supplying
  five logical selectors;
- `RSP2`: one shared checkpoint and one 25-line legacy event log supplying the
  sixth selector.

The six logical restart identities, their fixed order, and the two-container
grammar remain owned by PRO18. PREF26 may verify that the real bytes, NPZ
structure, selected payloads, and source-specific journal semantics agree with
that frozen map. It must not replace the map, normalize the history, or infer
a physical outcome from it.

## Read-only source contract

The binder has a fixed maximum source budget of 27,053,743 bytes: 5,942,233
bytes for the two NPZ containers and 21,111,510 bytes for the two event logs.
It also fixes 53 total NPZ members, at most 28,083,115 uncompressed NPZ bytes,
and 74 total journal lines. The bounds make archive expansion and journal
parsing explicit, bounded operations rather than hidden environmental input.

The completed execution performed all of the following without mutation:

1. hash the four declared source files and match their frozen byte counts;
2. validate the archive layouts and extract only the six PRO18-selected
   restart/state identities;
3. replay each source’s documented legacy journal grammar, order, terminal
   identity, and checkpoint-embedded-history linkage;
4. construct compact, acyclic AUTH1 input evidence for the later distinct
   `FGC-1-PRO18-AUTH1` stage. PREF26 does not yet construct the final PROTO17
   `GenesisSpec`, because the future production runtime and runner source pins
   belong to AUTH1.

It does **not** claim to re-execute all historical PDE steps. The records were
authenticated and semantically checked as legacy evidence, not rerun as new
physical trajectories.

The result safely retained and hashed one byte image of every source before
parsing it, rejected path and same-inode mutation races, validated all `53` NPZ
members and `28,083,115` uncompressed bytes, reconstructed the six selected
physical/restart identities, and replayed `49 + 25 = 74` journal records. The
selected PROTO12 event and RSP2 endpoint were cross-linked to the recomputed
state hashes and checkpoint counters. The focused adversarial suite additionally
rejects unsafe paths, hostile ZIP flags and NPY versions, nonfinite metadata,
source/resource drift, result tampering, and authority or claim promotion.

## Temporal and authority boundary

PREF26 must not inspect, list, create, or otherwise access either future
`PROTO17` output root. The absent-root observation belongs only to the PRO18
freeze boundary and cannot be repeated as though it were current evidence.
`FGC-1-HLT15-GEN1` instead owns the later same-invocation authority/import
recheck immediately before atomic generation zero.

The complete successor order remains:

```text
FGC-1-PRO18-PREF26  actual read-only source evidence
  -> FGC-1-PRO18-AUTH1  separately committed authority result
  -> FGC-1-HLT15-GEN1  isolated authority recheck and atomic genesis
  -> FGC-1-PRO18-PREF27  actual reuse and foreign-store rejection
```

PREF26 does not authorize `AUTH1`, launch, genesis, or any later stage. Its
successful result authorizes only the separate implementation/design work needed
to construct and test the future AUTH1 certificate. A
missing, malformed, oversized, changed, or semantically inconsistent source is
an invalid binder premise, not evidence about FGC-QR or a reason to repair the
historical record.

## Duty transition and nonclaims

The completed PREF26 execution discharges exactly two HLT14 duties: actual
raw-bundle byte/hash recomputation and semantic historical-journal payload
replay. Live external Git/blob/import authority
remains incomplete until `FGC-1-PRO18-AUTH1`; production namespace reuse and
foreign-store rejection remain incomplete until `FGC-1-PRO18-PREF27`.

The canonical result earns `PREF26_completed=true`,
`actual_raw_bundle_byte_and_hash_recomputed=true`,
`actual_legacy_history_semantically_replayed=true`, and
`AUTH1_input_evidence_derived=true`. AUTH1 itself, the final PROTO17
`GenesisSpec`, HLT15 genesis,
PREF27, output-root observation or creation,
pretrajectory operations, fresh GR-0 calibration, SGB-L, FGC-QR, DEF1,
retained-EFT promotion, transition, singularity resolution, child topology,
dark-sector, varying-`c`, and every physical claim remain false.
