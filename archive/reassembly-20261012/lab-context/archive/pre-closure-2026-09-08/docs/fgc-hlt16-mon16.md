# FGC-1-HLT16-MON16: durable-runtime implementation inventory

**Status:** sealed compact Phase-2 runtime qualification. It documents and
hash-binds the implementation boundary; it is not a run result, calibration,
or physical authorization.

MON16 binds the frozen PRO19 first-edge contract and the sealed PREF27
generation-zero binder, then inventories the HLT16 sources, tests, status
owner, runner owner, and manifest by repository-relative path and SHA-256.
Every inventory hash is sealed only after the implementation, tests, status,
runner, and documentation stop moving. The independent PRO19 launch manifest
then binds the complete committed image; MON16 itself cannot grant execution.

The pinned AUTH1 compatibility contract is Python `3.14.3`, NumPy `2.5.1`,
Darwin, and `arm64`. MON16 records that contract as provenance; it does not
import NumPy, inspect arrays, open the live run store, construct a state, or
run a member.

The qualified ownership boundary is explicit. Every post-generation-zero
payload, journal, checkpoint, or terminal-lock mutation requires a live
in-process writer capability that is revalidated against the on-disk lease
and the authenticated checkpoint lineage. Released or superseded leases are
atomically moved to content-addressed, independently validated retired records
and are never deleted by pathname. Historical raw arrays may be opened only
at the exact authenticated generation-zero boundary. The public low-level
state store rejects evolved-state writes; those writes exist only behind the
campaign capability boundary. Before the first checkpoint, the complete six-
member payload/descriptor/journal/checkpoint bridge is derived in memory and
every existing final leaf must be an exact prefix of those bytes before replay
may mutate the store. A two-link publication stage is replayable only when its
final hard link and exact bytes prove the inode pair. A one-link prepublication
stage is retained and rejected as ownership-ambiguous: a deterministic name is
not an ownership credential. In that early interruption case, recovery uses a
fresh verified copy of the immutable generation-zero authority while the
ambiguous tree is preserved as evidence. Every published HLT16 checkpoint then
makes subsequent recovery descriptor/payload-authoritative.
The production command is pinned to this repository, the canonical launch
manifest, and the authenticated calibration store.

The local integrity threat model is explicit. `bootstrap.guard` and the later
checkpoint-bound writer lease serialize every authorized project writer, and
the adversarial qualification covers crashes, stale or foreign pre-existing
trees, unsafe filesystem nodes, forks, and byte mutations. It does not claim
to withstand an actively malicious, non-cooperating process running as the
same macOS account while publication is in progress; such a process can also
rewrite the repository and every local evidence file. That stronger claim
would require an operating-system or external signing trust root not supplied
by this local campaign. Any concurrent suffix that is observed by the normal
validators remains an invalid store, never an alternative history.

## Exact claims

The artifact may state only:

- durable runtime implementation exists;
- the implementation inventory and copied-generation-zero qualification are
  sealed; and
- an external committed-image authority may assess exactly one first-event
  preflight without treating MON16 as the authority itself.

It explicitly does **not** authorize a one-event preflight, production state
advance, GR-0 calibration completion, candidate execution through SGB-L or
FGC-QR, or any physical claim. The result remains a finite, source-text
provenance record for **FGC-2-SF1-PROTO18**.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt16_mon16.py --verify
```

Use `--write` only while producing the immutable integration checkpoint after
the final inventory refresh. The ordinary verifier consumes the compact result
only; it never treats MON16 as permission to read, create, or advance the
calibration namespace.
