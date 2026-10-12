# FGC-1-PRO18-PREF27: generation-zero store binder and refusal controls

**Status:** sealed read-only binding of the HLT15 calibration generation-zero
store and temporary no-adoption controls. It is not a calibration trajectory.

PREF27 authenticates the real calibration-only `FGC-1-HLT15-GEN1` store against
the committed AUTH1 authority tuple and reconstructed `GenesisSpec`. It checks
the canonical generation-zero receipt, checkpoint, six inherited state objects,
and exact eight-leaf store grammar. The real store is read-only throughout.

The binder separately exercises the following controls in temporary stores:

- refusal of partial, file-shaped, extra-leaf, and symlink namespace reuse;
- refusal to adopt an internally self-consistent store built from a different
  GenesisSpec; and
- an arriving foreign root at exclusive rename; and
- no staging residue after the byte-identical-store, internally consistent
  foreign-store, and arriving-foreign-root controls.

The original real-store tree is rechecked after those temporary controls. This
establishes only that the authenticated generation-zero store is not silently
reused or replaced by the tested malformed or foreign alternatives.

## Portability and verifier boundary

The raw calibration store is ignored runtime evidence. A completely absent
store in a clean clone is permitted: the compact canonical PREF27 result
records its authenticated hashes and control outcome. A partially present store
is an error at the historical seal and fails closed. Ordinary repository
verification checks only the tracked compact result; it does not reopen the
mutable live campaign store. The explicitly named historical live audit may
inspect store names and entry types, and after an authorized successor has
legitimately progressed the campaign it is expected to reject that evolved
store as no longer being the original eight-leaf generation-zero tree. Raw
reimport, temporal root observation, generation-zero materialization, and the
temporary refusal/race controls run only through the explicit focused
historical target.

## Exact nonclaim boundary

PREF27 performs no state advance. It proves no accepted macro-step, common
event, retry, debit, terminal result, calibration completion, GR-0 eligibility,
candidate execution, SGB-L/FGC-QR result, DEF1 result, retained-EFT promotion,
transition, or physical conclusion. The calibration store remains nonterminal
and its six members remain `FRESH_READY`.
