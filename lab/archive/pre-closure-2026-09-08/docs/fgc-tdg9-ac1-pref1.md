# FGC-1-TDG9-AC1-PREF1: independent AC1 terminal binder

`FGC-1-TDG9-AC1-PREF1` independently binds the completed two-leaf AC1
production-classifier terminal produced under immutable authority commit
`ec4195d5f09788bad3e08936b9506716a9eafa3d`, whose single parent is
`660f369365ab8417f4c68eded46acb9bb2b05a5b`. The one-time construction path
uses no-following, race-checked reads for the canonical raw manifest and
terminal, verifies their exact SHA-256 identities, authenticates the authority
lineage and six committed AC1 blobs, and checks the scientific campaign store
before and after against its sealed 115-leaf identity.

The binder imports no AC1 runner serialize, validate, reduce, run, or status
code. It does not import TI2 shadow restore or prepare helpers and does not
reconstruct the seven SSPRK3 proposals. For the completed terminal it
independently decodes every exact binary64/rational interval, reclassifies all
eighteen channels through design-only `classify_tdg6_channel`, enforces
`TDG6_COMPLETE_STATE_CHANNELS` order, requires debit equal to the finest upper
bound, requires `failed_channels` if and only if the admission flags fail, and
requires the raw terminal class if and only if that independent reduction.

The bounded outcome is exact. Seventeen of eighteen channels admit. The only
nonadmitted channel is `u:R`, with class `order_inconclusive`. This confirms
that one frozen SSPRK3-tableau-on-inherited-SBP4 state and width has one
unresolved `u:R` channel. It licenses only a separately frozen diagnosis of
`u:R` or a sharper discriminator. It does not earn a production method,
selected remedy, state advance, common event, calibration, candidate,
mechanism, EFT, transition, or physics.

Ordinary verification reads only the tracked compact certificate. The raw
namespace and campaign store are opened only by the explicit one-time live
construction command.

```bash
make fgc-tdg9-ac1-pref1
make verify-fgc-tdg9-ac1-pref1

# Explicit one-time construction only; not an ordinary verifier:
python3 scripts/reproduce_fgc_tdg9_ac1_pref1.py --live
```
