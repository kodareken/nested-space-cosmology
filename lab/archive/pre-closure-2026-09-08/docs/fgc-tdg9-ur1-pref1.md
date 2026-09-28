# FGC-1-TDG9-UR1-PREF1: independent UR1 envelope-owner binder

`FGC-1-TDG9-UR1-PREF1` independently binds the completed two-leaf UR1
`u:R` envelope-owner terminal produced under immutable authority commit
`17d8c4cf728ffd18ce6e137874501e48b35fbb32`, whose single parent is
`22f7dc5e181af45eb43f6f19dee23fbc29f89c0c`. The one-time construction path
uses no-following, race-checked reads for the canonical raw manifest and
terminal, verifies their exact SHA-256 identities, authenticates the authority
lineage and all sixteen committed UR1 delta blobs, and checks the scientific
campaign store before and after against its sealed 115-leaf identity.

The binder imports no UR1 or AC1 runner, authority, or serialize code. It does
not import TI2 runner or binder helpers, the runtime envelope implementation,
or LOC1 row reconstruction, and it does not reconstruct the seven SSPRK3
proposals. For the completed terminal it independently decodes every published
`u:R` exact binary64/rational envelope field, reconstructs the lower bound as
`max(0, Fraction(raw) - Fraction(construction debit) - Fraction(arithmetic debit))`
followed by downward binary64 rounding, records Bernstein certification slack
without including it in that clip, names both zero lowers
`coefficient_plus_arithmetic_debit_clips_positive_raw_maximum`, and
reclassifies D01/D12 through design-only `classify_tdg6_channel`.

The bounded outcome is exact. The 2,044-row stream matches
`b4a48c7bf72e015f3f3d9a340ea550d608514fb8cd886562265158c33c0b3969`. D01 is
`[0x0.0p+0, 0x1.2d198e246e459p-38]` and D12 is
`[0x0.0p+0, 0x1.2d1c31bb91376p-38]`. Construction/raw ratios are D01
`654019454794891/25865933815808` (~25.285) and D12
`163504863789045/5012309316608` (~32.621). Both D01/D12 lowers are zero under
the current global aggregation. The current globally aggregated
coefficient-construction debit alone is sufficient on this frozen state/width;
the arithmetic debit is negligible and arithmetic alone is not sufficient.
The production classifier remains `order_inconclusive`. This licenses only
later, separately prospectively frozen envelope or admission design. It does
not earn a production method, selected remedy, state advance, common event,
calibration, candidate, mechanism, EFT, transition, or physics.

UR1 does not prove that raw and construction maxima are co-located, and it
does not prove that candidatewise/local pairing repairs admission. That is a
separate prospective numerical question. FGC-QR remains unopened.

Ordinary verification reads only the tracked compact certificate. The raw
namespace and campaign store were opened only by the explicit one-time live
construction command. The resulting compact certificate is now frozen at
SHA-256 `27c2aef9b033092285b868f583441842de99f61941ff57aa7e62749cc8d91300`.

```bash
python3 scripts/reproduce_fgc_tdg9_ur1_pref1.py

# Explicit one-time construction only; not an ordinary verifier:
python3 scripts/reproduce_fgc_tdg9_ur1_pref1.py --live
```
