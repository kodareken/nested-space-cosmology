# FGC-1-TDG9-TI2-FRZ1: fingerprint-only recovery authority

## Result

**FGC-1-TDG9-TI2-FRZ1** prospectively authorizes one bounded recovery of the
same-state tableau counterfactual that TI1 could not start. TI1 remains the
immutable invalid attempt at commit `8a32636654dbd6479fcb91a88d4557a8090041d3`:
its replay fingerprint rejected the legitimate finite binary64 value
`monitor.last_accepted_time = 0x1.78554de5a30e0p+0` before the first shadow
proposal was constructed. The operator-observed stdout digest is retained with
the explicit limitation that it was not written to a timestamp-authenticated
raw namespace.

TI2 changes only the representation used inside replay fingerprints. A finite
built-in Python `float` is encoded as its exact binary64 identity:

```json
{"binary64_hex":"0x1.78554de5a30e0p+0"}
```

NaN, infinity, NumPy scalars, foreign numeric types, and extended numeric types
remain invalid. Fractions, booleans, integers, strings, nulls, lists, and
mappings retain their prior exact semantics. The general scientific JSON
serializer still rejects floating-point values; this recovery is confined to
the replay fingerprint.

Before authorization, TI2 restores the three exact RK4-2049 retry predecessors
and requires transaction fingerprints:

| Retry | Exact transaction SHA-256 |
|---:|---|
| 3 | `7a90163d2eb4252fa7a1bbdf55d12f92a9127ed7cba5a37c28456e72d732b22b` |
| 4 | `abab435c8222de0a86368e9562a2949578de38fe822c6032319410b49e4c6f6c` |
| 5 | `69fd069034fdd3e2a80a5f3fa989b639b6659105a8ff3e8742e32e78739fd461` |

This prelaunch restoration constructs no compositor proposal and leaves the
sealed 115-leaf campaign store unchanged.

## Scientific question and unchanged scope

If the authority is later committed and explicitly run, it asks only:

> On the exact three restored RK4-2049 retry predecessors already bound by
> AR1/LOC2, does replacing only the RK4 time tableau with SSPRK3 clear the ten
> sealed exact complete-C failures?

The state, accepted time, 2049-point grid, SBP4 right-hand-side operator,
projector, transaction, tracers, temporal ledger, retry widths 3/4/5, ten
occurrences, `p >= 3/2` predicate, exact primary/v2 routes, depths, ceilings,
and work budgets are unchanged. The combined runtime constant is used only to
select the SSPRK3 tableau. The actual spatial operator remains inherited
RK4-2049 SBP4.

This is not production SSPRK3, an independent-method comparison, temporal
admission, a fine-path commit, continuation, calibration, a candidate branch,
a mechanism result, or physics. No outcome selects a remedy.

## Bounded classifications

The diagnostic may later report only:

- `completed_all_ten_complete_failures_clear_under_tableau_shadow`;
- `completed_one_or_more_complete_failures_persist_under_tableau_shadow`;
- `bounded_localization_or_resource_inconclusive`;
- `shadow_proposal_premise_stop`; or
- `invalid_provenance_or_implementation`.

A changed ownership class alone is not clearance. No fourth width, resource
escalation, state advance, campaign write, PDE/candidate execution, or physical
claim is authorized.

## Verification

Store-blind compact verification reads only the tracked certificate and never
opens the live campaign store or restores a predecessor:

```bash
make fgc-tdg9-ti2-frz1
python scripts/check_repo.py --only-tdg9-ti2-frz1
```

The stronger bounded prelaunch proof deliberately performs read-only live-store
snapshot and namespace checks and restores all three predecessor fingerprints,
while constructing no shadow proposal and performing no write:

```bash
make verify-fgc-tdg9-ti2-prelaunch
```

The explicit shadow run target remains separate and is not invoked by either
verification path.
