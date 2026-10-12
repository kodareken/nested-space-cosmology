# FGC-1-TDG2-PREF17: absolute-tail history diagnosis

`FGC-1-TDG2-PREF17` binds the outcome-neutral calculation authorized by
FGC-1-TDG2-FRZ1 at immutable commit `a1b6b99…`. It restores the exact CAL10
terminal checkpoint read-only, classifies all six method-owned history ladders,
advances no state, resumes no trajectory, and reads no SGB-L or FGC-QR data.

## What TDG2 actually establishes

The absolute-tail outcome is mixed. Across `576` method/tracer/field ladders,
TDG2 records `1,152` separate field- and derivative-tail classifications:

| Method and measure | Contracts toward zero | Nonzero and resolution-independent within the conditional enclosure | Enclosure dominated |
|---|---:|---:|---:|
| RK4 field | 35 | 94 | 159 |
| RK4 derivative | 35 | 94 | 159 |
| SSPRK3 field | 11 | 91 | 186 |
| SSPRK3 derivative | 11 | 90 | 187 |

No ladder enters `converges_to_nonzero_resolved_power` or
`unresolved_absolute_tail_behavior`. That does **not** mean every ladder is
resolved: the third class explicitly says the frozen instrument cannot
separate the alternatives with the present enclosure.

## Relation to the historical PROTO13 failures

The historical temporal rule has `152` finest-resolution individual failures.
For field-tail power, TDG2 classifies them as:

```text
46  contract toward zero at the required power order
 5  are nonzero and resolution-independent within the conditional enclosure
101 remain floor- or interpolation-enclosure dominated
```

For derivative-tail power the corresponding counts are `46`, `4`, and `102`.
The same result can therefore neither preserve all failures as resolved
physical content nor dismiss all of them as truncation noise.

The dominant unresolved owner is the conditional interpolation debit. Of the
finest legacy failures, `47` RK4 field ladders and `37` SSPRK3 field ladders
have a finest interval reaching zero because of interpolation; derivative
power gives `47` and `38`. A further `14` field/derivative cases lie entirely
below the frozen DFT-amplitude floor, and three RK4 field/derivative cases have
cross-resolution intervals that overlap without directional separation.

The independently audited binary64 DFT arithmetic owns no finest legacy
failure. TDG2 therefore narrows the remaining instrument question from a
generic spectral or floating-point suspicion to the time-grid alignment and
its pointwise enclosure.

## Cross-method evidence

The two methods agree on `259/288` field classifications and `258/288`
derivative classifications. They jointly classify `10` signals as
zero-contracting, `91` field (`90` derivative) signals as nonzero and
resolution-independent, and `158` as enclosure dominated.

Agreement is supporting evidence, not a replacement admission. The methods
have different temporal and spatial discretizations, but they consume the
same frozen compact-history construction and the same interpolation rule.
Shared enclosure domination can therefore reflect a shared diagnostic
limitation.

## Scientific interpretation

TDG2 earns three bounded conclusions:

1. Some PROTO13 individual failures are consistent with a high-frequency tail
   that contracts toward zero at the required rate.
2. A small subset remains nonzero across the present resolution ladders under
   the declared conditional enclosure.
3. Most finest failures cannot yet be classified because the round-trip
   interpolation debit reaches zero or otherwise prevents directional
   separation.

The third conclusion prevents a replacement temporal admission. The next
permitted task is a prospective **TDG3 interpolation-ownership discriminator**
that avoids or independently bounds the common-grid resampling step. It must
be frozen before reading any TDG3 outcome. PROTO14 cannot yet be designed.

## Immutable and clean-clone boundary

PREF17 verifies the TDG2 freeze commit and every frozen tracked blob, hashes
the four-file CAL10 raw bundle, restores the terminal histories through the
independent CAL10 loader, recomputes every ladder with the frozen TDG2 module,
and attacks history, method-ladder, lineage, count, interpolation, replacement-
gate, and eligibility mutations.

When the ignored raw bundle is present, verification reproduces the complete
diagnosis. A clean clone may verify the compact canonical record and its
tracked hashes without that bundle. A partially present raw bundle fails
closed.

## Claim boundary

```text
TDG2_actual_histories_diagnosed = true
TDG2_absolute_tail_outcome_is_mixed = true
individual_temporal_budget_failure_cause_fully_derived = false
replacement_temporal_admission_defined = false
TDG3_interpolation_ownership_discriminator_design_may_begin = true
TDG3_frozen = false
PROTO14_frozen = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

Reproduce or verify the result with:

```bash
make fgc-tdg2-pref17
```
