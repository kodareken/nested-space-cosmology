# FGC-1-TDG11-C1R1 — exact C1 representation optimization (prospective)

This owner is a prospective performance implementation of the selected C1
mathematical object. It changes neither that object nor a scientific gate.
The sealed TDG11-IMP1 implementation remains unchanged and is the reference.
No C1R1 production use, full-size measurement or new trajectory is authorized
by this document.

Small 9/33/129-node profiles identify exact Fraction construction,
normalization, five-point Bernstein evaluation and coefficient hashing as the
dominant synthetic cost. The provisional full-size timings are not a formal
qualification and do not establish physical-RHS or worst-case fallback cost.

## Fixed equivalence contract

Preserve all original RK4/SSPRK3 stages and endpoints, exact accumulation
reconstruction, all eighteen complete-state channels, p>=3/2, zero/inconclusive
cases, dual-localizer fallback, resource-stop behavior and full non-cancelling
public debit. Only the original binary64 fine endpoint is committable.

Binary64-derived direct controls lie in the dyadic module `(1/3)Z[1/2]`.
That set is not a ring closed under arbitrary products: `(1/3)*(1/3)=1/9`.
The C1R1 file and helper names keep the historical "ring" identifier. An
integer-plus-binary-exponent representation can evaluate their dyadic linear
operations without repeated Fraction gcd work. Specialized exact evaluations
at 0,1/4,1/2,3/4,1 and dyadic de Casteljau restrictions must preserve the
reference mathematical values, heap/tie decisions, coefficient/row hashes
and public wire encoding. Inputs outside a supported fast domain use the
unchanged reference path on the original logical stream; they are not
silently discarded, misclassified or given a relaxed limit. One-pass
iterators are replayed from a prefix of at most `expected_row_count + 1`
rows, matching IMP1's extra-row detection budget.

The implementation map is new `evolution/tdg11_c1r1_ring.py` and
`evolution/tdg11_c1r1_enclosure.py`, plus focused equivalence/adversarial tests.
Runtime integration, if needed, lives in a new `tdg11_c1r1_runtime.py` adapter.
There is no monkey-patch of IMP1. Any use of its receipt/ledger formats is an
explicit reference-wire compatibility claim, not a claim that the old
implementation ran. Production authority must name the actual C1R1 owner and
implementation hashes. No swapped method or corrected endpoint is permitted.

## Implemented APIs

Mathematical owner: MSEL1-PREF1 / IMP1 C1 object
`exact_accumulation_reconstruction_with_debit`.

Implementation owner:

- `C1R1_IMPLEMENTATION_ID = tdg11_c1r1_integer_exponent_direct_ring_v1`
- `ARTIFACT_ID = FGC-1-TDG11-C1R1`

Reference-wire format: IMP1 enclosure dataclasses, IMP1 hash domains,
`tdg11_imp1_exact_bernstein_with_dual_fallback_v1`, and the IMP1 ledger
rejection event type.

| Module | Public API | Role |
|---|---|---|
| `tdg11_c1r1_ring.py` | `RingElement`, `from_input`, `hermite_bernstein`, `split_bernstein`, `bernstein_samples`, `sample_lower`, `hull_upper` | Immutable integer/exponent representatives of the dyadic module `(1/3)Z[1/2]`; `exp3=2` only for squared witnesses in `(1/9)Z[1/2]`. |
| `tdg11_c1r1_enclosure.py` | `bound_c1r1_cubics`, `assess_c1r1_rows`, `assess_c1r1_channel`, `assess_c1r1_family` | Certify-once Bernstein enclosure returning IMP1-wire assessments. |
| `tdg11_c1r1_runtime.py` | `prepare_c1r1_initial_runtime`, `require_c1r1_admission`, `replan_c1r1_retry`, `prepare_c1r1_retry_runtime`, `commit_c1r1_runtime` | C1R1-named adapter. Assesses with `assess_c1r1_family`; does not call sealed IMP1 `_prepare` and replace the answer. |

Non-ring exact rationals (odd denominator primes other than 3, or a 3 in a
dyadic-only slot) take `assess_imp1_rows` / `bound_exact_cubics` unchanged.
Fallback classification still uses the sealed
`assess_rational_complete_c_rows` and `classify_tdg6_channel` instruments.
Localizer intersections may leave `(1/3)Z[1/2]`; those gate/extra values stay
ordinary Fractions, matching IMP1. Public debit and accumulation remain in
the direct ring.

## Private compatibility dependencies

These are named so they cannot be mistaken for C1R1-owned implementations:

- `_imp1_reference_assess_rows` / `_channel` / `_family` / `_bound_exact_cubics`
  — unsupported-fast-domain fallback only.
- `_imp1_reference_shadow_path`, `_clone_tracers`, `_restore_tracers`,
  `_tracer_snapshot`, `_validate_shadow_paths`, `_check_current_boundary`,
  `_plan`, `_wire`, `_digest`, `_coordinates`, `_state_hash` — sealed
  shadow/guard/lattice/receipt helpers.
- `accept_imp1_step`, `append_imp1_rejection`, `TDG11IMP1TemporalRejection`
  — IMP1 ledger wire, including `rejected_TDG11_IMP1_temporal_admission`.
- IMP1 Bernstein hash domains, used so coefficient/row/assessment SHA-256
  values equal the reference wire.

C1R1 prepared receipts include both `implementation_id` and
`reference_wire_artifact_id`. Semantic receipt equivalence is equality of
family hash, IMP1-wire assessment mapping, shadow endpoints and public
debit, not equality of the C1R1 preparation SHA-256. Retaining the IMP1
rejection event type is explicit ledger-wire compatibility; C1R1 receipts
must still bind `C1R1_IMPLEMENTATION_ID` so later source authority names
the implementation that ran. Private IMP1 guard helper messages may remain
attributed compatibility text.

## Proof before any replacement

Require exact equality to IMP1 for bounds, classifications, debit, sample and
polynomial accounting, hashes and semantic receipts on its full focused
enclosure corpus and both methods' 9/33/129-node families. Include zero rows,
large-base cancellation, signs, dyadic/third denominators, non-ring fallbacks,
interior extrema, threshold straddles, tie ordering, forced fallback, malformed
data and every affected resource limit. Use independent exact point controls
as well as reference comparisons. Do not earn speed by omitting validation or
channels.

Small profiling may guide implementation, but a speedup is not assumed.
Any later production-size measurement needs its own prospective execution
record and coordinator-owned commit/output handling. Permission failure is a
stop, not permission to use alternate Git stores or a different measurement
driver. C1R1 remains unqualified for production until this proof, HLT17
integration and the full source/resource/authority gates are complete.

## Remaining unsupported / unclaimed paths

- Production method selection, HLT17 integration and campaign authority.
- Physical or raw RHS, 2049-point or worst-case fallback measurement.
- Non-ring inputs as a fast path; they are delegated on the original
  logical stream, and IMP1 still rejects them where the reference rejects
  them.
- Using a C1R1 prepared object as an IMP1 prepared object, or swapping a
  C1R1 assessment onto an IMP1 receipt.
- Any claim that identical C1R1/IMP1 assessment hashes mean IMP1 executed.
