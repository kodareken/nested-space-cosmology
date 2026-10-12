# FGC-1-TDG7-IMP3: synthetic stage-safe runtime-repair qualification

`FGC-1-TDG7-IMP3` is the implementation and synthetic qualification layer
authorized by the sealed `FGC-1-TDG7-PREF23` binder at commit
`e9289ce479a2cd720d6c58f2ca3d653f573d8eec`. It binds the exact PREF23
artifact and independent-theorem bytes from both that commit and the worktree,
then binds the live TDG7 stage-safe runtime and its direct core test by
SHA-256. A byte change to any bound input invalidates reproduction.

## What executes

The reproducer constructs a small deterministic synthetic six-field state and
uses the existing TDG6 transaction, normal-flow tracer, and temporal-ledger
interfaces. It executes both real method paths through the TDG7 adapter:

| Method | 1/2/4 shadow proposals | Total stage records | Fine stage records |
| --- | ---: | ---: | ---: |
| RK4 | 7 | 35 | 20 |
| SSPRK3 | 7 | 28 | 16 |

For each method, one immutable TDG7 binary64 plan supplies the exact
coordinates for the one-, two-, and four-substep paths. The result records the
five plan boundaries, the three geometric stage coordinates for every path,
and every actual RK4/SSPRK3 stage name/time as binary64 hexadecimal. The
adapter revalidates that record before TDG6 admission or commit. Only the
complete fine four-quarter path can atomically become the accepted state;
outer and medium paths remain evidence.

The controls also exercise every typed coordinate-lattice stop with poison
objects in all runtime-hook positions. Each stop is required to occur before
accessing a shadow, source, projector, tracer, monitor, state, coordinate
array, or temporal ledger. The certificate records zero hook accesses for:

- non-finite coordinates;
- a non-forward target;
- an out-of-envelope request;
- a non-`Q`-aligned endpoint;
- a target not aligned to the `8Q` stage lattice;
- a cap with no positive aligned width; and
- a width below TDG6's frozen minimum.

## TDG6 is inherited, not rewritten

TDG7 owns coordinate selection and shared-plan validation. TDG6 remains the
sole owner of all 18-channel continuous admission, temporal retry evidence,
and the fine-path atomic transaction. Synthetic interval controls retain:

- an all-18 exact-zero pass;
- one resolved-order-failure channel vetoing the complete all-of admission;
- a large convergent pass retaining an outward debit greater than `9` in every
  fine component; and
- no absolute, physical-signal, or outcome-scaled magnitude tolerance.

Late synthetic failures at a fine-only stage confirm distinct ownership. A
source-only failure remains owned by PROTO7; a non-source failure is terminal
evidence rather than a source or temporal retry. Both preserve the real state,
monitor, causal ledger, tracer, and TDG6 ledger. A deliberately failing real
tracer commit confirms rollback of monitor, causal, and complete tracer
boundaries.

A failed 18-channel admission is first written to a synthetic durable sink,
then replanned through TDG7 using TDG6's exact half-cap. The resulting
immutable `TDG7RetrySuccessor` binds the predecessor preparation, smaller
shared plan, complete TDG6 retry evidence, and its exact updated ledger.
`prepare_tdg7_retry_successor` is the only tested resumption route and consumes
that exact ledger. The executed fresh-entrypoint control passes it an **active
retry ledger**, which it rejects. This is not a claim that the stateless adapter
can detect a discarded external retry receipt from later identical zero-ledger
inputs. Separately, the retry-successor control rejects a mismatched method,
over-limit count, forged serialized rejection, or forged full rejection prefix
before another shadow path is made. The successor check retains both the exact
full durable-rejection prefix and TDG6's
`last_accepted_macro_step_temporal_retry_count` field.

The normal minimum-width control starts at `2 * TDG6_MINIMUM_MACRO_STEP`,
replans once to exactly `TDG6_MINIMUM_MACRO_STEP`, and then receives TDG6's
typed, nonphysical `minimum_macro_step` exhaustion on the next failed
admission. TDG7's `TDG7RetryLatticeExhausted` is not the ordinary minimum-step
result: it is exercised only defensively with an impossible forged
`TDG6TemporalRetryRequired` wrapper around TDG6's exhaustion evidence. The
runtime also freezes its minimum-width lock so callers cannot override it.

The adapter intentionally has no fresh-versus-retry protocol discriminator:
`fresh_vs_retry_protocol_discriminator_implemented` is false. A stateless
adapter cannot distinguish identical rewound inputs after an external retry
receipt has been discarded. A future protocol cursor must own that lifecycle
decision before any calibration can be considered. **Only prospective design
of a replacement protocol freeze with that durable fresh-versus-retry cursor is
authorized.** No replacement protocol is frozen, and this artifact creates no
protocol file, runtime, run, output namespace, calibration, or candidate.

## Scope boundary

This layer does not load CAL11 or PROTO14 history arrays, does not load or
resume a campaign checkpoint, does not create a runs namespace, and executes
no production trajectory. It preserves PREF23's conclusion that CAL11's
terminal result is an invalid runtime outcome, not physics.

The following remain false: historical PROTO14 mutation or resume; any new
protocol freeze, protocol implementation, or output namespace; calibration
completion or eligibility; SGB-L,
FGC-QR, or DEF1 execution; retained-EFT evolution; and every physical,
transition, singularity-resolution, topology, dark-sector, or varying-local-
`c` claim.

## Reproduction

```bash
python3 scripts/reproduce_fgc_tdg7_imp3.py --verify
python3 -m unittest tests/test_fgc_tdg7_imp3_reproduction.py -v
```

The first command rebuilds from `build(load_config())` and requires the
tracked JSON bytes to match exactly. The second includes adversarial lineage,
hash, stage-plan, active-retry-ledger, retry-lineage, over-limit, and nonpromotion
attacks.
