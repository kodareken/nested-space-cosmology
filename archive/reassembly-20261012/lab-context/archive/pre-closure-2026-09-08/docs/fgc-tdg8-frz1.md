# FGC-1-TDG8-FRZ1 — one bounded replay-evidence diagnosis

PREF28 sealed the generation-nine invalid terminal without inventing a root
cause that the durable record did not contain. TDG8 now freezes the only next
question: does the exact generation-eight retry replay disagree semantically,
or does native tuple/list representation alone make equal JSON evidence compare
unequal?

The diagnosis is read-only. It may reconstruct the affected `RK4-2049`
generation-eight state and replay evidence in memory, but it may not mutate or
resume the terminal campaign. Three outcomes were fixed before that read:

1. canonical JSON equality with only tuple/list container differences;
2. a scalar, key, length, or channel-order difference;
3. a different replay classification or runtime stop.

Only the first outcome permits one surgical repair: normalize the replayed
mapping through the existing persistence normalizer at the equality boundary.
The durable validator must continue to require JSON lists, all 18 channels in
their frozen order, unchanged scalar values, and the original TDG6 thresholds,
retry counts, debit, and TDG7 lattice.

Any later execution must start a fresh GR-0-only campaign at
`runs/fgc-2-sf1/proto19/calibration` from authenticated generation zero. The
locked `proto17` campaign remains immutable evidence. This freeze neither
applies the repair nor authorizes the new campaign, and it contains no common
event, calibration, candidate, or physical result.

## Frozen-question result

The later read-only `FGC-1-TDG8-DIAG1` discriminator returned the first
predeclared outcome. Raw native/Python equality was false, while canonical JSON
equality was true. The complete type-difference set was exactly:

```text
$.channel_admissions: tuple -> list
$.failed_channels:    tuple -> list
```

There were no scalar, key, length, or channel-order differences. The replayed
ledger and TDG7 successor plan both matched the durable generation-eight
objects, and the locked terminal tree was unchanged before and after the
diagnosis.

The conditional repair therefore normalizes only the recomputed mapping at the
replay equality boundary. A second canonical-byte comparison preserves scalar
type identity in addition to array order and values. Focused regression tests
reject scalar-value drift, integer/float drift, reordered channel admissions,
changed failed channels, and zero or duplicate sink emissions. TDG6/TDG7
classification, thresholds, retry counts, lattice construction, health rules,
and physical rules are unchanged.

This reproduces a lifecycle mechanism capable of causing the observed generic
invalid terminal; it does not rewrite PREF28 into a terminal-message claim.
The repair still authorizes no execution. A separately committed image and a
fresh campaign remain mandatory.
