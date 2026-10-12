# FGC-1-PRO19-PREF28 — generation-nine invalid-terminal binder

PREF28 is the outcome-neutral, read-only binder for the first SID3-authorized
attempt to continue the bounded GR-0 common event. The attempt did not complete
the event. It appended one fail-closed invalid terminal record and installed a
matching terminal lock while preserving every accepted member state from
generation eight.

The exact durable edge is:

```text
generation 8 / sequence 8
  checkpoint 6d2b6e37...
  nonterminal, event 23, target 3/2
        |
        | one terminal-lock record; no accepted-state record
        v
generation 9 / sequence 9
  checkpoint 770728fa...
  invalid_terminal, event 23, target 3/2
  terminal record 9e8b8a58...
  terminal lock 1726b50a...
```

The persisted terminal evidence says only that `RK4-2049` stopped under the
generic invalid owner with exception type `HLT16AttemptError`. It contains no
new TDG6 rejection record, no source or CFL exhaustion, and no scientific
stop. The exception message was not serialized. PREF28 therefore does not
infer or announce a root cause from the terminal record.

## What is independently bound

The binder cross-checks immutable SID3 authority commit
`b8fea44384be62918670ccbab5bd1a2d92057224`, its implementation parent, and
its exact authority configuration, execution closure, result, HLT16 compact
configuration/result, and consumed launch-manifest bytes. The last three are
historical inputs after terminal closure; they are not refreshed against later
repository code. The binder then traverses the terminal campaign without
following symlinks and verifies:

- the exact 56-leaf terminal store;
- the exact four-leaf append-only suffix over SID2's certified 52-leaf
  generation-eight store;
- the canonical and raw identities of generation eight, generation nine,
  sequence nine, and the terminal lock;
- the generation-eight parent and sequence-eight journal parent links;
- the complete terminal schema and `invalid_terminal` disposition;
- absence of an active writer, staging residue, orphan, fork, symlink, or
  special file;
- exact equality of the generation-eight and generation-nine member maps.

The four new leaves are one checkpoint, one journal record, one retired writer
lease, and one terminal lock. There is no new state descriptor or payload.
All six accepted descriptors, cursors, ledgers, counters, and accepted times
are unchanged. In particular, `RK4-2049` remains at finite time
`0x1.78554de5a30e0p+0`; the other five members remain at `23/16`.

## Interpretation

This is a trustworthy invalid-instrument terminal, not a scientific negative.
It proves that the runtime failed closed and that its last accepted finite
state is recoverable as evidence. It does not prove why the in-memory adapter
raised `HLT16AttemptError`; that question belongs to a separately frozen,
read-only diagnostic.

The locked campaign cannot be resumed, rolled back, reclassified, or repaired
in place. Any repair and any successor campaign require new prospective
authority and a new namespace. The unchanged threshold, retry, health,
constraint, and physical rules remain inherited; PREF28 weakens none of them.

No common event, GR-0 calibration, trapped interval, SGB-L or FGC-QR
trajectory, activation, affine-null defocusing result, retained-EFT result,
transition, or physical claim is contained in this artifact.
