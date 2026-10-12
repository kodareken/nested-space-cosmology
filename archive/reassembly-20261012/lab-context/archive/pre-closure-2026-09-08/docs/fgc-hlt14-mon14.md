# FGC-1-HLT14-MON14: synthetic PROTO17 adapter and runtime qualification

**Status:** completed synthetic-only implementation and qualification. HLT14
implements the production-shaped authority, raw-bundle, and durable-runtime
pathways frozen by PROTO17/PREF25, but exercises them only with temporary test
Git repositories, generated NPZ bundles, and temporary lifecycle stores. It
does not open a real campaign bundle, inspect historical campaign payloads, or
create a production `runs/` namespace.

## Result

The bounded HLT14 corpus contains 72 unique executed cases:

```text
4 nominal construction/recovery cases
6 all-of raw-bundle admissions
18 mutations
44 injected exception windows
--------------------------------------
72 total
```

Their sorted-name digest is
`445765193eddba247b7452635ca382388950a55a3a2bfb128e0b125f530b6baf`.
Every case passed its frozen classification.

This establishes only the following synthetic statement:

> The implemented HLT14 pathway can authenticate a complete temporary
> six-member PROTO17 authority, recompute every generated NPZ and array hash,
> materialize a content-addressed temporary genesis, replay the frozen
> six-record accepted prefix, construct the predecessor-only common-event
> edge, reject the named semantic forgeries, and recover only the declared
> durable images under the tested in-process exception model.

## Two evidence vectors that must not be conflated

HLT14 deliberately runs two different vectors.

The **sealed abstract vector** uses the unchanged PREF25 fixture. The runtime
reproduces its exact byte-level construction hashes:

```text
generation zero  01e8214dcbbef025fe7883670043659d3ec66bb47e9b3e9dd4b8089b8ba5f8de
common receipt   17f925e8489341392ac5df4fba1cc47b66a27153e64eca4455c3f701f37cbc17
successor        c4c31f53a6a97632addafd54d6b3ec08bdcffe4aa4321269f8c8bca957e5380b
```

The **real-NPZ synthetic vector** generates seven finite little-endian
`float64` arrays for each of the six canonical members. It recomputes the raw
archive hash, individual C-order byte hashes, physical-state hash, restart
payload hash, metadata binding, state-object address, and resulting genesis.
Those hashes are correctly input-dependent and therefore are not the three
PREF25 fixture hashes. Reporting the two vectors separately prevents a real
bundle check from being mistaken for reproduction of an abstract descriptor.

## Authority and raw-input boundary

The temporary Git authority test requires a full lowercase commit, one
canonical committed authority blob, a complete normalized `GenesisSpec`, and
exactly three pinned source roles: runtime, adapter, and runner. It compares
the committed blobs with the live files and checks the current Python import
origins.

That import-origin check is a **preflight**. It does not close a later time-of-
check/time-of-use gap and is not authentication of an actual production
launch. The raw adapter similarly exercises production-shaped checks on
generated inputs; it does not validate the real historical campaign bytes.

The all-of bundle gate rejects raw-file drift, array-hash or aggregate-hash
drift, wrong frozen shapes, nonfinite or non-little-endian arrays, metadata
semantic drift, noncanonical or duplicate-key metadata JSON, unsafe ZIP
members, symlinks, malformed manifests, and incomplete six-member admission.

## Semantic and durability boundary

The runtime materializes generation zero through a uniquely allocated sibling
staging directory, fsyncs its files and directories, and atomically adopts the
temporary root. It then admits only the exact synthetic six-record accepted
prefix and the exact PROTO17 common-event construction. Fully rehashed changes
to an executed plan or cursor/high-water state are rejected by semantic replay,
not merely by local hash mismatch.

The fault corpus fixes the expected image for each injected exception window:

- a common receipt interrupted through `before_replace` preserves event 23;
- a replaced common receipt recovers the exact event-24 successor;
- every common-event checkpoint interruption recovers that same successor;
- an accepted record made visible without its complete six-member checkpoint
  is an invalid prefix, never a recoverable partial event;
- a fully replaced accepted-prefix checkpoint opens exactly generation one at
  event 23;
- generation-zero staging faults leave no adopted partial namespace.

These are deterministic **in-process injected exceptions**. Python cleanup
handlers run after the injected exception. This is not evidence of power-loss,
kernel-crash, filesystem, hardware, or production restart durability.

## Mutation boundary

The executed mutations cover foreign or malformed authority, pinned-source
drift, raw and aggregate hashes, frozen array shape, metadata semantic and
canonical encodings, duplicate JSON keys, bundle symlinks, incomplete and
extra accepted suffixes, pre-existing roots, foreign stores, the explicit
production-path guard, and two fully rehashed semantic forgeries. Separate
reproducer tests mutate the actual config and result authorization maps; a
pretrajectory, namespace, candidate, or physical-claim promotion is rejected.

## Four production duties remain open

All four production-instance requirements remain exactly
`required=true, completed=false`:

1. recompute the actual raw production bundle bytes and hashes;
2. validate the actual external Git authority blob and live import at launch;
3. semantically replay the real historical journal payloads;
4. reject reuse or a foreign store at the actual production namespace.

HLT14 therefore authorizes only the design and freeze of a separate production
prelaunch successor. It does not authorize that successor's execution,
pretrajectory operation, production namespace creation, fresh GR-0
calibration, GR-0 eligibility, SGB-L, FGC-QR, DEF1, retained-EFT promotion, a
transition, singularity resolution, a child domain, a dark-sector mechanism,
varying locally measured light speed, or any physical conclusion.

## Immutable inputs

HLT14 binds PREF25 and PROTO17 compact inputs at immutable commit
`4da668fad42d4006d3e09f511eafa0ef0e69d7ed`. Each declared blob must match both
the exact `git show` bytes at that commit and the current worktree bytes. The
canonical result also binds this document, the HLT14 config, reproducer, three
implementation modules, and three focused test modules.
