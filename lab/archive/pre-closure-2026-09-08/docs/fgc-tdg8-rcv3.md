# FGC-1-TDG8-RCV3-FRZ1: generation-nine projection freeze

RCV3 preserves the failed source campaign and starts a separately named
recovery fork from its last valid boundary. FRZ1 is the prospective,
read-only authority for that bootstrap. It does not create a store or run a
trajectory.

The source terminal tree contains 59 regular leaves. FRZ1 hashes every leaf,
including generation 10, journal sequence 11, the terminal lock, four
quarantined writer leases, and the two empty coordination guards. The selected
projection contains exactly 50 leaves: checkpoints 0--9, journals 0--10, all
11 payloads, all 17 state descriptors, and the generation-zero receipt.

All nine excluded leaves have a reason:

- generation 10 and sequence 11 are the invalid terminal edge that the fork
  must not inherit;
- `terminal.lock` enforces that source terminal and belongs only to it;
- four quarantined active-writer locks are retired process leases, not accepted
  state;
- `bootstrap.guard` and `writer.guard` are coordination primitives, not
  append-only lineage.

The projected store must therefore begin with a newly created empty `locks/`
directory. Internal campaign identifiers remain byte-for-byte unchanged so no
cursor, checkpoint, or journal body is rewritten. The external fork identity
is the tuple `FGC-1-TDG8-RCV3-PROJECTION-1`, the fixed destination
`runs/fgc-2-sf1/tdg8-rcv3/calibration`, generation-nine checkpoint
`eb6fddc...`, sequence-ten journal tip `5b533eb...`, and the eventual immutable
FRZ1 commit.

FRZ1 also binds the two execution sources repaired at commit `cde091d...` and
the exact pending temporal retry: depth 2, cap
`0x1.aaa9612df8000p-11`, accepted descriptor `778471...`, and unchanged
physical state `3cd9f5...`.

The only bootstrap image it permits is the projection-only runtime
`tdg8_rcv3_fork_runtime.py` at SHA-256 `9c7e40dd...` and
`bootstrap_fgc_tdg8_rcv3.py` at SHA-256 `0e3bbfc4...`. Before calling that
fixed installer, the script requires an exact committed `HEAD`, no tracked
worktree drift, compact FRZ1 validity, and committed/live equality for the
FRZ1 config, result, runtime, and script. It constructs the external authority
tuple from that immutable commit and those hashes. It imports no evolution or
candidate runtime, and its executable entry point requires `python -I -B`.

The destination-absence observation is intentionally one-time. The compact
verifier later checks the tracked certificate without reopening the live
namespace, so a valid subsequent bootstrap cannot retroactively make FRZ1
fail.

Passing FRZ1 means only that the exact prefix is frozen and bootstrap is
permitted. Bootstrap completion, independent projection authentication,
bounded execution, common-event completion, GR-0 calibration, SGB-L, FGC-QR,
DEF1, and every physical claim remain false.
