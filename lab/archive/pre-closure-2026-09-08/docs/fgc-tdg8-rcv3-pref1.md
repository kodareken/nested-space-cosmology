# FGC-1-TDG8-RCV3-PREF1: installed projection binder

PREF1 is the independent, read-only post-installation check for the RCV3
generation-nine recovery projection. It authenticates a finite restart
boundary; it neither runs the bootstrap again nor advances a PDE state.

The binder opens the terminal source and installed destination through
non-following descriptor-relative traversal, rejects links and special nodes,
and requires stable repeated snapshots. It verifies the complete 59-leaf
source terminal and the exact 50-leaf destination projection byte for byte.
The nine excluded source leaves remain absent from the destination: generation
10, journal sequence 11, the terminal lock, four quarantined writer leases,
and the two coordination guards. The destination lock directory is empty and
there is no terminal or foreign suffix.

The external bootstrap receipt is checked as two different identities:

- content address `59e92883...`;
- raw canonical-file SHA-256 `ed6b9c7f...`.

Its authority tuple binds immutable installation commit `ec24cd25...`, the
FRZ1 result, the projection runtime, and the bootstrap script. The binder reads
those bytes from Git itself and does not import or call the bootstrap writer or
its decision code.

Inside the projected store, PREF1 reconstructs checkpoints 0--9, journals
0--10, all six legacy states, all 11 evolved descriptors, and all 11 payloads.
Every JSON object is duplicate-key resistant and canonical. Every NPZ archive
has the exact ordered nine-member inventory; each NPY header, shape, dtype,
layout, byte count, finite value set, semantic digest, raw digest, and
descriptor relation is independently checked with `allow_pickle=False`.

The authenticated generation-nine anchor is checkpoint `eb6fddc4...` with
sequence-ten tip `5b533eb7...`. All six member identities are recomputed from
their persisted arrays. `RK4-2049` remains `RETRY_PENDING` under the temporal
owner at exact retry depth 2 and cap `0x1.aaa9612df8000p-11`; the other five
members remain `FRESH_READY`.

Ordinary repository verification consumes only the tracked compact result and
does not open either run store. The explicit live binder is separate and
read-only.

Passing PREF1 means the installed copy is exactly the frozen finite prefix and
is independently restart-authenticated. It does not authorize bounded event
execution, complete a common event or GR-0 calibration, open SGB-L or FGC-QR,
measure activation, trapping, or DEF1, validate a retained EFT, establish a
transition, or earn any physical result.
