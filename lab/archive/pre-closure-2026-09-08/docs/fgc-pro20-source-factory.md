# FGC-1-PRO20-SRCF1 — physical-source origin rebase implementation

This is implementation code awaiting prospective source/environment/resource
qualification. It is not a compact result, live runner, campaign authority or
physical result.

`evolution/pro20_source_factory.py` consumes a byte-revalidated HLT15 capture,
the four exact static input byte strings and an external source-closure
reference. It then:

1. calls the preserved `build_static_gr0_shells`, including its time-zero
   physical-source validation;
2. restores each retained generation-one HLT16 descriptor/payload into the
   corresponding fresh historical-type shell;
3. wraps the supplied Proto12 source/projector with SRCB1 configuration
   binding;
4. seeds a fresh IMP1/C1-reference-wire ledger without resetting physical
   step/serial/source/CFL counters or the inherited TDG6 sibling;
5. constructs the new HLT17 live member and an explicit PROTO17-to-PROTO19
   rebase projection using the real descriptor and full old cursor hashes;
6. creates a generation-zero in-memory HLT17 checkpoint at `t=23/16` with
   sole target `t=3/2`.

The descriptor address is deliberately distinct from the origin-receipt hash.
Matching the supplied source-closure digest does not authenticate it: HLT17's
external reference and every local source-binding authority flag remain false.

Ordinary tests exercise only refusal and nonexecution boundaries. The success
path must be exercised exactly once under a separately committed SRCQ1 freeze
that pins implementation bytes, exact environment, input/capture identities,
all six members, resource limits and no-write/no-state-advance controls. That
qualification result must be independently reproduced before PRO20 authority.
