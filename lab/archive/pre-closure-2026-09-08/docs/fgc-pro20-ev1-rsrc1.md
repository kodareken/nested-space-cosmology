# FGC-1-PRO20-EV1-RSRC1 — resource-layer design only

This document owns a prospective resource-isolation implementation. It is not
a freeze, not a live target, and not an authority.

[PREF1](fgc-pro20-ev1-pref1.md) independently bound the one authorized
PRO20 first event as a typed `resource_exhausted` stop with partial
accepted state and no common event. SSPRK3-8193 stopped because
`max RSS exceeded`. That result does not license resume, takeover, a new
namespace, another production run, calibration, SGB-L, DEF1, holdout,
candidate, mechanism, or physics.

## Status

The isolated implementation core is present. It introduces a fresh child for
each seed member and each scheduled member attempt, a one-member physical
reconstruction, bounded canonical handoffs, a complete six-member in-memory
seed-cohort gate, and parent-side byte reducers. A separately named mechanical
schema sibling, atomic store, and parent-only scheduler/publisher are now also
implemented under the fresh prospective namespace
`runs/fgc-2-sf1/pro20-rsrc1-event1/event`. They cannot parse, open, append, or
take over the closed `pro20-event1` store. No authority module or live target
exists, so `may_execute_state_change` remains false.

The new member identity is
`FGC-2-SF1-PRO20-RSRC1-EVENT1`; it cannot alias the closed PRO20 campaign.
The physical protocol, action, amplitude, HLT15 origin, C1R1/IMP1 arithmetic,
CFL formula, retry caps, event target, and accepted-state semantics are
unchanged.

## Scope

RSRC1 fixes the following operational isolation policy:

- one child process constructs exactly one scheduled member;
- six independently returned seed bundles must form the exact canonical cohort
  before a future parent can publish anything; a resource-stopped or missing
  seed leaves no namespace;
- the child peak-RSS ceiling remains 4 GiB, rather than fitting a larger
  ceiling to the failed run;
- the store-owning parent is limited to 1 GiB current RSS;
- 4 GiB remains unallocated as host reserve, so the pre-child available-memory
  premise is exactly 9 GiB;
- one child receives the inherited 86,400-second wall allowance;
- each canonical request or response is bounded by 256 MiB and stderr by
  1 MiB;
- the child imports no PRO20 production store, reads no production namespace,
  and writes no endpoint or campaign state;
- the parent reauthenticates every descriptor, payload and cursor and derives
  the protocol kind from predecessor/successor bytes rather than accepting a
  child label;
- an explicit clean child resource response can only return the unchanged
  accepted bundle; a signal, timeout, nonzero exit, stderr, malformed response,
  or response-hash mismatch is forensic uncertainty and cannot be published as
  a typed resource or scientific terminal.
- only the parent imports the fresh RSRC1 store; it publishes the whole seed
  container atomically, owns the exclusive writer/session, authenticates the
  complete chain before every schedule decision, and appends one validated
  generation at a time;
- the new root, journal, checkpoint, generation, session, terminal-lock, and
  exclusion schemas bind `FGC-1-PRO20-EV1-RSRC1-FRZ1` prospectively and use
  store kind `production_rsrc1_event1`; constructing them in tests is not
  authority.

RSRC1 may not:

- resume or take over the closed `pro20-event1` store;
- mutate a terminal generation, journal, checkpoint, or lock;
- authorize a new production namespace;
- rerun or continue the consumed FRZ1 event;
- open PRO21, Wave 2, calibration, SGB/DEF, holdout, candidate,
  mechanism, or physics.

## Next freeze

Before a freeze, the remaining owners are:

1. independent real-origin parity for all six one-member seed bundles;
2. current/available-memory and largest-grid physical child qualification,
   exact environment and process
   preflight, and synthetic crash/kill/partial-publication controls;
3. the exact-delta `RSRC1-FRZ1` authority and one-shot status/run routing;
4. an independent binder seam that does not trust child or parent labels.

Any new measurement still needs a separately committed
`FGC-1-PRO20-EV1-RSRC1-FRZ1`. Until that freeze exists, no status or run target
is authorized. Focused component tests are read-only:

```bash
python3.14 -B -m unittest \
  tests/test_fgc_pro20_rsrc1_member.py \
  tests/test_fgc_pro20_rsrc1_attempt.py \
  tests/test_fgc_pro20_rsrc1_isolation.py \
  tests/test_fgc_pro20_rsrc1_seed.py
```
