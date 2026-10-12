# FGC-1-PRO16-FRZ1: common-event commit and trusted-genesis freeze

**Status:** prospective contract only. PROTO16 corrects one identified gap in
PROTO15: after all six members reach a common target, the successful event
must itself become a durable, recoverable transition. This artifact implements
no runtime, adapter, runner, checkpoint, namespace, or trajectory.

## Frozen successful edge

`COMMON_EVENT_COMMIT` is the only new successful edge. It may occur only when
all six canonical members are `FRESH_READY`, each accepted boundary equals the
active common-event target, and every current TDG6 retry count is zero. The
canonical journal-record payload is strictly predecessor/evidence-only: it
contains the protocol/campaign identity, prior checkpoint and event identity,
completed and next target times, member keys, and prior
cursor/ledger/state-set hashes. It explicitly excludes its own hash and every
successor cursor, set, and checkpoint hash, because those values depend on the
completed record.

The outer `COMMON_EVENT_COMMIT` journal record is written and fsynced first;
its SHA-256 is derived from the canonical envelope and is not embedded in its
own payload. Only after that record hash exists may the six successor cursors
be derived with the record as journal tip. The new checkpoint binds that
derived cursor-set hash while retaining ledger-set and state-set hashes
unchanged. It increments the event index, advances the exact target by `1/16`,
and issues fresh cursor generations and attempt identities; each new cursor
names its preceding cursor as parent.

Only a lone durable `COMMON_EVENT_COMMIT` suffix may later be reconstructed.
Any other uncheckpointed suffix is invalid. Record/checkpoint fault injection
is a later HLT14 runtime obligation, not a result here.

## Trusted genesis

HLT13 proved bounded local persistence integrity but explicitly did not
authenticate wholesale replacement of a self-consistent store. Before any
PROTO16 production namespace is created, future **FGC-1-HLT14-MON14** must be
committed and embed a `GenesisSpec` binding protocol/freeze/HLT13, run plan,
runtime, adapter, runner, environment, campaign, namespace, member descriptor,
physical bundle, checkpoint, journal-root, and cursor/ledger/state-set hashes.
The spec includes the complete canonical six-member descriptor collection;
its aggregate descriptor hash is a checksum, not a substitute for those visible
inputs. This makes generation zero independently derivable without trusting a
value created inside the mutable production namespace. At launch the runner
records a separate authority tuple: authorization commit,
authorization-result path and hash, plus the `GenesisSpec` hash. The result
must not self-reference its own commit or hash; instead the runner proves the
result blob is tracked at the supplied immutable authorization commit. That
future authority is not created by this freeze.

The six restart descriptors are inherited from PROTO15 by its exact compact
hash. All physical, numerical, TDG6, and TDG7 rules are inherited unchanged.
No raw run, historical checkpoint, SGB-L outcome, FGC-QR outcome, calibration,
or mechanism state is read.

## Nonclaims

PROTO16 does not authorize implementation, pretrajectory execution, namespace
creation, GR-0 calibration, eligibility, SGB-L, FGC-QR, DEF1, retained-EFT
promotion, transition, or any physical claim. The only true claims are that
the PROTO16 freeze, common-event-commit contract, and trusted-genesis contract
are frozen.

Reproduce this compact artifact with:

```bash
python3 scripts/reproduce_fgc_pro16_frz1.py --verify
python3 -m unittest tests/test_fgc_protocol_v16.py \
  tests/test_fgc_pro16_frz1_reproduction.py -v
```
