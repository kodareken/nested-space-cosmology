# FGC-1-TDG8-RCV3-REC1-AUTH1: exact recovery, then same-event continuation

## Purpose

`FGC-1-TDG8-RCV3-REC1-AUTH1` is a prospective, premise-only successor to the
first RCV3 execution authority. It owns one already-durable publication cut:

```text
complete checkpoint: generation 9 / journal sequence 10
uncheckpointed suffix: sequence 11 TDG6 rejection
                       sequence 12 cursor transition
writer:                exact same-host verified-stale lease
predicted closure:     generation 10
event:                 23, target 3/2
branch:                GR-0, amplitude 3
```

The first RCV3 attempt accepted no new physical state. It durably wrote the
third temporal rejection and its retry cursor, then failed while validating
the successor checkpoint. The checkpoint body was not the defect: runtime and
recovery both construct checkpoint
`6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187`.
The in-memory rejection record retained tuple-valued `failed_channels` and
`channel_admissions`, while the validator required the canonical JSON list
representation. Persistence had already normalized the exact records to JSON
arrays, so read-only recovery succeeds. Commit
`ee888156b226c8210a902ed733fd144e0ba24598` closes that representation split at
the journal-schema boundary without changing the action, equations, state,
threshold, retry policy, record bytes, or predicted checkpoint.

REC1 binds that diagnosis as an executable gate. It does not infer physics
from it.

## Exact entry and predicted closure

The compact authority fixes:

- predecessor RCV3 AUTH1 commit
  `afe855be31fe983e940479bffa669e9e60a26344`;
- tuple/list normalization commit
  `ee888156b226c8210a902ed733fd144e0ba24598`;
- generation-nine checkpoint
  `eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56`;
- sequence-eleven rejection
  `341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a`;
- sequence-twelve transition
  `d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4`;
- the complete bounded entry-tree manifest and stale-lease identity;
- predicted generation-ten checkpoint
  `6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187`.

The predicted checkpoint has generation nine as parent, sequence twelve as its
journal tip, event `23`, and target `3/2`. All six accepted descriptors and
physical `u,p,q` identities are unchanged. Only `RK4-2049` receives the
already-durable cursor/ledger transition: retry depth `3` and pending cap
`0x1.aaa9612df8000p-12`. No accepted macro-step or accumulated TDG6 debit is
added by recovery.

## Two-phase runtime rule

The explicit REC1 runner is the only mutation path. It must:

1. authenticate the committed REC1 image and exact live entry;
2. recheck that the bound lease is still same-host and verified stale;
3. take over the lease and reconcile only sequence eleven/twelve;
4. require the resulting complete checkpoint to equal the frozen generation-
   ten identity and invariants;
5. only then hand the still-open event `23` to the existing bounded GR-0
   event runner.

No numerical proposal is permitted before step 4. A mismatch, live or foreign
lease, additional suffix, staging object, orphan state, lost ancestor, or
changed generation-ten body stops invalid. An interruption may resume only
from the exact entry, exact clean generation ten, or a lawful event-23
descendant containing generation ten as an immutable ancestor.

The terminal boundary remains exactly one of:

- event `24` complete with next target `25/16`;
- a typed event-`23` scientific or invalid terminal;
- a recoverable interruption under the same authority.

## Why no intermediate binder is required

Generation nine to ten is a deterministic metadata-only closure, not a new
physical sample. The prospective authority freezes its exact result and the
runner must compare that result before proposing another step. The later
outcome-neutral event binder will independently reconstruct the full
generation-nine-to-ten edge and every later record/state. A separate committed
post-recovery binder would add a pause point but no independent fact.

If either the runtime comparison or the final independent reconstruction were
absent, this compression would not be safe.

## Claim boundary

REC1 does **not** establish:

- a completed common event or GR-0 calibration;
- an activation, trapped interval, or Raychaudhuri result;
- an SGB-L or FGC-QR trajectory;
- a candidate-action success or failure;
- retained-EFT validity, singularity resolution, transition, child domain, or
  any physical mechanism.

It changes no physics and opens no candidate branch. It authorizes only exact
metadata recovery followed by the previously frozen GR-0 event-`23`
continuation.

## Reproduction boundary

Ordinary verification is compact and store-blind:

```bash
make fgc-tdg8-rcv3-rec1-auth1
make verify-fgc-tdg8-rcv3-rec1-prelaunch
```

Read-only live status and explicit execution are separate:

```bash
make status-fgc-tdg8-rcv3-rec1
make run-fgc-tdg8-rcv3-rec1
```

The ordinary verifier never opens the ignored run store, tests lease
liveness, publishes generation ten, or executes a PDE proposal.
