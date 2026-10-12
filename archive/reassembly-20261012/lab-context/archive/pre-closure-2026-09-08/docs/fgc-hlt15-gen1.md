# FGC-1-HLT15-GEN1: authenticated calibration generation zero

**Status:** generation zero has been materialized in the calibration namespace;
this is not a numerical evolution result.

HLT15 consumed an externally supplied authority tuple for committed
**FGC-1-PRO18-AUTH1**, then rechecked the pinned authority bytes, production
source closure and import origins, frozen evolution environment, and actual
PROTO12/RSP2 source reimport in the same launch invocation. Only after those
checks did it make the first observation of the calibration output root and
atomically materialize one `FGC-2-SF1-PROTO17` generation-zero store.

The resulting store has exactly six restored state objects, one generation-zero
checkpoint, and one receipt. The receipt identifies generation `0`, records the
external AUTH1 tuple, and binds the checkpoint and state-object set. Each member
is `FRESH_READY` at the inherited common event `23` and coordinate time
`23/16`.

## Exact nonclaim boundary

Generation zero installs authorized starting state; it does not execute a
macro-step. The store is nonterminal: it has no accepted new macro-step, no
new common event, no retry, no journal record, no debit, and no terminal
classification. It does not establish calibration completion or eligibility,
authorize candidate execution, open SGB-L or FGC-QR, establish a DEF1
observable, promote retained-EFT validity, or support a transition or a
physical claim.

The sealed `FGC-1-PRO18-PREF27` record authenticates the real eight-leaf store
and records refusal of namespace reuse and
foreign-store adoption, including stale
sibling HLT15 staging-residue rejection and three explicit no-residue control
receipts. That later binder may close the final HLT14 production-duty category,
but it does not change this generation-zero/non-evolution boundary.
