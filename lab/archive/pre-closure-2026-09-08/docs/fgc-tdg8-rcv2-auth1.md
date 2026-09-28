# FGC-1-TDG8-RCV2-AUTH1: bounded retry-chain authority

RCV2 is a premise-only execution authority for the existing GR-0 TDG8
successor campaign. It replaces the need to freeze a new authority after every
ordinary TDG6 temporal halving.

The authority starts from the authenticated generation-eight checkpoint and
its exact sequence-nine rejection plus sequence-ten cursor transition. That
entry deterministically implies generation nine. After that entry is closed,
RCV2 permits any descendant produced by the original immutable progression
plan, provided the complete HLT16 store validates and remains inside event 23
with target `3/2`, or reaches only the event-24 `25/16` completion boundary or
a typed terminal.

The bounded transition relation is not open ended:

- every TDG6 rejection halves the exact durable plan;
- the retry count is at most 32 for each macro step;
- no successor below `0x1.0000000000000p-30` is executable;
- a rejection preserves the accepted descriptor, accepted time, and physical
  `u,p,q` identity;
- every publication cut must be one already recognized by the HLT16 store;
- a stale writer may be replaced only after same-host dead-PID verification;
- every restart must retain the exact generation-eight ancestor, original plan,
  campaign, and one-event scope.

The compact verifier reads only tracked config/result bytes. The live ignored
campaign is read once while the prospective result is created, and again only
by the explicit status or run command. RCV2 never reimports GEN0, rechecks the
historical destination-absence fact, or opens SGB-L or FGC-QR.

`physical_state_advanced_by_recovery = false` refers specifically to rejection
reconciliation. A later fine step may advance the GR-0 state only after the
unchanged numerical, health, constraint, causal, and persistence gates admit
it under the original plan.

Permitted outputs remain outcome neutral: event completion, a typed scientific
or invalid terminal, or a recoverable interruption. None is a calibration,
trapping, DEF1, mechanism, retained-EFT, transition, or physical result.
