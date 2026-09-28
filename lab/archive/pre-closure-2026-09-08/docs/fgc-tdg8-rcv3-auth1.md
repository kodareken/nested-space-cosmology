# FGC-1-TDG8-RCV3-AUTH1

## Result

AUTH1 prospectively authorizes one bounded GR-0 continuation on the installed
RCV3 projection. It consumes immutable PREF1 commit `46abf807...`, binds its
exact compact result, the external projection receipt, and the unchanged
internal HLT16 campaign identity at generation `9` / journal sequence `10`.

The permitted transition is only:

```text
event 23, target 3/2
        -> existing source/CFL/TDG6/TDG7/scientific-stop relation
        -> event 24 boundary, next target 25/16
```

It preserves the `RK4-2049` depth-two temporal retry, pending cap
`0x1.aaa9612df8000p-11`, retry ceiling `32`, and minimum macro-step
`0x1p-30`. The original internal execution commit, plan, campaign ID, physical
inputs, thresholds, and numerical engines are unchanged.

## Verification boundary

Ordinary verification is compact and store-blind. Immediately before status
or mutation, the runner independently reauthenticates the external receipt,
generation-nine ancestor, current lineage, persisted payloads, retry state,
one-event scope, and writer/restart status. The status path is read-only.
The compact implementation inventory explicitly hash-binds both predecessor
authority modules whose frozen-policy and repository helpers AUTH1 executes.

The runner can only restore persisted state, construct the existing static
GR-0 shells, and call the existing bounded one-event engine. It has no
generation-zero import, bootstrap, SGB-L, FGC-QR, calibration promotion, or
physical-result path.

## Nonclaims

AUTH1 has not executed a PDE proposal or completed a common event. It earns no
trapping calibration, candidate, mechanism, retained-EFT, transition, or
physical result. Retry exhaustion, provenance rejection, or runtime failure
remains `invalid_implementation_or_nonconverged_run`, never candidate physics.

## Commands

```bash
make fgc-tdg8-rcv3-auth1
make verify-fgc-tdg8-rcv3-event-prelaunch

# Only after the committed AUTH1 SHA exists:
python3 scripts/run_fgc_tdg8_rcv3_event.py \
  --authority-commit <full-auth1-commit> --status
python3 scripts/run_fgc_tdg8_rcv3_event.py \
  --authority-commit <full-auth1-commit> --run
```

The run command is explicit and is never part of ordinary verification.
