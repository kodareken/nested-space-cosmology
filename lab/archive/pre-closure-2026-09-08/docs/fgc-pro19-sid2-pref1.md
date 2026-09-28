# FGC-1-PRO19-SID2-PREF1 — post-recovery generation-eight binder

SID2 is an independent, read-only binder for the one metadata-only recovery
authorized by SID1. It authenticates the generation-seven TDG6 rejection,
the sequence-eight cursor transition, and the complete generation-eight
checkpoint as one exact parented chain. It also proves that all six persisted
state descriptors are unchanged across that recovery and independently
reopens the affected `RK4-2049` payload to bind its semantic archive, raw
archive, and `u,p,q` identities.

The recovery changed no accepted finite field. It converted an authenticated
uncheckpointed TDG6 rejection into a durable `RETRY_PENDING` cursor owned by
the temporal controller. The stored cap is
`0x1.aaa9612df9000p-10`; the already-derived lattice-safe successor macro
width is `0x1.aaa9612df8000p-10`. The rejected predecessor width
`0x1.aaa9612df9000p-9` may not be replayed.

Recovery appended the already frozen sequence-eight cursor transition and
complete generation-eight checkpoint. It did not execute or replay a PDE
proposal, change any accepted descriptor or boundary time, complete event 23,
or authorize trajectory resume.

The live audit uses an independent non-following traversal and rejects
symlinks, special files, foreign leaves, staging residue, forks, terminal
state, or an active writer. Ordinary verification consumes only SID2's
canonical compact evidence after this temporal boundary has been sealed; it
must not later demand that a legitimately progressed store still equal the
generation-eight snapshot.

SID2 does **not** authorize numerical resume by itself. A later committed-image
authority must bind the execution-source closure, exclude untracked import
shadows, reconstruct the original campaign plan rather than a new plan, and
pass a read-only resume preflight against this exact ancestor. Until that
separate gate closes, `safe_to_restart=true` means only that the durable store
is recoverable; `safe_to_resume_trajectory` remains false.

No common event, GR-0 calibration, SGB-L or FGC-QR trajectory, activation,
trapped interval, defocusing result, retained-EFT conclusion, transition, or
physical claim is contained in SID2.
