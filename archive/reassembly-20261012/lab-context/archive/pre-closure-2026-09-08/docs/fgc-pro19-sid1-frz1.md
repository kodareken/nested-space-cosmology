# FGC-1-PRO19-SID1-FRZ1 — generation-seven suffix identity freeze

SID1 is a read-only, premise-only authority for one already-written GR-0
campaign suffix. It binds generation-seven checkpoint
`c0ebfe09ef4058e117e46d34424edc8139419bfb4680cec979a1f5ab64d37be9`,
the checkpoint-owned CFL journal record at sequence six
`78e65c08380a4fc979269e21db314232fd3fdc0e428c7e259f2031e66bd30bf7`,
and the later uncheckpointed TDG6 rejection at sequence seven
`e23d81705f1038f8240249658d0160cd5ca8ea8e33969cb4b32de3f03d426473`.

Those are three different hashes with three different meanings. The
generation-seven checkpoint binds the sequence-six CFL record as its journal
tip. The sequence-seven TDG6 record has that CFL record as its parent and is
only a crash suffix; it is not a checkpoint and cannot be substituted for one.

The authority also binds the preserved `RK4-2049` descriptor, raw payload,
semantic payload, and independently hashed evolution arrays `u,p,q`. It
recursively reconstructs `Binary64CubicEnvelope`, then the TDG6 magnitude
intervals, channel admissions, and retry evidence in memory. It does not run
a PDE proposal, take a writer lease, reconcile the suffix, mutate the store,
or open SGB-L or FGC-QR.

The exact successor was subsequently published after immutable SID1 commit
`c83ef958...` and is independently bound by
[FGC-1-PRO19-SID2-PREF1](fgc-pro19-sid2-pref1.md). That later fact does not
rewrite SID1: this artifact itself remains a pre-recovery freeze and earns no
calibration, trapping, activation, defocusing, EFT, or physical claim. SID2
also does not unlock numerical resume; a separate resume authority and
read-only preflight remain mandatory.
