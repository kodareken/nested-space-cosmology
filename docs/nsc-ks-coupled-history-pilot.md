# Coupled-history source evaluator pilot

This bounded control exercises both local history functions with the
original positive and negative energy sources for both angular signs of
group14. Saved upstream columns and their preparation digests are restored
by `RetainedUpstreamArchive`; preparation is not repeated.

The initial geometric guess uses `w(center)=w''(center)=0`,
`w'=-baseline_beta/F(0)` and
`U=-(baseline_N+C*(w')²)/A(0)`. These are solver controls only. The state
is then newly evolved from the original upstream data through this history,
and the same-action source-cutoff edge and its retarded derivative enter
both constraints. No frozen-C0 root is relabelled as physical.

Both functions have eight Chebyshev coefficients. The64-point spatial cell
and degree32 energy interpolant on[0,160] are a throughput control, not a
numerical-accuracy certificate. Every operator carries all16 retarded
directions and uses the common union of eight solve and17 verification
nodes. Two original angular operators are allowed, within30 CPU seconds.

Only group14 matter changes are included. Other retained groups, finite
tails, source and field accuracy, and between-node errors remain OPEN.
The pilot does not run a physical optimizer or claim a local gate result.
Its purpose is to measure the new candidate's computation and retain its
actual coupled response before selecting a budget for a full source solve.

Saved operators replay the complete source reconstruction and constraint
assembly without repeating field evolution. Original source values and
historical records are unchanged.

```sh
python scripts/derive_nsc_ks_coupled_history_pilot.py --check
```
