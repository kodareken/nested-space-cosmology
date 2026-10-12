# CF4 temporal accuracy on the actual local prepared-state control

The midpoint control has a measured local60/120-step matter difference
(1.66130e-7,1.40888e-7). Reuse its exact history alpha=.001,w=T1*plateau,U=0,
the resolved7601 initial columns, four signed energies, all three coherent
source fibers and the original channel/multiplicity. Only the time method
changes to the cited CF4 owner; no initial/spatial/source resampling occurs.

The bounded comparison runs60 and120 steps on[0,.18]. A240-step case is
allowed only if their maximum raw N,beta matter difference on I exceeds
1e-11. Every case carries the actual retarded tangent and PDE F_z,dF_z.
Maximum three solves and360 CPU seconds, with checkpoints after each.
No automatic fourth case follows a failed indicator or CPU cap.

Compare nested time nodes directly. Record the raw matter and tangent
changes, fixed-preparation equality, source-phase and flux algebra. The
partial constraint diagnostic still includes the unmodified approximate
baseline and the actual geometric/state changes; no numerical drift is
removed. A time-comparison PASS is not a continuum error bound, a complete
source, or a local constraint root. Those physical gates stay OPEN.

```sh
python3 scripts/derive_nsc_local_cf4_response.py --run
python3 scripts/derive_nsc_local_cf4_response.py --check
```

The check path replays saved restrictions and provenance without propagation.
