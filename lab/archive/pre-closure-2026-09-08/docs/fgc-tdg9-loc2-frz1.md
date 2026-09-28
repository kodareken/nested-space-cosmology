# FGC-1-TDG9-LOC2-FRZ1: scale-invariant independent localization freeze

`FGC-1-TDG9-LOC2-FRZ1` is the direct, still-unexecuted successor to LOC1
authority commit `14e26c43edbc88003d68cbcf24681bb090d98534`. It consumes the
premise-only [`FGC-1-TDG9-LOC1-ERR1`](fgc-tdg9-loc1-err1.md) diagnosis and
changes only the independent stationary-root instrument.

The v2 route normalizes every derivative to primitive integer
`A t^2 + B t + C`, with deterministic positive leading sign. Degenerate and
linear cases are exact. Roots at `t=0` or `t=1` are exactly deflated before the
remaining rational factor is solved. A perfect-square integer discriminant is
solved exactly. Only a genuinely nonsquare discriminant enters the un-clipped
dyadic enclosure route.

For a primitive height bit-length `h`, the per-polynomial proof cap is
`B_proof=2h+8`. Refinement begins at 192 bits, doubles deterministically, and
ends with the exact `B_proof` value when it is not already in the schedule.
The global authorization ceiling is 9,216 bits. Exact Sturm variation counts
over `(0,1)` must match the strictly interior, ordered, disjoint root
enclosures. Proof-cap exhaustion, location, separation, or root-count
uncertainty is a distinct typed invalid diagnostic; no boundary-straddling
interval is clipped into an interior root.

The prospective audit covered the unchanged `122,640` base and `367,920`
`V/S/C` component cubics. It found maximum primitive height `h=96` and maximum
proof cap `200`, far below 9,216. The raw coefficient, derivative, and
discriminant numerator/denominator maxima are recorded in the compact
authority. This audit did not localize extrema or read a LOC2 outcome.

LOC2 retains the ten PREF1 failures, retries `3..5`, exact TDG6 pass `<=` and
strict-failure `>` predicates, `490,560` per-family and `1,471,680` aggregate
candidate ceilings, explanatory-only ULPs, two exact routes, full survivor
interval comparison, and an independently generated per-polynomial stationary
count digest. Before localization, the runner must match ERR1's authenticated
115-leaf construction snapshot; after localization it must match the same
snapshot again. The campaign store is read-only. A successful diagnostic may
publish only canonical `manifest.json` and `terminal.json` under
`runs/fgc-2-sf1/tdg9-loc2/rcv3-retries-3-5`.

No fourth width, Stage-2 source decomposition, SSPRK3 comparison, PDE state
commit, continuation, calibration, candidate branch, mechanism result, or
physical result is authorized.

```bash
make fgc-tdg9-loc2-frz1
make verify-fgc-tdg9-loc2-prelaunch
python3 scripts/run_fgc_tdg9_loc2.py --authority-commit <committed-authority>
```
