# FGC-1-RSP1-PREF12: checkpoint-recomputed RSP1 result

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-run numerical result; the targeted
amplitude-three veto is cleared for successor design, with no calibration,
candidate, or physical-mechanism promotion

```text
RSP1_study_terminated_normally = true
RSP1_all_six_members_reached_endpoint = true
RSP1_terminal_checkpoint_recomputed = true
RSP1_source_retry_total_zero = true
amplitude_three_spectral_veto_cleared_for_successor_design = true
RSP1_generic_all_field_raw_spectral_admission_passed = false
PROTO12_design_may_begin = true
PROTO12_frozen = false
fresh_GR0_dynamic_calibration_completed = false
FGCQR_holdout_execution_authorized = false
```

## Result

The prospectively frozen RSP1 study completed from immutable commit
`7f448bf6e5f0736e025c484deb664ba7887269da` under study identity
`0de794592ecafb7b93aff5b316d5b393a4007008c877535baf826a3c8d45f53f`.
It advanced amplitude `3` on exactly the declared
`4097 -> 8193 -> 16385` grids with the independent RK4 and SSPRK3 methods to
the only readable endpoint, `t=1/16`. The run took about 7,572 wall seconds.

All six members reached that endpoint. Their accepted steps were
`38, 73, 134` for RK4 and `20, 39, 77` for SSPRK3. Every member used the
prospectively frozen SRC3 source evaluator and recorded zero source retries.
The nonzero RK4 CFL trial counts remain separately serialized; they are
adaptive step-selection events, not source failures or scientific stops.

The unchanged target is the direct adjacent-grid ratio of the
`phi/Lambda` top-band derivative-weighted power fraction. The ceiling is
strictly `<1/4`; no saturation, epsilon floor, fitted asymptote, or threshold
change is permitted. The evolved values are:

```text
method    4097 -> 8193        8193 -> 16385
RK4       0.22812049822366257 0.18656614804836605
SSPRK3    0.2142551901915683  0.18877527466403243
```

All four ratios pass. By contrast, the corresponding public time-zero values
were `(0.8797163045076191, 0.3459398350356097)` for RK4 and
`(0.8846382581520159, 0.3459231990125775)` for SSPRK3, so both pairs failed
before evolution. The endpoint result is therefore neither a restatement of
the initial condition nor a threshold convenience.

Both methods also pass the complete target premises:

- the common-event constraint admissions pass, with minimum observed finite
  finest-pair orders about `2.3063` for RK4 and `1.9440` for SSPRK3;
- every complete field profile contracts in both infinity and RMS norms and
  through the declared spatial round trip;
- the medium/fine individual absolute spectral budgets pass; and
- the causal boundary remains excluded by a large positive margin for every
  member.

This clears CAL8's amplitude-three direct `phi/Lambda` spectral veto for the
design of a successor protocol. It also supplies evolved evidence that SRC3
removes the previously diagnosed source-evaluation wall on this six-member,
one-endpoint envelope without mutating the continuum equations or raw source
threshold.

## Independent checkpoint reduction

The tracked certificate does not trust only the runner's terminal JSON. Its
reproducer:

1. binds the raw manifest, two-record canonical event log, terminal checkpoint,
   and result by SHA-256;
2. verifies the immutable Git commit and every implementation blob named by
   the launch manifest;
3. reproduces the prospective RSP1 authorization and its six initial-state
   identities;
4. restores all terminal `u`, `p`, `q`, transaction, causal, and tracer state
   from the checkpoint;
5. recomputes the complete common endpoint from those arrays; and
6. requires exact serialized equality with both the terminal event and the
   study result.

The ignored raw checkpoint remains outside Git. The compact tracked result
contains its hash, the six terminal state hashes, the bounded diagnosis, and
the exact nonclaim boundary.

## Important all-field caveat

RSP1 was frozen to answer the one unresolved CAL8 target. It does not claim
that every raw nested-tail ratio for every field passes on the new ladder. In
fact, the generic all-field raw spectral aggregate is false for both methods at
the endpoint. That fact is serialized and remains public.

This does not contradict the RSP1 pass: its all-of decision contains the two
target ratios plus constraints, target absolute budgets, complete-profile
contraction, source health, and causal/boundary premises. It does mean that a
future PROTO12 cannot advertise blanket spectral convergence or silently
replace the repository's existing ownership and guarded-resolution rules.
Those rules must be inherited or justified prospectively.

## Successor boundary

The result authorizes design—not execution—of **FGC-1-PRO12-FRZ1**. Before a
fresh GR-0 trajectory is read, PROTO12 must prospectively freeze:

- the amplitude set and order;
- one common resolution ladder for every compared member, or an explicit
  predeclared reason for any branch-specific ladder;
- SRC3 as the sole source-evaluation graph;
- event times and the calibration decision rule;
- the existing constraint, health, boundary, spectrum, retry, rollback, and
  observable contracts; and
- a fresh output namespace.

RSP1 itself is not a completed GR-0 calibration and does not make amplitude
`3` an eligible case. It does not inspect SGB-L or FGC-QR, classify collapse or
trapping, test regulator activation, calculate an affine-null defocusing
margin, validate the retained EFT, resolve a singularity, derive a child
domain or dark sector, or vary a locally measured speed of light.

## Reproduction

With the immutable ignored raw namespace present:

```bash
python3 scripts/reproduce_fgc_rsp1_pref12.py
python3 scripts/reproduce_fgc_rsp1_pref12.py --check
```

The canonical configuration is
[`configs/fgc/fgc-1-rsp1-pref12.toml`](../configs/fgc/fgc-1-rsp1-pref12.toml),
and the compact result is
[`results/fgc-1-rsp1-pref12.json`](../results/fgc-1-rsp1-pref12.json).
