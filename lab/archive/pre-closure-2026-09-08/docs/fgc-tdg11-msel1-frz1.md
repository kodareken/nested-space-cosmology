# FGC-1-TDG11-MSEL1-FRZ1 — prospective temporal-method selection

## Status and scope

This owner freezes a prospective, one-shot numerical selection diagnostic.
It is not an executed result or a production admission authority. The run
requires the complete config, source closure, controls, and compact freeze to
be committed as the exact direct successor of the implementation base below.
The independent result binder must precede any `TDG11-IMP1` permission.

The implementation base is
`dbfe92df40bf18fd287663e0b9997a812fe4ccb2`. The scientific predecessor is
`FGC-1-TDG10-QA2-PREF1`, compact SHA-256
`08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6`.
QA2 rejected only the tested SSPRK3-on-inherited-SBP4 exact-complete-C remedy.
It selected no production successor.

Two origins must not be conflated:

- The measurement points are the authenticated retry-3/4/5 predecessors of
  the affected `RK4-2049` member, all at accepted time
  `0x1.78554de5a30e0p+0`, physical-state SHA-256
  `3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a`.
- The synchronized HLT15 `t=23/16` six-member origin is a separately pinned
  future production restart. It is not a replacement measurement state and
  no diagnostic endpoint may be transplanted into it.

| Retry | Restored generation | Checkpoint journal tip | Historical rejection | Exact attempted binary64 width |
|---|---:|---:|---:|---|
| 3 | 9 | 10 | 11 | `0x1.aaa9612df8000p-11` |
| 4 | 10 | 12 | 13 | `0x1.aaa9612df8000p-12` |
| 5 | 11 | 14 | 15 | `0x1.aaa9612df0000p-13` |

Later rejection records are evidence, never executable predecessors. The
retry-5 width is not replaced by an informal floating-point halving of retry
4. The original stage-safe lattice owns the exact boundaries.

The diagnostic measures the affected RK4/SBP4 construction. A separately
implemented compensated SSPRK3 tableau is covered by outcome-blind controls,
but no SSPRK3-on-SBP4 hybrid is relabeled as a production comparator. Full
two-method/three-resolution qualification remains HLT17 and the fresh-event
owner in the forward roadmap.

## Numerical objects and fixed selection order

Every candidate retains the eighteen ordered `u,p,q` channels, the owned-row
slice excluding the centre and four projector rows, and the exact
`8 * fine_upper**2 <= outer_lower**2` form of `p >= 3/2`, including equality.
No absolute state tolerance, signal normalization, threshold fitting, new
physics, or post-result candidate modification is permitted.

The candidates are evaluated in this fixed precedence:

1. exact-accumulation reconstruction with an explicit roundoff debit;
2. compensated coherent RK4/SSPRK3 accumulation;
3. embedded RK quadrature-defect estimation.

The first candidate passing its controls and every frozen width/channel is
selected for a separately owned `TDG11-IMP1`. A name in this list does not
establish a theorem or a measurement result.

### Candidate 1 — theorem-backed accumulation saturation

This route changes the *declared numerical reconstruction used for admission*,
not the accepted state, the integrator, or the equations. The raw binary64
complete-C assessment remains visible alongside it.

For one recorded RK step, let `y_j` be the actual binary64 endpoint, `h_j`
the exact rational value of the executed binary64 width, `k_ji` the exact
rational values of the recorded RHS arrays, and `b_i` the exact rational
Butcher weights. Define

```text
Delta_j = h_j * sum_i b_i k_ji
delta_j = y_(j+1) - y_j - Delta_j
z_0 = y_0
z_(j+1) = z_j + Delta_j.
```

Before `delta_j` may be called an accumulation defect, the implementation must
replay the original binary64 operation schedule from the recorded stage RHS
and initial state and match every owned internal-stage and endpoint byte.
An arbitrary endpoint perturbation, changed coefficient, projector change,
or unrecognized arithmetic schedule is invalid; it is not subtracted away.

The corrected piecewise cubic Hermite reconstruction uses `z_j,z_(j+1)` and
the same fresh endpoint RHS slopes as the recorded path. All values and
coefficients are exact rationals. The following finite-algebra facts are
proved and mutation-tested before measurement:

1. `y_j-z_j = sum_(i<j) delta_i` exactly.
2. The difference between the recorded and corrected Hermite segment has
   zero endpoint slopes. Its absolute maximum is bounded by the larger
   absolute endpoint defect, because the cubic smooth-step weights are
   nonnegative and sum to one on `[0,1]`.
   The public path bound is the row maximum of `sum_j abs(delta_j)`, not a
   signed sum: it bounds every cumulative endpoint defect without allowing
   roundoff cancellation to reduce the debit.
3. Corrected outer/medium/fine complete-C differences are enclosed by both
   independent rational localizers, and the unchanged TDG6 classifier is
   applied to their intersection.
4. The public fine-path debit is the corrected finest-pair upper bound plus
   the complete fine-path accumulation-defect bound. No correction is hidden
   or cancelled from the debit.

An outcome-blind large-base constant-RHS RK4 control makes this distinction
explicit. At `y0=2^52` and macro width `1/4`, the recorded binary64 paths have
complete-C maxima `D01=1/36`, `D12=1/72`, which fail the order gate. All three
corrected paths are the same exact affine polynomial, but the fine-path
roundoff debit is still `1/4`, not zero. Agreement of a reconstruction is not
permission to discard the actual path's error.

This is a theorem about a finite recorded RK data object. It is not a bound
on RHS evaluation error, a nonlinear ODE stability theorem, or a rigorous
global PDE error bound. Source, constraint, health, method agreement, spatial
convergence, and the later DEF1 stability map remain independent burdens.
The route must not claim that the original binary64 path acquired a resolved
order when only the corrected reconstruction did.

### Candidate 2 — compensated coherent accumulation

Create a new proposer rather than editing `numerical_engine.py` or any sealed
runner. It retains the exact RK4 and SSPRK3 mathematical tableaux, stage
abscissae, RHS, projector, and guard ordering, with a separate arithmetic
identity in every receipt.

RK4 stages remain the usual three base-plus-increment stages. SSPRK3 is
written in coherent Butcher form:

```text
y1 = y0 + h*k0
y2 = y0 + h*(k0+k1)/4
y3 = y0 + h*(k0+k1+4*k2)/6.
```

These are exact real-arithmetic identities of the existing Shu–Osher form;
they avoid repeatedly weighting the large base state. Every weighted update
uses a fixed compensated sum/product schedule, including a prospectively
fixed twofold representation of non-dyadic rational coefficients. Overflow,
unsupported underflow, nonfinite intermediates, or an arithmetic-contract
mismatch is a typed premise/resource stop, not an alternative precision run.

The implementation uses binary64 arrays only and records its exact operation
schedule. Independent rational controls verify the error-free primitives and
the coefficient residual enclosure; exact order-condition controls verify
both tableaux. No claim of universally improved accuracy is inferred from a
few controls.

The fixed arithmetic ID is
`tdg11_binary64_twofold_dot2_coherent_rk_v1`. Use Knuth's six-operation
`TwoSum` and Dekker's `TwoProduct` with splitter `2^27+1`. Every elementary
operation must remain finite and normal-or-zero; a nonzero product flushed to
zero is unsupported underflow. For exact coefficient `c`, set
`hi=RN(c)`, `lo=RN(c-hi)`, and retain the exact rational remainder
`c-hi-lo` in the coefficient controls. A twofold does not exactly represent
`1/3`. For each update start `acc=base`, `correction=0`, then, in the frozen
RHS order and `hi`-then-`lo` order, compute

```text
(product, product_error) = TwoProduct(coefficient_part, RHS)
(acc, sum_error) = TwoSum(acc, product)
correction = RN(correction + RN(sum_error + product_error))
answer = RN(acc + correction).
```

The last line runs once after all terms. This is a compensated fixed schedule,
not a promise of correctly rounded dot products. Its supported normal-domain
arithmetic and failures are separate from its measured temporal contraction.

The measurement uses the resulting *actual* binary64 endpoints and fresh RHS
records in the exact complete-C classifier. There is no admission override
for this candidate. A failure in one frozen channel rejects this candidate on
the measured neighborhood.

### Candidate 3 — embedded quadrature defect

This is a method-owned finite-data estimator, not a claim that finitely many
stages identify an unrestricted continuous PDE history.

For RK4, use the existing fresh endpoint RHS as a fifth stage and the
predeclared third-order embedded weights obtained by replacing the `k4`
weight with the endpoint weight. Its exact local quadrature defect is
`h*(k4-k_endpoint)/6`. For SSPRK3 controls, use the predeclared second-order
Heun weights on its first two stages; the defect is
`h*(-k0-k1+2*k2)/3`. The main/embedded order conditions and the first failed
higher-order condition are checked exactly, not inferred from measurements.

At each refinement level, accumulate absolute local defects without
cancellation, componentwise on every owned row, then take the row maximum for
each of the eighteen channels. Require the unchanged `p >= 3/2` test for
both outer-to-medium and medium-to-fine estimator contractions. Exact-zero
semantics remain explicit. The fine estimator and full fine accumulation
defect are retained as the public debit.

The embedded estimator is conditional numerical evidence, not a rigorous
global state-error enclosure. A true residual-to-solution bound requires an
applicable stability estimate; no such estimate is borrowed from a
principal-symbol or source-residual certificate. The later DEF1 error owner
must still close its independent stability burden.

## File and code solution map

The implementation is split by responsibility:

| New owner | Required contents |
|---|---|
| `tdg11_msel1_contract.py` | Exact candidate IDs/precedence, width and predecessor records, eighteen-channel order, resource ceilings, terminal enum, nonclaims, strict config parser. |
| `tdg11_msel1_reconstruction.py` | Exact rational RK accumulation, schedule-coherence checks, corrected Hermite rows, endpoint-defect bound, embedded quadrature defects, finite-algebra controls. No I/O or selection reduction. |
| `tdg11_compensated_rk.py` | Fixed compensated primitives and coherent RK4/SSPRK3 proposer, separately identified arithmetic, finite/underflow guards, no state commit. |
| `tdg11_rational_complete_c.py` | Rational cubic-family adapter to the two existing independent localizers; new hash domains, exact bounds/counts/digests, unchanged contraction predicate. No campaign access. |
| `tdg11_msel1_runtime.py` | Shadow-only guarded path construction, complete record collection, all-of per-width/channel assessment, explicit candidate status. No historical store writes. |
| `tdg11_msel1_authority.py` | Compact predecessor/origin pins, clean direct-successor image and exact delta, no-follow source/store snapshot, environment and namespace/process preflight. |
| `scripts/reproduce_fgc_tdg11_msel1_frz1.py` | Outcome-blind compact freeze and control reproduction; no proposal or raw terminal creation. |
| `scripts/run_fgc_tdg11_msel1.py` | Isolated one-shot status/run entry, authority recheck, bounded shadow execution, atomic diagnostic terminal only. |

All modules live under `src/recursive_horizons/fgc/evolution/` except the
named scripts. The exact config is
`configs/fgc/fgc-1-tdg11-msel1-frz1.toml`; the compact freeze result is
`results/fgc-1-tdg11-msel1-frz1.json`. Focused tests mirror each owner and
attack cross-owner substitutions.

Reuse is limited to the unchanged state/RHS types, spatial operator,
projector, stage health/causal/constraint machinery, authenticated persisted
state decoder, stage-safe lattice, exact rational localizers, and neutral
future evidence I/O. Historical private transaction helpers may be isolated
behind one explicitly pinned adapter only where no public equivalent exists.
No `run_fgc_gr0_calibration_v*` import is introduced. No historical source or
compact result is rewritten.

The current Make fragment receives the new freeze, focused verifier, status,
and one-shot run targets. The artifact catalog receives exactly one active
prospective authority for the new run target; all 21 historical mutation
tombstones remain closed. The full repository checker gains a focused TDG11
compact route without making compact verification raw/store/Git dependent.

The pure reconstruction boundary accepts three tuples of recorded
`StepProposal` objects, of lengths `(1,2,4)`, plus the method and arithmetic
IDs. It verifies every stage name/time, byte-coherent owned state, shared
endpoint/RHS, exact uniform subdivision, and common initial state before
forming any cubic. It exposes one channel at a time as a stream of rational
Hermite segments `(left_value,left_rhs,right_value,right_rhs,width)`; corrected
segments additionally accumulate the endpoint defect without cancellation.
The localizer adapter accepts only these exact-rational row streams and
returns exact bounds, count/digest receipts, and the inherited TDG6 channel
decision. It has no runner, raw-path, restoration, or selection dependency.

The runtime owns the only proposer-substitution adapter. It composes the
pinned PROTO7 accepted-source precheck, complete-proposal failure preview,
tracer preview, and acceptance guard in their existing order, with fresh
cloned transaction and tracer objects for each path. It never assigns a
replacement to `proto7_runtime.propose_step` or any historical module global.
The existing `TDG6ShadowPath` is a compatible in-memory record; its old
floating envelope assessment and production compositor are not invoked.

Restoration uses `HLT16CampaignStore.authenticated_checkpoint_at_generation`
and `restore_member_with_overlay`, with repairs disabled and no writer
capability. `build_static_gr0_shells` is called once to create its six pinned
time-zero templates; only the affected member is restored for measurement.
The six template source prechecks are separate from the 42 diagnostic
accepted-source prechecks and 210 retained stage/endpoint records. No new
module imports an executable historical runner merely to restore a member.

## Execution and independent binding

One original and one compensated RK4/SBP4 path family are constructed for
each of the three widths. Each family has one outer, two medium, and four
fine proposals. The measurement budget is therefore 42 proposals and 210
stage-plus-endpoint records, plus the separately counted accepted-state
source prechecks. Candidates 1 and 3 consume the original family; candidate
2 consumes the compensated family. Reusing recorded inputs does not merge
their classification decisions.

The order is fixed: construct the three original-RK4 families at retries
3, 4, and 5 first; then the three compensated-RK4 families in that same retry
order. On each original family assess the raw baseline, corrected complete-C,
and embedded defect, in that order, each over the eighteen ordered channels.
On a compensated family assess its actual complete-C. A stopped group is
never rerun or retuned. A localizer/premise stop stays visible while other
predeclared groups may still be evaluated; a global resource stop ends the
remaining schedule. Provenance, implementation-integrity, or store drift
prevents publication of a scientific terminal.

Each candidate passes only with all 54 width/channel gates complete and
passing. One valid sufficient nonpass rejects that candidate's all-of claim
even if other channels are unresolved; an unresolved channel is never itself
called a sufficient failure. Without either a complete pass or a valid
counterexample the candidate is inconclusive. The first passing candidate is
recorded as an **unbound diagnostic selection**; the raw terminal explicitly
keeps `licenses_only_separate_TDG11_IMP1=false` and requires PREF1 first.

The wall budget is 14,400 seconds, the raw bundle is at most 16 MiB, and exact
rational serialization is capped at 32,768 bits per numerator/denominator.
These are resource ceilings, not scientific tolerances. Git inventory reads
are capped at 16,384 entries; directory publication retains the independent
1,024-entry infrastructure bound and this diagnostic permits only two leaves.

The new raw namespace is `runs/fgc-2-sf1/tdg11-msel1`. It will contain only a
canonical manifest and terminal, published as one exclusive directory under
the already existing `runs/fgc-2-sf1` parent. It
serializes no diagnostic endpoint, payload, cursor, journal, checkpoint, or
accepted state. The historical 115-leaf store must remain byte-identical.
The manifest binds the authority commit, exact environment, implementation
hashes, source-store manifest, predecessor fingerprints, counts, and terminal
hash. Publication is exclusive and atomic; ambiguous publication stops
without informal retry.

After the raw terminal, `FGC-1-TDG11-MSEL1-PREF1` is a separate owner. Its
binder authenticates the raw leaves and authority, reconstructs the path
families, independently recomputes rational reconstruction/defect reductions
and every candidate width/channel classification, and independently applies
the fixed precedence. It imports no runner selection or publication decision.
The numerical proposer may be replayed as the hash-bound implementation under
test; independent exact controls, not shared labels, establish its algebra.

Ordinary compact verification of PREF1 is raw-, store-, shadow-, and
Git-blind. Live construction is an explicit one-time operation. The raw
baseline exact-C labels remain visible and cannot be replaced by the selected
candidate's labels.

## Outcome branches

- A candidate passes only if its controls and all three widths/all eighteen
  channels pass. Select the first passing candidate in the frozen order and
  license only the separate `TDG11-IMP1` implementation/qualification.
- Valid nonpasses of all three candidates close this tested numerical branch
  and require a new solver/formulation plan. They do not reject FGC-QR physics.
- A premise/resource/provenance interruption is typed and independently
  bound. If no candidate has a valid pass and an unresolved candidate remains,
  the terminal is inconclusive rather than a scientific nonpass. Another
  measurement requires a new prospective freeze.

No branch opens SGB-L, FGC-QR, holdout, calibration, retained-EFT evolution,
mechanism, physical transition, publication, or push.

## Literature boundary

Error-free sum/product transforms provide the arithmetic building blocks,
not a PDE convergence theorem; see
[Ogita, Rump, and Oishi (2005)](https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf).
Coefficient representation and update accumulation can both affect RK
roundoff; the long-time Hamiltonian experiments of
[Hairer, McLachlan, and Razakarivony](https://www.unige.ch/~hairer/preprints/roundoff.pdf)
are context, not evidence for this collapse problem. A residual becomes a
solution-error bound only through an applicable stability estimate, as made
explicit by
[Ranocha and Giesselmann](https://arxiv.org/html/2307.12677v2). None of these
sources licenses a post-outcome tolerance or a probabilistic roundoff floor.

## Finish-line proof

The first claimable result is the independently bound one-shot selection
terminal, not production advancement. Its proof route is:

1. prospective compact freeze/controls and exact committed authority;
2. read-only status with `safe_to_run=true`;
3. one bounded run and unchanged source-store manifest;
4. independent PREF1 reconstruction/reduction and adversarial compact tests;
5. the current-development and canonical repository verification gates,
   followed by a clean task-owned Git boundary.

## Failure forecast

The most likely failures are that a corrected reconstruction still fails a
width/channel (distinguished by preserved raw and corrected exact bounds),
that compensation encounters an unsupported arithmetic case or fails to
remove the observed obstruction (distinguished by its typed arithmetic stop
versus exact complete-C nonpass), or that an estimator/localizer resource or
premise is insufficient (distinguished by a typed inconclusive terminal and
unchanged state/store fingerprints). No branch permits threshold fitting,
silent method relabeling, or converting a numerical failure into physics.
