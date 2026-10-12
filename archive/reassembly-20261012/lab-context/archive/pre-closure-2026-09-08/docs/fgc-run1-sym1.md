# FGC-1-RUN1-SYM1: scoped spherical run authorization

## Decision

`FGC-1-RUN1-SYM1` is the fail-closed authorization contract for one bounded
classical spherical diagnostic of the frozen FGC-QR equations. Its answer is
**no**. PROTO3 passes the frozen structural schema,
`FGC-1-SRC1-NL1` closes the nonlinear-source predicate, and
`FGC-1-CON4-PHY1` closes the conditional boundary-free physical/gauge/reduction
constraint predicate. `FGC-1-CTR1-REG1` and `FGC-1-ID1-FAM1` now jointly close
the regular-centre/data-family predicate. `FGC-1-DOM4-RUN1` and
`FGC-1-HYP2-MD1` now jointly close the classical run-domain and
multidirectional early-kill predicate. `FGC-1-BND2-CP1` closes the scoped
boundary row through strict causal isolation, and `FGC-1-HLT1-MON1` closes the
transactional classical-health-monitor row on direct and injected controls.
`FGC-1-NUM1-VAL1` closes the independent numerical and affine-measurement row
using two discretizations, unchanged-source cross-checks, centre, restart,
transaction, boundary, and physical-null controls without opening an FGC-QR
trajectory. In this preserved composition seven of eight structural predicates
pass and the resolved outcome-neutral holdout manifest is absent. The later
[`FGC-1-CAL0-PREF2`](fgc-cal0-pref2.md) preflight proves that this manifest
cannot honestly be resolved under PROTO3 because the declared compact pulse
fails its own literal spectral rule and calibration lacks a branch-specific
stop-applicability map. Consequently:

```text
classical_spherical_diagnostic_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
FGCQR_holdout_execution_authorized = false
```

This result remains the immutable PROTO3-era authorization record.
**FGC-2-SF1-PROTO4** and
[`FGC-1-PRO4-FRZ1`](fgc-pro4-frz1.md) now freeze the permitted premise-only
successor, but they do not retroactively change RUN1 or make its missing row
true. `FGC-1-HLT2-MON2` now implements the successor admission formulas, and
`FGC-1-ID2-ALL1` now closes the outcome-neutral all-amplitude/all-case static
ledger with `[5/2, 3]` as the only fresh-calibration candidates. Neither is a
trajectory or a successor authorization result. Fresh nested GR-0 calibration,
the SGB-L-owned health route, the DEF1 stability map, and PRO4-HLD1 must still
be independently constructed before a successor RUN1 composition can be
evaluated.

This is not a second version of `FGC-1-EFT1-OPEN1`. EFT1 asks whether the
candidate is authorized as a retained Wilsonian EFT. RUN1 asks the narrower
question of whether one symmetry-reduced calculation of the declared
classical equations has enough continuum, constraint, health, numerical, and
measurement control to test the proposed mechanism without confusing a code
artifact with physics.

The logical direction is one way:

```text
physical-transition claim
        => retained-EFT authorization
        => scoped classical diagnostic prerequisites
```

No converse is assumed. A future scoped calculation may become authorized
while retained-EFT and physical-transition authorization remain false.

## Why this gate exists

The post-COMP1 ladder established a real local/radial analytic foundation, but
it did not establish an initial-boundary-value problem or a collapse solution.
Starting a solver directly from that ladder would silently turn numerical
choices into missing continuum assumptions. Requiring all of EFT1's open
Wilsonian and nonspherical programme before asking whether the classical
mechanism occurs at all would instead couple a bounded falsification to a much
larger physical-promotion problem.

RUN1 separates those questions without weakening either one. It permits only
an outcome-neutral mechanism diagnostic after every scoped premise has its own
hash-bound result artifact. A configuration boolean is never evidence. Missing
evidence is false.

## Frozen study protocol

The active source of truth is
[`configs/fgc/fgc-2-sf1-protocol-v3.toml`](../configs/fgc/fgc-2-sf1-protocol-v3.toml),
artifact **`FGC-2-SF1-PROTO3`**. The original
[`FGC-2-SF1-PROTO1`](../configs/fgc/fgc-2-sf1-protocol.toml) remains unchanged
as the failed-premise record. ID0 proved before any candidate outcome that
PROTO1 omitted the regulator unit-normal momentum and that its amplitude box
could not reach its own compactness floor under the minimal completion.
**FGC-2-SF1-PROTO2** therefore makes only three authorized changes: it freezes
`Pi_phi=0` in the future-slice unit-normal frame, replaces the amplitude list
with the GR-0-only preflight scan `[2, 5/2, 3, 7/2, 4, 9/2, 5]`, and writes to
`runs/fgc-2-sf1/proto2/holdout`. Every other scientific and numerical rule is
validated by exact reduction to the sealed PROTO1 schema.

ID1 then exposed a second premise contradiction before evolution: the already
declared `9/8` width leaves the exact positive centre-vacuum buffer
`3-(1/2)(9/8)=39/16 L0`, while PROTO2 required the larger round floor `5/2 L0`.
PROTO3 changes only that floor to `39/16 L0`, versions the future manifest and
output namespace, and requires a static all-case initial-premise ledger before
GR-0 dynamic calibration may qualify. The immutable PROTO2 and ID1 diagnosis
blobs at commit `fae2c5b76ba53680f9e70188f79784b160826734` are checked directly;
the widest support retains 78 complete coarsest-grid intervals of exact vacuum.
No SGB-L or FGC-QR evolution outcome was inspected. The PROTO3 amendment and
its provenance enter the active semantic hash. The protocol freezes the scientific
choices before the FGC-QR holdout is opened:

- branch order `GR-0`, `SGB-L`, then `FGC-QR`;
- the unredefined ACT1/VAR1 equations as the physical equations;
- the modified-harmonic terms as a formulation that vanishes on the
  metric-defined gauge surface, not a new physical source;
- the exact compact-bump profiles for `chi` and `phi`, their centre and vacuum
  buffers, and the same physical `chi` profile and normal derivative across
  branches before each branch-specific constraint solve;
- one ordered GR-0-only amplitude calibration: the first candidate that forms
  a converged trapped sphere without initial trapping and passes the comparator
  method is selected; if no candidate qualifies, the calibration terminates as
  `calibration_failed_no_eligible_GR0_case` and FGC-QR execution stays closed;
- disjoint development, calibration, and held-out cases;
- an explicit nonzero `phi` seed of amplitude `1/131072` and explicit
  unit-normal momentum `Pi_phi=0`, with activation by numerical roundoff
  forbidden;
- three nested grids, a primary SBP4/RK4 method, and an independent
  SBP2/SSPRK3 comparator;
- typed branch, hyperbolicity, constraint, scale-control, boundary, and
  affine/measurement stops;
- one fixed future-null orientation and a physical-metric null generator
  launched from `t=0`, labelled by its initial areal radius, normalized by
  `-g_ab k_+^a n^b=1`, assigned affine origin `lambda_+=0`, and transported by
  `k_+^b nabla_b k_+^a=0`, with a frozen residual bound and cross-resolution
  trajectory alignment;
- the complete spherical Raychaudhuri expression, rather than the Ricci term
  alone;
- a positive margin greater than four times the combined Richardson,
  constraint, initial-data, nonlinear-source, gauge, affine, null,
  trajectory-alignment, interpolation, extraction, and boundary error
  estimate;
- a dynamically updated all-cone boundary domain-of-dependence check,
  including stencil reach and a required outer-boundary displacement test;
- a rational nonzero physical robustness box separated from numerical and
  gauge diagnostics, with an interval or continuity certificate required
  before finite samples may be called an open neighborhood; and
- a revision rule that reclassifies post-unblinding changes as exploratory.

The numerical values are study coordinates, not a fit to nature. In
particular, the cutoff is externally declared and is not inferred from the
same local fixture it is meant to control.

The protocol itself does not resolve the GR-0 calibration outcome. A separate
`FGC-1-PRO3-HLD1` manifest must convert the deterministic selection rule into
fully specified, hashed holdout inputs before any FGC-QR outcome is inspected.
That manifest must serialize every ordered calibration candidate, its config
and result hashes, eligibility decision, the selected amplitude, every expanded
holdout case hash, an initially empty output inventory, and an append-only
execution-event-log hash. Before GR-0 dynamic eligibility, it must additionally
serialize an ordered static ledger proving constraint compatibility, regular
centre, finite mass, compactness-window membership, no initial trapping, and
the exact centre buffer for all five expanded FGC-QR inputs without inspecting
an FGC-QR evolution outcome. The runner must refuse overwrite or pre-existing
holdout output.

This is a machine-verifiable **repository sequencing** claim, not a claim that
private human observation can be proved absent. The protocol says so
explicitly. Each future certificate must bind the exact RUN1 config, semantic
protocol hash, ACT1 and VAR1 identities, action/variation source hashes, branch,
and common run-envelope hash. All simultaneously present future certificates
must agree on that envelope. A result for a different action, protocol, or
run box therefore cannot accidentally unlock this study.

## Inherited evidence and its boundary

RUN1 consumes the following records without restating their proofs:

| Record | Positive input | Boundary retained by RUN1 |
|---|---|---|
| `FGC-1-ACT1` and `FGC-1-VAR1` | Frozen action and unredefined covariant bulk equations | No nonlinear collapse solution |
| `FGC-1-HYP1-DOM3-UHYP1` | Compact nonflat radial branch graph, real eigenframe, and symmetrizer | No all-covector or nonspherical theorem |
| `FGC-1-HYP1-CON2-MPROP1` | Metric-defined gauge subsidiary and kinematic reduction identities | No complete physical constraint system |
| `FGC-1-HYP1-CON3-CAU1` | Conditional boundary-free gauge preservation | Assumes a sufficiently smooth compatible full solution |
| `FGC-1-HYP1-BND1-MD1` | Frozen annular radial main-system dissipation | No regular centre or nonlinear constraint-preserving IBVP |
| `FGC-1-DEF0-OBS1` | Exact affine physical-null observable and sign controls | No trajectory or defocusing result |
| `FGC-1-EFT1-OPEN1` | Explicit all-of retained-EFT audit | Remains negative at `3/12`; evolution is not authorized as retained EFT |
| `FGC-1-SRC1-NL1` | Same-evaluator floating REF1 residual/Jacobian and safeguarded QIFT1-local acceleration solve | No run-domain existence theorem, constraints, centre, evolution, or holdout outcome |
| `FGC-1-CON4-PHY1` | Unredefined Hamiltonian/momentum projections, exact physical-to-normal-gauge map, composed boundary-free subsidiary theorem, reduction closure, and normalized monitors | Conditional on a smooth solution before boundary arrival; no initial-data family, regular centre, constraint-preserving boundary map, or evolution |
| `FGC-1-CTR1-REG1` | Exact `v=rV`, `R=rA` centre variables, parity/elementary-flatness enforcement, rederived reference limits, finite REF1 Taylor system, and exact curvature/defect/first-grid controls | No compatible finite-mass data family, centre existence theorem, IBVP, evolution, or holdout outcome |
| `FGC-1-ID0-PREF1` | Exact pre-holdout PROTO1 obstruction plus verification of its canonical record and the prior RUN1 `2/8` stop at immutable Git checkpoint `d4f0cc8f58408ef4ee12fb32fe3619916e231795` | Rejects only PROTO1 as written; it contains no SGB-L or FGC-QR outcome and cannot authorize ID1 or evolution |
| `FGC-1-ID1-FAM1` | Complete specialized FGC-QR Hamiltonian/momentum solve, exact full-evaluator and REF1-gauge controls, regular centre, finite-mass vacuum exterior, two-method convergence, and a nonzero open local data family | Initial hypersurfaces only; no evolution, holdout outcome, box-wide interval enclosure, or qualification of PROTO3's nine-eighths-width case |
| `FGC-1-DOM4-RUN1` | Nonzero parameter/state container, complete ID1-to-REF1 initial-slice bridge, strict branch-continuous initial acceleration witnesses, and a separately hashed quantitative domain definition | No global trajectory tube, all-case initial-data admission, nonlinear existence time, evolution, or holdout outcome |
| `FGC-1-HYP2-MD1` | Full 12-field covariant principal symbol and coefficient-complete all-spatial-covector sufficient health ball, with exact spherical restriction and strict initial witnesses | Classical principal health only; no retained-EFT authorization, trajectory containment, nonlinear existence theorem, boundary result, evolution, or nonspherical robustness |
| `FGC-1-BND2-CP1` | Exact local-to-coordinate all-cone speed conversion, strict immutable causal ledger, finite-difference extraction pad, and grid-aligned outer-boundary variants | Causal isolation of the retained measurement interval only; no nonlinear constraint-preserving boundary theorem, evolved boundary-invariance result, or trajectory |
| `FGC-1-HLT1-MON1` | Canonical dimensionless monitors and an immutable pre-acceptance transaction with all 26 typed failures independently injected | Synthetic monitor correctness only; no accepted trajectory, numerical validation, retained-EFT control, collapse, or defocusing |
| `FGC-1-NUM1-VAL1` | Independent SBP4/RK4 and SBP2/SSPRK3 convergence; exact/scalar/batched REF1 source agreement; centre-limit, atomic-stop, restart, boundary, null, affine, and complete Raychaudhuri controls | Outcome-blind implementation validation only; no FGC-QR trajectory, collapse, activation, trapping, defocusing, retained-EFT control, or transition |

The inherited results are necessary inputs. None substitutes for the future
evidence below.

## The eight all-of predicates

The scoped diagnostic becomes authorized only if every row passes through its
declared result gate.

| Predicate | Required future artifacts | What failure means |
|---|---|---|
| Frozen outcome-neutral protocol and resolved holdout hash | `FGC-1-PRO3-HLD1` | FGC-QR inputs and their static initial premises were not fully sealed before outcome access |
| Classical run domain, branch, and multidirectional early-kill envelope | `FGC-1-DOM4-RUN1`, `FGC-1-HYP2-MD1` | The intended run may leave its healthy branch or fail in an untested covector direction |
| Nonlinear unredefined REF1 source and acceleration solver | `FGC-1-SRC1-NL1` | A computed trajectory is not verified as a solution of the declared equations |
| Physical, gauge, and reduction constraint system | `FGC-1-CON4-PHY1` | Constraint error could imitate the proposed mechanism |
| Regular centre and finite-mass compatible data family | `FGC-1-CTR1-REG1`, `FGC-1-ID1-FAM1` | There is no regular nonzero-width collapse problem |
| Boundary or domain-of-dependence control | `FGC-1-BND2-CP1` | Incoming or reflected error could contaminate the measured interval |
| Classical health and typed stop monitoring | `FGC-1-HLT1-MON1` | A failed premise could be crossed without terminating interpretation |
| Independent solver validation and affine measurement contract | `FGC-1-NUM1-VAL1` | Implementation or measurement error could be mistaken for a scientific outcome |

One false predicate stops the FGC-QR holdout. A future result file without its
declared source configuration, or either member of a partial config/result
pair, is rejected rather than interpreted as absent evidence.

The preserved PROTO3 composition is therefore `7/8`: the nonlinear-source row, the
physical/gauge/reduction constraint-system row, the combined
regular-centre/finite-mass-data row, and the combined classical run-domain/
multidirectional early-kill row pass, together with the BND2 causal-isolation
row, HLT1 transactional-health row, and NUM1 numerical/measurement row. DOM4 supplies a nonzero declared
container and strict initial-slice source/branch witnesses; HYP2 supplies a
coefficient-complete sufficient principal-health bound for every spatial
covector. Neither result proves that a future trajectory remains in that
container; HLT1 can only stop a future stage that reports the required
diagnostics and does not prove such a stage exists. BND2 similarly proves the
declared measurement region can be isolated from the outer boundary under its
strict debit rule, not that a future solver will satisfy the budget. ID1 obtains the centre/data
premise from a smooth radial ODE with strict central margins plus nearby
two-method stress points; it does not claim an interval enclosure of its entire
sampled box. PROTO3 resolves the contradictory round-number buffer requirement,
but ID1 still does not qualify the expanded nine-eighths-width slice. ID2 now
supplies the required outcome-blind ledger through independent constraint
solves for all seven amplitudes and all five future FGC-QR modifiers. It finds
31/35 FGC-QR slices statically admissible and freezes `[5/2, 3]` as the only
fresh-calibration candidates; that later evidence does not rewrite this
PROTO3-era RUN1 decision. NUM1 establishes that the declared machinery
passes its independent controls; it does not establish that PROTO3's declared
pulse passes HLT1, or that an FGC-QR stage exists or succeeds. CAL0 closes that
missing composition and forbids resolving PRO3-HLD1 by reinterpretation. A
PROTO4 freezes branch applicability and convergent spectral/constraint
admission, HLT2 implements those exact formulas, and ID2 closes their static
pre-calibration application. Fresh calibration and a new manifest remain
required.

## Outcome semantics

RUN1 authorization is permission to perform the held-out calculation, not its
scientific result. Subsequent run records must distinguish:

```text
invalid_implementation_or_nonconverged_run
stopped_branch_loss
stopped_hyperbolicity_loss
stopped_constraint_loss
stopped_scale_control
stopped_boundary_contamination
stopped_affine_or_measurement_contract_loss
completed_no_defocusing
completed_defocusing_candidate
```

A crash, failed Newton iteration, or unconverged trajectory is an invalid run.
A single typed stop is an obstruction candidate, not a publication-grade
negative result. A scientific negative requires either completion of the
predeclared admissible family with no defocusing or a converged obstruction
whose location and type persist under the relevant neighborhood, resolution,
and method checks. A branch-rejection certificate must name the exact
parameter-box x data-family x formulation scope, use an invariant or canonical
dimensionless obstruction quantity, converge under refinement, survive both
numerical methods plus gauge/domain/boundary perturbations, and be supported by
an analytic inequality or an independent solver formulation. Otherwise the
strongest permitted label is a stopped or invalid run.

Likewise, `completed_defocusing_candidate` is deliberately not called a
transition. It must still pass DEF1 and ROB1: a finite trapped affine interval,
agreement of the direct and complete Raychaudhuri routes, an error-separated
positive margin, three-resolution convergence, independent-method agreement,
no prior stop, and persistence on a nonzero parameter neighborhood.

The five frozen holdout cases are finite probes, not a proof of an open set.
They may support `completed_no_defocusing` only on their declared finite cases,
time window, and measurement domain. ROB1 may claim a nonzero neighborhood only
after an interval or continuity certificate covers the frozen rational
physical box; its numerical, gauge-cone, dissipation, CFL, and outer-boundary
variations remain separate artifact diagnostics rather than physical axes.

## What a future positive result could say

If RUN1 later authorizes the study and DEF1/ROB1 pass, the strongest permitted
conclusion is:

> The declared classical action admits a controlled spherical interval of
> collapse-induced metric-null defocusing under the tested conditions.

That statement would be a mechanism result. It would not by itself establish
a retained EFT, applicability to nature, global singularity resolution, a
child universe, a dark-sector origin, or varying locally measured light speed.
RUN1 permanently hard-locks `physical_transition_claim_authorized=false`.
Changing that value belongs to a later global authorization that must add an
invariant transition layer, causal continuation, exterior and flux control,
an entropy ledger, and an upper closure or finite handoff; neither RUN1 nor
even a positive spherical DEF1/ROB1 certificate can supply those facts.

## What a future negative result could say

A robust negative may reject the frozen FGC-QR branch on the declared study
domain, or a precisely stated larger mechanism class if the obstruction proof
actually quantifies that class. It cannot reject “gradients,” thermodynamic
phase transitions, variable cones, or nested domains in general merely because
one coupling choice or numerical implementation fails.

This scope is essential to the programme's truth criterion. Failure carries
scientific value only in proportion to the assumptions it genuinely closes.

## Reproduction

Regenerate the current record without rerunning the expensive DOM3 interval
proof:

```bash
python3 scripts/reproduce_fgc_run1_sym1.py \
  --output results/fgc-1-run1-sym1.json
```

The reproducer verifies the hashes binding every inherited result to its direct
configuration, validates the strict PROTO3 schema, the immutable PROTO2/ID1
diagnosis at `fae2c5b…`, and the preserved PROTO1/ID0 lineage,
records absent future
evidence as false, and emits canonical unique-key JSON. A future artifact is
consumable only when its artifact-specific generator and tests establish its
payload and the common reproducer verifies the exact schema, source,
predecessor, implementation, protocol, action, variation, and run-envelope
bindings. ID0's anti-circular lineage check requires the repository's history
to contain checkpoints `d4f0cc8f58408ef4ee12fb32fe3619916e231795` and
`fae2c5b76ba53680f9e70188f79784b160826734`; a shallow
archive lacking that object fails closed rather than treating chronology as a
prose claim. A self-asserted `passed=true` field is never sufficient. The result is
[`results/fgc-1-run1-sym1.json`](../results/fgc-1-run1-sym1.json).

## Fail-closed nonclaims

RUN1 does not define a complete Wilsonian EFT or UV completion. It does not
authorize retained-EFT evolution or a physical transition claim. It derives no
collapse solution, affine-null defocusing interval, invariant transition
surface, singularity resolution, child topology, dark-sector mechanism, or
variation of locally measured `c`.
