# FGC-1-DOM4-RUN1: classical spherical run domain

## Decision

`FGC-1-DOM4-RUN1` defines and verifies a nonzero classical state domain for
the scoped FGC-QR mechanism diagnostic. It is deliberately not a global tube
containing a future collapse trajectory. It proves the smaller prerequisite
that the protocol coordinates fit inside one declared container and that the
complete REF1 acceleration equations possess strict, branch-continuous,
nonsingular initial-slice witnesses inside it.

The certificate passes only this gate:

```text
nonzero_classical_spherical_run_envelope_passed = true
```

It does not open the FGC-QR holdout by itself. `FGC-1-HYP2-MD1` now closes the
all-covector principal-health prerequisite and `FGC-1-HLT1-MON1` now verifies
the transactional stop formulas on direct and injected controls.
`FGC-1-PRO3-HLD1` must still solve and hash the static premises for the five
resolved inputs, and NUM1 must validate that the future solver calls HLT1
before every accepted stage and stops before leaving the domain.

## Domain that is actually declared

In code units with `L0=4`, DOM4 contains

```text
2 <= A_chi <= 5,
7/4 <= w_chi <= 9/4,
1/2 <= s_phi <= 2,
phi_peak = s_phi/131072,
0 <= r <= 128,
0 <= t <= 32,
0 <= r_measure <= 24.
```

The three parameter axes have exact dimensionless coordinate volume `9/4`.
If the physical seed amplitude is used instead of its dimensionless factor,
the volume is `9/524288`. Every PROTO3 amplitude candidate and every declared
holdout width/seed coordinate lies in this container.

Container membership is not initial-data admission. A particular member must
still satisfy the physical constraints, regular-centre and finite-mass
conditions, the initial compactness window, and absence of initial trapping.
PROTO3 intentionally assigns that all-case decision to PRO3-HLD1. Thus DOM4
does not conceal the fact that some arbitrary corners of the rectangular
coordinate container can miss an initial-premise inequality.

## Complete initial-slice bridge

ID1 supplies `lambda`, the angular extrinsic curvature `k`, and their first
radial derivatives on a constraint-solved maximal polar--areal slice. The
complete REF1 source additionally needs their second radial derivatives and
the six coordinate-time accelerations.

The numerical bridge differentiates the stored first derivatives with a
centred five-point operator. A centred three-point value is retained as an
independent local error indicator. The remaining second-jet entries follow by
differentiating the frozen identities

```text
R=r,
partial_t R=-r k,
partial_t(lambda^2)=4 lambda^2 k,
C^t=C^r=0.
```

The six accelerations are then obtained from the same complete REF1 tensor
residual used by SRC1. Its forward-AD acceleration Jacobian supplies the first
linear predictor; safeguarded Newton acceptance still enforces residual,
condition-number, branch-displacement, acceleration-box, backtracking, and
monotonic-decrease limits. No second implementation of the field equations is
introduced.

## Why the local state set is open

At seed factors `0`, `1`, and `2`, DOM4 continues the centre-support REF1 root
for `A_chi=2`, `w_chi=2`. Each root has a strict residual margin and a finite,
nonsingular acceleration Jacobian. The coefficients of the finite-dimensional
REF1 system are smooth on this regular chart. The ordinary implicit-function
theorem therefore supplies a nonzero open state neighborhood around each
strict witness.

A second stress slice uses `A_chi=5`, `w_chi=9/4`, and seed factor `2`. It is
regular, finite mass, initially untrapped, and inside the protocol compactness
window. Both its profile centre and peak-compactness support points pass the
same complete source/branch checks. This is an adversarial initial-slice
fixture, not evidence that the later trajectory remains healthy.

The seed-zero point remains on the frozen FGC-QR action branch; it is the
zero-regulator start of the continuation, not a relabelled GR evolution.

## What DOM4 consumes and what it does not promote

DOM4 consumes:

- DOM3's compact nonflat radial eigenframe and symmetrizer;
- SRC1's same-evaluator residual, acceleration Jacobian, and typed Newton
  stops;
- CON4's physical/gauge/reduction constraint identities and monitors;
- CTR1's exact regular-centre descendants; and
- ID1's nonzero finite-mass constraint-compatible family.

DOM3 remains radial. DOM4 does not turn it into an all-direction theorem.
HYP2 must independently analyze the full ten-metric-plus-two-scalar principal
system. Likewise, strict initial witnesses plus the implicit-function theorem
do not provide a quasilinear existence time or a global run tube. The later
monitor must evaluate the canonical domain inequalities at every accepted
Runge--Kutta stage and stop before the first crossing.

## Hash boundary

The RUN1 scope binding continues to use PROTO3's semantic holdout-contract hash
as the common declared-run-envelope identity, matching SRC1, CON4, CTR1, and
ID1. DOM4 additionally hashes its narrower domain definition—scope,
parameters, spacetime window, formulation, and source/branch limits—inside
the quantitative evidence. HYP2 must reproduce that exact secondary hash.

## Reproduction

```bash
python3 scripts/reproduce_fgc_dom4_run1.py \
  --output results/fgc-1-dom4-run1.json
```

The result is canonical unique-key JSON and binds its direct configuration,
owner document, implementation, predecessors, ACT1/VAR1 identity, active
PROTO3 semantic hash, and RUN1 config.

## Nonclaims

DOM4 performs no time evolution and reads no FGC-QR holdout outcome. It does
not prove retained-EFT validity, global hyperbolicity, trajectory containment,
collapse, defocusing, a transition, singularity resolution, child topology, a
dark-sector mechanism, or variation of locally measured `c`.
