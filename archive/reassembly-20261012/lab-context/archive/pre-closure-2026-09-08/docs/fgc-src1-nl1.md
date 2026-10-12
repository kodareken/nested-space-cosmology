# FGC-1-SRC1-NL1: nonlinear REF1 source and branch solver

## Decision

`FGC-1-SRC1-NL1` closes one implementation premise for the scoped spherical
mechanism diagnostic:

> The six complete REF1 physical rows can be evaluated in finite binary64
> arithmetic, differentiated with respect to all six coordinate-time
> accelerations through the same tensor code, and solved on the already-proved
> compact QIFT1 branch by a safeguarded, fail-closed Newton method.

This is a source/solver certificate, not an evolution. Its certified domain is
the small annular QIFT1 local branch box at coordinate radius `r=4`. It does not
show that a collapse trajectory remains in that box, include the regular
centre, close the physical constraints, prove a multidirectional health
envelope, authorize the FGC-QR holdout, or establish defocusing.

## One equation implementation

The numerical residual is not a hand-transcribed approximation to the field
equations. The floating adapter supplies the existing `Jet2` interface to the
same function used by REF1:

```text
modified_harmonic_full_residuals
        (
          FloatJet2 state,
          fixed reference connection,
          physical ACT1/VAR1 couplings
        )
             -> six complete REF1 rows
```

`FloatJet2` subclasses the established exact `Jet2` type and implements the
same two-jet arithmetic with finite binary64 scalars. The exact evaluator and
its source files are unchanged. Consequently the floating backend reuses the
same curvature, scalar, and gauge-extension assembly rather than maintaining a
second set of field equations.

This architectural identity is necessary but not sufficient. Floating
arithmetic can still introduce implementation or conditioning error, so the
artifact cross-checks it against exact rational fixtures and an independent
finite-difference derivative.

## Analytic acceleration Jacobian

For fixed

```text
z = (u, p, q, p_r, q_r),
a = p_t,
```

the solver evaluates

```text
R(a; z) = 0.
```

Each acceleration column is seeded with a binary64 first tangent. Six calls
through the same REF1 evaluator produce

```text
J_a = partial R / partial a.
```

The seeded primals must agree within a bound derived from binary64 machine
epsilon. Nonfinite values, a malformed matrix, singular conditioning, or a
primal inconsistency terminates the solve with a typed failure. The analytic
Jacobian is also compared with a central finite-difference Jacobian at a
non-rational control point. Finite differencing is a cross-check only; it is
not the production derivative.

## Safeguarded branch solve

At iteration `n`, the primary step is

```text
J_a(a_n; z) delta_a_n = -R(a_n; z).
```

The implementation then applies a deterministic halving line search. A step is
accepted only if:

- the complete REF1 infinity norm strictly decreases;
- the trial acceleration remains inside the open QIFT1 acceleration box;
- the displacement from the supplied warm start remains below the frozen
  branch-continuity limit; and
- every residual, Jacobian, condition estimate, direction, and trial value is
  finite.

The solver stops before crossing the first failed premise. It never returns the
last iterate as a physical state after failure. Its frozen limits agree with
active `FGC-2-SF1-PROTO3` where that protocol owns the value; those solver
limits are unchanged through PROTO1, PROTO2, and PROTO3:

| Quantity | Bound |
|---|---:|
| Complete REF1 residual infinity norm | `1e-12` |
| Newton iterations | `16` |
| Kinetic-block infinity condition number | `1e10` |
| Warm-start branch displacement | `1/16` |
| QIFT1 parameter half-width | `1/65536` |
| QIFT1 acceleration half-width | `1/128` |

The condition and component norms are local implementation coordinates, not
invariant physical observables. DOM4 now defines the classical run envelope
and HLT1 now defines the canonical scaling and transactional stops, but NUM1
must still validate their integration into a solver before any holdout
evolution is authorized.

## Controls

The machine certificate requires all of the following:

1. **IMP1 flat root.** The floating residual is zero and the floating analytic
   acceleration Jacobian agrees with the exact rational `6 x 6` matrix.
2. **COMP1 at zero acceleration.** The nonzero exact rational residual and
   acceleration Jacobian agree with the floating path.
3. **COMP1 exact root.** The floating path agrees with the exact rational root
   and complete zero residual.
4. **Non-rational derivative control.** A deterministic non-rational local
   state and acceleration compare first-tangent and central finite-difference
   Jacobians under frozen absolute and scaled tolerances.
5. **Warm-start continuation path.** Five predeclared interior points are
   solved in order, with each root seeding the next point; every accepted step
   must strictly reduce the complete residual.
6. **Typed failure injection.** Nonfinite input, parameter-box exit,
   acceleration-box exit, excessive condition number, and iteration exhaustion
   must each return the declared reason before a result can be interpreted.

These controls verify the floating adapter and the safeguarded solver in the
domain already covered by QIFT1. They do not sample the FGC-QR physics holdout.

## RUN1 consequence

SRC1 fills only the `nonlinear_source` evidence slot in
`FGC-1-RUN1-SYM1`. The seven other complete predicates remain false, so
`classical_spherical_diagnostic_authorized`,
`FGCQR_holdout_execution_authorized`, `retained_EFT_evolution_authorized`, and
`physical_transition_claim_authorized` all remain false.

The shared run-envelope hash currently identifies the frozen semantic
`FGC-2-SF1-PROTO2` contract. SRC1 itself certifies only its QIFT1-local
source/solver premise. A later DOM4 result must prove that the actual proposed
run domain exists and is contained in a healthy branch; the shared hash cannot
substitute for that theorem.

## Reproduction

```bash
python3 scripts/reproduce_fgc_src1_nl1.py \
  --output results/fgc-1-src1-nl1.json
python3 scripts/reproduce_fgc_run1_sym1.py \
  --output results/fgc-1-run1-sym1.json
```

The generated SRC1 record binds the exact action, variation, REF1, FO1,
QIFT1, COMP1, RUN1, and protocol identities; the reused evaluator sources; the
floating adapter and solver; this derivation document; and the canonical
configuration.

## Nonclaims

SRC1 derives no physical constraint solution, regular-centre system,
quasilinear initial-boundary problem, multidirectional hyperbolicity theorem,
retained-EFT envelope, collapse solution, finite affine-null defocusing
interval, invariant transition layer, singularity resolution, child domain,
dark-sector mechanism, or varying locally measured speed of light.
