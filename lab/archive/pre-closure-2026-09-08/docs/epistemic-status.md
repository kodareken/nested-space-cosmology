# Epistemic status

This project separates a proposition’s meaning from its evidentiary status. A
beautiful idea can be unproved; a negative result can be useful; an established
equation can appear inside an unsupported inference.

Current claim assignments are owned by the [claim ledger](claim-ledger.md).
This document owns the labels and review rules, not the artifact chronology.

## Labels

### E — Established input

A result supported by a mature theory in its stated domain, a standard
derivation, and/or replicated observation. “Established” never means
unrestricted. An E claim states its domain and primary source.

Examples include local null propagation in a Lorentzian metric, the
Bekenstein–Hawking entropy of a stationary semiclassical horizon, and the
Kantowski–Sachs form of the Schwarzschild interior.

### O — Observation or published constraint

A result reported by an experiment or survey, with its release, estimator,
model, and uncertainty. Observational inferences such as CMB-derived `H0` remain
conditional on their cosmological model. An O label does not identify a unique
microscopic ontology.

### D — Repository derivation

A consequence reproduced from stated premises:

- **D-I: identity** — follows from definitions or algebraic substitution;
- **D-C: conditional consequence** — follows if explicit model postulates hold;
- **D-N: numerical model result** — output of a specified numerical system
  with convergence and error evidence.

A D claim is not an independent prediction when the target data fixed its
inputs. A file hash establishes byte identity, not the physical conclusion.

Public-data calculations name their layer: raw reanalysis, compressed
likelihood, posterior comparison, metadata reconstruction, or arithmetic
reproduction. Completing a subproblem does not promote the enclosing project.

### H — Hypothesis or model postulate

A new physical proposition introduced by Recursive Horizons. It may be
motivated and compatible with known physics but is not established by those
ingredients. Every H claim needs a derivation or falsification route.

### A — Analogy or interpretation

A picture useful for reasoning but not itself an action, stress tensor,
solution, channel, or likelihood. Pressure cookers, trampolines, resonance,
simulation language, “sink holes,” and philosophical or theological readings
belong here until translated into operational variables.

### N — Negative, null, or invalidated result

A declared route that failed, an advertised number contradicted by its own
calculation, or an invalid estimator. Negative results are retained and are not
renamed “insensitive” merely to protect a preferred conclusion.

### R — Retired claim

A proposition removed because it lacks a mechanism, conflicts with established
physics, depends on invalid analysis, or is unnecessary. Revival requires a new
derivation and evidence record.

### Q — Open question

A required piece not yet derived or measured. Q is the work programme, not a
weakness to hide.

## Organizing law versus implementation

Finite Gradient Closure is an **H/Q organizing law**. It orders calculations
and failure conditions; it is not raised to D by a successful toy control or
local candidate certificate.

Artifact IDs do not transfer evidence upward. A pointwise principal-symbol
proof, source identity, exact local root, finite initial-data family, numerical
engine validation, or campaign binder establishes only its stated subproblem.
The enclosing action still needs the remaining constraints, health, global
solution, error, robustness, and observational gates.

No candidate owns FGC. Failure of Einstein–Cartan, one scalar action, one
temporal method, or one numerical formulation narrows that branch. It does not
become success for another candidate and does not reject FGC unless the result
actually quantifies the wider class.

The current artifact-level assignments are intentionally absent from this
guide. Use the [claim ledger](claim-ledger.md) and
[active code map](active-code-map.md).

## Prediction, inference, calibration, and consistency

These words are not interchangeable:

- **Prediction:** parameters are fixed independently of the target data and the
  result could have disagreed.
- **Inference:** target data estimate a model parameter.
- **Calibration:** known data deliberately adjust a parameter or transfer.
- **Consistency check:** independently specified quantities are compared with
  uncertainties and covariance.
- **Identity check:** software confirms an algebraic relationship already
  guaranteed by definitions.
- **Sensitivity analysis:** published inputs vary to show how a conditional
  result moves.

For example, `Lambda=3/r_s^2` is D-C after the optional saturation postulate.
Solving it for `r_s` using measured `Lambda` is inference. Substitution back
into the same equation is an identity check, not a prediction.

## Local, coordinate, global, and cross-domain statements

The repository distinguishes:

- **local invariant:** measured by a local orthonormal observer or scalar;
- **coordinate quantity:** chart dependent, such as Schwarzschild `dr/dt`;
- **global causal statement:** reachability, horizon structure, topology, or
  geodesic completeness; and
- **cross-domain statement:** requires a physical transition map and a
  dimensionless comparison convention.

A Schwarzschild-coordinate speed approaching zero is not a locally measured
change in `c`. “Time ends” must name proper time, geodesic completeness,
conformal boundary, phase boundary, or finite handoff. A finite conformal
diagram does not imply instantaneous information transfer.

A singular endpoint can be a wall in current knowledge. It is not evidence of
a demonstrated transition or destination. Repeated decimal digits are not
physical labels for nested domains.

## Physical settings, units, and mathematics

- A **dimensionful quantity** changes numerical value with units; a variation
  claim must become an operational ratio of clocks, rods, masses, or cones.
- A **dimensionless physical setting** can be compared once its measurement
  and dynamics are defined.
- An **effective parameter** may depend on state or scale; calling it emergent
  does not prove observable spacetime variation.
- A **mathematical constant** is not a material coupling. A new physical
  structure may map mathematics differently, but `pi` and `e` do not vary as
  domain substances.

Temperature and Planck length are different categories. `T=0` does not mean
zero energy or no quantum fluctuations. Cross-domain frequency or causal
claims require observer-defined proper frequency and dimensionless reference
ratios.

## Background equivalence is not physical identification

Matching `a(t)`, `H(z)`, or an equation of state does not identify fields or
mechanisms. A stationary scalar offset can be vacuum-like without being the
observed dark energy. A rapidly oscillating scalar can average matter-like
without satisfying abundance, clustering, lensing, structure, and longevity.
An effective transfer `Q` is physical only after its stress tensors and flux
law are derived.

The same applies to particle ontology. QFT already describes excitations by
states, poles, resonances, correlators, amplitudes, and detector events. A
deeper inherited-band interpretation must still recover Standard Model masses,
spin/statistics, charges, unitarity, Lorentz behavior, and measured rates and
must predict a distinguishing feature.

## Open-minded without madness

Dark matter and dark energy name constrained gravitational and expansion
phenomenology; the names do not uniquely settle microscopic ontology. The
phenomena are not optional. A replacement must confront the joint ΛCDM-plus-GR
evidence package rather than one residual.

The canonical rule is **no unaccountable explanatory variable**, not “no
unobserved variable.” For every proposed degree of freedom, record:

```text
status and units
invariant observable or transition map
symmetry/gauge redundancy and prior domain
calibration inputs and held-out target
falsifier and removal condition.
```

If two descriptions predict the same distributions for every declared
experiment, current evidence does not select between their ontologies.

Antimatter is kept separate because positrons, antiprotons, and antihydrogen
are directly produced and measured. The open questions concern cosmic
matter–antimatter asymmetry and higher-precision CPT/equivalence tests.

“Relational” is used operationally: observables need invariant quantities or
explicit clocks, rods, detectors, and transition maps. No quantum
interpretation supplies missing FGC dynamics.

## A source certificate is not a spacetime or quantum state

A local source, interaction, torsion, constraint, or principal-symbol identity
can close one decisive question without being a solution. It does not by
itself supply complete stress, regular-centre data, exterior, curvature
evolution, Raychaudhuri history, topology, stability, or a quantum many-body
state.

Classical c-number spinors are not automatically a fermion gas or quantum
singlet. Signature/gamma translations must carry tetrad, adjoint, gamma-five,
orientation, index lowering, action, curvature, and stress conventions
together.

## Review rule

When a reviewer identifies a conflict, ask:

1. Which label was wrong?
2. Which premise generated the conclusion?
3. What survives if that premise is removed?
4. What equation or observation distinguishes the revised model?
5. Is the output dimensionless and operational?
6. Is it a local, coordinate, global, or cross-domain statement?
7. Is the data layer named correctly?
8. Were parameters fixed before the target data?
9. Does the new degree of freedom have a measurement map and falsifier?
10. Does a source certificate actually derive the complete source and
    evolution being claimed?
11. Does the restoring sector have the required sign and scaling before health
    or EFT failure?
12. Are lower and upper closure, quasi-local flux, and entropy all addressed?
13. Has an unbounded finite process been confused with a completed infinity or
    perpetual motion?
14. Does a dark-sector claim recover background, perturbations, lensing,
    structure, early-universe, and local evidence?
15. Does an alleged robust result survive the declared methods, resolutions,
    gauge/domain/boundary diagnostics, and a nonzero neighborhood?

The purpose is not to defend a sentence. It is to leave the surviving map of
reality more accurate than before.

## Current routing

The exact current scientific boundary is intentionally maintained elsewhere:

- [claim ledger](claim-ledger.md) for proposition status;
- [active code map](active-code-map.md) for the current checkpoint;
- [research roadmap](research-roadmap.md) for dependencies and finish line;
- [runtime matrix](fgc-runtime-matrix.md) for execution ownership; and
- [dated review](../archive/reviews/review-summary-2026-08-21.md) for the
  historical corpus audit.

A status example in Git history remains evidence of its own boundary. It is not
kept current by appending every later artifact to this methodological guide.
