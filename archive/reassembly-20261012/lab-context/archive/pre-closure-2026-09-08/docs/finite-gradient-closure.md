# Finite Gradient Closure

## Status

**Finite Gradient Closure (FGC)** is the provisional organizing law of the
Recursive Horizons programme. It is a hypothesis to implement and test—not a
result established by a black-hole analogy, the closed de Sitter control,
Einstein–Cartan torsion, thermodynamics, or any repository artifact.

Current implementation status belongs to the [claim ledger](claim-ledger.md),
[active code map](active-code-map.md), and [research roadmap](research-roadmap.md).
This document owns the law, its operational meaning, and its failure gates.

## Provisional law

> **Finite Gradient Closure.** Whenever gravitational localization steepens
> beyond a dynamically determined threshold, a healthy higher-order geometric
> or matter-gradient sector activates and dominates the lower-order focusing
> channel before the retained description loses predictivity. The resulting
> finite transition must map finite predecessor state to finite successor state
> while preserving constraints, causal orientation, quasi-local flux, and the
> appropriate generalized-entropy law. Every realized causal domain is finite,
> although their succession need not terminate.

The proposed sequence is

```text
gradient
  -> inward flux and localization
  -> steeper invariant structure
  -> faster-growing restoring response
  -> finite transition
  -> expanding continuation
  -> later upper closure or finite handoff.
```

The infinity, if the mechanism works, belongs only to the absence of a final
iteration. No one realized state may require infinite energy, density,
curvature, bandwidth, entropy, or state data.

## Completed infinity versus unbounded continuation

A **completed physical infinity** is a realized state whose extent or required
description has no finite closure. An **unbounded process** is a sequence with
no final step even though each realized state and transition is finite.

FGC rejects neither mathematical infinity nor useful asymptotic models. It
requires an operational completion only when a physical claim relies on a
singular endpoint or an actually completed infinite state.

“Multiverse” in this programme therefore means a possible derived genealogy of
finite domains—gradients within gradients—not an unconstrained inventory of
disconnected universes.

## Divergence is a diagnostic, not a transition law

A divergent expression may mean:

1. a coordinate or variable choice failed while invariants remain regular;
2. invariant curvature, density, or tidal data diverge and the classical
   solution is incomplete;
3. the retained EFT, gauge, formulation, or numerical method left its validity
   domain; or
4. new finite dynamics supplies a physical handoff.

The divergence alone does not choose among these possibilities. A failed
extension, geodesic incompleteness, or infinite decimal expansion does not
encode a daughter domain.

A candidate handoff surface must supply:

- finite or explicitly controlled limiting state and curvature data;
- a closed map for physical, gauge, and reduction constraints;
- finite quasi-local energy and boundary flux;
- a declared generalized-entropy law;
- causal data sufficient for a unique or probabilistically defined successor;
- a characterization independent of a removable coordinate singularity; and
- evidence that the handoff is not merely lost hyperbolicity or EFT control.

Without those items, “the mathematics reached a wall” is a useful diagnosis,
not evidence that another spacetime was created.

## What a gradient cascade must mean

For a covariant state variable `Phi`, derivative levels include

```text
nabla_mu Phi
nabla_mu nabla_nu Phi
Box Phi.
```

Higher derivatives do not automatically regularize a collapse. They may make
a divergence worse. FGC requires both:

1. the activated response has the restoring sign in the physical focusing
   channel; and
2. it grows faster under localization than the term driving collapse.

The implementation could use an auxiliary regulator, torsion or another
connection degree of freedom, a proved-degenerate higher-gradient theory, or
a controlled quantum transition. No candidate owns the organizing law.

## Finite-scale controls are not gravity theories

The spatial energy

```text
E[Phi] = integral sqrt(h) [
  V(Phi) - alpha (D Phi)^2/2 + beta (D^2 Phi)^2/2
],

alpha > 0, beta > 0
```

has quadratic mode kernel

```text
-alpha k^2 + beta k^4
```

and selects

```text
k_*^2 = alpha/(2 beta),
ell_* = sqrt(2 beta/alpha).
```

This demonstrates finite scale selection only. It does not prove a healthy
time-dependent principal symbol, bounded Hamiltonian, constraint propagation,
bounce, or child spacetime.

A two-sided radius control,

```text
V_eff(R) = -A/R + B/R^n + C R^m,
A,B,C>0, n>1, m>0,
```

can bound finite-energy motion between `R_min` and `R_max`. Its terms still
need a covariant derivation. The exact closed de Sitter solution
`a=L cosh(t/L)` has a finite lower scale but no upper closure and therefore is
not an FGC mechanism.

## Local geometric gate

For an affinely parametrized, irrotational null congruence,

```text
d theta/d lambda
  = -theta^2/2 - sigma_ab sigma^ab - R_ab k^a k^b.
```

Inside a trapped region, a successful metric-null restoring sector must give a
finite interval where

```text
-R_ab k^a k^b > theta^2/2 + sigma_ab sigma^ab,
```

or derive a precise reason the metric-null Raychaudhuri equation is not the
complete causal description in that phase.

A coordinate scale-factor minimum is not a substitute. The current exact
observable contract is only a pointwise spherical measurement prerequisite;
the finite trajectory and error-separated interval remain later gates.

## Activation scale

For a Schwarzschild control,

```text
K = R_abcd R^abcd = 12 r_s^2/R^6.
```

An illustrative trigger `K ~ 1/ell_*^4` gives

```text
R_* ~ 12^(1/6) (r_s ell_*^2)^(1/3).
```

This is an ansatz, not an inheritance law. A physical trigger may require
several dimensionless invariants, expansion, shear, density, and temporal
Hessian data. A regular centre can have zero first spatial gradient while its
density or curvature is maximal.

## Local regularity is not a child universe

Avoiding a local singularity does not select a global outcome. The same core
physics might yield a bouncing star, remnant, oscillating interior,
white-hole-like release, or child cosmology.

A child-domain claim needs one global causal solution containing

```text
theta_+ < 0 and theta_- < 0
  -> finite regular transition
  -> expanding branch with the required anti-trapped region,
```

with a consistent exterior, no unaccounted shell, finite invariants, and a
derived flux law. Topology change and inner-horizon stability carry additional
burdens.

## Energy, flux, and entropy

Gradients store or expose free energy; they do not create unrestricted energy.
The implemented ledger must use constraints, quasi-local mass, boundary work,
and flux because gravitational energy has no unique local tensor in general
relativity.

A schematic nested balance is

```text
dot(E_i) = J_(i-1->i) - J_(i->i+1) + W_i,
```

but every term must be defined on the actual interface. The corresponding
matter/horizon/generalized entropy must obey its derived law. A descendant
cannot exert an ordinary causal pull through a permanent classical event
horizon; any handoff requires a dynamical/global/quantum construction.

An “eternal engine” means only an unbounded succession of finite
transformations. It never means energy creation or perpetual work extraction
from equilibrium.

## Falsifiable requirements

An FGC implementation fails at the claimed level unless it supplies:

1. a covariant action and declared dynamical variables;
2. a calculable invariant activation threshold and nonzero scale;
3. a healthy principal symbol, causal propagation, and no ghost instability;
4. the restoring sign and magnitude in the physical Raychaudhuri/source gate;
5. finite active fields, stress, curvature, density, and torsion invariants;
6. constraint propagation and finite quasi-local mass;
7. a regular trapped-to-expanding global continuation if a child is claimed;
8. a lower closure and upper closure or finite handoff;
9. a complete energy, boundary-flux, and generalized-entropy ledger;
10. perturbative stability against shear, rotation, and inhomogeneity;
11. recovery of tested parent black-hole and child local physics; and
12. an independently fixed observable rather than a fitted identity.

Failure of one action narrows that action’s implementation space. It does not
become evidence for another candidate and does not reject all gradients unless
the proof actually quantifies that wider class.

## Current programme routing

FGC-1 owns action/formulation/health prerequisites. FGC-2-SF1 owns the local
constraint-controlled spherical mechanism test. FGC-3 would own global causal
continuation, exterior and flux closure, entropy, upper handoff, rotation, and
shear.

The live artifact frontier changes too quickly to duplicate here. Use:

- [active code map](active-code-map.md) for the current checkpoint and allowed
  next design;
- [claim ledger](claim-ledger.md) for exact artifact propositions;
- [research roadmap](research-roadmap.md) for the forward dependency graph;
- [FGC runtime matrix](fgc-runtime-matrix.md) for arithmetic and command
  ownership; and
- [results index](../results/README.md) for compact evidence discovery.

The companion [Nested Gradient Spectrum](nested-gradient-spectrum.md) begins
only after a controlled transition supplies a mode/operator/observer map. Dark
components, causal calibration, particle inheritance, complexity, and any
recurrence remain outputs to derive, not premises that can reinterpret an SF1
terminal.
