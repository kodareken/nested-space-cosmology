# Whole-realization translation control

[The bounded consumer](../src/recursive_horizons/nsc_discovery_translation.py)
reads the saved coupled nf128/nf256 states at $T=0.3$ and $T=1$, resolves
the existing FFT carrier, and translates the entire realization by $s=5.5$
on the period-8 carrier. The translated field is $f_s(x)=f(x-s)$ and
physical observer worldlines move to $x+s$ modulo 8.

Geometry and momenta use integer Fourier modes. Both spinor blocks, the
original source columns and the reference columns use antiperiodic
half-integer modes. Translation by a full period gives minus the spinor,
while its probability and bilinears are periodic. A naive periodic spinor
roll is not this transformation.

The columns of $W$ are translated as spatial basis functions. Their
canonical geometry/momentum coefficients, parent/detail indices and
orthogonality are preserved. Thus the basis interpretation follows the
physical windows; it is not reassigned to an unshifted coordinate cut.
The original state, basis and source objects remain unchanged.

The default offset lies on both saved fermion lattices and all associated
quadrature lattices. Products therefore shift by an exact quadrature-node
permutation, while geometry translation uses its odd Fourier band.
Other offsets for the full-equation control must also lie on the fermion
lattice. Standalone Fourier translation supports arbitrary shifts in its
own band; this control does not claim arbitrary-shift covariance of aliased
nonlinear collocation products.

The child $[1,3]$ moves to $[6.5,8]\cup[0,0.5]$, and the parent $[0,4]$
to $[5.5,8]\cup[0,1.5]$. The diagnostic integrator splits wrapped
intervals before using the owning finite-interpolant integral. Translated
pair interval metadata uses unwrapped endpoints; callers should use this
wrapped diagnostic helper rather than `model.metrics` on those endpoints.
Clock worldlines at $1,2,3$ move to $6.5,7.5,0.5$. Their saved accumulated
proper clocks follow those same physical worldlines without a new clock
integration.

The record compares physical reconstruction and canonical rates, AP
prolongation, full source forces, Hamilton/momentum constraint profiles,
actual projected time jets, actual radial/angular tidal components, $R_h$,
$R_4$ and Weyl scalar profiles. It also compares wrapped proper lengths,
probabilities, signed normal energies, clock rates and source/reference
amplitudes. Constraint covariance is scaled by the matter/geometric term
magnitudes because the residual subtracts large terms. Absolute gaps,
scales and the engineering comparison tolerance are retained in every check.
This comparison does not impose a scientific refinement gate.

## Reproduction

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_translation.py --check
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_translation.py \
  --write results/development/nsc-discovery-translation-v1/translation-shift5.5.json
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_translation.py -q
```

`--check` reads saved inputs and prints an assessment without writing files.
`--write` exclusively creates a new JSON record under the new translation
result directory; existing paths are refused before computation. Input
and consumer/physics hashes are bound, and input hashes are checked again
at the end. Execution is bounded by 30 CPU seconds and a 64 MiB output cap.
`--step 0.0005` optionally compares one owned RK4 step in both charts;
steps above $0.001$ are refused. No trajectory campaign is started.

Passing comparisons support covariance at the saved states when the periodic
seam and coordinate origin move by the stated discrete transformation.
This control does not vary ambient
extent, test delayed echoes, establish a physical regeneration mechanism,
or certify a continuum limit. A passing covariance identity alone makes
no physical-mechanism claim.
