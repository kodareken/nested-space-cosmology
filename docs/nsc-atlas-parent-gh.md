# Nonlinear evolution of the nonspherical inherited parent

The [finite initial-data owner](nsc-atlas-finite-parent-data.md) supplies
the unchanged real scalar, internal vacuum and positive constrained
geometry. A full spatially varying Einstein--two-scalar evolution
implementation now exists in the same package. Its finite-amplitude
physical results remain under qualification; a completed coarse run
with large Einstein residuals is not admitted as a prediction.

## The implemented connection

The evolution uses the covariant second-order GH equation24 of
[Lindblom--Szilagyi--Taylor](https://arxiv.org/html/1312.0701v2), with
static round reference geometry and the canonical fields
S=sqrt(84)(s-s_star), chi. The potential is
W=V(s)+mu^2(s)chi^2/2, including the constant vacuum once.
All metric components and both scalar waves evolve. Constraint damping
vanishes on an Einstein solution; no physical mode is suppressed.

The global unit-S3 frame has [e_i,e_j]=2 epsilon_ijk e_k. Real normalized
Wigner harmonics provide sparse exact coefficient derivative generators,
with Casimir sum D_i^2=-l(l+2). The physical reference derivatives have
the additional factor1/a0. Covariant tensor derivatives retain the frame
connection and curvature; raw ordered second derivatives are not
symmetrized as though the frame commuted.

The preserved proper paired-flip group gives diagonal characters for
each tensor component and scalar. It reduces computation without a
spherical average. The larger evolution sector retains modes absent
from the kernel-free initial-constraint symmetry sector.

Hopf coordinates separate each complex Wigner mode into a radial
factor and two circle phases. Two Fourier transforms and radial
quadrature evaluate/project it without the dense values matrix.
The fast and dense representations agree on independent controls.
This is computational use of circle fibers, not an assumed physical
velocity field or a particle interpretation.

## Implementation owners

| Role | Component |
| --- | --- |
| Real modes and regular frame | [s3_frame](../src/recursive_horizons/nsc_atlas_s3_frame.py) |
| Fast Hopf transform | [s3_hopf](../src/recursive_horizons/nsc_atlas_s3_hopf.py) |
| Parallel local GH/Klein--Gordon equations | [gh_pointwise](../src/recursive_horizons/nsc_atlas_gh_pointwise.py) |
| Independent physical ADM constraints | [gh_constraints](../src/recursive_horizons/nsc_atlas_gh_constraints.py) |
| Nonlinear driver and geodesic observer | [parent_gh](../src/recursive_horizons/nsc_atlas_parent_gh.py) |
| Phase/amplitude-gradient diagnostics | [gh_wave_diagnostics](../src/recursive_horizons/nsc_atlas_gh_wave_diagnostics.py) |

The ADM calculator derives physical spatial connection and curvature
directly from Koszul/Riemann formulas, independently of the GH kernel.
The driver checks Hamiltonian, momentum and GH constraints, accepted
metric/lapse validity, spectral-shell power and unresolved accelerations.

## Coordinates and the clock protocol

An observer begins at the original scalar maximum, initially normal to
the inherited slice. Its geodesic, proper duration and field derivative
are evolved with the geometry. Phase is read from the resolved field's
Z=chi-i u(chi)/omega_ref; it does not drive the PDE time-step controller.
The current sampled unwrapping flag is a heuristic, not a certified
winding count. A wrapped increment cannot exclude an extra full turn;
accepted physical comparisons must use directly reconstructed phase or
a separately resolved counting protocol. Exact interference nodes are
not spacetime singularities.

The original prescribed vacuum clock H0=3H tanh(Ht) has a demonstrable
coordinate limitation. For homogeneous lapse N and scale a,

$$
\mathcal C_0=H_0+\partial_t\ln N-3\partial_t\ln a=0,
\qquad N=(a/a_{\rm init})^3\operatorname{sech}^3(Ht).
$$

Thus

$$
\int_0^{\tau(t)}(a_{\rm init}/a)^3d\tau
=\frac{\operatorname{sech}(Ht)\tanh(Ht)+\arctan(\sinh Ht)}{2H}
\leq\frac{\pi}{4H}.
$$

A regular static geometry has a finite proper-time window in that chart.
The zero-H harmonic gauge instead gives
t=integral (a_init/a)^3 d tau, removing that imposed proper-time ceiling.
For de Sitter, N=cosh(H tau)^3 and the finite target proper time remains
reachable. These are coordinate statements, not a change of action or
a conclusion about the finite source's fate.

## Qualification and current limits

The local kernel has independent de Sitter, nonlinear flat-coordinate,
anisotropic and matter controls. Frame commutators, orthonormality,
Laplacian spectrum, fast-transform normalization and derivative ordering
are tested separately. Integrated vacuum evolution is compared in proper
duration, not by equating different gauge-coordinate times.

The nonlinear driver evolves a spectral projection of Psi^4 and its gauge
data. It does not preserve the original conformal initial constraint
exactly at finite degree. Initial as well as final ADM/GH errors must
shrink. The full-amplitude coarse diagnostic was underresolved: its late
shells and unresolved acceleration contained substantial power. Higher
resolution and phase/amplitude diagnosis are therefore required before
claiming physical concentration, trapping or a clock prediction.

The gradient diagnostic splits |DZ|^2 into amplitude and phase-gradient
contributions where Z is nonzero. It reports undefined nodes and their
excluded contribution separately. Those are quadrature-gradient powers,
not conserved matter-energy fractions or baryon charge.

The phase and amplitude both contribute in the coarse diagnostic. A
carrier-only or homogeneous replacement is therefore not silently used
to obtain a more convenient finite result. The active target remains
nonlinear evolution of the original nonspherical finite preparation.

## Running and preserving work

```powershell
python scripts/run_python.py scripts/evolve_nsc_atlas_parent_gh.py --degree 12 --epsilon 1 --threads 8 --budget-seconds 1500 --state-cache .build/gh-L12-state.npz --output .build/gh-L12-summary.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_s3_frame.py tests/test_nsc_atlas_s3_hopf.py tests/test_nsc_atlas_gh_pointwise.py tests/test_nsc_atlas_gh_constraints.py tests/test_nsc_atlas_parent_gh.py tests/test_nsc_atlas_gh_wave_diagnostics.py
```

Small checkpoints stay under .build for diagnosis; compact source-bound
records belong in the result archive, including failed or unresolved
comparisons under their explicit status. A source hash change during a run
is explicitly flagged. Incomplete endpoints retain their status and
cannot satisfy the finite-evolution target. The current action, physical
source, SI-scale gap and charge-selection gap remain unchanged.
