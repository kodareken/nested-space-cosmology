# EC-1: homogeneous Einstein--Cartan spin-fluid collapse cap

EC-1 is a controlled reduced benchmark motivated by Einstein--Cartan gravity.
It implements the homogeneous random-spin ultrarelativistic effective-fluid
reduction used in a published collapse model. It is **not** full
Dirac/tetrad/connection evolution, a smooth vacuum-boundary calculation, a
global black-hole-to-child spacetime, or GMF-1B.

The executable record is
[`results/einstein-cartan-collapse.json`](../results/einstein-cartan-collapse.json),
reproduced with:

```bash
python3 scripts/reproduce_einstein_cartan.py
```

## Exact reduced system

With `c=hbar=1`, `kappa=8 pi G`, and a closed-FLRW cap,


\[
ds^2=-dt^2+a^2(t)
\left[d\chi^2+\sin^2\chi\,d\Omega_2^2\right],
\qquad 0<\chi_b<\frac{\pi}{2},
\]

the benchmark takes

\[
F(a)\equiv\dot a^2=-1+\frac{A}{a^2}-\frac{B}{a^4},
\qquad A>0,\quad B>0,\quad A^2>4B.
\]

The exact turning points in `x=a^2` are

\[
x_\pm=\frac{A\pm\sqrt{A^2-4B}}{2},
\qquad a_{\min}=\sqrt{x_-},\quad a_{\max}=\sqrt{x_+},
\]

The code evaluates the small root as `x_min=B/x_max`, rather than subtracting
nearly equal numbers, and evaluates `F` in its factorized form
`(x_max-a^2)(a^2-x_min)/a^4`. This preserves the finite positive bounce root
in the high-dynamic-range `B << A^2` regime.

and

\[
\ddot a=-\frac{A}{a^3}+\frac{2B}{a^5}.
\]

The effective total density and pressure are

\[
\rho=\frac{3}{\kappa}\left(\frac{A}{a^4}-\frac{B}{a^6}\right),
\qquad
p=\frac{1}{\kappa}\left(\frac{A}{a^4}-\frac{3B}{a^6}\right).
\]

They split into positive radiation (`w=1/3`, `rho proportional to a^-4`)
and a negative spin--torsion term (`w=1`, `rho proportional to -a^-6`). The
reduced no-particle-production system has `Q_internal=0` and exactly obeys

\[
\dot\rho+3H(\rho+p)=0.
\]

The effective null and strong combinations are

\[
\rho+p=\frac{1}{\kappa}\left(\frac{4A}{a^4}-\frac{6B}{a^6}\right),
\qquad
\rho+3p=\frac{1}{\kappa}\left(\frac{6A}{a^4}-\frac{12B}{a^6}\right).
\]

For the benchmark both are negative at `a_min`. That is the effective
null/strong-condition violation and defocusing which permits the reduced
bounce. It is a property of the assumed averaged spin-fluid closure. It does
**not** prove microscopic Dirac stability, an EFT cutoff, or control of the
full torsion/fermion system.

The benchmark does **not** derive a late `w=-1` component, dark energy, or an
external exchange law.

## GMF-1B-PF1: the thermal perfect-fluid closure fails its radial preflight

The homogeneous bounce is not automatically an admissible radial material.
For this particular averaged `T^4-T^6` perfect-fluid closure, define

\[
z=\frac{B}{Aa^2}.
\]

Its equation-of-state ratio and barotropic adiabatic derivative are exactly

\[
w=\frac{1-3z}{3(1-z)},\qquad
c_s^2\equiv\frac{dp}{d\rho}=\frac{2-9z}{3(2-3z)}.
\]

The relevant thresholds are `rho+p=0` and the sound-speed pole at `z=2/3`,
`c_s^2=0` at `z=2/9`, and `rho+3p=0` at `z=1/2`; `rho=0` only at `z=1`.
Thus `c_s^2<0` for `2/9<z<2/3`, while it is greater than one for `z>2/3`.
The lower turning point has

\[
z_{\rm bounce}=\frac{1+\sqrt{1-4B/A^2}}{2},\qquad \frac12<z_{\rm bounce}<1.
\]

Any evolution from the radiation-like large-`a` regime (`z` near zero) to the
bounce therefore enters a negative-sound-speed interval and, when it reaches
`z>2/3`, also a pole followed by a superluminal barotropic derivative. This
is a rejection of this **naive averaged perfect-fluid closure** as a causal,
regular hyperbolic radial GMF-1B matter system. `dp/d rho` is an adiabatic
fluid diagnostic rather than a microscopic characteristic proof; that is why
the conclusion does not reject Einstein--Cartan theory, full EC--Dirac
dynamics, or another action-derived tapered closure. It does say that a
`1+1` solver may not simply promote this homogeneous closure to radial matter.

## Finite reduced bounce and causal signs

On its real branch all displayed invariants are finite:

\[
\mathcal R=\frac{6B}{a^6},
\qquad
R_{\mu\nu}R^{\mu\nu}=\kappa^2(\rho^2+3p^2),
\]

\[
R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}
=12\left[\left(\frac{\ddot a}{a}\right)^2+
\left(\frac{F+1}{a^2}\right)^2\right],
\qquad C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma}=0.
\]

For future radial null normals normalized by `k_plus . k_minus = -2`,

\[
\theta_\pm=2\left(H\pm\frac{\cot\chi}{a}\right).
\]

Marginal radii at a fixed `chi` solve `F=cot^2(chi)`:

\[
a^2_{\rm marg,\pm}=
\frac{A\pm\sqrt{A^2-4B(1+\cot^2\chi)}}
{2(1+\cot^2\chi)}.
\]

For the deterministic benchmark `A=10`, `B=1`, and `chi_b=pi/3`, EC-1
checks one **interior metric** along the sequence

```text
normal at maximum -> marginal -> trapped -> marginal -> normal at bounce
-> marginal -> anti-trapped -> marginal -> normal at maximum.
```

This establishes only a homogeneous-cap causal sequence. It does not match
those spheres to a parent exterior or prove a complete child topology.

## Finite-cap mass/work diagnostic

The cap boundary has

\[
R_b=a\sin\chi_b,
\qquad
M_{\rm MS,b}=
\frac{\sin^3\chi_b}{2G}
\left(\frac{A}{a}-\frac{B}{a^3}\right).
\]

The reduced conservation identity is exactly

\[
\frac{dM_{\rm MS,b}}{dt}+4\pi pR_b^2\dot R_b=0.
\]

The same boundary mass differs between `a_max` and `a_min`. Therefore a
single fixed Schwarzschild/Kottler vacuum mass cannot remain smoothly matched
to this moving pressured cap without a flux and/or resolved boundary layer.
EC-1 calculates that mass drift but supplies neither the layer nor the flux.
For `A=10`, `B=1`, and `chi_b=pi/3`, the executable ratio
`M_MS,bounce/M_MS,a_max` is about `0.1010205`. The record labels the resulting
fixed-parent crossing only as a **comparator**; it is not an event-horizon
calculation or a global matching event.

The half-cycle proper time is integrated without an endpoint singularity by
putting

\[
x=a^2=x_-+(x_+-x_-)\sin^2\vartheta,
\qquad
dt=\sqrt{x(\vartheta)}\,d\vartheta,
\qquad 0\le\vartheta\le\frac{\pi}{2},
\]

and applying refined Simpson quadrature.

## Exact nonclaims

The result record deliberately retains all of the following as `false`:

- `global_solution_constructed`;
- `smooth_static_vacuum_boundary_proven`;
- `child_topology_proven`;
- `derived_external_Q`;
- `late_dark_energy_derived`;
- `variable_c_derived`;
- `radial_stability_proven`.

It sets `homogeneous_spin_fluid_closure_assumed=true`. The benchmark is a
reduced action-motivated fluid closure, not a proof that an Einstein--Cartan
Dirac field has produced this cap, and not evidence for black-hole ancestry.

## GMF-1B successor contract

EC-1 unlocks, but does not substitute for, a 1+1 spherical GMF-1B calculation.
Two smaller source gates now sit between this homogeneous benchmark and that
solver. GMF-1B-ECD-SYM1 identifies a classical spherical two-spinor ansatz with
a surviving minimal-ECD axial-current channel and a metric-null Dirac principal
cone. GMF-1B-ECD-INT1 fixes the canonical-to-numerical signature bridge and
closes the reduced contact-action/direct-cubic interaction identity. Neither
derives the complete ECD stress or solves a gravitational constraint.

The next task, GMF-1B-ECD-ID1, is to derive the **complete causal tapered
closure** before solving the `1+1` equations. It must vary the kinetic,
spin-connection, contact, and tetrad dependence in one convention; derive the
full effective stress and nonlinear radial equations; have a principal
symbol/characteristic analysis; and construct compact support or a taper with
the appropriate radial pressure and torsion source vanishing at its vacuum
boundary. It must then solve horizon-regular Einstein--Cartan constraints with
one common null orientation and a Misner--Sharp unified-first-law residual.
It must show finite Ricci, Ricci-squared, Kretschmann, torsion, and stress
invariants and perform radial convergence and stability checks. If its global
development contains an inner or Cauchy horizon, it must also evolve
counter-streaming perturbations and demonstrate bounded mass/invariants under
refinement; an inner-horizon shortcut is not a stability proof.

The conservative physical target is an expanding closed-FLRW-like **cap**:
`0<=chi<=chi_b` has finite `B^3` spatial topology and a boundary. It is not a
claimed boundaryless `S^3` child. Under smooth connected globally hyperbolic
fixed-topology evolution, a noncompact finite-mass parent Cauchy slice cannot
become a later complete `S^3` Cauchy slice. A future complete-child proposal
must name the changed assumption (for example a failure of global
hyperbolicity, disconnectedness, a degenerate/boundary phase, or a specified
quantum topology change).

Primary sources: [N. Popławski, *Gravitational Collapse with Torsion and
Universe in a Black Hole* (2023)](https://arxiv.org/abs/2307.12190), especially
the closed reduction and its thermal spin-fluid assumptions; and [A. H. Ziaie,
P. V. Moniz, A. Ranjbar, and H. R. Sepangi, *Einstein--Cartan gravitational
collapse of a homogeneous Weyssenhoff fluid* (2014)](https://doi.org/10.1140/epjc/s10052-014-3154-2). The
restricted focusing, topology, and inner-horizon gates are documented in the
[GMF-1 note](global-matching-flux-closure.md) with their primary sources.
