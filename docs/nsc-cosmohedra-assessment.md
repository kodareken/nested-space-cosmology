# Cosmohedra: value, a checked seed and the inheritance question

## Assessment

The subsequent [actual-mode interaction](nsc-atlas-cosmohedron-interaction.md)
now supplies derived vertices and a nonzero selected exchange on the sourced
two-wall geometry, under an explicit probe extension and conditional state.
This advances the low-point seed to a physical mode calculation while
retaining the distinction from the full Einstein wavefunctional.

Cosmohedra are worth a bounded investigation as an interacting-wavefunction
and geometric-composition entry in the atlas. Their connection to the current
quadratic boundary response is not an established identification. The useful
next calculation is a nonlinear coefficient from a declared action and state,
not a forced match of an interacting amplitude to a fixed finite matrix.

Figueiredo's award recognizes important related work and motivates close
reading. Mathematical value for NSC comes from a reproducible connection.
Credit her and her collaborators for their constructions; a citation does
not imply their participation in or endorsement of this project.

The direct sources are [Arkani-Hamed--Figueiredo--Vazao, Cosmohedra](https://arxiv.org/html/2412.19881v1)
and [Ardila-Mantilla et al., Combinatorics of the cosmohedron](https://arxiv.org/html/2603.03425v3).
Read here: the former's action and wavefunction definition, Sections2/4.1--4.2
and Appendix A through the low-point examples; the latter's introduction,
physical motivation and Matryoshka definition. The full fan/polytope proofs
and general loop extraction were not independently reproduced. Grok supplied
an independent comparison against the existing operator/field owners.

## What is actually shared

The cosmohedron organizes coefficients of a specified cubic scalar
wavefunction. Its faces encode nested subpolygons of a momentum polygon,
with factorization and recursion. For a single cut, its smaller geometries
and graph-associahedron factor are specified by equation18 of the original
paper. This is useful mathematical nesting with an explicit observable.

The existing [gravitational entry](nsc-atlas-gravitational-response.md) instead
composes a quadratic axial boundary action through physical interfaces.
Its Schur complement returns a linear response. The
[vortex calculation](nsc-atlas-vortex-transfer.md) retains nonlinear momentum
covariance in a classical fluid initial-value problem. These are distinct
objects to connect; the word nesting does not define the connecting map.

There is a constructive shared ingredient. A Gaussian internal field with a
cubic source J=g_L phi1 phi2+g_R phi3 phi4 can be eliminated to give a term
proportional to J G J. The propagator G and its state/boundary conditions
matter. Evaluating it on external oscillator profiles produces the simplest
cosmohedron graph kernel. A frequency-independent static Schur matrix does
not already contain these time integrations or the wavefunction state.

## An independent low-point check

Take five oscillator coordinates (q1,q2,q3,q4,s), with positive frequencies
(w1,w2,w3,w4,y). The formal Hamiltonian is

$$
H=\frac12\sum_i(-\partial_{q_i}^2+\omega_i^2q_i^2)
+g_Lq_1q_2s+g_Rq_3q_4s.
$$

This is a formal perturbative comparison about the Gaussian sector. A pure
real cubic potential is not bounded below, so this does not assert the
existence of a global normalizable ground state of that polynomial model.

For Psi=exp(-S), put S0=(1/2)sum omega_i q_i^2. The first-order terms are

$$
S_L=\frac{g_Lq_1q_2s}{E_L},\quad
S_R=\frac{g_Rq_3q_4s}{E_R},\quad
E_L=w_1+w_2+y,\quad E_R=w_3+w_4+y.
$$

The mixed connected tree term is

$$
S_{LR}=-g_Lg_Rq_1q_2q_3q_4 K_4,\qquad
K_4=\frac1{E_t E_L E_R},\quad E_t=w_1+w_2+w_3+w_4.
$$

Inserting these terms into the stationary Schrodinger equation and extracting
only g_L g_R gives zero residual exactly. The other second-order sectors
are not asserted to be solved by this truncated ansatz. The coefficient in
log(Psi) is +g_L g_R K4, up to the stated field/coupling conventions.
This independently reproduces the energy denominator of equation49 in the
cosmohedra paper. For frequencies (1,2,3,4,2), the stripped K4 is1/450.
That number is a formal example, not an astronomical prediction.

The [two symbolic controls](../tests/test_nsc_cosmohedra_seed.py) preserve
the Schrodinger residual and the scaling check. They need no PDE grids,
new dependency or state fitted to a target result; allow10 CPU seconds.

## The meaning of T*D=D in this comparison

For a common positive energy scale Omega, the stripped kernel obeys
Omega^3 K4(Omega omega)=K4(omega). That is homogeneity. In the oscillator
family, a normalized Hamiltonian pullback can also be defined with

$$
\omega_i'=\Omega\omega_i,\quad
q_i'\mapsto\Omega^{-1/2}q_i,\quad
p_i'\mapsto\Omega^{1/2}p_i,\quad
g_{L,R}'=\Omega^{5/2}g_{L,R}.
$$

Then Omega^(-1) U* H' U=H formally. Both the cubic term and the mixed
four-point term of S are preserved after the canonical coordinate scaling.
This is a compatible algebraic example of the inherited-law form. It does
not derive the coupling transformation from our spacetime: the actual
mode overlaps and state would have to produce it. Normalized scale covariance
and physical descent retain the distinction already made in the
[finite inherited-law owner](nsc-nested-qualities.md).

Fluid scaling supplies another exact covariance, with its own clock law:
u_a(x,t)=a^-1 u(x/a,t/a^2), p_a=a^-2 p(x/a,t/a^2), fixed viscosity, and
force f_a=a^-3 f(x/a,t/a^2). All Navier--Stokes residual terms scale as
a^-3. The [gravitational response](nsc-atlas-gravitational-response.md) uses
a relativistic length/time scaling instead, preserving dimensionless
reflection when mass/radii scale together and frequency scales inversely.
Identifying these transformations requires a physical clock/state map;
they are not one already identified T.

## Oscillations, horizons and the next calculation

A suitable fluid flow can produce inward or outward wave horizons. Its
normal velocity, wave speed, forcing and medium determine them. A black/white
time reversal additionally changes causal boundary conditions; a viscous
fluid time reversal at fixed positive viscosity needs a changed forcing or
energy history. The common local law and the selected response/state are
separate parts of the inherited description.

Cosmohedra can potentially organize the interaction coefficients of such
fluctuations after their physical action and state have been selected.
Their polygon coordinates are not yet the atlas's horizon radii, flow
velocities or proper clocks. The checked K4 seed makes a bounded pilot
worthwhile; it does not supply a black/white-hole transition or a Lambda value.

The recommended pilot was to derive one cubic coupling between actual modes of
a specified atlas region, compute its connected four-point boundary
coefficient with a declared quantum state, and compare its energy poles
and cut recursion with K4. The [subsequent realization](nsc-atlas-cosmohedron-interaction.md) now supplies
that selected mode/current kernel. A full physical map still needs the
remaining constrained terms, broader modes and a state beyond the conditional
finite domain; the earlier fluid/gravity results retain their own scopes.

```powershell
python scripts/run_python.py -m pytest -q tests/test_nsc_cosmohedra_seed.py
```
