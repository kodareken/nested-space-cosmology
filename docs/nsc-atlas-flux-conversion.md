# From a flux transition to vacuum and an evolving component

The subsequent [coupled preparation](nsc-atlas-scalar-preparation.md)
solves the necessary scalar layer, gravity and membrane together under
regular Euclidean boundary data. Its continued child state is far less
excited than the rapid transfer assumed here. The present calculation
remains the conditional energy account, not the preparation of that branch.

## The decision and the scope

The [P002 construction](nsc-atlas-flux-transition.md) established a finite
charged geometry. Its immediate successor need not be a full nucleation
calculation: rates cannot create a vacuum absent from its flux sector, and
a constant four-form supplies no additional local propagating field.
The useful connection is whether an explicit transfer can populate an
evolving component while also fixing its residual vacuum.

We compare three routes. P002 formation amplitudes require a finite massive
parent and its state/boundary conditions. Observer fitting inside a constant
flux region tests ordinary Einstein gravity with that effective Lambda.
A scalar/four-form coupling instead supplies an explicit energy-conversion
law. We calculate the last route conditionally, after a short exact test of
the first route's selection capacity.

This adds one scalar and two physical parameters: its total frequency and
mixing ratio. The bare vacuum and displacement flux remain owned inputs.
The initial state is specified by a rapid-transfer assumption, not fitted
independently for each output. SciPy, SymPy and the existing FluxLaw suffice;
no new runtime, PDE engine or scan is introduced. The original P002 solution
and its manuscript remain unchanged.

**Admission.** Borrow the field elimination in
[Kaloper--Sorbo, 0811.1989v2, equations (1)--(3)](https://arxiv.org/html/0811.1989),
then explicitly add a positive quadratic scalar potential. Do not transfer
their axion symmetry or radiative-protection claims to this added term.
Use [Hirano, 1804.09985, equations (2.37), (2.45)](https://arxiv.org/html/1804.09985)
only for isolated zero-mass bubble controls. Calculate eight symbolic
identities, three homogeneous ODE comparisons with at most 30000 RHS calls
each, one specified two-e-fold interval and at most 60 CPU seconds.
The stopping test is energy conservation plus agreement between proper-time
and expansion-time evolution. No observed abundance or SI scale is fitted.

## What transition rates can and cannot select

For fixed charge e>0 the pure four-form has

$$
f_n=f_0-ne,\qquad
\Lambda_n=\Lambda_b+\frac{\kappa}{2}(f_0-ne)^2.
$$

The residue of f_0 modulo e is unchanged by any sequence of signed
integer-charge jumps. Allowing all integer n, the kinematic minimum is

$$
\Lambda_{\min}=\Lambda_b+\frac{\kappa e^2}{2}
\operatorname{dist}(f_0/e,\mathbb Z)^2.
$$

For n>=0 only, minimize over the nonnegative integers; a negative initial
f_0 then has its minimum at n=0. These are lower bounds on a reachable
set, not proofs that dynamics reaches its minimum. Rates can exclude or
reweight members, but cannot create a level outside the set.

In the original normalized benchmark f_0=1.5, e=.5 and Lambda_b=0,
the minimum is zero and the smallest strictly positive value is 1/8.
The open interval (0,1/8) contains no level. Changing the supplied
f_0 to 1.55, with the same charge and bare vacuum, gives a floor 1/800.
These exact values expose dependence on the initial offset. A zero
endpoint in the aligned example is not universal scale selection.
This is the single-flux spacing issue discussed by
[Bousso--Polchinski, Sections 2--3](https://arxiv.org/html/hep-th/0004134),
applied to the actual chosen data.

For the ordinary positive-cap, zero-Kottler-mass bounce, set
k=kappa sigma/2, delta=Lambda_out-Lambda_in and
B=delta/(6k)+k/2. When Lambda_out>Lambda_in>=0 and B-k>0,

$$
R_*^{-2}=\Lambda_{\rm in}/3+B^2,\quad
\beta_{\rm in}=BR_*,\quad \beta_{\rm out}=(B-k)R_*,
$$
$$
S_{\rm bounce}=
\frac{4\pi^2 kR_*^3}
{\kappa(1+\beta_{\rm in})(1+\beta_{\rm out})}.
$$

This stable expression includes the Lambda_in=0 limit. Independent
variation of the off-shell wall-plus-volume action checks its stationary
radius and value.

| Pure-flux step | R_*/L_ref | Normalized bounce action |
| --- | ---: | ---: |
| 1.5 ->1 | .8579969015643 | .6956583643046 |
| 1 ->.5 | 1.4180593547073 | 3.1684103810633 |
| .5 ->0 | 3.8709677419355 | 70.9872066808842 |

The physical exponent multiplies these action coefficients by
L_ref^2/kappa_phys. That factor, determinant prefactors, thin-wall
control and a semiclassical regime are not selected by the classical
benchmark. In particular, these numbers are not P002 transition rates:
P002 has nonzero masses, different initial radii and a finite dust exterior.
The second isolated nucleation radius even exceeds its initial outer wall.
We therefore do not launch a rate campaign by substituting these formulas
into that different geometry.

## A shared potential gives two linked components

Let q be the displacement flux conserved between charged transitions.
With the scalar mixing, the physical four-form is proportional to
q+mu phi. Eliminating it with the proper boundary ensemble gives

$$
U(\phi,q)=\rho_b+\frac{m^2\phi^2}{2}+\frac{(q+\mu\phi)^2}{2},
\quad
\omega^2=m^2+\mu^2,\quad \epsilon=m^2/\omega^2.
$$

The canonical scalar has its usual kinetic term. Completing the square,

$$
\phi_*(q)=-\frac{\mu q}{\omega^2},\qquad
U=\rho_{\rm vac}(q)+\frac{\omega^2}{2}(\phi-\phi_*)^2,\qquad
\rho_{\rm vac}(q)=\rho_b+\frac{\epsilon q^2}{2}.
$$

Assume the parent field is at its minimum and a rapid local transfer changes
q_p to q_c=q_p-e, inheriting continuous phi and zero initial velocity.
The new displacement and excitation are then fixed:

$$
\chi_0=-\frac{\mu e}{\omega^2},\qquad
E_{\rm osc,0}=\frac{1-\epsilon}{2}e^2,\qquad
\frac{E_{\rm osc,0}}{\rho_{\rm vac}(q_c)-\rho_b}
=\frac{1-\epsilon}{\epsilon}\frac{e^2}{q_c^2}.
$$

The ratio requires epsilon>0 and q_c!=0. The actual energy available
to the transfer, per local physical volume in this preparation, is

$$
W_{\rm release}
=\rho_{\rm vac}(q_p)-\rho_{\rm vac}(q_c)-E_{\rm osc,0}
=\epsilon e q_p-\frac{e^2}{2}.
$$

For e>0 no external energy injection is required by this local balance
only if e<=2 epsilon q_p. This is a necessary energy condition for the
specified rapid transfer, not a universal condition on tunnelling. A
negative value cannot be hidden by calling the excitation reheating.
For epsilon=.01, q_p=1.5 and e=.5, the required input is .1175.
For epsilon=0 the relaxed flux vacua all equal rho_b, but the same
frozen-field transfer requires e^2/2 of work.

There is a useful consequence for the ground-state preparation. For
rho_b=0, q_p>0 and 0<epsilon<1/2, the no-injection condition implies
e/q_p<=2epsilon<1. Since r/(1-r) increases for0<r<1,

$$
\frac{E_{\rm osc,0}}{\rho_{\rm vac}(q_c)}
\leq \frac{4\epsilon(1-\epsilon)}{(1-2\epsilon)^2}.
$$

Strong suppression epsilon<<1 therefore does not simultaneously produce
a matter-dominated child by this ground-state rapid transfer. This bound
is specific to the preparation; it is not a rejection of excited parents,
other transition histories or additional sources.

For a parent displacement delta from its minimum and inherited velocity v,
the same law instead gives the state map and energy account

$$
\chi_c=\delta-\frac{\mu e}{\omega^2},\qquad v_c=v,\qquad
E_c=\frac12\left[v^2+\omega^2\chi_c^2\right],
$$
$$
W_{\rm release}=e(\epsilon q_p+\mu\delta-e/2),\qquad
\rho_{\rm parent,total}=\rho_{\rm vac}(q_c)+E_c+W_{\rm release}.
$$

These are two explicit parent-state inputs, not extra couplings or a
reset of the child state. Their values must come from the parent's
actual history or stated calibration before claiming a prediction.

**The actual boundary force.** With no direct scalar coupling on the
membrane, finite-energy scalar matching requires continuous phi and
normal derivative. Tangential derivatives then agree as well. The
kinetic pieces of T_nn cancel across the boundary, leaving

$$
[T_{nn}]=-[U]=-e(q_{\rm out}+\mu\phi_\Sigma-e/2),\qquad
\sigma\overline K=e(q_{\rm out}+\mu\phi_\Sigma-e/2).
$$

The tangential flux jump also vanishes. Thus the earlier normal
conservation law extends with the *actual* wall field, not separately
relaxed vacuum values. At the inherited parent minimum this force is
exactly W_release. Replacing it by the relaxed vacuum drop would
overestimate it by E_osc,0. This is the source-to-boundary connection
which a full evolving solution must retain.

This local pressure identity is not an already integrated membrane
energy or a solved moving wall. The boundary must provide the work and
scalar gradients; an instantaneous change of q on a whole spacelike
slice is not the causal motion of a charged wall.

## One conditional homogeneous evolution

As a separate post-transfer Cauchy problem, supply a homogeneous flat
region with chi=chi_0, chi_dot=0 and expanding H fixed by its energy.
No intervening wall, spatial curvature, thermal bath or density
perturbations are imported from P002. For constant child q,

$$
\ddot\chi+3H\dot\chi+\omega^2\chi=0,\qquad
3H^2=\kappa\left[\rho_{\rm vac}
+\tfrac12(\dot\chi^2+\omega^2\chi^2)\right],
\qquad
\dot H=-\frac{\kappa}{2}\dot\chi^2.
$$

The exact conservation law is rho_dot+3H chi_dot^2=0.
The Friedmann constraint propagates. For the oscillating component,

$$
p_{\rm osc}=\tfrac12(\dot\chi^2-\omega^2\chi^2),\qquad
\frac{d}{dt}(a^3\chi\dot\chi)=2a^3p_{\rm osc}.
$$

A full cycle whose endpoints have chi chi_dot=0 has zero integrated
comoving pressure. This is not exactly zero instantaneous pressure or
exact conservation of rho_osc a^3. When H varies little over a cycle,
it supplies the pressureless-matter approximation.

Choose epsilon=1/4, omega=20, q_p=1.5, e=.5 and rho_b=0 in normalized
variables. The exact initial account is

| Quantity | Normalized value |
| --- | ---: |
| Parent vacuum energy | .28125 |
| Child residual vacuum | .125 |
| Initial coherent excitation | .09375 |
| Energy released to the transfer | .0625 |
| Excitation / residual vacuum | .75 |

All four energies follow one potential and the supplied preparation.
After two expansion e-folds, a=e^2=7.38905610,

$$
t=9.31040100558,\quad H=.204316778764,\quad
\rho_{\rm osc}=.000236038253365,\quad
\frac{\rho_{\rm osc}a^3}{E_{\rm osc,0}}=1.01572936291.
$$

The endpoint differs from ideal dust dilution by 1.573%. The scalar was
evolved directly, not replaced by dust in the equations.
The reference starts with less oscillation energy than vacuum energy;
it does not contain a matter-dominated era.
Its initial omega/H is 74.07. This has not established cold-dark-matter
clustering, an observational abundance, reheating or ordinary matter.
The added mass and mixing are specified theory parameters, not derived
from geometry, and the bare vacuum remains unselected.

Here phi_hat=sqrt(kappa_phys) phi, q_hat=sqrt(kappa_phys)L_ref q,
omega_hat=L_ref omega and rho_hat=kappa_phys L_ref^2 rho.
Arithmetic kappa=1 therefore selects no physical SI ruler.

## Evidence and the next physical question

The [compact record](../results/summaries/nsc-atlas-flux-conversion-v1.json)
binds the producer and tests. Eight symbolic identities vanish exactly.
Proper-time integration evolves H with Raychaudhuri; an independent
expansion-time integration enforces Friedmann algebraically. Endpoint
clock/field differences are below 2.6e-12; maximum Friedmann and integrated
energy residuals are below 5.8e-13. Two tolerances differ by 2.7e-10 or less.
Those checks are numerical diagnostics, not observational error bars.

```
python scripts/run_python.py scripts/check_nsc_atlas_flux_conversion.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_flux_conversion.py -q
```

The priority is now the physical preparation, including scalar continuity,
normal stress and boundary work, followed by perturbations of this same
source. The homogeneous calculation alone does not close a parent-child
spacetime or justify an observational prediction-register entry.
It identifies the energy relation that the next construction must preserve.
Do not silently insert these numbers into P002 or adjust epsilon separately
for vacuum and clustering.
