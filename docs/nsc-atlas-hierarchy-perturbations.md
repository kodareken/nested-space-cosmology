# Radial breathing and quadrupole response on the same evolving core

The [finite hierarchy evolution](nsc-atlas-hierarchy-evolution.md) retained
a localized core while its parent expanded, but its radius oscillations
increased. This connection tests how neighbouring states respond, including
departures from spherical symmetry. It transports perturbations on that
actual evolving source, with the reservoir and parent expansion retained.

The useful claim is a finite response bound for a declared family, not
asymptotic stability for every possible perturbation. A field-phase change
and a density deformation are different outputs and are measured separately.

## Calculation and mathematical inputs

The choices were extending the radial trajectory, solving a stationary
eigenproblem, and linearizing the actual nonstationary trajectory. Select
the third: a stationary solution or a longer radius trace would not test
the present time-dependent source's response to angular distortions.
The source stepper is imported unchanged and its final field hashes must
match the earlier results exactly at every retained resolution.

The new state is a small set of tangent fields. Each supplies its own
linearized gravitational potential; it does not alter the base solution.
The implementation uses the existing spectral radial kinetic step for l=0
and a unitary centrifugal radial step for l=2. Both are checked against
independent analytic calculations. Five resolution/time/domain cases and
two nonlinear radial derivative checks fit a 4500 CPU-second budget, with
independent runs parallelized and completed fields retained only locally.

[Harrison, Moroz and Tod](https://arxiv.org/abs/math-ph/0208045) study
stationary Schrodinger--Newton stability; that stationary conclusion is
not assumed for this expanding, breathing source.
[Zagorac et al.](https://arxiv.org/abs/2109.01920) study soliton breathing,
quadrupolar distortions and time-dependent perturbations. Their abstracts
and bibliographic records were checked. The present calculation derives
and evaluates the response on the already retained NSC preparation rather
than substituting another profile or potential.

## Real-linear perturbation equations

Let psi0(x,tau) and V0(x,tau) be the retained total envelope and potential.
Use real angular functions with angular mean square one; for an axisymmetric
quadrupole Y20=sqrt(5) P2. Write

\[
\psi=\psi_0+\varepsilon\eta_lY_{lm}+O(\varepsilon^2),\qquad
V=V_0+\varepsilon v_lY_{lm}+O(\varepsilon^2).
\]

The linearized equations are

\[
i\partial_\tau\eta_l=
\left[-\frac{D_l}{2A^2}+V_0\right]\eta_l+v_l\psi_0,\qquad
D_lv_l=\frac{\operatorname{Re}(\psi_0^*\eta_l)}{A},
\]

\[
D_l=\partial_x^2+\frac2x\partial_x-\frac{l(l+1)}{x^2}.
\]

The response is linear over real coefficients, not over arbitrary complex
coefficients: the gravitational source contains the real density variation
d_l=2 Re(psi0* eta_l). Phase directions must therefore be included explicitly
in the initial real basis.

Regularity at the origin and decaying exterior potential give

\[
v_l(x)=-\frac{1}{2A(2l+1)}\left[
\frac{1}{x^{l+1}}\int_0^x s^{l+2}d_l(s)ds
+x^l\int_x^{L}s^{1-l}d_l(s)ds\right].
\]

The finite domain has zero additional source outside L. The field boundary
is homogeneous, as in the base experiment; larger-domain comparison tests
its influence. The l=2 coefficient has no monopole mass variation because
its angular mean is zero. Its radial coefficient integral must not be
misreported as a total mass change.

## Exact tangent potential flow and a conserved bilinear form

During a potential-only substep of duration h, the base density is fixed.
Its first variation is also fixed: differentiating 2 Re(psi0* eta) along
the joint base/tangent potential flow gives zero. The base field itself
rotates; freezing that field would define a different linear equation.
Consequently the exact tangent substep is

\[
\eta_l^{\rm new}=e^{-ihV_0}
\left(\eta_l-ihv_l\psi_0\right).
\]

In q=x eta variables the kinetic operator is
-d_xx/(2A^2)+l(l+1)/(2A^2x^2). The l=0 step differentiates the exact same
discrete kinetic flow as the base. The l=2 Crank--Nicolson step is unitary;
a free solid-harmonic Gaussian checks its centrifugal coefficient, regular
origin and refinement independently.

For two tangent solutions, the real symplectic form

\[
S_{jk}=\operatorname{Im}\int q_j^*q_k\,dx
\]

is conserved even though A and the base state vary in time. The kinetic
and V0 terms cancel by self-adjointness. The remaining derivative is a
difference of two density/potential pairings, which cancels because the
multipole Green kernel is symmetric. This is an exact identity of the
linearized equations; conservation does not by itself imply stability.
The code monitors it, and a passive quadrupole with v_l suppressed checks
unitarity separately.

## A response bound for all coefficients in the declared family

Choose two smooth compact radial shapes proportional to
(A_seed*x)^l exp[-(A_seed*x)^2/(2w^2)], with w=2 and 4 and a taper from
physical scaled radius 15 to 20. Include each shape and its imaginary phase
partner. The radial family is projected to preserve total mass to first
order; the quadrupole already has zero angular mean. Real Gram
orthonormalization supplies four independent canonical directions q_j.

For real coefficients c with sum c_j^2=1, write a perturbation of fractional
initial core-field norm epsilon as

\[
\delta q(\tau)=\varepsilon N_c\sum_{j=1}^4c_jq_j(\tau),
\qquad N_c=\|u_{\rm core}(\tau_{\rm seed})\|_2.
\]

The real field Gram matrix G_jk=Re integral q_j* q_k then proves

\[
\frac{\|\delta q(\tau)\|_2}{\varepsilon N_c}
\leq\sqrt{\lambda_{\max}(G(\tau))}.
\]

For density, form the real matrix from d_j=2 Re(psi0* q_j/x) inside the
fixed physical core cutoff. Divide its largest eigenvalue's square root
by the base excess-density norm and multiply by Nc. This gives the maximum
fractional density response per fractional initial field perturbation.
This coefficient is a density response per initial field-norm perturbation,
not a density-to-density growth factor. It includes initially imaginary
directions whose initial density signal
is zero; no division by that zero is made.

Spherical symmetry makes the radial propagator identical for all five
quadrupole components. Angular orthogonality gives the direct-sum Gram
matrix I5 tensor G. Its largest eigenvalue is the same one. Thus the
four-column calculation bounds every combination in the corresponding
twenty-dimensional real quadrupole span, not just one chosen orientation.
It does not cover every radial profile, higher angular multipoles or
nonlinear perturbations of finite size.

The source is sampled each unit of slow time. Numerical maxima refer to
those saved events, not a rigorous continuous-time enclosure. A second
quadrupole family is initialized on the evolved state at tau 480 and followed
to 600. This tests a late breathing phase independently of the initial seeds.

## Physical deformation and an independent nonlinear check

The normalized second-moment quadrupole uses an additional x^2 weight.
For the axisymmetric component its first-order coefficient is

\[
\frac{\int R^2P_2\,\delta\rho\,d^3R}
{\int R^2\rho_{\rm excess}\,d^3R}
=\frac{\varepsilon N_c}{\sqrt5}
\frac{\sum_jc_j\int x^4d_jdx}{\int x^4d_0dx}.
\]

The angular average of P2*Y20 is 1/sqrt(5), which fixes the factor here.
This integrated P2 moment is different from the coefficient of P2 or a
Cartesian trace-free tensor component. Direct angular quadrature of the
full perturbed density independently checks the normalization.
Its maximum follows directly from the norm of the four response coefficients.
This is a core-shape diagnostic, not an already computed lens image.
The density response of the direction with largest field-norm growth is
also retained, preventing phase sensitivity from being silently called
geometric disruption.

The worst final radial density direction is checked by two independent
centred differences of the original nonlinear spherical envelope flow,
using initial amplitudes1e-4 and5e-5. Convergence toward the tangent result
checks the full radial linearization. A separate multipole Green-function
test uses an analytic radial Laplacian, and the free centrifugal propagator
is checked against a solid-harmonic Gaussian. The nonlinear radial check
does not substitute for a nonlinear three-dimensional evolution.

## Calculated finite response

The [compact response record](../results/summaries/nsc-atlas-hierarchy-perturbations-v1.json)
contains all five controls and the nonlinear radial checks. Every spherical
source replay matches the preceding retained field hash exactly.

| Family | Largest sampled field-norm amplification | Largest sampled density response per initial field perturbation |
| --- | ---: | ---: |
| Radial, seeded at 0 | 6.756 +/- .162 | 3.592 +/- .033 |
| Quadrupole, seeded at 0 | 1.08113 +/- .00025 | 1.15675 +/- .00007 |
| Quadrupole, seeded at 480 | 1.0820 +/- .0011 | 1.1720 +/- .0049 |

Thus a 1% fractional field perturbation in the quadrupole span seeded at
time zero has a calculated first-order core-density response of at most about
1.157% at the saved events, with the stated numerical comparison range.
This is not a nonlinear guarantee for a finite 1% perturbation, nor a
rigorous upper bound at unsampled times. The maxima for different outputs
need not be reached by the same initial coefficient vector.

The radial response is more sensitive to refinement. At the final event,
for the direction with largest field norm, 99.4% of its squared norm lies
inside the core. Of that core part, 99.5% is locally phase-like relative to
psi0. This is a measured decomposition of that direction, not a statement
that all radial disturbances are phase changes. Spatial phase gradients
can affect later flow. The radial density response remains a separate,
nonzero quantity, and asymptotic radial stability is not established.

The maximum symplectic discrepancy is 1.15e-11 in the radial channel and
9.76e-14 in the quadrupole. The passive quadrupole norm stays at one within
5e-14. For the independent nonlinear radial flow, halving the initial
amplitude from 1e-4 to 5e-5 reduces the complex-field derivative error from
3.29e-6 to 8.24e-7 and the density derivative error from 2.18e-6 to 5.46e-7.
This is the expected centred-difference ordering. The underlying campaign
used about 4073 CPU seconds within its 4500-second budget.

These results provide finite evidence against rapid quadrupole disruption
in this tested family, despite the base core's growing breathing. They do
not classify every radial profile, reservoir-supported seed, higher angular
multipole, nonlinear configuration or full relativistic solution.

## The exterior probe exposed a boundary condition, not a new prediction

A further response was evaluated at fixed physical radius 2Rd, outside the
initial reservoir. The full multipole kernel was used, including any
perturbed matter outside that radius. For initial coefficient norm one,
Cauchy--Schwarz bounds the angular-RMS potential by the norm of its four
response coefficients. The addition identity sum_m Y_2m(n)^2=5 gives a
pointwise all-orientation bound with an additional sqrt(5).

The original Dirichlet-field probe is **unresolved**: changing the outer
domain at fixed spacing changes its coefficient from 3.33e-14 to 9.28e-15.
That datum must not be reported as an accepted nonzero tidal prediction.
The failure motivated one new, bounded calculation rather than another
unexplained refinement campaign.

In a homogeneous ambient region with no incoming density disturbance,
let v2=C(tau)/r^3. Since D2 v2=0 and psi0=b there, the linearized equation
has the exact particular solution

\[
\eta_2(r,\tau)=-i b r^{-3}\int C(s)ds,\qquad
2\operatorname{Re}(b\eta_2)=0.
\]

It is a harmonic phase tail, with no first-order density change. For
q=r eta this gives

\[
q\propto r^{-2},\qquad q_r=-\frac{2q}{r}.
\]

Forcing q=0 at a finite outer radius suppresses this driven phase and can
create spurious density near the boundary. It can contaminate a small
exterior tidal signal even while the core response remains insensitive.
The statement is exact for this ideal homogeneous harmonic exterior;
arbitrary radiative or Jeans disturbances need a more complete boundary
response.

The separate `nsc_atlas_hierarchy_tidal_tail.py` implements the leading
harmonic continuation as a real change to the last kinetic-Hamiltonian
diagonal. An exact rank-one correction of the existing Cayley map performs
that step; the operator remains self-adjoint and the step unitary. Dense
matrix comparison and the discrete harmonic residual check it independently.
The frozen Dirichlet experiment is preserved as a different boundary model.

The [tidal-tail record](../results/summaries/nsc-atlas-hierarchy-tidal-tail-v1.json)
contains three targeted comparisons under a separately declared 800 CPU-second
budget, admitted after the boundary diagnosis; they used about 615 CPU seconds.
The outer-domain change drops from 2.40e-14 to 7.26e-17, a factor of about 331.
The final core-density response changes by less than 1.8e-8, so this correction
does not explain away the resolved core result.

The exterior amplitude is still not resolved. Its coarse-to-medium change
is 8.80e-15; twice that change is larger than the medium coefficient 9.34e-15.
The record therefore explicitly leaves `outward_probe_resolved` false.
This targeted continuation establishes the necessary tidal-phase matching
and removes the main domain sensitivity, but does not deliver a trustworthy
nonzero external amplitude or a clock/lens prediction. Further computation
is not automatically admitted by this diagnosis.

## Verification and remaining scope

The exact inputs are the retained evolving source, two compact radial
profiles and their phase partners, their seed ages, and a declared outer
boundary model. Mathematical identities establish the real Gram bound,
orientation degeneracy, symplectic transport and harmonic tail. Numerical
comparisons evaluate those identities and responses with finite error
estimates. The core response is useful within that scope; the unresolved
exterior amplitude and universal stability remain open.

A Grok review checked the linearization, multipole kernel, real Gram argument
and centrifugal solver. Interpretation questions were resolved by writing
the actual observable: the joint base/tangent potential flow rotates the
base field, and an integrated P2 moment differs from an angular coefficient.
Nonlinear differentiation and direct angular quadrature independently test
those definitions.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_perturbations.py --output results/summaries/nsc-atlas-hierarchy-perturbations-v1.json --figure
python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_tidal_tail.py --output results/summaries/nsc-atlas-hierarchy-tidal-tail-v1.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_hierarchy_perturbations.py tests/test_nsc_atlas_hierarchy_tidal_tail.py
```
