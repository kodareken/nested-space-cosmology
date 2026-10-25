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
two nonlinear radial derivative checks fit a4500 CPU-second budget, with
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
(A_seed*x)^l exp[-(A_seed*x)^2/(2w^2)], with w=2 and4 and a taper from
physical scaled radius15 to20. Include each shape and its imaginary phase
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
quadrupole orientations. Angular orthogonality gives the direct-sum Gram
matrix I5 tensor G. Its largest eigenvalue is the same one. Thus the
four-column calculation bounds every combination in the corresponding
twenty-dimensional real quadrupole span, not just one chosen orientation.
It does not cover every radial profile, higher angular multipoles or
nonlinear perturbations of finite size.

The source is sampled each unit of slow time. Numerical maxima refer to
those saved events, not a rigorous continuous-time enclosure. A second
quadrupole family is initialized on the evolved state at tau480 and followed
to600. This tests a late breathing phase independently of the initial seeds.

## Physical deformation and an independent nonlinear check

The normalized second-moment quadrupole uses an additional x^2 weight.
For the axisymmetric component its first-order coefficient is

\[
\frac{\int R^2P_2\,\delta\rho\,d^3R}
{\int R^2\rho_{\rm excess}\,d^3R}
=\frac{\varepsilon N_c}{\sqrt5}
\frac{\sum_jc_j\int x^4d_jdx}{\int x^4d_0dx}.
\]

The angular average of P2*Y20 is1/sqrt(5), which fixes the factor here.
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

## Current calculation

The equations and analytic tests are implemented. The finite response
campaign is active; the first coarse result distinguishes radial sensitivity
from modest quadrupole response. Final values and numerical ranges await
the other resolutions and the nonlinear derivative check. A source-hash
mismatch, failed symplectic identity, unresolved derivative or significant
domain sensitivity would stop that interpretation.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_perturbations.py --output results/summaries/nsc-atlas-hierarchy-perturbations-v1.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_hierarchy_perturbations.py
```
