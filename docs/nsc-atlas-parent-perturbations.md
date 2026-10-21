# Gravitational structure from the inherited parent scalar

The [shared scalar atlas](nsc-atlas-shared-scalar-slice.md) now has a
constraint-compatible perturbation supported inside its actual parent.
Keeping the action and parent background fixed, the coupled linear
Einstein-scalar evolution raises its central comoving density contrast by
15.188159% through one e-fold. A fixed-metric field control instead gives
a 1.606587% decrease. Longer and shorter modes have different responses,
and an independently variable initial scalar velocity changes the result.

This establishes conditional gravitational structure growth, beyond the
previous homogeneous dilution calculation. It does not establish nonlinear
halos, an inherited primordial spectrum, observed abundance or unique
evidence for nesting. The complete perturbed wall future is not assumed.

## Decision and borrowed mathematics

The alternatives were full coupled wall evolution, an averaged dust/Jeans
estimate, and the unaveraged scalar/metric response in the protected parent.
The third gives a direct structure test while respecting the unresolved
wall branch. Its perturbation variables are the field displacement, its
velocity and the metric potential. They remove the ambiguity of treating
background dilution as clustering, without adding a physical coupling.

The general curved-background Einstein-scalar perturbation equations are
established mathematics: see [Hwang and Noh, Sections 2 and 4](https://arxiv.org/html/astro-ph/0103244).
Their conventions are translated below and the resulting constraints,
spatial Einstein equation and energy conservation are checked symbolically.
The contribution here is their connection to the actual shared-parent
data, a perturbation that preserves the initial membrane joins, its finite
transfer and the scale condition on this inherited scalar.

[Hwang and Noh, equations 23--27](https://arxiv.org/html/0902.4738) supply
the averaged massive-scalar Jeans comparison. Their leading oscillation
approximation and flat-background assumptions are not imported as the
parent evolution. The main calculation retains curvature and every scalar
oscillation. [Noh, Hwang and Barrow](https://arxiv.org/html/2003.14108v2)
also discuss the special lowest harmonics and nonlinear integrability in
closed cosmologies. Here n=0 and n=1 are excluded by an explicit compensated
preparation, not dismissed as universally unphysical modes. No nonlinear
integrability or finite-amplitude error theorem is claimed.

The admitted calculation uses four diagnostic modes, a compact packet at
three spectral cutoffs, two time formulations and a fixed-metric control.
It has a 120 CPU-second budget and stores only a compact record, producer
hashes and coefficient hashes. Failure of constraint propagation or the
declared comparisons stops the interpretation. There is no wall-grid scan.

## Same action, regular perturbation equations

In the parent, q is fixed and the scalar potential is exactly

\[
U=\rho_v+\frac{\omega^2\chi^2}{2},\qquad
\rho_v=9/32,\quad \omega=20,\quad \kappa=1.
\]

Use the existing background with a0=6.626077076100, H0=.442688940038,
chi0=.0433012701892 and Pi0=0. Its excitation E0=3/8 came from placing the
common initial field at the grandchild minimum. In general, the two equal
charge steps give E0=2(1-epsilon)e^2; this preparation's E0/rho_v=4/3.
The packet introduced below is additional specified initial data, not a
state-selection result obtained from that homogeneous inheritance rule.

In zero-shear/Newtonian gauge write

\[
ds^2=-(1+2\Phi)dt^2+a^2(1-2\Phi)\gamma_{ij}dx^idx^j,
\qquad {}^{(3)}K_\gamma=1.
\]

There is no linear anisotropic stress for this canonical scalar on a
homogeneous background, so the two metric potentials agree. For a scalar
harmonic with Laplacian eigenvalue -lambda, lambda=n(n+2), n>=2, set
X=delta chi and Y=dX/dt. Then

\[
\dot X=Y,\qquad
\dot Y=-3HY-(\omega^2+\lambda/a^2)X
             +4\Pi\dot\Phi-2\omega^2\chi\Phi,
\]

\[
\dot\Phi=-H\Phi+\frac{\kappa}{2}\Pi X.
\]

The energy constraint is

\[
C=\left(\frac{\lambda-3}{a^2}-\frac{\kappa\Pi^2}{2}\right)\Phi
 +\frac{\kappa}{2}\left[\Pi Y+(\omega^2\chi+3H\Pi)X\right]=0.
\]

Direct differentiation using these equations and the existing background
gives the exact identity

\[
\dot C=-3HC.
\]

The spatial Einstein equation follows independently:

\[
\ddot\Phi+4H\dot\Phi+(2\dot H+3H^2-a^{-2})\Phi
=\frac{\kappa}{2}\delta p,
\qquad \delta p=\Pi Y-\Pi^2\Phi-\omega^2\chi X.
\]

No equation divides by Pi or by the time-varying coefficient multiplying
Phi in C. The integration passes eighteen scalar-velocity sign changes.
Solving a constraint by dividing by Pi would fail at precisely these
ordinary turning points; no physical singularity is inferred from them.

## A density and metric relation without a gauge singularity

The directly evaluated perturbed density is
delta rho=Pi Y-Pi^2 Phi+omega^2 chi X. Its comoving combination remains
well-defined even at Pi=0:

\[
\delta\rho_c=\delta\rho+3H\Pi X,
\qquad \Delta=\frac{\delta\rho_c}{E}
=-\frac{2(\lambda-3)\Phi}{\kappa a^2E}.
\]

Under an infinitesimal time shift T, delta rho gains 3H Pi^2 T while
X loses Pi T, so delta rho_c is unchanged. E is the strictly positive
background excitation, not the total vacuum-plus-scalar density.

The same field thus supplies density and metric transfer; the local
lensing-potential combination is Phi+Psi=2Phi. This is a linear local
relation, not a ray-traced lensing survey prediction. The observable
geometry and normalization remain those of the protected parent.

## A packet that preserves the actual atlas joins

Let theta be angle from the regular parent pole and choose theta_b=1.4.
For x=cos(theta) and c=cos(theta_b), define weight
w(x)=(x-c)^8 sqrt(1-x^2), c<=x<=1. Let mu, v and m3 be its normalized
mean, variance and third central moment. The exact initial potential is

\[
\Phi_0(\theta)=A(x-c)^8
\left[(x-\mu)^2-\frac{m_3}{v}(x-\mu)-v\right],\quad \theta<\theta_b,
\]

and zero elsewhere. Choose A so Delta0 at the pole equals one; a physical
linear forecast multiplies the whole response by a small epsilon, here
1e-5. The exact potential is C7 at its support boundary and regular at the
pole. Its n=0 and n=1 projections vanish by the two defining moments.
These are preparation choices, not gauge claims or membrane eigenmodes.

At Pi0=0 choose Y0=0 and set

\[
X_0=\frac{2(\nabla_\gamma^2+3)\Phi_0}
          {\kappa a_0^2\omega^2\chi_0},\qquad
\dot\Phi_0=-H_0\Phi_0.
\]

The field profile is C5. The linear spatial metric changes by
delta h=-2Phi0 h, while the mixed slice extrinsic curvature remains
unchanged. Its curvature variation is
delta R3=4(Delta_gamma+3)Phi0/a0^2, exactly 2 kappa delta rho0.
The momentum constraint holds because the background scalar momentum
vanishes and its initial field is spatially uniform. All perturbation data
vanish on an open collar of both walls, so their initial scalar, charge,
surface-stress and Israel data remain unchanged. This is an explicit
linear perturbation of the full compact initial atlas.

The existing total conformal bound is eta_infinity<=.515418396612,
and the parent wall starts at theta_w=2.850764359770. The margin

\[
\theta_w-\theta_b-2\eta_{\rm bound}=.419927566546>0
\]

keeps the packet's entire causal cone inside the protected parent domain.
The smooth-parent perturbation can therefore be calculated without a
reflecting condition at the material wall. This statement concerns the
linear response on the existing background domain; it does not construct
the entire nonlinear membrane spacetime.

We use radial S3 harmonics Q_n=sin((n+1)theta)/[(n+1)sin(theta)] as a
computational extension, with their polynomial continuation at the poles.
The infinite expansion reconstructs the compact profile. A finite expansion
has small tails; these are measured rather than declared identically zero.
At cutoffs 64,128,256, maximum sampled initial density reconstruction errors
are 2.79e-7,3.13e-9,1.55e-11. At the endpoint, tails beyond the causal support
are 8.94e-9,1.23e-10,1.71e-12 respectively.

## Calculated growth and the remaining initial-state freedom

At one e-fold, proper time2.88884551314, the central normalized density is

\[
\Delta_{\rm center}=1.15188159122.
\]

For the example amplitude epsilon=1e-5, this is the linear forecast
1.15188159122e-5. The central contrast oscillates during the interval; the
reported increase is an endpoint change, not a monotonic growth theorem.

Removing the metric response while retaining the same background and
initial scalar displacement yields .983934130329. This test-field control
deliberately violates the perturbed Einstein constraints; it isolates the
role of gravitational response and is not a second coupled atlas solution.
The difference in the two central response coefficients is .167947460890.

For each diagnostic harmonic we retain both independent initial columns:
density seed (Delta0=1,Y0=0) and phase seed (Delta0=Phi0=X0=0,
Y0=omega chi0). Final linear density is their linear combination.

| n | Initial physical wavenumber | Density-seed response | Phase-seed response |
| --- | ---: | ---: | ---: |
| 2 | .426863 | 1.15091875 | .00702776 |
| 8 | 1.349859 | 1.15226198 | .10801274 |
| 32 | 4.978035 | .92380409 | 1.36724676 |
| 128 | 19.467948 | -.98302417 | -.98772566 |

The shorter density seeds are suppressed or change sign over the interval,
showing oscillatory response rather than uniform dust growth. A sign change
alone is not a theorem of stability at every later time. The phase column
also proves that equal initial density and metric do not specify a unique
later density. The scalar velocity or its statistical preparation is needed.

![Inherited parent perturbation growth and scale dependence](figures/nsc-atlas-parent-perturbations.png)

## The vacuum and structure scales must work together

The leading averaged, short-curvature-scale Jeans estimate is

\[
k_{J,\rm physical}^4=2\kappa E\omega^2,\qquad
\ell_\Lambda=\sqrt{3/\Lambda},\qquad \Lambda=\kappa\rho_v.
\]

For this same inherited preparation, its initial values give
kJ=4.161791450288 and lambdaJ/ellLambda=.462258805666. The background
H/omega=.022134447 and curvature/(a^2 kJ^2)=.001315002 quantify two
small comparison parameters; they are not a certified error bound.
Neither this estimate nor its dust approximation was used to evolve the
field. The exact linear calculation above supplies the scale dependence.

Because E0/rho_v=4/3 here, the estimate implies

\[
\frac{\omega}{\sqrt{\Lambda}}
=\frac{\pi^2\sqrt{2/3}}{(\lambda_J/\ell_\Lambda)^2}.
\]

This is a concrete additional requirement on a physical identification:
the present benchmark's initial pressure scale is a substantial fraction
of its vacuum radius. It cannot simply be described as cold on every
subcosmological scale. A much smaller cutoff requires a correspondingly
larger omega/sqrt(Lambda), together with a compatible history. The current
inheritance rule does not select that ratio or an observed power spectrum.

The local equations are established scalar-gravity physics. An unnested
model with the same local initial data would have the same response in this
protected domain. Nesting contributes the composed geometry and allowed
inheritance/preparation relations; this growth alone is not unique evidence
of an unseen surrounding nest.

## Verification and next connection

Eight symbolic identities include constraint propagation, the independent
spatial Einstein equation, energy conservation and gauge invariance. Two
clock formulations differ by at most3.01e-9 in the diagnostic transfer
coefficients. The final packet changes by2.94e-9 from cutoff128 to256;
its independent-clock maximum sampled difference is8.07e-11. The retained
comparison tolerance is1e-7 in normalized linear response. These are finite
numerical consistency checks, not a uniform spectral enclosure or a
finite-amplitude nonlinear remainder.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_parent_perturbations.py
python scripts/run_python.py scripts/plot_nsc_atlas_parent_perturbations.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_parent_perturbations.py
```

The [record](../results/summaries/nsc-atlas-parent-perturbations-v1.json)
contains the initial state, transfer basis, compact profile prescription,
cutoff comparisons and hashes. No wave history or large coefficient array
is committed. The background rates are shared with the preceding parent
calculation rather than maintained as a second evolution implementation.

Next: connect the prepared perturbation state or covariance to the parent
formation/inheritance mechanism, and use the linked density/metric transfer
for observer and lensing calculations. Physical scale selection is a
required test before calling this benchmark a model of measured dark matter.
