# Use the full spectral source at the child curvature

The stored development geometry reaches curvatures for which the local
heat expansion is not controlled. Applying an existing four-sphere Dirac
spectrum to the already computed compact weight gives a usable full
Euclidean endpoint response. At the child curvature its radius derivative
has the opposite sign to the local approximation. A stationary radius of
that approximation is not stationary in the full determinant.

This decides the next source formulation: use the original spectral
functional at these curvatures. The local curvature EFT remains applicable
within a demonstrated perturbative regime. Neither calculation yet supplies
the complete quantum functional or the global parent-to-child state.

## Reused evidence and the missing connection

The [clock and child geometry](../results/nsc-5-clock-horizon.json) already
establish the asymptotic directional rate H=sqrt(3 pi) in the benchmark's
length units. The [horizon source](nsc-horizon-source.md) supplies the
recorded horizon. The [finite-term owner](nsc-finite-terms.md) supplies the
curvature contractions, and the [compact matching](nsc-compact-matching.md)
supplies the full weight h and its local coefficients. Their generators
are not rerun.

For the endpoint, import the sphere spectrum from
[Camporesi and Higuchi, equations 3.51–3.54](https://arxiv.org/html/gr-qc/9505009v1).
Their anti-Hermitian eigenvalues become real after the conventional factor
of i. On a four-sphere of radius a, the self-adjoint D4 has eigenvalues
plus/minus m/a, m=2,3,..., with 2(m^3-m)/3 states per sign. Squaring doubles
that multiplicity. The sphere volume is 8 pi^2 a^4/3. These are established
results; the new calculation composes them with the project's h.

The decision is whether the local action can be used for the intended
source solve, or whether the computed spectral remainder changes it. No
new geometry solver, collapse calculation, cutoff fit or vacuum coefficient
is introduced for this comparison.

## Curvature and the actual eliminated scales

For ds^2=A dt^2-d rho^2/A-r^2 dOmega^2, substitute the existing curvature
owner into the quasi-global chart:

\[
a_1=A''/2,\quad b_1=A'r'/(2r),\quad
c_1=A r''/r+b_1,\quad d_1=(1-A(r')^2)/r^2.
\]

Its scalar contractions extend through A=0. No static observer frame is
assigned to the trapped interior. Use the largest absolute sectional
entry as one local curvature scale; it is not a general proof controlling
all covariant derivatives, gauge fields, states or nonlocal effects.

The stored parameters are Lambda=2, ell=2, warp path 1 and matching nu=1.
The first nonzero compact eigenvalue at zero four-momentum is approximately
2.496676. They are development inputs, not identified physical constants.

| Location | Largest sectional magnitude | Divided by Lambda^2 | Divided by nu^2 | Divided by first compact eigenvalue |
|---|---:|---:|---:|---:|
| Recorded horizon | 0.216796 | 0.0542 | 0.2168 | 0.0868 |
| Areal-radius neck | 4.712389 | 1.1781 | 4.7124 | 1.8875 |
| Asymptotic child | 9.424778 | 2.3562 | 9.4248 | 3.7749 |

Thus this parameter choice supplies no small curvature hierarchy at the
neck or child limit. The smaller horizon ratios alone do not certify a
complete local approximation there. Under a physical dilation of the
profile by L, curvature scales as L^-2: these comparisons depend on
Lambda L, nu L and the actual compact scale, not a change of units. They
do not determine L, zeta, Omega or the full cW/A ratio.

## Full endpoint response versus the local action

The existing two-copy normalization gives

\[
\Gamma_5(a)=\frac23\sum_{m=2}^{\infty}(m^3-m)
 h_\Lambda(m^2/a^2),
\qquad
\partial_{\ln a}\Gamma_5
=-\frac43\sum_{m=2}^{\infty}(m^3-m)\frac{m^2}{a^2}
 h_\Lambda'(m^2/a^2).
\]

There is no additional copy or spin multiplicity. The mode sum already
includes the volume dependence. The dimensionless derivative holds cutoff,
compact size and warp fixed: it is a partial metric equation, not a
simultaneous scale solution.

The comparison retains the canonical light field exactly and approximates
only the matched complement through a4:

\[
\Gamma_{\rm local}(a)=\Gamma_{\rm light,\nu}(a)
 +\frac{8\pi^2a^4}{3}\left(V_H-\frac{12 A_H}{a^2}\right)
 +64\pi^2 C_{E,H}.
\]

On S4 the Weyl tensor and box-R vanish; the integrated Euler term is
constant. The free Dirac complement has no independent R-squared term.
This does not fix finite coefficients of the completed theory. Both
expressions use the original physical metric, before the optional EFT
field redefinition. The induced Einstein contribution is included once.

| Sphere radius a | Full log-radius derivative | Local log-radius derivative |
|---|---:|---:|
| 0.325735, matching the stored child H | 0.001239600 | -0.164273768 |
| 0.43742254, local stationary control | 0.084326048 | approximately 0 |
| 1 | 14.8033344 | 14.8393135 |
| 2 | 284.848997 | 284.854519 |
| 4 | 4754.593055 | 4754.594317 |

The local stationary control is found without fitting an observation.
It is a root of the stated approximation, used to test that approximation.
The nonzero full derivative rejects that root as a stationary point of
the same full Dirac modulus. The near agreement at larger a supports the
expected low-curvature use of the local coefficients. The calculation
does not prove absence of every stationary branch of h, or of the complete
action with its remaining state, compensator and interactions.

For a homogeneous Euclidean variation, define
Pi_a=(partial_log_a Gamma)/(4 Vol(S4)). At the child-radius control this is
approximately 0.001045912 in the declared units. It is a scalar projection
of the full Dirac metric variation, which already includes its induced
geometric terms. It must not be added as a second matter stress alongside
those same terms. Lorentzian stress identification requires the causal
state and continuation below.

## Error control and physical domain

The generalized compact eigenvalue derivative is compared with an
independent finite difference of the total action using two step sizes.
Compact mode counts 18,26,34, overlap quadrature 112,144, and eight extra
sphere levels are varied separately. At the local stationary control,
the last compact change in the full derivative is about 9.7e-10; its value
is 0.084326. The decision is well separated from the measured errors.

For s=Lambda a and last retained sphere index M beyond the decreasing
envelopes, the established warp min-max comparison gives

\[
|\Gamma_{m>M}|\le\frac{C s^4}{3}e^{-M^2/s^2},\qquad
|\partial_{\ln a}\Gamma_{m>M}|
\le\frac{2e^{\alpha}C}{3}s^2(M^2+s^2)e^{-M^2/s^2},
\quad C=1+\frac{\ell\Lambda}{\sqrt\pi}.
\]

Here alpha is the maximum absolute compact warp. The derivative bound
uses Cauchy–Schwarz on the existing generalized quadratic form:
|partial_y epsilon| <= exp(alpha) sqrt(epsilon/y), with epsilon>=y.
Integral comparison of the decreasing m exp(-m^2/s^2) and
m^3 exp(-m^2/s^2) envelopes then gives the displayed bounds. These
continuum sphere-tail bounds are evaluated in floating point; they do
not bound errors of the retained Galerkin eigenvalues. The latter have
measured convergence, not a proved remainder bound. A direct unwarped
double sum checks the new sphere/compact composition independently.

The round S4 operator is untwisted, with no gauge flux, link mass or
spacetime zero mode. A regular Euclidean sphere can serve as an invariant
de Sitter vacuum construction only with its Lorentzian continuation and
state prescription. It is not the global continuation of the black-universe
chart: that interior has Kantowski–Sachs spatial topology, and its finite
time shear, magnetic flux, inherited correlations and throat domain do not
automatically extend to the round sphere. Local decay of magnetic density
does not prove a global bundle or state identification.

The existing [Allen–Lütken spinor two-point construction](https://pure.mpg.de/view/item_153497)
is a reference for that state matching. A bounded Grok investigation
independently confirmed the Camporesi–Higuchi multiplicity and identified
this state reference. The state identification is not counted as a
completed physical result here.

The full positive spectral modulus is defined here. Its determinant phase,
compensator, finite normalization terms, quantum metric measure and global
in-in state have not been supplied by this calculation. In particular it
does not derive a cosmological energy-transfer rate or a physical dark
fraction. Working directly in the original variables avoids introducing
the optional field-redefinition Jacobian into this Dirac-only comparison;
it does not solve the remaining quantum measure problem.

Evidence: [immutable record](../results/development/spectral-endpoint.json),
reproducer scripts/check_nsc_spectral_endpoint.py --check. Comparison covers
all scientific fields, formulas, scope statements and source/input hashes.
Previous generators remain unused.
