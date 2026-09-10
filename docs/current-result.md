# What the calculations show

The project starts from one simple claim: a wave, a particle, an unseen
boundary response, and a changing geometry can be different measurements of
one inherited spectrum. This page shows the equations and numbers behind that
claim.

| Plain-language idea | Calculated connection |
|---|---|
| A particle is a stable note of a field | One scalar link gives the exact massive Dirac dispersion and propagator self-energy |
| Another room can be unseen but still act here | Direct propagation and Schur elimination give the same curved boundary response |
| A boundary stores and transfers energy | The real and imaginary parts of one coherence give link energy and occupation transfer |
| Geometry can turn into particles | Work from a changing radius equals the energy of the created Dirac pairs |
| An inherited state can push geometry apart | Both computed radial null contractions are negative at the smooth neck |
| Empty uniform space defines the zero | The recursive $a_0$ projector fixes $V_{\rm full}=0$ and preserves all gradient terms |
| Charged radius and spectral coefficients meet at one scale | The first cutoff-resolved checked branch gives $\Omega=3.9730743688$ and $\zeta=15.7853199397$ |
| The direct MMP-LLL substitution is decided | Its recursive power is fixed, but both neck null components have the wrong sign |

The construction connects the locally resolved field, the response across its
boundary, and the stress acting on its geometry through the same operator.
General research status is stated in the [README](../README.md#research-status).

The [September 10 development snapshot](development-update-2026-09-10.md),
through laboratory commit `7ee7447`, now contains fifteen post-preprint records.
Its latest gates bind the charged field, parity, AP state and normalization to
the imported MMP throat; implement canonical CTP as the causal owner; and
apply the recursive zero-tadpole law that fixes the homogeneous unlinked
$V_{\rm full}=0$. The charged radius equation now supplies a fixed-$q$ scale
candidate. Applying its simplest inherited transparent LLL state to the actual
black-universe neck gives positive rather than required negative null sources,
so the missing owner is the full charged angular/compact CTP tensor on that
domain rather than another scale scan.

The subsequent [charged CTP neck calculation](nsc-charged-ctp-neck-source.md)
returns the retained child-frame tensor
$(\rho,T_{01},p_\parallel,p_\perp)=(0.0468138,0.00122387,-0.106079,0.0302511)$
and negative radial null components $-0.0568175$ and $-0.0617129$ at the
locked scale. This passes the retained sign test, while the hard full gate
remains open for the first positive compact levels' nonlocal CTP covariance
and four ADM variations; no cosmological projection is accepted yet.

## 1. Shared mass and visible response

In the invariant scalar-link sector,

$$
[H_8,\Pi_-]=0,\qquad
E^2=\lvert\mathbf p\rvert^2+\Phi^2,\qquad
\Sigma_p(E)=\Phi^2(E+\boldsymbol\alpha\cdot\mathbf p)^{-1}.
$$

The same $\Phi$ determines the mass shell and visible self-energy.
[The Dirac reduction](nsc-observable-bridge.md) specifies the projector, spin
frame and charge conventions. The complementary sector is spectrally
identical. The evaluated massless throat kernel has vector and axial
Clifford content; the scalar coupling and sector selection are action/domain
requirements.

## 2. Evaluated curved boundary maps

Parent and child maps are obtained from the smooth Dirac geometry with
oriented normal data. Direct inversion, Schur elimination and transfer agree
below $2\times10^{-15}$ relative error in the recorded finite system.
Independent radial integration differs by 0.00947%–0.03645% at the stated
resolution. The [full spinor calculation](nsc-chiral-boundary.md) retains both
angular sectors and the common spin frame.

At common dimensional energy, the recursive response is

$$
\Gamma_p(x)=K_p(x)-\Omega^{-1}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.
$$

The explicit $\Omega^{-1}$ is the first-order normalization factor.
[Finite recursion](nsc-regulated-recursion.md) and
[endpoint criteria](nsc-tail-limit.md) supply the recorded controls.

## 3. The state completes the energy account

For the finite quadratic system,

$$
z_B=\mathrm{Tr}_p(B C_{cp}),\qquad
E_{\rm link}=2\Re z_B,\qquad \dot N_p=2\Im z_B,
$$

$$
\dot E_{\rm tot}=\mathrm{Tr}(C\dot H).
$$

[The common-source derivation](nsc-common-source-derivation.md) includes
energy stored in the link. [The boundary-state implementation](nsc-boundary-state.md)
retains occupied and empty kernels and mixed initial correlations.
Its finite covariance reconstruction agrees with full evolution to
$1.1\times10^{-15}$.

A prescribed changing radius creates pairs with excitation energy equal to
the supplied work. In the representative $\kappa=1$ channel,
$\epsilon=0.02$, $\tau=0.5$ gives approximately
$1.40181\times10^{-5}$ pairs and energy $3.46640\times10^{-5}$
in throat units. [The vacuum-work record](nsc-vacuum-work.md) and
[normalized influence functional](nsc-influence.md) connect this response
to its geometric fluctuation weight.

## 4. Parent data produce a definite child stress

Affine-horizon data and the actual complex exterior reflection determine
the canonical massless state on the expanding child. In the child
orthonormal frame, with $\hbar=c=L_{\rm throat}=1$, $\mu=1$:

| Quantity | Recorded neck value |
|---|---:|
| Density $\rho$ | $-0.0033068962$ |
| Longitudinal pressure $p_\parallel$ | $-0.0494945966$ |
| Angular pressure $p_\perp$ | $-0.0429957713$ |
| Radial null $+$ | $-0.0528075900$ |
| Radial null $-$ | $-0.0527953956$ |
| Outward parent Killing power | $0.0001422206795$ |

Stress has dimensions $L^{-4}$, and power $L^{-2}$. Angular refinement
from 16 to 32 after the recorded tail correction changes tensor components
by less than $4\times10^{-6}$ across the sampled radii. This measures
numerical sensitivity under that asymptotic-tail approximation.
[Source and state](nsc-unruh-state.md).

The positive-Einstein static neck requires

$$
(\rho,p_\parallel,p_\perp,T_{\hat T\hat z})
=(2a,2a(1-3\pi),2a(1-3\pi),0),\quad a>0.
$$

The computed null signs are negative; density, pressure anisotropy and flux
leave explicit componentwise residuals. The [ADM frame map](nsc-adm-neck-source-map.md)
places those quantities in the metric equations. Under the stated asymptotic
Einstein normalization, the outgoing channel contributes
$dM_B/du=-P_0$, selecting a backreaction calculation.
[Flux and mass balance](nsc-parent-backreaction-gate.md).

## 5. One normalization and one metric variation

The [canonical–spectral bridge](nsc-canonical-spectral-bridge.md) partitions
the compact determinant into its physical tower, cutoff conversion and warp
contribution. The [Gaussian measure](nsc-gaussian-cutoff-measure.md) gives an
exact finite representation of the adopted proper-time action.
[Compact anomaly matching](nsc-compact-anomaly-bridge.md) fixes the domain and
the scope of the imported anomaly result.

The [ADM source equations](nsc-adm-source-constraints.md) retain lapse, shift,
radial metric and sphere radius through variation. Their Ward identities
bind stress, energy, flux and geometric work. The
[local warp source](nsc-warp-local-neck-source.md) evaluates an identified
summand on the actual neck; its nonlocal remainder has no small-error bound.

## 6. The causal and vacuum source now have explicit owners

Canonical Lorentzian evolution supplies the causal state-dependent response:

$$
\Gamma_{\rm one}^{\rm CTP}
=S_{\rm induced}[+]-S_{\rm induced}[-]
-i\log\det(I-C_0+C_0U_-^\dagger U_+).
$$

The [causal common functional](nsc-causal-common-functional.md) returns the
four metric forces, stress projections, boundary power, retarded response,
noise, and discrete Ward residuals from one covariance. Equal histories,
unitarity, covariance trace, and energy/work balance agree at approximately
$10^{-16}$ in its finite control.

The homogeneous gravitational zero is fixed separately by the same-action
relational law

$$
\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one},
\qquad
(V,A,C)\mapsto(0,A,C).
$$

This gives $V_{\rm full}=\lambda_4=\Xi=0$ while preserving Einstein, gauge,
curvature, Casimir, and link response. The
[normalization derivation](nsc-relational-vacuum-normalization.md) proves
idempotence, additivity, CTP normalization, and common room-scale covariance.

The physical geometry is the common stationary point:

$$
\left.\frac{\delta\Gamma_{\rm one}^{\rm CTP}}
{\delta g_\Delta^{\mu\nu}}\right|_{g_\Delta=0}=0.
$$

The same solution determines the following physical links:

| Target | Quantity to extract from the common solution |
|---|---|
| Constants and resolution | Characteristic cones, gravitational coupling, physical thermal state |
| Finite continuation | Evolved geometry, constraints, stress and boundary flux |
| Dark response | Background pressure and density, anisotropic stress, growth and lensing |
| Particles and measurement | Poles, residues, interactions, detector statistics |
| Mass and antimatter | Chiral interaction, sector selection, gauge and conjugate representations |
| Recursive infinity | State marginals, clock transfer and limiting observables |

## Inspect the record

The [preprint manifest](../results/manifest.json) contains 100 records;
[the later snapshot](../results/development-snapshot.json) contains fifteen.
[The development index](development-update-2026-09-10.md) links all fifteen
new derivation notes. The canonical paper retains its technical appendices,
attribution and historical scale corrections.

The default demonstration displays stored evidence. Use
[the reproduction guide](reproducing.md) for an explicitly requested replay.
