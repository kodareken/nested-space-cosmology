# A scalar layer, a charged boundary and a prepared expanding child

The [conversion calculation](nsc-atlas-flux-conversion.md) supplied a
potential and an energy account, but assumed the state after a rapid
transfer. Here the scalar profile, gravitational geometry and charged
boundary are solved together. Regular Euclidean cap data select a
candidate preparation; the Lorentzian continuation gives an expanding
wall and open child cosmology. This is a different prepared geometry
from P002's two-wall Kottler/dust initial data.

The result combines exact matching/reduction identities with a converged
nonlinear boundary-value calculation. It is not a validated interval
existence theorem, a proof of uniqueness, or a tunnelling probability.

![Scalar profile, coupled geometry and continued wall](figures/nsc-atlas-scalar-preparation.png)

## The necessary spatial connection

Suppose a homogeneous child phi_c(t) is joined directly to a constant
parent phi_p across a timelike wall without a scalar surface coupling.
Scalar continuity gives phi_c(t_W(tau))=phi_p. In the FLRW child,
timelikeness gives dt_W/dtau>=1. Differentiating continuity therefore
forces phi_dot_c=0 along the whole matched interval. The scalar field
equation then requires U_phi(phi_p,q_c)=0. But at the parent minimum
U_phi(phi_p,q_p)=0 and

$$
U_\phi(\phi_p,q_c)-U_\phi(\phi_p,q_p)=-\mu e\ne0.
$$

Thus the direct homogeneous construction cannot satisfy the stated
equations. A spatial scalar layer, a changed parent field or an additional
surface interaction is necessary. We keep the action and add the layer;
we do not insert a scalar impulse at the membrane.

## Inputs and borrowed geometry

Retain

$$
U(\phi,q)=\rho_b+\frac{m^2\phi^2}{2}+\frac{(q+\mu\phi)^2}{2},
\quad \omega^2=m^2+\mu^2,\quad \epsilon=m^2/\omega^2.
$$

The numerical inputs are epsilon=1/4, omega=20, q_p=1.5, q_c=1,
rho_b=0 and constant membrane tension sigma=.2. The normalized arithmetic
coupling is kappa=1, as in the conversion owner; no SI ruler is selected.
The two relaxed densities are rho_i=.125 and rho_o=.28125.
The radius and both scalar pole values are outputs, not supplied wall data.

The scalar/four-form action is the attributed
[Kaloper--Sorbo construction](https://arxiv.org/html/0811.1989), with our
explicit additional positive scalar mass. The field/gravity equations
and continuation reuse [Coleman--De Luccia](https://doi.org/10.1103/PhysRevD.21.3305)
and [Dong--Harlow, Sections 2.1--2.2](https://arxiv.org/html/1109.0011).
The charged-wall convention remains the one checked against
[Hirano, Section 2.3](https://arxiv.org/html/1804.09985).
No supersymmetric or field-dependent membrane tension is imported.

Closely related combined scalar/membrane/gravity preparations are already
published in [Bandos, Blanco-Pillado, Sousa and Urkiola, Section5](https://arxiv.org/html/2306.09412).
Their field-dependent membrane tension also sources a scalar derivative
jump. Our constant-tension quadratic example retains continuous normal
derivatives. The general combination is prior art; the result here is
the stated model's checked profile, boundary reduction and child state.

**Finite question and budget.** Construct one regular two-cap solution
with the declared action, then determine its field transfer and child
excitation. Compare adaptive collocation tolerances 1e-7 and 1e-9,
with at most 12000 nodes and 3000 vector RHS calls per solve. Use separate
logarithmic transfer equations and direct action quadrature as checks.
The evidence creation has a 120 CPU-second budget. There is no parameter
sweep, P001 grid campaign or independently chosen incoming waveform.
Failed constraints stop this preparation candidate; success closes this
branch rather than inviting finer certification of its tiny tail.

## Coupled cap equations and complete wall conditions

Each cap has its own coordinate xi_j from its regular pole to the wall:

$$
ds_E^2=d\xi_j^2+a_j(\xi_j)^2d\Omega_3^2,\qquad
\chi_j=\phi_j-\phi_*(q_j),\qquad
\phi_*(q)=-\mu q/\omega^2.
$$

For j=i,o the equations are

$$
\chi_j''+3\frac{a_j'}{a_j}\chi_j'=\omega^2\chi_j,\qquad
a_j''=-\frac{\kappa a_j}{3}
\left[\chi_j'^2+\rho_j+\frac{\omega^2\chi_j^2}{2}\right].
\tag{S1}
$$
$$
{\cal C}_j=a_j'^2-1-\frac{\kappa a_j^2}{3}
\left[\frac{\chi_j'^2}{2}-\rho_j-\frac{\omega^2\chi_j^2}{2}\right]=0.
\tag{S2}
$$

Direct differentiation gives C_j'=0 on S1. Regular pole data
a_j(0)=0, a_j'(0)=1, chi_j'(0)=0 therefore impose S2.
The near-pole series, with phi_0 the computed pole value, is

$$
a_j=\xi_j-\frac{\kappa U(\phi_0,q_j)}{18}\xi_j^3+O(\xi_j^5),
\qquad
\chi_j=\chi_j(0)\left[1+\frac{\omega^2\xi_j^2}{8}
+O(\xi_j^4)\right].
$$

The unknown pole-to-wall lengths are L_i and L_o. At the wall,

$$
a_i=a_o=R_*,\quad
\chi_i-\chi_o=-\Delta,\quad
\chi_i'+\chi_o'=0,\quad
a_i'+a_o'=kR_*,
\tag{S3}
$$
$$
\Delta=\phi_*(q_c)-\phi_*(q_p)=\frac{\mu e}{\omega^2},
\qquad k=\kappa\sigma/2.
$$

The signs matter: the physical outward derivative on the parent side is
minus its pole-to-wall derivative. Hence K_i=3a_i'/R_* and
K_o=-3a_o'/R_*. The ordinary branch has a_i'>0>a_o' at the wall.
Eight first-order equations, two unknown lengths, six pole conditions
and four wall conditions define the boundary-value problem.

S2 and S3 imply the full normal force:

$$
[K^a{}_b]=-k\delta^a_b,\qquad
\sigma\overline K=U(\phi_W,q_p)-U(\phi_W,q_c)
=e(q_p+\mu\phi_W-e/2).
\tag{S4}
$$

The scalar derivative is continuous in the common normal, and all
tangential derivatives vanish in the symmetric wall frame. The
tangential flux jump is zero. These are the distributional conservation
conditions of the common action, not a pressure fitted afterward.

## A physical scalar boundary reduction

At the solved geometry, define regular-cap maps
s_i=chi_i'(L_i)/chi_i(L_i) and s_o=chi_o'(L_o)/chi_o(L_o).
For a positive cap scale factor, (a^3 chi')'=a^3 omega^2 chi shows
that a nonzero regular mode keeps its sign and grows in absolute value
from the pole. It has no zero on the cap, and both coefficients are positive.
These are Dirichlet-to-Neumann coefficients of the actual scalar operator

$$
\mathcal D_j=-a_j^{-3}\partial_{\xi_j}
(a_j^3\partial_{\xi_j})+\omega^2.
$$

More explicitly, D_j acts as
-a_j^-3 d(a_j^3 d chi/d xi_j)/d xi_j + omega^2 chi.
The maps give

$$
\chi_i(L_i)=-\frac{\Delta s_o}{s_i+s_o},\qquad
\chi_o(L_o)=\frac{\Delta s_i}{s_i+s_o}.
\tag{S5}
$$

The bulk equation gives the exact boundary identity
(a^3 chi chi')'=a^3(chi'^2+omega^2 chi^2). Consequently the scalar
excitation part of the Euclidean action, evaluated at this fixed geometry, is

$$
S_{\chi,\rm bulk}
=\pi^2R_*^3\left(s_i\chi_i^2+s_o\chi_o^2\right)
=\pi^2R_*^3\Delta^2\frac{s_i s_o}{s_i+s_o}.
\tag{S6}
$$

This is an explicit physical Schur elimination. It is the regular
Euclidean saddle sector, not a retarded Lorentzian memory kernel or
an additional bulk Lambda. The gravitational and vacuum parts of the
action remain separate; variations of the geometry also vary the maps.

The finite calculation gives

| Quantity | Value |
| --- | ---: |
| R_* | 2.70634531321651 |
| L_i | 2.86731015671418 |
| L_o | 7.07049444470919 |
| phi_W | -.0543355679824020 |
| Common normal scalar derivative | -.215595888249873 |
| Wall pressure sigma K_bar | .154440177981835 |
| s_i | 19.5387048900498 |
| s_o | 20.3079350371445 |
| s_i s_o/(s_i+s_o) | 9.95794752937923 |
| Scalar excitation action | .913188612281953 |

Direct bulk quadrature differs from S6 by 2.8e-13 or less. The scalar
layer's energy per membrane area on the symmetric initial slice is
.00233741846805, computed using
R_*^-2 sum_j integral a_j^2(chi_j'^2+omega^2 chi_j^2)/2 dxi_j.
It is distinct from the bare membrane tension .2.

## Lorentzian continuation and the child state

Continue a sphere angle to pi/2+i eta in each open cap:

$$
ds_L^2=d\xi_j^2+a_j(\xi_j)^2
\left[-d\eta^2+\cosh^2\eta\,d\Omega_2^2\right].
\tag{S7}
$$

The bulk equations are still S1--S2, with a Lorentzian de Sitter
three-slice. The same surface satisfies S3--S4. Its proper time is
tau_W=R_* eta and its areal radius is

$$
R_W(\tau_W)=R_*\cosh(\tau_W/R_*).
$$

The eta=0 slice is an additional direct check of the physical connection.
It has induced metric dxi^2+a^2 dOmega_2^2 and vanishing extrinsic
curvature in time. S1--S2 give

$$
{}^{(3)}R=-\frac{4a''}{a}+\frac{2(1-a'^2)}{a^2}
=2\kappa\left(\frac{\phi'^2}{2}+U\right).
$$

The wall derivative jump [a']=-kR_* contributes
-4[a']/R_*=4k=2kappa sigma to the curvature delta. The Hamilton
constraint therefore includes exactly the rest membrane energy.
All scalar and wall momentum densities vanish on this slice. Thus
the preparation supplies compatible Lorentzian initial data, including
the surface source, rather than only formal Euclidean profiles.

For eta=1, tau_W=2.70634531321651 and R_W=4.17610904394739.
The wall is timelike; its areal-radius derivative exceeding one is not
a local superluminal velocity. In its rest frame the scalar energy
flux is zero for this symmetric prepared solution. The spatial layer
and its stresses are present in the geometry.
Likewise d^2R_W/dtau_W^2 is an areal-radius second derivative, not the
invariant proper acceleration; the two side accelerations are K_i/3
and K_o/3 in the adopted orientation.

Near the child pole, xi=i tau gives an open FLRW chart,

$$
ds^2=-d\tau^2+b(\tau)^2dH_3^2,\quad
\ddot\chi+3\frac{\dot b}{b}\dot\chi+\omega^2\chi=0,\quad
\dot b^2=1+\frac{\kappa b^2}{3}(\rho_i+E),
\tag{S8}
$$
$$
E=\frac12(\dot\chi^2+\omega^2\chi^2),\quad
b(0)=0,\quad\dot b(0)=1,\quad
\chi(0)=\chi_i(0),\quad\dot\chi(0)=0.
$$

The coordinate origin b=0 is regular, with
b=tau+kappa U_0 tau^3/18+... and
chi=chi_i(0)(1-omega^2 tau^2/8+...).
It is not a divergent-curvature big bang.

For an expanding child, E_dot=-3H chi_dot^2<=0. Set
h_i^2=kappa rho_i/3 and h_+^2=kappa(rho_i+E_0)/3.
Comparison of S8 gives

$$
\frac{\sinh(h_i\tau)}{h_i}
\le b(\tau)\le
\frac{\sinh(h_+\tau)}{h_+}.
\tag{S9}
$$

Together with bounded chi and chi_dot this establishes future
continuation of the homogeneous chart for every finite tau, conditional
on the regular prepared solution. No global geodesic-completeness or
quantum-formation probability is inferred.

The pole amplitude is extracted from an independent regular mode:
u=chi'/chi obeys u'=omega^2-3(a'/a)u-u^2, with
u=omega^2 xi/4+... near the pole. Integrating its logarithmic gain avoids
subtracting nearly equal field minima. It yields

$$
\chi_i(0)=-6.9370429853\,10^{-25},\qquad
E_0=9.6245130760\,10^{-47}.
$$

The gain from child pole to wall is 51.1210048481. The large homogeneous
excitation .09375 from the earlier frozen transfer is therefore not
realized by this preparation. E_0/.09375=1.0266e-45, and S9 limits
the child excitation to E/rho_i<=7.70e-46 throughout this homogeneous
future branch. This diagnoses the physical consequence, not a demand
for further tiny-tail certification.

For an illustrative local response at tau=1, evolve the normalized
scalar on the open de Sitter background. The actual Friedmann correction
at fixed b is bounded by kappa E_0/3=3.21e-47. The computed excitation
is 2.84e-50; the local free-fall tidal-trace increment is -4.44e-50
against vacuum trace .125. These are approximation-controlled reference
values, not observational error bars or a measured dark-matter abundance.

## Verification and what is still open

The [compact record](../results/summaries/nsc-atlas-scalar-preparation-v1.json)
binds source, producer, tests and figure. Six identities vanish exactly.
The two collocation tolerances change R_* by 3.81e-12 and the logarithmic
pole-transfer amplitude by 9.0e-11 relatively. The independent Hamilton
constraint is below 6.9e-13, the sampled scalar equation below 2.2e-9,
and the weighted tangential Einstein equation below 1.2e-9. Wall matching
and force residuals are below 1e-14. These are controlled numerical checks,
not an interval-certified nonlinear existence or uniqueness proof.
The zero-mixing control recovers the imported pure-flux cap solution.

```
python scripts/run_python.py scripts/check_nsc_atlas_scalar_preparation.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_scalar_preparation.py -q
python scripts/run_python.py scripts/plot_nsc_atlas_scalar_preparation.py
```

Regular Euclidean preparation is a stated state-selection assumption.
The negative-mode count, contour, determinant and dominant transition
probability have not been established. A black-hole collapse history,
P002's dust exterior and a second nested generation are not part of
this solved boundary problem. Absolute vacuum scale and couplings remain
inputs, and the infinite inherited operator is not identified by S5 alone.

The useful next question is a physically specified excited-parent or
nonequilibrium preparation if a substantial child matter component is
claimed. This regular branch should remain a completed reference.
Its tiny homogeneous excitation must not be replaced by the earlier
instant-transfer number, nor should another refinement of this branch
substitute for a different preparation mechanism.
