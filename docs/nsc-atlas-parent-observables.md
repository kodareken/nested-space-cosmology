# Surrounding geometry, an expanding parent and observable gravity

This atlas entry connects the existing nested body to an expanding
four-dimensional parent, then compares two published ways a surrounding
five-dimensional geometry produces local gravity. A fourth calculation
tests whether the existing wall can select the cap curvature by stable
equilibrium. These are geometries with specified laws, not interchangeable
descriptions of one solution.

The book's purpose is to make these constructions legible together:
draw the geometry, name its source, exhibit the connecting equation, and
follow an observer. Published solutions are building blocks credited to
their authors. The project contribution is their checked attachment or
comparison with the owned atlas, and the resulting inference questions.

![Calculated comparison of surrounding geometries and wall balance](figures/nsc-atlas-parent-observables.png)

The four panels are calculated from the equations below. A and C use the
current four-dimensional family. B and D use distinct five-dimensional
backgrounds; their successes have not been combined in one model.

## A. Replace the isolated end by a cosmological parent

**Borrowed construction.** Einstein--Straus matching joins a spherical
vacuum region to comoving dust without a material surface layer.
[Mars, Mena and Vera, Section 2](https://arxiv.org/html/1307.4371#S2)
provide the matching data and its generalizations. Here the new
attachment uses the original nested body's increasing-radius right
Schwarzschild exterior, beyond its throat. The cap, both material walls,
host and bridge retain their checked local histories.

This is an alternative completion of the earlier AE end. It does not
turn the [compact S3 completion](nsc-atlas-compact-parent.md) into an open
space, or replace its causal structure during evolution.

Set G=c=1, keep M=.4 ell and the original ell=2/k, and put the new join
at R_b(0)=4 ell. The parent is a flat dust FLRW solution with Lambda=0.
The chosen radius and matter law are additional initial/model data.
Write its metric and the collar metric as

$$
ds_F^2=-d\tau^2+a(\tau)^2(d\chi^2+\chi^2d\Omega^2),\qquad
ds_S^2=-f\,dT^2+f^{-1}dR^2+R^2d\Omega^2,\quad f=1-\frac{2M}{R}.
\tag{P1}
$$

At chi=chi_b=4 ell, let R_b=a chi_b. Equality of the angular extrinsic
curvature, with the normal directed toward increasing R, gives

$$
\frac{\sqrt{\dot R_b^2+f_b}}{R_b}=\frac1{R_b},
\qquad
\dot R_b^2=\frac{2M}{R_b},\qquad \dot T_b=\frac1{f_b}.
\tag{P2}
$$

The boundary is geodesic: R_b''=-M/R_b^2. Its temporal extrinsic
curvature is zero on both sides. The induced metrics are both
-d tau^2+R_b^2 dOmega^2. Thus both fundamental forms agree, not only
the angular equation. Normal pressure and flux vanish at this join.
There is no added shell stress.

The parent Friedmann equation then implies the mass-compensation identity

$$
M=\frac{4\pi}{3}\rho_d R_b^3,\qquad
H_0^2=\frac{2M}{R_b(0)^3},\qquad
a(\tau)=\left(1+\frac32H_0\tau\right)^{2/3}.
\tag{P3}
$$

M is the already measured exterior mass, including the gravitational
effect of its attached interior. It replaces the excised dust mass.
The proper integral of energy inside a nontrivial nested geometry is
not an additional exterior mass to add to P3. Positive energy, binding,
wall stresses and the chosen hypersurface all enter that distinction.

For the fixed geometry,

$$
H_0\ell=.111803398875,\quad
\rho_d(0)\ell^2=.00149207759149,\quad
\dot R_b(0)=.447213595500.
\tag{P4}
$$

At tau=.01 ell, R_b=4.004470886886 ell. At tau=4 ell,
R_b=5.632222928894 ell. The exterior/cosmological equations are exact
on this interval. Their continuation does not claim that the separately
studied inner throat remains regular for four ell of its own normal time.

### A compatible initial surface, with its clock choice exposed

The old right collar has T=constant. A cosmic-time surface at the new
join instead has dT/dR=v/f, where v=sqrt(2M/R). Connect them by
dT/dR=w(R)v/f: w=0 on R<=ell, w=1 on R>=3.5 ell, and
w=10s^3-15s^4+6s^5 between, s=(R-ell)/(2.5 ell).
This C2 splice has induced radial metric

$$
g_{RR}^{\Sigma}=\frac{1-w^2v^2}{f}>0,\qquad
V_\Sigma=\frac{\sqrt f\,wv}{\sqrt{1-w^2v^2}}.
\tag{P5}
$$

It preserves the old inner collar, becomes the flat cosmic slice near
the join, and remains spacelike because R>2M and 0<=w<=1.
The endpoint derivatives agree; no source is inserted by choosing this
surface. Its induced metric and extrinsic curvature satisfy the vacuum
Hamiltonian constraint. Choosing the integration constant fixes one
relative chart origin. This changes the outer slicing, so the earlier
whole-slice normal averages are not reused as averages of this parent.

If instead the old zero-velocity normal were imposed at the join itself,
an expanding flat dust boundary would be incompatible with it. A
closed-dust turnaround is a different possible match. P5 explicitly
changes that outer initial slice while preserving the interior collar;
it does not identify a static normal with a moving comoving boundary.

### An observer consequence of this same join

A static emitter at R_e=2 ell emits at the collar Killing time of the
join's tau=0 event. This is a declared synchronization convention.
Reception occurs on the outward-moving join. Both events and the ray
are on the parent side of the bridge; this is not a message from the
causally separated cap.

For rstar=R+2M log(R/(2M)-1), the reception equation is

$$
\int_0^{\tau_r}\frac{d\tau}{1+v_b(\tau)}
=r_*(R_b(0))-r_*(R_e).
\tag{P6}
$$

The constant photon Killing energy and the receiver's unit tangent give

$$
\frac{\omega_r}{\omega_e}
=\frac{\sqrt{f_e}}{1+v_b(\tau_r)},\qquad
\tau_r/\ell=3.921887982097,\quad 1+z=1.778824923298.
\tag{P7}
$$

A separate null ODE gives the same arrival within 1.3e-12 ell; the
32/64 quadrature difference is 4.4e-16 ell. Differentiating reception
time with respect to the emitter's proper time checks the frequency
ratio independently. These are numerical checks of a conditional
instrument prediction, not observational error bars or a fit to data.
The physical seconds require an independently supplied ell/c.

This parent decelerates with q=1/2. Its cosmological Lambda contribution
from the attachment is exactly zero. In contrast its local Ricci scalar
divided by four initially equals .009375/ell^2, from dust. This example
again shows why a nonzero curvature trace alone does not identify vacuum.

## B. A surrounding bulk can induce local curvature

**Separate candidate and imported equations.** A timelike four-dimensional
brane in a five-dimensional Einstein bulk differs from the existing
two-sphere wall in four-dimensional spacetime. The bulk has geometric
Lambda_5=-6/L_5^2, Z2 reflection, tension lambda and confined stress T.
The [Shiromizu--Maeda--Sasaki projection, equations 15--22](https://arxiv.org/html/gr-qc/9910076v3)
and [Maartens--Koyama, Section 3](https://arxiv.org/html/1004.3962)
give

$$
G_{\mu\nu}=-\Lambda_4 q_{\mu\nu}+\kappa_4^2 T_{\mu\nu}
 +\kappa_5^4\pi_{\mu\nu}-E_{\mu\nu},\qquad
\kappa_4^2=\frac{\kappa_5^4\lambda}{6},\qquad
\Lambda_4=\frac{\Lambda_5}{2}+\frac{\kappa_5^4\lambda^2}{12}.
\tag{P8}
$$

Here E is the projected bulk Weyl tensor, with zero trace.
The quadratic term is

$$
\pi_{\mu\nu}=-\frac14T_{\mu\alpha}T_\nu{}^\alpha
+\frac1{12}T T_{\mu\nu}
+\frac18q_{\mu\nu}T_{\alpha\beta}T^{\alpha\beta}
-\frac1{24}q_{\mu\nu}T^2 .
\tag{P9}
$$

The calculator substitutes the Israel curvature directly into the
contracted Gauss equation and compares every tensor component with P8,
including momentum and anisotropic stress. Conservation retains
D^mu E_munu=kappa_5^4 D^mu pi_munu when confined matter is conserved
and the bulk contains only its cosmological term. Arbitrary choices of
T and E in an algebra check do not construct their bulk evolution.

An explicit vacuum example closes that missing geometry:

$$
ds_5^2=dy^2+A(y)^2ds_{\mathrm{dS}_4}^2,\qquad
A(y)=\cosh(y/L_5)-s\sinh(|y|/L_5),\qquad
s=\frac{\kappa_5^2\lambda L_5}{6}.
\tag{P10}
$$

For s>1, use 0<=|y|<L_5 atanh(1/s) and A(0)=1.
Direct differentiation gives A''=A/L_5^2 and
A'^2-A^2/L_5^2=(s^2-1)/L_5^2. The normal and tangential bulk
sectional curvatures both equal -1/L_5^2. Israel fixes
A'(0+)/A(0)=-kappa_5^2 lambda/6. Thus

$$
\boxed{\Lambda_4 L_5^2=3(s^2-1).}
\tag{P11}
$$

The fixed demonstration s=1.1 gives Lambda_4 L_5^2=.63 and
H_4 L_5=.458257569496. The bulk coordinate horizon is at
1.522261218862 L_5. At s=1 the induced geometry is flat.
This is a positive local-curvature construction with a full bulk metric,
not a claimed selection of s or L_5 from NSC inheritance.
The reference length L_5 is distinct from the original membrane ell.

This decomposition makes a useful distinction calculable. In a homogeneous
and isotropic brane the Weyl contribution, if separately conserved,
has p_E=rho_E/3 and rho_E proportional to a^-4. It supplies no nonzero
pure vacuum tensor. For confined dust the quadratic correction has
rho_quad=p_quad=rho^2/(2 lambda). The constant term in this example
comes from bulk curvature and tension, not either of those components.

## C. Matter outside the observed region: a lensing test

[Garriga--Tanaka, equations 26--28](https://arxiv.org/html/hep-th/9911055v3)
supply a different explicit example: two flat branes, one with positive
and one with negative tension, with an unstabilized massless radion.
In their weak, static, zero-mode limit, a nonrelativistic source on
the other brane contributes through the tensor
T_munu-q_munu T/3. Its time component is 2rho/3, each spatial diagonal
component rho/3. With

$$
ds^2=-(1+2\Phi)dt^2+(1-2\Psi)d{\bf x}^2,
\qquad \gamma=\Psi/\Phi,
\tag{P12}
$$

this gives gamma_shadow=1/2. Integrating the weak null deflection along
a straight ray at impact b gives

$$
\alpha=\frac{2(1+\gamma)G_4M_{\rm dyn}}{bc^2},\qquad
\frac{\alpha_{\rm shadow}}{\alpha_{\rm GR}}=\frac34
\quad\hbox{at equal dynamical mass}.
\tag{P13}
$$

The familiar 25% difference is attributed to Garriga--Tanaka, not
claimed as an NSC discovery. It is valuable because it converts
"outside matter gravitates here" into two distinct measurements.

For an observer on the positive brane, let f be the shadow fraction
of the dynamical potential and assume common radial profiles. With
d/L_5=5,

$$
\gamma_v=\frac{3-e^{-2d/L_5}}{3+e^{-2d/L_5}},\qquad
\gamma=(1-f)\gamma_v+\frac f2,\qquad
\frac{\alpha}{\alpha_{\rm GR}}=\frac{1+\gamma}{2}.
\tag{P14}
$$

The calculated lensing ratios at f=0,.5,1 are .999984866919,
.874992433460 and .75. The fraction f is not Omega_c. Stabilizing
the radion, adding massive bulk modes or changing the background
requires its own response calculation. The flat two-brane solution
and the detuned dS brane in B have not been combined into one geometry.

## D. Does the actual wall select its own cap scale?

First test the original massless de Sitter cap joined to Schwarzschild.
The [owned tension wall](../src/recursive_horizons/nsc_atlas_nested_slice.py)
uses a constant positive sigma and surface pressure p_s=-sigma,
with k=4 pi G sigma in c=1 units. Reusing the published
[vacuum-bubble construction](https://arxiv.org/html/0711.2198v1#S3),
its signed angular junction gives

$$
\beta_-=\frac{M}{kR^2}+\frac{(k^2-H^2)R}{2k},\qquad
\beta_+=\frac{M}{kR^2}-\frac{(k^2+H^2)R}{2k},
\tag{P15}
$$

$$
\dot R^2+V(R)=0,\qquad
V=1-\frac{M^2}{k^2R^4}
-\frac{M(1-H^2/k^2)}{R}
-\frac{(k^2+H^2)^2R^2}{4k^2}.
\tag{P16}
$$

This is the existing equation, checked directly against its implementation.
The proposed stable-equilibrium closure has a short exact obstruction.
Write V=1-A/R^4-B/R-CR^2, where A=M^2/k^2>0 and
C=(k^2+H^2)^2/(4k^2)>0. At every stationary point,
B/R^3=2C-4A/R^6. Consequently

$$
\left.V''\right|_{V'=0}=-\frac{12A}{R^6}-6C<0.
\tag{P17}
$$

No static solution of this constant-tension law is a radial stable
minimum. This is an identity, not a failed grid search. Near a static
solution delta R''=-V'' delta R/2 has an exponentially growing mode.
At fixed zero energy, that unstable displacement is allowed by the
constraint. The result does not exclude stability for a different
material law.

The entire static family can also be stated without scanning. Put
u=H^2/k^2, x=kR, mu=kM and q=mu/x^3. Then

$$
q=\frac{u-1+\sqrt{9u^2+14u+9}}8,\qquad
x^2=\frac1{3[(1+u)^2/4-q^2]},\qquad \mu_{\rm crit}=q x^3.
\tag{P18}
$$

For u>=0 these solutions have beta_->0 and beta_+<0. They satisfy
the signed temporal junction too. Thus their wormhole-side exterior
orientation is retained; they are not ordinary gravastar joins.
Along the critical curve,
V_u=x beta_+<0 and V_mu=-2 beta_-/x^2<0, so
d mu_crit/du=-V_u/V_mu<0. At u=0, mu_crit=16/27.
The existing state has mu=.8>16/27: at its fixed M and k it has
no static solution for any nonnegative H^2.
The stationary equation has only one positive radius, and its value is
V_max=1-(mu/mu_crit(u))^(2/3). Thus above the critical mass the
potential is negative at every radius, not merely at the original wall.

The actual nested outer wall has inner Kottler mass m=.01 ell.
Its potential has the same form with
A=(M-m)^2/k^2 and B=M+m-(M-m)u; C is unchanged.
The instability identity P17 therefore still applies. Its dimensionless
potential equals the massless-cap expression at mu=k(M-m)=.78,
minus 2km/x=.04/x. The massless expression is already everywhere
negative for .78>16/27, so this actual outer wall also has no static
root for any u>=0. The massless critical curve P18 is used as a
comparison bound, not misidentified as the nested wall's critical curve.

For illustration, keeping u=.25 but changing the mass to the critical
one gives x=1.126996043952, mu=.512485700729 and V_xx=-3.881941479272.
Its Lambda ell^2=12u=3 is still an input ratio, and the wall is unstable.

The scale map

$$
(R,M,k,H)\longmapsto(aR,aM,k/a,H/a)
\tag{P19}
$$

preserves the same dimensionless wall law for every positive a.
It leaves kM and H/k unchanged, without selecting either value.
This exhibits the freedom left by scale-compatible inheritance of
these equations. It does not prove that the full operator relation
T*D=D can never constrain them; the physical map from that operator
and inherited state to these material parameters remains to be supplied.
The relevant mass invariant is kM, not the dimensionful ratio M/sigma.

### Why gravastars remain useful building blocks

[Mazur--Mottola](https://arxiv.org/html/gr-qc/0109035) use a de Sitter
core and Schwarzschild exterior with a material transition layer.
[Visser--Wiltshire](https://arxiv.org/html/gr-qc/0310107), Sections 2--5,
make radial stability calculable for a thin shell with an equation of
state. Their shell is not restricted to p_s=-sigma. In a no-flux
surface description its conservation law is

$$
\frac{d\sigma}{dR}=-\frac{2(\sigma+p_s)}R.
\tag{P20}
$$

An inherited material/state law could therefore change P16 and its
stability. Choosing that law solely to enforce a desired potential
would supply the answer as input. The useful next connection is to
derive its pressure response from inherited degrees of freedom, then
use the existing gravastar stability machinery.

The throat null-stress condition also has established antecedents in
[Morris--Thorne](https://doi.org/10.1119/1.15620). The earlier compact
note applied the radial equations at a specified dynamical turn; it
did not discover null-energy violation for traversable throats.
[Easson--Brandenberger](https://arxiv.org/abs/hep-th/0103019) supplies
a black-hole-interior cosmology comparison. Those constructions do not
by themselves determine the present wall's inherited equation of state.

## E. Growing black-hole masses and the missing energy equation

The [cosmological-coupling proposal](https://arxiv.org/abs/2302.07878)
is useful as a test of source accounting. For fixed comoving number,
M proportional to a^n gives rho_B proportional to a^(n-3).
The conservation equation is

$$
\dot\rho_B+3H(\rho_B+p_B)=Q.
\tag{P21}
$$

If this is a separately conserved effective component, Q=0 implies
w_B=-n/3. If it is pressureless mass receiving energy, p_B=0 instead
gives Q=nH rho_B. The same mass scaling cannot decide between these
physical descriptions.

A finite explicit counterexample takes n=3, rho_B=.1 and a pressureless
reservoir rho_X=a^-3-.1, in units where the initial total density is one.
For 1<=a<=2 the reservoir is positive, the two transfers cancel,
rho_total=a^-3 and q_total=1/2. Black-hole density stays constant
without total vacuum pressure or acceleration. This homogeneous
conservation example is not a microscopic black-hole accretion solution.
It shows exactly which further source law must be derived.
Published [binary-merger tests](https://arxiv.org/abs/2306.08199) supply
additional independent consequences of particular coupling laws.

## The next connection and reproduction

Completed here: a cosmological attachment to the retained exterior,
its clock prediction, explicit bulk-curvature and shadow-gravity
comparisons, and the algebraic decision on constant-tension equilibrium.
The next NSC-specific unknown is an inherited stress/state law or
gravitational reduction that selects dimensionless parameters and changes
local observables. Another quiet throat integral or a finer grid cannot
select that law.

No Omega_Lambda, Omega_c, Omega_b or baryon asymmetry has been predicted
by these comparisons. Those are distinct future outputs of a common
cosmological and thermal history. The supplied cap H, the new dust
parent density and the separate bulk detuning are all recorded inputs.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_parent_observables.py
python scripts/run_python.py scripts/plot_nsc_atlas_parent_observables.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_parent_observables.py -q
```

The [calculator](../src/recursive_horizons/nsc_atlas_parent_observables.py),
[23 checks](../tests/test_nsc_atlas_parent_observables.py) and
[compact scalar record](../results/summaries/nsc-atlas-parent-observables-v1.json)
retain equations, domains, source hashes and numerical indicators.
No field grid is saved. The calculation uses 882 parent-state evaluations
and completes below a tenth of a process-CPU second; the recorded Windows
CPU timer can round this short run to zero. This analytic workload benefits
from parallel independent review, not numerical grid expansion.
