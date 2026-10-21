# A finite parent scalar source with its gravity and wall response

This entry constructs a gravitationally active preparation using the
same scalar, four-form and charged membrane as the
[regular reference](nsc-atlas-scalar-preparation.md). The
[transport map](nsc-atlas-parent-transport.md) remains a causal control;
the source here is no longer a test body. The calculation solves a
time-symmetric Cauchy slice, not a finite interval of its evolution.

For a specified compact parent field displacement of amplitude .02,
the coupled solution gives wall radius 2.764905208 L_ref, compared
with 2.706345313 L_ref without that displacement. Its initial areal
acceleration increases by 4.569211 percent. The proper scalar excess
is 1.112671383 in the normalized energy units defined below.

For a second-order smooth initial corner, the mathematical connection
is a scalar-layer correction. Solving gravity while freezing the old
layer fails that classical compatibility condition. The second-normal correction
is determined by the conformal geometry, not an independently fitted
impulse. After including it, both sides agree on the wall's scalar
acceleration and initial energy-flux derivative.

A second preparation keeps the child initially at its minimum while
storing excitation uniformly in the finite parent cap. Its round
gravitational constraints close exactly. Allowing a finite-energy weak
initial corner gives a conditional inward-flux coefficient of 18.75
and a linked pressure-change coefficient of the same value. This
short-time result is derived below; a full coupled evolution is not
established by the tangent calculation.

## Rethinking decision and scope

A full moving-wall evolution, an additive collision preparation and
a finite constrained scalar preparation were compared. Full evolution
first requires compatible data. Additive free passage in the present
quadratic model lands on another minimum. More test-body calculations
do not include the source's gravity. The chosen preparation supplies
that missing initial-data connection with one positive conformal
factor and one determined wall coefficient.

The free state data are a packet amplitude, its position and width,
a wall-collar width, and time symmetry. The action and its couplings
stay fixed. The conformal factor removes the ambiguity of adding
matter to an unchanged metric; the wall coefficient removes the
incompatibility introduced by freezing the scalar layer. They are
solved together. The shape and amplitude of this preparation are
supplied examples, not a derivation of the parent's previous history.

The bulk method reuses the scalar Hamiltonian constraint and spatial
curvature transformation in
[Choquet-Bruhat, Isenberg and Pollack, equations 4 and 8](https://arxiv.org/html/gr-qc/0610045).
The material-wall equations require the additional derivation below.
Their smooth compact-manifold existence results are not a theorem for
this distributional membrane. Israel's full tensor condition, including
the embedding, is also the starting point of
[Mukohyama's junction analysis, Section II](https://arxiv.org/html/hep-th/0006146).

The calculation admits amplitudes 0, .01 and .02, one rejected frozen-layer
control, and one coarse/fine comparison. Each boundary solve has at most
10000 nodes and 3000 vector RHS calls; the producer has a 120 CPU-second
budget. These are adaptive constraint solves, not a spacetime grid scan.
Failed constraints or matching stop this preparation. Successful initial
data become input to evolution; they do not justify a new precision campaign.

## The finite state and its Hamiltonian constraint

Let the original time-symmetric slice be

$$
h_0=d\xi^2+a^2d\Omega_2^2,\qquad
K_{ij}=0,\qquad \Pi=u^\mu\nabla_\mu\phi=0.
\tag{C1}
$$

Both caps retain their own pole-to-wall coordinates. The common wall
normal points from child to parent; it is positive child xi and negative
parent xi. The original scalar is phi_0 and the physical spatial metric
to be solved is

$$
h=\Psi^4h_0,\qquad \phi=\phi_0+f+g.
\tag{C2}
$$

The child has f=g=0. In the parent let d=L_o-xi_o. The prescribed
packet is

$$
f=A(1-z^2)^6,\quad z=\frac{d-.3}{.12},\quad |z|<1,
\qquad f=0\ \hbox{otherwise}.
\tag{C3}
$$

It is C5 at its compact-support endpoints and initially has zero
conjugate momentum. It is therefore a time-symmetric excitation, not
a prepared purely inward travelling wave. Its support is within the
reference's accessible parent collar. New causal arrival times must
nevertheless be computed in the evolved sourced geometry.

Write U_0=U(phi_0,q) and U=U(phi,q), with the unchanged quadratic
potential. The exact bulk constraint is

$$
\Delta_0\Psi=
\frac{\kappa}{8}\left[(\phi_0')^2-(\phi')^2\right]\Psi
+\frac{\kappa U_0}{4}\Psi-\frac{\kappa U}{4}\Psi^5,
\quad \Delta_0\Psi=\Psi''+\frac{2a'}{a}\Psi'.
\tag{C4}
$$

This follows directly from R[h]=Psi^-5(-8 Delta_0 Psi+R[h_0] Psi)
and R[h]=kappa(|D phi|_h^2+2U). The momentum constraint is exactly
zero because K_ij and Pi both vanish. At each pole Psi'=0.

The physical areal radius is R=Psi^2 a and proper radial distance
obeys dl=Psi^2 dxi. Thus beta=dR/dl=a'+2a Psi'/Psi. Keeping the
same material tension requires

$$
[\Psi]=0,\qquad
[\Psi']=\frac{k}{2}(\Psi-\Psi^3),\qquad
k=\frac{\kappa\sigma}{2}.
\tag{C5}
$$

Brackets use the common normal. In stored cap coordinates the second
condition is Psi_i'+Psi_o'=k(Psi^3-Psi)/2. This retains the actual
wall energy; setting both derivatives continuous would change the source.

## Scalar matching for a second-order smooth initial corner

Keep phi and its first seed-normal derivative at the wall equal to
the original values. Write s=phi_0,n=-.21559588825 there. Continuity
of phi and the physical phi_n=s/Psi_W^2 is then preserved. However,
the wave equation at a timelike interface has the normal form

$$
\phi_{;nn}+K\phi_n+\Box_\Sigma\phi=U_\phi.
\tag{C6}
$$

The wall's intrinsic scalar is single-valued. With [K]=-3k and
[U_phi]=mu e, its second-order compatibility condition is

$$
[\phi_{;nn}]-3k\phi_n-\mu e=0.
\tag{C7}
$$

For a frozen layer, the conformal transformation changes the physical
normal derivatives but does not change the charge jump. The resulting
residual is

$$
\mu e(\Psi_W^{-4}-1)
+2ks(\Psi_W^{-4}-\Psi_W^{-2}).
\tag{C8}
$$

Multiplying by Psi_W^4 factors this as
-(Psi_W^2-1)[mu e(Psi_W^2+1)+2ks]. For the fixed couplings, the
second factor is positive for every positive Psi_W. Hence a frozen
layer can satisfy C7 only at Psi_W=1. This is an algebraic obstruction,
not a discretization error.

To supply the missing derivative without changing phi or phi_n at the
wall, use a parent collar of seed width delta=.08:

$$
g(d)=\frac{\lambda}{2}d^2(1-d/\delta)^6,
\quad 0\le d<\delta,\qquad g=0\ \hbox{otherwise}.
\tag{C9}
$$

It has g(0)=g'(0)=0 and g''(0)=lambda. C7 determines

$$
\lambda=\mu e(\Psi_W^4-1)+2ks(\Psi_W^2-1).
\tag{C10}
$$

C4, C5 and C10 are solved simultaneously. The finite primary value is
lambda=.3779028621. The collar shape and width are part of this
declared preparation; C10 fixes its second derivative, not a unique
global preparation mechanism.

## Instantaneous gravitational and scalar response

At the momentarily resting wall, let B_j=beta_j/R and let A_j denote
the wall's normal proper acceleration evaluated from side j, always
with the common normal. The actual-field pressure is
P=e(q_p+mu phi_W-e/2)=.154440177981835. Angular Israel gives
B_o-B_i=-k. Full temporal Israel and normal force then determine

$$
\overline{A}=\frac{P}{\sigma}-(B_i+B_o),\qquad
A_i=\overline{A}+\frac{k}{2},\qquad
A_o=\overline{A}-\frac{k}{2}.
\tag{C11}
$$

The scalar equation and the radial Einstein equation determine
two independently evaluable initial wall derivatives:

$$
\frac{d^2R}{d\tau_W^2}=
-\frac{m_j}{R^2}-\frac{\kappa R p_{r,j}}{2}+A_j\beta_j,
\quad m_j=\frac{R}{2}(1-\beta_j^2),
\tag{C12}
$$
$$
\frac{d^2\phi_W}{d\tau_W^2}=
\phi_{;nn,j}+\frac{2\beta_j}{R}\phi_n-U_{\phi,j}+A_j\phi_n.
\tag{C13}
$$

For C12, the spherical warped-product Einstein identity gives
G_nn=-2 R_{;uu}/R-2m/R^3. Adding A_j beta_j converts the timelike
Hessian into the wall proper-time derivative. Equivalently, the right
side follows from R_{;uu}=-m/R^2-kappa R p_r/2. It is an areal-radius acceleration,
distinct from the invariant proper acceleration A_j.

The two C12 values agree exactly when the angular, temporal and
force conditions hold. The difference of the two C13 values is
precisely C7. This supplies independent checks of the coupled
geometric and scalar matching rather than a scalar-only probe.

Energy flux in the outward common-normal direction is
F=-T(n,u)=-phi_n d phi_W/dtau_W. Initially F=0, and

$$
\frac{dF}{d\tau_W}(0)=
-\phi_n\frac{d^2\phi_W}{d\tau_W^2}(0).
\tag{C14}
$$

Both sides give the same value, consistent with constant membrane
tension and continuous tangential energy exchange. For the primary
source it is +.03803680229: the initial flux slope points from child
to parent. This state has not demonstrated inward matter deposition.
Later arrival of the specified parent packet requires actual evolution.

## Finite results and independent checks

All quantities below use the existing normalized variables. Energy
is the proper integral of the difference between the new scalar and
phi_0 on the **same sourced metric**, including interference:

$$
E_{\rm excess}=4\pi\int a^2\Psi^6\left[
\frac{\delta\phi'^2+2\phi_0'\delta\phi'}{2\Psi^4}
+\frac{\omega^2}{2}(\delta\phi^2+2\chi_0\delta\phi)
\right]d\xi_o.
\tag{C15}
$$

Here delta phi=f+g and chi_0=phi_0-phi_*(q_p). This is neither an
ADM energy nor a conserved total energy for a closed universe.
For the primary source its positive quadratic part is 1.112605397
and its interference part is .000065985551. No SI ruler is selected.

| Supplied amplitude | Proper scalar excess | R/L_ref | Initial L_ref R''/c^2 | Initial outward flux slope |
| --- | --- | --- | --- | --- |
| 0 | 0 | 2.706345313 | .369501998 | 0 |
| .01 | .270536650 | 2.721705563 | .374020762 | .010375353 |
| .02 | 1.112671383 | 2.764905208 | .386385323 | .038036802 |

For A=.02, Psi_W=1.010761098, and the conformal factor stays
between .979049929 and 1.012986587 on the checked slice. The wall
radius increases by 2.163800 percent. C12 gives .38638532338659
from each side; C13 gives .18024389453205 from each side. All
reported matching discrepancies are below 3e-15.

The frozen-layer control satisfies the gravitational constraint but
has scalar compatibility residual -.362064088. Adding C10 removes
that incompatibility. It does not introduce a new wall force or
change the membrane's scalar surface coupling.

An independent integrated mass relation,

$$
\frac{dm}{dl}=\frac{\kappa}{2}R^2\rho\frac{dR}{dl},
\tag{C16}
$$

recovers wall masses .44414721424740 and .96826983393285 from
their respective regular poles, with discrepancies below 7e-13.
Direct spatial-curvature checks give Hamilton residuals below
3e-8. Changing both background and source tolerances from 1e-7
to 1e-9 changes radius by 4.36e-10 and the determined wall coefficient
by 2.82e-9. The retained numerical comparison tolerance is 2e-7.
Five symbolic identities, the zero-source limit and the rejected
frozen-layer control accompany the numerical checks. These are not
interval-certified existence, uniqueness or error-enclosure theorems.

![Finite parent source and its gravitational geometry](figures/nsc-atlas-parent-source.png)

## A weak initial corner can carry a real inward pulse

Second-order smoothness at the initial interface is a regularity choice,
not an additional physical law. C8 rejects a frozen layer as C2 corner
data; it does not exclude a causal finite-energy weak development with
an initial transient. Treating smoothness as compulsory would discard
a different useful preparation before examining its energy transfer.

Set Pi=K_ij=0 and choose one uniform field on the entire compact slice,
phi=phi_*(q_c). The child starts exactly at its minimum. The parent
is displaced from its own minimum by mu e/omega^2, so its excitation
density is 3/32 with the unchanged benchmark couplings. The initial
densities are rho_i=1/8 and rho_p=3/8. The latter is an excited
initial density, not the parent's future cosmological constant.

The Hamiltonian constraint gives round caps with h_i^2=1/24 and
h_p^2=1/8. The angular and force conditions close exactly:

$$
B_i=\frac{7}{15},\qquad B_o=\frac{11}{30},\qquad
R_*^2=\frac{1800}{467},\qquad
R_*=1.96325997301,\qquad P_0=\frac14.
\tag{C17}
$$

Here B_j=beta_j/R_* and the initial A_j equal B_j. The initial
areal acceleration is 1/R_*=.509356893; this is a turning instant,
not static equilibrium. Both poles are regular. The parent excitation
has finite proper integral 38.30744838 in the declared units. Scalar
value and first normal derivative match, but the initial C7 residual
is -mu e. This seed has a weak initial corner and is distinct from
the second-order smooth packet above.

**Conditional local transfer statement.** Suppose these data have a
causal local Einstein/scalar/membrane development with continuous scalar
value and normal derivative, bounded one-sided geometric coefficients
near the wall, and no additional initial scalar wave. Existence of that
full weak development is not proved here. In local wall proper-time
and normal coordinates (tau,x), parent x>0 and child x<0, put
S=mu e=5 sqrt(3). The rescaled leading scalar equation is

$$
(\partial_\tau^2-\partial_x^2)\chi=-S H(x),\qquad
\chi(0,x)=\partial_\tau\chi(0,x)=0.
\tag{C18}
$$

Here H is the step function. Retarded integration gives, within the
two sides of the initial cone,

$$
\chi_i=-\frac{S}{4}(\tau+x)^2\quad(-\tau<x<0),\qquad
\chi_p=-\frac{S}{4}(\tau^2+2\tau x-x^2)\quad(0<x<\tau).
\tag{C19}
$$

It is zero ahead of the child front; beyond the parent front the
leading displacement is -S tau^2/2. The value and first derivatives
match at x=0 and at the fronts. The initial scalar momentum is zero;
the different potential slopes generate these derivatives dynamically.

The leading statement follows from a local scaling, not a flat-space
replacement of the full geometry. On a cone of size tau, bounded
one-sided connection coefficients multiply field derivatives of order
tau; the mass term multiplies a displacement of order tau^2. Deviations
of the principal metric from its tangent value are of order tau.
Retarded integration therefore puts these corrections at order tau^3
in the field and order tau^2 in its first derivatives, subject to the
stated regularity assumptions. Finite wall acceleration and source
backreaction cannot change the displayed leading coefficient within
that class. No numerical bound on the remainder is supplied.

At the wall this gives actual inward energy flow, conditional on the
development just specified:

$$
F_{\rm in}=T(n,u)=\frac{S^2}{4}\tau^2+O(\tau^3),\qquad
\frac{E_{\rm child}}{\mathcal A_0}
=\frac{S^2}{12}\tau^3+O(\tau^4).
\tag{C20}
$$

The energy is the excitation above the initially quiet child, integrated
through its initial inward layer; A_0=4 pi R_*^2. Geometry and the
moving area's corrections enter beyond the leading term. Numerically,

$$
\lim_{\tau\to0^+}\frac{F_{\rm in}}{\tau^2}=18.75,\qquad
\lim_{\tau\to0^+}\frac{E_{\rm child}/\mathcal A_0}{\tau^3}=6.25,
\qquad
\lim_{\tau\to0^+}\frac{E_{\rm child}}{\tau^3}=302.72306084.
\tag{C21}
$$

These coefficients are not finite-time energy predictions or observed
abundances. They are fixed by the same charge and scalar coupling that
determine the membrane's actual-field pressure. In particular,

$$
P(\tau)=P_0+S\chi_W
=\frac14-18.75\tau^2+O(\tau^3),\qquad
\sigma\overline{K}=P(\tau),
\tag{C22}
$$
$$
\lim_{\tau\to0^+}
\frac{P_0-\sigma\overline{K}(\tau)}{F_{\rm in}(\tau)}=1.
\tag{C23}
$$

This is a conditional relation between geometric wall response and
transmitted scalar energy from one supplied state. It does not require
fitting those two outputs independently. It also does not select the
source's history or the absolute vacuum scale.

The energy has an explicit origin. On the affected parent layer, the
leading change in density is S^2(x^2-tau x)/2. Its integral from
x=0 to x=tau is -S^2 tau^3/12, exactly the negative of the child
gain per area. Thus the leading pulse redistributes the initial parent
energy; it does not add an independent source of energy.

Within the initially quiet child, domain of dependence gives the earliest
central front time

$$
\tau_{\rm front}=
\frac{1}{h_i}\operatorname{arsinh}\left(\frac{h_i}{B_i}\right)
=2.07981585635\,L_{\rm ref}/c.
\tag{C24}
$$

This is the leading edge in a causal development, not a detector threshold
or a waveform amplitude. The child remains the initial de Sitter solution
ahead of that front. The parent is not continued as a vacuum spacetime.

**Independent massive-field control.** On a planar fixed interface the
full retarded massive equation is exactly solvable at the wall:

$$
\chi_W=-\frac{S}{2\omega^2}(1-\cos\omega t),\quad
\partial_t\chi_W=-\frac{S}{2\omega}\sin\omega t,\quad
\partial_x\chi_W=-\frac{S}{2}\int_0^t J_0(\omega u)du.
\tag{C25}
$$

At t=.01, omega=20, it gives inward flux .001856325865 and
child energy per area .000006212608434. A separate volume integral
of the retarded field and its derivatives agrees with the integrated
boundary flux; the parent loses the same amount. This checks the source normalization, sign and energy
account. Those finite-time values belong to the planar control, not
to the curved coupled spacetime. A finite wall thickness could also
change the sub-crossing-time transient; C18--C23 use the stated sharp
membrane model.

## What this closes and what follows

The preparation now includes finite scalar energy, its spatial gravity,
the wall's material stress, initial scalar compatibility through second
order, and the resulting initial wall response. The input family is
explicit and the failed frozen-layer assumption has a calculable repair.
The bulk vacuum values remain .125 and .28125; this calculation does
not select a new cosmological Lambda.

These are separate prepared Cauchy data. Their different initial wall
geometries do not describe an instantaneous influence travelling from
the distant packet. Nor does smooth bulk Cauchy theory by itself prove
well-posed evolution of this independent distributional membrane.
All-order corner compatibility, finite-time existence, transmitted
energy, stability, isotropy and clustering remain unestablished.

The next calculation should evolve the uniform excited-parent seed in
the coupled scalar/metric/membrane system and test C20--C23 at finite
times with a controlled remainder. The smooth source remains a separate
initial-data control. A scalar probe, an assumed
child density, or further refinement of the regular pole tail cannot
substitute for that connection. The three likely failures are higher
interface compatibility, constraint drift during evolution, and a
transmitted state whose stress cannot provide the required clustering.

Evidence: [compact source record](../results/summaries/nsc-atlas-parent-source-v1.json).

```
python scripts/run_python.py scripts/check_nsc_atlas_parent_source.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_parent_source.py -q
python scripts/run_python.py scripts/plot_nsc_atlas_parent_source.py
```
