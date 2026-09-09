# The compact chiral domain supplies the Einstein boundary term

The declared compact domain determines the leading geometric boundary terms
of the Euclidean determinant magnitude. The Einstein boundary coefficient
has the required ratio two to the bulk coefficient; no independent boundary
weight is needed at this order. This replaces the conditional GHY comparison
in the [vacuum matching record](nsc-vacuum-charge-matching.md) for this
specific reflecting interval realization.

The boundary is the compact endpoint, not the radial parent–child throat.
The Lorentzian state, determinant phase, transmitting interaction and full
finite action are still separate parts of the common-functional problem.

## Match the adjoint before importing a heat formula

Reuse the [conformal compact reduction](nsc-compact-mass-map.md) and its
paired charged projector. In a smooth Euclidean development collar,

\[
g_5=e^{2\sigma(Y)}(g_4+dY^2),\quad
P=\tfrac12(1-\eta\tau_3\otimes\gamma^5),\quad \eta=\pm1,
\qquad (1-P)\Psi|_{\partial M}=0.
\]

Here P describes allowed boundary values. The two choices remain admissible.
The Cauchy-current matrix is I2 tensor i beta gamma5, whereas the Euclidean
normal gamma is plus or minus I2 tensor gamma5. Consequently

\[
P\alpha_Y P=0,\qquad P\Gamma_n P\ne0,\qquad
(1-P)\Gamma_n P=0.
\]

The Lorentzian Hamiltonian domain is isotropic for its current. It does not
make the Euclidean D_E self-adjoint in the measure sqrt(g5) Psi-dagger Psi.
The latter's adjoint domain has complementary allowed boundary values 1-P.
This distinction is needed when applying the boundary Green-form criterion
in [Vassilevich, equation3.34](https://arxiv.org/html/hep-th/0306138v3).

The positive operator for the determinant magnitude is D_E-dagger D_E.
Its domain requires both (1-P)Psi=0 and P D_E Psi=0. On the conformal collar,
the projected tangential term vanishes. Using the existing conformal Dirac
law, the remaining condition is

\[
(\partial_Y+2\sigma')P\Psi=0,
\qquad
(\nabla_{n_{\rm in}}+S)P\Psi=0,
\qquad S=-\tfrac12K P.
\]

K=4 n_out^Y exp(-sigma) sigma' is the outward mean curvature; it equals
Vassilevich's L_aa, defined with the inward normal in the second fundamental
form. The sign works at both endpoints. The local Dirichlet/Robin problem
is strongly elliptic and its nonnegative quadratic form is ||D_E Psi||².
Zero modes must be treated separately when defining the determinant.
Using this modulus does not determine its phase or prove a continuation
through a trapped region.

## Apply the existing mixed heat coefficients

Use the universal formulas in
[Vassilevich, equations5.30–5.32](https://arxiv.org/html/hep-th/0306138v3).
For unit smearing, the torsionless geometric endomorphism is E=-R/4.
The rank is r=8 for two bulk Dirac copies. With chi=2P-1, equal projector
ranks and tr(chi_:a chi_:a)=r K_ab K_ab give

\[
a_1=0,\qquad
a_2=-\frac{r}{12(4\pi)^{5/2}}
 \left(\int_M\sqrt g R+2\int_{\partial M}\sqrt h K\right),
\]
\[
a_3=\frac{r}{128(4\pi)^2}
 \int_{\partial M}\sqrt h\,(K^2-2K_{ab}K^{ab}).
\]

Thus this domain produces no order-Lambda4 boundary-volume term. The a2
coefficient yields the usual Einstein/GHY ratio. On the actual conformal
collar K_ab=(K/4)h_ab, the displayed a3 density is nonnegative. These are
applications of known heat formulas, not new universal heat coefficients.
The runner checks the projector and curvature contractions directly with
the project's spinor matrices, for both complementary orientations.

## Consequence for the source coefficient

The same proper-time window weights a2 by (Lambda3-nu3)/3 and a3 by
(Lambda2-nu2)/2. The induced boundary term cancels the endpoint part of the
previously computed bulk Einstein warp contribution:

\[
V_{\rm warp,bulk}=A_5(8B-12J),\quad
V_{\rm boundary,a2}=-8A_5B,\quad
V_{\rm warp,total}=-12A_5J.
\]

A5 includes the same fermion copies; B and J retain the definitions in the
authenticated vacuum-matching record. Its already computed Gaussian
integrals are reused. The remaining negative correction is bounded by
0.000697052/zeta relative to the retained positive volume term in the
proper-time window. The displayed a3 contribution on this static collar is

\[
V_{\partial,a3}=\frac{N(\Lambda^2-\nu^2)}{8(4\pi)^2}
\sum_{Y=\pm L_\star} e^{2\sigma}(\sigma')^2\ge0,\qquad N=2.
\]

For these retained channels and zeta>=1 the earlier charged-seed combination
therefore obeys Xi>=2.9979088 q². The necessary charged-extremal range is
still Xi<=1/4. A controlled heat approximation additionally requires
nu_match L_star much larger than one and small curvature/gradient scales.
Higher curvature, the actual link and gauge backgrounds, finite boundary
terms and the remaining quantum contributions have not been evaluated here.
This is not a bound on the complete action.

## The remaining metric equation is an interface balance

A derived GHY coefficient does not hold an unfixed boundary in place. Reuse
the metric variation in
[Krishnan–Raju, equations4–5](https://arxiv.org/html/1605.01603v1):
pi^ab=-A5 sqrt(|h|) epsilon (K^ab-K h^ab). Here epsilon is the squared
outward-normal sign. If the boundary metric is varied freely and no other
terms are retained, its stationarity condition is K_ab-K h_ab=0.

For the declared Gaussian warp, each endpoint instead has

\[
K_{ab}-Kh_{ab}=\frac{6a e^a}{L_\star}h_{ab},
\qquad a=\frac{18}{1015}.
\]

This nonzero residual is checked at both ends. It is an equation for the
remaining source, not permission to insert a brane tension. For a common
metric on a joined interface, the action requires

\[
\boxed{\pi_p^{ab}+\pi_c^{ab}
      +\frac{\delta\Gamma_{\rm rest}}{\delta h_{ab}}=0.}
\]

Gamma_rest includes the remaining boundary/quantum/field contributions
and excludes the geometric terms already used in pi. Shared spin/metric
frames, domains and states must be matched before evaluating it. Holding
the endpoint metric fixed is a valid development boundary condition; it
does not establish the self-sourced interface. The higher-order and
transmitting metric variations have not been solved here.

## Decision and remaining owner

The positive bulk result made another cutoff-profile or throat scan
uninformative. Three different completions were considered: the complete
quantum/constraint contribution, a spectral volume constraint, and the
boundary terms fixed by the existing domain. The boundary map was chosen
because it supplies an actual missing coefficient with no additional field
or weight. The full quantum/state problem remains necessary.

Spectral volume quantization is a different, existing construction.
[Chamseddine–Connes–Mukhanov, equations4–16](https://arxiv.org/html/1409.2471v2)
adds coordinate-map relations and constrains compact Euclidean four-volume.
It removes the fixed-volume contribution from local metric variations but
leaves cosmological curvature as an integration constant. Those extra
relations are not consequences established by the NSC recursion, and do not
select the charged-source coefficient here. No volume constraint, mimetic
clock or vacuum subtraction has been adopted.

Next evaluate this interface balance with the determinant phase, the
field/link-dependent terms and common quantum measure in the actual
transmission domain. The compact
geometric a1–a3 terms above are now reusable inputs; they do not require a
second derivation or a new numerical throat solution.

Reproduction:

    python scripts/check_nsc_compact_boundary_action.py --check

The [record](../results/development/compact-boundary-action.json) compares
every field and authenticates its dependency closure. Its float tolerance is
3e-13 absolute and relative, with no exceptions. It runs only the new domain
and coefficient checks. No prior generator or old warp integration is run.
