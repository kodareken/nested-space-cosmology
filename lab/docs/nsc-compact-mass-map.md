# The existing compact carrier and four-dimensional Dirac masses

The missing source input is a physical four-dimensional mass, not another
spatial gap. Applying the established conformal Dirac and Kaluza–Klein (KK)
reductions to the declared free carrier resolves its classical mass operator:
on its length-2 interval with the stated chiral endpoint conditions, there is
one massless Weyl mode and a massive Dirac tower. The warp does not shift this
tower in the four-dimensional metric frame defined below. Compact size,
domain selection, interactions and occupations remain physical inputs to
determine before using the [quantum source](nsc-torsion-source-compatibility.md).

This is an application of prior mathematics, not a claim of a new mass
mechanism. Reuse [Fischmann, section3](https://arxiv.org/html/1311.4182v1)
for semi-Riemannian conformal Dirac covariance and
[Grossman–Neubert, section2](https://arxiv.org/html/hep-ph/9912408v3) for
five-dimensional spinors, chiral boundary choices and KK reduction. Their
Randall–Sundrum metric and phenomenological mass assignments are not imported.

## Metric, units and field normalization

The laboratory's section20.39 specifies a conformal factor multiplying the
**entire** five-dimensional spacetime. Restore the chart length L_star and
write Y=L_star y:

\[
g_5=e^{2\sigma(Y)}(g_4-dY^2),\qquad
\sigma(Y)=-\frac{18}{1015}\frac{Y^2}{L_\star^2},\qquad
-L_\star\leq Y\leq L_\star.
\]

Here g4 is the y-independent four-dimensional horizon-penetrating metric.
The cancellation below also holds for any y-independent Lorentzian g4.
The unwarped PG lapse is 1; the full lapse is exp(sigma), and its radial
shift is the same as that of g4. Use signature +----, natural units, one
massless torsionless five-dimensional Dirac field, and no added gauge,
mass or boundary interaction. The fifth tangent gamma matrix is
Gamma^Y=i gamma5, with gamma5 the physical four-dimensional chirality.

Under the canonical identification of the spin frames, the known conformal
law gives

\[
\mathscr D_5=e^{-3\sigma}
  (i\gamma^\mu\nabla_\mu^{(4)}-\gamma^5\partial_Y)e^{2\sigma},
\qquad \chi=e^{2\sigma}\Psi.
\]

The volume factor exp(5 sigma), the two spinor factors and the Dirac factor
cancel in the classical quadratic action:

\[
S_{5,\mathrm{free}}=\int d^4x\sqrt{|g_4|}\,dY\;
\bar\chi(i\gamma^\mu\nabla_\mu^{(4)}-\gamma^5\partial_Y)\chi.
\]

This rescaling also preserves the Cauchy norm: the spatial volume transforms
by exp(4 sigma), while Psi-dagger Psi transforms by exp(-4 sigma). After
the established radial half-density transformation, the Hamiltonian is

\[
H_5=H_{4,0}-i\alpha_Y\partial_Y,\qquad
\alpha_Y=i\beta_D\gamma^5.
\]

H4,0 includes the already derived PG shift and spin connection. Its
massless part commutes with physical chirality. This equality concerns free
evolution in the g4 time frame; it does not assert invariance of a regulated
quantum determinant or its stress under a metric change.

## The declared interval domain yields the usual massive tower

First take the candidate domain with P_R chi=0 at **both** endpoints, where
P_L,R=(1 ∓ gamma5)/2. The endpoint current form vanishes because
P_L alpha_Y P_L=0. This uses full physical chirality; the node/edge labels
of the old spatial discretization do not by themselves derive that domain.
Its nonzero modes are

\[
m_n=\frac{n\pi}{2L_\star},\qquad
f_{L,n}=\frac{\cos[m_n(Y+L_\star)]}{\sqrt{L_\star}},\qquad
f_{R,n}=\frac{\sin[m_n(Y+L_\star)]}{\sqrt{L_\star}},\quad n\geq1.
\]

They satisfy f_L'=-m_n f_R and f_R'=m_n f_L. With
F_n=f_L,n P_L+f_R,n P_R, the exact embedding identity is

\[
H_5 F_n=F_n(H_{4,0}+m_n\beta_D),\qquad
\int_{-L_\star}^{L_\star}F_n^\dagger F_n\,dY=I_4.
\]

Thus each nonzero level contributes one ordinary four-component Dirac field
with action bar(psi_n)(i slash(nabla_4)-m_n)psi_n. This is a scalar mass
coefficient in the four-dimensional action. It is independent of the radial
potential kappa/r and requires no inserted five-dimensional bare mass. On a
curved background this action-level mass is well defined; a global flat-space
momentum pole is not assumed.

The two chiral profiles are parts of one five-dimensional field. Identifying
their compact coupling with the parent–child sheet field Phi is a further
boundary/interaction question, not an equality supplied by the KK reduction.

The zero mode is F0=P_L/sqrt(2 L_star). Its integrated Gram matrix is P_L,
not I4: it supplies a single Weyl field. It cannot be passed to a massive
four-component source calculation as though it contained both chiralities.
The opposite endpoint chirality supplies the mirror zero mode and the same
massive tower. Self-adjointness alone does not choose between them; the two
choices are also explicit in Grossman–Neubert's section2.

For the retained dimensional convention, the mass input is

\[
\frac{m_n}{\Lambda}=\frac{n\pi}{2\sqrt\zeta},\qquad
\zeta=(\Lambda L_\star)^2.
\]

Neither zeta nor a mode's identification as an observed species is fixed by
this relation. In ordinary units the rest energy is n pi hbar c/(2 L_star).
L_star refers to the g4/conformal-coordinate normalization. The proper
compact length is the integral of exp(sigma) dY and is a different quantity.

## What the historical parity label does not determine

The implemented `warped_compact_operator` has coordinate interval [-1,1]
and length 2. The formulas above use that same interval as the fundamental
physical interval. A reflected covering circle would then have circumference
4 L_star. If instead [-L_star,L_star] were the covering circle and Y and -Y
were identified inside it, its quotient interval would have length L_star
and a different tower. The label S1/Z2 alone cannot distinguish those choices.

The even condition belongs to chi, not to the original warped spinor:

\[
\partial_Y\chi_L=0
\quad\Longleftrightarrow\quad
(\partial_Y+2\sigma')\Psi_L=0
\quad\hbox{at the endpoints}.
\]

The first-order domain imposes the chiral projector; this derivative condition
is its compatible mode/squared-operator condition. Imposing ordinary Neumann
conditions directly on Psi would change the domain.

The quadratic warp is smooth on the closed interval, but sigma' is
+36/(1015 L_star) at the left endpoint and -36/(1015 L_star) at the right.
Even reflection across those endpoints has derivative jumps. A smooth global
extension or an interface action has to account for them before treating the
carrier as a self-sourced orbifold. This compact issue is separate from the
already smooth radial throat.

## Why the old spatial scan is not the mass input

The [spatial calculation](nsc-boundary-response.md) uses
Q_y=exp(-sigma/2) partial_y exp(-sigma/2) and a projected warp multiplication.
Its compact singular values and leakage describe the intrinsic spatial
operator. Restoring the lapse and the five-dimensional spinor normalization
gives the different classical operator above. No repeat of that spatial scan
can decide which four-dimensional mass to insert into the backreaction
equations.

## Next physical use and its conditions

The free part of the imported four-dimensional quantum source can now use
the KK mass operator, conditional on this interval realization. The source
still needs the interacting reduction, compact boundary action, scale and
state. A chosen massive level is not generally a closed interacting sector:
mode overlaps can couple it to the Weyl zero mode and the rest of the tower.
The conditional torsion coefficient |xi_eff|=1/3 was a four-dimensional
same-source comparison; carrying it to a reduced five-dimensional tower
requires the interaction and gravitational normalizations to be matched.

There is an existing interaction reduction to reuse too:
[Castillo-Felisola et al., arXiv:1405.0397v1](https://arxiv.org/html/1405.0397v1),
equations13 and18–20, give the higher-dimensional torsion contact and its
four-dimensional axial-vector/axial-tensor channels. Their application keeps
only zero KK modes with stipulated profiles in a Randall–Sundrum metric.
Reuse the contact structure; replace those profiles with the actual NSC
mode embeddings and retain their overlap integrals. Their zero-mode
truncation and removal of the tensor channel do not establish a truncation
for our massive tower. No generic higher-dimensional torsion derivation is
needed to start that comparison.

Keep this transverse compact carrier distinct from the axial P/AP cylinder
used in the vacuum controls and from the unwrapped recursive room chain.
Changing the compact size in spacetime also requires the additional geometric
terms of that varying-size reduction. No static mass formula determines
occupations, a cosmological injection rate, or a stationary geometry.

The executable [convention check](../scripts/check_nsc_compact_mass_map.py)
uses the existing full spinor matrices. It verifies the mass embedding,
current form, warp derivative cancellation, endpoint profiles and mode
normalization exactly for integer n≥1. Run
`python3 scripts/check_nsc_compact_mass_map.py`.
It neither repeats a numerical spectrum nor adds a new numerical record to
the 81-record v0.3.0 collection.
