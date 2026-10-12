# Physical fixed-transfer condition on incoming first-normal functions

For `Eo=E+omega/2`, `Ei=E-omega/2`, and each fixed finite real omega, the
actual prepared covariance response has endpoint limits

\[
\boxed{\begin{aligned}
\lim_{E\to\infty}E^2\delta\widehat C(E_o,E_i)_{01}
&=\frac14\left(\frac{\ell a_1\widehat A}{r_1}
-\frac{\ell a_1^2\widehat R}{r_1^2}-im a_1\widehat A\right),\\
\lim_{E\to\infty}E^2\delta\widehat C(E_o,E_i)_{10}
&=\frac14\left(\frac{\ell a_1\widehat A}{r_1}
-\frac{\ell a_1^2\widehat R}{r_1^2}+im a_1\widehat A\right).
\end{aligned}}
\]

Here Ahat and Rhat are the Fourier transforms, with positive exponent,
of delta a_T(z) and delta r_T(z). They are generally **complex**. If matching
holds at every pair and m*ell is nonzero, subtracting both equations forces
Ahat=0 and then Rhat=0 for every omega. Fourier injectivity on smooth compact
axial tangents implies the necessary pointwise conditions

\[
\boxed{\delta a_T(z)=\delta r_T(z)=0.}
\]

All intrinsic deltaN,delta beta,deltaa,deltar must vanish **as functions of z**
at Sigma. Their pure spatial derivatives then vanish too, and the incoming
normal/half-density restriction has no parameter derivative. A nonzero
intrinsic field, including shift, is a different leading-order case.

## Reuse, decision and stopping condition

Reuse the physical affine-vacuum P16 error and its transport defect, the
[endpoint proof](nsc-incoming-endpoint-obstruction.md), exact midpoint Weyl
vertices and the same retarded source preparation. The gap is the physical
fixed-transfer response of zero-mean first-normal functions. Four local
inverse-gap terms close it without propagation or an energy scan.

The claim uses the existing **unbounded signed-real-energy continuum** and
unchanged coherent horizon/infinity source law. It is a first-order necessary
condition for the stated smooth compact-axial class. It does not impose a
new Hadamard gate, invent a cutoff, classify higher jets, or exclude alternate
preparations, noncompact tangents, finite-amplitude histories or NSC generally.
Stop with this coefficient/physical-remainder connection; root separately
owns any intersection with the constraint equations.

## Exact physical pair equation

Let Po,Pi be the exact affine-vacuum projectors at the two energies. On the
unchanged compact radial slab,

\[
Q'=G_oQ-QG_i+D P_i-P_oD,\qquad Q(\rho_u)=0,
\]
\[
P_oQ+QP_i=Q,\qquad
G_E=\frac{i}{a_0}\left(-m\sigma_1+\frac\ell{r_0}\sigma_2-\frac E{a_0}\sigma_3\right).
\]

For Fourier metric directions, D is `i/a0` times the actual raw Hamiltonian
vertex at the **Weyl midpoint E=(Eo+Ei)/2**. The lapse, scale, radius and shift
vertices are respectively
`-m sigma1+ell sigma2/r-E sigma3/a`, `E sigma3/a²`,
`-ell sigma2/r²` and `+E I`. For real metric functions their Fourier
coefficients may be complex; no Hermiticity of an individual D(omega) is
assumed. The background left/right transports are unitary.

The exact projector constraint is propagated by this equation: for
T=D Pi-Po D,
`Po T+T Pi-T=D(Pi²-Pi)-(Po²-Po)D=0`.
The zero upstream tangent follows from the unchanged past/operator and
source law. The Fourier-pair response is finite; no delta(0), energy
quadrature weight or angular multiplicity enters it.

## Mathematical frames and four inverse-gap terms

Use the existing normalized mathematical Riccati projector P16 at each
energy. With z16=-iaS16, the explicit frame

\[
V=\frac1{\sqrt{1+|z_{16}|^2}}
\begin{pmatrix}1&-\bar z_{16}\\z_{16}&1\end{pmatrix},\qquad
P16=VP_0V^\dagger,\quad P_0=\operatorname{diag}(1,0),
\]

is exactly SU(2). This frame construction does not normalize an archived
column or replace the physical state. The transformed generator is
`L=V†GV−V†V'`. Its off-diagonal part is controlled by the exact identity

\[
V^\dagger(P16'-[G,P16])V=-[L,P_0].
\]

The existing finite Riccati defect is O(E^-16), pointwise with the needed
radial derivatives on this compact slab. Its positive full-collar integral
also gives `Paff−P16=O(E^-16)` uniformly. These are reused background
estimates; no derivative of a numerical error inequality is taken.

Set `Dtilde=Vo†DVi`. Its transformed source
`F=Dtilde P0−P0 Dtilde` is exactly off-diagonal. Since
`D=E D1+D0` with D1 diagonal, V=I+O(E^-1) and its positive-order radial
derivatives O(E^-1) make F and its required radial derivatives O(1). These
frame estimates follow from the finite Riccati coefficient functions on the
compact slab. The opposite-band gaps satisfy

\[
\lambda_{01}=(L_o)_{00}-(L_i)_{11}=-2iE/a_0^2+O(E^{-1}),
\quad
\lambda_{10}=+2iE/a_0^2+O(E^{-1}).
\]

They are uniformly bounded away from zero by a positive multiple of E
for sufficiently large E, at each fixed omega. The gap concerns opposite
bands, not the small difference Eo−Ei.

For either off-diagonal entry f of F, construct

\[
q_0=-f/\lambda,\quad q_1=q_0'/\lambda,\quad
q_2=q_1'/\lambda,\quad q_3=q_2'/\lambda.
\]

The residual telescopes exactly:
`(q0+q1+q2+q3)'-lambda*(q0+q1+q2+q3)-f=q3'`.
The four terms have orders E^-1 through E^-4; smooth compact coefficient
derivatives give residual O(E^-4). All terms vanish near the upstream
endpoint. No convergence of an infinite formal expansion is required.

Let Z contain these off-diagonal sums and zero diagonal. It obeys the
projector-pair constraint for the approximants exactly. Ignoring the
transformed off-diagonal generators contributes only
`O(E^-16)*O(E^-1)=O(E^-17)`. Replacing physical projectors in the source
by P16 contributes `O(E)*O(E^-16)=O(E^-15)`. Returning to the original
frame and applying exact unitary left/right variation of constants gives

\[
\|Q_{\rm physical}-V_o ZV_i^\dagger\|=O(E^{-4})
\]

uniformly on the slab. This supplies the actual prepared-state connection;
a formal reference projector alone would not establish it.

## Endpoint extraction and intrinsic-data assumptions

Put `G=E A+B`, `A=-i sigma3/a²`, and write Q=Q1/E+Q2/E²+... . The checked
leading coefficient is off-diagonal:

\[
(Q_1)_{01}=\tfrac12[(m+i\ell/r)\delta a-i a\ell\delta r/r^2],
\quad
(Q_1)_{10}=\tfrac12[(m-i\ell/r)\delta a+i a\ell\delta r/r^2].
\]

Lapse cancels and shift adds no term to Q1. Intrinsic directions vanish at
Sigma, so Q1=0 there. The next source terms also vanish there, and
`{A,Q1}=0`. The projector constraint fixes the diagonal of Q2 to zero;
the off-diagonal equation is `[A,Q2]=Q1'`. Substituting
`delta a_rho=−Ahat/a1`, `delta r_rho=−Rhat/a1` gives the boxed coefficients.
No fixed-omega, lapse-first-normal or shift-first-normal contribution enters
this E^-2 coefficient under the stated assumptions. This does not say that
those directions have zero response at other orders or transfers.

At Sigma f=0 exactly. The first inverse-gap term therefore vanishes there;
its next term is `−f'/lambda²`, with the same coefficient. The remaining
terms and frame corrections affect E^-3 and beyond. The physical O(E^-4)
error is smaller still.

## Full source and Fourier injectivity

Keep the complete coherent source covariance. By the unchanged source law
and canonical coisometry,
`norm(C(E)-Pv(E))<=max(sqrt(f_H(E)),n_in(E))`, which is exponentially small
for positive large E. The full-minus-vacuum pair equation has source
`D(Ci-Pi)-(Co-Po)D` and zero upstream data. The same left/right unitary
estimate gives O(E exp(-cE)) at fixed transfer, without dropping coherence,
reflection phases or source-complement terms. The exact physical covariance
therefore has the stated leading limits.

For nonzero m,ell the complex two-entry coefficient map has determinant
`i*m*ell*a1³/(8*r1²)`. Both entries, not the real/imaginary parts of one entry
with an assumed real Fourier amplitude, force Ahat=Rhat=0. Reality of the
original functions instead gives Ahat(-omega)=conj(Ahat(omega)) and the
corresponding pair-adjoint relation; the module checks that convention.

Requiring matching for **all** pairs supplies these conditions for every
fixed omega. Fourier injectivity then applies to the smooth compact axial
first-normal functions. No uniform bound as abs(omega) tends to infinity
is needed for this logical implication, and none is claimed. A finite set
of tested transfers would not supply the pointwise conclusion.

## Verification and scope

The module checks frame/projector identities, four-term telescoping, the
pair transport algebra, Weyl midpoint, endpoint equations and the complex
coefficient map. Uniform Big-O estimates and Fourier injectivity are analytic
arguments, with their hypotheses stated above; symbolic zeros are not
numerical error bounds. The authenticated
[record](../results/development/nsc-incoming-fixed-transfer.json) supplies
no finite-energy threshold, propagation, energy scan or arbitrary-data
four-dimensional regularity gate.

This is a necessary first-normal condition in the declared fixed-intrinsic,
unchanged-source, unbounded-energy matching problem. Higher normal jets,
noncompact tangents, finite-amplitude preparation, global constraint/history
solutions, other source prescriptions and full NSC remain outside the result.

```sh
python3 scripts/derive_nsc_incoming_fixed_transfer.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_fixed_transfer.py
```
