# A canonical geometry owner for the local source action

The local two-derivative spherical part of the common action now has an
explicit first-order formulation. Its coupling to the calculated charge
mode closes the classical constraint algebra, including the area-dependent
interaction. This identifies an action-specific canonical measure and ghost
operator for that reduced theory. It replaces a guessed scale weight with
the variables and constraints of an existing gravity formulation.

This is a conditional low-energy interface, not a quantization of the full
nonlocal spectral action. The complete coefficients, matter measure, quantum
state, boundary matching and higher-derivative terms remain to be supplied.
No new geometry is fitted or evolved in this calculation.

## Reuse the spherical gravity formulation

Use [Grumiller–Kummer–Vassilevich (GKV)](https://arxiv.org/pdf/hep-th/0204253v9),
sections2.1–2.3 and7, for spherical reduction, first-order gravity and canonical
quantization. These are established methods. The new checks concern the
normalization inherited from the NSC coefficient map, its interacting
source, and the domain at a finite-radius minimum.

The [matched Dirac coefficients](nsc-compact-matching.md) and
[projected gauge interaction](nsc-gauge-source.md) are reused without running
their generators. A4,V4,C_F below denote the full effective coefficients;
the known Dirac contributions are not substituted as their completed values.
Keep A4>0, C_F>0, V4 unrestricted and magnetic flux 2 pi q.

## Normalize the area field and the action

The laboratory's Lorentzian convention is +---, with R_E=-R_L. The local
four-dimensional action is

\[
S_4=\int\sqrt{|g_4|}\,[-A_4 R_L-V_4-C_F F^2].
\]

For ds4^2=ds2^2-r^2 dOmega^2, define the dimensionless area field

\[
X=8\pi A_4r^2>0.
\]

Spherical reduction and the standard integration by parts give

\[
S_2=\int\sqrt{|g_2|}
 \left[-\frac X2 R_2-\frac{(\nabla X)^2}{4X}
       +8\pi A_4-\frac{V_4}{2A_4}X
       -\frac{16\pi^2A_4C_Fq^2}{X}\right].
\]

With GKV's U,V notation this is minus their equation2.9, using

\[
U(X)=-\frac1{2X},\qquad
V(X)=-8\pi A_4+\frac{V_4}{2A_4}X
                 +\frac{16\pi^2A_4C_Fq^2}{X}.
\]

The Ricci orientation agrees with the laboratory convention; the sign
difference is the action convention. The implementation uses GKV's
chapter7 canonical component convention, including the sign explained
in footnote57. A direct component-action check, including its boundary
derivative, fixes the relative sign with the positive-energy matter action.

X is the area already in the metric. The first-order X+,X- and connection
omega are auxiliary representations of this spherical geometry. Their
auxiliary torsion is not the independent five-dimensional skew torsion
in the earlier Einstein–Cartan candidate.

The spherical integration-by-parts term is the boundary integral of
epsilon n.dot(grad X). It cancels the radius-derivative part of the same
four-dimensional GHY term, leaving

\[
S_{\rm GHY,2}=-\int_{\partial M_2}\sqrt{|h_1|}\,\epsilon X K_2.
\]

Here epsilon is the squared unit-normal sign. The canonical-to-bulk
spatial boundary derivative is [X omega_0+X+ e^-_0+X- e^+_0] at the ends.
Corner, null-boundary and transmitting-interface terms require their actual
domain. This local conversion does not supply all those boundary data.

## The same charge interaction closes the constraints

For the full Maxwell convention g4^2=1/(4 C_F), the existing collective
charge scale becomes

\[
\mu^2(X)=\frac{A_4|q|}{2\pi C_F X},\qquad
W_m(X,\chi)=\frac12\mu^2(X)\chi^2.
\]

Its area derivative reproduces the Maxwell angular pressure. The displayed
matter field is the charged collective mode. At |q|>1 the additional neutral
CFT sector must also be included; the preceding numerical source example
uses |q|=1. No new fundamental scalar is added.

Take canonical coordinates q_i=(omega_1,e^-_1,e^+_1) and momenta
p_i=(X,X+,X-), with scalar momentum pi_chi. Write
calV=U X+X-+V+W_m. In the chapter7 convention the constraints are

\[
G_1=X'+X^-e^+_1-X^+e^-_1,
\]
\[
G_2=(X^+)'+\omega_1X^+-e^+_1\mathcal V
       +\frac{(\chi'-\pi_\chi)^2}{4e^-_1},
\]
\[
G_3=(X^-)'-\omega_1X^-+e^-_1\mathcal V
       -\frac{(\chi'+\pi_\chi)^2}{4e^+_1}.
\]

The oriented canonical patch uses e^+_1>0 and e^-_1<0, so constant-x0
slices are spacelike. Other coordinate patches and global gauge coverage
require their own domain check. The canonical x0 is not automatically the
parent Killing time tau used for the earlier outgoing-power calculation.

The gravitational Poisson matrix has entries P12=-X+, P13=X-,
P23=-(U X+X-+V). It satisfies Jacobi. Including W_m in the corresponding
structure functions gives the smeared classical algebra

\[
\{G_i[a],G_j[b]\}=\int ab\,
 \frac{\partial P^{ij}_{\rm matter}}{\partial p_k}G_k.
\]

All three brackets were differentiated as functionals, including their
spatial derivative terms and integration by parts. The residuals vanish
exactly. This matters because the 1/X mass with a constant kinetic coefficient
is not the factorized matter potential printed in GKV equation7.1. Its
area-force term must be retained. An independent Legendre transform of the
covariant scalar action fixes its sign in these same constraints. This is a classical consistency check;
it does not prove the full regulated quantum constraint algebra anomaly-free.

## The measure has an explicit operator owner

The canonical phase-space measure contains Dq Dp Dchi Dpi_chi, constraint
multipliers and BRST ghosts. In the local temporal gauge
(omega_0,e^-_0,e^+_0)=(0,1,0), the ghost operator has the structure

\[
M=\begin{pmatrix}
\partial_0&-1&0\\
0&\partial_0&0\\
\partial_X\mathcal V&X^-U&\partial_0+X^+U
\end{pmatrix}.
\]

This is the established construction of GKV equations7.27–7.32, with the
actual matter-dependent entry checked here. Its formal determinant is
(det partial_0)^2 det(partial_0+X+ U). The record also checks a finite-block
case with noncommuting time and multiplication matrices. A continuum
factorization needs common domains and zero-mode prescriptions; separately
regulated determinants cannot simply be substituted.

The gauge-fixed canonical bulk equations for a fixed charge history include

\[
\partial_0X=X^+,\qquad
\partial_0X^+=-(\partial_0\chi)^2,\qquad
\partial_0X^-=-UX^+X^--V-W_m.
\]

Their full linearized Jacobian in p_i is the same M. This explains the
known geometric ghost/Jacobian cancellation in the applicable canonical
integration. It does not make the matter measure disappear. The standard
covariant scalar configuration measure retains a factor sqrt(e^+_1) in
this gauge; its origin in bosonization and its finite regulator/zero-mode
terms still have to match the common functional. It is not a new prior to
choose for the resolution field.

The coframe equation gives partial_0 e^+_1=U X+ e^+_1, hence
e^+_1/e^+_{1,initial}=sqrt(X_initial/X) on this local gauge patch when the
auxiliary generating sources are zero; the charge history remains. This is
a consequence of the coupled geometry equations. It does
not justify replacing the complete metric/ghost/matter measure by an
independently chosen one-dimensional weight proportional to X^(-1/4).

## Energy feedback and a regular field at the neck

The known primitive w'=exp(Q)V, Q=-log(X)/2, gives the geometric Casimir

\[
\mathcal C=\frac{X^+X^-}{\sqrt X}+w(X).
\]

The new source has the canonical balance

\[
\partial_0\mathcal C
 =-\frac{X^+W_m+X^-(\partial_0\chi)^2}{\sqrt X}.
\]

Thus the geometry is linked to the same charge history and area-dependent
potential. Quantum use requires consistently renormalized source expectations,
the remaining matter measure and a specified state. The classical equation
is not substituted for those missing quantum inputs.

Matching the known vacuum solution fixes G=1/(16 pi A4),
lambda4=V4/(2 A4), rQ^2=C_F q^2/(4 A4) and
calC=-32 pi A4 sqrt(8 pi A4) M_geom. More generally, it agrees with the
charged, vacuum-subtracted spherical mass function. In the asymptotically
flat parent limit with unit-normalized Killing time and the stated reference
subtraction, energy is M_geom/G. A complete outgoing source then obeys
dM_geom/dtau=-G P_infinity with fixed magnetic/vacuum reference terms.
This is not an identification of every gauge time with tau or of a de Sitter
boundary with asymptotically flat infinity. No vacuum metric or collapse
is re-solved here.

At a positive-radius minimum the field transformation r to X is invertible,
even though X'=0. A PG null coframe has determinant1 and remains regular
both at the horizon and at that minimum. In contrast, taking X itself as
a spacetime coordinate fails where its gradient vanishes. Also, X+=X-=0
does not alone imply a Poisson rank drop: the rank is2 if V(X) is nonzero,
and drops to0 only if V(X)=0. Constant-area and sourced branches must be
distinguished rather than discarded as coordinate failures.

## Boundary of this result

The local canonical geometry owner, source coupling and ghost operator are
now explicit. The known curvature-squared coefficients and nonlocal h(D4^2)
remainder are not contained in this first-order truncation. Their smallness
has not been established on the imposed neck. Neither the complete 5D
metric measure nor the common/relative resolution-field measure follows
from this spherical construction. Global gauge coverage, boundaries,
initial state, finite matching terms and physical couplings remain open.

The executable record is `results/development/spherical-action.json`.
`scripts/check_nsc_spherical_action.py --check` authenticates its inputs and
compares every exact identity, formula, scope field and source hash. The
known source, spectrum and horizon generators are not run. The next work
must carry the required matter and nonlocal contributions through this
action/measure interface, rather than repeat the spherical reduction or
choose a convergence weight to make the source stationary.
