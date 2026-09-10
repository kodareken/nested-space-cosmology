# From the common operator to its source equations

This derivation joins the existing NSC Hamiltonian, quantum state, boundary
response and metric variation. Its useful connection is explicit: one complex
parent–child coherence determines link energy and transfer, and the variation
of that same Hamiltonian supplies its geometric force. The source is not a
second fluid assigned to the geometry.

The calculation reuses the [canonical influence](nsc-influence.md),
[mass/charge identification](nsc-observable-bridge.md),
[covariant source](nsc-covariant-source.md) and
[compact matching](nsc-compact-matching.md). No spectrum, angular tower,
thermal comparison or collapse solution is recalculated. The finite
quadratic fermion identities below are exact. Their absolute continuum
normalization must be supplied by the same matched quantum functional.

## 1. One state and one operator

Work in natural units, hbar=c=1, and a common canonical spin frame. Let X
collect the lapse, shift, spatial metric, sphere radius, gauge connection
and existing link data. With the domain pulled back to a fixed Hilbert
space, the parent–child Hamiltonian and covariance are

\[
H[X]=\begin{pmatrix}H_p&B\\B^\dagger&H_c\end{pmatrix},\qquad
C_{ij}=\langle c_j^\dagger c_i\rangle,\qquad
0\le C\le I,\qquad i\dot C=[H,C].                         \tag{1}
\]

The covariance includes C_pc and C_cp: eliminating a room does not erase
the correlations in the prepared state. Changing the spin frame sends
B to U_p B U_c^dagger and C_cp to U_c C_cp U_p^dagger. Thus

\[
\boxed{z_B=\operatorname{Tr}_p(B C_{cp})}                 \tag{2}
\]

is invariant under those frame changes. It is the relevant complex
contraction, rather than a sign assigned to a bare off-diagonal matrix entry.

The wave/particle reduction already established in the project is retained:
for a constant scalar link in its invariant chiral sector,

\[
H_{\rm eff}=\boldsymbol\alpha\cdot\mathbf p+\Phi\beta,
\qquad E^2=|\mathbf p|^2+|\Phi|^2.                       \tag{3}
\]

That relation is conditional on the scalar Clifford structure and invariant
domain; it does not identify every geometric link B with a mass. For a
mode written a(x)u(x) exp[-i S(x)], the leading eikonal equation gives
(gamma^mu pi_mu-m)u=0, with pi_mu=partial_mu S-e A_mu and pi^2=m^2.
This connects the phase gradient, frequency and mass shell in the scalar
sector. A real sine is a projection of that complex phase. Opposite charge
orientation is represented by the charge-conjugate spinor and gauge
representation; the two signs of a sine alone differ by a phase shift.

## 2. The real and imaginary parts have definite physical meanings

The interaction energy and regional occupation numbers are

\[
E_{\rm link}=2\Re z_B,\qquad
N_p=\operatorname{Tr}C_{pp},\qquad N_c=\operatorname{Tr}C_{cc}.
\]

Taking the pp block of (1) gives

\[
\dot N_p=-i\operatorname{Tr}_p(B C_{cp}-C_{pc}B^\dagger)
=2\Im z_B=-\dot N_c.                                  \tag{4}
\]

For equal charge e on the two sectors, the corresponding charge current
into the parent is e dot N_p. More generally, use the actual conserved
charge representation, which must intertwine with B.

A phase coordinate on the existing link, B(theta)=exp(i theta)B, obeys,
at fixed instantaneous C,

\[
\boxed{E_{\rm link}=2\Re z_B,\qquad
\dot N_p=2\Im z_B=-\left.\partial_\theta E_{\rm link}
\right|_C.}                                            \tag{5}
\]

This does not introduce a scalar field or a freely chosen current. It
relates two components of the same coherence. A simultaneous gauge/frame
change of the link and state leaves z_B unchanged; such a change is not a
physical phase deformation. A variation of the link amplitude likewise
differentiates the same E_link, rather than adding a new source energy.

Energy transfer carries Hamiltonian weights. Define

\[
\begin{aligned}
E_p&=\operatorname{Tr}(C_{pp}H_p),&
J_p&=2\Im\operatorname{Tr}_p(H_p B C_{cp}),\\
E_c&=\operatorname{Tr}(C_{cc}H_c),&
J_c&=-2\Im\operatorname{Tr}_p(B H_c C_{cp}).
\end{aligned}
\]

Equation (1) and cyclicity of the trace give the complete account:

\[
\begin{aligned}
\dot E_p&=\operatorname{Tr}(C_{pp}\dot H_p)+J_p,\\
\dot E_c&=\operatorname{Tr}(C_{cc}\dot H_c)+J_c,\\
\dot E_{\rm link}&=2\Re\operatorname{Tr}_p(\dot B C_{cp})-J_p-J_c,\\
\dot E_{\rm tot}&=\operatorname{Tr}(C\dot H),\qquad
E_{\rm tot}=E_p+E_c+E_{\rm link}.                       \tag{6}
\end{aligned}
\]

The last identity also follows immediately from
Tr([H,C]H)=0. For a time-independent complete Hamiltonian, total energy is
conserved. For evolving geometry, the right side is its work on the field.
The link stores energy, so J_p and J_c need not cancel separately. External
ports or additional rooms must be included when applying this finite
account to a stationary open region. Cosmological density supply further
requires the proper clock, volume and deposition map already identified
in the [observable bridge](nsc-observable-bridge.md).

## 3. The same interaction supplies geometric stress

For prescribed histories X_plus and X_minus, let U_plus and U_minus be
generated by (1), with one prepared covariance C_0. The existing
normalized Gaussian functional is

\[
Z_F[X_+,X_-;C_0]=\det(I-C_0+C_0U_-^\dagger U_+),\qquad
\Gamma_F=-i\log Z_F.                                  \tag{7}
\]

This uses [Klich's Fock-space trace identity](https://arxiv.org/abs/cond-mat/0209642).
Use its existing project implementation; there is no new determinant
algorithm here. The determinant form applies to the quadratic fermion
sector. With genuine many-body interactions, the corresponding Fock-space
trace remains the definition, but a one-particle covariance alone no longer
determines the state or its evolution.

On the physical history, the first variation is already known:

\[
\left.\frac{\delta\Gamma_F}{\delta X_\Delta^A(t)}\right|_{X_+=X_-}
=-\operatorname{Tr}\left(C(t)\frac{\delta H}{\delta X^A(t)}\right).
                                                               \tag{8}
\]

Here X_Delta=X_plus-X_minus and C_0 is held fixed in the canonical frame.
If the preparation or domain itself is varied, its initial-surface
contribution must be included. Re-preparing a different vacuum for every
metric is a different variation.

For a metric variation, the link part of (8) contains

\[
-2\Re\operatorname{Tr}_p\left(\frac{\delta B}{\delta g^{\mu\nu}}C_{cp}\right).
                                                               \tag{9}
\]

Thus link stress, link energy and link transfer depend on the same B and
C, with the prescribed measure, frame and domain variations retained.
For the project's +--- convention,

\[
T^F_{\mu\nu}=\frac{2}{\sqrt{-g}}
\left.\frac{\delta\Gamma_F}{\delta g_\Delta^{\mu\nu}}\right|_{\Delta=0}.
                                                               \tag{10}
\]

This is the finite source identity. In the continuum it denotes the
covariantly normalized expectation value, including its assigned local
terms, rather than an unrenormalized coincident trace. The existing static
source formulas for lapse, radial metric and sphere radius are its
coordinate variations. The influence-functional organization of metric
backreaction is established prior work, for example
[Martín–Verdaguer](https://arxiv.org/abs/gr-qc/9904021); its application here
retains the project's Dirac operator and state.

## 4. Eliminating the child leaves its metric variation in the action

Write the complete contour kernel, including preparation data, as

\[
\mathscr K=\begin{pmatrix}K_p&-V\\-W&K_c\end{pmatrix},\qquad
G_c=K_c^{-1},\qquad S=K_p-VG_cW.                       \tag{11}
\]

In the differential bulk V=B and W=B^dagger; initial cross-correlations
can also contribute at the contour's preparation surface. In a finite
representation with the required inverses, block elimination gives

\[
\log\det\mathscr K=\log\det K_c+\log\det S,
\]
\[
\delta S=\delta K_p-\delta V G_cW-VG_c\delta W
                 +VG_c\delta K_cG_cW,                 \tag{12}
\]
\[
\boxed{\delta\log\det\mathscr K
=\operatorname{Tr}(G_c\delta K_c)
 +\operatorname{Tr}(S^{-1}\delta S).}                 \tag{13}
\]

Use a consistent local logarithm branch and the same normalization on both
sides. The child's own determinant occurs once. All link variations are
already inside delta S; adding another copy would double-count them.
Continuum domain and determinant issues are not replaced by this finite
identity. The normalized first-order recursive equation in the roadmap is
the scale-converted version of (11), with its stated 1/Omega and ordered
matrix products retained.

The [explicit boundary-state construction](nsc-boundary-state.md) now
supplies the finite quadratic state's two-time lesser/greater kernels,
mixed initial-surface terms and normalized endpoint action. It removes
that finite state kernel from the list of placeholders. Its room partition
does not by itself perform the compact/light vacuum matching in (16).

## 5. The coupled equation and the remaining specified input

Let S_b,Theta denote the real branch contribution actually obtained by
matching the remaining terms of the same action, with the retained
canonical determinant counted once. When such a branch representation is
available, its finite Gaussian assembly is

\[
\Gamma_\Theta[X_+,X_-;C_0]
=S_{b,\Theta}[X_+]-S_{b,\Theta}[X_-]+\Gamma_F[X_+,X_-;C_0].
                                                               \tag{14}
\]

Equations (1) and (8) then give the self-consistent construction:

\[
\boxed{
i\dot C=[H[X],C],\qquad
\frac{\delta S_{b,\Theta}}{\delta X^A}
=\operatorname{Tr}\left(C\frac{\delta H[X]}{\delta X^A}\right).
}                                                              \tag{15}
\]

In particular, defining E^b_munu=-2(delta S_b/delta g^munu)/sqrt(-g)
gives E^b_munu=T^F_munu. Metric, link and relative-scale equations are
different variations of (14); no independent dark-fluid function is added.
If sectors have been integrated out with memory and fluctuations, their
full mixed contour functional must replace the simple branch-action term
in (14). A real branch action cannot stand in for an unspecified quantum
influence kernel.

The stored Euclidean allocation is

\[
\Gamma_5^E=\Gamma_{\rm light,ren,\mu}^E
             +\Gamma_{H,\nu}^E+\mathcal C_{\nu,\mu}^E.          \tag{16}
\]

It supplies definite spectral coefficients and a Euclidean remainder.
The [canonical/spectral bridge](nsc-canonical-spectral-bridge.md) now
resolves that remainder into its canonical massive-field contribution,
light endpoint conversion and warp/covariant conversion. The real-time
conversion to the absolute canonical source remains the next owner.
The [explicit Gaussian measure](nsc-gaussian-cutoff-measure.md) supplies
its finite Euclidean operator representation. At fixed cutoff and
normalization, its ordered warp Jacobian has the variation

\[
\delta S_{\rm measure}
=\operatorname{Tr}\delta\sigma
 +\tfrac12\operatorname{Tr}[(I-e^{-L_s/\Lambda^2})L_s^{-1}\delta L_s],
\qquad L_s=D_s^2.
\]

The apparent inverse is removable at L_s=0, with limit 1/Lambda².
For the pointwise Weyl tangent delta D_s=-{delta sigma,D_s}/2 at fixed
canonical D0, cyclicity combines the two terms into
delta S_measure=Tr[delta sigma exp(-L_s/Lambda²)]. This recovers the
retained regulated warp source before taking a continuum limit, without
an additional pressure coefficient. The finite measure is now explicit;
its physical time-contour and preparation prescription is still required.
The [vacuum-matched CTP development](nsc-vacuum-matched-ctp.md) now states
one explicit prescription and evaluates its real geometric vertex through
quadratic order in the radius probe. Its finite state factor is canonical;
the complete constrained metric source and continuum continuation remain
to be constructed. This is not an assertion of uniqueness or a physical
stability result.
It has not supplied the complete physical-contour version of the last
two terms, its measure/phase contributions, and their common state.
Consequently (14) is a derived assembly rule, not an evaluated S_b,Theta
for the full NSC action. Equation (15) makes the exact missing input visible:
the matched vacuum/boundary functional whose variations enter its left side.
The [state-regulator result](nsc-state-regulator.md) already establishes why
one cannot simply substitute the raw finite heat determinant for (7).
That comparison is not repeated here.

There are two distinct constructions: retaining (16) as a full functional
requires its actual contour completion; retaining its local spectral
coefficients as matching data for a canonical low-energy action defines
an effective realization with a stated truncation. The latter does not
establish equality with the entire finite heat determinant. It also needs
control of the omitted terms before application at the stored neck.

For the existing development geometry the [computed parent state](nsc-unruh-state.md)
has (rho,p_parallel,p_sphere,T_Tz) approximately
(-0.00330689624,-0.04949459659,-0.04299577129,-0.000003048593).
These are reused data. Their negative null contraction is a useful source
contribution, while the complete componentwise equation (15) remains to
be solved. Defining the missing action derivative to equal a desired
residual would choose the answer. Its derivation must instead come from
the common operator, normalization and state, and the geometry can change
when its source equations are solved.

The [ADM source construction](nsc-adm-source-constraints.md) now makes
the independent lapse, shift, radial and sphere variations explicit.
The stored parent-state tensor has also been
[transported to the PG normal frame](nsc-adm-neck-source-map.md), giving
the canonical contribution to all four equations while preserving the
same Killing power. Its Einstein coefficient remains symbolic; the
remaining forces must come from the common action rather than from
selecting a componentwise fit.

## What this construction establishes

The [parent flux condition](nsc-parent-backreaction-gate.md) now selects
evolution for the current empty-incoming free-source branch. An analytic
geometric remainder with the stated stationary symmetry cannot cancel the
recorded outward power. A different directed state/domain contribution
must be explicit; no global time reflection of the Lorentzian horizon
patch is assumed. The corresponding initial mass balance is conditional
on the matched asymptotic Einstein normalization.

The [matched local warp contribution](nsc-warp-local-neck-source.md)
now supplies a nonzero part of the remaining force vector on the actual
geometry. Its vacuum and Weyl-squared sources are counted on the source
side; its Einstein increment belongs in the full geometric coefficient.
The uncomputed nonlocal and other-sector contributions remain in the
same componentwise matching equation.

The exact binding is (2), (5), (6), (8) and (13): common coherence, transfer,
work, geometric force and child elimination form one consistent finite
quadratic system. These use established quantum mechanics and operator
identities; they are not claimed as newly discovered laws. Equation (15)
is the corresponding self-sourcing equation once its matched functional
is specified. A full NSC solution additionally requires that functional,
its state/domain and a solution of its independent equations with the
same Theta. The recursive identity remains the organizing condition.

This derivation adds no simulation or stress-test campaign. Its next
construction task is the missing physical-contour contribution in (16),
including its link and metric dependence, rather than another refinement
of the completed canonical source calculation.
