# The paired compact Dirac operator and its universal boundary anomaly

The universal parity-even conformal-anomaly contribution is now evaluated
for the existing free conformal compact carrier. It vanishes on this
geometry, as does the first metric variation of its universal Wess–Zumino
representative. The previously computed finite-cutoff warp contribution
remains nonzero. These are different parts of the quantum source.

This application reuses the [compact adjoint/domain construction](nsc-compact-boundary-action.md)
and [canonical/spectral allocation](nsc-canonical-spectral-bridge.md).
The new connection is an exact map to the published doubled-Dirac boundary
problem, followed by substitution of the actual compact geometry into its
known invariants. No heat coefficient, spectrum or source integral is
recomputed. The [record](../results/development/compact-anomaly-bridge.json)
and [derivation script](../scripts/derive_nsc_compact_anomaly.py) retain
the exact import algebra and its assumptions.

## Match the closed operator, not just the projector's appearance

For the two identical torsionless copies, the original Euclidean
differential expression is D=I2 tensor D0, with allowed boundary values

\[
P=\tfrac12(I-\eta\tau_3\otimes\gamma^5),\qquad\eta=\pm1.
\]

Its adjoint has the complementary allowed domain. Thus the original
positive determinant operator is L=D^dagger D, not D squared on P.
Let U=tau1 tensor I4 and define the new representative

\[
\widetilde D=UD,\qquad\widetilde\Gamma^A=U\Gamma^A.
\]

U is a parallel, constant unitary; UPU=I-P and it commutes with the
identical-copy bulk differential expression. Consequently

\[
\operatorname{Dom}\widetilde D^\dagger
=U\operatorname{Dom}D^\dagger=\operatorname{Dom}D,
\qquad \widetilde D^\dagger=\widetilde D,
\qquad\boxed{\widetilde D^2=D^\dagger D.}
\]

For the last equality the domain is also equal: UD psi must have allowed
values P exactly when D psi has the adjoint's allowed values I-P.
This is left multiplication, not a similarity transformation. It preserves
the modulus and its heat coefficients. Left unitarity alone does not fix
the original determinant phase.

Put chi=2P-I and Gamma_star=tau2 tensor I4. At outward normal sign
epsilon=±1 the published boundary form is recovered by

\[
\{\chi,\widetilde\Gamma_n\}=0,\qquad
[\chi,\widetilde\Gamma_a]=0,\qquad
\chi=(-\eta\epsilon)i\Gamma_\star\widetilde\Gamma_n.
\]

The opposite component signs preserve the original coordinate-chirality
choice at both endpoints. They do not change the parity-even conformal
charges. The compatible Robin coefficient remains -K/2. Both endpoint
orientations and both allowed eta values are retained in the exact
matrix calculation. If the Lorentzian time gamma is transformed together
with D, its factors cancel when extracting the canonical Hamiltonian.
Different copy couplings or domains require their own map.

## Apply the published coefficients

[Faraji Astaneh–Solodukhin, PRD 108, 085015 (2023)](https://scoap3-prod-backend.s3.cern.ch/media/files/81024/10.1103/PhysRevD.108.085015.pdf),
Eq.80 and section VI.A.3, give the doubled-field gravitational anomaly
in the normalization 1/[5760(4 pi)^2], with

\[
(a,c_1,\ldots,c_8)
=(0,-429/4,621,0,0,270,990,-210,-360).
\]

Use their invariant definitions, Eqs.81–89. On the project geometry

\[
g_5=e^{2\sigma(Y)}(g_4+dY^2),\qquad\partial_Yg_4=0,
\]

the endpoints satisfy K_ab=(K/4)h_ab, so the trace-free shape K_hat
vanishes. The product geometry has W_nabc=0, which remains zero under
the conformal transformation. Therefore I1, I2, I5, I6, I7 and I8 vanish.
The remaining intrinsic-curvature invariants E4, I3 and I4 can be nonzero,
but their coefficients are zero. No conformal-flatness assumption on g4
is used.

The gauge term in the same paper, Eq.96, is proportional to F_an F^an.
It vanishes for F_aY=0, including the retained Y-independent tangential
gauge backgrounds. Thus

\[
\boxed{\mathcal A_{\partial,\mathrm{universal}}=0}
\]

for this paired free sector and geometry. This does not assert cancellation
of all chiral, gauge or global anomalies of an interacting realization.

## The universal anomaly action and its source

The nontrivial anomaly coefficients above determine the universal
Wess–Zumino representative. Since a=0, its gravitational terms are the
boundary conformal invariants weighted by the Weyl parameter. Sigma(Y)
is constant along each boundary component.

At the unwarped product, K_hat=0, W_nabc=0 and the normal derivative
of W_anbn vanishes. The nonzero-coefficient invariants start at quadratic
or higher order in deviations: I1/I2 are quartic in K_hat; I5/I6 are
quadratic in K_hat; I7 is quadratic in W_nabc; I8 contains products of
vanishing shape/normal-variation data and their tangential derivatives.
The normal gauge term is also quadratic about F_an=0. Their first
variations therefore vanish, including variations away from the product
ansatz. Conformal invariance transports this statement to the warped
background; constant boundary sigma introduces no tangential weight terms.

Consequently the universal representative has

\[
\mathcal J_{\mathrm{WZ,univ}}=0,\qquad
\delta_g\mathcal J_{\mathrm{WZ,univ}}=0
\]

on this class. Its corresponding local real-time branch difference is
also zero. This fixes this part of the source; it does not remove the
conformally invariant integration functional, scheme-dependent trivial
anomalies, state dependence or finite-cutoff terms. No finite coefficient
has been set to zero by a renormalization choice here. A zero trace anomaly
also does not mean that vacuum energy or the entire stress tensor vanishes.

## A separate statement about the paired phase

The same representative has a domain-preserving spectral grading
G=tau3 tensor I4:

\[
[G,P]=0,\qquad\{G,\widetilde D\}=0.
\]

Its nonzero eigenvalues pair as +lambda and -lambda. Therefore its eta
function vanishes where the spectral sum is defined, and so does its
analytic continuation. For a discrete spectrum with no zero modes,
closed tangential sections and no additional boundary terms, the present
integrated a5=0 also gives zeta_L(0)=0.

For the spectral cut with arg(-lambda)=pi, pairing gives

\[
\zeta_{\widetilde D}(s)
=\tfrac12(1+e^{-i\pi s})\zeta_L(s/2),\qquad
\det_\zeta\widetilde D
=(\det_\zeta L)^{1/2}e^{i\pi\zeta_L(0)/2}.
\]

The representative's phase is then fixed to one in this convention.
Transporting a paired convention through U is a specified phase choice,
not a consequence of left unitarity alone. This does not settle zero-mode
occupations, noncompact relative eta invariants, changing topology or the
full Lorentzian transition determinant. No phase from those problems is
silently discarded.

## Preserve the finite-cutoff result

The [canonical/spectral bridge](nsc-canonical-spectral-bridge.md) already
records nonzero warp-path derivatives of Q1 and Q2, approximately -0.1037363
and -0.407662 at its development point. They are read unchanged. The raw
finite-endpoint conversion J_Lambda is not the universal logarithmic
anomaly: it includes the explicit cutoff and covariant-measure response.
It must not be set to zero using the result above.

The remaining real-time source is therefore the finite conversion C5+J,
with the physical cutoff/state/measure prescription specified, plus any
interactions and genuine transmission contributions of the same action.
An unspecified universal anomaly coefficient can no longer serve as its
owner in this free conformal compact sector. The full self-sourcing
equations remain unsolved.

## Evidence and reuse

```sh
python -B scripts/derive_nsc_compact_anomaly.py --check
```

The command authenticates stored inputs and reproduces the new domain
algebra and invariant substitution across every record field. The
closed-domain and first-variation arguments are stated above; a matrix
commutator alone is not substituted for those arguments. The published
charges are imported facts, not a new heat-kernel calculation or a claim
of original universal anomaly coefficients.
