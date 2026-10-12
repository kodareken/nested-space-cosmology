# Prepared fixed-C0 tangents fix all twenty incoming normal slots

For the declared smooth compact history class, full prepared fixed-C0
matching at first order requires

\[
\boxed{\delta\partial_T^k a(z)=\delta\partial_T^k r(z)=0,
\qquad k=1,2,3,4.}
\]

Their spatial derivatives therefore vanish, fixing all twenty permitted
normal changes with `k>=1` and `k+j<=4` in
[`IncomingNormalJetChange`](../src/recursive_horizons/nsc_incoming_cauchy_jets.py).
This is necessity for a tangent about the reference, not sufficiency for
matching and not a finite-amplitude exclusion.

## Hypotheses and reused physical result

Reuse the [physical fixed-transfer proof](nsc-incoming-fixed-transfer.md),
including the authenticated physical P16 affine-projector bound, mathematical
SU(2) frames and full coherent source estimate. No physical generator,
source integral, energy scan or background coefficient calculation is rerun.
The new gap is the induction from first-normal functions through the existing
order-four normal-jet domain. Stop after these four orders and twenty slots.

The metric history tangent is real, C-infinity, energy independent, compact
in the axial coordinate and supported on a fixed compact radial slab away
from the horizon. It vanishes in an open upstream neighborhood, retaining
the original preparation and source law. Intrinsic variations of N,beta,a,r
vanish **as functions of z along Sigma**. Normal variations of N and beta
through order four also vanish as functions along Sigma. A restriction at
one point alone would not satisfy these hypotheses.

The smooth history supplies bounded derivatives through the six-term
residual, including forcing derivatives through order six and the necessary
background frame/gap derivatives. This is not a statement about a bare
twenty-number germ. Higher derivatives need not vanish and are not classified.
Use the existing unbounded signed-real-energy continuum and one retained
block with m*ell nonzero, with `Eo=E+omega/2,Ei=E-omega/2` at each fixed omega.

## Six inverse-gap terms and the physical remainder

In the owned mathematical P16 frames, each opposite-band entry has gap
`lambda=+/-2iE/a0²+O(E^-1)` and forcing f whose needed derivatives are O(1).
Set

\[
q_0=-f/\lambda,\qquad q_{j+1}=q_j'/\lambda\quad(j=0,\ldots,4).
\]

The exact telescoping residual of their six-term sum is q5_prime=O(E^-6).
The terms vanish near the upstream endpoint. The reused off-diagonal
generator defect contributes O(E^-17); replacing the physical projectors
in the forcing by P16 contributes O(E^-15). Exact unitary left/right
transport therefore gives a physical pair remainder O(E^-6). The full
thermal/coherent correction is O(E exp(-cE)), including its correlations.
No derivative of a numerical error bound or infinite-series convergence
assumption is used.

## Induction and complex endpoint coefficients

At step k, assume the lower normal a/r functions vanish. Together with the
fixed N/beta jets this gives f^(j)(Sigma)=0 for j<k. Consequently the lower
q terms vanish exactly and `qk(Sigma)=-f^(k)(Sigma)/lambda^(k+1)`.
Gap-derivative terms multiply vanished lower forcing derivatives. In the
fixed axial Fourier coordinate, the normal chain rule is
`partial_rho^k delta_g_hat=(-1/a1)^k delta_g_T^k_hat`; the other chain terms
multiply lower normal derivatives and vanish. The inherited leading Q1
coefficient then gives, with complex Fourier amplitudes A_k and R_k,

\[
\begin{aligned}
\lim E^{k+1}Q_{01}
&=-\frac{a_1^k}{(2i)^{k+1}}
\left[\left(\frac\ell{r_1}-im\right)A_k-
\frac{\ell a_1}{r_1^2}R_k\right],\\
\lim E^{k+1}Q_{10}
&=\frac{(-1)^ka_1^k}{(2i)^{k+1}}
\left[\left(\frac\ell{r_1}+im\right)A_k-
\frac{\ell a_1}{r_1^2}R_k\right].
\end{aligned}
\]

Endpoint frame corrections and higher q terms are O(E^-k-2). The physical
O(E^-6) remainder is smaller than every displayed leading term through k=4.
The two-entry coefficient determinant is

\[
\frac{i m\ell a_1^{2k+1}}{2^{2k+1}r_1^2}\ne0.
\]

Thus both entries vanishing forces A_k=R_k=0 without assuming real Fourier
amplitudes. Matching at **every** fixed transfer and Fourier injectivity
give the normal functions zero pointwise. No uniform large-omega estimate
is needed. Inducting through four orders then permits differentiation in z,
giving4+3+2+1 slots for each of a and r.

## Verification and scope

The [record](../results/development/nsc-incoming-jet-rigidity.json) authenticates
the reused physical bridge and stores exact six-term telescoping, endpoint
quotient/chain identities, all four complex maps and determinants, conjugacy,
and the twenty-slot inventory. Symbolic zeros verify finite identities;
the uniform remainder and injectivity are analytic deductions under the
stated hypotheses, not numerical error bounds or finite-energy thresholds.

The principal risks are a missed lower-jet product term, insufficient history
regularity, and treating complex Fourier amplitudes as real. The explicit
endpoint identities and hypotheses address each. Fixed lapse/shift jets,
intrinsic data and unchanged preparation are substantive restrictions.
The result does not classify varying lapse/shift jets, other preparations,
noncompact histories or finite-amplitude matching. Constraint intersection
is a separate owner; global constraints and extended stationarity remain
OPEN. No new cutoff, action term, state prescription or metric evolution is
introduced.

```sh
python3 scripts/derive_nsc_incoming_jet_rigidity.py --write
python3 scripts/derive_nsc_incoming_jet_rigidity.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_jet_rigidity.py
```
