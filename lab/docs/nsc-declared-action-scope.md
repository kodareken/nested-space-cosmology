# Declared action and the transmitting history domain

Douglas Ek clarified the action inventory on 2026-09-16. This clarification
supersedes the earlier request for an additional definition of
$\Gamma_{\rm rest}$. It introduces no new term and changes no scientific
record, coefficient, quantum state or boundary condition.

## Count each contribution once

In the interface balance

$$
\pi_p^{ab}+\pi_c^{ab}
+\frac{\delta\Gamma_{\rm rest}}{\delta h_{ab}}=0,
$$

$\Gamma_{\rm rest}$ denotes remaining same-action quantum/field contributions
after excluding the geometry already assigned to $\pi$ and local endpoint
terms already counted. It is a bookkeeping category, not a requirement to
introduce a further independently adjustable functional.

| Contribution | Existing owner and treatment |
|---|---|
| Einstein–GHY and geometric $\pi$ | The [compact boundary action](nsc-compact-boundary-action.md) fixes its coefficient. The [smooth-seam variation](nsc-smooth-seam-variation.md) cancels the two oriented contributions on the shared smooth cut. |
| Weyl, Euler and $\Box R$ endpoint primitives | The [Weyl boundary jets](nsc-weyl-boundary-jets.md) and smooth-seam record own them. They are not added again. |
| Full Gaussian CTP $\Gamma_G=-i\log\det Q$ | [CausalCommonFunctional](../src/recursive_horizons/nsc_causal_common.py) and the [transmitting CTP formulation](nsc-transmitting-ctp-resolvent.md) retain the field, initial correlations and both Schur factors. The complementary factor is not another interface force. |
| Reference subtraction and induced local terms | The existing reference and local-history owners retain their allocation. Their derivatives must be included once in the same metric variation. |
| $93.54264532195464$ | A diagnostic discrete time-node derivative, not a force or counterterm at the physical seam. |

The inspected canonical owner sums its Gaussian matter force and the supplied
same-spectrum induced force. It declares no additional independently weighted
boundary functional. Consequently completion requires the variation of this
declared action,

$$
\frac{\delta\Gamma_{\rm one}^{\rm CTP}[g,C[g]]}{\delta g^A}=0,
$$

with its physical history, reference and state derivatives. An unevaluated or
nonzero bulk variation remains an OPEN or FAIL residual of that equation.
Neither geometric cancellation nor the absence of an additional term sets
the quantum source or the total variation to zero.

## B1d: distinguish metric coordinates from Cauchy time

The [PG-to-KS owner](nsc-pg-ks-metric-pullback.md) fixes the reference chart

$$
dT=-\frac{d\rho}{a_0(\rho)},\qquad
dz=d\tau+\frac{\beta_0(\rho)}{a_0(\rho)^2}\,d\rho.
$$

Thus $T=T(\rho)$, while $z=\tau+S(\rho)$. A homogeneous KS history $g_K(T)$
is not a history $g_K(\tau)$. The transmitting surface $\rho=0$ is one
KS-time slice, with $z$ as its spatial coordinate. Its two oriented sides
are not the first and last time nodes of the old diagnostic history.

The already computed source-mode derivative acts on the full PG field. The
raw KS metric Jacobian and fixed normal variation are also computed. What
the existing `EndpointBranchJets` consumer additionally needs is the chain
rule for its *particular* endpoint coordinates:

$$
\frac{\partial U}{\partial g^B_{K,e}}
=\sum_A\int d\tau\,d\rho\;
\frac{\delta U}{\delta g_P^A(\tau,\rho)}
\frac{\partial g_P^A(\tau,\rho)}{\partial g^B_{K,e}},
$$

with the Cauchy restriction/half-density derivatives if those Cauchy
coordinates vary. The smooth common-jet condition must be imposed on the
variation as well as on the background.

The current endpoint object names the two diagnostic KS time nodes; its
`transmitting_embedding_to_KS_nodes` field is unevaluated. The published
PG/KS owner supplies a compact metric-profile Jacobian, but explicitly does
not identify that profile with these two endpoint directions. Substituting
the profile, identifying $T$ with $\tau$, or treating a same-surface spin
frame change as unitary transport would not implement this chain rule.

An unselected **family** of metric histories can be differentiated; a
stationary solution need not already exist. The missing item is the
definition of the endpoint variation family on the transmitting domain,
not a chosen physical duration and not a new action term. Until that family
is specified or the full history variation is evaluated directly in its
own domain, this eight-coordinate endpoint consumer remains OPEN.

## Current execution status

- The author's action inventory is adopted; no additional
  $\Gamma_{\rm rest}$ is requested or introduced.
- B1d is OPEN on the physical history/endpoint chain rule above. This
  clarification supplies no numerical transmitting endpoint jets.
- B2, $V_c$ and extended stationarity remain unevaluated; there is no new
  existence or non-existence certificate for the extended class.
- Previously computed geometric and Gaussian residuals retain their
  recorded scope. No scientific generator was rerun for this clarification.
- $A,q,\Omega,\zeta,V_{\rm full}$, seeds and PDF v0.26.0 are unchanged.
  No finite stress, selected history or metric timestep is produced.
