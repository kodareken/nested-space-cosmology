# Adopted source-fixed incoming history problem

Douglas Ek explicitly adopted this formulation on 2026-09-18: "same quantum
state" means the same upstream/source preparation, not a covariance pinned
to the incoming surface independently of the metric. The definition question
is resolved. This document owns that decision and supersedes the earlier
pending-choice wording.

The [finite preparation theorem](nsc-evolved-incoming-state.md)
excludes simultaneous included incoming constraint closure in the declared
smooth compact history class with exact incoming C0 fixed. More amplitudes,
normal jets or finite-pair compensators in that same class cannot remove
the inherited shift residual. The larger transmitting theory remains open.

The fixed-C0 theorem is retained as a regression certificate. It does not
constrain the new search to preserve C0, and compensators or extra jets to
restore that frozen covariance are no longer the active task.

## Source-fixed formulation

Keep the authenticated upstream covariance `C_up`, horizon/infinity source
occupations, phases, channel inventory, couplings and intrinsic incoming
identification. For each allowed metric history define

\[
C_\Sigma[g]=U_g(T_\Sigma,T_u)\,C_{\rm up}\,
U_g(T_\Sigma,T_u)^\dagger,
\]

in the same canonical KS variables. Equivalently use the existing full
source columns and fixed restriction map. This formula determines the
state; it does not fit or prescribe an independent incoming stress.

The incoming equations are

\[
\mathcal E_B^{\rm local}[j_\Sigma g]
+\mathcal E_B^{\rm reference}[j_\Sigma g]
+\mathcal E_B^{\rm matter}[g_\Sigma,C_\Sigma[g]]=0,
\qquad B=N,\beta.
\]

Their history derivative includes the already owned retarded
`delta C_Sigma[g]`. The old fixed-C0 polynomial coefficients remain useful
local/reference responses, but their constant matter insertion is no
longer the complete source functional for this problem.
No old fixed-C0 root would become physical merely by relabeling it.

This formulation retains the same action and source preparation while
changing the frozen incoming-state condition. No new stress, coupling,
source occupation, phase, channel inventory or physical duration is fitted.
The intrinsic incoming identification and fixed rho0 seam remain unchanged.

## Implementation and proof requirements

Use the existing evolved source columns and their retarded tangents as the
computational representation. Retain complete source fibers and coherent
cross-position kernels. A supplied frozen covariance or a legacy fixed-C0
matter dictionary is not an evolved incoming state and must be rejected by
the new constraint interface. Source-fixed means delta C_src=0 explicitly;
it does not mean delta C_Sigma=0.

For nonstationary columns, momentum vertices act through axial derivatives
with the owned symmetric ordering. The source energy remains a source label;
substituting k=-E as the outgoing momentum would silently restore the old
stationary assumption. Quadrature measures and multiplicities are retained
once, and finite sampled controls do not establish full spectral accuracy.

First wire the state and its derivative into both incoming constraints.
Then close or bound the low-energy/subgap and changed-history source errors
needed for the actual residual claim, and solve under this state law. A
result is EXISTENCE with actual history/state residuals, a necessary scoped
NON-EXISTENCE in this class, or OPEN with the exact remaining error or input.
Old fixed-C0 roots are never promoted by relabeling. Extended stationarity
and metric evolution retain their existing existence gates; publication
requires Douglas's request.

## Retained upstream compensation lemmas

At the reference $N_K=1,\beta_K=0,k=-E_*$, the actual lapse, axial-scale
and radius vertices have Pauli coefficient matrix

$$
V=\begin{pmatrix}-m&0&0\\\ell/r&0&-\ell/r^2\\
-E_*/a&E_*/a^2&0\end{pmatrix},\qquad
\det V=-\frac{m\ell E_*}{a^2r^2}\ne0.
$$

A common nonnegative, nonzero smooth upstream bump gives
$K_{\rm comp}=-i\mu_b(\overline V_b\alpha)\cdot\sigma$,
$\mu_b=\widehat w(0)\int b/a_0>0$, where $\overline V_b$ is the positive
weighted mean of the transported vertices. The reference coordinate factor
$a_0$ remains fixed under raw metric variation. For group `14_1` at
`E=0.5562120090641313`, every such bump supported in
$[203/200,41/40]$ has $\sigma_{\min}(\overline V_b)>51/200$;
the physical margin is $\mu_b,51/200$, with no bump-independent lower
bound on $\mu_b$. The rational argument uses $\rho_c=51/50$, halfwidth
$1/200$, $\sigma_{\min}(V_c)>3/8$ and variation bound $3/25$.

These lemmas permit cancellation of one selected diagonal traceless
response while leaving incoming metric jets unchanged. They choose no
amplitudes or finite history. A physical metric must use the same functions
and coefficients at every energy and channel, respect the union-of-supports
retarded guards, and retain positive fields. One-pair rank is neither full
$C_0$ preservation nor constraint closure. Under the adopted source-fixed
law, these are reusable regression mathematics rather than a requirement to
restore the old frozen covariance.
