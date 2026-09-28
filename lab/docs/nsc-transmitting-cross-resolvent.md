# Domain-correct transmitting cross-response and static metric jets

This record chooses **Option B: the whole-line resolvent**, followed by
normalized source/measurement projections. It uses Z1's inherited self-adjoint
Dirac operator and T's transparent seam. No given link matrix, reflecting
wall, seed-state remapping or metric trajectory is introduced.

## 1. The cross object and its finite measurements

Keep the convention

$$
R(z)=(H_D-z)^{-1},\qquad G_{pc}(z)=P_pR(z)P_c.
$$

The new numerical matrix is the compression
$\mathcal R(z)=J^\dagger R(z)J$, where $J$ contains two spin components of
each of the normalized spatial packets

$$
f_p(\rho)=\sqrt{30}\,\rho(1-\rho)\quad(0<\rho<1),
\qquad f_c(\rho)=f_p(-\rho)\quad(-1<\rho<0).
$$

Both vanish outside their support. They are independent $L^2$ probes in
$D(H)$, not four new physical modes or an initial quantum state. The ports
$-1,0,+1$ and the two nonreal energies are reused from the existing boundary
record. The parent probe lies in the trapped part of the parent region;
this calculation does not evaluate an exterior source beyond the horizon.

## 2. Retarded and advanced restrictions without finite walls

Both characteristic velocities are negative throughout this probe cell.
For $\mathrm{Im}\,z>0$, the whole-line resolvent is the forward causal
Laplace transform $R(z)=i\int_0^\infty e^{izt}e^{-iH_Dt}dt$. Consequently
the response is zero above the source support. For the lower half-plane the
direction reverses. This uses the already established propagation domain.

With $U_z(\rho)=R(z)J$ and the inherited multiplication term $M$, evaluate

$$
U_z'=v^{-1}\big[i(zI-M)-v'/2\big]U_z+i v^{-1}J(\rho).
$$

The upper-half-plane calculation starts with $U_z(1)=0$ and integrates toward
$-1$; the lower-half-plane calculation starts with $U_z(-1)=0$ and integrates
toward $+1$. These are causal restrictions of the unique whole-line
resolvent, not boundary conditions on a reflecting box. The outflow endpoint
is not assigned a vanishing value. The two integrations meet continuously
at $\rho=0$.

For the upper half-plane, $\mathcal R_{pc}=0$ while
$\mathcal R_{cp}$ is the evaluated parent-to-child response. The lower
half-plane satisfies $\mathcal R(\bar z)=\mathcal R(z)^\dagger$.
Zero reverse response is a causal statement; it is not a declaration $B=0$.

The individual projected terms $iv(0)U_z(0)\delta(\rho)$ are generally
nonzero. Their opposite parent/child contributions cancel because the full
solution, and its metric sensitivities, share the same seam trace. They are
never inserted as a new source or renamed as a link matrix.

## 3. A generating function instead of an instantaneous B

Define the inverse projected response

$$
\mathcal K(z)=\mathcal R(z)^{-1}.
$$

Its off-diagonal block $\mathcal K_{cp}(z)$ is an energy-dependent cross-link
generating kernel in the declared packet basis. The full projected response,
including its diagonal blocks, determines it. A cross block alone would not
fix this normalization.

The standard Feshbach identity explains its domain-correct meaning. Since
the finite probe functions belong to $D(H)$, their projector $P_J=JJ^\dagger$
preserves the operator domain; $Q_J=I-P_J$ does too. Eliminating the actual
remaining bulk gives

$$
\mathcal K(z)=J^\dagger H_DJ-zI
-J^\dagger H_DQ_J(Q_JH_DQ_J-z)^{-1}Q_JH_DJ.
$$

The complement is existing bulk field content, not an added bath or particle
sector. Its contribution is retained through the whole-line resolvent. This
finite packet projector is distinct from the domain-unsafe sharp half-line
projectors. No instantaneous parent/child $B$ is claimed, and
$\mathcal R$ is not assumed to be the resolvent of a closed four-mode
Hamiltonian.

## 4. Metric derivatives of the same response

For fixed canonical probes and the fixed coordinate cut, use the imported
identity $\delta R=-R\delta H_D R$. The implementation computes its static
functional kernel for compactly supported variations of
$(\log N,\beta,\log q_{\rm PG},\log r)$ inside $(-1,1)$:

$$
\delta\mathcal R=\sum_A\int s_A(\rho)\mathcal J_A(\rho;z)\,d\rho,
\qquad
\delta\mathcal K=-\mathcal K\,\delta\mathcal R\,\mathcal K.
$$

With $U=R(z)J$, $L=R(\bar z)J$, $V_A=\partial v/\partial g_A$ and
$M_A=\partial M/\partial g_A$, integration by parts in the owned Dirac
vertex gives

$$
\mathcal J_A=\frac{i}{2}(L^\dagger V_AU'-L'^\dagger V_AU)-L^\dagger M_AU.
$$

The compact variations vanish at the outer ends; the internal seam terms
cancel by continuity. The test profile is a smooth compact spatial source,
not a prescribed metric history or an observational fit.

The linearized radial solve is compared with the independent weak-resolvent
quadrature, and with finite differences for one existing channel from each
family. Static resolvent jets are stored as such. They cannot be inserted
directly into `EndpointBranchJets`, which requires the full unitary
time-history derivative. That conversion still needs the spectral/time
reconstruction and the physical state/memory data.

## 5. Focused verification and scope

The verifier checks the adjoint relation, seam continuity, positive imaginary
part, the resolvent norm bound and the true resolvent identity

$$
\mathcal R(z)-\mathcal R(w)
=(z-w)\int U_{\bar z}(\rho)^\dagger U_w(\rho)\,d\rho.
$$

The advanced and retarded tails have disjoint support outside $[-1,1]$, so
this integral is the whole-line product for these probes. It is not replaced
by the generally incorrect product of two compressed matrices.

**Z2a passes for the declared cross-response/generating function. Z2b passes
for its compact static metric kernels.** A constant Hamiltonian link and
time-history CTP jets are not claimed. Moving-neck jets remain out of scope.
$\Gamma_{\rm rest}$, the physical Weyl mismatch, a selected recursive state
and metric stationarity remain open. All locked coefficients and earlier
certificates are preserved.

The [result record](../results/development/nsc-transmitting-cross-resolvent.json)
authenticates its input owners and the numerical payload.

```sh
python3 scripts/derive_nsc_transmitting_resolvent.py --check
python3 -m pytest -q tests/test_nsc_transmitting_resolvent.py
```
