# Independent bulk spaces on the common PG time slice

**Z1 defines the parent and child bulk subspaces and their projectors.** The
definition uses the inherited whole-line PG time $\tau$ and the 33 retained
reduced Dirac channels. T remains the transmission condition on the seam;
its two complete traces are not reused as bulk canonical variables.

The [record](../results/development/nsc-common-time-bulk-split.json) verifies
the independent support split, its compatibility with T, and the subsequent
operator-domain test for $P_pH_DP_c$.

## 1. One common time and independent bulk supports

On a $\tau=$ constant slice, the canonical field is
$\chi=r\sqrt{q_{\rm PG}}\,\psi$. Per retained channel, define

$$
\mathcal H_\tau=L^2(\mathbb R_\rho,d\rho;\mathbb C^2),\qquad
\mathcal H_p=L^2((0,\infty);\mathbb C^2),\quad
\mathcal H_c=L^2((-\infty,0);\mathbb C^2).
$$

The full retained space is the direct sum over the existing channels and
multiplicity labels. The bounded orthogonal projectors are

$$
(P_p\chi)(\rho)=\mathbf1_{\rho>0}\chi(\rho),\qquad
(P_c\chi)(\rho)=\mathbf1_{\rho<0}\chi(\rho).
$$

The single point $\rho=0$ has zero $L^2$ measure; its one-sided values belong
to operator-domain traces. Hence $P_pP_c=0$ and $P_p+P_c=I$ on the Hilbert
space. Canonical fields with disjoint bulk supports have the imported CAR
relations with zero cross anticommutator. This independence does not prescribe
a factorized quantum state.

The time foliation is already owned by the
[Dirac tetrad result](../results/nsc-4-dirac-tetrad.json):
$g^{\tau\tau}=1/N^2>0$, including inside the horizon. The whole spacelike
three-surface $\Sigma_0:\rho=0$ is different from its two-sphere intersection
with a common-time slice. Enforcing transparent traces at that intersection
for every $\tau$ recovers T's seam condition.

## 2. Operator and transmitting domain

Use the existing symmetric Dirac principal part in the current convention
fixed by T:

$$
H_c[g]=-i\big(v\partial_\rho+\tfrac12v'\big)+M_c,
\qquad v=\frac{N}{q_{\rm PG}}\sigma_2-\beta I,
$$

$$
M_c=-Nm_c\sigma_1-N\frac{\lambda_c}{r}\sigma_3.
$$

The finite $m_c,\lambda_c$ are the retained compact and angular labels.
A constant current-preserving phase rotates the massless term into the
previous $+\lambda_c\sigma_1/r$ convention. The operator is applied in the
common PG foliation; the stored frequency covariance on $\Sigma_0$ is not
silently reinterpreted as a spatial covariance on a $\tau$ slice.

On the locked background, take the inherited self-adjoint graph closure of
$C_c^\infty(\mathbb R;\mathbb C^2)$. Reuse the existing complete-flow and
bounded-perturbation argument in [the tetrad note](nsc-dirac-tetrad.md): the
principal velocities have at most linear growth and complete flows. Adding
the retained multiplication matrices preserves this realization because
$r\ge1$ and $\|M_c\|\le\sqrt{m_c^2+\lambda_c^2}$ for every retained channel.
No new boundary condition at infinity, reflecting cut or periodic carrier is
introduced. Smooth compactly supported coefficient variations keep the same
asymptotic realization; no metric solution is selected.

Near the noncharacteristic seam, domain representatives are piecewise $H^1$
and obey $\chi_p(0)=\chi_c(0)$ in the common smooth frame. This is an
operator-domain condition, not equality of their independent bulk variables.
The Hilbert direct sum is valid while the operator domain is not a product of
independent half-line operator domains.

T's current matrix is consistent with this common-time flux:

$$
G_{T,j}=-q_{\rm PG}r^2w_jv(0).
$$

The record checks this against the authenticated T artifact, without
reconstructing that artifact.

## 3. Finite projector/CAR verification

The executable `CommonTimeBulkSplit` acts on spatial samples on both sides
of the cut. The check uses midpoint quadrature in the previously used
$[-4,4]$ window, solely as a rank/norm witness. No finite-box Hamiltonian or
outer boundary condition is assigned. Per-channel sample counts reuse the
numerical budget; frequency nodes and their covariance are not mapped to
these spatial nodes.

On this finite witness the parent and child ranks are 1,904 each, for a total
independent rank of 3,808. Their extraction/reassembly is a permutation, with
CAR equal to the identity on that rank. This is an actual spatial split of
one set of bulk variables. The locked duplicated-trace mismatch of one is
preserved as a different, failed representation.

## 4. Z2: sharp Hilbert projectors do not preserve the Dirac domain

The imported distributional product rule applied to this NSC principal
matrix gives

$$
H(P_c\chi)-P_cH\chi=i\,v(0)\chi(0)\delta(\rho).
$$

At the inherited seam, $v(0)=\sigma_2-\sqrt{3\pi/2}\,I$ is invertible.
Therefore a nonzero transmitting trace produces a delta distribution when
zero-extended to only one half-line. That zero extension is not in $D(H)$.
This is a projection-domain defect; no physical delta source is added to the
unsplit field. The parent and child distributional terms cancel when the
matching field is assembled.

**Z2 remains OPEN as an ordinary Hamiltonian block.** The product $P_pH_DP_c$
is not defined on general child $L^2$ data through the transmitting operator
domain. Restricting to zero-trace functions gives zero cross action but loses
the transmitting condition; it cannot justify a physical $B=0$. The derived
$iv(0)$ coefficient is not renamed as a finite `LinkHistory` matrix.

A domain-respecting weak/Galerkin or boundary-resolvent realization is needed
before identifying a finite matrix block and its jets. The bounded response

$$
G_{pc}(z)=P_p(H_D-z)^{-1}P_c,\qquad \mathrm{Im}\,z\ne0,
$$

is already well-defined by the whole-line self-adjoint operator and bounded
projectors; it is not evaluated in this definition gate. This supplies a
precise response object without an invented hopping term.

The split keeps the inherited coordinate cut fixed. Moving-surface jets are
out of scope for this truncation. The physical boundary remainder, Weyl
mismatch, selected recursive kernel and stationary geometry remain open.
All prior certificates and scales are unchanged; no metric timestep is run.

```sh
python3 scripts/define_nsc_common_time_bulk_split.py --check
python3 -m pytest -q tests/test_nsc_common_time_bulk_split.py
```
