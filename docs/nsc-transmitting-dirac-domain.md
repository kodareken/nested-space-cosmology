# A transmitting Dirac trace domain at the stored spacelike seam

This definition binds all 33 stored reduced channels to parent and child
boundary data on **the same** hypersurface. It uses the already declared
transparent condition, the serialized Cauchy basis and the authenticated
child/PG coframe map. It accepts no supplied link matrix.

The [record](../results/development/nsc-transmitting-dirac-domain.json) and
`TransmittingDiracSeamDomain` define a finite retained spectral trace domain.
They do not claim a complete continuum domain or a selected metric history.

## 1. The surface and its normal

The inherited PG coframe is

$$
\theta^0=N\,d\tau,\qquad
\theta^1=q_{\rm PG}(d\rho+\beta\,d\tau).
$$

Define the seam by $X_{\Sigma_0}(\tau,\theta,\phi)=(\tau,0,\theta,\phi)$.
For the stored geometry, $N=q_{\rm PG}=r=1$ and
$\beta=\sqrt{3\pi/2}$. The induced axial scale is

$$
a_\Sigma=\sqrt{q_{\rm PG}^2\beta^2-N^2}=1.9267560770332839.
$$

Thus $\Sigma_0$ is a spacelike three-surface. It is distinct from its
two-sphere intersection with a constant-PG-time slice. On the seam the
child spatial coordinate is $z=\tau$, while child future evolution proceeds
toward decreasing $\rho$. These facts come from the existing
[coframe map](../results/development/adm-neck-source-map.json).

The domain extends algebraically to positive $N,q_{\rm PG},r$ with
$q_{\rm PG}\beta>N$, keeping the coordinate cut fixed. Its future normal is

$$
n^\mu=\left(\frac{q_{\rm PG}\beta}{Na_\Sigma},
             -\frac{a_\Sigma}{Nq_{\rm PG}}\right)_{(\tau,\rho)}.
$$

The parent region $\rho>0$ has this future outward normal at the cut;
the child region $\rho<0$ has its negative. Null and timelike cuts are
rejected by this particular domain. A moving throat or the selected
$\Sigma_L$ requires its own embedding and Cauchy transport.

## 2. Stored channels to boundary traces

The stored Pauli convention is
$H_c=h_x\sigma_1+h_y\sigma_2+h_z\sigma_3$, with axial current $\sigma_3$.
The existing PG/current representation uses $\sigma_2$. The explicit
unitary basis map is

$$
R=\frac{I+i\sigma_1}{\sqrt2},\qquad
R\sigma_3R^\dagger=\sigma_2,\quad
R\sigma_2R^\dagger=-\sigma_3.
$$

The entire reduced operator and covariance transform together. In particular,
$R H_cR^\dagger=h_x\sigma_1-h_y\sigma_3+h_z\sigma_2$. This is a basis
convention within the reduced channels, not a selection of a four-dimensional
chiral sector or a new particle species.

The coframe conversion is separate. On this same seam its normal-current
matrix and rapidity are

$$
J_\Sigma=\frac{q_{\rm PG}\beta I-N\sigma_2}{a_\Sigma}
=e^{-\xi\sigma_2},\qquad
\tanh\xi=\frac{N}{q_{\rm PG}\beta}.
$$

$\xi$ is the inherited child/PG coframe rapidity; it is not the small
Landau rapidity associated with the other Cauchy normal.

For each stored frequency weight $w_j>0$, define the map from an independent
canonical coefficient $a_j$ to its coordinate spinor trace by

$$
\psi_j=T_j a_j,\qquad
T_j=\frac{e^{+\xi\sigma_2/2}R}{r\sqrt{a_\Sigma w_j}},\qquad
G_j=a_\Sigma r^2w_jJ_\Sigma.
$$

It obeys $T_j^\dagger G_jT_j=I$. The raw trace covariance is
$T_jC_jT_j^\dagger$; pulling it back recovers the stored canonical $C_j$.
CAR is assessed in that canonical space, not by applying an ordinary
Euclidean eigenvalue bound to unnormalized coordinate traces.

The stored frequency is already on the parent-neck coordinate axis. The
payload verifies $a_\Sigma h_{z,j}=-\omega_j$, so the seam spatial momentum
is $k=-\omega$. The original quadrature weights, incoming/outgoing occupations,
angular and compact labels, multiplicities and source factors are retained.
No second $\Omega$ rescaling is applied. Degeneracy-compressed channels stay
compressed; unresolved angular spinor functions are not manufactured.

## 3. The transmitting domain

Reuse the existing common-frame condition

$$
\psi_p|_{\Sigma_0}=\psi_c|_{\Sigma_0}.
$$

In the retained spectral representation, allowed pairs are the graph

$$
E_j a_j=(T_j a_j,T_j a_j),\qquad
E_j^\dagger\begin{pmatrix}G_j&0\\0&-G_j\end{pmatrix}E_j=0.
$$

The graph is maximal isotropic for the oriented finite Dirac boundary form.
Across the 1,904 nodes it has 3,808 independent complex coefficients inside
a 7,616-dimensional trace-pair space. The two traces are restrictions of
one field. They are not two independent canonical particle copies; their
joint trace covariance must not be passed as a doubled CAR state.

This applies the declared same-field transmission rule to the actual stored
channel basis. It does not choose a free dynamical $V_c$, identify a pointwise
boost with $U_{L0}$, or propagate the seed covariance onto $r_\star$.

## 4. Definition gate and what follows

**T — DEFINITION PASS in the finite retained seed-seam domain.** The verifier
evaluates the spin/current map, geometric normal and measure, all 1,904 trace
normalizations and covariance recoveries, and the matching graph. It rejects
bare `LinkHistory.values`, an attempted conversion of the spacelike sewing
relation to an instantaneous hopping block, and a copy of the seed state to
the selected radius. Source and map bytes are authenticated.

**U remains OPEN.** A spacelike same-field sewing domain does not define the
independent simultaneous-room decomposition assumed by the given-$B$
Hamiltonian control. A physical covariant parent/child representation and
its conversion to that Hamiltonian must be specified before the sewing data
can determine $B[g,X_\Sigma]$ and its jets. The boundary norm $G_j$ is not
renamed as that energy-valued link.

The fixed-seam normal now has an explicit metric formula. Moving-/selected-
surface jets, $\Gamma_{\rm rest}$'s metric and normal-jet boundary dependence,
the physical two-sided Weyl mismatch and full stationarity remain unevaluated.
The eight Weyl coefficients and their $93.54264532195464$ diagnostic maximum
are unchanged. The previous non-existence, CTP differential and bulk-jet
certificates are retained; no metric evolution is started.

```sh
python3 scripts/define_nsc_transmitting_dirac_domain.py --check
python3 -m pytest -q tests/test_nsc_transmitting_dirac_domain.py
```
