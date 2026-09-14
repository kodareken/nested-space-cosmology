# A physical PG Cauchy preparation for the transparent LLL channel

The massless lowest-Landau channel now has an explicit covariance on a common
PG-time slice, including correlations with the bulk outside the original
packets. It uses the already declared affine-horizon state and inherited
incoming occupation. The stored seed covariance is a verification input,
not the source of a copied spatial matrix.

**Scope:** one of the 33 retained channel groups, with its locked magnetic
multiplicity. The 32 massive angular/compact groups are not assigned this
massless preparation. Full C1b and C2--C4 remain open.

## 1. The global domain and the missing horizon partner

Import the Z1 whole-line self-adjoint PG operator. In the LLL, $m=\lambda=0$,
and the two characteristic velocities are $a_\pm=\pm1-\beta(\rho)$.
The horizon zero of $a_+$ separates its interior and exterior flow sectors;
$a_-$ crosses the whole radial line. Thus the same two-component field has
three independent characteristic sectors:

$$
P_{e+}=\mathbf{1}_{\rho>\rho_h}P_+,
\qquad P_{i+}=\mathbf{1}_{\rho<\rho_h}P_+,
\qquad P_-=\frac{I-\sigma_2}{2},
\qquad P_+=\frac{I+\sigma_2}{2}.
$$

They sum to the identity; no field or particle copy is added. The original
parent/child cut remains $\rho=0$. The spacelike seed cylinder is sufficient
for trapped-domain propagation but does not intersect the exterior outgoing
characteristics. Its two covariance slots alone cannot determine that third
sector. The existing horizon preparation supplies the missing partner and
its correlations, rather than an arbitrary occupation being inserted there.

## 2. State and unitary spatial map

Use the imported characteristic coordinates and half-density transformation:

$$
\frac{dx_s}{d\rho}=\frac{1}{a_s(\rho)},\qquad
\chi_s(\rho)=\frac{\widetilde\chi_s(x_s(\rho))}{\sqrt{|a_s(\rho)|}}.
$$

Each sector becomes the already known generator $-i\partial_{x_s}$ on a
complete flow line. This is a unitary change of Cauchy representation, not
a Lorentz boost substituted for evolution. The exterior/interior $x_+$
coordinates share the logarithmic finite part at the same horizon. A common
translation is immaterial. This preserves the relative affine phase of the
existing horizon covariance and introduces no transport duration.

For all signed Killing frequencies $E$, keep

$$
f(E)=\frac{1}{1+e^{\beta_H E}},\qquad
n(E)=\frac{1}{1+e^{\beta_{\mathrm{in}}E}},\qquad
\beta_H=\frac{2\pi}{\kappa_h},\qquad \beta_{\mathrm{in}}=\frac{\beta_H}{\Omega}.
$$

In the order $(e+,i+,-)$ the already declared preparation is

$$
\mathsf C(E)=
\begin{pmatrix}
f&-i\sqrt{f(1-f)}&0\\
i\sqrt{f(1-f)}&1-f&0\\
0&0&n
\end{pmatrix}.
$$

The horizon block is the original $C_H$ from the
[parent-matched state](nsc-unruh-state.md); the incoming occupation is the
one used in the [charged source](nsc-charged-ctp-neck-source.md).
Its eigenvalues are $0,1,n(E)$, so the full direct-integral covariance obeys
$0\le C\le I$. The unitary flow map carries this bound onto the PG slice.
No projection onto negative Killing energy is used to invent a global
vacuum: the inside and outside outgoing occupations are complementary.

At the seed, restriction to the interior outgoing and throughgoing incoming
sectors gives $\mathrm{diag}(1-f,n)$. The existing T coframe/half-density map
then recovers every stored LLL covariance block. This comparison does not
assign the seed frequency matrix to a spatial packet.

## 3. Spatial covariance and probe/bulk correlations

The global covariance kernel is the Fourier transform of $\mathsf C(E)$,
with the sector half-density factors. The standard thermal chiral Fourier
kernels are used with the inherited sign convention: the diagonal Fermi
block has contact term $\delta(x-y)/2$ and regular part
$-i/[2\beta\sinh(\pi(x-y)/\beta)]$; its complement has the opposite regular
sign. The horizon cross block is
$-i/[2\beta_H\cosh(\pi(x_e-x_i)/\beta_H)]$.
These are imported Fourier identities; the NSC calculation is their
normalization, affine matching and projection on the actual PG geometry.

The contact term is retained in every smeared covariance. No finite
frequency grid is treated as the complete CAR kernel. The code refuses a
coincident-point value; it does not assign a finite vacuum stress there.

The original four packet coordinates in $J$ are unchanged. Three additional
orthogonal test coordinates inspect the eliminated bulk: two spin components
of an odd polynomial packet on $(-1,0)$ orthogonal to the child packet, and
one outgoing exterior packet on $(3,4)$. These are observables of the same
field, not new source channels. Their seven-dimensional compression records
$C_{JJ}$, $C_{Jb}$ and $C_{bb}$ without factorizing the state. Those three
extra coordinates do not exhaust the eliminated bulk; the global kernel
defines the rest.

At equal PG time this supplies actual state-dependent contour data:

$$
G^<=iC,\qquad G^>=-i(I-C),\qquad G^K=-i(I-2C).
$$

The independent characteristic integral of the retarded propagator is
compared with the authenticated LLL packet resolvent. This binds the
preparation to the same memory operator used in C1. The projected covariance
is not inserted as a seven-mode closed-system state into `influence`; its
complement and memory remain part of the full Gaussian construction.

## 4. Remaining gate

**LLL C1b preparation PASS; full retained C1b OPEN.** The physical PG state,
its horizon partner and packet/bulk correlations are evaluated for this exact
channel. The massive groups require their own global mode reconstruction,
including the signed-frequency spin/angular map and original horizon/parent
preparation. Their multiplicities do not define that map.

Full transmitting `EndpointBranchJets`, the KS endpoint pullback,
$\Gamma_{\rm rest}$, the two-sided Weyl mismatch and extended stationarity
remain unevaluated. The eight Weyl coefficients and previous certificates
are unchanged. This construction selects no metric history and starts no
metric timestep; moving-neck Z3 remains outside scope.

The [record](../results/development/nsc-pg-lll-preparation.json) serializes the
new covariance, contour matrices, correlations, residuals and source hashes.

```sh
python3 scripts/derive_nsc_pg_lll_preparation.py --check
python3 -m pytest -q tests/test_nsc_pg_lll_preparation.py
```
