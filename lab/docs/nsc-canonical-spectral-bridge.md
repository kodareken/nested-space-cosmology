# Canonical compact fields and their covariant spectral conversion

The stored compact spectral source now has an explicit allocation to its
canonical field content, light cutoff conversion and warp/covariant-measure
conversion. This identifies the actual physical field split to which the
[real-time Gaussian construction](nsc-boundary-state.md) can apply. It does
not treat a difference of spectral weights as a new bath.

The [record](../results/development/canonical-spectral-bridge.json) is
assembled from the authenticated
[compact-matching record](../results/development/compact-matching.json).
The owner is [nsc_canonical_spectral_bridge.py](../src/recursive_horizons/nsc_canonical_spectral_bridge.py).
No compact eigenvalue, momentum integral or angular source is recalculated.

## The physical fields are already determined for the free domain

Reuse the [canonical compact reduction](nsc-compact-mass-map.md). For the
two complementary five-dimensional Dirac copies, with the fixed length-ell
chiral interval, the free four-dimensional content is

\[
\mathcal H_{5,\mathrm{can}}
=\mathcal H_{D,0}\oplus
\bigoplus_{n\ge1}(\mathcal H_{D,m_n}\oplus\mathcal H_{D,m_n}),
\qquad m_n=\frac{n\pi}{\ell},\qquad \ell=2L_\star.
\]

The two Weyl zero modes form one Dirac field; every positive level has two
Dirac copies. These are fields of the declared free realization, not an
identification with observed particle species. The masses use the conformal
coordinate interval ell, not the warped proper length ell J.

The conformal field normalization chi=exp(2s sigma)Psi removes the warp
from the classical free evolution in the g4 frame. For static general
lapse N, the canonical mass term is m_n N beta; it is m_n beta in the
unit-lapse PG frame. A constant-size variation at fixed g4 gives
partial_ell H_n=-(m_n/ell)N beta. A time- or space-dependent size also
changes the reduction and must retain its additional geometric terms.

For a mode-block-diagonal Gaussian preparation and the free block-diagonal
Hamiltonian, the real-time normalized functional factorizes at finite
mode number:

\[
\Gamma_{F,5}^{\mathrm{CTP}}[C]
=\Gamma_{F,0}^{\mathrm{CTP}}[C_0]
 +\sum_{n\ge1}\sum_{a=1}^2\Gamma_{F,m_n}^{\mathrm{CTP}}[C_{n,a}].
\]

The infinite limit needs the common covariant normalization. Intermode
initial correlations require the full covariance determinant, and the
already derived interactions can mix levels. No truncation of those
interactions is inferred here. For this free, fixed-size product domain
the light-heavy fermion link vanishes. Heavy fields still contribute to
the metric determinant and its stress.

## Separate the three contributions exactly

The existing Euclidean weight at zero warp is

\[
h_{\Lambda,0}(y)=E_1(y/\Lambda^2)
 +2\sum_{n\ge1}E_1((y+m_n^2)/\Lambda^2).
\]

Let h_Lambda,s be the already computed warped weight. Define

\[
\begin{aligned}
k_\Lambda(y)&=2\sum_{n\ge1}E_1((y+m_n^2)/\Lambda^2),\\
d_{\Lambda,\nu}(y)&=E_1(y/\Lambda^2)-E_1(y/\nu^2),\\
j_{\Lambda,s}(y)&=h_{\Lambda,s}(y)-h_{\Lambda,0}(y).
\end{aligned}
\]

Then the previously computed complement obeys

\[
\boxed{H_\nu(y)=k_\Lambda(y)+d_{\Lambda,\nu}(y)+j_{\Lambda,s}(y).}
\]

The three terms are the finite proper-time weight of the canonical massive
fields, the light endpoint conversion, and the change in the covariant
modulus under the warp. The latter includes the declared measure and
fixed-cutoff prescription; it is not a separately adjustable source or a
regulator-independent Jacobian asserted without a definition.

The warp conversion is a specified Euclidean functional,

\[
\mathcal J_{\Lambda,s}^E[g_4]
=\tfrac12\operatorname{Tr}_4 j_{\Lambda,s}(D_4^2)
=\Gamma_{5,\Lambda,s}^E-\Gamma_{5,\Lambda,0}^E.
\]

Its full metric variation uses the existing spectral calculus, including
measure/domain transport:

\[
\delta\mathcal J^E=\operatorname{Tr}_4
[D_4 j'_{\Lambda,s}(D_4^2)\,\delta D_4].
\]

Its warp variation is the already derived dual heat trace, not a new
coefficient. With L_s=D_s^dagger D_s and R_s=D_s D_s^dagger,

\[
\partial_s\mathcal J^E
=\tfrac12\operatorname{Tr}\sigma
(e^{-L_s/\Lambda^2}+e^{-R_s/\Lambda^2}-P_{\ker L_s}-P_{\ker R_s}).
\]

The normalization, adjoint domains and kernel terms are those of the
[warped-source owner](nsc-warped-source.md). They are not discarded by
calling this a conformal anomaly. The underlying heat/Dirac and conformal
methods are established, for example in
[Vassilevich, sections 3.3 and 7.1](https://arxiv.org/html/hep-th/0306138v3).

## The existing coefficients have definite origins

Linearity of Q0, Q1 and Q2 gives the same three-part allocation for the
stored vacuum, Einstein, gauge and curvature-squared terms. In particular,

\[
Q_0[d]=2\log(\Lambda/\nu),\quad
Q_1[d]=\Lambda^2-\nu^2,\quad
Q_2[d]=(\Lambda^4-\nu^4)/2.
\]

The flat-tower Q1 and Q2 are read from the previous analytic calculation.
The previous zero-limit formula gives Q0[k] directly; it requires no
Galerkin operator. The warp moments are the stored warped-minus-flat
moments. The same linear coefficient map is applied to each column.

For the stored development values Lambda=2, ell=2, s=1, nu=1:

| Origin | Vacuum coefficient V | Einstein coefficient A | Gauge coefficient C_F |
|---|---:|---:|---:|
| Massive canonical tower, finite proper-time weight | 0.08283595387 | 0.00244705046 | 0.00196700137 |
| Light endpoint conversion | 0.09498860966 | 0.00316628699 | 0.00292627053 |
| Warp/covariant conversion | -0.00528370024 | -0.00011105300 | -0.00005602121 |
| Sum: stored complement | 0.17254086329 | 0.00550228445 | 0.00483725070 |

These are coefficients of V-A R_E+C_F F²+..., not separately measured
cosmic energy densities or physical G and gauge couplings. The allocation
also preserves all the stored Weyl-squared, Euler and total-derivative
coefficients. No additional finite R² coefficient is selected. A local
a4 truncation remains conditional on small external curvatures and momenta.

The record includes every column at the three previously recorded matching
cutoffs, the stored target fields and all reconstruction residuals. The
warp and canonical-tower columns are independent of nu. Their assembled
moments agree with the stored complement within 5e-11 and their coefficients
within 5e-13; the small differences retain the previous refined-versus-base
moment integration differences. No new convergence claim is made.

## Join the canonical state without double-counting the vacuum

Let Gamma_can,5,Lambda^E denote the unwarped canonical tower with the same
finite proper-time prescription. Before any infinite-mode limit, define
its conversion to a specified renormalized canonical determinant by

\[
\mathcal C_{5,\Lambda,\mu}^E
=\Gamma_{\mathrm{can},5,\Lambda}^E
 -\Gamma_{\mathrm{can},5,\mathrm{ren},\mu}^E.
\]

The complete identity is therefore

\[
\boxed{\Gamma_{5,\Lambda,s}^E
=\Gamma_{\mathrm{can},5,\mathrm{ren},\mu}^E
 +\mathcal C_{5,\Lambda,\mu}^E+\mathcal J_{\Lambda,s}^E.}
\]

The auxiliary light cutoff nu has disappeared. This joins the previously
separate canonical field map and warped modulus without a second gravity
weight or a new matter field. On the free fixed-ell, fixed-g4 orbit, the
canonical Hamiltonian and a fixed canonical preparation have no s
dependence; the retained spectral warp source belongs to J.

For a physical real-time source, the canonical state determinant is now
specified once its preparation is given. The remaining work is the common
real-time prescription for C5+J, including its measure, phase and boundary
data. A Euclidean subtraction identity does not by itself define that
prescription. The previous state-regulator result is retained: its
state-dependent conversion cannot be replaced by a local vacuum constant.
This record neither replaces the full action by a low-energy truncation
nor claims that the source equations or self-sourced geometry are solved.

## Reproduction

```sh
python -B scripts/check_nsc_canonical_spectral_bridge.py --check
```

This authenticates the existing dependency closure, performs only the new
allocation, and compares every record field with the declared policy. It
does not call any old spectrum or source generator. New output paths refuse
overwrite. Full physical closure remains the goal of the common source
equation, with this decomposition as a completed input.
