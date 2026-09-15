# What the Weyl endpoint derivative actually varies

The existing Weyl action now has an explicit boundary-jet covector. It
separates the physical kinds of boundary variation from the frozen discrete
end-node gradient whose largest magnitude is $93.54264532195464$. Every
stored coefficient is retained and reconstructed; no compensating force is
inserted.

## Apply the existing density

Use the same $C_W$, KS fields and Weyl density as
`GeneralKSLocalInducedHistory`, per unit homogeneous axial coordinate:

$$
L_W=k,N a r^2 F^2,\qquad k=-\frac{16\pi}{3}C_W,
$$

$$
F=\frac{\ddot a-\dot a\dot N/N}{N^2a}
-\frac{\dot a\dot r}{N^2ar}
-\frac{\ddot r-\dot r\dot N/N}{N^2r}
+\frac{1+\dot r^2/N^2}{r^2}.
$$

Import the higher-derivative first-variation identity. For this density its
boundary part is

$$
\left[P_N\delta N+P_a\delta a+P_r\delta r
+Q_a\delta\dot a+Q_r\delta\dot r\right]_0^L,
$$

where $Q_a=\partial L_W/\partial\ddot a$,
$Q_r=\partial L_W/\partial\ddot r$ and
$P_i=\partial L_W/\partial\dot g_i-dQ_i/dT$; $Q_N=Q_\beta=P_\beta=0$.
The new owner evaluates these independent jet derivatives, rather than
identifying them with a displacement of one time-grid node.

## The normal rate is shear

The same density gives

$$
Q_a=\frac{2kr^2F}{N},\qquad
Q_r=-\frac{2karF}{N},\qquad
aQ_a+rQ_r=0,
$$

$$
N P_N+\dot a Q_a+\dot r Q_r=0.
$$

Changing to the proper normal rates
$H_a=\dot a/(Na)$ and $H_r=\dot r/(Nr)$ therefore expresses the boundary
form as

$$
\left[\overline P_a\delta a+\overline P_r\delta r
+\Pi_\sigma\delta(H_a-H_r)\right]_0^L,
\qquad \Pi_\sigma=NaQ_a.
$$

The independent rate is the difference between axial and spherical expansion.
The lapse boundary coefficient cancels under this change of variables by
the displayed identity. This is an application of the existing Weyl action,
not a new shear source or physical prediction.

On the stored diagnostic jets, $\Pi_\sigma$ is approximately $-0.84788331$
at the initial end and $0.85747641$ at the final end, before applying the
final-minus-initial boundary orientation. The history remains the previously
specified control, not a selected duration or spacetime solution.

## Exact relation to the saved eight numbers

The local owner uses trapezoidal weights $W$ and the `edge_order=2`
first-derivative matrix $D$. Define $B=WD+D^TW$ and
$S=\mathrm{diag}(-1,0,\ldots,0,1)$. Its discrete action gradient decomposes
exactly as

$$
G=W E+SP+D^TSQ+(B-S)P+D^T(B-S)Q.
$$

The terms are the bulk chain-rule contribution, canonical endpoint-jet
contribution and the stencil's integration-by-parts defect. Here $E$ is
bookkeeping for the existing density variation, not a newly solved metric
equation. The existing stencil has $\lVert B-S\rVert_2=0.86237244$;
setting this defect to zero would change the discrete calculation.

For the final radius node the three contributions are approximately

$$
0.68615308-87.27311859-6.95567981=-93.54264532.
$$

The [record](../results/development/nsc-weyl-boundary-jets.json) contains the
full precision decomposition at every node. Its reconstruction differs from
the immutable old gradient by at most $9.13\times10^{-12}$, below the declared
$3\times10^{-11}$ tolerance. The original eight-coefficient binding remains
exactly zero. The highest-jet and proper-lapse identities hold below
$3.14\times10^{-16}$.

## What this changes in the closure calculation

Fixing the intrinsic endpoint metric alone does not fix its normal shear.
Fixing both metric and normal rates makes the canonical Weyl boundary form
vanish on those allowed variations; this is a boundary-data condition, not
evidence that an unfixed interface self-sources. A free or glued physical
interface requires the actual matching of metric, normal shear and their
oriented conjugate data. The beginning and end of the diagnostic history
are not automatically the two sides of the fixed $\rho=0$ interface.

Consequently, matching an unspecified extra force directly to $93.5426$ is
not a substitute for stating that variational domain. The new boundary-jet
object supplies the formerly missing local first-jet coefficients while
preserving the old nodal basis and all its numbers.

The Gaussian Schur determinant is a factorization of the already counted
full Gaussian CTP action. It cannot be added again as $\Gamma_{\rm rest}$.
Euler and $\Box R$ endpoint terms are already in the local action. The
declared remaining transmitting functional and its permitted boundary
variations still need their physical identification; this calculation assigns
neither an extra term nor zero to them.

**Boundary-jet extraction PASS; physical B1/B2 stationarity OPEN.** No
physical two-sided mismatch, new stress, selected history or metric timestep
is asserted. The active work is the actual intrinsic-jet gluing/endpoint
domain and remaining same-action derivative, rather than more adapters for
the diagnostic end-node vector.

```sh
python3 scripts/derive_nsc_weyl_boundary_jets.py --check
python3 -m pytest -q tests/test_nsc_weyl_boundary_jets.py
```
