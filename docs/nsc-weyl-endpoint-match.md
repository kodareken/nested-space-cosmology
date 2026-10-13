# Weyl endpoint: a common covector, with the boundary derivative OPEN

The [endpoint record](../results/development/nsc-weyl-endpoint-match.json)
extracts the stored local action gradient without recomputing the local
history. In the raw nodal coordinates $(N,\beta,q_{\rm ADM},r)$, the first
and last nodes of its 33-node control give:

| Coordinate | Initial-node derivative | Final-node derivative |
|---|---:|---:|
| $N$ | -1.6709564486074024 | -1.7122278287564878 |
| $\beta$ | 0 | 0 |
| $q_{\rm ADM}$ | -43.58985597013112 | 45.868233687373596 |
| $r$ | 83.87051420812537 | -93.54264532195464 |

Extraction agrees exactly with the authenticated local record: residual zero,
tolerance zero. The full four-field nodal gradient is retained alongside
these eight components. The signs are those of $dS/dg$, with the original
time orientation; an extra normal-sign flip is not applied.

The maximum $93.54264532195464$ is the existing **diagnostic end-node local
gradient**. It depends on the control history, node coordinates and
discretization. It is not yet an evaluated physical two-sided boundary
mismatch, and by itself cannot establish non-existence of an extended
stationary history. The original record and this nonzero value are preserved.

## One endpoint-variation object

`EndpointVariation` identifies the history by the local record hash, the
four coordinates, the two node indices, units and action-gradient convention.
A transmitting boundary derivative can be added only after it is pulled
back to that same basis:

$$
R_{\partial,A}
=w_{\partial,A}+b_{\partial,A},\qquad
b_{\partial,A}
=\left[\mathrm{pullback}\!left(
\frac{\delta\Gamma_{\rm rest}}{\delta h}\right)\right]_A.
$$

Here $w$ denotes the recorded Weyl projection; $b$ must identify the residual
action allocation being matched, so already counted local terms are not
added again. A numerical covector sum is only this projected subcheck, not
the complete metric equation. The boundary action must also specify its
embedding and admissible endpoint variations, including the derivative jets
required by the higher-derivative action. Fixing endpoint metric values alone
does not supply those data. A temporal end node is not automatically the
transmitting parent-child hypersurface.

The original endpoint-match receipt marked independent first-jet
derivatives unavailable, with `null` entries for all four fields at both
ends. Its nodal arrays include finite-difference dependence. The later local
boundary-jet extraction below supplies continuum covectors in their own
jet basis; the transmitting functional and its pullback remain missing.

**OPEN:** neither the required transmitting metric derivative nor its
pullback is supplied by the existing action inventory. Accordingly the
boundary vector and the post-matching residual are both `null`, with target
tolerance $3\times10^{-11}$. The code rejects mismatched endpoint bases and
does not interpret an absent derivative as zero or as minus the local vector.

This binds Step B to the same missing action identified in
[Step A](https://doi.org/10.5281/zenodo.23090015). It introduces no endpoint
completion or compensator and selects no control duration as physical.


## Continuum boundary jets and the discrete end-node decomposition

For the same local Weyl density $L_W=kNar^2F^2$,
$k=-16\pi C_W/3$, with

$$F=\frac{\ddot a-\dot a\dot N/N}{N^2a}
-\frac{\dot a\dot r}{N^2ar}
-\frac{\ddot r-\dot r\dot N/N}{N^2r}
+\frac{1+\dot r^2/N^2}{r^2},$$

the independent continuum boundary covector is
$[P_N\delta N+P_a\delta a+P_r\delta r+
Q_a\delta\dot a+Q_r\delta\dot r]_0^L$, with
$Q_i=\partial L_W/\partial\ddot g_i$ and
$P_i=\partial L_W/\partial\dot g_i-dQ_i/dT$.
The exact identities are

$$Q_a=2kr^2F/N,\quad Q_r=-2karF/N,\quad aQ_a+rQ_r=0,
\qquad NP_N+\dot aQ_a+\dot rQ_r=0.$$

With proper rates $H_a=\dot a/(Na)$, $H_r=\dot r/(Nr)$ this becomes
$[\overline P_a\delta a+\overline P_r\delta r+
\Pi_\sigma\delta(H_a-H_r)]_0^L$, $\Pi_\sigma=NaQ_a$.
The lapse boundary coefficient cancels, while normal shear remains
independent of the intrinsic metric. Before final-minus-initial orientation,
the stored diagnostic shear coefficients are approximately -0.84788331
and0.85747641 at the initial/final ends.

For the original trapezoid weights W and edge_order2 derivative D, put
$B=WD+D^TW$ and $S=\operatorname{diag}(-1,0,\ldots,0,1)$.
The full saved nodal gradient reconstructs as

$$G=WE+SP+D^TSQ+(B-S)P+D^T(B-S)Q.$$

The stencil has $\|B-S\|_2=0.86237244$, so its integration-by-parts
defect cannot be discarded. The final radius node's bulk, canonical
boundary and stencil contributions are approximately
$0.68615308-87.27311859-6.95567981=-93.54264532$.
The retired extractor authenticated the33-node history and exact eight
old coefficients, preserved all nodal components, and compared their
reconstruction to the unchanged gradient within
$9.13\times10^{-12}<3\times10^{-11}$; highest-jet/proper-lapse
identities were below $3.14\times10^{-16}$. These are boundary-jet and
discrete-bookkeeping controls, not physical two-sided matching.

Fixing metric and normal rates kills the canonical form on those allowed
variations. Fixing the metric alone leaves free shear. A free/glued
interface needs the actual intrinsic jet, normal shear and oriented
conjugate matching; temporal control endpoints are not automatically its
two sides. Local boundary-jet extraction therefore passes while full
B1/B2 stationarity remains OPEN. The Gaussian Schur determinant factors
the already counted Gaussian CTP action; Euler/BoxR endpoints also already
belong to the local action. None supplies an extra compensator, an assumed
zero transmitting remainder or a physical counterforce to93.5426.
