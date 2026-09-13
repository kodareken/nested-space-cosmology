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

The record explicitly marks the independent first-jet derivatives as
unavailable, with `null` entries for all four fields at both ends. The stored
nodal arrays already include their finite-difference dependence; they do not
separately identify these continuum boundary-jet covectors.

**OPEN:** neither the required transmitting metric derivative nor its
pullback is supplied by the existing action inventory. Accordingly the
boundary vector and the post-matching residual are both `null`, with target
tolerance $3\times10^{-11}$. The code rejects mismatched endpoint bases and
does not interpret an absent derivative as zero or as minus the local vector.

This binds Step B to the same missing action identified in
[Step A](nsc-transmitting-boundary-selection.md). It introduces no endpoint
completion or compensator and selects no control duration as physical.
