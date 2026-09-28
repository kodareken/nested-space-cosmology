# Conditional continuous constraint and assembly bounds

The final certificate cannot replace a continuous supremum by denser nodes.
This primitive accepts a directed bound on each cell's full residual derivative
and intersects the two endpoint Lipschitz cones. It returns the continuous N
and beta maxima and the additional between-node remainder. The derivative
bound must come from the final continuous field/source owner; until then the
gate component remains `None`.

The same module supplies exact-rational `gamma_n` arithmetic for the final
addition of baseline, geometry, matter and edge terms. It deliberately excludes
matrix-product and field-contraction roundoff, which require interval replay of
the final arrays. It therefore establishes the reusable interface without
prematurely filling either budget component.

The production successor is the
[257-node / 256-cell enclosure](nsc-ks-continuous-constraint-enclosure.md)
together with remaining assembly arithmetic and
[Gate Budget v5](nsc-ks-gate-budget-v5.md). The Lipschitz primitive now also
exposes the exact 257/256 verification-grid restriction. Production between-node
and arithmetic components stay `null` until the directed inputs exist.
