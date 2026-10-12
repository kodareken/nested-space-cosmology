# Stable coherent matter difference

The local constraint split needs `G[C_Sigma[g]]-G[C_ref]`. For a prepared
state from the joint reference/difference owner, this contraction can be
evaluated directly from `F=F_ref+dF` without subtracting two large background
values. Here `dF` denotes the finite field difference; the retarded parameter
tangent is retained separately by the prepared-state object.

The density difference includes both cross terms and the quadratic term:

\[
\Delta D=dF C_{src}F_{ref}^\dagger+F_{ref}C_{src}dF^\dagger
          +dF C_{src}dF^\dagger.
\]

Define

\[
K=dF_z C_{src}F_{ref}^\dagger+F_{ref,z}C_{src}dF^\dagger
  +dF_z C_{src}dF^\dagger,\qquad \Delta J=(K-K^\dagger)/(2i).
\]

The existing KS vertices contract these matrices with the original signed
family multiplicity. Source column weights enter once. The full source
matrix, including coherence, is preserved. The reference is the homogeneous
solution from the identical upstream preparation; no free reference stress
or frozen incoming covariance is accepted. Its history tangent is zero.

The KS constraint adapter selects the independently cached homogeneous
reference from the same preparation. In that case the finite difference
is evaluated as `D + (A_joint - A_cached)`. The small difference between
the two numerical reference integrations is retained, so the stable
contraction remains the same raw current-minus-reference assembly; no
measured drift is discarded.

The returned derivative is the full retarded derivative of the evolved
matter term. The tests compare both algebraic forms at a finite nonzero
history, check the derivative, and check exact observable cancellation of a
common constant column phase. The latter is an algebraic example, not a
claim that the phase is a radius-history solution.

This change avoids numerical cancellation; it supplies no UV subtraction,
source error certificate, new term in the action or physical PASS. All
spectral coverage and field error requirements remain explicit.
