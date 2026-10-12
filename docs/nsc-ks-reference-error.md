# Reference error from a causally unaffected point

The homogeneous reference A_ref is evolved jointly with the difference D.
Its numerical error must be included when bounding a matter difference.
It is not removed as a measured stress drift.

Choose the numerical period's origin z0, outside the compact axial support.
On the owned slab, `a>=4/5`, so the total axial characteristic distance is
at most `(rho_up-1)/(4/5)^2`. If both distances from z0 to the neighbouring
periodic support copies exceed this bound, the exact inhomogeneous envelope
there is exactly the homogeneous reference solution. This is finite
propagation of the same Dirac operator, not a periodic physical-state claim.

If the completed field proof gives a uniform row error e0, then

\[
\|A_{ref,true}-A_{ref,num}\|_{row}
\le e_0+\|D_{num}(z_0)\|_{row}.
\]

The numerical difference at z0 is retained on the right-hand side. It is
neither forced to zero nor subtracted from the physical field. Multiplication
by sqrt(2) gives a weighted Frobenius field bound. The homogeneous reference
has only the source carrier, so its axial error is at most E_max times that
field bound. Both current and reference contributions can then enter the
same coherent matter-error inequality.

The initial source approximation and spectral coverage remain separate
requirements. This argument supplies a numerical reference-error bound
only after the full-history, whole-period field-error bound is complete.
