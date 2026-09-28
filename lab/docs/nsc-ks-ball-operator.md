# Time and spatial coefficients of the same radius history

`AnalyticRadiusFamily` restricts validation to the approved common normal
window with owned `LocalAxialFunction` profiles. Amplitudes, binary
Chebyshev maps and cutoffs are unchanged. Generic callbacks cannot be
identified from sample agreement and are rejected by this proof path.

For each rho cell the background and normal-window jet owners enclose
`1/a^2`, `1/a`, `1/(a*r_ref)` and the coefficients multiplying
`w(z)^p U(z)^q` in the reciprocal-radius expansion. All scalar coefficients
come from that geometry. The field residual does not accept a new physical
potential or source term.

A Taylor polynomial of degree d at the cell centre is accompanied by

\[
|f-P_d|\le
\sup_{\rho\in J}\frac{|f^{(d+1)}(\rho)|}{(d+1)!}
\left(\frac{|J|}{2}\right)^{d+1}.
\]

The derivative coefficient is enclosed over the whole cell by interval
jets, including the flat-cutoff boundary treatment. Polynomial coefficients
are converted to the actual step fraction. A cell requiring subdivision
does not receive a fabricated remainder. At rho=1 the enclosing chart ball
uses exact dyadic radius input; the default radius constructor can round
outward across the restricted chart boundary.

Spatial product jets are formed from the same w,U functions. These enclose
ordinary derivatives, with factorials restored from the normalized Taylor
coefficients. The time Taylor remainder and reciprocal-radius series
remainder are distinct contributions and both must enter the final
continuous residual bound.
