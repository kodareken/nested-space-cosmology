# Certified profile reduction for residual validation

The enclosed axial Fourier table extends to index 8192. A smaller
convolution can represent exactly the same physical profile when every
discarded coefficient is added to the error budget. If index K is retained,
the derivative-j remainder is increased by

\[
\sum_{K<|k|\le K_{old}} |2\pi k/L|^j\,|c_k|_{upper}.
\]

The prior continuum tail beyond K_old remains present. This uses the
certified coefficient balls and therefore does not fit a decay curve,
drop a physical frequency, or change the profile. All derivatives needed
by the field residual are accounted for.

The performance control repeats the already owned continuous residual-cell
calculation with retained profile index 1024. It reuses the same compute
code object in a private namespace whose convolution wrapper applies this
bound-preserving reduction. Process globals are unchanged. All time and
radius remainders, source weights and history bindings remain in the proof.
The original cell remains an immutable comparison record. No field or
source equation is evolved in this control.

This improves the cost of a complete history validation. It does not by
itself turn a one-cell certificate into a full-history or physical gate
certificate.
