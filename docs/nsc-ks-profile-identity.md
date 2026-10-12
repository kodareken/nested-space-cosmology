# Analytic profile identity across numerical resolutions

Profile identities serialize the exact binary amplitudes, all w/U
Chebyshev coefficients, centres, axial windows and coordinate maps. The
history identity also includes normal windows. Zero-amplitude directions
remain present because their retarded variations need the same basis.

Numerical meshes are absent from this identity: changing resolution does
not change an analytic history. A different profile with the same amplitudes
does change the identity. The axial-only identity permits reuse of a
geometric Fourier table across normal-window changes, while the complete
state/history binding must still reject such a change.

This identity accompanies the existing source/preparation and sampled
geometry bindings; it does not replace their checks or certify a field.
