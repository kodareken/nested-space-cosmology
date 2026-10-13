# RN/PG wall, current domain and local scales

The paired-angular NSC-8 boundary response carries all four spinor components
at every implicit-midpoint trace. Its PG current and allowed wall domain
must be derived in that representation. NSC-3 uses a staggered node/edge
regulator, so its domain cannot be substituted for this one. A finite
spatial boundary solve does not replace the Lorentzian trapped-region
initial-value problem; its inverse coincident Green trace is a port kernel.

For the projected magnetic Dirac sector, QED2 bosonization gives

$$
\mu^2(r)=\frac{|q|g_4^2}{4\pi^2r^2},\qquad
V_{\rm charge}(r)=A(r)\mu^2(r).
$$

The exterior metric and coupling are inputs. This collective charge-mode
mass introduces neither a fermion rest mass nor an independently added
scalar, and it does not complete a self-sourced four-dimensional reduction.

Recorded orthonormal pressures define the crossing locator

$$
B(r)=\frac{|p_r(r)|}{|p_r(r)|+|p_\Omega(r)|}.
$$

$B=1/2$ locates equal pressure magnitudes. It does not establish equilibrium,
a child boundary, exterior decoupling or a regular centre. The curvature
comparison uses $R_4$ reconstructed from the same recorded PG state; the
analytic source-free `Ricci2` reference is excluded from that pressure ratio.
Inner-horizon reach and extrapolated pressure slopes do not settle centre
regularity.

The two local ends use the same state and grid:

$$
\omega_{\rm zero}=\frac{N}{r},\qquad
\omega_{\rm cutoff}=\frac{N+|\beta|}{\Delta r},\qquad
Q=\frac{\omega_{\rm zero}}{\omega_{\rm cutoff}}
=\frac{N\Delta r}{r(N+|\beta|)}.
$$

This is the implemented per-room ratio, measured from finite recorded
metrics; it is invariant under common rescaling of $r$ and $\Delta r$ at
fixed dimensionless lapse and shift. It does not derive Planck's constant
or an absolute zero. The geometry map describes one space as its free
parameter changes and tests the ordering of its features. That map is a
description, not a new geometry discovery. Original diagnostics are at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).
