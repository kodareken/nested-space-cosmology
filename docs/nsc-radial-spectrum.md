# What the smooth throat does and does not supply

The [exact calculation](../scripts/check_nsc_radial_spectrum.py) and
[record](../results/nsc-3-radial-spectrum.json) test the unwarped spatial radial
Dirac operator already declared by the project:

\[
D=-i\sigma_2\partial_\rho+\sigma_1\frac{\kappa}{\sqrt{1+\rho^2}},
\quad \mathcal H=L^2(\mathbb R,d\rho)\otimes\mathbb C^2,
\quad\operatorname{Dom}D=H^1(\mathbb R)\otimes\mathbb C^2.
\]

The Hermitian potential is bounded, so this is a bounded perturbation of the
self-adjoint free operator. A constructive Weyl sequence gives a direct
answer about its mass gap, independently of a finite-box eigenvalue search.

Take `f(s)=sqrt(15/16)(1-s²)` on `[-1,1]`, zero outside, and
`v_plus=(1,i)/sqrt(2)`. This bump is in H¹ and has norm one;
`||f'||²=5/2`. For any real energy E, define

\[
\psi_n(\rho)=n^{-1/2}f((\rho-n^2)/n)e^{iE\rho}v_+.
\]

These normalized states escape to infinity, converge weakly to zero, and obey

\[
\|(D-E)\psi_n\|
\le\frac{\sqrt{5/2}}{n}
+\frac{|\kappa|}{\sqrt{1+(n^2-n)^2}}\longrightarrow0.
\]

Thus every real E is in the essential spectrum. In particular, the complete
radial throat operator has **no nonzero spectral mass gap**. This statement is
stronger than observing decreasing finite-box levels and does not depend on
which finite-box boundary condition was used in the earlier scan.

The zero-energy equations are also explicit:
`u=C exp(-kappa asinh(rho))`, `v=C exp(+kappa asinh(rho))`. Each nonzero solution
grows at one end (or is constant for kappa=0), so no global L² zero mode exists.
Zero belongs to the continuum; it is not a normalizable particle state.
The record does not assert a complete exclusion of embedded eigenvalues.

This clarifies the missing physical connection. A smooth throat supports an
energy-dependent scattering/boundary response. It does not, by gluing alone,
supply the constant nonzero Phi in the elementary two-sheet massive symbol.
That Phi requires an actual dynamical internal/boundary sector and its
stationary equations. The result does not cover the full compact warped,
interacting, recursive or Lorentzian operator, and it does not reject the
nested-space architecture. The Weyl-sequence method is established spectral
theory; this is its application to the project's chosen radial realization.

Reproduce: `python3 scripts/check_nsc_radial_spectrum.py --check`.
