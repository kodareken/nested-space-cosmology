# Boundary stiffness is not a mass gap

The [threshold calculation](../scripts/check_nsc_threshold_response.py) gives
an explicit counterexample to assigning a constant particle mass from the
zero-energy Dirichlet-to-Neumann jump alone. Its [record](../results/nsc-3-threshold-response.json)
uses the same isolated spatial radial throat whose essential spectrum is
[proved to be the whole real line](nsc-radial-spectrum.md).

For kappa=1 let `w=1/sqrt(1+rho²)` and `S=asinh(rho)`. The upper component of
`D²` is `A†A=-d²+w²-w'`. Solve its zero-energy equation on each half with
`u(0)=1` and `u(+-R)=0`. Reduction of order gives the solution proportional to
`exp(-S)` times an integral of `exp(2S)`.

Writing `t=asinh(R)`, the two positive integrals are

\[
I_p=e^{3t}/6+e^t/2-2/3,\qquad
I_c=2/3-e^{-t}/2-e^{-3t}/6.
\]

The parent is rho>0 and its outward normal at the cut is minus partial_rho.
The child is rho<0 and its outward normal is plus partial_rho. Therefore

\[
N_p=1+I_p^{-1},\qquad N_c=-1+I_c^{-1},\qquad
\lim_{R\to\infty}(N_p+N_c)=\frac32.
\]

The units are inverse throat length. Independent second-order ODE integration
checks the maps at R=4,8,16,32. This positive threshold jump coexists with a
gapless full radial Dirac spectrum. It is a local boundary stiffness, not a
particle mass of 3/2. To find a mass one must solve the full spectral pole
problem with its residues and physical state identification.

For nonzero energy the first-order ratio m(E)=v/u determines the covariant
radial derivative through `(d+w)u/u=E*m(E)`. At zero, m can be singular; the
product must be evaluated by its limit or the squared-operator boundary
problem. Substituting m(0) as a finite constant would miss this distinction.

The result does not determine the full compact warped, interacting or
recursive response. It is a reproducible project-specific clarification
using established boundary-value mathematics, not a new physical law.

Reproduce: `python3 scripts/check_nsc_threshold_response.py --check`.
