# Source-preserving finite-energy KS propagator

`nsc_ks_energy_propagator.py` reuses the existing joint reference/difference
Dirac evolution on two basis columns per energy. The resulting operator can
be applied to several disjoint source batches without repeating their field
evolution. The physical incoming law remains
`C_Sigma[g] = U_g C_up U_g†`.

For the canonical envelope, let `W_g(E,z)` evolve from the identity matrix
on the fixed upstream slice. Linearity gives

\[
 X_g(E,z)=W_g(E,z)A_{\rm up}(E),\qquad
 F_g(E,z)=e^{-iEz}W_g(E,z)A_{\rm up}(E).
\]

The identity is an operator basis. The original physical `A_up`, `C_src`,
source energies and quadrature weights are retained on application. Their
values are never interpolated. All three coherent columns at each original
energy are reconstructed, and their covariance is contracted by the existing
matter owner. The baseline and channel multiplicities still belong to the
constraint accumulator.

`evolve_energy_propagator(interval, degree, family, grid, target_z, mass,
angular, rho_up, **solver_options)` returns a `KSEnergyPropagator` containing
the reference, evolved difference, its axial derivative, and retarded history
tangents at the operator's energy nodes. The existing evolution is used without
changing its generator. Its fingerprint binds the sampled history and mesh;
the operator digest additionally binds all operator values, target positions,
energy nodes and channel parameters.

`operator.apply(source, A_up)` interpolates those operators and constructs a
`KSDifferenceIncoming` bound to that exact source and initial field. It forms

\[
 F_z=e^{-iEz}(W_z-iEW)A_{\rm up},\qquad
 \delta F_z=e^{-iEz}(\delta W_z-iE\delta W)A_{\rm up}.
\]

The input energy remains a label. The `W_z` terms are retained. Retarded
variations use the same evolution and fixed `A_up`; no source derivative is
introduced. Positive and negative energies can be handled separately by the
existing signed-state map.

The numerical interpolant passes through the actual rounded Chebyshev nodes.
Barycentric weights are calculated from those nodes. Extrapolation outside
the declared finite energy interval is rejected. Ideal Chebyshev remainder,
node evolution error, node displacement, interpolation arithmetic and source
preparation errors are separate requirements. Recovery of unphased tangents
from the existing field interface is also floating-point arithmetic and is
included in the unresolved reconstruction error. No formal polynomial order
or requested ODE tolerance supplies those bounds.

Focused synthetic tests compare all four fields and both matter insertions,
including the full retarded derivatives, against independent direct-column
evolution. They retain a non-diagonal source, exercise original interval
endpoints and detect dropped coherence. Different upstream columns reuse the
same operator while yielding different physical states. These are numerical
implementation controls. Full finite-source accuracy, the changed-history
UV remainder and the physical local incoming gate remain OPEN.
