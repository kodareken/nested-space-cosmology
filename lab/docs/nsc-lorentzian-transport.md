# Physical Dirac propagation through the fixed horizon geometry

The previous boundary maps described the intrinsic spatial Dirac operator.
This calculation implements the actual coordinate-time Hamiltonian on the
project's unwarped four-dimensional black-universe benchmark. Its derivation
is independently checked in [the tetrad calculation](nsc-dirac-tetrad.md).

\[
H_\kappa=-i(\sigma_2-\beta I)\partial_\rho
+\frac{i}{2}\beta'I+\frac{\kappa}{r}\sigma_1,
\qquad j=u^\dagger(\sigma_2-\beta I)u.
\]

The shift changes both the operator and the current. The spin connection fixes
the beta-prime term. Omitting it violates the local number-conservation law.
The state is a compact Cauchy packet in one massless Dirac angular sector;
no vacuum, particle mass or gravitational stationary state is inferred.

## The correct domain changes the calculation

The characteristic coordinate speeds are `1-beta` and `-1-beta`. At the throat
they are about -1.1708 and -3.1708: both move toward the child side, even when
the local normal-frame light direction is outward. The horizon is at
`rho=1.90069160547`.

At the left trapped truncation both fields are outflow. At the right parent
truncation only the minus field is incoming, and its data are set to zero.
There is no reflecting condition at the inner cut. A finite reflecting
Hermitian matrix would change this physical boundary problem. The complete
radial line instead admits norm-preserving evolution through complete
characteristic half-density flows and a bounded angular potential.

## Numerical intervention and independent proof

The [runner](../scripts/check_nsc_lorentzian.py) uses a summation-by-parts
split-form spatial operator and an inflow penalty. The exact semidiscrete
Green identity gives a norm budget consisting of physical outward flux and
the separately recorded inflow-penalty debit. Time evolution and flux
integration use the same RK4 stages. The [record](../results/nsc-4-lorentzian-transport.json)
compares every recorded quantity, including sampled complex fields.

A packet initially supported inside the horizon at positive rho is evolved
through the throat. At the finest 3,201-point grid and time 1.2, roughly
0.94575 of its initial normalized number remains on the child side and
0.05425 has left through the inner outflow boundary. Number and throat-transfer
budgets agree to rounding-level accuracy. These are conserved Dirac-norm
budgets, not a calculation of the complete gravitating energy ledger.

Independent controls include exact characteristic transport with the angular
potential removed, decreasing-grid profile differences, a separate packet
leaving the inner boundary, and causal support bounded by the two null flows.
The formal kappa=0 control is not a spinor-sphere eigenvalue. Numerical
dispersion outside the characteristic support decreases under refinement;
it is not asserted to vanish at finite grid spacing. Profile convergence for
the narrow packet is pre-asymptotic, so no rigorous continuum error bar is
assigned from the three grids.

An independent time refinement at fixed 801-point space resolution uses CFL
0.45, 0.225 and 0.1125. Successive field L² differences are approximately
8.01e-10 and 5.00e-11, giving fourth-order time convergence. The record reports
that order to two decimal places because it is estimated from tiny field
differences. The appreciably larger spatial errors remain explicit.

## What advances, and what remains open

This supplies a causal, conserved-number fermionic propagation sector on the
same explicit horizon geometry, without applying elliptic heat methods to
its nonelliptic spatial Hamiltonian. The metric remains an imported fixed
background. Its source, quantum state, retarded metric response and common
covariant stationary action are still required for physical closure.
The generator of coordinate time need not have positive Killing energy
inside the trapped region; that is not a ghost or a negative-norm inference.

Reproduce: `python3 scripts/check_nsc_lorentzian.py --check`.
