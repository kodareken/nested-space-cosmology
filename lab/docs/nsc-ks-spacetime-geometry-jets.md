# Mixed geometry jets for the exact finite-history defect

The exact Weyl-defect estimate needs higher mixed derivatives of the existing
history. `spacetime_geometry_jets` supplies interval Taylor coefficients in
actual canonical coordinate time T and axial coordinate z, through total
order9 by default. It reuses the owned background, plateau, Chebyshev-profile
and analytic-history identity implementations.

The metric remains

\[
N=1,\quad\beta_K=0,\quad a=a_{\rm ref}(T),\qquad
r=r_{\rm ref}(T)+\chi(s)[s\,w(z)+s^3U(z)/6],\quad s=T-T_\Sigma.
\]

Every normal and mixed derivative comes from these same functions. The
module introduces no independent higher jets, compensators, fitted scales
or new state data.

The background owner supplies rho-derivatives. To obtain T-derivatives,
solve the finite Taylor recurrence `rho_T=-a(rho)` and compose the existing
rho-jets with its increment. The increment `rho(T0+tau)-rho(T0)` has exact
zero constant even when the common basepoint ranges over a real ball. This
structural zero is imposed before composition; subtracting two independent
interval copies of the basepoint would lose that dependency.

The normal profile is then evaluated directly in T: `s=s0+tau`, where
`s0=-integral_1^rho d rho/a`. The two scalar time series `chi(s)*s` and
`chi(s)*s^3/6` multiply the existing axial series for w and U. Coefficients
are derivatives divided by `t! z!`. The reference radius and axial scale
have no axial derivatives. All arithmetic uses directed balls.

A rho/z input box encloses the derivatives at every basepoint in that box.
A box that cannot be classified against a plateau boundary requests
subdivision. Radius positivity must be proved by its interval value; no
positive floor is inserted. The owned preparation slab remains
`rho in [1,33/32]`.

Tests check the coordinate-time chain rule, the incoming identities
`delta r_T=w` and `delta r_TTT=U` with mixed derivative factorials,
independent high-precision differentiation of the nonzero normal-window
transition, and containment of mixed derivatives by a geometry box. A
nonzero `d_T^4 d_z^5 r` is retained in the transition. No field or source
evolution is run.

These are geometric coefficient enclosures. They do not yet evaluate the
higher projector derivatives, the exact Weyl-defect seminorm or the physical
incoming constraints. The local gate remains OPEN.
