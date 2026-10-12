# A geometry-matched initializer for the same affine vacuum

At the true exact-profile horizon distance `delta=1e-10`, the new numerical
ratio is

$$
w(\delta)=\sqrt\delta\,(w_0+w_1\delta),\qquad
P_{\rm match}=\frac{(1,w)^T(1,w)^\dagger}{1+|w|^2}.
$$

This is a newly constructed rankone vacuum approximation. No archived
column is normalized or modified. Its limiting projector is the same
`diag(1,0)` as the existing affine-horizon vacuum. The
[certificate](../results/development/nsc-incoming-matched-horizon-initializer.json)
encloses the group32, both-sign, `E in [24,32]` mathematical endpoint error.
Physical propagation and the actual source integral remain subsequent work.

## Decision and bounded reuse

Reuse the exact fixed profile, existing horizon enclosure, the
[first omitted profile coefficient](nsc-incoming-horizon-endpoint-bound.md),
and the scalar/projector defect identity from the vacuum-tail owner. The
missing connection is a sufficiently accurate finite initializer for the
same affine vacuum. Include only the first `delta^(3/2)` profile correction,
then enclose its remainder once. Stop when the mathematical endpoint lapse
bound fits the contextual remaining `2.10969e-11` budget or report OPEN.
Do not propagate a field, integrate a source, repeat the previous loose
endpoint bound or revise a physical parameter during this step.

The material risks are losing the geometric indicial exponent, mistaking a
new normalized approximation for normalization of old columns, and claiming
a uniform arithmetic result from a few point controls. The owner makes each
distinction explicit.

## Exact-profile expansion and remainder

Let `B(delta)=W(q_h-delta)/delta` for the unchanged exact profile, and use
`kappa_g=W_rho(rho_h)/2` only for the geometric indicial exponent. With

$$
B=b_0+b_1\delta+\delta^2 B_2,\qquad
K=\frac{-mr+i\lambda}{\sqrt B}=K_0+K_1\delta+\delta^2 K_2,\qquad
D=E/B=D_0+D_1\delta+\delta^2D_2,
$$

the two coefficients are

$$
w_0=\frac{-iK_0}{1/2+2iD_0},\qquad
w_1=\frac{-iK_1-2iD_1w_0+i\bar K_0w_0^2}{3/2+2iD_0}.
$$

Exact algebra cancels the first two orders of the declared Riccati equation

$$
\delta w_\delta=-i\sqrt\delta K-2iDw+i\sqrt\delta\bar K w^2.
$$

The existing fixed profile supplies interval bounds on its second and third
derivatives over the entire tiny horizon collar. The integral Taylor form of
`B=W/delta` avoids division by zero: `B'=integral_0^1 t W'' dt` and
`B''=-integral_0^1 t² W''' dt`. Interval Taylor bounds for `r` and `B^(-1/2)`
give the directed `K2` enclosure. An exact reciprocal remainder identity
gives

$$
D_2=\frac{E(b_1^2-b_0B_2+b_1B_2\delta)}{b_0^2B}.
$$

These are whole-collar enclosures, not sampled residual maxima. The owner
proves `B>0` and obtains `abs(Riccati residual)<=C delta^(5/2)` uniformly in
energy and each actual angular sign. The imported projector identity and
unitary residual estimate yield

$$
\|P_{\rm affine}(\delta_0)-P_{\rm match}(\delta_0)\|
\le\frac25 C\delta_0^{5/2}.
$$

The following incoming-energy/lapse conversion is conditional on exact
subsequent unitary transport. It does not certify an actual numerical ODE.
Multiplicity and both angular signs enter once.

## Source identity and numerical representation

The exact-profile coordinate removes the old stored-horizon coordinate
bridge. Geometric `kappa_g` sets the field equation's indicial exponent;
the binary source occupation `surface_gravity` is recorded unchanged and is
not inserted into a new thermal law. For this vacuum pilot, the limiting
diagonal projector is independent of that occupation parameter. No new
state definition is needed. This statement does not extend to a low-energy
coherent horizon pair without its required phase/sewing correspondence.

`represent_initializer` returns a directed ratio enclosure at a supplied
point energy, a binary ratio and a bound for its normalized mathematical
projector. The exact distance between two such projectors is at most the
ratio distance, so this includes parameter and binary representation error
without assuming a floating matrix normalization is exact. Six point
certificates demonstrate the API; their largest value is **not** promoted
to a uniform continuum arithmetic bound. A future validated propagator can
request the directed initial data at each needed energy.

Archived modes, sources, initializers and their receipts remain intact.
Physical thermal, ODE and source quadrature errors stay separate. No physical
IV, constraint root, stationarity result or metric timestep is claimed.

```sh
python3 scripts/derive_nsc_incoming_matched_horizon_initializer.py --write
python3 scripts/derive_nsc_incoming_matched_horizon_initializer.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_matched_horizon_initializer.py
```
