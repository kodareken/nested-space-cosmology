# Finite-collar vacuum endpoint: the exact-profile coordinate matters

For group32, both angular signs and `24 <= E <= 32`, the existing finite
initializer has the directed mathematical error estimate

$$
\|C_{\rm frame}(\delta_{\rm num})-P_{\rm affine}(q_{h,\rm num}-\delta_{\rm num})\|
\le D_{\rm transport}+D_{\rm coordinate}.
$$

The [receipt](../results/development/nsc-incoming-horizon-endpoint-bound.json)
keeps these terms, their incoming-lapse projection and finite-series norm
defect separate. A broad upper bound is not a measured source error or a
non-existence result. No archived mode, source value or initial matrix changes.

## Decision and stopping condition

Reuse `PairedHorizonSeedMap.working_frame`, its ten-term `horizon_frame`, the
exact fixed profile, the existing horizon enclosure and the unitary
variation-of-constants estimate already used in the vacuum-tail owner.
The missing connection is the finite-offset initial vacuum versus its
affine-horizon limit. Evaluate one uniform endpoint enclosure on group32
`[24,32]` before any validated field propagation. Stop with this receipt,
PASS or OPEN against the contextual remaining lapse budget `2.10969e-11`;
other LOW errors are still outside that budget. No radial or source solve
is part of this check.

The three material risks are confusing the stored horizon coordinate with
the exact profile root, conflating the geometry slope with the frozen source
occupation parameter, and interpreting a finite-difference frame residual
as a uniform error bound. The proof treats them separately.

## Same finite formula, phase-free projector residual

Write the phase-free Frobenius column as `v=(a,b)` before its constant unitary
angular rotation, with `z=mu sqrt(2 delta/kappa_source)` and
`nu=E/kappa_source`. Its common tortoise phase cancels exactly in `v v†`.
The finite series obey

$$
\partial_y a+izb/2=izb_{\rm last}/2,\qquad
\partial_y b+iza/2+i\nu b=0,\qquad y=\log\delta.
$$

The only truncation residual is order `delta^10`. Its integrated covariance
residual also bounds the finite-series norm defect; it is counted once.
No column or covariance is normalized to enforce a desired result.

The exact profile has `W=2 kappa_geometry delta (1+epsilon)`, where a directed
bound on its second derivative gives `abs(epsilon)<=J delta`. Compare the
exact generator with the leading Rindler generator used by the frame. A
transverse generator difference is bounded using `2 norm(v)^2`; its diagonal
part uses the sharper `norm([sigma3,v v†])=2 abs(a b)`. Since `b=O(sqrt(delta))`,
the remaining profile terms integrate as powers `delta^(1/2)` and
`delta^(3/2)`, not a bound on the rapidly oscillating full-frame phase.

The old `interior_W` evaluator retains a fourth-order polynomial near the
horizon. The receipt encloses the exact profile's fifth-order Taylor
remainder, but never substitutes that approximate evaluator for the exact
generator. Its stored-coordinate and coefficient rounding differences are
not declared zero, nor is the old ODE arithmetic certified by this estimate.

## Coordinate and parameter accounting

The comparison point is the original mathematical chart point
`q=q_h_stored-delta_stored`. Its distance from the exact-profile horizon is
`delta_true=delta_stored+q_h_exact-q_h_stored`. The existing root enclosure
gives `q_h_exact-q_h_stored` about `1.49026e-15`, which is material compared
with the numerical collar `delta_stored≈1e-10`.

The coordinate bridge is bounded by the phase-free projector derivative on
the small interval between these two distances. Geometry uses the derivative
of the unchanged exact profile. The source and initializer retain the frozen
binary `surface_gravity` value. Replacing it by the enclosed geometric slope
inside a source occupation would change the supplied numerical source and
is not performed. Exact defining masses/angular labels and their stored
values are enclosed together in the geometric comparison.

The bound concerns the **mathematical finite initializer formula**. The
sampled comparison with the actual floating `working_frame` is a convention
and arithmetic indicator only. Uniform floating-initializer and subsequent
`DOP853` propagation errors remain explicit `null` fields. The future energy
projection uses preservation of operator distance under **exact subsequent
unitary transport** and the generic two-dimensional vertex bound. It does
not certify the archived numerical ODE or source quadrature and supplies no
evolved field.

## Next constructive connection

An aligned chart removes the coordinate bridge but does not by itself make
the present absolute-residual estimate sharp. The first omitted geometric
coefficient is the `delta^(3/2)` term of the vacuum component ratio. Write

$$
W=2\kappa_g\delta+w_2\delta^2+\cdots,\quad
K(\delta)=\frac{-mr(\delta)+i\lambda}{\sqrt{W/\delta}}=K_0+K_1\delta+\cdots,
\quad D(\delta)=E\delta/W=D_0+D_1\delta+\cdots.
$$

For `b/a=sqrt(delta)(w0+w1 delta+...)`, the already declared Dirac generator
gives

$$
w_0=\frac{-iK_0}{1/2+iE/\kappa_g},\qquad
w_1=\frac{-iK_1-2iD_1w_0+i\bar K_0 w_0^2}{3/2+iE/\kappa_g}.
$$

The leading Rindler frame contains the nonlinear last term but omits the
`K1,D1` profile corrections. Including those corrections in a consistent
exact-profile coordinate is the next representation to certify, rather
than repeating this loose residual estimate. The geometric indicial exponent
uses `kappa_g`, while the source occupation parameter stays frozen.

For this high-energy **vacuum-only pilot**, the affine limiting projector is
the same diagonal vacuum irrespective of the stored thermal kappa. A
geometry-matched frame can therefore improve its numerical representation
without changing that vacuum state. This does not authorize an analogous
phase change in the genuinely low coherent source: its horizon-pair phases,
sewing and source covariance must be matched together. Positive thermal
remainders stay with their separate owner.

```sh
python3 scripts/derive_nsc_incoming_horizon_endpoint_bound.py --write
python3 scripts/derive_nsc_incoming_horizon_endpoint_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_horizon_endpoint_bound.py
```
