# Directed mixed-normal shift-response coefficient

On the unchanged incoming surface, the retained action gives

$$
\boxed{16.575975622058497\le c_v\le16.575975622058518},
\qquad \mathcal E_\beta=S_\beta+c_v\,\delta r_{Tz}.
$$

This encloses a fixed action response, not an adjustable source/coupling or
selected physical initial datum. The physical source cancels in this
coefficient calculation.

## Decision and stopping condition

Reuse the established incoming plane `u=delta r_TTT`, `v=delta r_Tz`, the
fourth-order reference projector with its spatial Weyl product, and the full
two-dimensional local action Euler operator. The fixed physical source
cancels from these response coefficients. The missing result is a directed
nonzero enclosure of `c_v` in `E_beta=S_beta+c_v*v`; the closed `c_u` proof
is not rerun.

First verify, rather than assume, the proposed first-spatial-gradient trace
identity from the existing `star_order` convention. If valid, extract only
the reference orders needed for the linear `v` response, integrate the
coefficient functions using standard full-line moments, and combine with
the exact local Euler coefficient. Compare against authenticated stored
reference nodes and local channel responses as controls. Stop after one
directed coefficient certificate, or a precise obstruction. Do not rerun
probe/reference generators, fit coefficients, or select physical initial data.
The quadratic lapse coefficient is outside this bounded task unless it
comes without a separate derivation.

Material risks are the sign/factor of the Weyl trace, omission of mixed
coordinate Euler derivatives, and promoting stored float response samples
to an exact integral certificate. Keep the full spatial convention and
vary beta before imposing the background gauge; use the stored data only
as independent checks. All geometry, action and source parameters remain
the existing declared inputs.

## Verified first-spatial-gradient trace

For the homogeneous formal projector `Q=(I+b.sigma)/2`, write the first
spatial correction as `R=(t I+x.sigma)/2`. The actual Weyl convention gives

$$
(Q\star Q)_{(1)}=-\frac\epsilon4(b_z\times b_k)\cdot\sigma,
\qquad
QR+RQ-R=\frac12[(b\cdot x)I+t\,b\cdot\sigma].
$$

Star-idempotency and formal `b.b=1` therefore give
`Tr R = epsilon*b.(b_z cross b_k)/2`. Exact Pauli algebra verifies the sign
and factor; a direct control uses the original `star_order` function.
The inherited Bloch recurrence verifies normalization through the required
order. On `r(T,z)=r_background(T)+v*T*z`, second and higher spatial gradients
are at least quadratic in v. Hence this first-gradient identity supplies
the complete linear-v shift response through reference order four.

The raw shift vertex is `-k I`. Its reference order3 term need not vanish
pointwise: it is odd in momentum and integrates to zero on the full line.
Order1 vanishes; orders2 and4 remain. No symmetry of the physical C0 state
is imposed.

## Closed reference coefficient

Define `L=ell/r`, `M2=m^2+L^2`, `Ha=a_T/a`, `Hr=r_T/r`,
`A2=a_TT/a`, `R2=r_TT/r`, and

$$S=A2+R2+2H_r^2-5H_aH_r+2H_a^2.$$

With `d=copy_count*degeneracy`, including the actual angular signs once,
the full-line integral with `dk=a dp` gives

$$
c_{v,2,g}=\frac{d a L^2}{12\pi r M2},
$$

$$
c_{v,4,g}=\frac{d a L^2}{120\pi r}
\left[\frac{S}{M2^2}
+\frac{2m^2(H_r^2+4R2)}{M2^3}
-\frac{20m^4H_r^2}{M2^4}\right].
$$

These follow by coefficient extraction from the existing homogeneous Bloch
recursion plus the verified Weyl trace. The elementary full-line moments are
`2/(3 M^2)` for `p^2/omega^5`, and `32/(315 M^8)`, `16/(315 M^6)`,
`4/(63 M^4)` for `p^2,p^4,p^6` over `omega^11`. They are evaluated as exact
beta/gamma-function identities, not fitted to saved nodes. Zero-angular
reference responses, including the separately owned LLL, are zero.
The retained sums are approximately `19.12120967609074` and
`0.2531553470423731`.

## Full two-dimensional local Euler contribution

The exact beta response is extracted from the original spherical curvature
contractions and light/compact densities with N=1 and a=a(T), retaining
general radius spatial jets until the Euler derivatives are applied:

$$
\mathcal E_\beta=L_\beta-\partial_TL_{\beta_T}
-\partial_zL_{\beta_z}+\partial_T\partial_zL_{\beta_{Tz}}.
$$

The mixed derivative appears once. Only after this operation are spatial
jets set to the established plane. Exact algebra yields:

| Existing channel | Contribution to c_v |
|---|---|
| Einstein | `-16*pi*A*a*r` |
| Compact Weyl | `32*pi*C_Weyl*a*r*S/3` |
| Barred radial R2 squared | `-h_q*a*r*S/(30*pi)` |
| LLL geometry | `a*q/(12*pi*r)` |
| WZ Weyl | `-a*r*(S*log(r)-R2+3*Hr^2+1/r^2-2*Ha*Hr+A2)/(30*pi)` |
| WZ Euler | `11*a*r*(Hr^2+1/r^2)/(180*pi)` |
| WZ BoxR | `a*r*(R2-2*Hr^2+Ha*Hr+A2-2*Ha^2)/(60*pi)` |

The constant Euler bulk, Maxwell, cylinder and WZ gauge contributions are
exactly zero. This is not the homogeneous shift-zero calculation: mixed
coordinate Euler derivatives are essential. New algebra-only density
derivative controls compare the formulas with the unchanged original
`local_action_densities` implementation; no old contour/probe generator runs.

## Inputs, checks and scope

Directed hulls include the exact baseline geometric/harmonic definitions and
stored evaluations. `A` and `C_Weyl` remain the declared binary ledger values;
no uncertainty for an unprovided alternative exact ledger is invented.
Stored n24 reference nodes, formal-order integrals, local channel responses
and the independent mixed control verify the new formulas without preparing
new reference nodes. The closed c_u proof is not called.

This result establishes `c_v != 0` on the existing retained incoming plane.
Together with the earlier positive c_u interval, it removes the coefficient
nonsingularity condition in that plane's conditional relation. Source
finiteness/accuracy, a surface solution, parent preparation, endpoints and
extended stationarity remain separate. c_v2 is not newly certified here.

```sh
python3 scripts/derive_nsc_incoming_shift_coefficient.py --prepare
python3 scripts/derive_nsc_incoming_shift_coefficient.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_shift_coefficient.py
```
