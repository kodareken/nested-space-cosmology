# One local history from two coherent axial functions

The local incoming search uses one actual smooth radius history,

$$
\delta r(s,z)=\chi(s)\left[s\,w(z)+s^3U(z)/6\right],
\qquad s=T-T_\Sigma.
$$

`LocalIncomingFamily` supplies w and U as Chebyshev series multiplied by the
same C-infinity plateau. Its interval is I=S(1)+[.12,.18], where the plateau
equals one; its support is inside S(1)+[.09,.21]. The normal cutoff retains
the existing .007/.03 radii. These are coordinate domains, not elapsed
cosmological time. The coefficient arrays have shape(2,n). The base owner accepts n=8,16,32,64.
n=128 is the same (w, U) class through `LocalIncomingFamily128` in the
high-mode owner; the base allowlist is unchanged so earlier source hashes
stay valid. Amplitude order is all w coefficients followed by all U
coefficients.

Every derivative through order four is the derivative of that same function.
The metric adapter reuses `CompatibleRadiusDirection` and
`CompatibleIncomingMetric`; it supplies the six mutually compatible incoming
slots w,w_z,w_zz,w_zzz,U,U_z. The geometry never accepts six independent
pointwise jets as a surrogate for the surface functions.

On either transition interval, with u=(abs(z-center)-inner)/(outer-inner),
the window is 1/(1+exp(1/(1-u)-1/u)). Both logistic factors are evaluated
independently to retain derivative accuracy near the inner edge. The value
is exactly one/zero in the constant regions, where all higher derivatives
are zero. This is the same smooth plateau construction as the normal window.

Since abs(T_j)<=1 on the supported Chebyshev coordinate, the coefficient
one-norms W and V give a conservative global radius bound

$$
r_g\ge 1-s_{\max}W-s_{\max}^3V/6.
$$

The implementation cushions the nonnegative sums/products upward and the
final subtraction downward in binary64. A nonpositive bound signals that
this sufficient positivity check fails; it does not prove negative radius.
This bound says nothing about matter accuracy or constraint closure.

Four focused controls compare window derivatives to independent60-digit
differentiation, check exact support and interior polynomial derivatives,
check coherent surface slots/amplitude derivatives and expose invalid inputs
or an unsafe radius bound. They do not evolve a field or alter source data.

The returned geometry description keeps the physical gate OPEN. Both N,beta
still require the actual C_Sigma[g] evolved from the fixed source and a full
error enclosure on I. No state, stress, physical duration or root is supplied
by this geometry object; the rho0 transmitting seam stays fixed.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_local_incoming_family.py
```
