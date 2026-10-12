# Conditional incoming geometry on the same Cauchy surface

The [incoming constraint gate](nsc-incoming-constraint-gate.md) excludes the
compact pulse's frozen incoming jets as a physical solution. The next
executable domain keeps the **same** intrinsic KS surface at $\rho=1$, with
$N=1$, $\beta=0$, $a^2=3\pi/2-4$, $r^2=2$, and the same normal/spin frame.
It frees the normal derivatives of $a$ and $r$ rather than copying a tensor
onto a different radius or tilted slice.

For each of $f=a,r$, the domain admits

$$
\delta(\partial_T^n\partial_z^p f)|_\Sigma,
\qquad n\ge1,\quad n+p\le4.
$$

There are ten entries per field. Every purely intrinsic/spatial derivative
is unchanged; the baseline normal derivatives use the existing reference
chart $d\rho/dT=-a_0$. Derivative values enter the common Taylor algebra as
coefficients divided by $n!p!$, so mixed derivatives commute. The origin
$T=0$ labels this same incoming surface and is not an elapsed time.
These are off-shell Taylor slots, not twenty additional physical degrees
of freedom; the constraints, evolution equations and gauge determine their
admissible relations.

The [owner](../src/recursive_horizons/nsc_incoming_cauchy_jets.py) exposes the
changed entries, $k_\parallel=n(a)/a$, $k_\perp=n(r)/r$, and their spatial
derivatives. It supplies no preferred nonzero values. Holding lapse and shift
fixed here is a background gauge choice; their constraint equations remain
independent variations and must still be imposed.

## State identification and reference check

Because the intrinsic metric, normal and canonical frame are unchanged,
the abstract canonical density operator $C_0$ lives on the same Hilbert
space. Holding it fixed is an explicit **conditional initial-data
experiment**, not an invented inter-slice isometry. This does not establish
that the altered normal jets have the old smooth Bronnikov past or reproduce
its global horizon preparation.

The fourth-order reference must be recomputed for the changed jets. The
[record](../results/development/nsc-incoming-cauchy-data.json) applies the
existing reference owner to a labeled derivative probe in the LLL, angular,
compact-only and mixed sectors. It compares both momentum signs at
$|k|=16,32,64,128$ and all four raw metric vertices. Initial $H$, pointwise
vertices and $P_0$ remain exactly unchanged. The LLL projector is constant;
its separately owned geometric allocation is not duplicated here.

The sampled massive reference-vertex differences decrease in this check.
This is a necessary high-momentum screen, not a proof of the complete
Hadamard condition, integrated stress convergence, or global preparation
compatibility. The band/boundary remainder and local-action change still
belong to the same source calculation.

The next physical operation is the **incoming lapse and shift constraint
balance** on this domain, using the same-slice energy/current source and the
recomputed reference/local allocation. It must determine compatible data or
reject this domain. The probe's normal-derivative values are not such data;
no metric timestep, duration, new state, parameter refit or physical
stationarity is claimed.

```sh
python3 scripts/derive_nsc_incoming_cauchy_data.py --check
```
