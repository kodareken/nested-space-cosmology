# Conservative scalar evolution on the retained wedge

The [scalar forcing derivation](nsc-atlas-scalar-forcing.md) supplies the
equation, source pullback and wall lift. The
[profile-moment successor](nsc-atlas-profile-moment.md) supplies economical
boundary queries. This note specifies their numerical connection to the
first exterior scalar correction, with the same source, action and finite
reference geometry. It introduces no finite amplitude or preparation law.

## Fixed coordinates and the moving wall

Let $y$ be the areal scalar correction pulled back to the original null
coordinates, $g_1(V)$ its source pullback, and $B_U=\mathfrak b$ its wall
lift. Set $z=y-g_1-B$. On the original wedge,

$$
z_V=F(U,V)-F(U_0,V)
-\int_{U_0}^{U}\{K_0[z+g_1+B]+G\}\,dU.
$$

The source and lower boundary have $z=0$, and $z_U|_W=0$. Keep a fixed
$U$ grid. Clip inactive nodes to $U_W(V)$ and give them the terminal
value of $z$; evaluate $B$ at that clipped coordinate too. The material
rate of this terminal state equals its partial $V$ rate because its
wall derivative is zero. No wall-advection term belongs in the $z$ ODE.

The reconstructed wall signal has a different derivative:

$$
\frac{dy_W}{dV}=g_1'(V)+z_V|_W+U_W'(V)\mathfrak b(U_W(V)).
$$

The last term is the wall pullback of the lift. Omitting it loses the
material derivative; inserting it into the $z$ equation counts it twice.
With normalized time $t=(V-V_a)/\Delta V$, multiply the entire ODE
right-hand side by $\Delta V$ exactly once. Radial integration still uses
the original $U$ measure.

## Cache the known forcing

At each cache time, use
$\eta=(U-U_0)/(U_W(V)-U_0)$ to place the requested forcing samples on the
retained radial interval. The existing cumulative evaluator integrates
across those samples and its original field cells. Its final requested
grid entry may be set to $U_b$; the evaluator clips that entry to the
actual wall. All other requested entries already lie inside the wall.

Interpolate $F(\eta,t)$ and $G(\eta,t)$ separately within the source-time
blocks $[0,1/4]$, $[1/4,3/4]$, $[3/4,1]$. The initial implementation uses
cubic interpolation in time and linear interpolation in $\eta$. Interpolate
$g_1$ with its supplied first derivative, converting that derivative to
normalized time for the Hermite fit. Preserve its post-pulse tail.

At an ODE query, map the original $U$ coordinates to $\eta$ and evaluate
these known coefficient interpolants. The evolved unknown remains in
$(U,V)$: this composition introduces no advection or Jacobian in the PDE.
The same source value of $F$ supplies both the array endpoint and the
subtraction. Check $F_W=0$ and the quiet corner directly. Evaluate $K_0$
from the reference geometry rather than interpolating it in time.

Training times include both endpoints and source edges. Withheld slices
stay outside every fit. Include radial midpoints on withheld slices:
radial-knot discrepancies assess time interpolation, while midpoint
discrepancies assess the combined time/radial interpolation. Compare
$F,G,g_1$, the cell-width-scaled source
slope, and the source/wall identities. Time interpolation, radial cache
interpolation, leading-field resolution and metric quadrature are distinct
numerical approximations. Their refinements must remain separately visible.

## Evolution and independent controls

Use the existing clipped Volterra topology and cumulative trapezoidal
integration, with DOP853 for the normalized-time ODE. Split the evolution
at source edges, grid/wall births, cache knots and lift knots. Share the
preparation's process-CPU and scalar-kernel budget throughout; cache queries
inside the ODE perform no new physical field or metric integration.

The independent software controls include a nonconstant source and lift
with zero bulk forcing, and a manufactured solution with a nonlinear
moving wall, nonzero source flux, potential and density. Construct the
manufactured density from the differentiated exact field so that the
radial integrand has an independently integrable polynomial form. These
checks resolve source-subtraction signs, normalized-time scaling, newly
uncovered grid cells and the total wall derivative. They do not establish
accuracy for the physical input.

Numerical comparisons must distinguish cache-time resolution, cache-radial
resolution, correction-grid resolution and ODE tolerance. The existing
conditional operator estimate with $Q=0.13882099234357653$ may propagate
genuine uniform input bounds. Substituting sampled discrepancies gives a
sensitivity indicator only, not a uniform theorem or a detector-derivative
bound. The next cap reconstruction must still include its forcing, regular
reflection and observer-clock pullback.

## Implementation

The [forcing cache](../src/recursive_horizons/nsc_atlas_forcing_cache.py)
and its [controls](../tests/test_nsc_atlas_forcing_cache.py) own coefficient
interpolation. The [scalar solver](../src/recursive_horizons/nsc_atlas_scalar_evolution.py)
and its [controls](../tests/test_nsc_atlas_scalar_evolution.py) own fixed-$U$
evolution and the wall readout. The [admission](prior-art-and-open-claim.md)
and [current result](current-result.md) record the finite calculation and
its actual outcome.
