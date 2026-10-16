# What the calculations show

Current presentation, 16 October 2026: [the living atlas manuscript](../paper/nested-space-cosmology.md) and [overview](../paper/figures/geometry-overview.png). Dated entries below preserve their original numerical domains and chronology. Former current-paper labels describe those historical editions; [the catalog](../paper/catalog.json) owns publication status.

## First nonlinear exterior rectangle, 16 October 2026

The [four-grid record](../results/development/nsc-atlas-nonlinear-rectangle-v1/rectangle.json)
now propagates the exact-model source across the original exterior rectangle.
It compares zero amplitude and `alpha=0.00012983127142551443` on 16 and 32
cells per direction, with the original compact pulse and quiet bottom.
All 2,560 cells complete in 5.71875 process CPU seconds and 23,237 charged
evaluations. Each cell uses two full Newton updates; no tube change,
subdivision or alternative seed is needed. The field before the first source
front is exactly zero.

The zero-amplitude numerical geometry agrees with analytic Schwarzschild:
at 32 cells, the maximum radius error is `3.4444e-9` and the maximum
`log(C)` error is `5.9628e-7`, both about four times smaller than at 16.
The following comparisons use the same 225 interior points on both grids:

| Quantity | 16-cell value | 32-cell value at common points |
| --- | ---: | ---: |
| Maximum raw finite-amplitude `abs(E_U)` | `1.11663e-5` | `2.79078e-6` |
| Maximum raw finite-amplitude `abs(E_x)` | `0.2074424` | `0.0591010` |

The finite-amplitude state differences at those points are `1.02824e-8`
in radius, `1.37207e-6` in `log(C)` and `6.50796e-5` in the normalized
scalar. Comparing the finite-minus-zero response between grids gives
`1.79285e-9`, `2.69272e-7` and `6.01643e-7`, respectively. These paired
differences retain the raw residuals above; they do not subtract away a
constraint error. Source transverse scalar derivatives retain the nonzero
post-pulse value, with the maximum finite-amplitude defect falling from
`7.29062e-6` to `1.61650e-6` over each grid's post-pulse source points.

The [march and diagnostics](../src/recursive_horizons/nsc_atlas_nonlinear_rectangle.py)
and [creation-only recorder](../scripts/record_nsc_atlas_nonlinear_rectangle.py)
pass 31 focused implementation/authentication tests and independent review.
Authentication of the actual saved result passes with zero source
integrations or cell solves. The [nonlinear owner](nsc-atlas-finite-amplitude.md#original-exterior-rectangle-and-independent-diagnostics)
specifies the equations and independently differenced diagnostics.

The saved source identifies the main derivative issue: its finite-difference
`E_x` maximum is `0.208647` at 16 cells and `0.0594446` at 32, even though
the source's analytic constraint is below `1.06e-15`. These nearly reproduce
the propagated residuals. The sampled transport difference
`r*E_x(Ubar,x)-r_source*E_x(0,x)` is much smaller: `3.75545e-4` and
`1.14297e-4` at common interior points. The steep known source is therefore
an important part of the raw differentiation error. This diagnosis does
not remove the remaining propagated discretization error.

The next connection uses the exact source plus an interpolated propagated
increment to supply values and controlled first jets on `Ubar=Ubar_b`.
That leg feeds new exterior cells adjoining the moving wall; it is not
itself the wall. The regular polar cap and hybrid interface then supply
the coupled continuation. This result is a numerical exterior region
with measured discretization changes, not a full nested spacetime,
formation/stability result or uniform finite-amplitude remainder bound.

## Exact source initialization and characteristic cells, 16 October 2026

The [source integrator](../src/recursive_horizons/nsc_atlas_nonlinear_source.py)
now constructs the same model's finite-amplitude initial leg and all its
first jets, with the original quiet bottom and corner. Continuous polynomial
profile bounds select the initialization control
$\alpha=1.2983127142551443\times10^{-4}$, half the source-only continuation
threshold. Tightening both ODE tolerances changes the sampled Einstein
radius by at most $1.44\times10^{-14}$ and geometric mass by
$1.71\times10^{-13}$, using 3,063 shared reference/profile evaluations.
The [38 source controls](../tests/test_nsc_atlas_nonlinear_source.py) pass in
1.71875 process CPU seconds. They include independent vacuum field
integrals, quadratic source moments, frame identities and corner matching.

The [cell solver](../src/recursive_horizons/nsc_atlas_characteristic_cell.py)
now advances the nonlinear mixed equations through one implicit-midpoint
characteristic cell. Its analytic Jacobian, fixed error scales, compact
state bounds and finite iteration limits make convergence or failure
explicit. The [35 controls](../tests/test_nsc_atlas_characteristic_cell.py)
pass in 0.59375 process CPU seconds. Exact-corner Schwarzschild and cap
tests approach the derived fourth-order local error under balanced cell
halving. This is compatible with second-order accuracy for a complete
smooth grid; no regional order has yet been measured.

These implementations preserve the source and governing equations. The
source threshold applies to initialization, and its algebraic constraint
identities are not an ODE error enclosure. The next calculation assembles
these cells into the original prepared-source rectangle and measures its
propagated constraints and refinement. The moving wall, regular cap and
full observer signal remain downstream of that regional step.

## Regional equations and mixed-chart closure, 16 October 2026

The [normalized bulk operator](../src/recursive_horizons/nsc_atlas_nonlinear_bulk.py)
now evaluates the exact mixed metric and scalar equations with a continuous
zero-amplitude limit. It keeps the potential increment separately and
distinguishes mass rates on the null-constraint surface from derivatives of
arbitrary input jets. Independent canonical-field differentiation and exact
vacuum checks include the trapped Schwarzschild region used by this model.

The [regular cap derivation](nsc-atlas-finite-amplitude.md#regular-polar-areal-cap)
supplies mass and lapse through radial constraints, with Jordan central
proper time fixing the lapse constant. Radius-squared variables and a
regular mass-density integral remove explicit division by the axis radius;
the inherited pulse fronts still require their stated piecewise regularity.
The same owner now specifies every source first jet through the exact
source ODE and three mixed transport equations, with the original quiet
corner fixing their constants.

A [mixed-chart wall construction](nsc-atlas-finite-amplitude.md#closing-a-polar-cap-against-the-null-exterior)
joins that polar cap to the null exterior. Angular matching solves the
velocity, and transmission supplies exactly one incoming cap scalar
characteristic. The cap wall lapse comes from the interior, so it is not
independently prescribed a second time. On the branch `A_cap>0,H>0`, the
exterior acceleration equation and cap constraints imply the remaining
wall laws and incoming exterior constraint. This implication has two
independent derivations and stays regular at the zero of `beta_plus`.

The [hybrid boundary helper](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
implements this map. The combined bulk, hybrid and existing boundary/clock
controls pass 102 tests in 0.734375 process CPU seconds, with independent
mathematical and code review. These developments specify a coherent coupled
system and its pointwise operators. Exact source integration, the characteristic
interior update, regular cap march and moving interface must still be
assembled into a trajectory. No finite-amplitude spacetime or uniform
remainder estimate is claimed by these algebraic controls.

## Exact boundary and central-clock helpers, 16 October 2026

The [implemented boundary operator](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
now converts outgoing radial, scalar and lapse derivatives into the six
incoming derivatives and wall state rates. It imposes the same action's
angular and temporal jumps, scalar transmission and normal force. Its
normalized field formulation remains regular at zero amplitude, zero
radial speed and the reference zero of `beta_plus`.

A companion centre helper converts regular null-coordinate data to Jordan
proper-clock rate and the normalized central scalar. Both outputs respect
null-coordinate rescaling, and the clock is anchored by the caller at the
same quiet event. The [independent controls](../tests/test_nsc_atlas_nonlinear_boundary.py)
pass 32 tests, including canonical-action identities at nonzero amplitude,
vacuum recovery and observer conversion. No field or spacetime evolution
is run by these helpers.

Compatible regional lapse evolution and a regular central solution remain
inputs to connect. This implements the boundary algebra needed for the
coupled solver; it does not yet provide a nonlinear trajectory or its
amplitude/error interval.

## Reference scalar bounds and constraint transport, 16 October 2026

The [nonlinear owner](nsc-atlas-finite-amplitude.md#a-bounded-vacuum-reference-collar)
now supplies explicit reference-geometry bounds on a normal collar of
half-width $\ell/20$ over $0.45\le R/\ell\le0.55$. Its lapse stays
above $2/3$, and its one-sided coefficient derivatives have rational
bounds. These close the scalar energy inequality with amplification
coefficient $119/8$, retaining the outer-boundary flux and the prescribed
forcing norm. Homogeneous value and normal-derivative transmission cancel
the internal wall flux. Independent derivations check the bounds and signs;
no field evolution or fitting enters this result.

Commuting once with wall time closes a second energy inequality, with
coefficient $1549/8$ and both forcing norms explicit. A separate regular
cap strip gives an exact incoming-characteristic formula for the central
field error and a first-commuted energy trace bound. These estimates use
compatible initial and outer data; the outer traces are not new sources.

Direct differentiation also gives
$\partial_V(rE_U)=\partial_U(rE_V)=0$ for the null Einstein-constraint
residuals. It proves their propagation in each smooth region. At the wall,
the full junction laws transfer the tangential residual, while an
independent normal-force condition removes the remaining incoming
constraint mode. The earlier temporal-junction shortcut already assumes
the constraints and cannot supply their boundary closure circularly.

A noncircular boundary choice now imposes the incoming exterior constraint,
temporal junction and normal force independently. It propagates initial
angular matching and forces the cap's incoming constraint to vanish. This
is a consistency implication for a regular solution. Constructing that
coupled solution map and bounding its metric, field and residual remain
the next connection to a nonzero nonlinear amplitude interval.

## Central scalar coefficient, 16 October 2026

The [central reconstruction](nsc-atlas-central-scalar-response.md) now
connects the exterior correction to the interior observer. It uses the
exterior wall derivative, current clock pullback, outgoing cap response
and reflected incoming forcing. The
[first comparison](../results/development/nsc-atlas-central-scalar-v1/central.json)
holds ray order 16 fixed and varies incoming radial order 16/32 at five
arrival labels strictly inside the finite window.

The fixed-central-clock field is
$\ell\phi_C/\varepsilon=p_0+\alpha p_{1,C}+O(\alpha^2)$.
Write $T=\tau_C/\ell$ for dimensionless central proper time. At radial
and ray orders both 32 the retained coefficients are:

| Central time $T$ | $p_0$ | $p_{1,C}$ |
| --- | ---: | ---: |
| -0.54363418 | 0 | $3.49\times10^{-137}$ |
| -0.46546305 | 42.21372894 | 3395.219 |
| -0.42507962 | -0.265835912 | -34413.869 |
| -0.38358521 | -40.14663419 | 4369.585 |
| -0.29640937 | 0.00127167605 | 1454.355 |

The maximum radial change in $p_{1,C}$ is $1.09781\times10^{-5}$.
Propagating the existing exterior cache/grid/ODE comparisons through the
same cap algebra changes it by at most $1.15477\times10^{-4}$,
$1.66666\times10^{-4}$ and $5.95797\times10^{-8}$, respectively. These
are separate sampled sensitivities, not a total uncertainty enclosure.
The independent [incoming-ray 16/32 comparison](../results/development/nsc-atlas-central-ray-v1/ray.json)
changes the central coefficient by at most $5.28497\times10^{-5}$,
the optical time integral by $1.33819\times10^{-6}$ and the reflected
forcing integral by $3.60495\times10^{-7}$. It reuses the saved coarse
case and computes five new rays with one leading-field preparation,
taking 50.046875 process-CPU seconds and 26,350,673 scalar kernels.

One leading-field preparation and ten incoming rays take 27.9375
process-CPU seconds and 14,498,108 scalar kernels. No wall, cap-clock or
exterior correction is evolved again. At radial-32/ray-16, the travel
integral reproduces the earlier cap result within $1.43\times10^{-14}$;
fixed-material-time arrival agrees within $2.85\times10^{-14}$. Ray order
32 changes each by $5.21\times10^{-10}$. Quiet recovery is below
$3.49\times10^{-137}$. Two positive-radius samples on each ray retain
the computed centre-near potential separately from its analytic zero
at the centre. The frame cancellation closes below $9.10\times10^{-13}$
as an algebraic identity; it is not a second independently evolved solution.

The large cubic coefficient is a result of this preparation and geometry.
In particular, a leading field near zero does not make its absolute
correction small. No finite $\alpha$ or uniform relative-error claim is
selected from this table. The nonlinear amplitude interval and remainder
need their own calculation. The [exact nonlinear owner](nsc-atlas-finite-amplitude.md)
derives the source constraint, smooth bulk equations and transformed wall
laws from the same action. The invariant reduced field and its normal
derivative transmit continuously. On a branch with positive Einstein-radius
expansion, the temporal wall condition follows from angular matching and
the scalar jump without dividing by the vanishing exterior normal-radius
factor. A source-only continuation bound keeps radius,
slope and conformal factor positive for sufficiently small nonzero
amplitude. The coupled bulk/interface/centre remainder remains unresolved.

## Exterior scalar correction, 16 October 2026

The [first propagated correction](../results/development/nsc-atlas-scalar-evolution-v1/evolution.json)
evolves the conservative scalar equation with the original prescribed
reduced-field source and the resolved 32-panel clock lift. It restores the
saved forcing and boundary inputs, with no new field, wall or metric solve.
Four evolutions separate cache-time resolution, correction-grid resolution
and ODE tolerance. Each retains 17 wall readout cuts.

| Independent comparison | Maximum wall-value change | Maximum material-derivative change | Value change / sampled peak |
| --- | ---: | ---: | ---: |
| Forcing cache, 17/33 times; 129 correction nodes | $1.75501\times10^{-3}$ | $1.32631\times10^{-4}$ | 0.654856% |
| Correction grid, 65/129 nodes; 33 cache times | $6.98421\times10^{-6}$ | $9.31679\times10^{-5}$ | 0.00260606% |
| ODE tolerance, $2\times10^{-6}/2\times10^{-7}$; fine cache/grid | $3.72805\times10^{-9}$ | $3.95592\times10^{-8}$ | 0.00000139107% |

The finer sampled peak is $0.2679995064$ in the dimensionless exterior
areal-gauge coefficient $y_W$. Its terminal value is $0.2194607242$;
its terminal material derivative with respect to $V$ is $-0.5827560530$.
These multiply $\alpha$ in $\chi/\varepsilon=c_0+\alpha c_1+O(\alpha^2)$;
they are not finite-amplitude field predictions by themselves.

The complete calculation takes 5.9375 process-CPU seconds and 950,926
scalar kernels. The four RHS counts are 1,874, 1,987, 1,536 and 2,677,
all within their unchanged 6,000-call ceilings. Source-value, source-rate
and wall-flux defects vanish, as do the saved/reference geometry
differences. The tighter ODE comparison changes the wall value by only
$2.13\times10^{-6}$ of the cache comparison's change.

The retained interior arrays allow the discrete wall equation and material
derivative to be checked without another evolution. The comparisons measure
finite numerical sensitivity, not a uniform error enclosure. The central
observer still requires the incoming cap forcing, regular reflection and
clock pullback; its analytic reconstruction is the next connection.

## Clock-resolved scalar boundary, 16 October 2026

The [clock comparison](../results/development/nsc-atlas-lift-clock-v1/clock.json)
resolves the interior interpolation plateau identified below. It holds the
leading field, geometry, 41 boundary knots and 40 withheld queries fixed,
reuses the 81 authenticated outgoing cap jets, and changes only the wall
clock integration from 16 to 32 panels. Both calculations retain the
continuum clock rates as their supplied derivatives.

| Diagnostic | 16 clock panels | 32 clock panels |
| --- | ---: | ---: |
| Maximum withheld boundary-value discrepancy | $1.73136\times10^{-6}$ | $5.74819\times10^{-8}$ |
| Maximum cell-width-scaled slope discrepancy | $2.16504\times10^{-6}$ | $4.08960\times10^{-8}$ |
| Maximum clock $V_{1,s}$ polynomial/continuum gap | $8.86470\times10^{-4}$ | $6.29210\times10^{-5}$ |

The finer boundary passes the unchanged $2.2795554\times10^{-7}$ sampled
target. The maximum lift difference is $3.24913\times10^{-10}$, and the
clock derivative-gap identity closes below $1.92\times10^{-12}$. The pair
takes 35.234375 process-CPU seconds and 5,002,056 scalar kernels, including
one leading-field preparation and two wall solves. It uses no new outgoing
cap quadrature. This establishes a sufficiently resolved sampled boundary
input for the first exterior propagation; it is not a uniform error bound.

The first saved-input scalar propagation attempt stopped at its 6,000-RHS
ceiling in the first case, before producing a solution record. A no-solve
diagnosis found 2,609 distinct crossings of radial cache knots. Splitting
at every crossing would require at least 41,760 RHS calls before adaptive
refinement. The completed comparison above therefore keeps the original
cache and splits, uses tolerances appropriate to the measured input
accuracy, and measures the effect of a tighter tolerance separately.

## Exterior forcing cache, 16 October 2026

The [conservative evolution implementation](nsc-atlas-scalar-evolution.md)
now caches the known exterior forcing in normalized radius and evaluates
it back at the original null coordinates. The
[17/33-time comparison](../results/development/nsc-atlas-forcing-cache-v1/cache.json)
shares 65 direct slices: 33 training times on 129 radial knots and 32
withheld times on those knots plus their 128 radial midpoints. Both fits
retain the same field, geometry and radial quadrature.

| Maximum withheld discrepancy | 17 training times | 33 training times |
| --- | ---: | ---: |
| $F$, at radial training knots | 0.0655104 | 0.00700171 |
| $G$, at radial training knots | 7.91622 | 0.0908529 |
| $F$, at radial midpoints | 0.0624395 | 0.00862133 |
| $G$, at radial midpoints | 7.89536 | 0.0893962 |
| Source pullback $g_1$ | 0.000120579 | 0.0000141481 |
| Time-cell-width-scaled source slope | 0.000656558 | 0.00000352091 |

The knot comparisons isolate time interpolation at those coordinates;
midpoint comparisons measure combined time and radial interpolation.
The complete pair takes 22.84375 process-CPU seconds and 4,431,777 scalar
kernels. The 49 repeated slice requests are reused. Direct source value
and slope discrepancies remain below $4.56\times10^{-14}$ and
$1.78\times10^{-12}$; wall flux is exactly zero. The independent wall-mass
discrepancy remains $6.50\times10^{-7}$, separate from cache interpolation.

These inputs refine substantially, but sampled discrepancies are not
uniform bounds. Their effect on the scalar wall signal must be compared
using the same correction grid, separately from grid and ODE refinement.
The new scalar solver passes independent moving-wall manufactured
solutions, source subtraction and material-derivative controls. No physical
scalar propagation is claimed by this input record.

## Profile moment and local lift, 16 October 2026

The [pre-arrival moment](nsc-atlas-profile-moment.md) now evaluates the cap
mass from a cumulative first-profile-derivative integral. The
[comparison record](../results/development/nsc-atlas-profile-moment-v1/moment.json)
holds the prepared field, wall and geometry fixed, comparing moment Gauss
orders 4/8 against the original radial calculation. It includes nine
boundary cuts and 24 independent mass points on the outgoing cap rays.

| Quantity | Maximum moment 4/8 difference | Maximum moment-8 / radial difference |
| --- | ---: | ---: |
| Cap Einstein mass coefficient | $1.07\times10^{-14}$ | $2.50\times10^{-9}$ |
| Outgoing scalar coefficient $J_-$ | $8.89\times10^{-16}$ | $5.87\times10^{-11}$ |
| Material derivative $(J_-)_s$ | $5.69\times10^{-14}$ | $4.94\times10^{-9}$ |
| Wall boundary value $\mathfrak b$ | $2.67\times10^{-15}$ | $1.52\times10^{-10}$ |
| Wall slope $\mathfrak b_U$ | $1.82\times10^{-12}$ | $7.18\times10^{-8}$ |

The whole comparison takes 12.25 process-CPU seconds and 2,637,474 scalar
kernels, including one field and wall preparation and the original radial
reference. The new integration is stable under its quadrature refinement.
The larger radial/moment difference remains a separate numerical
representation discrepancy: the dense field and its continuum first jet
need not be exact derivatives of one floating interpolant, and the original
radial reference also has quadrature error. These nine cuts do not identify
each contribution separately or establish a uniform bound.

The maximum boundary-value discrepancy is over 1,500 times smaller than the
existing $2.2795554\times10^{-7}$ interpolation target. Even multiplying the
slope discrepancy by the entire retained wall-$U$ width gives less than
$1.8\times10^{-9}$. This supports using moment order 8 for the admitted local
lift comparison. Its five-node, fixed-per-cell Hermite construction and
exact antiderivative pass independent polynomial and continuity tests.

The [local comparison](../results/development/nsc-atlas-local-lift-v1/lift.json)
then shares 81 direct jets between 25/41 training knots and 40 withheld
queries. It takes 17.9375 process-CPU seconds and 3,737,156 scalar kernels.

| Diagnostic | 16 active cells | 32 active cells |
| --- | ---: | ---: |
| Maximum withheld boundary-value discrepancy | $3.95548\times10^{-5}$ | $2.90055\times10^{-5}$ |
| Maximum cell-width-scaled slope discrepancy | $7.95539\times10^{-5}$ | $9.98697\times10^{-5}$ |
| Terminal lift | 0.0117052349145 | 0.0117052398354 |

The maximum sampled lift difference is $2.47058\times10^{-8}$, while
boundary values still miss the original target. Endpoint value jumps are
below $5.33\times10^{-15}$; exact polynomial differentiation and wall
transport agree below $1.78\times10^{-15}$ and $1.52\times10^{-16}$,
respectively. Thus assembly roundoff is not the observed limitation.
The post-pulse value error is below $2.99\times10^{-10}$.

The largest fine-fit value discrepancies occur at $t=0.5078125$ and
$0.4921875$, through the pulse interior. The near-edge discrepancy has
substantially improved, but further knot refinement alone does not explain
this interior plateau. The inherited wall clock uses eight panels and
four-node polynomial antiderivatives for its values, while the supplied
clock slopes use the continuum rate. The earlier nine-cut comparisons lie
at these panel boundaries. A saved-jet decomposition into cap backscatter,
material radius and clock pullback is the next diagnostic before changing
the wall-time resolution. The unresolved boundary is not yet used for
physical scalar propagation. The independent forcing-cache implementation
and manufactured evolution controls can proceed meanwhile.

The [saved-jet decomposition](../results/development/nsc-atlas-lift-components-v1.json)
isolates the fine-grid plateau: maximum value
errors are $5.99\times10^{-10}$ for backscatter, $1.67\times10^{-8}$ for
material radius and $2.901\times10^{-5}$ for the clock. Their direct values
sum to within $3.56\times10^{-15}$ and their fitted errors to within
$7.44\times10^{-15}$. The first two components now meet the original
target; the clock does not. Its authenticated recomputation takes 0.546875
process-CPU seconds and no physical solve. Refining the wall clock's panels while retaining
its continuum derivative is the next physical comparison.

## Scalar lift and exterior forcing, 16 October 2026

The [scalar-input calculation](nsc-atlas-scalar-forcing.md) now evaluates
the boundary's analytic first derivative and the conservative exterior
forcing. Both use the same prepared field, quiet anchors and computed wall.
The comparison holds 129 field nodes and 257 geometry knots fixed, retaining
cap radial/ray order 16 and the existing wall settings.

The [9/17-knot lift comparison](../results/development/nsc-atlas-scalar-inputs-v1/lift/inputs.json)
shares 33 direct boundary queries. Sixteen midpoints remain outside both
fits. The boundary is $\mathfrak b(V)=r_Uc_{1,r}|_W$; its Hermite
antiderivative supplies the lift exactly for that interpolant.

| Interpolation diagnostic | 9 knots | 17 knots |
| --- | ---: | ---: |
| Maximum withheld boundary-value discrepancy | 0.5546024823 | 0.0610876282 |
| Maximum local-width-scaled slope discrepancy | 3.1990105788 | 0.0596985619 |
| Terminal lift | 0.0139211082 | 0.0117209086 |

The maximum lift difference is 0.0022004862. Neither fit meets the
$2.2795554\times10^{-7}$ numerical target, chosen as one tenth of the
sampled 129/257 field sensitivity in $\mathfrak b$, with a $10^{-11}$
floor. This is a diagnostic target, not a uniform mathematical enclosure.
The derivative/potential identities agree within $2.88\times10^{-13}$;
arrival and wall-coordinate rate discrepancies are below
$6.01\times10^{-10}$. The unresolved part is the interpolation accuracy,
so this coarse lift is not yet used for a physical scalar evolution.

The [exterior forcing pair](../results/development/nsc-atlas-scalar-inputs-v1/forcing/inputs.json)
compares radial Gauss orders 2/4 on nine time slices. A single cumulative
mass/lapse integration supplies all 129 requested coordinates per slice.

| Quantity | Maximum 2/4 change |
| --- | ---: |
| Einstein mass coefficient $\mu_E$ | $3.10\times10^{-11}$ |
| Einstein lapse coefficient $\zeta_E$ | $2.04\times10^{-10}$ |
| Conservative scalar flux $F$ | $7.20\times10^{-14}$ |
| Scalar forcing density $G$ | $4.45\times10^{-10}$ |

The independent timelike mass comparison has maximum discrepancy
$6.50\times10^{-7}$ at both radial orders. The smaller quadrature changes
do not remove that separate field/transport discrepancy. Source field and
tangent-derivative discrepancies are below $4.56\times10^{-14}$ and
$1.39\times10^{-12}$, respectively; the leading wall derivative and $F_W$
remain zero. The corrected source pullback persists after the prescribed
leading pulse ends: at the final cut $g_1=-1.3969369088$ and
$g_1'=-10.3608466089$. Setting that correction to zero after the leading
input would change the declared preparation.

The complete lift and forcing pairs take 17.05 and 9.52 process-CPU seconds,
using 4.89 million and 1.20 million scalar kernels respectively, including
their field and wall preparations. Both records bind their producer and
input bytes. The next work is the boundary interpolation, checked time
interpolation of the exterior forcing, and propagation of the scalar
correction. Neither these inputs nor the mass agreement supplies the
central corrected waveform or a finite-amplitude remainder.

An additional analytic identity expresses the pre-arrival cap mass as
$\nu_E=\tfrac12\int_{q_a}^{q}[f'(z)]^2\,dz-(1-r^2)f(q)^2/(4r)$.
It follows from the full radial
constraint and quiet regular centre. The current evaluator still uses the
original radial metric; the identity gives a route to cheaper independent
boundary refinement, with its numerical comparison still to be done.

The [saved-data assessment](../results/development/nsc-atlas-scalar-lift-assessment-v1/assessment.json)
then tests local higher-order Hermite stencils on the same 17 training knots,
retaining all 16 withheld queries and splitting at the source edges. It
requires 0.078125 process-CPU seconds and no new field or wall solve.

| Requested local training nodes | Maximum withheld value error | Maximum width-scaled slope error |
| --- | ---: | ---: |
| 4 | 0.0123173 | 0.0286763 |
| 5 | 0.00956095 | 0.0220658 |
| 7 | 0.0849837 | 0.335182 |
| 9 | 0.765149 | 3.46095 |

All four miss the original target. The post-pulse block, which has five
available training nodes, reaches $3.01\times10^{-10}$ value error. The
active block carries the failure, with peak discrepancies at its near-edge
withheld cuts $t=0.28125$ or $0.71875$. Increasing the degree on these
unchanged knots is therefore not the next refinement. Additional active-edge
samples, after comparing the cheaper mass-moment evaluator with the existing
radial metric, are the concrete next connection. This assessment does not
construct a new global lift or propagate a scalar field.

## Scalar backscatter boundary, 16 October 2026

The [scalar-metric calculation](nsc-atlas-scalar-metric-response.md) now
evaluates a boundary input for the next scalar response. The computed cap
metric has a nonzero scattering potential. On the certified no-return
interval, its outgoing correction at the wall is a known forcing integral
$J_-$, determined from the same leading field and metric. It requires no
additional outgoing waveform. Combining that integral with the existing
wall displacement and exterior clock gives the correction's exterior
areal derivative $c_{1,r}$.

The continuation keeps the original prescribed reduced-field input exact.
At the final reference radius $x=0.55$, with 257 field/geometry nodes and
16-point cap radial and outgoing-ray rules:

| Dimensionless coefficient | Value |
| --- | ---: |
| Outgoing cap derivative $J_-$ | -1.8223975002 |
| Backscatter contribution to $c_{1,r}$ | -2.9736230746 |
| Wall-displacement contribution | -0.5062222284 |
| Exterior areal derivative $c_{1,r}$ | -3.4798453030 |

These are coefficients in $\chi/\varepsilon=c_0+\alpha c_1+O(\alpha^2)$,
where $\alpha=\kappa_4\varepsilon^2/\ell^2$. The actual correction is cubic
in the source amplitude; no finite amplitude has been selected.

| Independent comparison | Maximum change in $J_-$ | Maximum change in $c_{1,r}$ |
| --- | ---: | ---: |
| [Field 129/257](../results/development/nsc-atlas-scalar-boundary-v1/field/boundary.json) | $3.16\times10^{-7}$ | $1.73\times10^{-6}$ |
| [Radial order 8/16](../results/development/nsc-atlas-scalar-boundary-v1/radial/boundary.json) | $2.44\times10^{-7}$ | $3.99\times10^{-7}$ |
| [Ray order 8/16](../results/development/nsc-atlas-scalar-boundary-v1/ray/boundary.json) | $2.00\times10^{-3}$ | $3.26\times10^{-3}$ |
| [Radial order 16/32](../results/development/nsc-atlas-scalar-boundary-v1/radial_fine/boundary.json), field 129 | $7.07\times10^{-11}$ | $1.17\times10^{-10}$ |
| [Ray order 16/32](../results/development/nsc-atlas-scalar-boundary-v1/ray_fine/boundary.json), field 129 | $7.90\times10^{-11}$ | $1.31\times10^{-10}$ |

Each comparison changes one setting at nine common wall cuts. The larger
coarse-ray difference motivated the further refinement. The primary
incoming-label and radius quadratures use independently placed nodes and
agree in $J_-$ within $5.59\times10^{-11}$. General and quiet-profile
potential formulas agree within $2.95\times10^{-13}$. The algebraic wall
flux residual is at most $2.23\times10^{-16}$; this reconstruction check
is not an independent propagated-field error estimate.

At the final cut, the signed forcing integrals are $0.0384206688$ and
$-1.8608181690$. Every queried ray remains in the established cap, with
outgoing no-return margin at least 0.9052343761 and central-time quiet
margin at least 0.2868530625. No central second-derivative query is needed.
The five complete pairs take 19.06--34.61 process-CPU seconds and at most
8.45 million scalar kernels each, including preparation and wall evolution.
All records bind the current sources and are authenticated without replay.

Forty-nine focused tests pass: fourteen earlier operator/coordinate checks,
nine additional flux/lift/clock identities, sixteen manufactured boundary
checks and ten record/budget checks. The next connection is an interpolant
with measured numerical error, followed by the forced exterior and cap
waveforms. The current nine samples do not provide that interpolant, a
finite-amplitude remainder or a physical preparation mechanism. Refinement
differences are sensitivity evidence rather than a rigorous total enclosure.

## Smooth cap and central arrival clock, 15 October 2026

The [cap calculation](nsc-atlas-cap-response.md) reconstructs the smooth
interior metric from the same prepared scalar that drives the exterior
and moving wall. It keeps the reflected field, transformed cap potential
and regular centre. The lapse is normalized so the cap time coordinate
is Jordan proper time at the central observer. No independent cap source,
mass constant or observer-clock history is added.

The radial Einstein mass agrees with the independently integrated timelike
wall mass. The metric and existing wall displacement determine the wall's
cap-time correction and the inward null-arrival correction. Converting
from fixed wall proper time to a fixed exterior input label includes the
already computed exterior embedding shift.

At the final reference event $x=0.55$, the primary setting has 257 field
nodes, 257 geometry knots, 16-point radial/ray rules and 16 cap clock panels:

| Dimensionless quadratic coefficient | Value |
| --- | ---: |
| Cap Einstein excess mass $\nu_E$ | 17.7590647726 |
| Cap Jordan excess mass $\nu_J$ | 17.7590784210 |
| Cap lapse $\delta_E$ | 103.4370285720 |
| Wall cap-time shift $T_1$ | -8.4727761614 |
| Central arrival at fixed proper wall time $\Delta_C$ | -128.0670033577 |
| Central arrival at fixed exterior label | -130.3664151052 |

These are coefficients in expansions such as
$T_C(V)=q_0(V)+\alpha\left.\Delta_C\right|_V+O(\alpha^2)$, with
$\alpha=\kappa_4\varepsilon^2/\ell^2$. They are not finite-amplitude
time shifts. The final conversion between the two clock comparisons is
2.2994117476 in the same coefficient units.

The [8/16 radial comparison](../results/development/nsc-atlas-cap-response-v1/radial/response.json)
reduces the maximum cap mass-path discrepancy from
$1.49\times10^{-4}$ to $2.52\times10^{-9}$. This maximum includes
the cap clock quadrature nodes and nine saved wall cuts. The direct
induced-metric residual falls to $1.33\times10^{-8}$; the derivation
shows that it is the same mass mismatch multiplied by $2/[xA_0(x)]$,
not an additional independent error estimate. The radial comparison
changes the arrival coefficient by at most $6.30\times10^{-3}$.
The [16/32 radial refinement](../results/development/nsc-atlas-cap-response-v1/radial_fine/response.json),
with 129 field nodes fixed, changes it by $1.97\times10^{-9}$.
This checks the interior radial integrals along the arrival rays, where
most of the coarser change arose. Its wall mass-path discrepancy falls
from $2.48\times10^{-9}$ to $1.54\times10^{-9}$.

The [129/257 field comparison](../results/development/nsc-atlas-cap-response-v1/field/response.json)
changes the fixed-label arrival coefficient by at most
$6.42\times10^{-7}$. The [8/16 ray comparison](../results/development/nsc-atlas-cap-response-v1/ray/response.json)
changes it by $2.31\times10^{-5}$; the independent
[16/32 ray refinement](../results/development/nsc-atlas-cap-response-v1/ray_fine/response.json),
with 129 field nodes held fixed, changes it by only $1.84\times10^{-9}$.
Changing the [regular-centre threshold](../results/development/nsc-atlas-cap-response-v1/centre_radius/response.json)
from $10^{-4}$ to $10^{-3}$ changes it by at most $9.53\times10^{-10}$.
The primary small-radius choice reduces derivative-query cost while
retaining agreement in the reported coefficients.

With 129 field nodes fixed, the [129/257 geometry comparison](../results/development/nsc-atlas-cap-response-v1/geometry/response.json)
changes the arrival coefficient by at most $1.06\times10^{-8}$.
The [2/4 derivative-quadrature comparison](../results/development/nsc-atlas-cap-response-v1/derivative/response.json)
leaves the saved coefficients equal at floating-point precision; the
[8/16 centre-quadrature comparison](../results/development/nsc-atlas-cap-response-v1/centre_order/response.json)
changes them by at most $7.11\times10^{-15}$. These are sensitivities of
this finite calculation, not exact continuum identities.

The [8/16 cap-time-panel comparison](../results/development/nsc-atlas-cap-response-v1/time/response.json)
changes the arrival coefficient by at most $1.92\times10^{-5}$. Its
second solution reuses the identical reference-ray integrals while
recomputing the wall-time integral; fresh and reused nodes are reported
separately. This time-panel sensitivity exceeds the field and finer
radial/ray sensitivities, and is retained explicitly in the accuracy account.

Signed quadrature contributions locate the negative ray integral inside
the cap. At the final primary cut, the positive contribution is
$6.47\times10^{-9}$ and the negative contribution is
$-112.0145129776$. Their sum is $-112.0145129711$.
The centre-normalized lapse can outweigh the positive mass term in the
optical integrand. The recorded extrema are sampled values, not bounds
over the continuous ray.

The first five completed records retain their exact source context at
`19983a7b4730d1fdc8e59ae4ba2f3efd560913a6`. The geometry comparison
exposed an overly strict validator: interpolated arrival endpoints differed
by at most $2.22\times10^{-16}$. The corrected validator compares domains
within an absolute $10^{-13}$ tolerance and keeps the actual coordinates.
It does not alter any field, metric equation or saved result. The
[reproduction guide](reproducing.md) names the source context explicitly.
The geometry, derivative and centre-order records retain their corrected
producer at `5e58c9fd2866f8adcc107ee22f13f389e7f15c77`.
The time and further-radial comparisons use the current producer. All ten
complete pairs finished in 35.98--56.31 process-CPU seconds and used fewer
than 24.94 million scalar-kernel evaluations each, within the unchanged
120-second and 30-million budgets. Dense exploratory arrays are not part
of the evidence package.

This connects the formal sourced exterior, moving wall and regular cap
to a central null-clock correction. The [next scalar derivation](nsc-atlas-scalar-metric-response.md)
now specifies its forced wave equation, source continuation and wall
pullbacks, with fourteen independent analytic checks. The cap's outgoing
correction on the original no-return interval is a known forcing integral;
the scalar waveform has not yet been evolved. A nonzero finite-amplitude range, remainder estimate,
formation and stability remain to be established. The refinements here
measure numerical sensitivity; they do not supply those stronger results.

## Moving wall and exterior clock, 15 October 2026

The [wall calculation](nsc-atlas-wall-evolution.md) evaluates the same
prepared field across the full retained exterior wedge. Radial mass and
lapse integration reaches the moving wall. Its Jordan mass is compared
independently with the timelike functional; both one-sided mass laws then
drive the quiet-anchored wall displacement and exterior null clock.

At the final reference radius $x=0.55$, with 257 field nodes, 257 geometry
knots, radial Gauss order 4 and 16 clock panels:

| Dimensionless quadratic coefficient | Value |
| --- | ---: |
| Exterior mass $\mu_+$ | -28.8834223256 |
| Cap excess mass $\nu_-$ | 17.7590784205 |
| Wall displacement $\rho$ | -5.2868390417 |
| Exterior clock ratio $\eta=V_1/V_0$ | -0.8178125677 |

These multiply $\alpha=\kappa_4\varepsilon^2/\ell^2$ in their respective
expansions. The cap excess subtracts the fixed vacuum mass at the actual
radius; its total pulled-back mass coefficient is 15.3601752054 here.
The calculation selects no finite amplitude.

The [129/257 field comparison](../results/development/nsc-atlas-wall-evolution-v1/field/response.json)
reduces the maximum radial/timelike mass discrepancy from
$6.49\times10^{-7}$ to $1.63\times10^{-7}$. These maxima include the
clock quadrature nodes as well as nine reported time cuts. At the common
cuts, the changes in $\mu_+,\nu_-,\rho,\eta$ are at most
$3.38\times10^{-7}$, $2.09\times10^{-7}$, $3.64\times10^{-8}$ and
$3.57\times10^{-9}$ respectively. Angular matching, the two displacement
equations, frame conversion and clock normalization have residuals below
$4.3\times10^{-14}$; the reference wall-spline clock normalization has
its separate $9.85\times10^{-10}$ diagnostic at 257 geometry knots.

The [8/16 clock-panel comparison](../results/development/nsc-atlas-wall-evolution-v1/time/response.json)
changes $\eta$ by at most $2.41\times10^{-8}$, while the
[2/4 radial-rule comparison](../results/development/nsc-atlas-wall-evolution-v1/radial/response.json)
changes it by less than $1.35\times10^{-13}$. Doubling the
[field derivative quadrature](../results/development/nsc-atlas-wall-evolution-v1/derivative/response.json)
leaves the eight recorded comparison quantities identical in this run.
The [geometry interpolation](../results/development/nsc-atlas-wall-evolution-v1/geometry/response.json),
[source quadrature](../results/development/nsc-atlas-wall-evolution-v1/source/response.json)
and [wall ODE tolerance](../results/development/nsc-atlas-wall-evolution-v1/ode/response.json)
comparisons change these quantities by at most $5.64\times10^{-11}$,
$2.85\times10^{-14}$ and $8.25\times10^{-13}$ respectively. Each
complete two-setting comparison used 50.50--99.29 process-CPU seconds
and at most 12,761,511 scalar-kernel evaluations within its fixed budget.
These are measured refinement differences, not a continuum enclosure.

An [independent overlap calculation](../results/development/nsc-atlas-wall-overlap-v1/response.json)
compares the radial method with the earlier mixed metric solver at $U=U_a$
on nine cuts. Their mass coefficients differ by at most $4.08\times10^{-7}$,
consistent with the earlier metric-grid discrepancy scale. The lapse
coefficients differ by at most $1.30\times10^{-6}$. This comparison uses
a different radial quadrature, exposing spatial error that a constraint
computed with the same trapezoidal rule can cancel. The integrated
constraint residual alone is therefore not a total metric error estimate.
This check used 11.24 process-CPU seconds and 4,108,105 kernel values.

The formal boundary now has computed displacement and exterior embedding
coefficients, with no extra adjustable wall datum. The next connection is
the smooth cap metric with its regular centre and fixed central clock,
followed by the corrected arrival map and next scalar response. A
finite-amplitude remainder, formation and stability remain to be calculated.
The dated entries below retain the scope of their original results.

## Smooth exterior metric response, 15 October 2026

The [bulk calculation](nsc-atlas-bulk-response.md) evolves the radius and
lapse coefficients from the original prepared scalar and gravitational
initial data on $[U_0,U_a]\times[V_a,V_b]$. It retains the original null
coordinates and fixed Jordan source lapse, including the source's conformal
and affine-parameter correction. Both scalar derivatives come from the
field's integral equations; first metric derivatives are evolved alongside
the metric.

The [33/65 metric comparison](../results/development/nsc-atlas-bulk-response-v1/metric/response.json)
reduces the largest sampled disagreement between geometric mass and its
independent flux expressions from $1.65\times10^{-6}$ to
$4.31\times10^{-7}$. The two integrated null-constraint defects decrease
from $1.98\times10^{-6}$ and $1.01\times10^{-6}$ to
$4.94\times10^{-7}$ and $2.61\times10^{-7}$. The areal-null lapse
diagnostic is the same first constraint in different variables.

At $(U_a,V_b)$ the 65-node calculation with 257 field nodes gives:

| Dimensionless quadratic coefficient | Value |
| --- | ---: |
| Jordan radius $h_J$ | -5.7131513827 |
| Jordan logarithmic lapse $c_J$ | -1.5461215575 |
| Jordan geometric mass $\mu_J$ | -28.8833375991 |

On common sampled interior nodes and nine time cuts, the metric-grid
differences in these three quantities are at most $1.15\times10^{-7}$,
$3.82\times10^{-6}$ and $1.23\times10^{-6}$. The
[129/257 field comparison](../results/development/nsc-atlas-bulk-response-v1/field/response.json)
changes them by at most $2.14\times10^{-8}$, $9.55\times10^{-8}$ and
$5.78\times10^{-7}$. Doubling the
[derivative quadrature](../results/development/nsc-atlas-bulk-response-v1/derivative/response.json)
changes them by less than $7.2\times10^{-15}$. The separate
[geometry interpolation](../results/development/nsc-atlas-bulk-response-v1/geometry/response.json)
and [source quadrature](../results/development/nsc-atlas-bulk-response-v1/source/response.json)
comparisons change these three coefficients by less than
$7.2\times10^{-15}$ and $1.5\times10^{-13}$ respectively. Tightening
the [time tolerance](../results/development/nsc-atlas-bulk-response-v1/time/response.json)
tenfold changes them by less than $1.21\times10^{-11}$. All six
comparisons finished within their individual budgets (17.16--21.91
process-CPU seconds; at most 10,448,594 scalar-kernel evaluations).
These are separate refinement indicators, not a continuum or
finite-amplitude error enclosure.

This supplies a computed transverse characteristic leg at $U=U_a$.
The [same owner](nsc-atlas-bulk-response.md#continuation-to-the-wall-an-unevaluated-analytic-route)
derives a route to the emerging wall strip using radial mass/lapse
reconstruction and the already derived wall-displacement law. That strip,
the smooth cap, observer clock and higher-order field correction remain
to be calculated. The coefficients above do not select a finite physical
amplitude or establish the complete sourced configuration.

## Cap mass and shared wall displacement, 15 October 2026

The [two-sided boundary calculation](nsc-atlas-cap-wall-response.md) derives
the cap's vacuum-subtracted mass from the same leading field as the exterior
functional. The cap vacuum is evaluated at the actual moving radius. Its
frame transformation supplies an additional potential contribution that
the exterior's Ricci-flat expression does not contain.

Both mass coefficients determine a quiet-anchored linear equation for the
quadratic wall displacement. It divides by the positive wall speed, stays
regular at the signed exterior-normal crossing, and satisfies the angular,
temporal and normal junction conditions with the shared stress transport.
The [record](../results/development/nsc-atlas-cap-wall-response-v1/cap-wall.json)
binds the exact algebra and the
[common model specification](../results/development/nsc-atlas-cap-wall-response-v1/model-specification.json).
The next dynamical input is the compatible bulk radius and lapse response,
followed by corrected observer clocks and field propagation. Numerical
displacement values and a finite-amplitude remainder bound have not yet
been calculated.

## Exterior wall-mass closure, 13 October 2026

The [analytic functional](nsc-atlas-wall-mass-closure.md) connects the
already known wall field to a conditional quadratic geometric mass on
the smooth exterior side. Both null mass fluxes contribute along the
timelike wall. Their integral and the full conformal endpoint term use
only the fixed leading trace, its first derivative and reference geometry.
The quiet corner fixes the integration constant. Original-frame transport
with the improved stress independently agrees, and the expression remains
regular where the signed exterior radial normal vanishes.

The [record](../results/development/nsc-atlas-wall-mass-closure-v1/wall-mass.json)
binds exact algebra and unchanged prerequisites. No new field solve or
numerical mass curve is supplied. The wall's field tail does not justify
the initial source leg's mass plateau or a late Jordan mass sign. This is
a necessary coefficient for a compatible perturbative family, not proof
that the coupled cap/wall problem has a solution. The next decision is
controlled evaluation of the functional or closure of a named remaining
cap/wall input. Observer corrections, finite-amplitude validity and other
atlas routes retain their separate questions.

## Conformal scalar normal-force compatibility, 13 October 2026

The [normal-force check](nsc-atlas-scalar-normal-force.md) closes the
regular normal-stress jump left open by the earlier surface calculation.
The full canonical conformal stress, scalar transmission law and Gauss
identity give zero additional jump at quadratic order on the fixed wall.
The cap's Einstein-tensor contribution is essential to the cancellation.

The [record](../results/development/nsc-atlas-scalar-normal-force-v1/normal-force.json)
checks the direct projection, geometric identity and force sign with exact
algebra, without a new waveform or pressure curve. With fixed vacuum terms
and material tension, a later compatible perturbation must have zero
variation of the mean extrinsic-curvature trace. This variation includes
metric and wall displacement; it does not set either correction to zero.
The analytic successor above now closes exterior-side geometric mass
transport from the existing leading wall traces. Coupled cap/wall evolution,
observer corrections and finite-amplitude validity remain open. Earlier
dated entries retain their original scope.

## Geometric mass on the prepared null surface, 13 October 2026

The [mass calculation](nsc-atlas-source-mass.md) closes the transverse
radius equation on the initialized source leg, using the same compact
conformal field preparation, corner data and fixed null lapse. It
determines a coordinate-invariant geometric Misner--Sharp mass coefficient.
The Jordan and Einstein frame coefficients generally differ during the
input; the full conformal transformation includes both scalar derivatives.
They agree after the compact trace, where the coefficient is negative and
constant on the remaining finite leg.

The [record](../results/development/nsc-atlas-source-mass-v1/mass.json)
compares frame-converted flux, mixed-radius geometry and the full improved
Jordan stress. The sign follows analytically from the source support and
white-hole-sector orientation. This is not a vacuum future, global
mass-loss history or asymptotic energy. No finite amplitude, wall motion,
corrected interior waveform or stability result is supplied. The next
connection is to carry compatible sourced geometry toward the moving wall
and test its dynamical junction conditions. Other atlas routes remain open.

## Characteristic gravitational initialization, 13 October 2026

The [null constraint](nsc-atlas-characteristic-initialization.md) now determines
the initial areal-radius correction required by the same prepared conformal
probe. With the corner geometry and affine normalization fixed, the
quadratic coefficient $h$ and its affine slope $h'$ are negative after the
compact input. The unchanged lapse and radius cannot both remain reference
functions on that sourced leg. Independent integrated and smooth-exterior
Einstein-frame expressions check the result.

The [record](../results/development/nsc-atlas-characteristic-initialization-v1/initialization.json)
retains the unchanged field preparation and the declared gravitational data.
This is an intrinsic initial-leg constraint, not a mass shift, wall
displacement or evolved backreaction. Zero $T_{VV}$ after the trace does
not remove transverse scalar derivatives or establish a vacuum future.
Coupled bulk evolution, normal-force balance, finite-amplitude error bounds
and stability remain open. The next useful connection is a mixed or
transverse equation linking these initialized data to an invariant exterior
response or to the wall. Other atlas routes remain available.

## Conformal probe stress and shell compatibility, 13 October 2026

The [analytic check](nsc-atlas-scalar-wall-consistency.md) fixes the scalar's
canonical normalization and connects its quadratic stress to the existing
wall. The full improved stress has zero singular surface coefficient and
zero regular tangential flux jump. The published scalar-tensor junction
law independently gives the unchanged material tension. Keeping only the
canonical kinetic flux would miss the cancellation.

This removes those two local obstructions for the original smooth compact
preparation; no wave evolution was rerun. The
[record](../results/development/nsc-atlas-scalar-wall-consistency-v1/consistency.json)
binds exact algebra and unchanged finite prerequisites. Normal-force balance,
bulk Einstein response, an amplitude-validity bound and stability remain
uncomputed. The next useful connection is the pulse's smooth bulk stress and
the geometric correction it requires, under a separate bounded admission.
Physical source selection and the other atlas routes remain open.

## Fixed-model observer inverse and preparation exclusions, 12 October 2026

The [proof](nsc-atlas-observer-inference.md) and
[record](../results/development/nsc-atlas-observer-inference-v1/inference.json)
connect a complete compatible signed central history to its unique source
on the same retained v7/v8 wedge. The result borrows resolvent mathematics;
direct differentiation and an independent analytic review check the PDE,
wall, initial data, regularity at the corner and central normalization.
There is no numerical inverse solve or new waveform in this cycle.

The inverse covers initially quiet C2 signed preparations, with unrestricted
later shape. Requiring a nonzero nonnegative source to end at the original
cutoff imposes a positive post-input field and negative accumulated signed
field. A smooth entirely positive observer pulse requires continuing late
input and is excluded by that restricted class. This is a conditional
within-model test; no second physical geometry has been shown equivalent.

The next unlock is an independently justified source restriction or a
source-linked additional observable before comparing geometries. Formation,
remote preparation, quantum state and the other atlas routes remain open.

## Prepared interior response and finite positive tail, 12 October 2026

The eighth atlas evaluates the v7 finite exterior preparation operator with
its geometry, conformal scalar, compact source and initial cap held fixed.
The [scientific owner](nsc-atlas-prepared-response.md) and
[response record](../results/development/nsc-cosmological-atlas-v8/prepared-response.json)
connect the retained Schwarzschild-IV data to the moving-wall trace and
central observer. This closes the numerical `wall-trace` and
`observer-response` connection for those independently specified data.

The normalized central field Y=ell phi_C/epsilon has sampled extrema
47.29833694 and -45.69595941. After the exterior input ends, a positive
field tail reaches the centre over
-0.340772595<=T/ell<=-0.250230701, with terminal Y=0.00129049806.
Its sign follows from an analytic positive-kernel bound. Three independent
reference-integral brackets contain the computed trace; spatial refinement
reduces the maximum waveform difference from 1.11704e-6 to 2.79286e-7, and
the post-input derivative-trace difference is 3.70648e-9. These are numerical
diagnostics supported by analytic bounds, not rigorous interval enclosures,
energy fractions or a claim about later persistence.

An initial unsegmented solve stopped at its 6,000-call cap. The admitted
event-aligned correction completed with at most 4,352 calls per solve,
5,597,264 kernel values and the unchanged 60-CPU-second guard. The owner
preserves both outcomes and the unchanged theory, source and tolerances.
The right-parent P exclusion remains; this signal is prepared in retained
IV, not driven from the ordinary right exterior.

The next named unlock is
[observer identifiability under preparation freedom](nsc-cosmological-atlas.md#observer-identifiability):
which exterior features can an interior signal distinguish when different
input preparations are allowed? The analytic successor above resolves the
fixed-model inverse and source-class restrictions. Remote generation, quantum state, backreaction, formation and
the many-body future remain open. All other atlas routes and historical
records below are retained.

## Finite exterior preparation and programme, 12 October 2026

The [seventh atlas](nsc-cosmological-atlas-v7.md) connects a finite
characteristic preparation to the existing interior response operator. The
fixed conformal scalar and wall admit a retained white-hole-side strip with
0.4 <= R/ell <= 0.55, wall endpoints 0.45 and 0.55, independent compact
characteristic data and zero earlier cap data. The wall's outgoing condition
follows from the absence of a reflected central signal in this window.

The [derivation](nsc-atlas-field-preparation.md) supplies a wedge integral
equation and an analytic operator-norm bound q* = 0.1388209923 < 1. This
establishes a unique finite conditional solution and a reduced-field bound,
without a PDE solver. The [record](../results/development/nsc-cosmological-atlas-v7/preparation-access.json)
keeps the causal, domain, coefficient and unit checks. Its central arrival
window is approximately -0.581759 <= T/ell <= -0.250231. The numerical waveform
from these data is unevaluated. The old pulse's late zero field is not a
condition imposed on this successor.

The designated right-exterior observer P still cannot send this signal.
Remote generation of the IV data, formation, stability, quantum state and
backreaction remain open. The [programme graph](../results/development/nsc-cosmological-atlas-v7/programme.json)
records the prerequisite, decision and downstream unlock for each candidate
calculation. It is distinct from physical causal arrows. All earlier scenarios,
multiple-offspring, communicating-mouth and exterior-response alternatives remain.

## Exterior-first atlas, 10 October 2026

The [sixth atlas](nsc-atlas-field-transfer.md) adds an explicit
classical field and initial-boundary state to the earlier causal map. A
massless conformal scalar in the spherical mode crosses the fixed wall with
its curvature-induced jump included. For the reduced field R phi the shell
jump cancels, while the smooth exterior scattering problem remains open.
The cap begins with zero field and time derivative; one early incoming
characteristic is prescribed on the exterior wall side. Its remote source
and global quantum state are not supplied.

The exact conditional relation phi_C[T(tau_W)]=I(tau_W)/T'(tau_W) gives a
finite central waveform, with the same clock map controlling duration and
amplitude. The compact pulse is fixed in wall proper time around R/ell=0.5,
with half-width 0.03ell. The
[waveform figure](../results/development/nsc-cosmological-atlas-v6/field-response.png)
and [field record](../results/development/nsc-cosmological-atlas-v6/field-transfer.json)
record its input, received support and independent checks. The central field
returns to zero after this pulse through the chosen t_C=0 observation limit.
Zero incoming trace outside the pulse is part of the state specification;
this does not rule out uncomputed exterior tails or memory in another field.

The [v6 paper](https://doi.org/10.5281/zenodo.23090015) connects this calculation to
the [common research thread](papers.md#the-common-research-thread): boundaries,
geometry and state shape local measurements. The focused companion's
mode-space complement has not become a geometric parent merely by analogy.
Local spherical regions need not form a globally spherical multi-hole parent.
All v1-v5 evidence stays frozen, with the multiple-daughter, communicating-mouth,
minimally coupled scalar and exterior-state routes still open.

The [fifth atlas](nsc-atlas-causal-classification.md) supplies a causal map
around the central interior observer, using the published isolated completion
of the fixed de Sitter/Schwarzschild wall. Its
[spacetime diagram](../results/development/nsc-cosmological-atlas-v5/spacetime.png)
matches interior and parent charts through the same wall events. The
[observer access table](../results/development/nsc-cosmological-atlas-v5/observer-access.md)
uses separate central, wall and static-parent proper clocks.

Earlier centre-to-wall, wall-to-centre and centre-to-parent geometric light
paths exist. Emissions at or after the chosen slice cannot follow those paths
in this reference. Earlier signals can still arrive later, including
arbitrarily late central reception of wall emissions approaching its de Sitter
horizon crossing from below. The central escape threshold to the right parent
is t_C/ell=-1.2312944525 in the recorded clock convention; equality follows a
horizon generator and has no finite parent reception. The parent escape
horizon continues through the interior and differs from the centre's de Sitter
horizon. The frozen v5 record left field transmission and amplitudes
unspecified; v6 supplies the conditional scalar connection above.

This reference has a singular past. It supplies neither regular formation nor
the global future of the glued two-body parent. The
[illustrated manuscript v5](https://doi.org/10.5281/zenodo.23090015) preserves those
limits and the independent multiple-daughter, communicating-mouth and
exterior-response routes. That edition identified an early field/state imprint
as the next connection, supplied conditionally in v6. Alternative attachments
still need their own causal channels. All v1-v4 evidence stays frozen; the v5
calculation applies known junction and causality mathematics rather than
claiming a new spacetime classification.

The [fourth atlas](nsc-atlas-connection-slice.md) supplies one explicit
sourced initial slice: a regular flat de Sitter cap, the fixed-tension v3 wall,
a vacuum Schwarzschild slicing transition, a regular bridge and an exact AE
end. Bulk constraints and shell energy/momentum jumps pass; independent metric
derivatives check the vacuum calculation. The source is held fixed. The cap is
a finite spatial ball, not an assumed closed S3 universe.

The published many-body theorem acts only on the smooth vacuum end. Exact open
overlaps preserve the sourced remainder. This justifies a common-parent
initial-data existence family for body A, M_A/ell=0.4, and an ordinary vacuum
Schwarzschild body B, M_B/ell=0.6. B's other end remains: two total AE ends, one
designated common parent end. B has no established cosmological child. The
connector metric and numerical separation threshold are not supplied, and the
[v2 Brill-Lindquist slice](nsc-cosmological-atlas.md) is unchanged.

For an assembly retaining A's right vacuum data through R_cut/ell=4 and an open
overlap, a static observer at R_obs/ell=1.6 sees the same local vacuum geometry
as any other assembly with those retained data for
0 <= tau_obs/ell < 2.4812627896. This is an open local proper-time window and a
demonstrated geometric degeneracy, not an arrival prediction or shared cosmic
age. Field and quantum states require their own initial data.

The [illustrated manuscript v4](https://doi.org/10.5281/zenodo.23090015) and
[source-bound records](https://doi.org/10.5281/zenodo.23090015)
develop the connection and its limits. Formation, global shell evolution and
stability remain open. The v3 conservation obstruction for a moving
radiation/fixed-tension/vacuum pairing and the independent branching, mouth and
response routes remain in the map. All earlier atlas bundles and sources stay
frozen. This closes one project interface using borrowed mathematics, without
a new-to-world mechanism claim.

The active objective is to map possible parent regions, black-hole transitions,
children, siblings and descendants, then investigate their local consequences.
The [finite atlas](nsc-cosmological-atlas.md) records four distinct
connections: descent, attachment, causal influence and correlation. Its first
map has two separated black-hole interfaces, unequal possible children and a
further generation. A separate scenario imports multiple daughters from one
hole. These are conditional constructions; the missing physical junctions,
sources and clock conversions stay visible.

The [research inventory](research-inventory.md) includes both
worktrees, child-and-neck research and charged, influence, constants and clock
investigations. File discovery is distinguished from verification of a
derivation. Saved evidence and stopped campaigns retain their historical roles.

The [reuse/admission owner](prior-art-and-open-claim.md) now admits finite
compatibility, response and inheritance questions. The universal-action
programme is one longer-term candidate. The existing relation
$\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta$ and
$\Gamma_{\mathrm{one}}$ realization are reused within their domains. They
are not compulsory assumptions for every atlas branch. No cosmological dark
fraction, antimatter identification or cross-universe age is a result of the
map. The [literature cards](arxiv-studies.md#atlas-reuse-library-10-october-2026)
state what was actually read and what each source supplies.

## RN radial parent, 4 October 2026

The scientific starting checkpoint is
`5cb775a751add5f72eb06cb681ca16c24a3527d2`, which produced the
recorded NF256 confirmation below. The external two-room import
is the historical commit
`d95175c25a9989aea681967d6089aa7b2bf27706`. It is not an RN result
and not the live branch identity. Live branch commands are in the
[instructions](instructions.md).
The periodic parent below is historical input. Further periodic
autonomous exchange runs stay off this programme. This retained candidate
asks for a nonperiodic horizon-regular Painlevé–Gullstrand parent exterior
and a distinct near-throat child on its own areal radius. Exact
magnetic Reissner–Nordström calibration comes before any sourced case.

The source-free core is independently confirmed. `_measure_case`
builds $R_4$ and Ricci squared from the reconstructed PG jets and
compares them with the reference. `injected_analytic_jets` is false.
At $\mu\in\{1.002,1.01,1.04\}$ the worst residuals are
$R_4=1.2732490227040927\times10^{-5}$, Ricci squared
$4.3332330528444984\times10^{-5}$, $\beta_{rr}=1.480009\times10^{-5}$
and horizon gap $3.198843\times10^{-6}$. They pass $10^{-4}$. No
calibration record is saved.

The reviewed ingoing angular coupling is $+\alpha$. Classical flux
$q=1$ uses $C_F=0.18007744731304431$ and $G_N=0.441906928009764$,
rescaled from the audited flux-$4$ sector so $P^2=1$. That quantized
sector stays distinct. Rank $5$ carried an outer replica of $0.435$.
The rank-$9$ preparation has dominant-lobe peak $12.03$, mean $12.29$
and RMS $2.94$. Its declared finite frequency is not the requested
width $2$.

The driver, episode note and driver tests now include
`sat_norm_ledger` with twice the raw SAT coefficient. The saved
paired run is archived at
[rn-first-infall-pilot-2026-10-04](https://doi.org/10.5281/zenodo.23090015)
with [receipt.json](https://doi.org/10.5281/zenodo.23090015),
sha256 `e3e8226f7d748ea0a62068c4efa37a7ec9e3d5ff6bd3f4923a9536dd4d77f130`.
The wrapper `nsc_rn_pilot_march.py` is
`ef45c5eba4c66f922be0ec04159f208861e5e5e7299d319ce8dd10ea8483ee47`.
The receipt states that the driver hash was not captured. Commit
`89f0fadfb886ab73077c6c39d1ac7cd6613c8979` is the retrospective source
match for the archived PG, source, observables and reference inputs.
The archived inputs match those recorded source objects. That commit
did not produce the run, and it is not the live branch identity.
At $T=5$ and $T=10$ both arms have corrected probability closure below
$9.6\times10^{-13}$. At $T=10$ the coupled inner mass is
$1.040990280$ against the fixed value $1.04$, and the horizon moves
from $1.325655540$ to $1.330239387$. The frozen arm is an initial
sourced-metric counterfactual. Probability accounting is reconciled.
The saved frames do not store derived energies or curvature. The
driver command `--assess-archive` reads the rank-$9$ directory and
takes no time step. The official writer then returned `existing
outputs are sealed` and wrote no bytes. The accepted body is the
creation-only file
[nsc-rn-pilot-assessment-v1.json](../results/development/nsc-rn-pilot-assessment-v1.json),
bound to frozen producer `64e2c36fa2bd0688dd02a0b1bcda714c4767c294`
after the prior code commit. That commit is not the live source now
being edited, and the raw JSON is left as written. Over $41$ frames,
$\Delta E_H=-0.0022440861341528895$ with quarter-step gap
$-8.364479832214181\times10^{-7}$ ($0.0373\%$ of the change) and
$\Delta E_N=-0.0015357675638723422$ with quarter-step gap
$-8.965835266506236\times10^{-7}$ ($0.0584\%$). The weighted SBP
boundary-plus-SAT identity matches the direct derivative to
$1.5\times10^{-12}$. The outer normal spike $7.44\times10^{-5}$ is
the outer SAT. The measured inner bulk residual is $2.98\times10^{-5}$
at $t=8.5$, and the retained commutator is $3.82\times10^{-6}$. These
are not continuous bounds. Dynamic $R_4$ reaches
$7.688\times10^{-5}$, with $3.69\times10^{-5}$ on the full $t=10$
sample and $9.18\times10^{-7}$ in the ten-point curvature bulk. The
source Ricci-squared stress check is about $2.04\times10^{-6}$ in the
bulk, and the full source vacuum departure peaks at $0.001525$ at
$t=7.75$.

The finite continuation from that archived $t=10$ frame is
[nsc-rn-neutral-continuation-v1](../results/development/nsc-rn-neutral-continuation-v1/summary.json).
Its handoff payload is `t10.0000.npz`, sha256
`7586472db0e6a29d2b4dea6a02429bf5c2c28cecc6812f735e02647bc134ef74`.
The run uses $\mu=1.04$, $\varepsilon=0.001$, rank $9$ and
$\Delta t=0.001$, with no reset. It records $564.993$ CPU seconds
and $205504512$ bytes of memory. At $T=20$ the coupled excision
capture is $0.9318536299678606$, the outer outflow is
$0.021424558047058832$, and the remainder is
$0.04671377650607762$. The coupled inner mass is
$1.0410282460948723$, the horizon is $1.3303665716496194$, and
$E_{\mathrm{coordinate}}=2.0207173180768176\times10^{-5}$. The norm
interval closure is $-1.539302149964174\times10^{-13}$. The saved
CAR eigenvalue $6.999561444355473\times10^{-5}$ is the load-time
coupled $t=10$ cache on both arms. The recomputed current values
are $2.6763819258953356\times10^{-5}$ coupled and
$2.681557087002335\times10^{-5}$ fixed. The states are valid PSD
matrices. The cached number is not the current CAR, and the old JSON
is not rewritten. The current code computes a fresh CAR and leaves that saved JSON
unchanged. Refinement and a held prediction remain pending. This
finite segment is not a completed neutral programme and not a
production seal.

The finished radial operator stores the same $D(2,1)$ coefficients
in a dense reference and a CSR production matrix. Their step-$0$
difference is zero. The recorded speedups are $2.055\times$ at $869$
points, $10.744\times$ at $1750$, and $30.308\times$ at $3921$. All
six $T=20$ forecasts fit in six CPU-hours, and the current estimate
is about three hours. The default production point factor is one
backend choice. It is not a continuum validation and not new physics.
The writer now allows a new development file and still refuses the
archive, historical results, and any existing output. The current
reproducible assessment is
[nsc-rn-pilot-assessment-v2.json](../results/development/nsc-rn-pilot-assessment-v2.json),
current sha256 `e52bd6e96929d4b2b5f48b2fa114f3364da7bd02312346972f7d9cdc192469d0`,
$78506$ bytes after normalizing historical path labels; its numerical values
remain $41$ authenticated frames and CPU $1.013725$ seconds.
Its producing commit is
`ab03859ae553d197f99430146985598fdeee3382`, the recorded source of
that file, not a live branch label. The manual v1 body remains the
historical writer failure pinned to `64e2c36`. A separate numerical
family, `nsc-rn-neutral-family-v1`, was earlier recorded as still
running within six CPU-hours on at most six one-thread processes. It
is not running. No case has written a `summary.json`; four stopped on
the CPU projection at $t=15.98$, $16.20$, $18.83$ and $19.99$, and the
two $\mu=1.002$ cases stopped at $t=8.00$ with no failure file. Its
physical results are not claimed.

<a id="two-force-gradient-balance-4-october-2026"></a>

### RN pressure diagnostics: semantic correction, 10 October 2026

The [owning note](nsc-rn-gradient-bubble.md) corrects the
interpretation of the recorded radial and transverse orthonormal pressures.
Their equal-magnitude crossings and local logarithmic slopes remain finite
state diagnostics. They do not establish opposing parent/child forces,
mechanical equilibrium or a regular centre. The live v6 consumer withdraws
those interpretations; the v1–v5 records remain unchanged. No new family
reduction or evolution was run for this correction.

The historical v5 reduction used six cases and 177 recorded stations.
Equal-magnitude crossings occurred at 49.7% of stations, with the recorded
frequency falling from about 65% to 31% across its parameter choices. Those
counts do not measure the persistence of a mechanically balanced surface.
Its $4\pi\int |p|r^2dr$ is a coordinate radial pressure integral, not a
surface force or an invariant proper-volume integral.

The attributed Hayward condition $g\geq-1$ was unsupported. For an actual
power law $r^g$, pointwise finiteness at zero requires $g\geq0$; convergence
of an integral is a different question with a different measure. Hayward's
regular-core metric and finite stress are imported examples, not a universal
slope test. Fitting pressure slopes on a positive-radius interval cannot
establish the centre's limit. Even reaching an RN inner horizon $r_->0$
would not establish regularity at $r=0$; curvature, stress and extendibility
of the proposed central geometry require their own checks.

The retained [v5 record](../results/development/nsc-rn-gradient-bubble-v5.json)
also contains the stencil audit. Exact rational weights agree with the
recorded coefficients to $5.1\times10^{-16}$ across 120 audited stencils;
the earlier polynomial-moment guard was ill-conditioned. That numerical
correction does not validate the withdrawn physical interpretations. The
six-case numerical family stopped without a summary in any case; two cases
stopped at $t=8$ without a failure record.

The [two-worktree audit](research-inventory.md) also identifies a
different clock error in the unmerged charged/influence extensions. Their
claimed radial infall quadrature uses $1/(1+\sqrt f)$, while their own PG
reference uses $\beta=\sqrt{1-f}$ and inward radial null coordinate speed
$dr/dt=-(1+\beta)$. Thus those saved quadratures and downstream ratios cannot
currently be imported as physical infall times. A coordinate transit time
would still not be an observer's proper age. The source paths and exact
reading boundaries are recorded in the inventory; the second worktree and
its raw evidence remain unchanged. Its empirical constants table is not a
derivation of universal constants, and the toy influence closure does not
eliminate all possible exterior explanations of local gravitational effects.

The first shell is neutral massless Dirac, $\kappa=1$, a spherical
closed Gaussian, multiplicity 4 counted once, fixed flux integer $q=1$,
sea-free, $V_4=0$, and a positive exterior static frequency. The parent
outer boundary is absorbing at $R_{\mathrm{out}}=32 r_m$. Excision
between the horizons is all outflow. $N$ and $\beta$ stay coupled from
the action. The primary readout is neutral Dirac infall, the mass,
horizon and RN deviation, and the held local response. Renewal lies
outside the acceptance criterion. The charged lowest-Landau Gauss
branch is a later experiment. The MMP negative-Casimir study runs in
parallel and leaves the classical throat unblocked.

Registered owners are
[nsc_rn_reference](nsc-rn-reference.md),
`nsc_rn_pg`, `nsc_rn_source`,
[nsc_rn_observables](../src/recursive_horizons/nsc_rn_observables.py)
and [derive_nsc_rn_parent](../scripts/derive_nsc_rn_parent.py).
The episode note is
[nsc-rn-parent-episode.md](nsc-rn-parent-episode.md).
The later charged branch is
[nsc-rn-charged-source.md](nsc-rn-pg-action.md).
The parallel MMP study is
[nsc-rn-mmp-source-study.md](nsc-rn-pg-action.md).
Root freezes completed producer files at integration.

The historical parent's bounded $n_f=256$ confirmation is recorded by
that same commit. Both arms stop on the matched centre clock, target
$\tau=2.8568786296681035$, at $2.856878629668112$ and
$2.8568786296680826$. The column-0 child effect is
$0.0004291836505783686$ beside the NF128 effect $0.0004292078493$.
Total child content is $0.0028419435306693414$ beside $0.0028420021143$.
Child proper length is $0.004767169795066806$ beside $0.0047672394727$.
CPU time is $140.664864$ seconds. The record is
[confirmation.json](../results/development/nsc-discovery-parent-cut-confirmation-v1/confirmation.json)
and [confirmation.npz](https://doi.org/10.5281/zenodo.23090015).
The payload is $128332$ bytes,
`32c2b22dd5f95ce6b9bb51c9ee865be6a4465b9c28bc6acbd0147b9d99dac1ef`.
Eight input hashes and forty-two producer hashes match the current
bytes. The finite-band comparison is an indicator. The record stores
no $R_4$. Source-free and matter-jet assessment is future Stage 0.
The written metadata stays as produced.

## Responsive parent, 3 October 2026

Historical input, started at `eb4798c`. Current routing is the
[exterior-first atlas](#exterior-first-atlas-10-october-2026). Primary finite
preparation of the magnetic branch is achieved, and the frozen,
empty-source and strong-offset controls below are recorded. The NF128
exterior-to-cut response is measured, and its bounded $n_f=256$
confirmation is recorded above. The renewal programme is not complete.
No autonomous renewal is asserted.
The
[2 October discovery programme](#active-discovery-programme-2-october-2026)
remains preserved from `ea6c699`. The incoming gate stays paused.

Owners are the [preparation](https://doi.org/10.5281/zenodo.23090015),
[episode](https://doi.org/10.5281/zenodo.23090015),
[observer](https://doi.org/10.5281/zenodo.23090015),
[response](https://doi.org/10.5281/zenodo.23090015),
[step control](https://doi.org/10.5281/zenodo.23090015)
and [Dirac contact](https://doi.org/10.5281/zenodo.23090015).
The response adapter remains one admitted step.

The magnetic branch replaces the retained initial radius by
$r_{\mathrm{mag}}=\sqrt{\mathfrak m/g}$. That radius stays within $0.321\%$
of the original collar candidate. $Q$, the canonical columns, the population
weights and the observer are unchanged. The balanced width-$1.4$ source
selects the common initial kinetic value

$k^*=4g/(3r_{\mathrm{mag}})\sqrt{(\int\rho^2/Q)/(\int Q)}=0.002621213673327129$.

That choice is an initial condition, not a dynamical law.
Producer `4392350e` records it in
[feasibility-magnetic-v5](https://doi.org/10.5281/zenodo.23090015).
The $n_f=64$ and $n_f=128$ solves meet the retained momentum, hard anchor
and CAR conditions with no fallback. Their physical radial-rate maxima are
$0.1173$ and $0.11632$. The $n_f=256$ initial confirmation is valid, and its
curvature comparison is less settled. The six source records and selection
ledger are
[six-magnetic-nf128-v1](https://doi.org/10.5281/zenodo.23090015).

The episode
[parent-episode-magnetic-t3-v1](../results/development/nsc-discovery-parent-episode-magnetic-t3-v1/manifest.json)
is pinned at `8c47d0c17fa5962f1d2210f55db3a96c6457a82e`. All six cases run
continuously through $T=1$ and $T=3$, use $186.65$ CPU, and record no resets
or clamps. CHECK passed on that commit. The manifest's top-level pure and
pool flags are the inherited plan. The NPZ states, step counts and CPU
ledger are the run. Those bytes stay as written.

Every case loses child content. At $T=3$ the child fractions are $0.0221$
and $0.0231$ for child-heavy, $0.0593$ and $0.0605$ for balanced, and
$0.1062$ and $0.1088$ for parent-heavy, plus then minus. Proper child length
contracts by $97$–$98\%$. Mean radius on the two signs stays within about
$2\%$ of the magnetic radius. Renewal is not observed.

On the parent-heavy branch the child fraction is about $0.174$ at $T=1.25$
and about $0.251$ at $T=2.25$. At $T=3$ that child content is $99.41\%$
(plus) and $98.99\%$ (minus) the original annulus column, so the return is
not shown to be outgoing child recycling.

The balanced confirmation
[nf256-balanced-v1](../results/development/nsc-discovery-parent-episode-magnetic-nf256-balanced-v1/manifest.json)
uses $49.42$ CPU. Child fraction moves from $0.059278$ to $0.059058$, child
length from $0.016201316$ to $0.016201667$, and centre proper time from
$1.4723228$ to $1.4722193$. Depletion and contraction remain. Global $R_4$
at $T=3$ moves from $4741$ to $983$, and the weighted Euler defect from
$0.257$ to $0.050$. Projection curvature stays unresolved. This is not a
continuum instability, and the raw curvature reference is not a replacement.
No $n_f=512$ episode is recorded.

Total energy closes near $10^{-13}$, and CAR near $10^{-12}$. For the
balanced child, the exchange against the normal ledger is
$\Delta E=+0.01565$, with boundary $-0.04919$, pressure $+0.06382$ and
lapse $+0.00099$. The time-quadrature indicator is $0.19\%$.

Completed controls, measured from the saved arrays, are
[frozen parent-heavy](../results/development/nsc-discovery-parent-frozen-parent-heavy-v1/manifest.json),
[weak-$k$ empty source](../results/development/nsc-discovery-parent-empty-source-v1/manifest.json),
[strong $T=1$](../results/development/nsc-discovery-parent-strong-t1-v1/manifest.json)
and [strong $T=3$](../results/development/nsc-discovery-parent-strong-t3-v1/manifest.json).
The strong-$T=3$ record is pinned at `cc238577`. The authoritative weak
offset is $k=0.002621213673327129$ and the strong offset is
$k=0.45388971484879426$. The legacy episode label `common_k` near
$10^{-8}$ is not the weak offset.

On the strong branch the minus case moves from child length $1.76097$ and
fraction $0.34039$ at $T=1$ to length $0.14933276$ and fraction $0.090496$
at $T=3$. The plus case at $T=3$ has length $0.00053105$ and fraction
$0.059670$. The early length is an expansion transient. Renewal is not
demonstrated. Intermediate strong-$T=3$ observation headers still carry
the $T=1$ text, so these values are the array measurements. The
[array-authenticated assessment](https://doi.org/10.5281/zenodo.23090015)
reconstructs the measurements without rewriting those headers. The
[assessment successor](https://doi.org/10.5281/zenodo.23090015)
preserves the earlier measurements and binds the current consumer bytes.

Frozen parent-heavy child fraction at $T=2.25$ is $0.33464$, against
$0.25065$ for the coupled parent-heavy case. Empty-source evolution at the
weak $k$ has child length $0.01633475$ at $T=3$, against $0.0162013$ for
the coupled balanced case. These comparisons separate transport, geometric
clocks and background. They are not a universal verdict.

The [matched strong-$k$ empty-source control](https://doi.org/10.5281/zenodo.23090015)
reaches $T=3$ on both signs at producer `5ce6bd02`, using $33.805365$
aggregate CPU seconds. Both preparations preserve the actual offset and
initial $Q,r$, and solve their own empty-source momentum constraints.
The [matched-clock assessment](https://doi.org/10.5281/zenodo.23090015)
authenticates those states and reuses the frozen field columns by exact
occupation reweighting, with no additional field integration.

On the minus branch, occupied and empty centre proper times reach
$3.0507895$ and $1.9788562$. At the occupied $T=1$ proper clock,
the coupled child fraction is $0.3403922$, compared with $0.3887899$
in the reweighted frozen control. Its coarsened-sampling indicator is
$0.00037945$, rather than an interpolation bound. Late minus clocks
exceed the saved control ranges and are not extrapolated. These results
separate geometric response from maintained child content; the source
changes the history without establishing renewal in this preparation.
The [owning comparison note](https://doi.org/10.5281/zenodo.23090015)
records the remaining constraint, curvature and interpolation limits.

Checkpoint summaries are repaired for new records at `5ce6bd02`; old
records stay unchanged. The completed exterior-population response is
the [cut-response batch](https://doi.org/10.5281/zenodo.23090015)
at producer `7004a791`, assessed in
[nsc-discovery-parent-cut-assessment-v1.json](https://doi.org/10.5281/zenodo.23090015)
and [its note](https://doi.org/10.5281/zenodo.23090015).
At the matched centre clock $\tau=2.8568786296681035$, with
$k=0.45388971484879426$, fixed $c_0$ and the authenticated initial
$Q$, $r$ and $\Phi$, a $+5\%$ parent occupation changes original
child-column content from $0.0385173686504$ by $0.000429207849261$,
which is $1.114\%$ of that baseline. The held $\pm5\%$ forecasts for
that column match within $0.3629\%$ and $0.3709\%$ of the predicted
effect. Total child content on the $+5\%$ arm moves by
$0.00284200211428$ against the locked $0.00285160805944$. Renewal is
not asserted. The bounded $n_f=256$ confirmation is the record cited
in the [RN radial parent](#rn-radial-parent-4-october-2026) section.
The same inheritance
law covers states whose physical frequencies may differ. A thermal
compression reading is conditional, and this recorded state is
non-thermal. The ancestry share $\eta$ is not a new gate. The force
here remains the full state–geometry force; a fixed-$J$ Schur
reduction is not substituted for it. Earlier papers and the paused
incoming gate stay preserved.

Earlier preparations stay in their own records:
[first pilot](../results/development/nsc-discovery-parent-v1/feasibility-nf64-128/pilot.json),
[hard-anchor v2](https://doi.org/10.5281/zenodo.23090015)
and
[width-1.4 hard anchor](https://doi.org/10.5281/zenodo.23090015).
The original $n_f=128$ convergence flag remains invalid there.
[Dirac contact](https://doi.org/10.5281/zenodo.23090015)
remains the off-shell mean-versus-connected reference.
The contour binding
[nsc-discovery-regions-v2-binding.json](https://doi.org/10.5281/zenodo.23090015)
stays pinned at `4f3d5fc`.

<a id="active-discovery-programme-2-october-2026"></a>

## Discovery programme, 2 October 2026

Preserved domain, started at `ea6c699`. Not the active cursor.

Douglas approved the **NSC Discovery Programme: Source-Selected Regions,
Regeneration and Inherited Local Response**, starting from `ea6c699`.
The accepted plan in that conversation was the task authority for this
preserved programme. The first
implementation wave repairs evidence-writing, supplies an equivalent FFT
backend and resumable episode runner, and adds regional observables,
stability diagnostics and a coupled-response specification. Its first
scientific batch continued the saved $T=0.3$ states through $T=1,3$ with
matched controls. Its measured outcome is recorded below.

Implementation owners are the [episode runner](https://doi.org/10.5281/zenodo.23090015),
[physical observer](https://doi.org/10.5281/zenodo.23090015),
[coupled response](https://doi.org/10.5281/zenodo.23090015), and
[FFT adapter](https://doi.org/10.5281/zenodo.23090015).
The integration review corrected work/flux measure and frozen-control
jet defects in the new observer. Dense/FFT, ledger, frozen-control and
analytic directional-response checks pass. The pinned one-thread FFT
benchmark is 12.95 times faster end to end at $n_f=256$.
Intermediate states and regular scalar observations are retained; restart
continues the saved clocks and accumulated work/transport channels.
The first continuation batch ran from numerical checkpoint `1c9e770`.

Scientific batches use a six CPU-hour aggregate budget. Independent cases
run in parallel on the Mac; Codex owns integration and one campaign
executor owns scientific launches. Regular-interior results are imported
foundations, not a bounce-rediscovery campaign. The incoming gate stays
paused. The publication target is one unified article after the central
discovery chain; the existing papers and release remain preserved.

### First discovery continuation

All six cases reached $T=3$, using 277.186 aggregate child CPU seconds.
The [campaign](../results/development/nsc-discovery-episode-v1/manifest.json)
retains immutable station states and the exact $T=0.3$ hand-off. At
$n_f=256$, half-step coupled evolution gives:

| Measurement | $T=0.3$ | $T=1$ | $T=3$ |
|---|---:|---:|---:|
| Child proper radial length | 2.375218 | 1.488905 | 0.00102044 |
| Probability inside the fixed child collar | 1.276791 | 1.499940 | 0.00300920 |
| Normal clock at $x=2$ | saved upstream clock | 1.102546 | 1.351988 |

Radial contraction accompanies growing areal radius and field transport.
This is not a demonstrated regenerative regime. Spatial and timestep
comparisons resolve these large geometric and probability changes. Late
curvature is more sensitive: at $T=3$, $128\to256$ changes the Weyl invariant
by about 28 percent, whereas the child length changes by about 0.046 percent.
An independent full-rate analytic acceleration agrees with the existing
curvature helper; the discrepancy amplifies finite-band jet differences
when $Q$ becomes small. The selected
[$n_f=512$ confirmation](https://doi.org/10.5281/zenodo.23090015)
also reaches $T=3$; its Weyl maximum differs from $256$ by about 0.08 percent.
Four-metric tidal measurements provide the observer-defined readout of this
growth. The frozen
control holds geometry at its $T=0.3$ value; coordinate-time comparisons
are not equal-proper-time comparisons.

### Withheld exterior-source prediction

The [prediction record](../results/development/nsc-discovery-prediction-v1/prediction-run.json)
propagates the full state–geometry tangent from the same source-consistent
initial preparation and compares independent nonlinear cases at equal
proper time. The middle occupations remain $0.5,0.5$; exterior occupations
change from $0.75,0.75,0.25,0.25$ to $0.76,0.76,0.24,0.24$.
At $n_f=256$, the predicted child-content change is $-0.0009409095$ and the
withheld measurement is $-0.0009366253$. Their difference is 0.455 percent of
the predicted effect. Two smaller centred perturbations check the derivative;
the $128\to256$ movement of the measured change is about 0.52 percent.
This is a quantitative conditional response prediction on the declared
finite construction. The secondary mean-radius linear prediction differs
from the withheld result by about 2.1 percent; its finite nonlinearity is
retained. The proper clock is integrated through RK4 stages and the
nonlinear cases reach its target by bracketed steps, rather than by applying
the prediction's linear clock correction to the measured value.

### Source-family outcome and the remaining preparation assumption

The [15-member source study](https://doi.org/10.5281/zenodo.23090015)
varies five occupation imbalances and three packet widths, solving each
member's own initial radius. Fourteen members reach $T=3$; all show radial
contraction, areal expansion, and loss of the original child packet from
its prepared window. Final packet retention ranges from about 0.024 to
0.127 percent. Less than 0.30 percent remains inside the parent collar.
The reversed-imbalance, widest-packet member fails during the positive
initial-radius homotopy and remains unresolved. This is a preparation
failure, not an exclusion of that source.

These runs retain the supplied initial $Q_0$, auxiliary field and geometric
momenta. They test a family of instantaneous constraint-adjusted slices,
not source-selected holding geometries. Width comparisons also change the
prepared measurement window. Confirmation therefore includes a common
fixed child window, selected finer cases and their own frozen-geometry
controls. The [external run binding](https://doi.org/10.5281/zenodo.23090015)
authenticates the frozen `c4a22e5` producer and completed artifacts without
rewriting their original records. Total measured family CPU is 732.65 seconds.

The [finer successor](https://doi.org/10.5281/zenodo.23090015)
solves all five selected sources, including that formerly unresolved member,
and evolves them with three own-source frozen controls to $T=3$. Its total
measured CPU is 620.11 seconds. All selected coupled preparations share the
same qualitative fate. In the common fixed $(1,3)$ child window, the original
child pair retains 0.0388, 0.0593 and 0.0878 percent for narrow, nominal and
wide sources. The corresponding frozen controls retain 0.459, 0.515 and
0.547 percent. These are comparisons at the same coordinate time; their
proper clocks differ. Source width has a quantitative effect, but none of
these packets establishes a holding mechanism. The predecessor's preparation
failure remains preserved in its original record.

The next balance query uses the actual retained spectrum to prepare a
stationary source while leaving the geometric means free.
[Static source-selected control](https://doi.org/10.5281/zenodo.23090015).
The [source-free control](https://doi.org/10.5281/zenodo.23090015) provides
an exact same-action comparison with constant areal radius and tidal curvature
despite contracting radial proper lengths. The
[translation control](https://doi.org/10.5281/zenodo.23090015) moves the
complete state, geometry, observers and clocks across the periodic seam;
it does not replace a future ambient-extent comparison.

The [exact source-free record](https://doi.org/10.5281/zenodo.23090015)
reproduces constant areal radius and bounded tidal invariants while $Q$ and
proper lengths contract. The
[seam translation record](../results/development/nsc-discovery-translation-v1/translation-shift5.5.json)
passes after moving the complete field, geometry and observer across the seam.
Both are finite controls of the implemented action.

The [first static spectral-source attempt](https://doi.org/10.5281/zenodo.23090015)
preserves its best finite iterate after 1,000 optimizer evaluations and
101.68 CPU seconds. The radial seed satisfies its own equation, but the
remaining $\dot p_Q$ and lapse residual maxima are 3.87 and 5.79.
The unsatisfied coupled balance is retained; the optimizer stop does not
exclude the declared static class. Subsequent seed selection uses an
integrated action balance rather than repeating that Jacobian search.

The next [grandchild experiment](https://doi.org/10.5281/zenodo.23090015)
starts from the transported equal-proper-time state, retains nested observer
details and the ambient response, and locks a forecast before generating
the withheld future local measurement. Its physical regional readout is
separate from its modal observer covariance.

### First carrier circuit and the choice of the next local measurement

All six unchanged [continuations](../results/development/nsc-discovery-crossing-v1/manifest.json)
reach $T=8$ using 431.73 aggregate child CPU seconds. The fine half-step child
probability returns to 0.999687; its frozen control reaches 0.994718.
The coupled child proper length is $9.14\times10^{-13}$, while the frozen
length stays 2.3752. This is a carrier-period return candidate with severe
geometric contraction, not a maintained physical region. An ambient-extent
comparison will distinguish the carrier dependence, retaining the changed
initial constraint preparation on a longer domain.

Actual radial and angular tidal components survive the spatial and timestep
comparisons; cancellation-sensitive $R_h$ and contracted Weyl channels do
not. The [crossing assessment](https://doi.org/10.5281/zenodo.23090015)
keeps those verdicts separate and records the cadence gap after $T=5.8$.
The clock at $x=2$ reaches 1.3521155, advancing only about $1.27\times10^{-4}$
between coordinate times three and eight.

The [short grandchild record](https://doi.org/10.5281/zenodo.23090015)
preserves a locked forecast followed by an independent continuation. Its
regional-content change is only $-1.68\times10^{-11}$; the linear prediction
differs by about 5.12 percent of that tiny effect. This interval supplies
too little regional response for a headline prediction. The successor moves
the readout to $\Delta\tau=0.75$, near coordinate time one, rather than
refining the early tail. The separate
[width response](https://doi.org/10.5281/zenodo.23090015) changes the
physical preparation on the same chart and predicts a withheld width;
it is not a coordinate rescaling or an assumed universal power law.

The [later grandchild measurement](https://doi.org/10.5281/zenodo.23090015)
reaches the declared proper clock 1.10785546 from the transported opening
state. Its locked linear prediction is $-0.000440453$ and the measured
regional-content change is $-0.000450768$, a 2.29 percent remainder relative
to the measured effect. Direct and sequential causal routes reconstruct the
same field and covariance at numerical precision, on one shared generated
geometry. The [numerical confirmation](https://doi.org/10.5281/zenodo.23090015)
uses 169.56 CPU seconds. Timestep halving moves the measured effect by
0.000046 percent; the $128\to256$ spatial comparison moves it by 0.0773
percent. Both indicators are much smaller than the prediction remainder.
They do not bound the opening prefix or isolate finite-perturbation error.
This is a measured inherited response, rather than a universal scaling law.

The [physical width experiment](https://doi.org/10.5281/zenodo.23090015)
changes the source width from 1 to 1.05 with its own source-consistent radius
and evaluates at the same proper clock. The measured change is $-0.0895641$,
compared with the first-order forecast $-0.105968$. Their difference is about
18.3 percent of the measured effect, far above the available frame and event
indicators. The finite-width response needs a nonlinear prediction. The
existing record remains a completed effect measurement and a limit of the
linear approximation; spatial and timestep refinement remain unmeasured.

The [nonlinear successor](https://doi.org/10.5281/zenodo.23090015)
derives a second proper-clock response coefficient of 9.55010 from four
neighboring preparations, without fitting the earlier width-1.05 result.
For its new withheld width 1.03, the measured content change is
$-0.0579681$. The linear discrepancy is 9.68 percent of that effect;
the quadratic discrepancy is 2.27 percent. It uses 88.50 CPU seconds and
preserves the earlier result. The subsequent
[confirmation](https://doi.org/10.5281/zenodo.23090015)
uses 53.16 CPU seconds. Timestep halving moves the effect by
$3.83\times10^{-10}$; the spatial comparison moves it by about 0.992
percent. The quadratic prediction remainder is about 2.27 percent; these
comparisons remain indicators rather than a continuum bound.

### Carrier dependence and an admissible dynamic preparation

The [period-eight/twelve intervention](https://doi.org/10.5281/zenodo.23090015)
continues the period-twelve state exactly after correcting a real periodic
integration diagnostic. Tagged probability returns near time eight in the
shorter carrier and near twelve in the longer one, in coupled and frozen
branches. At $T=8$ the coupled period-twelve child fraction is
$5.52\times10^{-10}$; at $T=12$ it is 0.997158, with child proper length
$1.80\times10^{-15}$. The changed initial radius/lapse is about twenty
percent, so this is an ambient intervention with changed constraint
preparation, rather than exact local Cauchy-data matching.

The [saved-state free comparison](https://doi.org/10.5281/zenodo.23090015)
provides a mechanism for the coupled return. Exact antiperiodic free
transport of the decoded $T=3$ columns reproduces the $T=8$ field to
$2.98\times10^{-6}$ relative weighted error and the fixed-window probability
to $9.96\times10^{-13}$. The conditional unitary comparison needs only
$\kappa\int\|Q\|_\infty dt$, without an exponential growth factor; that
integral has not been bounded from the saved endpoint samples. This return
is explained by carrier transport and does not establish renewed structure.

The [general positive-radius preparation](https://doi.org/10.5281/zenodo.23090015)
removes the constant-$Q$ initial-data restriction. With $\chi=0$, zero
momenta and the actual positive standing field modes, both a patterned
$Q$ and its uniform control solve the owned initial constraints. Retained
lapse residuals are about $1.9\times10^{-9}$ and full residuals about
$2.0\times10^{-8}$; the shift residual is about $3.2\times10^{-14}$.
The calculation uses 0.139 CPU seconds. The $Q$ pattern remains a declared
initial input; the radius is source-consistent, and subsequent motion must
be generated by the existing coupled equations.

The [Gaussian geometric influence](https://doi.org/10.5281/zenodo.23090015)
also has a nonzero connected third bilinear vertex,
$\kappa_3=1.14989\times10^{-5}$, although elementary fields remain Gaussian.
Its actual child/parent instantaneous response slope is $-0.0150914$.
The occupied/empty cut and an independent counting-field check reproduce
these relationships. This uses the existing mixed field/geometry coupling;
it inserts no cubic interaction and identifies no cosmohedron equivalence.
The representative-channel result does not yet determine total angular
stress noise.

The [action-balance selector](https://doi.org/10.5281/zenodo.23090015)
finds three curvature/source/flux-balanced shapes from 33 evaluations using
2.16 CPU seconds. Harmonics 8, 10 and 12 select cosine amplitudes 0.16416,
0.28722 and 0.34671. The local force and lapse equations remain unsatisfied.
The [local-critical query](../results/development/nsc-discovery-stationary-critical-v3.json)
then eliminates radial and auxiliary critical equations and tests the remaining
shape force on a named CAR-admissible occupation branch. Its 600 evaluations
use 31.26 CPU seconds and do not find a local critical state. The selected
shape approaches uniformity: mean $Q$ falls to $0.0009454$, while mean areal
radius grows to 1287.0. The normal-radius denominator approaches zero, yet
the retained $p_Q$ force remains 32.96 and the full lapse residual 35.07.
This is a noncompact, poorly conditioned search tendency, not a proof that
all stationary or regenerative states are absent. That search stops here;
a dynamic preparation must satisfy the initial constraints without also
being a stationary solution. [Saved-array diagnostic](https://doi.org/10.5281/zenodo.23090015).

The [independent limit assessment](../results/development/nsc-discovery-stationary-limit-v1.json)
derives the escape analytically: the mean conformal scale is proportional to
the vanishing shape amplitude, the areal radius is inversely proportional,
and their product tends to 1.216734. The predicted occupation limit 1.218715
matches the saved 1.218714. The old checker rejects a reordered gradient
identity at its absolute tolerance; that failure is retained. A stable
Einstein coefficient improves the affine arithmetic, while the actual local
force remains unsatisfied. Neither numerical conditioning nor this limit
is treated as a physical no-existence result.

The [source-consistent wave episode](https://doi.org/10.5281/zenodo.23090015)
reaches $T=3$ for all six matched cases using 159.25 CPU seconds. The coherent
source changes child probability by 0.100012; the timestep comparison differs
by $3.34\times10^{-14}$. Uniform and patterned stationary field preparations
retain their quarter-domain content. All three geometries contract radially
and expand in areal radius in this initial sector; spatial confirmation has
not been run.

The [initial-sector control](https://doi.org/10.5281/zenodo.23090015)
tests a consequential preparation choice. Constant $Q$, $\chi=0$ and zero
momenta force $\dot Q(0)=0$, $\ddot Q(0)=-Q^3$, regardless of the field source.
An admissible initial expansion/contraction rate changes extrinsic data;
$\chi=-2$ instead changes the initial invariant curvature. These controls
retain the same action, field source and observers and impose no rate after
initial preparation.

The [six initial-sector continuations](https://doi.org/10.5281/zenodo.23090015)
reach $T=3$ using 178.86 CPU seconds. Initial radial expansion, initial
contraction and the $\chi=-2$ control all end with smaller proper radial
lengths and larger areal radii. The same uniform source stays uniform;
these controls do not demonstrate a localized holding region.

The source also fixes a unique initial curvature when its initial auxiliary
acceleration is required to vanish. This is a conditional preparation
criterion, not a stationary-state or universal scale-selection law. The
[balanced successor](https://doi.org/10.5281/zenodo.23090015)
tests that value and fixed $\pm1$ percent controls using 137.43 CPU seconds.
The candidate contracts from child length 2.5497 to 0.0023929 by $T=3$.
The lower-curvature control reaches the positive-radius chart edge near
$T=1.482$; both timesteps retain the last admissible state and bracket the
event without clamping. The upper control contracts as well. Actual tides
are measured separately from the growing auxiliary field. These outcomes
stop this preparation-tuning route; they neither establish a singularity
nor exclude other action branches or nested constructions.

### Action-domain finding and the next coupled comparison

The [constrained curvature assessment](https://doi.org/10.5281/zenodo.23090015)
derives a nongauge extra mode from the owned local action, including both
linearized constraints. On the proper-flat patch its dispersion is
$\omega^2=k^2-m_*^2$, where $m_*^2=-A/(2C_W)=26.6667$.
Its extra pole has the opposite residue from the Einstein pole. The
implementation and declared Lorentzian action agree; this is not a sign
bug or an all-source instability theorem.

The [existing curvature-EFT branch](https://doi.org/10.5281/zenodo.23090015)
uses a perturbative local expansion without those extra initial data.
Treating the local four-derivative truncation as an exact Cauchy theory
therefore tests a different branch. Late tidal scales also leave the small
curvature regime; raw two-dimensional $R_h$ must not be used as physical
four-dimensional curvature in that check.

The [physical-unit assessment](https://doi.org/10.5281/zenodo.23090015)
remeasures 24 saved checkpoints. For the uniform preparation,
$|C_W|/A$ times the normal tidal scale rises from 0.00249 initially to
0.4295 at $T=1$ and $4.05\times10^6$ at $T=3$. Patterned and coherent
sources show the same loss of scale separation. The conditional balanced
preparation delays that growth, but does not retain small curvature at
$T=3$. Proper generator levels and the actual source-rate scale are reported
separately; high retained field energies are not treated as a sharp cutoff
failure by themselves.

The [independent work control](https://doi.org/10.5281/zenodo.23090015)
seals an analytic forecast before imposing one uniform geometric pulse.
The prepared full finite CAR state gives up $2.2678\times10^{-5}$ of
energy, within 0.085 percent of the forecast. Transition-probability work
and direct integrated work agree to $3.9\times10^{-18}$; preparation,
forecast and measurement together use 0.664 CPU seconds. The pulse is
prescribed, so this is a work-response control rather than autonomous
regeneration. A mathematical filled-negative reference absorbs energy.
Reference subtraction agrees for the mean, but the Gaussian characteristic
ratio and connected noise differ. Physical vacuum/influence matching is
not inferred from the equal mean.

The next finite comparison reuses the leading common action, unchanged
Dirac preparation, proper clocks and field–geometry feedback. It owns the
four leading canonical variables rather than deleting fields from the
singular auxiliary chart. The first-order Weyl contact for this Gaussian
Dirac source remains a separate derivation; the scalar specialization is
not copied into it. Earlier exact-branch measurements remain preserved in
their declared numerical domain.

### Completed leading comparison and source-dependent energy transfer

The [leading comparison](https://doi.org/10.5281/zenodo.23090015)
reaches T=3 in all eight cases. The original raw reaction norm mixed canonical
coordinate scales; the [numerical successor](https://doi.org/10.5281/zenodo.23090015)
changes only its timestep indicator, preserving the exact saved states and
all equations. Parent and successor ledgers charge 470.46 CPU seconds together,
within the original 600-second allowance.

At nf256 with half steps, the coherent case ends with child proper length
$4.91\times10^{-9}$, mean areal radius 0.274979, and child probability
0.798248. Its matched frozen control retains length 2.608961, radius 2.811734
and probability 0.769027. These are coordinate-time comparisons with separate
carried clocks. Spatial movement of the final angular tide is about 0.00098
percent; radial tide movement is 1.08 percent. The large contraction is
resolved, but no maintained regional size or regular continuation is inferred.
Late curvature is outside the established local-EFT separation.

The finite field's coordinate energy changes from +12.74674 to -0.543702
in the uniform case, and from +13.44904 to -0.910152 in the coherent case.
Geometry carries the opposite change. CAR admissibility remains intact;
negative instantaneous energy in this finite prepared model is not identified
with antimatter or a physical vacuum. The areal scale turns slightly while
proper radial length continues to contract. Its source-history relation is
the next locked prediction, rather than another bounce-rediscovery target.

The [exact static virial](../results/development/nsc-discovery-leading-virial-v1.json)
excludes only the static positive-energy periodic leading class: source energy
plus magnetic energy would have to vanish. On a finite region, the same action
instead supplies an explicit exterior gradient/traction term. This identifies
the parent response required for holding; it does not exclude dynamic or
oscillatory sectors, negative-energy states, or different boundary domains.

### Source-history-selected turn and independent forecast

The [homogeneous record](https://doi.org/10.5281/zenodo.23090015)
derives an invariant 16-real-state sector of the same leading equations;
it adds no force or law. Three momentum pairs and their spinors reproduce
the stored uniform full-field trajectory. Field energy crosses zero at
T=2.61301056, before the rooted areal minimum at T=2.61410338. The minimum
radius is 0.282440638. The constraint selects it through

$$
r_{\rm turn}^2=\frac{\mu}{8\pi A}
+\frac{E_f}{8\pi A\,\ell Q_{\rm turn}^2}.
$$

The scale depends on the actual transported state and geometry, as well as
couplings, flux and preparation. The T=2.6 sampled minimum is not substituted
for this event. Instantaneous energy-sector transfer is not identified with
antimatter or a physical filled-sea vacuum.

A locked tangent forecast changes the paired occupations from
(.75,.5,.25) to (.752,.5,.248), preserving trace three. The
[independent full-NF32 measurement](https://doi.org/10.5281/zenodo.23090015)
finds radius 0.2827044098 against the forecast 0.2827040102, an error of
0.152 percent of the predicted radius change. The original fine arm remains recorded at its CPU stop. The
[exact-state confirmation](../results/development/nsc-discovery-homogeneous-turn-confirmation-v1/confirmation.json)
resumes it without a source or radius reset, using 24.26 CPU seconds in a
separate 35-second allowance. Its radius is 0.2827044110; timestep movement
is $1.12\times10^{-9}$, or 0.00043 percent of the measured change. The
original record and its 30-second allowance remain unchanged. This is
a conditional finite-model scale prediction, not a new bounce mechanism,
regular-spacetime proof or completed regeneration. Strong-curvature EFT
validity and physical vacuum matching remain open.

### Source-responsive collars and required parent gradients, 3 October 2026

The [local collar calculation](https://doi.org/10.5281/zenodo.23090015)
uses the unchanged leading action with a real standing Dirac source on
$|x|\le1$, with an inner measurement at $|x|\le0.5$. Its center constraint
selects $r_0^2=\mathfrak m/g-Mw/(g\pi)$, $g=8\pi A$; the surrounding
radius, log-lapse and field are then solved together in space. The source
weight is a local amplitude parameter, not an assumed globally normalized
occupation. The actual collar covariance eigenvalue is measured.

Producer `b7aa5f5` seals the analytic first-order response before the new
nonlinear weight $w=0.0013$ is evaluated. The
[withheld measurement](../results/development/nsc-discovery-parent-traction-v1/measurement.json)
records:

| Source-induced change | Measurement | First-order prediction remainder / effect |
|---|---:|---:|
| Inner proper radial length | $-0.0004816671743$ | $0.0180\%$ |
| Larger-collar proper radial length | $+0.002264195085$ | $0.7933\%$ |
| Larger-collar boundary gradient charge | $+0.02671635197$ | $0.6256\%$ |

These are changes between static preparations, not temporal contraction or
expansion. The same geometry supports both length responses. The boundary
charge equals partial field energy plus the magnetic integral, including
its induced change. Tighter numerical settings move the length results by
about $2\times10^{-12}$ and $3\times10^{-12}$; the tighter static identity
residuals are below $4\times10^{-15}$. These remain numerical indicators.
The smaller raw inner-lapse effect has a 1.24 percent first-order remainder
and is preserved rather than retuned. Reported stage CPU is 0.089 seconds
within a 30-second allowance; read-only payload replay passes.

The exterior must supply the derived metric gradients and a compatible
continuation of the actual Dirac field. A globally normalized spectrum,
complete physical parent and autonomous renewal are not established by this
local branch. The [constraint continuation criterion](https://doi.org/10.5281/zenodo.23090015)
derives when a smooth dynamical exterior can carry real, regular momenta
without changing the leading action. Its nonnegative primitive and zero
regularity conditions must be checked on a proposed extended source; an
arbitrary smoothly glued parent is insufficient. No new parent simulation
or replacement roadmap is introduced.

### Coupled reduction with the same generated geometry

The [causal reduction experiment](../results/development/nsc-discovery-coupled-memory-v1/manifest.json)
continues the exact saved $T=0.3$ state to $T=1$. Its retained, exterior-drive
and memory columns reconstruct the full field, including their cross terms,
and supply the force that evolves the geometry. Independent full evolution
agrees within $1.6\times10^{-15}$ in field columns,
$6.1\times10^{-14}$ in source forces and $6.9\times10^{-13}$ in geometric
state components. The clocks follow the evolving metric in both descriptions.
This verifies coupled full-to-reduced force and geometry response on the
declared finite segment; the geometry is not prescribed to the reducer.

### Disjoint energy accounting and actual tidal growth

The [stage-integrated ledger](../results/development/nsc-discovery-balance-v1/coupled-T3.json)
uses disjoint spatial regions and the actual four RK4 stages. At $n_f=512$
on $[0.3,3]$, the child normal-energy change is $-5.42376$, comprising
boundary transport $-1614.30749$, pressure work $+1665.94192$ and
lapse-gradient work $-57.05819$. Its numerical accounting remainder is
$1.46\times10^{-6}$. The whole-carrier remainder is $7.12\times10^{-6}$;
halving the timestep reduces these remainders by roughly a factor of fifteen.
The measured projection contribution is retained separately. Normal energy
is signed in this sector; covariance admissibility does not establish a
positive local-energy theorem.

The [four-metric consumer](https://doi.org/10.5281/zenodo.23090015)
computes tidal components from the evolved metric and its actual projected
jets. At $x=2,T=3$, the radial and angular components are approximately
$-5.24\times10^5$ and $5.63\times10^6$. The angular component changes by
about $0.00012$ percent from $n_f=256$ to $512$. This resolves strong finite
tidal growth without relying on cancellation-sensitive scalar traces or
substituting the auxiliary field for curvature. It establishes neither a
future singularity nor a regenerative return.

## Completed finite parent–child pair, 2 October 2026

The [v2 confirmation](../results/development/nsc-nested-parent-child-confirmation-v2.json)
records `MEASURED_CONFIRMED_NESTED_PAIR`. A child collar $I_C=(1,3)$ lies
inside $I_P=(0,4)$; independent parent and child canonical geometry components
evolve as parts of one full metric under the same action. The inherited
statement is the common functional law under the declared coordinate/clock
pullbacks. Differing evolving profiles need not have identical instantaneous
normalized spectra. [Method and domain](nsc-nested-parent-child.md).

All four reciprocal physical responses meet their one-percent numerical
movement criterion through $T=0.3$:

| Physical response | Effect | Spatial movement / effect |
|---|---:|---:|
| Parent source to child probability | $0.01150720$ | $0.01348\%$ |
| Child source to parent-annulus probability | $0.00848896$ | $0.02622\%$ |
| Parent source to child mean radius | $0.0001759865$ | $0.22315\%$ |
| Child source to parent-annulus mean radius | $0.001109926$ | $0.01600\%$ |

These compare each control's change from its own initial slice with the
baseline change at common coordinate time. The first row uses the new
$n_f=256\to512$ confirmation and its finer effect denominator; the other
three retain the recorded $128\to256$ comparisons. V1's matched timestep
indicators also meet one percent. Reconstructed geometry stays positive, the Gaussian
covariance remains admissible, and the full source/link energy accounting
closes in the recorded numerical domain.

V2 also confirms the parent-to-child mean-radius curve at $256\to512$ with
movement $0.00033864\%$ of its finer effect. The figure uses that additional
comparison; the table retains the primary record's $0.22315\%$ indicator.

The same-trajectory reduction uses actual $\Phi(0.1)$ and fixed T0 observers
on $[0.1,0.3]$, retaining child, parent detail, and ambient source components.
The child occupation changes by $0.294235$; streamed/full discrepancy is
$0.04779\%$ of that change. Parent-detail drive and child/detail cross
omissions move occupation by $0.0331275$ and $0.0509377$, with their own
comparison discrepancies $0.20242\%$ and $0.16284\%$. All active control,
coherence, midpoint, and conditional/autonomous comparisons meet one percent.
These omissions are conditional diagnostics, not additive shares of a
physical cause. The exact cached linear-$Q$ operator changes neither the
source nor the generated geometry schedule.

This completes the declared finite construction; numerical movements are
refinement indicators, not continuum or total propagated state-error bounds.
The sealed [v1 record](../results/development/nsc-nested-parent-child-v1.json)
keeps its original `INCOMPLETE` status and twelve trajectories. V2 closes its
two named gaps without rewriting it. [Confirmation driver](../scripts/derive_nsc_nested_parent_child_confirmation.py),
[portable basis](../results/development/nsc-nested-parent-child-replay-basis-v1.json),
[independent replay](../scripts/check_nsc_nested_parent_child_confirmation.py),
and [figure](../results/development/nsc-nested-parent-child-figure.png).
The completed article and release scopes below stay unchanged.

## Present spherical loop

The conformal realization now measures the finite chain through $T=0.3$:
source differences, surface exchange, geometric response, maintained
localization, and a local response on the same generated trajectory.
The [owning note](https://doi.org/10.5281/zenodo.23090015) connects
the source controls, actual metric curvature, matched clock, and local
reduction. Reported numerical movements meet one percent of their own
claimed effects. They are refinement indicators, while stronger propagated
state-error and continuum certificates remain separate claims. Renewal is
false. The earlier prescribed-gauge records below keep their original domains.

The [14-page focused article](../paper/finite-regeneration.pdf), dated
1 October 2026, presents this measured finite result. Its
[source](../paper/finite-regeneration/main.tex) and
[manifest](../paper/finite-regeneration-manifest.json) bind the pinned
scientific evidence; human author review and the actual arXiv processor
check remain pending. The [paper guide](papers.md) preserves the foundation
and earlier focused edition.

The code owner is
[nsc_spherical_galerkin_coupling.py](https://doi.org/10.5281/zenodo.23090015).
The v5 replay,
[derive_nsc_spherical_galerkin_refinement_v5.py](https://doi.org/10.5281/zenodo.23090015),
remains the self-contained $T=0.005$ seed. The saved $T=0.05$ measurement
is [derive_nsc_spherical_feedback_episode.py](../scripts/derive_nsc_spherical_feedback_episode.py),
which reuses that preparation. The weak residual and the streamed reducer
are separate method owners, linked below.

`rk4_step` calls `rates` on each stage state, and `rates` calls
`compose_fine_hamiltonian`. That call prolongs the stage state, recomputes
the column source, and recomputes the geometric rates. With the matter force
included it feeds `force_Q` into the conjugate momentum,

$$
\dot p_Q \leftarrow \dot p_Q - F_Q/\Delta x_Q.
$$

$Q$, $r$, $\chi$, the momenta and the spinor columns evolve. $L$ and $\beta$
stay gauge controls. `solve_initial_radius`, and the v5 replay's full-source
Newton (`full_source_rho`), set the initial radius from that source,
$\rho = F_L/\Delta x_Q$.

The descriptions above refer to the historical default `prescribed`
gauge. The same owner also exposes the explicit conformal realization
below, whose dynamic lapse requires its additional chain-rule force and
$F_L\dot L$ work term.

<a id="same-action-conformal-episode-1-october-2026"></a>

### Same-action conformal episode, 1 October 2026

The [proof](https://doi.org/10.5281/zenodo.23090015) varies the same
action before choosing $L=Q$, $\beta=0$. The
[record](../results/development/nsc-spherical-conformal-episode-v1.json)
and [payload](https://doi.org/10.5281/zenodo.23090015)
store four runs from the unchanged v5 initial Cauchy arrays, with
$n_f\in\{256,512\}$ and $dt\in\{0.001,0.0005\}$. All reach $T=0.05$.
CPU is $120.662$ s and payload size about $4.07$ MB. The source,
preparation, initial-radius model, and eight initial fields are unchanged.

The finite constraint pair $h_c=QC$, $D$ satisfies
$R'\le\|f\|_2$ for its actual quadrature and projection forcing. On the
fine confirmation run, $R$ starts near $2.09546\times10^{-6}$ and ends
near $2.09543\times10^{-6}$. Its sampled forcing integral is
$7.30987\times10^{-5}$; no sampled budget violation occurs. The integral
is checked by timestep refinement and Simpson/trapezoid movement, with
`forcing_time_integral_certified=false`. It is a finite constraint budget,
not an observable-error certificate.

Field and gravitational energy exchange about $1.41447\times10^{-7}$,
with total-energy change about $-2.20\times10^{-12}$. Dynamic work includes
$F_L\dot L$, with $\dot L=\dot Q$. The normal clock is explicitly the
fixed worldline $x=1$, $d\tau=rQ\,dt$; it advances by about $0.0603183$.
The original fixed observer's occupation changes by about $-0.0842690$.
At equal proper time the coupled-minus-frozen occupation is about
$-2.16786\times10^{-6}$. Its largest timestep movement is about $0.379$
percent of that small effect. Field exchange moves by about $0.0477$
percent under timestep halving, and the largest spatial movement of the
compared radius/velocity changes is about $0.0287$ percent. The chart and
Gram stay admissible. Renewal remains false; the regional ledger and
constraint budget keep their own measurement domains.

### Finite conformal chain through $T=0.3$, 1 October 2026

The [continuation](../results/development/nsc-spherical-conformal-episode-v2.json)
reuses each exact $T=0.05$ final and reaches $T=0.3$ on all four cases.
The [surface assessment](../results/development/nsc-spherical-conformal-transport-v2.json)
resolves absolute integrated flows $0.03937795$, $0.07025487$, and
$0.06351828$ at $x=0,2,4$; their largest frame indicator is $0.095\%$.
Flow at $x=6$ remains unresolved in its own domain. Entry and exit remain
separate from their net contribution. Region 0 keeps shell and probability
shares above $0.65098$ and $0.65585$. Boundary flux, pressure and lapse work
explain its shell change, with a $0.1641\%$ trapezoid indicator.

The [source controls](https://doi.org/10.5281/zenodo.23090015)
prepare each changed occupation source through the owned radius/shift solve.
Reversal reverses signed $x=2$ flow and moves the maintained shell leader;
$\pm5\%$ imbalance changes it by approximately $\pm5.04\%$. Their
exploratory coarse and frame-comparison domains remain stated. The
[metric measurement](https://doi.org/10.5281/zenodo.23090015)
uses the actual projected rate Jacobian, independently checked, yielding
maximum Weyl invariant about $0.11331608$ at $T=0.3$. The
[analytic clock comparison](../results/development/nsc-spherical-conformal-clock-v2.json)
removes the old frozen-control interpolation: its small occupation separation
has timestep movement below $0.00079\%$ of the effect.

The [local production](nsc-coupled-local-response.md)
uses actual $\Phi(0.2)$ and the original $T=0$ observer on $[0.2,0.3]$.
Its occupation changes by $0.07342952$. Streamed/full discrepancy is
$6.3508\times10^{-6}$, or $0.00865\%$ of that effect. Memory, drive, and
actual cross omission move occupation by $0.01681854$, $0.04854120$, and
$0.09206143$, each tested against its own numerical comparison. No dense
exterior propagator is stored. Source dependence, transport, geometry,
maintained structure, and local response therefore coexist in one finite
realization. The [note](https://doi.org/10.5281/zenodo.23090015)
owns its domains and read-only replays. The tracked finite forcing budget
and initial-source bounds are reported; a certificate for an unknown
continuum or geometric path is a stronger claim, not the finite finish line.

### Recorded $n_f=512$, $T=0.005$ diagnostic

The [v5 record](../results/development/nsc-spherical-coupling-refinement-v5.json)
`nf512.window_dt_0_0005` was independently replayed. The stored window maxima
are

| Diagnostic | Recorded maximum |
|---|---:|
| Full Hamilton, `full_C_hamilton_max` | $9.981468739539423\times 10^{-6}$ |
| Full momentum, `full_D_momentum_max` | $1.3753934746951746\times 10^{-6}$ |
| Chart proper speed, `chart_proper_velocity_max` | $0.007950515372629119$ |
| Lifted proper speed, `lifted_normal_velocity_max` | $0.007950515381253666$ |

In short: full Hamilton about $9.98\times 10^{-6}$, momentum about
$1.38\times 10^{-6}$, maximum absolute proper radial velocity about
$0.0079505$. Both velocity entries are maxima of absolute values. Signed
radial motion is not inferred from them. The stored verdict is
`DIAGNOSTIC_MEASURED`, and `time_extension_T_0_05` is false. Renewal is
not demonstrated. The historical initial tolerance $10^{-8}$ is not an
automatic veto of this diagnostic. A weak residual alone would not be
sufficient; these full residuals stay reported.

<a id="recorded-t005-feedback-episode"></a>

### Recorded $T=0.05$ feedback episode

The [episode record](../results/development/nsc-spherical-feedback-episode-v1.json)
and [payload](https://doi.org/10.5281/zenodo.23090015)
store four runs, $n_f\in\{256,512\}$ and $dt\in\{5\times 10^{-4},2.5\times 10^{-4}\}$.
Each reaches $T=0.05$. Campaign CPU is $179.444027$ s. The verdict is
`MEASURED_FEEDBACK_UNRESOLVED_CONSTRAINT_CONTROL`. Renewal is false.
`constraint_consistent_solution` is false. Continuum error is unknown.

On the fine primary run, `nf512_dt_0_0005`, the final proper radial
velocity runs from $-0.06195632555218566$ to $0.0881593239017512$. Final
$Q_{\min}$ is $0.2195220548716644$. Final $\chi_{\max}$ is
$11.148261641496925$. Field energy changes by $-0.4123021259213786$ and
gravitational energy by $0.4123020914983755$. Their sum changes by
$-3.442300311462532\times 10^{-8}$, which is
$8.348975411586745\times 10^{-8}$ of the field exchange. That comparison
is against the exchange, not one percent of the near-zero balance.

All 62 time and space comparisons of the 31 physical observables meet the
one-percent test. The largest relative movement is
$5.546186987854776\times 10^{-4}$, on `observer_K_perp_min` in space.
The one-percent target applies to those claimed physical effects. The
assessment metadata was rewritten from the saved four-run arrays. The
evolution was not rerun, and the historical v5 JSON and NPZ bytes are
unchanged.

| Constraint diagnostic | $n_f=512$ final | $n_f=256$ final |
|---|---:|---:|
| Full Hamilton | $2.825670750915076\times 10^{-5}$ | $0.0068275893562912415$ |
| Full momentum | $4.294300980589583\times 10^{-6}$ | $0.001504909818827671$ |

Full, projected, and held-out residuals remain. The coarse Hamilton
residual stays far above the fine residual. No real bound validates a
constraint-consistent solution. Those residual changes do not veto the
physical effects. The historical initial tolerance $10^{-8}$ did not
abort the run and is not an automatic veto.

Owned tests are
[test_nsc_spherical_feedback_episode.py](https://doi.org/10.5281/zenodo.23090015).
Independent tests are
[test_nsc_spherical_feedback_episode_independent.py](https://doi.org/10.5281/zenodo.23090015).

The historical episode's stored Galerkin pin differs from the current file.
The current producer's strict `verify_saved` therefore refuses that replay.
The [historical reader](../scripts/check_nsc_spherical_feedback_episode_historical.py)
authenticates the archived source bytes, frozen inputs and saved arrays;
it does not import or re-evolve the old module because four transitive
imports were not part of the recorded frozen set. The scientific numbers
above retain their historical domain. This authentication is distinct from
a numerical replay under a complete historical execution context.

<a id="regeneration-episode-1-october-2026"></a>

### Regeneration episode, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015),
[record](../results/development/nsc-regeneration-episode-v1.json), and
[payload](https://doi.org/10.5281/zenodo.23090015) store
four continuations, $n_f\in\{256,512\}$ and
$dt_{\max}\in\{0.0005,0.00025\}$. Each begins with the predecessor's
`dt=0.0005` final Cauchy arrays, unchanged, at $T=0.05$. Each stops at
$T=0.085$ on `BRIDGE_ARC_ENTERED_PACKET`: the both-positive arc's left
edge reaches $x=3.96875$ inside the packet $[0,4)$. The arc still spans
the bridge. All four refinement pairs resolve the claimed physical
changes within one percent.

Window 0 stays the leader. Its share rises from about $0.60870$ to
$0.61653$ and its shell content by $0.08117$. The packet proper flux ends
near $-4.98756$, balanced by the reservoir partition. The normal-window
balance also includes pressure and lapse work; its integrated residual is
under about $10^{-5}$ of the exchange. Field energy changes by about
$-0.28608$, with a compensating gravitational change. Proper time advances
by about $0.10673$ and the leader clock by $0.06438$. The chart stays
positive, the column Gram gap stays near $10^{-9}$, and the largest
accepted $dt\,\omega$ is about $1.201$, below the declared $1.4$ cap.

`maintained_structure` is true and `renewed_structure` is false. The
joined trajectory is one drift with no reset. It carries frame $p_\chi$,
$\dot Q$, and a realized one-step rate increment. The
[assessment](https://doi.org/10.5281/zenodo.23090015) reads those
saved arrays for the metric/proxy comparison and proper-arc localization;
the increments and refinement movements remain numerical indicators.
The [independent audit](https://doi.org/10.5281/zenodo.23090015)
separates the handoff, covariance, ledger, geometry controls, and response.
No propagated observable-error bound is supplied. Full, projected, and
held-out constraints stay in the record; at $n_f=512$ the final full
Hamilton maximum is about $1.396\times10^{-4}$, and at $n_f=256$ about
$9.241\times10^{-3}$.

### Weak initial residual and initial-chart certificate

The [weak note](https://doi.org/10.5281/zenodo.23090015),
[helper](https://doi.org/10.5281/zenodo.23090015)
and [record](https://doi.org/10.5281/zenodo.23090015)
assess the saved seed. The sealed v1 status is `PARTIAL_TERMS_NO_TOTAL_BOUND`. No state
is stepped. Declared $G>0$ with `rho_independent` exactly true exposes
$\mu=\min(4,Q^2)=Q^2=0.0631642220827373$ without $y_*$. An omitted
independence flag does not establish the hypothesis. The sealed diagnostic
did not provide the full residual dual, a continuum positivity enclosure
of $G$, or a rounding enclosure.

On $n_f=512$ the dense owner maximum moves from
$9.981455056262689\times 10^{-6}$ at $n_q=2048$ to
$4.539313647278211\times 10^{-5}$ at $n_q=4096$. The Fourier product
$P=2rr''-(r')^2-Q^2 r^2+G$ stays from
$5.068517579063356\times 10^{-6}$ to
$5.400936559883778\times 10^{-6}$. That movement is a conditioning
indicator. `rounding_enclosure` is null. The v5 bytes are unchanged.
The saved record's `checkpoint_head` is `9a9090a`. `--check` and
`verify_saved` remeasure the scientific fields from the declared v5 bytes
and do not write. A live HEAD other than `9a9090a` is execution metadata.
The sealed generator hashes of the helper, the owned test, and the driver
differ from the current files, so replay reports those three limits. That
report is not a reexecution of the original producer. `total_certified`
stays false. The only write is an explicit `--output` path, and that path
cannot be the sealed JSON.

The separate [initial certificate](../results/development/nsc-spherical-cauchy-error-v1.json)
has status `FINITE_CONDITIONAL_INITIAL_LAPSE_BOUND`. On the saved fine
initial chart it encloses $G\ge0.0537330$ and a pointwise radius distance
of at most $4.06085\times10^{-6}$. The owning note proves existence and
uniqueness of the smooth positive critical radius under those hypotheses,
then a nearby exact initial state using the same-action shift profile.
The saved state retains $p_r=0$; that correction is not installed. The
initial radius comparison is about $0.2196$ percent of the smaller saved
radius effect. `evolution_error_bound` remains null. The helper's
constraint-transport and source-normalized residual budgets are separate
from a radius or occupation error estimate.

### Streamed retained-region reduction

[nsc_evolving_reduction.py](https://doi.org/10.5281/zenodo.23090015)
keeps the dense prescribed-control owner and adds `backend="streamed"`.
The streamed path matches the dense history and the independent exterior
sample on the same prescribed controls. It does not store a time-indexed
exterior propagator. Those prescribed-control checks do not consume a
coupled $H(g(t))$ trajectory, and the sine schedule remains a prescribed
$H(t)$. The conditional local response below calls the same backend on
stored episode frames and does not edit this reducer.

The toy $H$ and $B$ fractions are an optional witness.
Geometry-derived links are the target. Embedding v1 and the seam-regular v2
compression are an optional bridge. Neither that witness nor the historical
$10^{-8}$ tolerance sets the objective.

The earlier controls and consumers are recorded below in their original
domains. Renewal and continuum constraint control remain open. The completed checks have named
domains. The controls successor reads the stored $T=0.10$ endpoint and
leaves the v1 evolution unchanged. The local-boundary review checks stored
payloads, the projected endpoint, and a small causal ODE; the saved
$1024$-dimensional series remains that record. The same review binds the
immutable null-expansion and local-response files to their audited owners
at commit `5f10ecd`, replayed explicitly as historical sources.
The initial geometry-graded window is a completed measurement of the
initial packets. The records are
[regeneration controls](#regeneration-controls-1-october-2026),
[controls successor](#regeneration-controls-successor-1-october-2026),
[conditional local response](#conditional-local-response-1-october-2026),
[spherical null expansion](#spherical-null-expansion-1-october-2026),
[local-boundary review](#local-boundary-review-1-october-2026),
and
[initial geometry-graded window](#initial-geometry-graded-window-1-october-2026).

This loop is the present finite realization because its geometry and its
state forces already come from the same action. An invented rule for the
radius as a function of energy is unnecessary. The recorded outcome
guides the next model step. A candidate may be chosen creatively. The
computation and its error statement stay checkable.

The chart identity
$\dot Q=L p_\chi/(2F_\chi)+\partial_x(\beta Q)$ is now carried by the
continuation's saved frame rates. The metric assessment consumes its
realized increments; it does not substitute $\chi+2$ for an independent
$R_h$. The remaining uncertainty is an evolved error estimate at the scale
of the claimed geometry, transfer, and local-response changes. Constraint
smearings and refinement agreement provide diagnostics with named norms;
they do not by themselves bound trajectory or observable error.

The source-fixed incoming gate remains OPEN and its campaign remains paused.
It is not a prerequisite for this loop. This page does not change either
frozen PDF and does not assign the next release version. Frozen v0.26.0
remains the foundation.

The dated sections from finite circulation through the initial
geometry-graded window are evidence for this programme.
[Nested qualities](#nested-qualities-finite-mathematical-construction-27-september-2026)
records the finite theorem. From
[the incoming gate](#source-fixed-local-incoming-gate-open-campaign-paused)
through the Stage-2 notes, the text is that retained application.
[Earlier public checkpoints](#earlier-public-checkpoints) are history.

## Finite circulation and local response, 30 September 2026

A six-mode inherited finite example has nonzero channel currents, constant
regional contents, and a derived local filtered response. The
[note](https://doi.org/10.5281/zenodo.23090015) records it. The geometry and the
stationary preparation are fixed inputs. That prepared-current witness remains
valid for that fixed window. The present evolution is the spherical loop above.

## Incoherent imbalance mean circulation, 30 September 2026

On the same frozen six-mode window,
$C_0=I/2+(1/4)(P_0-P_2)$ starts with zero interregion channel current.
The executed Cesaro mean still carries inflow
$j_{(0,+),(1,+)}=0.0012511140479781876$, into $0+$ from $1+$. The
[note](https://doi.org/10.5281/zenodo.23090015) is the scientific owner.
The [record](../results/development/nsc-imbalance-turnover-v1.json)
replays 34/34 checks, with verdict `PASS_IMBALANCE_MEAN_CIRCULATION`.
Seven focused tests and the six existing finite-turnover tests passed, and
`--check` replayed the same four source bindings. An independent review
accepted the 34 conditions, the finite-time bound, and the forced
optional-failure semantics. Trace 3 and energy 99/16 are conserved;
instantaneous regional populations are not fixed. Cesaro averaging is not
relaxation. The channel current is not by itself heat or energy flow. The
filtered covariance is not a stress. Geometry remains an input in this
witness. The average does not by itself give autonomous regeneration.

## Spherical feedback action and direct conformal source, 30 September 2026

The [spherical feedback action](nsc-spherical-action.md)
is a checked first-order chart for the spherical Einstein–Maxwell–Weyl bulk
and boundary. Lapse and shift stay free. Seventeen tests passed. The module
sha256 is `708df11f3a9b69a32b28148298679098f2c04523edc329c2c4c836fe1ed0e24f`.
The check covers the auxiliary bulk, periodic $x$ with Dirichlet data
$Q$, $r$, and $\chi$ on the time caps, cancellation of the Euler and box
cap ledger, and an extra Weyl cap that does not double the Einstein–Hilbert
Gibbons–Hawking–York term. Six Hamilton right-hand sides match the chart.
Conditional nonuniform initial data solve the constraints. This is not a
coupled Cauchy problem and not a renewal. $C_W=0$ is rejected. The matter
source is not coupled.

The [direct conformal ADM source](https://doi.org/10.5281/zenodo.23090015)
is a separate owner. Ten tests passed. The module sha256 is
`7e935b981de6cb4876142819a1a15587b81574ef1c15bf11adc62a7d8cdaeaea`.
It fixes the finite-grid conformal Ward identity and leaves the historical
$\sqrt{N}$ source bytes unchanged. The result is raw nodal canonical forces
with one $4\kappa$ weight. It is not an absolute renormalized stress. Frame
mapping is absent.

Neither primitive embeds the frozen six-mode generator. A direct copy of
those sites into this chart is obstructed. A localized countercarrier
compression matches the onsite block $H$ and does not match the link $B$.
The full mapping is pending. The comparison residuals are not frozen data.

The vacuum-matched CTP branch remains the separate historical choice in
[its own note](nsc-local-observer-correspondence.md). These records are
equations and interfaces for the loop above.

## Provisional spherical coupling v1, 30 September 2026

The [note](https://doi.org/10.5281/zenodo.23090015) owns the failed v1 control.
The [record](../results/development/nsc-spherical-coupling-control-v1.json)
sha256 is
`88bea1c962478a2ec57e0c0f3844ab6ee9237d148dc1863f10122ee8fa5a9b60`.
The checkpoint that wrote it had module sha256
`8727a46aa7785a326e301eeeee370a77aeaa066ddb1dbaee022a255616887dbf`,
test sha256
`379fc239a1a3279a842db48e92d41c5be326ac85f616c773f420cd94ffa28083`,
and note sha256
`bd84c977965922694f8514db604476353375e1741e9966f5508e97d19a1d2ec7`.
The current source is the minus-column phase correction. Its module sha256
is `64e066b9b3210121461f902e748d4fc8a9cfc5978af2af47ab4b12fd8797e4f5`,
its test sha256 is
`da3f0b0a50e61487de017e04311acba0e90781d759ed41a4a174dbaf2a97f163`,
and the note sha256 is
`a0b7e8207c6ef6f0d17addd549d8cd4d8e5f02dcacd62449714d78b2ec700499`.
The saved verdict is `PROVISIONAL_CONSTRAINT_DRIFT`. Ten tests in the
current file passed. Seven require the record to report this failure.
Three check the column phase and the refusal to overwrite the v1 JSON.
They do not make the physical trajectory valid.

Initial N=64 constraints pass: the Hamilton residual maximum is
$9.556799795973347\times 10^{-12}$, and the field energy
$24.7195086$ cancels the gravitational energy, with total
$-1.4\times 10^{-14}$. By $T=0.005$ the N=64 Hamilton maximum is
$36.98886722265645$. Halving the timestep leaves
$36.988867326189634$. The same short interval gives
$2.330762882869135$ at N=128 and $0.575679312602035$ at N=256. N=256
still fails the validation tolerance $10^{-3}$. This is not a
self-consistent solution and not a renewal.

The calibration phase in the v1 JSON is on the odd spatial lobe. That
preparation error does not select the constraint drift. The current source
puts the recorded phase on the minus column and does not rewrite the v1
JSON. These v1 bytes stay checkpointed. Collocation constraint algebra,
sha256
`eefadb1bdf863da6c1e487304fe816bc1ae6514a6570b85431893b7d97b9db78`,
does not load the saved trajectory and does not certify the Galerkin
discretization.

The [evolving reduction](nsc-evolving-reduction.md) is a separate
prescribed-control owner. The module sha256 is
`17bd59363a25a8afba139b6cc6737f57684af5105a7ef148e68047316088e59d`.
The dense-test sha256 is
`2473425909c3fb67f786daa482e2cadb7ddf7052975976d51d1e81098fe68f97`.
The note sha256 is
`e679c08abd76278e08d2fa49ab41cb21cc93da25e16c74a670f7c72ad4f56f6e`.
The streamed-test sha256 is
`9debf8b7b11b676f83e2577939f95896926296f6f36a2ad910da6346329e3edf`.
Thirteen dense tests and nine independent tests, sha256
`bba311f6ee90acd2fecf631b0e5f5cb193a98b6726ae6cf4b1da2e56a67fe5d6`,
use one prescribed generator and an exterior integrator. Eleven streamed
tests check the same controls without a stored exterior propagator. The
sine schedule is a prescribed $H(t)$. It is not an autonomous trajectory,
not regeneration, and not a coupled geometry. No actual $H(g(t))$ was
consumed by these tests, and no coupled trajectory is accepted by them.
The conditional local response below calls the streamed backend on stored
episode frames and does not edit this reducer. Those tests do not
validate these failed coupling trajectories. The vacuum-matched CTP branch
remains a separate historical choice.

## Variational Fourier–Galerkin coupling v2, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015) and the
[record](../results/development/nsc-spherical-coupling-control-v2.json)
own the failed subspace control. The saved verdict is
`FAIL_HELD_OUT_CONSTRAINT`. Renewal is false. $T=0.05$ was not run.
A projected residual is not a physical pass. This is not a filter of the
v1 radius, not a constraint projection, and not a renewal. The v1 JSON
is unchanged. The Dirac compression of that window is recorded below and
is not an accepted coupled trajectory.

## Matched Galerkin refinement v3, 1 October 2026

The same [note](https://doi.org/10.5281/zenodo.23090015)
records the failed matched refinement. The
[driver](https://doi.org/10.5281/zenodo.23090015) and
the [record](../results/development/nsc-spherical-coupling-refinement-v3.json)
are immutable. The verdict is `FAIL_INITIAL_SOLVE`. Internal Newton
stopped at rho $=0.5$ with residual $1.2876060991167562\times 10^{-10}$.
No full-source result was produced. The v1 and v2 JSON bytes are
unchanged. Renewal is false. $T=0.05$ was not run for this refinement.

## Galerkin refinement v5, 1 October 2026

The [driver](https://doi.org/10.5281/zenodo.23090015),
[record](../results/development/nsc-spherical-coupling-refinement-v5.json),
[payload](https://doi.org/10.5281/zenodo.23090015)
and
[independent tests](https://doi.org/10.5281/zenodo.23090015)
are the completed $n_f=512$, $T=0.005$ geometric diagnostic summarized in
[the present loop](#present-spherical-loop). `compose_fine_hamiltonian`
recomputes the source and the rates and feeds $F_Q$ into $\dot p_Q$.
$Q$, $r$, $\chi$, the momenta and the spinor columns evolve. $L$ and
$\beta$ stay gauge. The initial radius uses the column source. The recorded
proper radial velocity is a maximum of the absolute value, about $0.0079505$.
Signed motion is not inferred from it. The saved $T=0.05$ episode and the
weak assessment are the present-loop sections above. This v5 record does
not itself contain $T=0.05$. The historical initial tolerance
$10^{-8}$ is not an automatic veto. A weak residual alone is not sufficient.

## Finite-window embedding v1, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015),
[module](https://doi.org/10.5281/zenodo.23090015),
[record](https://doi.org/10.5281/zenodo.23090015)
and [independent tests](../tests/test_nsc_finite_window_embedding_independent.py)
compress the frozen window onto the existing Dirac operator. Quadrature
matches the original $H$ and $B$ at about $4.45\times 10^{-14}$; a direct
analytic integral agrees at about $4.9\times 10^{-14}$. Modes are supported
on $(n,n+2)$. The subspace is not invariant, outside coupling is strong, and
the closed-window evolution is not transferred. Derivative jumps at the lobe
joins make the Fourier error decay as $N^{-2}$. This toy $H$/$B$ match is
an optional witness. Geometry-derived links are the target.

## Seam-regular embedding v2, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015),
[module](https://doi.org/10.5281/zenodo.23090015),
[tests](../tests/test_nsc_finite_window_embedding_smooth.py) and
[record](https://doi.org/10.5281/zenodo.23090015)
are the optional seam-regular bridge. The saved status is
`SEAM_REGULAR_COMPRESSION_NOT_INVARIANT`. Renewal is not claimed.

## Spherical Cauchy data, 1 October 2026

The [module](https://doi.org/10.5281/zenodo.23090015),
[note](https://doi.org/10.5281/zenodo.23090015),
[tests](https://doi.org/10.5281/zenodo.23090015) and
[record](../results/development/nsc-spherical-cauchy-data-v1.json)
provide a $p_Q$ antiderivative for $\chi=p_r=p_\chi=0$ and constant $Q$.
Three tests pass. The mean and held-out current residual stay explicit.

## Regional coordinate and normal-observer ledger, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015),
[module](schur-memory.md),
[record](../results/development/nsc-regional-energy-exchange-v1.json),
[tests](schur-memory.md) and
[independent tests](schur-memory.md)
are a diagnostic. Six owned tests and eight independent tests. Normal energy
is $F_L/r=V\rho\,dx$, and the proper pressures divide by $dx$. $F_Q\dot Q$
is coordinate metric work, not observer pressure work. Hamiltonian
derivatives are not pressures. Historical vacuum branch $B$ is not a
prerequisite. Renewal is not claimed.

## Regeneration controls, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015),
[module](https://doi.org/10.5281/zenodo.23090015),
[driver](../scripts/derive_nsc_regeneration_controls.py),
[record](../results/development/nsc-regeneration-controls-v1.json)
and [payload](https://doi.org/10.5281/zenodo.23090015)
are a thin consumer of the owned Galerkin step and the saved episode.
The column source, radius Newton, shift momentum, and `rk4_step` stay
with their owners. No radius force is added, and the episode driver is
not edited. Record status `MEASURED`. `maintained_structure` is true,
`renewed_structure` is false, and `one_drift` is true. `candidate_regime`
is false. That flag, and the note's comparison floors, belong to this
record. They are not programme requirements.

The saved baseline was not rerun. Frozen geometry keeps the saved initial
state and sets geometry velocities and applied fieldwork to zero. On that
control the $0.412$ field-energy exchange, the change in $\chi$, and the
proper-velocity change are absent, while the shell still rearranges. The
frozen trajectory drives the full Hamilton residual to about $24$ and is
not a constraint solution. Uniform occupations stay below a localized
share. Reversed occupations move the leader to the second window and
increase the field-energy exchange. $n_f=256$ reproduces those $T=0.05$
effects inside one percent.

The second episode starts from the saved `nf512_dt_0_0005` final arrays,
hash
`cd6b1f88a12635b137aaddc47444c06cc681c7568271ef989312f48f43a0649e`,
and runs to $T=0.10$ without a reset. Window 0 remains the leader. Its
share moves from $0.6087043074743482$ to $0.6198826640133578$, and its
shell content rises by $0.11695656049787839$. The joined series from
$T=0$ has reversal $0$. Packet proper flux ends at $-5.534007341562886$
and reservoir flux at $5.534007341562912$; their sum is
$2.5757174171303632\times 10^{-14}$. The sign does not reverse. Field
energy changes by $-0.395776510604783$ and gravitational energy by
$0.3957764909408823$. The same drift continues: $\chi_{\max}$ grows from
about $11.15$ to about $89.31$, $Q_{\min}$ falls from about $0.220$ to
about $0.125$, and $G_{\min}$ stays positive at about $0.0156$. The
$n_f=256$ leader-content change agrees to $6\times 10^{-10}$. That
agreement does not change the conclusion.

This is maintained localization of one finite drift. It is not renewed
structure and not autonomous regeneration closure. Continuum constraint
control remains unproved. The record leaves open whether $G$ or $Q$
reaches the chart boundary at a later finite time. It does not claim a
continuum limit or a cosmological regeneration.

<a id="regeneration-controls-successor-1-october-2026"></a>

## Regeneration-controls successor, 1 October 2026

[nsc-regeneration-controls-v2.json](../results/development/nsc-regeneration-controls-v2.json)
assesses the immutable v1 JSON and NPZ. Status `ASSESSED`.
`dynamics_rerun` is false. The independent checks are
[test_nsc_regeneration_uniform_independent.py](https://doi.org/10.5281/zenodo.23090015).
$G>0$ is the initial convex-radius ansatz for constant $Q$, vanishing
momenta, and $r=y^2$. The dynamic chart is `chart_failure` on positive
$r$, $Q$, and $L$. On the $n_f=14$, $n_q=64$ control the Galerkin
half-identity commutes with $H_G$. The same identity embedded in the
unprojected fine operator does not. At the stored $T=0.10$ endpoint,
$dt\,\omega=1.4994122493748205$, above the owned cap $1.4$ and below
the absolute RK4 limit $2\sqrt{2}$. The column Gram gap is
$1.3069428872469757\times 10^{-9}$. $r$, $Q$, and $L$ stay positive,
and `chart_failure` is empty. The owned stable timestep is
$4.668495940938623\times 10^{-4}$. A missed reversal proxy, a missed
leader-share proxy, and the cap label leave the maintained-drift reading
in place.

<a id="same-trajectory-local-response-successor-1-october-2026"></a>

## Same-trajectory local-response successor, 1 October 2026

The [v2 record](../results/development/nsc-coupled-local-response-v2.json)
and [series](https://doi.org/10.5281/zenodo.23090015)
consume the fine confirmation continuation from $T=0.05$ to $T=0.085$.
Status `MEASURED_SUCCESSOR_RESPONSE`. The fixed observer remains the two
original $T=0$ mode columns. At the handoff its retained eigenvalues are
about $0.41940$ and $0.64745$, and the actual initial cross block has
Frobenius norm $0.44224$. The initial correlations are carried into the
reduction.

The saved autonomous occupation changes by $0.26748893$. Conditional
full evolution on the stored linear $Q(t)$ differs from the autonomous
frames by at most $1.16662\times10^{-7}$ in occupation. The streamed
reduction differs from that conditional full series by $0.0002570304$,
about $0.0961$ percent of the occupation change. Output, exterior-step,
Hermite-schedule, and sibling-geometry movements stay below one percent
on their measured windows.

Omitting memory moves the occupation by $0.06952775$, omitting the
exterior initial drive by $0.15526034$, and omitting only the actual
initial cross correlations by $0.20434176$. Each exceeds the reduction
error. The saved full-window retained and cross-omitted covariance arrays
reproduce the last value directly. The local probability, proper clocks,
and normal-energy ledger keep their separate measures. No effective-action
kernel variation or stress is computed. The
[owning note](nsc-coupled-local-response.md) gives the source
bindings and the independently checked causal reduction.

The [independent supplement](../results/development/nsc-coupled-local-response-audit-v1.json)
propagates the same retained/exterior split with its own temporal and
covariance assembly. Cross omission moves occupation by $0.20437297$;
its streamed discrepancy is $0.0001195768$. Drive omission moves it by
$0.15497464$, with a separate streamed discrepancy of $0.0002856974$.
The memoryless recurrence reproduces the saved memory-off occupation to
roundoff. These errors are below one percent of their own effects.
The note also derives a conditional Duhamel relation from a future
Hamiltonian-path defect to local occupation error. It does not supply the
missing geometric-path bound.

## Conditional local response, 1 October 2026

The [note](nsc-coupled-local-response.md),
[module](https://doi.org/10.5281/zenodo.23090015),
[driver](../scripts/derive_nsc_conformal_local_response.py),
[record](../results/development/nsc-coupled-local-response-v1.json)
and [series](https://doi.org/10.5281/zenodo.23090015)
compare one streamed retained-region reduction with geometry the episode
already stored. The reducer is called with `backend="streamed"` and is
not edited. The episode is not run again. Record status
`MEASURED_CONDITIONAL_REDUCTION`; the stored-frame phase is
`MEASURED_STORED_FRAMES`. Renewal is false.

`phases.full_window` uses all eleven stored frames of `nf512_dt_0_0005`,
from $t=0$ through $t=0.05$, with piecewise-linear fine $Q$ and two
exterior substeps. The conditional full-band occupation of the fixed
region-0 observer changes by at most $0.3305919037025379$. The streamed
occupation differs from that series by at most
$0.0007353273811995242$. Saved autonomous spinors exist only at the
initial time and at $T=0.05$. Projected onto the same observer, the
conditional full-band endpoint differs from that saved projection by
$3.4123144243558556\times 10^{-7}$. The streamed endpoint differs from
the saved projection by $0.0007356686126419598$, the same gap as the
Volterra discrepancy against the conditional full evolution.

On that declared geometry, omitting memory separates the occupation from
the conditional full series by $0.36709975607817125$. Omitting the
exterior initial drive separates it by $0.028764901996227665$. Both
exceed the reduction error. A nodal $\dot Q$ series is not stored. No
stress is derived, because no effective-action variation was computed.
The run does not regenerate the spherical trajectory and does not close
the incoming gate. The pilot window $[0,0.035]$ remains the earlier
phase of the same record.

The [local-boundary review](#local-boundary-review-1-october-2026) checked
the stored payloads, the projected autonomous endpoint, and a small
independent causal ODE. The $1024$-dimensional series remains the saved
record. A nodal $\dot Q$ series is absent. Hermite interpolation is an
indicator of the midpoint geometry. No stress is computed. The live
observer label is initial mode columns 0 and 1, weights $0.75$ and $0.75$.
The stored v1 JSON keeps the observer string
`initial region-0 plus and minus source columns, unphased`.

## Spherical null expansion, 1 October 2026

The [note](nsc-spherical-null-expansion.md),
[module](https://doi.org/10.5281/zenodo.23090015),
[driver](../scripts/derive_nsc_spherical_null_expansion.py)
and [record](../results/development/nsc-spherical-null-expansion-v1.json)
are kinematic postprocessing of the stored episode. The production step
is not rerun. The measurement does not evolve a state, rebuild the
column source, or recompute inheritance, and it does not read $\chi$.
Status `MEASURED_SAMPLE_TRAPPING_CHARACTER_CHANGES`. The episode JSON
and NPZ bytes are unchanged.

At $T=0$ every quadrature sample is untrapped, and the expansion product
is negative on every node. The initial areal maximum and minimum already
lie in the chart. Where the stored slope vanishes, both expansions equal
$2v_n/r$, and those critical values sit inside the initial
normal-velocity noise, so the extrema are marginal. From the first
stored frame, $T=0.005$, through $T=0.05$, each run has one
both-negative arc at the areal maximum and one both-positive arc at the
areal minimum. Complementary nodes stay opposite in sign. On the primary
grid of 2048 nodes the final counts are 45 both-negative, 182
both-positive, and 1821 opposite-sign. At $T=0.005$ those same-sign
counts are 4 and 8. Halving $dt$ leaves the raw counts unchanged.
Shared $n_f=512$ and $n_f=256$ nodes disagree in sign at $0$ samples
on all eleven frames.

The sampled labels are a local property of the areal gradient. A global
event horizon, an interior, an exterior, black-hole birth, a child
region, and regeneration are outside the measurement. Still open: the
instant inside $(0,0.005)$ when the product first becomes positive; the
continuum coordinate of each crossing inside one grid spacing; the sign
of $2v_n/r$ at the initial critical points inside the normal-velocity
noise; any global causal reading; any time after $T=0.05$; any
identification of these expansions with $\chi$.

<a id="local-boundary-review-1-october-2026"></a>

## Local-boundary review, 1 October 2026

[nsc-local-boundary-review-v1.json](../results/development/nsc-local-boundary-review-v1.json)
stays the sealed correction record. Its `reference_unchanged` hash is
`68c2fa1559d705c3aae2d3b64d693af72311802e51cb5914432e1422063b4ea9`, and its
recorded null-expansion module hash is
`1aa47f0b3d680621d23ac9329e0d9e1e38070b7df42fc2ec48848ea367fd82d5`.
[nsc-local-boundary-review-v2.json](../results/development/nsc-local-boundary-review-v2.json)
binds those immutable v1 payloads and the owners audited at `5f10ecd`. The
independent test at that audit is
`d2d6afb8471481a74b816223d58e33066f112abf1026d15eae74b291222e44fa`.
The two recovered files are pinned in `lab/.source-history` on carrier
commit `594a11caeee760d172c4a4f73a7615b05cf30282`. That carrier is a
content object, not a historical laboratory commit. The original eighteen
objects and four evidence pins stay in place. The v2 owner sources, including the notes and the null-expansion test,
are authenticated at full commit
`5f10ecd365843d1616e50eb16a20d7acd8377e2c` for historical replay. The
current binding path remains strict for new requests and never substitutes
history for a changed current source. The in-memory measurement is the sample replay.
The null-expansion writer still replaces the v1 JSON, so it is not the
replay. The review domain remains the stored payloads, the projected
endpoint, a small causal ODE, the boost-return mutation, and the
nonvacuous initial extrema. The $1024$-dimensional series, a kernel
stress, and a negative rescaling inside `boosted_expansions` stay outside
that correction.

<a id="initial-geometry-graded-window-1-october-2026"></a>

## Initial geometry-graded window, 1 October 2026

The [note](https://doi.org/10.5281/zenodo.23090015),
[module](https://doi.org/10.5281/zenodo.23090015),
[driver](https://doi.org/10.5281/zenodo.23090015),
[tests](../tests/test_nsc_geometry_graded_window.py),
and
[record](https://doi.org/10.5281/zenodo.23090015)
measure the initial three-packet compression. Status
`INITIAL_GRADED_COMPRESSION_NOT_INVARIANT`. The measurement is complete.
No trajectory is integrated, and the old rational $B$ is not a target.
The record's writing CPU is $2.248243$ s. `--check` recomputes the
comparison and leaves the JSON unchanged.

On the initial $n_f=512$ packets the grading defects are about
$5\times 10^{-14}$. The $n_f=256$ and $n_f=512$ blocks differ by
$1.339442162796698\times 10^{-8}$. The corner Frobenius residual is
$1.2001004173602996\times 10^{-14}$. Projection leakage reaches
$0.9658068315009373$, so `closed_six_mode` is false. The six columns are
not the Dirac evolution.

The published Schur recurrence uses $B_{01}$ with the factor $1/\Omega$.
The congruence $\Omega^{-1/2}B_{01}$ is a different link. Substituting
that normalized link as the block coefficient moves the top resolvent
entry by about $1.17\times 10^{-2}$. On the saved $T=0.05$ frame the
action functional is the same `first_order_density`. Evolved $Q$ breaks
the kinetic grading: diagonal defect $0.012665526439756825$, link defect
$0.04059408339655854$. The instantaneous configuration is not the initial
graded window.

## Nested qualities: finite mathematical construction, 27 September 2026

The [connected proof](nsc-nested-qualities.md) and [focused article](../paper/local-incoming-gate-draft.pdf)
now establish four properties for a declared finite coupled-operator family:
regional composition, normalized inherited law, normalized spectra and the
local response of unresolved regions. Different regional states are allowed.
The proof holds at every finite depth; it is not a proof of infinite physical time.

Twelve exact controls and four focused tests accompany the construction.
The [source snapshot](../paper/nested-quality-evidence/snapshot.json) binds its
proof, implementation and control record to an immutable laboratory commit.
No cosmological abundance, physical $c$ or $\Lambda$ value, dark matter,
dark energy, antimatter, or eternal dynamics is derived by this result.
The incoming-gate interval, threshold, source and history below belong to
that application. This finite theorem does not use them. This is a
working-article update. It does not assign the next release version.
The former reservation of 0.27.0 for a closed gate result is dated history.

Read [the paper guide](papers.md). The opening of this page owns the present
programme. This section is the finite theorem.
[The incoming-gate record](#source-fixed-local-incoming-gate-open-campaign-paused)
through the Stage-2 notes is the retained paused application.
[Earlier public checkpoints](#earlier-public-checkpoints) are history and
keep their original numbers.

## Source-fixed local incoming gate (OPEN; campaign paused)

Retained application record. The certification campaign is paused. This
section is not the current regeneration programme, and it is not a
prerequisite for the spherical loop or for the finite theorem above.
The state law, interval, threshold and recorded seed below are unchanged
and apply only here. Missing error bounds still block a closed claim for
this application. Do not continue the expensive search from this page.

The laboratory owners are in `lab/`, indexed in the
[code map](active-code-map.md). The
[draft evidence manifest](../paper/local-gate-evidence/snapshot.json) is a frozen
selection supporting the existing article, not the full current lab.
Its state law and test domain are

$$
C_\Sigma[g]=U_g C_{\mathrm{up}}U_g^\dagger,
\qquad I=S(1)+[0.12,0.18],
\qquad \delta r=\chi(s)[s\,w(z)+s^3U(z)/6].
$$

**Status: OPEN.** The evolved state and its retarded derivative are implemented
in the laboratory. The earlier B1d implementation gaps described below are
historical checkpoints. Unfinished bounds for this application — not a live
campaign — are the rigorous whole-cone field, physical-source, ultraviolet,
continuous, and arithmetic bounds required before a nonlinear search on
this application could be certified.

| Quantity | Recorded meaning |
|---|---|
| Candidate nodal maxima, $N$ and $\beta$ | Approximately $9.64\times10^{-4}$ and $3.81\times10^{-4}$; solver seed, not a root |
| Local existence target | Each continuous residual plus its complete error bound at most $3\times10^{-11}$ |
| v4 error budget | Three finite reused components, six missing; partial sums about $9.00\times10^{-12}$ and $7.00\times10^{-14}$ |
| Field pilot | One spatial period in one trajectory cell; no whole-cone or physical-source certificate |
| Successor 257-node budget | All nine components require authenticated applicable bounds; 129-node phase coverage does not transfer automatically |
| Final outcome | Neither local EXISTENCE nor scoped NON-EXISTENCE is established |

Read the [focused PDF](../paper/local-incoming-gate-draft.pdf) and
[verification guide](local-gate-draft.md). The draft snapshot authenticates
selected original bytes. It is not a full scientific replay or a release of
all active laboratory code. The historical manifests and notebook below stay
unchanged. Frozen-state non-existence results apply to their own declared
classes and do not decide this switched-state gate.

## Shared UV coefficient transport, 29 September 2026

The [shared channel basis](https://doi.org/10.5281/zenodo.23090015) now implements
the existing angular/mass factorization. It reuses one geometry-jet march per
characteristic sign and assembles channel coefficients algebraically. Eleven
tests and independent review passed; the group14 center enclosures retain the
stored widths. This is not full-interval coverage or a higher UV remainder,
and it has not been inserted into the physical residual.

## Low-energy transport method, 29 September 2026

The independently reviewed [signed phase-transport successor](nsc-ks-source-operator-majorant-v7.md)
allows bounded growth in exterior cells where the old contraction-only proof
stopped. Original group14/low16_1 row15 completed through the existing horizon
frame and upstream Bloch transport, with signed source covariance errors about
1.05e-10. Fifteen tests and an independent pilot replay passed. This remains a
one-row result. Its [saved witness](https://doi.org/10.5281/zenodo.23090015) now
replays without either solver, with 13 structural and mutation tests passed.
It is included in the 868-row v3 source aggregate; the full upstream budget
remains OPEN.

## Expanded source coverage and initial vacuum columns, 29 September 2026

The complete middle window now contains 48 original positive rows and their
48 actual negative partners, saved and replayed without re-evolution. The
[v3 response bound](https://doi.org/10.5281/zenodo.23090015) replaces
the earlier single middle pair, rather than adding it twice. Together with
reused high-energy evidence and both verified low-row pairs it covers 868 signed
rows of family14_1. The continuous source-error contributions are bounded by
1.838e-13 for N and 1.461e-13 for beta. Another 300 signed rows in this family, other families,
quadrature and the other components remain missing. These are error bounds;
the candidate residuals above have not changed.

The independently reviewed [vacuum-column remainder](nsc-ks-source-operator-majorant-v7.md)
records eight original-channel cases and supplies phase-fixed initial
coefficients through A4 with an explicit omitted-term bound. It does not yet
supply the propagated changed-history remainder or finite-occupation errors.
Neither result closes the full upstream budget or the incoming gate.

## Continuous preparation-error contribution, 29 September 2026

The [source-response bound](nsc-ks-source-operator-majorant-v7.md)
now connects 772 covered preparation rows to the actual local test: their
combined error contributions are at most 4.064e-14 for N and 3.230e-14 for
beta, continuously on I for the current history. These are partial error
contributions, not new residual values. Remaining source rows/families,
quadrature and the other error components still leave the gate OPEN.

## Lower-middle source preparation, 29 September 2026

The [new correction proof](nsc-ks-source-operator-majorant-v7.md)
resolves the first retained group-14 middle-energy mode and its signed partner.
The archived preparation's covariance error is enclosed near 1.003e-11,
while the correction trajectory's Bloch error is below 4.076e-15. These
source-input quantities still need contraction through the actual N/beta
response. The original preparation remains intact and the gate remains OPEN.

## High-energy source preparation, 29 September 2026

The [new preparation remainder](nsc-ks-source-operator-majorant-v7.md)
bounds a covariance expansion using its equation defect and the exact metric.
Its comparison with the original archived data covers 768 signed rows in one
group-14 source family at |E|>=32. Their largest unweighted covariance error
is below 9.87e-13. The saved source arrays and occupations are unchanged.
Lower energies, other families, quadrature and contraction into the physical
constraints remain unfinished. This does not close the changed-history UV
remainder or any complete local-gate budget component.

## Archived upstream source bound, 29 September 2026

The [upstream successor](nsc-ks-source-operator-majorant-v7.md) bounds
the original saved source on the production upstream slice, including its
opposite-energy/opposite-angular partner. Both covariance errors are below
6.03e-12 for this one low-energy row. No saved preparation is replaced. The
proof replays a small stored curve and rejects changed payloads or a shifted
initial time. Full source coverage and insertion into the local constraint
budget remain unfinished.

## Massive UV coefficient extension, 29 September 2026

The [massive-field extension](nsc-ks-massive-cubic-uv.md) adds the
mass terms derived from the same operator and action, with a directed pilot
for the original group-14 pair at one incoming point. It also factors the
formal coefficient into three shared geometric terms, offering a way to
reduce later repeated calculations. The full high-energy remainder and
continuous interval coverage remain OPEN. No physical residual or budget
entry has been replaced by the formal coefficient.

## Original source covariance bound, 29 September 2026

The [source-preparation enclosure](nsc-ks-source-operator-majorant-v7.md)
now reaches the original rho=1 data for one low-energy massive mode, with
covariance operator error below 6.03e-12. It combines the exterior reflection,
matched horizon conventions and continuous interior transport without changing
the stored source. Remaining energies, angular partners, continuation to the
upstream slice and constraint contraction are still required. This source
sub-bound does not fill a complete error-budget component or change either PDF.

## Directed UV coefficient pilot, 29 September 2026

The [new enclosure](https://doi.org/10.5281/zenodo.23090015) gives rigorous
intervals for one original massless pair's formal cubic coefficient at one
incoming point. It uses analytic geometry jets and explicit integration
remainders. The leading term's integral and its uncertainty are recorded
separately. Continuous coverage of I, the remaining families, higher UV terms
and the physical error-budget component remain OPEN. The current candidate's
residuals and both papers are unchanged.

## Cubic UV method continuation, 29 September 2026

The [massless cubic-current calculation](https://doi.org/10.5281/zenodo.23090015)
now evaluates the generated third-order inverse-energy coefficient for the
original group-1 pair, retaining the linked transport terms and actual action
signs. Its finite quadrature remains unvalidated. C4, C_M and the physical UV
budget are still OPEN; sampled coefficients and refinement differences are
not error bounds.

## Horizon source-method continuation, 29 September 2026

The [full-profile horizon series](nsc-ks-source-operator-majorant-v7.md) now
has an analytic remainder bound, rather than a finite-collar comparison alone.
It connects the existing representative exterior phase enclosure to a horizon
reflection. The series tail at that endpoint is below 2.50e-22; the reflection
uncertainty is still dominated by the earlier phase bound. Production-frame
matching, interior transport and the saved physical source-column error remain
OPEN. No completed component has been added to the nine-component budget.

## Later Stage-2 continuation, 28 September 2026

A new [formal UV identity](nsc-ks-massive-cubic-uv.md) proves
that common scalar upstream constants cancel from the incoming C4 difference
when the higher coefficients are transformed consistently and the incoming
metrics match. It does not assume a zero upstream normalization error or give
a numerical C4. A [source phase continuation](nsc-ks-source-operator-majorant-v7.md)
now reaches a declared near-horizon endpoint, with a bound of about 9.04e-12
radians on one original energy label. Horizon sewing, the actual incoming
source-column error, C_M and the complete gate remain OPEN.

## Stage-2 method advances, 28 September 2026

The active lab now has a directed exterior subgap phase bound, new massless
minor-A3 jet bounds, and an enclosed field cell in the active cutoff transition.
See the [source method](nsc-ks-source-operator-majorant-v7.md),
[UV derivation](nsc-ks-massive-cubic-uv.md) and
[field pilot](https://doi.org/10.5281/zenodo.23090015).
These are sub-bound and method results. Radial/source transport, the remaining
UV constants, full field coverage and the complete local certificate stay OPEN.
The measured seed residuals and the nine-component physical budget are unchanged.

## Earlier public checkpoints

History. The sections below keep their original scope and recorded numbers.
They are background for the retained gate application above. They are not
the current regeneration programme and not a work queue.


**Actual geometric seam matching: PASS.** The [new boundary calculation](nsc-smooth-geometry.md)
uses the two sides of the same smooth $\rho=0$ surface. With opposite outward
normals, the Weyl and Einstein–GHY common boundary covectors cancel, as required
by that domain. The full declared quantum/field variation and bulk self-sourcing remain open;
no new interface force or zero quantum source is assigned. PDF v0.26.0 is unchanged.

On the actual reference seam, the parent Weyl coefficients are
$(\overline P_a,\overline P_r,\Pi_\sigma)=(0.41509246,-0.63006293,0.10900237)$.
The child has opposite metric covectors and the same shear-conjugate
coefficient, because its shear variation also reverses. All geometric
matching residuals vanish; the density/normal checks remain below
$6\times10^{-16}$ against tolerance $3\times10^{-11}$.

The full unit normal is $n=(\beta/a)\partial_\tau-a\partial_\rho$.
Only for stationary background scalars does its action reduce to
$-a\,d/d\rho$. This shorthand is not applied to time-dependent quantum fields.
The old $93.5426$ time-node diagnostic is preserved and is not identified
with either side of this surface.

**Action scope clarified by the author (2026-09-16):**
[$\Gamma_{\rm rest}$ counts only remaining declared contributions](nsc-declared-action-scope.md).
The geometric seam terms and full Gaussian CTP determinant, including its
Schur complement, are already owned. No extra independent boundary functional
is required or introduced. Stationarity must be evaluated from the declared
Gaussian, reference and local action with complete physical state/history jets.

**Historical B1d checkpoint (16 September 2026):** the physical
history-to-KS-endpoint chain rule was not yet implemented at that checkpoint.
Its two diagnostic time-node coordinates must not be identified
with the two sides of the single smooth throat. The existing chart has
$T=T(\rho)$ and $z=\tau+S(\rho)$; replacing KS time by PG time would give
a different derivative. B2, $V_c$ and extended stationarity remain unevaluated.
This clarification adds no numerical result; PDF v0.26.0 is unchanged.

**Weyl boundary variation clarified:** the [intrinsic boundary-jet calculation](nsc-weyl-endpoint-match.md)
reconstructs the stored $93.5426$ end-node derivative from its bulk, boundary
and numerical-stencil contributions. The physical boundary rate is the normal
shear $H_\parallel-H_\perp$. Its coefficient is now explicit; it is not an
extra force chosen to cancel the old number. Physical gluing and the remaining
same-action contribution remain to be specified. PDF v0.26.0 is unchanged.

For the final radius node, the exact stored derivative is recovered as
bulk $+0.68615308$, canonical boundary $-87.27311859$, and stencil contribution
$-6.95567981$. All old nodal derivatives reconstruct within
$9.13\times10^{-12}$ against tolerance $3\times10^{-11}$; the shared eight
original coefficients remain unchanged.

This changes the remaining gate: the physical interface must bind intrinsic
metric and normal-shear variations on one actual surface. The two diagnostic
time ends are not automatically the two sides of that surface. The Gaussian
Schur factor is already part of $\Gamma_G$ and cannot be counted again as a
new $\Gamma_{\rm rest}$ force. No additional term is set to zero.

**New CTP connection:** the [physical mode-history insertion](nsc-pg-ks-metric-pullback.md)
now acts in the same source-mode basis as the prepared PG covariance. Its
[raw KS metric pullback](nsc-pg-ks-metric-pullback.md) retains the full
spatial coordinate dependence. All 63 signed families pass the direct/weak
variation and CTP contraction checks. Complete history-endpoint jets remain
open; PDF v0.26.0 is unchanged.

| New component | Maximum residual | Tolerance |
|---|---:|---:|
| PG weak/direct Hamiltonian variation | $5.63\times10^{-10}$ | $3\times10^{-8}$ |
| Raw KS weak/direct variation | $2.18\times10^{-9}$ | $3\times10^{-8}$ |
| CTP tangent condition | $2.78\times10^{-17}$ | $3\times10^{-11}$ |
| Coordinate metric reconstruction | $1.78\times10^{-15}$ | $3\times10^{-11}$ |

The actual mode-space kernel determines the conditional time-ordered insertion.
Its raw KS conversion is integrated point by point; using only the neck's
Jacobian gives a different result. The fixed-surface normal derivative is
also evaluated. **Remaining B1:** the actual two-endpoint variation support
and Cauchy embedding, $\delta g_P^A(\tau,\rho)/\delta g^B_{K,e}$.
No source stress, physical history or metric timestep is assigned here.
Six focused tests and both all-field verifiers pass.

**Full retained PG preparation: C1b PASS.** The [common-time covariance](nsc-pg-retained-covariance.md)
now includes all **32 massive groups plus the LLL channel**, with both signed
angular families and correlations with the unresolved bulk. Its independently
integrated CAR residual is at most $1.40\times10^{-10}$, against the declared
$3\times10^{-9}$ numerical tolerance. Historical seeds and physical scales are
unchanged. The working PDF is now **v0.26.0**.

At this earlier checkpoint, the next calculation coupled the prepared state to
history jets and endpoint variation. The active source-fixed incoming-state
and constraint path is now implemented; its remaining validation work is
described in the opening section of this page.

| Computed quantity | Result |
|---|---:|
| Signed family evaluations, including LLL | 63 |
| Largest CAR residual | $1.40\times10^{-10}$ |
| Covariance eigenvalue range | $[0.0027381,0.9972619]$ |
| Largest measured low-energy covariance refinement | $2.89\times10^{-11}$ |
| Largest measured middle-band covariance refinement | $1.18\times10^{-14}$ |
| Largest covariance-tail representation bound | $9.73\times10^{-10}$ |

[Full record](../results/development/nsc-pg-retained-covariance.json) ·
[Join/error audit](https://doi.org/10.5281/zenodo.23090015).
The whole-field mode map retains $C_H$, incoming occupation and the probe
complement. No closed eight-mode state replaces the bulk. All-field replay
and three focused tests pass. $\Gamma_{\rm rest}$, the physical two-sided
Weyl mismatch and extended stationarity remain OPEN; no metric timestep is
started. Earlier entries below retain their dated checkpoint scope.

**Retained-state integration inputs: PASS.** Two further input collections
are now complete for the full covariance construction:

- [Compact subgap panels](nsc-pg-retained-covariance.md): all **38 signed families**
  in the 20 compact groups; largest summed quadrature estimate
  $3.14\times10^{-11}$.
- [High-energy tail inputs](nsc-pg-retained-covariance.md): all **61 signed families**
  for the remaining 31 groups. The largest covariance representation bound is
  $9.73\times10^{-10}$ and the Gram bound is $1.95\times10^{-9}$.

Five groups require the numerical transition at 320 instead of 160. This
changes the numerical integration split, with the same physical operator and
parameters. The [precision audit](https://doi.org/10.5281/zenodo.23090015) independently
recomputed all 29 new energy tables with 50-digit process-local precision;
their stored numerical values were identical.

**Full retained C1b remains OPEN at covariance assembly and its error checks.**
The completed group-13 result below remains the first full-energy massive
covariance. PDF v0.25.0 and the metric-evolution gate are unchanged.


**Earlier checkpoint — first full-energy massive PG covariance: PASS.**
The [new calculation](nsc-pg-retained-covariance.md) completes the retained
$m=\pi/2$, zero-angular-eigenvalue group on the common PG Cauchy surface.
The [record](../results/development/nsc-pg-group13-covariance.json) contains
the covariance, independently integrated CAR matrix, and lesser/greater/Keldysh
blocks, with the complete energy account and its numerical error checks.

| Quantity | Result |
|---|---:|
| CAR residual | $3.27\times10^{-11}$ |
| Covariance eigenvalue range | $[0.03570181,0.96429819]$ |
| Original-probe / bulk-probe correlation norm | $0.31816523$ |
| Finite-energy covariance refinement | $2.81\times10^{-15}$ |
| Tail mode-order change | $5.33\times10^{-15}$ |
| Tail spatial-grid change | $1.26\times10^{-12}$ |
| Fourier/normalization covariance-tail bound | $5.70\times10^{-10}$ |

The primary covariance uses the **calculated** CAR Gram; no identity correction
is used to fill a truncated frequency table. The eight probes inspect the full
mode-defined field and its unresolved-bulk correlations. They do not form a
closed eight-mode evolution model.

**Full retained C1b remains OPEN:** 31 massive groups remain. The locked action,
scales, seeds and all 209 earlier scientific files are unchanged. The existing
fourth-order reference/local allocation is preserved; no stress or metric
timestep is assigned here. **PDF v0.25.0 remains unchanged.**

Earlier entries below retain their checkpoint scope.


**Previous checkpoint — massive PG mode/Green construction: PASS.**
The [new derivation](nsc-pg-massive-mode-resolution.md) and
[record](../results/development/nsc-pg-massive-mode-resolution.json) evaluate
all **32 massive groups** on the complete existing node tables: **3,584
positive-energy/angular-sign combinations**, with signed-frequency partners
fixed by the existing operator map.

The input is the affine horizon covariance and inherited incoming occupation.
Complex Jost amplitudes give the exterior modes; the matched inner preparation
gives their transmitted components and independent interior horizon partner.
An open infinity channel remains present even when its occupation vanishes.

| Check | Maximum | Tolerance |
|---|---:|---:|
| Complex scattering isometry | $1.12\times10^{-12}$ | $3\times10^{-11}$ |
| Physical source law | $0$ | $3\times10^{-11}$ |
| Source CAR violation | $2.66\times10^{-17}$ | $3\times10^{-11}$ |
| Whole-line Green jump, relative | $8.29\times10^{-10}$ | $3\times10^{-9}$ |
| PG fields under outer-radius doubling | $7.79\times10^{-10}$ | $3\times10^{-9}$ |

Seven focused tests pass, including rejection of an omitted open infinity
channel. The old seed generations and all 201 previous scientific files are
unchanged. **Full C1b remains OPEN at the packet energy integral and its
contact/tail account.** A finite mode table is not used to fill the continuum
CAR identity. PDF v0.25.0 and the metric-evolution gate remain unchanged.

The earlier entries below retain their individual checkpoint scope.


**Previous checkpoint — full compact restart: H1 PASS.** The
[versioned restart](https://doi.org/10.5281/zenodo.23090015) now regenerates every
retained compact frequency from the original horizon/incoming data: 20 groups,
1,280 X blocks, 1,152 angular-partner Y blocks and their signed-frequency
partners. The historical seed is used only after construction, as a control.

| Start convention | Explicit control | Result |
|---|---|---|
| `historical_radial_r_h` | Immutable original compact seed | Maximum X residual `4.11e-14` |
| `matched_delta_q` | New `compact-seed-v2-matched-delta-q` checkpoint | Cold regeneration reproduces the stored artifact byte for byte |

The new matched control uses the corrected interior Frobenius argument and
shares the same inherited propagator, source and scattering with the historical
comparison. The two conventions are not made numerically equal. Their largest
X difference on the full grid is `2.73e-4`; it remains a diagnostic of the finite
start convention.

The [record](../results/development/nsc-compact-matched-restart.json) verifies
unitarity `5.97e-14`, matched modal normalization `5.99e-14`, T feedback
`4.37e-16`, and the prior independent matched modal control `9.20e-12`
(tolerance `3e-11`). The loader requires an explicit convention and hash.
The original seed, source file and previous result records are unchanged.

**H2 and full C1b remain OPEN.** The complete compact seed-frequency table is
now available under a consistent matched convention. It still needs the massive
global PG mode transform and current/Plancherel completeness before spatial
covariance assembly. No massive spatial C, transmitting endpoint jets, new
stress or metric evolution is inferred from this checkpoint. **PDF v0.25.0
is unchanged.**

**Preceding development — paired horizon construction and a compact normalization correction.**
The [new map](https://doi.org/10.5281/zenodo.23090015) normalizes the magnetic angular
spinors, their gauge transition and their signed pairing. Under the already
selected affine-horizon preparation, orthogonality and angular-sign conservation
fix $Z(E)=0$. The horizon-partner correlations inside each angular sector remain.
The other marginal $Y(E)$ is computed by its own massive mode propagation,
without using stored seed covariance as input.

The [record](../results/development/nsc-horizon-paired-pg-map.json) evaluates
one existing frequency in each of the 32 massive groups. Magnetic angular
residuals are below `1.16e-13`; matched frame and mode residuals are below
`6.44e-12` (tolerance `3e-11`). These are finite modal maps on the trapped
probe chart, not a completed spatial PG covariance or all-frequency state.

**Finite compact collar mismatch identified.** The stored interior uses
$\delta_q=q_{h,\mathrm{geom}}-q_{\mathrm{geom}}$, with
$q_{\mathrm{geom}}=\pi/2+\arctan\rho$; this coordinate is distinct from the
fixed magnetic flux. Since $|\rho-\rho_h|=r_h^2\delta_q+O(\delta_q^2)$,
the interior Frobenius distance is $\sqrt{2\delta_q/\kappa_h}$.
The historical compact initializer additionally divides by $r_h$.
That factor is correct for an exterior radial offset, but not for this interior
coordinate. The new owner uses the matched normalization and leaves the old
files untouched.

All twelve angular-only controls recover their stored X. For all twenty
compact controls, the legacy frame recovers X while the matched frame changes
it by norms between `3.88e-5` and `1.47e-4`. This isolates the discrepancy to
the finite initial frame. Both conventions have the same zero-collar limit;
no new physical interaction or fitted coefficient is introduced.

**C1b-H remains OPEN.** The compact retained-state checkpoint must be reconciled
under the matched collar convention, and the global signed PG spectral/current
resolution must be completed before spatial covariance assembly. No corrected
modal counterpart is silently substituted into the locked seed. C1b-M2 and
C2–C4 remain unexecuted. **PDF v0.25.0 is unchanged.**

**Preceding development — massive signed-spin map and data-sufficiency gate.**
The [new operator map](nsc-ks-signed-state.md) connects the imported
paired angular representation to T's current basis. Its coefficient identities
are exact; the maximum numerical residual is `6.06e-15` (tolerance `3e-11`).
T normal-current and signed-frequency checks have zero residual.

Fourteen groups admit simple radial signed-frequency maps: the twelve
angular-only groups and two compact-only groups. The eighteen groups with
both compact mass and angular potential need the paired representation for
this local signed-frequency symmetry. Their pointwise two-component axis flip
carries a required basis-connection term; dropping it would change the operator.

**C1b-M OPEN, with a verified information certificate.** Across all 1,152 stored
mixed-group seed blocks, two algebraic completions preserve the observed
marginal exactly while differing in the missing angular/signed partner by
Frobenius norm `sqrt(2)`. CAR violations from numerical roundoff stay below
`4.18e-14`. These are information controls, not physical state choices.
The [record](../results/development/nsc-massive-signed-preparation.json) names
the required opposite-angular marginal and cross block, to be generated by
the horizon/infinity-to-paired-PG modes in the current/gauge basis. The
operator involution does not by itself prove physical state symmetry.
No massive PG covariance or transmitting branch jet is assigned. Full C1b
and C2–C4 remain OPEN; **PDF v0.25.0 is unchanged**.

**Preceding development — physical LLL PG Cauchy preparation.** The
[LLL preparation](nsc-pg-lll-preparation.md) uses the original horizon pair and
inherited incoming occupation to construct a covariance on the same common
PG-time slice as the retarded memory operator. It retains the exterior outgoing
partner and packet/bulk correlations; stored seed blocks are comparison data,
not copied spatial inputs.

The [record](../results/development/nsc-pg-lll-preparation.json) gives seed
reconstruction `5.56e-17`, spatial CAR-Gram residual `7.57e-15`, and quadrature
refinement `4.34e-12` (tolerance `3e-11`). The independent characteristic
response matches the locked LLL resolvent to `1.53e-14`. Its seven-probe
covariance has eigenvalues between `0.11534` and `0.88466`; the packet-to-sampled-
bulk correlation norm is `0.50269`. These are quantum correlation matrices,
not stress or energy fractions.

**LLL C1b PASS; full C1b and C2–C4 OPEN.** The 32 massive groups still require
the signed-frequency spin/angular mode map in T's current basis, or the full
paired four-spinor PG description and its state correlations. Degeneracy and
copy count do not supply that map. Full transmitting `EndpointBranchJets`,
remaining boundary variation and physical Weyl mismatch stay unevaluated.
This development record follows the unchanged **PDF v0.25.0**; no new paper
version or metric evolution is claimed.

**Preceding result — causal memory and temporal metric transfer.** The
[new CTP connection](nsc-local-observer-correspondence.md) places the computed
packet response in the convention $G_J^R=-\mathcal R$, $D_J^R=-\mathcal K$.
It preserves the energy-dependent memory and calculates the four metric
vertices between distinct frequencies, using the same transmitted Dirac
operator and archived incoming fields.

All 33 channels pass: driven-versus-weak-form residual `5.57e-13`, adjoint
frequency-exchange residual `6.19e-13`, and inverse-memory Dyson residual
`5.14e-16` (tolerance `3e-11`). The zero-transfer limit reconstructs the
locked static derivatives; reverse retarded response and its metric transfer
remain zero in the trapped probe cell. The harmonic source is a linear-response
probe, not a selected metric history or physical frequency prediction.

**Retarded-memory C1 component PASS; full C1 and C2–C4 OPEN.** The
[record](../results/development/nsc-transmitting-ctp-resolvent.json) retains the
exact missing preparation: the physical covariance on a common PG-time slice,
including its packet/complement cross correlations. That state, full history
reconstruction and the KS endpoint pullback are required before the memory
vertices become `EndpointBranchJets`. The existing consumers leave the
remaining boundary derivative and physical Weyl mismatch unevaluated.

**Preceding result — domain-correct transmitting response: Z2a PASS.** The
[cross-resolvent calculation](nsc-transmitting-dirac-domain.md) evaluates
all 33 retained channels on the inherited whole-line Dirac domain. Continuous
seam data cancel the individual sharp-projector delta terms. Independent
parent/child spatial probes then give an energy-dependent inverse response
$\mathcal K(z)=[J^\dagger(H_D-z)^{-1}J]^{-1}$.

In the trapped probe cell, the retarded parent-to-child response is nonzero;
the reverse retarded response is zero by causal direction. This is a response
statement, not a choice $B=0$. The adjoint residual is `3.53e-14`, the true
resolvent-identity residual `1.18e-13`, and seam mismatch zero (tolerance
`3e-11`). The four static metric kernels also pass: independent weak-form
agreement is below `6.56e-13`; the finite-difference check is `1.01e-8`
(tolerance `3e-8`). The [record and payload](../results/development/nsc-transmitting-cross-resolvent.json)
store the response, its derivatives and uncompressed quadrature fields.

**Static Z2b PASS; time-history CTP bridge and Z4 OPEN.** The computed
resolvent derivatives need time/spectral reconstruction, physical state/memory
data and the endpoint-coordinate pullback before they become
`EndpointBranchJets`. The remaining boundary action and physical two-sided
Weyl mismatch are unevaluated. No stationary history or stress is assigned.
Moving-neck jets remain out of scope for the fixed coordinate seam.

**Preceding definition — common-time bulk split: Z1 PASS.** The
[new bulk domain](https://doi.org/10.5281/zenodo.23090015) uses the inherited PG time
and spatial supports $\rho>0$ and $\rho<0$. Its finite witness has parent and
child ranks 1,904 each, with zero projector/CAR residual and rank defect zero.
The current matching to T has residual `1.78e-15` (tolerance `3e-11`).
The [duplicated-trace FAIL](https://doi.org/10.5281/zenodo.23090015) is preserved.
The ordinary sharp-block Z2 route remains OPEN because projecting nonzero
transmitting trace data leaves the Dirac operator domain. No delta coefficient is promoted
to a Hamiltonian link or source.

**Preceding definition — transmitting Dirac seed-seam domain: PASS in its finite scope.**
The [new domain](nsc-transmitting-dirac-domain.md) identifies the 33 retained
channels with two oriented traces of one field on the stored spacelike
surface. All 1,904 blocks are mapped using the actual coframe, quadrature
and spin-current basis. Normalization residual is `6.66e-16`, covariance
recovery `5.55e-16`, and the oriented boundary-form residual `1.08e-16`
(tolerance `3e-11`). The same-surface domain is now explicit; an instantaneous
Hamiltonian link, selected Cauchy history and remaining boundary action are
still separate definitions. No physical propagator or stress is assigned.

**Preceding calculation — known bulk endpoint jets.** The
[Hamiltonian-to-CTP chain](https://doi.org/10.5281/zenodo.23090015) is now executable for
the already owned bulk vertices, using a fixed representative block in each
of the 33 channels on the frozen control. The derivative residual is
`1.0276e-10`; its contraction against the existing CTP action differs by
`4.9394e-11` (finite-difference tolerance `3e-8`). These are diagnostic bulk
jets, not complete physical transmitting jets. The
[J→M gate](https://doi.org/10.5281/zenodo.23090015) remains **OPEN** for the link,
embedding and remaining boundary derivatives.

**Preceding calculation — transmitting CTP first differential.** The
[new action differential](nsc-local-observer-correspondence.md) agrees with the
existing determinant in all 33 retained channels, with maximum error
`5.551121183780314e-12` (finite-difference tolerance `3e-8`). Common-branch
variation is exactly zero; relative variation supplies the conditional
Gaussian source differential. Its [endpoint chain rule](nsc-pg-ks-metric-pullback.md)
preserves the old covector and leaves the missing physical branch jets and
boundary remainder explicit. The [F→H decision](https://doi.org/10.5281/zenodo.23090015)
is **OPEN**; the old certificates are unchanged.

**Preceding binding result — full extended gate: OPEN.** The
[transmitting-action audit](https://doi.org/10.5281/zenodo.23090015) finds
norm-preserving directions that change the stored covariance in all 33 channels,
with maximum tangent norm `0.679080121964124` and isometry residual zero
(tolerance `3e-11`). The [Weyl binding](nsc-weyl-endpoint-match.md) extracts
all eight diagnostic end-node coefficients exactly. The physical boundary
kernel selector, endpoint derivative and its pullback remain unevaluated.
The value `93.54264532195464` is the local diagnostic gradient, not an
evaluated two-sided physical mismatch.

The project starts from one simple claim: a wave, a particle, an unseen
boundary response, and a changing geometry can be different measurements of
one inherited spectrum. This page shows the equations and numbers behind that
claim.

| Plain-language idea | Calculated connection |
|---|---|
| A particle is a stable note of a field | One scalar link gives the exact massive Dirac dispersion and propagator self-energy |
| Another room can be unseen but still act here | Direct propagation and Schur elimination give the same curved boundary response |
| A boundary stores and transfers energy | The real and imaginary parts of one coherence give link energy and occupation transfer |
| Geometry can turn into particles | Work from a changing radius equals the energy of the created Dirac pairs |
| An inherited state can push geometry apart | Both computed radial null contractions are negative at the smooth neck |
| The completed stress has a child-time energy ledger | Proper volume and the child clock fix its first conserved density derivative |
| Backreaction begins with a constraint gate | The locked tensor and neck fail the Hamiltonian and homogeneous momentum constraints at $T=0$ |
| The tensor can select valid replacement geometry | Its Landau frame and density determine a constraint-complete Kantowski--Sachs neck |
| Future pressure needs the full Gaussian state | Four integrated stress moments do not determine the mode commutators under the changed metric |
| The seed Gaussian state is now inspectable | All 1,904 retained covariance blocks are serialized and reconstruct the completed tensor |
| Cauchy propagation is history dependent | Conditional maps are unitary, while equal endpoint data produce different propagated states |
| The present joint BVP is non-unique | The seed cannot start and the endpoint problem lacks an executable same-action history selector |
| The homogeneous history class is decided | The complete shift ledger excludes every smooth frequency-diagonal no-interface history from the locked seed |
| The larger tilted class now has real component owners | Reference, local bulk and weighted interface gates run separately before stationarity composition |
| Empty uniform space defines the zero | The recursive $a_0$ projector fixes $V_{\rm full}=0$ and preserves all gradient terms |
| Charged radius and spectral coefficients meet at one scale | The first cutoff-resolved checked branch gives $\Omega=3.9730743688$ and $\zeta=15.7853199397$ |
| The direct MMP-LLL substitution is decided | Its recursive power is fixed, but both neck null components have the wrong sign |

The construction connects the locally resolved field, the response across its
boundary, and the stress acting on its geometry through the same operator.
General research status is stated in the [README](../README.md#research-status).

The [September 10 development snapshot](https://doi.org/10.5281/zenodo.23090015),
through laboratory commit `3ba9e01`, now contains forty-six post-preprint records.
Its latest gates bind the charged field, parity, AP state and normalization to
the imported MMP throat; implement canonical CTP as the causal owner; and
apply the recursive zero-tadpole law that fixes the homogeneous unlinked
$V_{\rm full}=0$. The charged radius equation now supplies a fixed-$q$ scale
candidate. Applying its simplest inherited transparent LLL state to the actual
black-universe neck gives positive rather than required negative null sources,
so the missing owner is the full charged angular/compact CTP tensor on that
domain rather than another scale scan.

The subsequent [charged CTP neck calculation](https://doi.org/10.5281/zenodo.23090015)
returns the retained child-frame tensor
$(\rho,T_{01},p_\parallel,p_\perp)=(0.0468138,0.00122387,-0.106079,0.0302511)$
and negative radial null components $-0.0568175$ and $-0.0617129$ at the
locked scale. This passes the retained sign test, while the hard full gate
remains open for the first positive compact levels' nonlocal CTP covariance
and four ADM variations; no cosmological projection is accepted yet.

That channel is now supplied by the
[positive compact CTP completion](https://doi.org/10.5281/zenodo.23090015). It reads the
previous angular tensor without rerunning it and adds the two positive compact
levels below the locked cutoff with massive exterior scattering and
fourth-order superadiabatic subtraction. The completed free-Gaussian tensor is

$$
(\rho,T_{01},p_\parallel,p_\perp)
=(0.10074836289,0.00122387016,-0.34978809527,0.12958962232),
$$

with $T_{++}=-0.24659199207$ and $T_{--}=-0.25148747269$. The hard neck sign
gate therefore passes in this realization. Its homogeneous neck projection is
recorded.

The [proper-volume/clock projection](https://doi.org/10.5281/zenodo.23090015) now places
that tensor in the stored child geometry. At the neck,
$V=24.21233094747$, $\dot V=12\pi$, and conservation gives

$$
\boxed{\dot\rho_0=-\frac{\dot V_0}{V_0}\rho_0
-H_\parallel p_\parallel=0.38776013531.}
$$

The volume term is $-0.15686733379$ and pressure work is
$+0.54462746910$. The recorded parent Killing power satisfies
$P_{\rm parent}=-a_\parallel VT_{01}$ with residual
$4.4\times10^{-16}$: inside the child it is the opposite of a conserved
spatial-momentum charge, not an additional energy-deposition rate. The free
compact blocks give regular bulk $Q=0$. This fixes a local Cauchy jet; the
pressure derivatives and global background require backreaction and evolution
of the same CTP state before any $H(z)$ calculation.

The [child backreaction gate](https://doi.org/10.5281/zenodo.23090015) applies the
existing ADM constraints before starting that evolution. With the locked
$A=0.04501936182826115$, the stored neck requires
$\rho=2A=0.0900387236565223$ and $T_{01}=0$. The completed tensor instead gives

$$
\boxed{\mathcal C_H=-0.118944813912720,
\qquad \mathcal C_M=-0.0135927088490713.}
$$

Both null components are still negative. The evolution stops at $T=0$
because stepping constraint-violating data would not solve the coupled
Einstein--CTP system. The required same-action completion is
$\Delta\rho=-0.0107096392302639$ and
$\Delta T_{01}=-0.00122387015580509$, or a source-derived initial geometry.
Neither repair is inserted in this result.

The [constraint-complete assembly](https://doi.org/10.5281/zenodo.23090015) then takes
the permitted geometry-from-source route. The tensor's unique subluminal
Landau boost,

$$
v=0.00491447570673408,
$$

sets the normal-frame momentum to zero without changing occupations. At a
temporal minimum, the Hamiltonian constraint fixes

$$
\boxed{r_\star=\sqrt{\frac{2A}{\rho_L}}
=0.945328394434129.}
$$

The maximum constraint and frame-invariance residual is
$2.2\times10^{-16}$, with
$T_{++,L}=T_{--,L}=-0.249027703022043$. This PASS geometry remains in the
existing Kantowski--Sachs class but replaces the stored unit-radius Bronnikov
profile. Joint evolution of the same covariance and this selected metric is
therefore the next executable owner.

The [coupled CTP/metric evolution gate](https://doi.org/10.5281/zenodo.23090015)
checks the first required stress re-evaluation. Because

$$
\frac{1}{r_\star}-1=0.0578334533139646,
\qquad
\left.\frac{\partial\rho}{\partial\log r}\right|_C
=-0.460687999769215,
$$

the unit-radius tensor cannot be reused unchanged on the selected geometry.
The published source records retain integrated stresses and covariance bounds,
but not the complex covariance for each compact/angular/frequency mode.
Consequently the exact break occurs at $T=0$, during same-covariance source
evaluation before the first metric step. The fixed-tensor constraints and null
signs still pass; future pressure becomes determined only after the recorded
`ModeResolvedCauchyState` contract is supplied.

The [mode-resolved state artifact](nsc-mode-resolved-cauchy-state.md) now
supplies the seed side of that contract. A deterministic 407,306-byte NPZ
contains 33 channels and 1,904 physical covariance blocks. Reload gives

$$
\lambda_{\min}(C)=-2.67\times10^{-16},
\qquad
\lambda_{\max}(C)=1+4.15\times10^{-14},
$$

and reconstructs all completed unit-radius tensor components with maximum
residual $9.11\times10^{-13}$. The physical Landau Cauchy isometry and the
general-KS fourth-order/local source history are not determined by endpoint
$(v,\eta,r_\star)$, so finite selected-surface stress and metric evolution
remain gated on those two explicit owners.

The [Landau Cauchy-isometry owner](nsc-landau-cauchy-isometry.md) now computes
the blockwise time-ordered Dirac exponential for a supplied ADM history. Two
smooth histories with identical seed and Landau endpoint data give

$$
\max\lvert U_{0.5}-U_{1.0}\rvert=1.99966699255,
$$

while each satisfies $\lVert U^\dagger U-I\rVert<1.2\times10^{-14}$ and
preserves the CAR spectrum. This makes the selection issue finite and exact:
the endpoint tuple does not determine the physical $U_{L0}$. Neither control
is promoted to a metric history, so a unique finite $r_\star$ stress and its
updated constraints remain undefined.

The [joint history/state BVP gate](https://doi.org/10.5281/zenodo.23090015) now decides
the two allowed solve routes. Dynamic route A stops on

$$
(\mathcal C_H,\mathcal C_M)
=(-0.118944813913,-0.0135927088491)
$$

at the seed. Two-point route B has at least two endpoint-compatible solutions,
with

$$
\max|C_A-C_B|=0.192089815976105.
$$

The chosen result is route C: non-uniqueness bound and stop. The missing
selector is the executable general-KS history variation of the existing CTP
action, including its fourth-order, local and boundary terms. No finite stress
or metric history is assigned until that owner exists.

The [general-KS same-action owner](https://doi.org/10.5281/zenodo.23090015) now
executes the independent shift variation in the minimal homogeneous class.
Every already owned local, fourth-order-reference and magnetic contribution
has zero homogeneous $T_{01}$, whereas the serialized state yields

$$
\boxed{\mathcal E_\beta(0)
=-\frac{T_{01}}{2A}
=-0.0135927088490713.}
$$

Smooth unitary evolution makes this nonzero source continuous, while the
homogeneous momentum constraint requires zero on every slice. Hence no smooth
stationary history exists in the frequency-diagonal no-interface KS
truncation. A tilted transmitting interface or spatially inhomogeneous
spherical history lies outside this certificate.

The extended class now has separate executable
[fourth-order reference](https://doi.org/10.5281/zenodo.23090015),
[node-wise local history](https://doi.org/10.5281/zenodo.23090015), and
[tilted interface](nsc-mode-resolved-cauchy-state.md) owners. Their composition
gives

$$
R_{\rm ref}\le 7.06\times10^{-12},\qquad
R_{\rm local,node}=2.08\times10^{-7},\qquad
R_{\rm interface}\le2.35\times10^{-14}.
$$

The [extended gate](https://doi.org/10.5281/zenodo.23090015) nevertheless remains
open: the diagnostic Weyl end-node gradient is $93.54264532195464$, and
the weighted-isometry equation leaves $447488$ real kernel directions.
The new [endpoint binding](nsc-weyl-endpoint-match.md) preserves that gradient
and leaves the physical mismatch unevaluated, with target tolerance
$3\times10^{-11}$. The same transmitting boundary action
must select one $V_c$ and complete that endpoint variation before
$\delta\Gamma$ is single-valued. No optimizer, finite stress or existence
claim is launched before then.

## 1. Shared mass and visible response

In the invariant scalar-link sector,

$$
[H_8,\Pi_-]=0,\qquad
E^2=\lvert\mathbf p\rvert^2+\Phi^2,\qquad
\Sigma_p(E)=\Phi^2(E+\boldsymbol\alpha\cdot\mathbf p)^{-1}.
$$

The same $\Phi$ determines the mass shell and visible self-energy.
[The Dirac reduction](nsc-observable-bridge.md) specifies the projector, spin
frame and charge conventions. The complementary sector is spectrally
identical. The evaluated massless throat kernel has vector and axial
Clifford content; the scalar coupling and sector selection are action/domain
requirements.

## 2. Evaluated curved boundary maps

Parent and child maps are obtained from the smooth Dirac geometry with
oriented normal data. Direct inversion, Schur elimination and transfer agree
below $2\times10^{-15}$ relative error in the recorded finite system.
Independent radial integration differs by 0.00947%–0.03645% at the stated
resolution. The [full spinor calculation](https://doi.org/10.5281/zenodo.23090015) retains both
angular sectors and the common spin frame.

At common dimensional energy, the recursive response is

$$
\Gamma_p(x)=K_p(x)-\Omega^{-1}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.
$$

The explicit $\Omega^{-1}$ is the first-order normalization factor.
[Finite recursion](nsc-regulated-recursion.md) and
[endpoint criteria](nsc-tail-limit.md) supply the recorded controls.

## 3. The state completes the energy account

For the finite quadratic system,

$$
z_B=\mathrm{Tr}_p(B C_{cp}),\qquad
E_{\rm link}=2\Re z_B,\qquad \dot N_p=2\Im z_B,
$$

$$
\dot E_{\rm tot}=\mathrm{Tr}(C\dot H).
$$

[The common-source derivation](nsc-common-source-derivation.md) includes
energy stored in the link. [The boundary-state implementation](https://doi.org/10.5281/zenodo.23090015)
retains occupied and empty kernels and mixed initial correlations.
Its finite covariance reconstruction agrees with full evolution to
$1.1\times10^{-15}$.

A prescribed changing radius creates pairs with excitation energy equal to
the supplied work. In the representative $\kappa=1$ channel,
$\epsilon=0.02$, $\tau=0.5$ gives approximately
$1.40181\times10^{-5}$ pairs and energy $3.46640\times10^{-5}$
in throat units. [The vacuum-work record](nsc-vacuum-work.md) and
[normalized influence functional](https://doi.org/10.5281/zenodo.23090015) connect this response
to its geometric fluctuation weight.

## 4. Parent data produce a definite child stress

Affine-horizon data and the actual complex exterior reflection determine
the canonical massless state on the expanding child. In the child
orthonormal frame, with $\hbar=c=L_{\rm throat}=1$, $\mu=1$:

| Quantity | Recorded neck value |
|---|---:|
| Density $\rho$ | $-0.0033068962$ |
| Longitudinal pressure $p_\parallel$ | $-0.0494945966$ |
| Angular pressure $p_\perp$ | $-0.0429957713$ |
| Radial null $+$ | $-0.0528075900$ |
| Radial null $-$ | $-0.0527953956$ |
| Outward parent Killing power | $0.0001422206795$ |

Stress has dimensions $L^{-4}$, and power $L^{-2}$. Angular refinement
from 16 to 32 after the recorded tail correction changes tensor components
by less than $4\times10^{-6}$ across the sampled radii. This measures
numerical sensitivity under that asymptotic-tail approximation.
[Source and state](nsc-unruh-state.md).

The positive-Einstein static neck requires

$$
(\rho,p_\parallel,p_\perp,T_{\hat T\hat z})
=(2a,2a(1-3\pi),2a(1-3\pi),0),\quad a>0.
$$

The computed null signs are negative; density, pressure anisotropy and flux
leave explicit componentwise residuals. The [ADM frame map](https://doi.org/10.5281/zenodo.23090015)
places those quantities in the metric equations. Under the stated asymptotic
Einstein normalization, the outgoing channel contributes
$dM_B/du=-P_0$, selecting a backreaction calculation.
[Flux and mass balance](nsc-parent-backreaction-gate.md).

## 5. One normalization and one metric variation

The [canonical–spectral bridge](https://doi.org/10.5281/zenodo.23090015) partitions
the compact determinant into its physical tower, cutoff conversion and warp
contribution. The [Gaussian measure](https://doi.org/10.5281/zenodo.23090015) gives an
exact finite representation of the adopted proper-time action.
[Compact anomaly matching](https://doi.org/10.5281/zenodo.23090015) fixes the domain and
the scope of the imported anomaly result.

The [ADM source equations](https://doi.org/10.5281/zenodo.23090015) retain lapse, shift,
radial metric and sphere radius through variation. Their Ward identities
bind stress, energy, flux and geometric work. The
[local warp source](nsc-parent-backreaction-gate.md) evaluates an identified
summand on the actual neck; its nonlocal remainder has no small-error bound.

## 6. The causal and vacuum source now have explicit owners

Canonical Lorentzian evolution supplies the causal state-dependent response:

$$
\Gamma_{\rm one}^{\rm CTP}
=S_{\rm induced}[+]-S_{\rm induced}[-]
-i\log\det(I-C_0+C_0U_-^\dagger U_+).
$$

The [causal common functional](https://doi.org/10.5281/zenodo.23090015) returns the
four metric forces, stress projections, boundary power, retarded response,
noise, and discrete Ward residuals from one covariance. Equal histories,
unitarity, covariance trace, and energy/work balance agree at approximately
$10^{-16}$ in its finite control.

The homogeneous gravitational zero is fixed separately by the same-action
relational law

$$
\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one},
\qquad
(V,A,C)\mapsto(0,A,C).
$$

This gives $V_{\rm full}=\lambda_4=\Xi=0$ while preserving Einstein, gauge,
curvature, Casimir, and link response. The
[normalization derivation](nsc-nested-qualities.md) proves
idempotence, additivity, CTP normalization, and common room-scale covariance.

The physical geometry is the common stationary point:

$$
\left.\frac{\delta\Gamma_{\rm one}^{\rm CTP}}
{\delta g_\Delta^{\mu\nu}}\right|_{g_\Delta=0}=0.
$$

The same solution determines the following physical links:

| Target | Quantity to extract from the common solution |
|---|---|
| Constants and resolution | Characteristic cones, gravitational coupling, physical thermal state |
| Finite continuation | Evolved geometry, constraints, stress and boundary flux |
| Dark response | Background pressure and density, anisotropic stress, growth and lensing |
| Particles and measurement | Poles, residues, interactions, detector statistics |
| Mass and antimatter | Chiral interaction, sector selection, gauge and conjugate representations |
| Recursive infinity | State marginals, clock transfer and limiting observables |

## Inspect the record

The [preprint manifest](../results/manifest.json) contains 100 records;
[the later snapshot](../results/development-snapshot.json) contains forty-six.
[The development index](https://doi.org/10.5281/zenodo.23090015) links their
derivation notes. The canonical paper retains its technical appendices,
attribution and historical scale corrections.

The default demonstration displays stored evidence. Use
[the reproduction guide](reproducing.md) for an explicitly requested replay.
