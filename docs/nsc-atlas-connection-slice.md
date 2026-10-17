# One expanding interior and an asymptotically Euclidean parent end

This successor connects the [fixed v3 moving wall](nsc-atlas-moving-boundaries.md)
to a regular de Sitter centre and a Schwarzschild bridge on one spacelike
initial slice. It supplies one body with an exact vacuum outer end. It does
not make the two-hole Brill-Lindquist data spherical or establish formation.

## Admission before execution

- **Borrowed:** the signed de Sitter/Schwarzschild wall already imported in
  v3; Kuchar [gr-qc/9403003v1](https://arxiv.org/html/gr-qc/9403003v1), Sections
  2.2, 3.1--3.2 and 4.1, equations (14)--(18), (23), (27)--(33), (81)--(82),
  for spherical slices, their constraints and invariant Schwarzschild mass;
  Hajicek [gr-qc/9708008v2](https://arxiv.org/html/gr-qc/9708008v2), Sections
  2 and 4.1, for shell stress projections and Lorentz-related frames.
- **Missing:** v3 had an admissible local worldtube but no slice connecting
  it to an exact AE end that could be retained in a parent construction.
- **Decision/unlock:** construct and check that slice, including the wall's
  distributional constraints. Its vacuum end can supply an independently
  justified exterior gluing route; that theorem's hypotheses remain separate.
- **Finite domain:** one fixed wall event Rw/ell=1.5, H ell=1, k ell=2,
  M/ell=0.4. A finite cap, finite vacuum transition, regular bridge, and a
  finite displayed portion of an analytically specified Schwarzschild end.
- **Budget/stop:** at most 60 CPU seconds and 20,000 scalar quadrature
  evaluations; no PDE solver, parameter fit, formation run or source changes.
  Wrong signs, a non-spacelike slice, failed bulk constraints or a failed
  shell delta coefficient stop this connection.
- **Independent checks:** induced metric in Schwarzschild coordinates,
  finite differences of radial geometry, direct Lorentz stress projection,
  Hamiltonian/momentum jump coefficients, mass invariance, quadrature order
  comparison, centre and throat regularity, and unit rescaling.
  This is an application connecting imported components, without a priority
  claim to new thin-shell dynamics or new constraint mathematics.

## Conventions and the slice at the shell

Use signature (-,+,+,+), future slice normal T, outward spatial radial unit
vector s, and K=-Lie_T(h)/2. Write h=dl^2+R(l)^2 dOmega^2 and V=T(R).
The de Sitter cap is a constant flat-slicing time surface:

$$
R'=1,\qquad V_-=HR,\qquad K^i{}_j=-H\delta^i{}_j,
\qquad \rho_-=3H^2/(8\pi),\qquad J_-=0.
$$

It has a regular centre, although its global past and formation are not
specified. At the chosen event the imported wall gives dot R and signed
beta_- and beta_+. With U=gamma(T+v s) and n=gamma(v T+s), its boost relative
to the cap is fixed, rather than fitted:

$$
\gamma={\beta_- -HR\dot R\over1-H^2R^2},\qquad
\gamma v={\dot R-HR\beta_-\over1-H^2R^2}.
$$

Using that same frame across the wall yields

$$
R'_+=\gamma\beta_+-\gamma v\dot R<0,\qquad
V_+=\gamma\dot R-\gamma v\beta_+ .
$$

Thus the slice starts at the centre, expands in areal radius to the wall,
then initially moves toward smaller Schwarzschild radius. It does not silently
substitute the other normal orientation for the imported wall.

## Vacuum transition and bridge

On the decreasing-radius Schwarzschild segment choose a smooth slicing
function V(R), then use the imported mass identity and vacuum momentum
constraint:

$$
R'=-\sqrt{1-2M/R+V^2},\qquad
K^\theta{}_\theta=-V/R,\qquad K^l{}_l=-dV/dR.
$$

Fix R0/ell=1.0 and let V(R)=V_+ S((R-R0)/(Rw-R0)), where S is a smooth step
from zero to one, flat to every order at both endpoints:

$$
S(x)={e^{-1/x}\over e^{-1/x}+e^{-1/(1-x)}}\quad(0<x<1),
\qquad S=0\ (x\leq0),\qquad S=1\ (x\geq1).
$$

This is coordinate slicing freedom in the same vacuum Schwarzschild geometry;
M, H and the physical wall tension remain fixed. It is not interpolation
between incompatible spatial metrics. Substitution into both constraints
shows that no matter is added in the transition.

Below R0, V and K vanish. Since R0>2M, the slice joins a time-symmetric
Schwarzschild neighbourhood before reaching the bridge's minimum sphere.
The regular bridge parameter u is

$$
R=2M\cosh^2u,\qquad
l-l_{\rm throat}=2M(u+\sinh u\cosh u),\qquad R'=\tanh u.
$$

The child-side branch has u<0; the parent-side branch has u>0. The throat at
u=0 is a regular minimum, R=2M, R'=0, R''=1/(4M). The parent side continues
analytically to its positive-mass AE end; a plotting cutoff does not truncate
the mathematical end. A spacelike bridge does not itself establish a causal
communication channel.

## Distributional shell constraints

Let brackets mean exterior minus cap at the shell. The signed wall law gives

$$
[R']=-4\pi R\sigma\gamma,\qquad [V]=4\pi R\sigma\gamma v.
$$

These jumps must be compared with the stress on this slice, not with the wall's
rest density alone. If eta is proper distance normal to the spacetime wall,
then at its intersection with the slice d eta/dl=gamma. Restricting its delta
function therefore gives delta(eta)=delta(l-lw)/gamma. Direct projection of
the pure-wall stress consequently gives integrated slice energy E= sigma gamma
and outward momentum J= sigma gamma v.

With A=K^l_l and B=K^theta_theta, the regular bulk constraints are

$$
{}^{(3)}R+4AB+2B^2=16\pi E_{\rm bulk},\qquad
{2R'\over R}(A-B)-2B'=8\pi J_{\rm bulk},
$$

where scalar curvature is -4R''/R+2(1-R'^2)/R^2. Their singular coefficients
are -4[R']/R and -2[B], respectively, exactly 16 pi E and 8 pi J. A change
in A across the wall does not add a delta to the Hamiltonian constraint;
the A derivatives cancel in the radial momentum constraint. The wall's
tangential pressure remains the fixed v3 value -sigma.

## Recorded connection and independent checks

At Rw/ell=1.5 the fixed wall has dot R=1.6503109509. Its speed relative to
the flat cap is v=0.1351659945, with gamma=1.0092620299. The exterior slice
starts with R'_+=-2.0277860898 and V_+=1.9092537179. Its integrated shell
energy and outward momentum are E ell=0.16062904094 and
J ell=0.02171158406. The corresponding jumps are [R']=-3.0277860898 and
[V]=0.4092537179. Their Hamiltonian and momentum delta residuals are below
1.78e-15 and 3.34e-16, respectively, in the recorded dimensionless units.

The vacuum transition has proper length 0.56555545198 ell. The bridge
minimum occurs at l/ell=2.8977385075 and R/ell=0.8; the displayed parent end
reaches R/ell=6 at l/ell=9.8158051461. These are spatial distances on this
chosen slice, not elapsed times or ages. In particular |R'| exceeds one near
the exterior wall, so a plot of R against l is an intrinsic radial profile,
not a Euclidean surface-of-revolution embedding.

Bulk Hamiltonian and momentum residuals are below 8.44e-15 and 3.56e-15.
Independent finite differences of the geometry and extrinsic curvature at
two radial step sizes give residuals below 9.56e-11 and 4.61e-11. The
Schwarzschild mass changes by less than 1.00e-15 ell across the sampled
vacuum geometry. A direct Schwarzschild metric check gives unit spacelike
tangent and unit timelike normal to 3.56e-15. The 32- and 64-point quadrature
evaluations agree at the stored floating-point precision; this is a numerical
comparison rather than a rigorous error bound. There are 33 cap points,
65 transition points and 97 regular bridge points, including its minimum.

## Scope of the new connection

The slice is regular away from its declared thin wall, and regular at its
centre and bridge minimum. It contains one connected cap, wall, bridge and
AE end. This is initial data in the already borrowed spherical wall model,
not a nonsingular formation history or a stability proof. The cap extends
beyond its de Sitter Hubble sphere and includes anti-trapped spheres; a
regular initial centre does not remove the literature's past-completeness
questions. No inherited constants, antimatter asymmetry or observer-age map
is supplied by the slicing choice.

Reading covered the stated Kuchar and Hajicek sections on 10 October 2026;
the latter's vacuum-shell canonical reduction is not imported wholesale into
the de Sitter bulk. Only its frame and surface-stress machinery is reused.
The focused source search covered spherical Schwarzschild initial data,
thin-shell Hamiltonian constraints and their boosts. It does not certify that
the same elementary connection has never appeared elsewhere.

The [adapter](../src/recursive_horizons/nsc_atlas_connection_slice.py) and
[tests](../tests/test_nsc_atlas_connection_slice.py) own the finite record.

## Retained v4 assembly

The historical v4 application used this sourced body A with $M_A/\ell=0.4$
and a complete vacuum Schwarzschild bridge B with $M_B/\ell=0.6$. Only
selected smooth vacuum ends enter the many-body existence theorem; exact
open overlaps reattach the original bodies. The resulting family has two
AE ends: the designated common parent and B's other end. That other end is
not a second established child. Sufficiently large separation is a theorem
hypothesis, without an explicit common connector or a certified BL distance.
The separate observer at $R_{\rm obs}/\ell=1.6$ is protected only while its
past domain stays inside the retained exterior through $R_{\rm cut}/\ell=4$
and the open overlap. The slice supplies the body-to-vacuum-end interface;
formation, the many-body future, stability and inherited parameters remain
open. Its preserved bundle is historical evidence, with replay selected
explicitly through the reproduction guide.

## Unequal sourced input, 18 October 2026

The [compact two-body record](../results/summaries/nsc-atlas-unequal-children-v1.json)
now supplies a second cap-wall-bridge body. Both bodies use the same
Einstein/conformal-scalar action, classical field `phi=0`, `H*ell=1`,
`k*ell=2`, `xi=1/6`, and wall radius `1.5*ell`. A has `M/ell=.4` and
stops its vacuum slicing taper at `R0/ell=1`; B has `M/ell=.6` and stops
at `R0/ell=1.3`. These stops lie strictly outside the respective throat
radii `.8*ell` and `1.2*ell`. Mass is object initial data; the taper selects
a slice of the same vacuum geometry. The supplied cap vacuum is unchanged.

For B the signed wall quantities have independent closed forms. Set
`s=sqrt(40801)`. In ell units,

$$
\beta_-={151\over120},\quad\beta_+=-{209\over120},\quad
\dot R={s\over120},\quad
\gamma={3s-302\over300},\quad\gamma v={453-2s\over300}.
$$

The measured numerical values are `Rdot=1.683271451008` and local normal-frame
three-speed `v=0.161244881293`. The first is an areal-radius derivative with
wall proper time, not a local three-speed. The checks reconstruct both signed
normals, retain the restricted wall-delta factor, and test the bulk and shell
constraints, mass, spacelike embedding and smooth throat. B's sampled analytic
Hamiltonian residual is below `1.49e-14/ell^2`. Its independent radial-difference
indicator decreases from `8.23e-9/ell^2` to `5.30e-10/ell^2` at the two fixed
steps. These comparisons are not rigorous error bounds.

Both bodies begin with a regular centre and retain the analytic vacuum bridge
only on `u in [u0,+infinity)`, with finite `u0`. The left asymptotic end of the
complete vacuum bridge is excluded by that domain. Each body therefore has
one AE end. End-only gluing of these two sourced bodies leaves one designated
parent end and no additional AE end. The historical v4 vacuum-bridge assembly
above retains its two ends; it is a different input set.

The [parent moments](schur-memory.md#unequal-sourced-inputs-for-the-parent-reduction)
reuse A's closed average and evaluate only B's new moments. Method `f904e0c`
was frozen before the target tests and application. All 42 focused checks
pass; the application used 0.078125 CPU seconds and 13,932 pointwise geometric
samples, without field evolution. The compact record preserves scalars and
source hashes rather than the spatial arrays. A particular shared connector,
formation history, physical inheritance and observational comparison remain
separate calculations.

## Two successive sourced boundaries

The next composition implements two nesting levels in one explicit radial
initial slice. It uses the quiet Einstein/conformal-scalar and positive
membrane action, with `xi=1/6`, `phi=0` and `k*ell=2` at both walls.
The vacuum values are the existing regional choices `3/ell^2` and `0`.
The inner membrane separates two regions of equal vacuum energy; it is
still a material tension source. Its formation is not supplied by that
choice. This is a new hierarchical state of the declared piecewise model,
distinct from P001's nonzero scalar pulse and from the two-sibling parent.

| Outward region | $M/\ell$ | $H\ell$ | Retained portion |
| --- | ---: | ---: | --- |
| Inner cap | 0 | 1 | Regular centre to inner wall at `R/ell=.8` |
| Host Kottler region | .01 | 1 | Inner wall on left static branch, black-hole throat, right static branch, expanding cosmological patch to outer wall at `R/ell=1.5` |
| Outer Schwarzschild region | .4 | 0 | Outer wall on left branch, throat at `R/ell=.8`, right AE end |

For either wall, put $f_\pm=1-2M_\pm/R-H_\pm^2R^2$. The same signed law is

$$
\beta_-={f_--f_++(kR)^2\over2kR},\qquad
\beta_+=\beta_--kR,\qquad \dot R^2=\beta_-^2-f_-.
$$

Choose the positive radial root, with `beta_minus>0,beta_plus<0` at both
initial cuts. The mass relation at each level is

$$
M_+=M_-+\tfrac12(H_-^2-H_+^2-k^2)R^3+kR^2\beta_-.
$$

It relates regional data; it does not determine them all from no inputs.
The inner/outer independent signed factors are respectively
`(517/640,-507/640)` and `(727/600,-1073/600)`. Differentiating the same
law supplies the temporal Israel jump and mean normal-force condition.

Each inner wall frame starts from a flat radial slice, `Rprime_minus=1`
and $V_-=\sqrt{2M_-/R+H_-^2R^2}$. The boost is fixed by
$a=(\beta_-+\dot R)/(1+V_-)$,
$\gamma=(a+a^{-1})/2$, $\gamma v=(a-a^{-1})/2$.
Use the same distributional jumps as above. In a smooth bulk segment an
assigned slicing function $V(R)$ gives

$$
R'=\varepsilon\sqrt{f+V^2},\quad
R''=M/R^2-H^2R+VV_R,\quad
K^l{}_l=-V_R,\quad K^\theta{}_\theta=-V/R.
$$

The host's left slicing taper stops at `R/ell=.1`, strictly above its
black-hole root. A static bridge then reaches the right radius `.5*ell`;
a smooth transition on `[.5,.9]*ell` reaches the expanding flat radial
profile $V=\sqrt{2M/R+H^2R^2}$. This profile passes through the host's
cosmological horizon with `Rprime=1` before reaching the outer wall.
The outer left taper stops at `R/ell=1`, reusing its exact vacuum bridge
and analytic AE continuation. At a simple black-hole root `rb`, use
`R=rb+y^2` with signed `y`; the proper-distance derivative tends to
`2/sqrt(f_prime(rb))`. This avoids treating a coordinate square-root edge
as a singular sphere.

The [bounded application](prior-art-and-open-claim.md#two-successive-sourced-boundaries-18-october-2026)
now has a [compact checked record](../results/summaries/nsc-atlas-nested-slice-v1.json).
It retains scalar milestones and integrated
moments, rather than sample arrays. A regular initial slice does not yet
provide the matched Kottler horizon charts or a common future clock.
The spatial chain alone implies neither a causal passage across the bridge
nor self-similarity of a physical operator. Establishing
$\mathcal T^*\mathbb D=\mathbb D$ for its field/metric response requires
the actual transition maps and domains. No fractal dimension or cosmological
Lambda is assigned by these inputs.

### Regular host charts and the clock convention

An explicit bulk chart is available before any field solve. Choose a simple
positive horizon root `C` of `1-2M/R-H^2 R^2`, and put

$$
\kappa_C={1-3H^2C^2\over2C},\quad
R=C+\kappa_Cuw,\quad
F={w^2\over R}\left({1\over C}+H^2(R-C)\right),
$$
$$
ds^2=Fdu^2+2du\,dw+R^2d\Omega^2.
$$

This rewrites [Lake, gr-qc/0507031v2, equations (5)--(10)](https://arxiv.org/html/gr-qc/0507031v2)
using `H^2=(C-2M)/C^3`. The two-dimensional determinant is `-1` and
the coefficients are regular for finite coordinates with `R>0`, including
`R=C`. Direct symbolic curvature gives
$R_{\mu\nu}=3H^2g_{\mu\nu}$, and
$|\nabla R|^2=1-2M/R-H^2R^2$; all four independent residual expressions
reduce to zero in a 0.5 CPU-second algebra check. No target field or
trajectory is evaluated by that check.

Use a black-hole-centred chart for the bridge and a cosmological-horizon
chart where the first chart reaches an infinite coordinate. On a static
overlap they match through the same areal radius and

$$
T=r_*(R)-\kappa_C^{-1}\log(|u|/L_0),\qquad r_*'=1/f.
$$

Branch signs and the arbitrary coordinate scale `L0` must be fixed by the
actual slice. They are not physical ages. In this single Killing coordinate
the initial oriented frame is

$$
s=(V/f,R'),\qquad N=(R'/f,V),\qquad R'^2-V^2=f.
$$

Thus `T_l=V/f` on both static branches. The left future has `N(T)<0`;
replacing its coordinate by a locally future time `T_L=-T` reverses both
time formulas. It does not change the geometry. On the right expanding
profile the outgoing coordinate obeys
`u_EF,l=-1/(V+Rprime)`, which stays finite at the cosmological horizon.
Matching the wall worldtubes and chart overlaps remains part of the shared
spacetime construction. The regular bulk formula alone supplies neither
global identifications nor a nonsingular replacement for `R=0`.

### A central test-flash protocol

For the specified inner wall, prescribe an ideal local frequency standard
and a radial inward test flash at its initial event. Use geometric optics
with negligible backreaction. Its receiver is the cap's central geodesic;
set cap flat-slicing time `t=0` at emission. In that cap, `t` is central
proper time and `Rwall=.8*ell`. The exact null path and wall boost give

$$
t_{\rm rec}=-H^{-1}\log(1-HR_w),\qquad
1+z={\gamma(1+v)\over1-HR_w}
={\beta_-+\dot R\over(1+HR_w)(1-HR_w)}.
$$

The Doppler factor belongs to the material wall, rather than to a comoving
cap emitter. Prescribing a frequency standard is a test protocol; the
tension action supplies no emission line. The host enters through the
signed junction factors, while `H` and masses remain supplied model data.

The photon stays inside the wall. In flat cap coordinates it obeys
`dr/dt=-exp(-H*t)`; a future timelike wall has strictly greater radial
coordinate speed. For this inner root in ell units,
`Rdot^2=2R^2-1+.01/R+.000025/R^4` is positive and increasing for `R>=.8`,
and `Rdot<=sqrt(2)*R`. The expanding wall persists to every finite cap
time and cannot recross that inward ray. This checks the reception domain
without a numerical wall evolution. It implies no outer-wall signal
through the host bridge. The finite application will record the conditional
arrival time and redshift with their inputs and numerical precision.

The implementation is [the nested slice](../src/recursive_horizons/nsc_atlas_nested_slice.py),
with [regular charts and the flash](../src/recursive_horizons/nsc_atlas_nested_spacetime.py).
Twenty focused controls pass. They include the independent symbolic Einstein
calculation, rational junctions, cap and Schwarzschild limits, horizon and
length controls. The [producer](../scripts/record_nsc_atlas_nested_slice.py)
requires frozen source bytes and a fresh output path. The
[figure producer](../scripts/plot_nsc_atlas_nested_slice.py) regenerates the
radius curve rather than retaining spatial arrays. These bulk charts and
cap-side observations do not yet assemble the joint two-wall clock atlas.

### First application and finite result

Frozen method `acc5a176ac074999ed23af770045003c4123cd57` passes the first
application in 0.078125 CPU seconds / 1,016 charged geometric points.
The result is 11,474 bytes. The figure regenerates 10,087 points in
0.515625 CPU seconds; neither command solves a field or evolves a grid.
The host black-hole and cosmological radii are
`0.020008009615388216*ell` and `0.989845863754293*ell`. The two walls
occur at proper distances `.8*ell` and `3.372632651397006*ell`.
The right retention sphere is at `l/ell=9.503059823421234`.

![Two successive sourced boundaries and the resolved inner throat](figures/nsc-atlas-nested-slice.png)

This is a proper-radius versus areal-radius profile of the initial geometry.
Its two minima are regular spatial spheres. It does not establish causal
travel across either bridge.

The central flash gives `t_rec/ell=1.6094379124341005` and
`1+z=3.7463933761061194`. At the supplied ruler `ell=10 km`, the central
proper interval is **53.6850701 microseconds**. Independent 80-digit
evaluation of the rational wall factors differs by `1.30e-16` in time/ell
and `1.45e-15` in `1+z`. These are roundoff indicators, not model or
experimental error bars; input uncertainties were not supplied.
The host changes the wall boost while the membrane law remains fixed.
This is a conditional observer prediction for a specified test source.

The [initial normal-clock average](schur-memory.md#two-successive-sourced-levels)
passes its global constraint check. Its geometric trace is negative and
its stress fails the vacuum-form test, so it does not generate a cosmological
Lambda. This finite calculation is closed. Further work concerns the parent
geometry, states and observer reduction, rather than a precision scan here.
