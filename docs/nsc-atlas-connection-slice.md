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
