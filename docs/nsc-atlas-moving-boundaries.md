# Moving boundaries: conservation and a separate vacuum-wall control

This successor preserves the [v2 junction screen](nsc-cosmological-atlas.md).
It asks whether boundary motion alone restores that rejected attachment, then
imports a different published moving wall to keep the surviving map visible.
Neither calculation is a claim to have discovered thin-shell dynamics.

## Admission recorded before execution

- **Borrowed:** Lobo and Crawford, [gr-qc/0507063v2](https://arxiv.org/html/gr-qc/0507063v2),
  Section 2, equations (3), (4), (8)--(10), and Section 3.2, equations
  (24)--(28): surface geometry, Israel stresses and the contracted Codazzi
  flux condition. Ansoldi, Guendelman and Shilon,
  [0711.2198v1](https://arxiv.org/html/0711.2198v1), Section 2, equation (2),
  and Section 3: signed de Sitter/Schwarzschild vacuum-wall dynamics.
- **Missing:** v2 imposed comoving radiation. Could allowing a noncomoving
  boundary repair it without changing its fixed pure tension or vacuum
  exterior? Which different imported bulk actually permits a moving wall?
- **Decision and unlock:** conservation can stop motion as a repair within
  the radiation candidate. A separate vacuum control supplies an admissible
  local moving worldtube, with explicit orientation; attachment to the
  Brill-Lindquist parent remains a separate missing connection.
- **Finite domain:** timelike spherical boundaries in classical four-dimensional
  Einstein gravity, geometric units. The radiation argument is local and
  algebraic for positive density and relative speed less than light speed.
  The control uses one fixed parameter tuple and 33 radii in a finite
  monotone interval, away from turning points, zero signed normal components
  and static-coordinate horizons. No formation trajectory is computed.
- **Budget and stop:** fewer than 10,000 scalar quadrature evaluations and
  60 CPU seconds; no parameter search or dynamics solver. A nonzero radiation
  flux retires motion alone as the repair. Any failure of the unsquared
  junction or the independent geometric check stops the control.
- **Verification and attribution:** direct Lorentz-frame stress contraction,
  signed angular and pressure junctions, metric normalization and independently
  differentiated acceleration, quadrature-order comparison, fixed parameters
  and unit rescaling. The result is a compatibility application and imported
  control. The bounded search cannot establish worldwide priority.

## Conservation decides the proposed radiation repair

Use signature (-,+,+,+), a shell tangent U, and one unit normal n directed
from the FLRW side to the exterior. A pure wall with constant positive
tension has S^a_b = -sigma delta^a_b and hence D_a S^a_b = 0. Contracted
Codazzi therefore requires [T(U,n)] = 0. This zero-jump statement is
independent of an overall convention for the normal or flux sign.

In the radiation rest frame let u be its unit timelike velocity and e the
outward radial unit vector. Write

$$
U=\gamma(u+v e),\qquad n=\gamma(vu+e),\qquad
\gamma=(1-v^2)^{-1/2}.
$$

The perfect-fluid stress gives

$$
T(U,n)=(\rho+p)\gamma^2v,
\qquad p=\rho/3,\quad \rho>0.
$$

Schwarzschild vacuum has zero flux. Thus constant pure tension requires
v=0. Normal reversal reverses the displayed flux, leaving that conclusion
unchanged. This is a necessary local condition, not a sufficient matching
condition. Combined with v2's comoving-wall obstruction, noncomoving motion
alone cannot repair that radiation/vacuum candidate on its nonconstant-radius
interval.

The word "transparency" sometimes denotes zero *net jump*, not a material
rule saying particles pass through a wall. Radiation passing into an otherwise
vacuum exterior changes the exterior stress tensor. Absorption, reflection,
surface fields or a radiating exterior require a separately specified stress
and exchange law; none is smuggled into this test. A non-perfect-fluid boundary
layer likewise lies outside this argument.

## Distinct imported moving control

Replace the radiation bulk by de Sitter vacuum, keeping a constant positive
wall tension and a Schwarzschild exterior. Now rho+p=0, so normal fluid flux
vanishes for any timelike boundary speed. This changes the candidate's bulk;
it does not rescue the radiation candidate while keeping its theory fixed.

For f_- = 1-H^2 R^2, f_+ = 1-2M/R, and k=4 pi sigma > 0, the published signed
junction is

$$
\beta_- - \beta_+ = kR,\qquad
\beta_\pm=\epsilon_\pm\sqrt{\dot R^2+f_\pm}.
$$

Its algebraic form retains the signs explicitly:

$$
\beta_-={M\over kR^2}+{k^2-H^2\over2k}R,\qquad
\beta_+={M\over kR^2}-{k^2+H^2\over2k}R,
\qquad \dot R^2=\beta_-^2-f_- .
$$

The adapter selects the expanding square root only after verifying the
unsquared equation. It evaluates proper time by integrating dR/sqrt(dot R^2).
Pressure is checked separately through
K^theta_theta = beta/R and K^tau_tau = (ddot R+f'/2)/beta.
Replacing signed beta by its absolute value is an explicitly tested error.
The squared radial speed has the form a/R^4+b/R+c R^2-1, with a,c positive.
Its unique positive stationary radius is its minimum; the adapter checks it
when it lies between requested endpoints. Quadrature therefore cannot skip a
narrow forbidden interval just because all its nodes happen to miss it.

The finite control fixes H ell=1, k ell=2, M/ell=0.4 and R/ell in [1.2,1.8].
Both signs stay fixed: epsilon_- = +1, epsilon_+ = -1. The outward normal
therefore points toward decreasing Schwarzschild areal radius. This is the
orientation that the published child-universe discussion associates with
growth on the far side of a throat. Our local segment alone does not establish
its global completion, formation, event horizons or communication.

The de Sitter interval lies beyond its static-patch horizon, so its static
coordinate t is spacelike. Proper shell time remains timelike. On the selected
Schwarzschild patch the signed time coordinate decreases along the chosen
future shell tangent; it must not be substituted for an observer's age. The
adapter checks the indefinite metric directly and does not integrate either
coordinate time. No cross-universe age or entropy arrow is inferred.

## Recorded finite result

For the fixed control, proper shell duration is 0.371870610298428 ell while
R/ell increases from 1.2 to 1.8. The endpoint radial rates are 1.2325948740
and 2.0574200012; these are areal-radius rates per proper time, not local
three-speeds. The fixed source is sigma ell=1/(2 pi)=0.15915494309189535
and P ell=-sigma ell. Signed beta_- runs from 1.0388888889 to 1.4117283951,
and beta_+ from -1.3611111111 to -2.1882716049.

The maximum unsquared angular residual is 4.45e-16 and the direct pressure-law
residual is below 2.78e-17 in dimensionless units. A separate five-point
differentiation of the tangent followed by the metric Christoffel contraction
checks the pressure to 2.86e-12 and K^tau_tau to 7.60e-11. Tangent normalization
agrees with -1 to 3.56e-15 even where the de Sitter static coordinate is
spacelike. The 32/64-point quadrature comparison differs by 5.56e-17 ell;
this agreement is a numerical check, not a rigorous quadrature error bound.
Discarding the beta signs instead leaves an angular residual of at least
2.7222, so the orientation is consequential.

## Reading depth and boundaries

The cited source sections and equations were read directly on 10 October
2026. The 0711.2198 discussion distinguishes singular expanding histories,
bounded histories and bounce histories, and stresses the sign data that a
potential alone loses. The original Blau--Guendelman--Guth paper is attribution
in that source, not independently reread here. No formula from the corrected
0803.0268 potential is used.

Aguirre and Johnson, [gr-qc/0508093v2](https://arxiv.org/html/gr-qc/0508093v2),
Sections II.1--II.3, equations (12)--(13), independently supplies the signed
junction and demonstrates why chart-time directions and global regions must
be recorded. Our beta means n^a partial_a R; coordinate-time signs depend on
the chosen chart convention. The paper's general nonspherical stability
investigation is not a stability certificate for this finite control.

Focused arXiv queries covered de Sitter/Schwarzschild signed thin shells and
FLRW/radiation noncomoving-wall flux. They also exposed other shell equations
of state and perturbation papers; these were not absorbed as solutions to the
present gap. No priority claim follows from that bounded coverage.

The [module](../src/recursive_horizons/nsc_atlas_moving_boundaries.py) and
[tests](../tests/test_nsc_atlas_moving_boundaries.py) own the adapter. Its record
keeps the radiation incompatibility and the vacuum control separate. Neither
is embedded in the shared two-hole initial slice. Formation, nonspherical
matching, perturbative stability, inherited constants and sibling influence
remain unresolved.
