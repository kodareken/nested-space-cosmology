# Local structure and the inherited parent need different length scales

The [fixed-observer image result](nsc-atlas-core-images.md) ties a local
delay to the parent's vacuum through one dimensional embedding. The next
question is whether the same field law admits local structures much smaller
than that parent. This calculation supplies an exact initial-constraint
family and a conditional asymptotic existence argument. Its large-frequency
members have not undergone the preceding finite evolution or lens experiment.

## The connection chosen

Three routes were considered: more image refinement, direct large-frequency
evolution, and a scale/energy-supply construction. The completed image pair
already has a stated numerical range. Direct evolution would hide a missing
preparation constraint behind a large scale ratio. The supply calculation
instead determines which initial configurations can be assembled at all.

Reuse the original canonical scalar/four-form potential, the nodeless core
shape, spherical gravitational constraints and protected-parent theorem.
Keep kappa=1, epsilon=1/4, fluxes (.5,1,1.5), charge .5, tension .2, wall
radii 1 and 1.9, and middle mass .12. Frequency omega and core amplitude zeta
are declared physical parameters. The reservoir's depletion a is supplied;
its outer radius is solved from zero total exterior mass offset. There is
no added outgoing pulse and the scalar normal momentum is zero initially.
This is a distinct preparation from the diffuse cloud used for the image pair.

`nsc_atlas_scale_hierarchy.py` imports the existing core and parent owners.
It introduces scaled core and reservoir coordinates, an independent
proper-distance Hamilton calculation and a solved reservoir radius. These
separate the two physical lengths instead of hiding them in one grid.
The producer uses thirteen independent cases in six processes with a 120
CPU-second budget. Failure of the capacity bound, chart positivity or
constraint matching has a specific diagnostic; no field evolution is run.

The existing [core owner](nsc-atlas-scalar-core.md) attributes the
Schrodinger--Newton and real-scalar envelope input, including
[Fodor, Forgacs and Mezei](https://arxiv.org/html/1007.0388). The spherical
mass law is the retained [Hayward construction](https://arxiv.org/html/gr-qc/9710089).
The new result is the exact scaled constraint, the supply obstruction and
the compatible hierarchy in this excited inherited parent.

## What remains unchanged when frequency changes

With q_g=.5 and mu=omega sqrt(1-epsilon), choose the common initial field
phi0=-sqrt(1-epsilon) q_g/omega. Then

\[
U(\phi_0,q)=\frac{\epsilon q^2}{2}
+\frac{(1-\epsilon)(q-q_g)^2}{2}.
\]

Thus the original initial stresses, wall forces, expansions and cap
curvature do not depend on omega. In the parent,

\[
\Lambda_v=\kappa U_v=\frac{9}{32},\qquad
U_0=\frac{21}{32},\qquad E_b=U_0-U_v=\frac{3}{8}.
\]

The reference radius stays 6.62607707610 and Lambda_v Lref^2 stays 12.34825239893.
Changing omega changes the physical action and subsequent scalar history.
It is not a unit conversion or identical parent evolution at all frequencies.
Every nested region in one member uses that member's common action.

## Exact core mass as three inverse powers of frequency

Use the existing compact core shape A(x), including its C2 taper at x=20:

\[
x=\omega\zeta R,\qquad
f=\sqrt{2/\kappa}\,\zeta^2 A(x),\qquad
\mathfrak{M}=\frac{\zeta}{\omega}N(x).
\]

Here mathfrak M is geometric excess mass above the excited homogeneous
parent. The background displacement from the parent minimum is d/omega,
where d=sqrt(1-epsilon)(q_parent-q_g)=sqrt(3)/2. Write C for the unchanged
initial parent curvature. The exact constraint reduces to

\[
N'=-\zeta^4x(A')^2N
+\frac{x^2}{2}\left[A^2+\zeta^2(A')^2\right]
+\frac{x^2\sqrt{\kappa/2}\,d A}{\omega\zeta^2}
-\frac{Cx^4(A')^2}{2\omega^2}.
\]

Its coefficient multiplying N is frequency independent. Linearity and
N(0)=0 therefore prove the **exact** decomposition, for this supplied shape,

\[
N(x;\omega)=N_0(x)+\frac{N_1(x)}{\omega}
+\frac{N_2(x)}{\omega^2},\qquad
\mathfrak{M}_c=\frac{\zeta N_0(20)}{\omega}
+\frac{\zeta N_1(20)}{\omega^2}
+\frac{\zeta N_2(20)}{\omega^3}.
\]

This is not a fit to a frequency scan. The first coefficient is positive
for a nonzero shape; the second carries the ambient/core interference and
the third the parent-curvature correction. At the retained zeta=.0826621713,
the calculated coefficients are approximately 1.756636312, 1118.302431 and
-.1349602551. Independent integration of the original areal-radius core
equation at omega 2000 agrees in geometric mass within 2.7e-16.

## A reservoir that scales with the core eventually fails

Let Rc=20/(omega zeta). In a reservoir with Pi=0 and positive radial chart B,
the exact mass equation is

\[
\mathfrak{M}_R=\frac{\kappa R^2}{2}
\left[\frac{B F^2}{2}+U-U_0\right]
\geq-\frac{\kappa E_bR^2}{2}.
\]

To reduce a positive core mass Mc to zero by radius Rd requires

\[
R_d^3-R_c^3\geq\frac{6\mathfrak{M}_c}{\kappa E_b}.
\]

Gradient energy only makes this requirement stronger. In an untrapped
region with B>=v^2 it also follows for arbitrary Pi by completing the square:

\[
\rho-v\Pi F=U+\frac{(\Pi-vF)^2}{2}
+\frac{(B-v^2)F^2}{2}\geq U_v.
\]

For this inherited charge ladder, kappa E_b=4 Lambda_v/3. The result is
the shape-independent, dimensionless bound

\[
\mathcal Q=\frac{\Lambda_v(R_d^3-R_c^3)}{\mathfrak{M}_c}>\frac{9}{2}.
\]

The strict inequality applies to the smooth finite-energy matches here.
Equality would require zero derivative cost and U=Uv almost everywhere
through the annulus. The scalar would then equal its parent minimum
throughout, contradicting its continuous matching to the excited ambient
value at the two ends. This lower bound is independent of the bump and its
depletion. The later value 18.330588 is a particular constructed family
above that floor. The bound is for the declared inherited excitation/vacuum
ratio; a different initial parent state changes its coefficient.

For Rd=b Rc with fixed b, the available excitation scales as omega^-3,
while the positive core mass scales as omega^-1. Such a fixed-proportion
preparation cannot persist to arbitrarily large frequency. For b=2 and
the retained zeta, the exact core polynomial gives the necessary ceiling
omega=6221.94. Above it even completely removing the available excitation
throughout that annulus is insufficient. This is an obstruction to that
preparation, not to the inherited field model or to larger reservoirs.

## Constructing the larger reservoir

Set Pi=0 and place a depleted scalar annulus between Rc and Rd:

\[
\chi(R)=\frac{d}{\omega}\,[1-a b(x)],\qquad
x=\frac{R-R_c}{R_d-R_c},\qquad
b(x)=1024x^5(1-x)^5.
\]

Take 0<a<=1. The compact bump is C4 at both joins. The scalar and its
derivatives recover the original parent at Rd. Solve the exact mass
constraint with Mc as its inner datum, selecting Rd by M(Rd)=0. Outside
Rd the fields, geometric mass and initial metric are exactly those of the
original parent. Both initial material walls consequently retain their
original matching data.

For clarity, put W=Rd-Rc and y=M/Mc. The reservoir equation is the affine
ODE

\[
y_x=-\frac{\kappa RF^2W}{2}y
+\frac{\kappa R^2W}{4\mathfrak{M}_c}
\left[(1-CR^2)F^2+d^2(a^2b^2-2ab)\right].
\]

Now fix zeta and a, write Mc=Cm/omega+O(omega^-2) with Cm>0, and try
Rd=z omega^(-1/3). The gradient source vanishes as omega^(-4/3), the
coefficient of y vanishes as omega^-2, and Rc/Rd vanishes as omega^(-2/3).
The endpoint equation converges to

\[
y(1)=1-\frac{\kappa d^2z^3}{4C_m}J(a),\qquad
J(a)=\int_0^1x^2[2ab(x)-a^2b(x)^2]dx>0.
\]

Its positive root has a nonzero derivative with respect to z. The
implicit-function theorem therefore continues it to sufficiently large
finite omega. The scaled affine ODE has a smooth limiting coefficient
system; this argument does not presume a numerical fit. Assume the
limiting core has a positive B margin. For sufficiently large omega the
reservoir's potential deficit dominates its gradient cost, M decreases
from Mc to zero, B remains positive, and Rd lies inside the material wall.
These conditions are directly checked for the finite examples.

Hence this initial-data family has

\[
R_c\sim\omega^{-1},\qquad
R_d\sim\omega^{-1/3},\qquad
\frac{R_d}{R_c}\sim\omega^{2/3},\qquad
R_d\ll L_{\rm ref}.
\]

The parent's causal energy-floor bound supplies a protected annulus for
each constructed member. Its future homogeneous scalar history depends
on omega; its inherited vacuum and the uniform causal bound do not.
This uses the existing domain-of-dependence result, not a completed future
of both material walls.

## A dimensionless mass--reservoir--vacuum relation

The hierarchy gives a joint relation that survives a choice of physical units:

\[
\frac{\Lambda_v R_d^3}{\mathfrak{M}_c}
\longrightarrow\frac{4\Lambda_v}{\kappa d^2J(a)}.
\]

For the declared half-depletion a=1/2,

\[
J(1/2)=\frac{5476736}{66927861},\qquad
\frac{\Lambda_vR_d^3}{\mathfrak{M}_c}
\longrightarrow\frac{200783583}{10953472}=18.3305880546.
\]

Here physical geometric mass means G M_phys/c^2. The coefficient depends
on the supplied depletion and shape. It is not a newly selected constant of
nature. Quarter- and three-quarter-depletion examples retain this freedom
explicitly; the same core can have different admissible supply radii.

The [compact result](../results/summaries/nsc-atlas-scale-hierarchy-v1.json)
contains thirteen examples and controls. At fixed retained zeta, increasing
omega from 2000 to 1e10 changes Rd/Rc from 1.5743 to 40578.7 while preserving
Lambda_v Lref^2=12.34825239893. The last vacuum-radius-mass product is
18.33027699, approaching the independently derived 18.33058805.

An illustrative small-amplitude member uses omega=1e9, zeta=.001, a=.5.
Assigning the **supplied reference** Lambda_v=1e-52 m^-2 gives:

| Quantity | Conditional dimensional value |
| --- | ---: |
| Core geometric excess mass | 6.32e10 nominal solar masses |
| Radius enclosing half the core geometric excess | 3.86 kpc |
| Compact core cutoff | 34.4 kpc |
| Compensating reservoir outer radius | 828 kpc |

The length and nominal mass conversions use the
[IAU nominal constants](https://www.iau.org/static/resolutions/IAU2015_English.pdf).
The reference Lambda and the field/preparation parameters were supplied,
not predicted or fitted to a galaxy. These numbers show that the constraint
construction admits widely separated local and cosmological scales. They do
not establish a stable galaxy, its rotation curve, or an observed abundance.

## Verification and next connection

The symbolic core-polynomial residual is exactly zero. Scaled and original
core masses agree; the scaled reservoir reproduces the original no-pulse
preparation's solved depletion. A proper-distance system independently
evolves the radial metric deficit from the Hamilton equation. Its exterior
mass residual is below 4.2e-13 of core mass across the retained examples.
Tightening integration tolerances changes the reported mass/radius quantities
by at most 1.1e-11 relatively. These are finite integration checks; the
input core shape retains its earlier numerical profile error and its
leading-envelope origin. No rigorous numerical interval or full evolution
error bound is asserted.

The next physical problem is evolving this **two-scale preparation** with
its supply region retained, and testing how much mass remains localized.
The previous omega 2000 diffuse lens is unchanged. Its image delay cannot
be transferred to this family by a unit conversion. The broad construction
has now identified what that evolution must include, rather than requesting
an unexplained source or a finer version of the old grid.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_scale_hierarchy.py --output results/summaries/nsc-atlas-scale-hierarchy-v1.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_scale_hierarchy.py
```

Grok independently confirmed the supply-floor argument, the failure of fixed
radius ratios and the larger-reservoir construction. Its leading weak-amplitude
mass estimate is distinguished from the exact three-coefficient constraint
used here, which retains the finite gradient and self-gravity corrections.
