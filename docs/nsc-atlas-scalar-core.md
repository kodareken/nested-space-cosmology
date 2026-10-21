# A self-gravitating scalar profile and its inherited environment

Self-gravity supplies a physical profile-selection mechanism beyond the
polynomial packet used in [NSC-P003](nsc-atlas-parent-clocks.md). Two
different localized clouds in a leading real-scalar envelope model shed
radiation and leave cores within about 4% of the independently calculated
nodeless profile after normalization by their remaining mass. A coherent
representative of one late core can then be admitted into the compact
atlas's exact initial-constraint equations, including both material walls.

These are connected but distinct results: a local vacuum-limit cooling
calculation and a sourced general-relativistic initial slice. Their full
dynamical matching in the excited parent remains open. The initial cloud,
its energy and the scalar frequency are specified; this does not derive
primordial structure or retroactively prepare P003's particular profile.

## Choice of route and mathematical attribution

Three preparation routes were compared. The protected parent cannot be
reached by its existing walls after the retained initial slice. A Gaussian
state or Euclidean boundary condition would add a state prescription, and
the earlier real O(4) saddle is a different geometry. Scalar self-gravity
instead gives an internal classical mechanism for selecting a localized
profile from already supplied cloud data.

[Seidel and Suen](https://arxiv.org/abs/gr-qc/9309015) studied gravitational
cooling of self-gravitating scalar configurations. The weak-field
stationary equations are also the Schrodinger--Newton system studied by
[Moroz, Penrose and Tod](https://doi.org/10.1088/0264-9381/15/9/019).
[Harrison, Moroz and Tod](https://arxiv.org/html/math-ph/0208046) give
time-dependent radial examples and discuss outgoing radiation and
rescaled residual cores. We independently solve the stationary and
time-dependent equations here.

For the relation to a real Klein--Gordon field with positive vacuum
curvature, use [Fodor, Forgacs and Mezei, Section III](https://arxiv.org/html/1007.0388).
The complex envelope encodes the amplitude and phase of one real scalar.
It is not a new charged field, a probability interpretation or a quantum
state-reduction law. A strictly isolated periodic lump is not inferred in
positive Lambda; late radiation and the full relativistic problem retain
their published limitations.

The new connection is the explicit hierarchy test, local finite cooling
with conserved-flux checks, and the admission of its mass-selected coherent
representative into this shared atlas. The ambient inherited field is
retained in that admission rather than added afterward as a separate fluid.

## An exact frequency family, not a change of units

Keep kappa=1, epsilon=1/4, charge e=1/2, tension sigma=1/5 and the three
fluxes (1/2,1,3/2). For each frequency omega, set mu=omega sqrt(1-epsilon)
and the common initial field to the grandchild minimum. Then

\[
\phi_0=-\frac{\sqrt{1-\epsilon}\,q_g}{\omega},\qquad
U(\phi_0,q)=\frac{\epsilon q^2}{2}
 +\frac{(1-\epsilon)(q-q_g)^2}{2}.
\]

The initial energies, wall forces, signed wall jets and original shared
geometry are independent of omega. Their later scalar dynamics are not.
This family changes an actual coupling; it is not a rescaling of units or
a derivation of that coupling from the initial geometry.
The wall radii remain 1 and 1.9, with common middle mass .12. Numerical
lengths and geometric masses use these same atlas coordinate units;
no SI ruler is selected.

Let zeta denote the weak-field core amplitude parameter, distinct from
the mixing epsilon. The local approximation needs zeta^2 small and

\[
h=\frac{H_\Lambda}{\omega\zeta^2}\ll1,
\qquad H_\Lambda=\sqrt{\kappa\rho_v/3}.
\]

For rho_v=9/32 and zeta=.1, the original omega=20 gives h=1.53093;
it is outside this small-h core approximation. The explicitly different
common-frequency member omega=2000 gives h=.0153093 and zeta^2=.01.
It is chosen to exhibit that hierarchy, not fitted to observations.
This comparison does not rule out relativistic structures at omega=20.
P003 remains the separate omega=20 conditional result.
The scale ratios are not a certified full-GR error estimate.

## The profile follows from the coupled field and gravity equations

In local coordinates rho=omega zeta r and tau=omega zeta^2 t,

\[
\chi=\sqrt{2/\kappa}\,\zeta^2
       \operatorname{Re}[\psi(\rho,\tau)e^{-i\omega t}],
\]

\[
i\partial_\tau\psi=-\frac12\Delta\psi
       +(V-h^2\rho^2/2)\psi,\qquad
\Delta V=\frac12|\psi|^2.
\]

These are leading slow-envelope equations in the local de Sitter vacuum
limit. They do not replace the full evolving parent by an exact vacuum.
For h=0, the stationary nodeless shape obeys

\[
S''+2S'/\rho+sS=0,\qquad
s''+2s'/\rho+S^2=0,\qquad s(\infty)=-1,
\]

with regular pole derivatives and S positive. Here V=-(s+1)/2 and the
envelope eigenvalue is -1/2. Independent quadrature and asymptotic mass
give

\[
S(0)=1.02149303621,\qquad
I=\int_0^\infty\rho^2S^2d\rho=3.50532970334.
\]

The half-mass and 99%-mass radii are 2.23950479521 and 5.67468871483 in rho.
The independent radial kinetic and gravitational energies are .584221617046
and -1.168443234123; their virial residual is below 3e-11.
Changing the outer solve interval from 25 to 20 changes the central amplitude
by 8.2e-11 and I by 5.2e-10.

The scaling S_lambda(rho)=lambda^2 S(lambda rho) gives mass proportional
to lambda and radius proportional to 1/lambda. Restoring the geometric mass
offset mathfrak M=G M in units c=1,

\[
\mathfrak M=\frac{\zeta I}{2\omega},\qquad
R_f=\frac{\rho_f}{\omega\zeta},\qquad
\mathfrak M R_f\omega^2=\frac{I\rho_f}{2}.
\]

The shape and dimensionless mass-radius relation are determined. The
frequency and the total amount of scalar material remain physical inputs.

There is an energetic answer to why the isolated reference is spherical.
At fixed envelope norm, removing phase gradients and taking a symmetric
decreasing rearrangement cannot increase the kinetic-plus-self-gravity
energy; see [Lieb's minimizing-solution result](https://doi.org/10.1002/sapm197757293).
For a norm-preserving trial dilation, the energy is b^2 T+b W, with T>0
and W<0. Its minimum is at b=-W/(2T), giving the virial relation and a
finite size. This attributed variational statement concerns the isolated
h=0 problem. Positive vacuum curvature adds an outward tidal term; it
does not give a global minimum of the same functional on all of space,
and our radial cooling runs do not establish general relativistic
nonspherical stability.

## Finite gravitational cooling, with the outflow accounted for

Use Gaussian and super-Gaussian clouds of width 6, both with initial radial
norm I and zero envelope phase gradient. Evolve to tau=600 with h=.0153093.
The kinetic and self-gravity steps preserve norm; an outer numerical
absorber starts at 70% of the finite radius. It acts only in the outgoing
region. Its removed norm and Hamilton energy are recorded explicitly.

For u=rho psi, the independent current is

\[
J=\operatorname{Im}(u^*u'),\qquad
\partial_\tau|u|^2+\partial_\rho J=0
\]

in the interior where the absorber vanishes. Integrating J at rho=20 checks that norm actually
leaves the central region. It is not merely disappearing at a boundary.
This norm is proportional to real-scalar mass at the retained envelope
order; it is not an exact conserved U(1) charge of the original real field.

| Calculation | Initial profile discrepancy | Late discrepancy | Norm remaining in finite domain |
| --- | ---: | ---: | ---: |
| Gaussian cloud | .81467 | .03767 | 2.90161 |
| Super-Gaussian cloud | .79474 | .04113 | 3.01613 |

The late discrepancy compares the square root of time-averaged density over
500<=tau<=600 with the nodeless shape. The remaining mass inside rho=10
fixes its scale, without fitting the peak or radial shape. The cores retain
roughly 83% and 86% of the initial norm and still have residual oscillations.
This is finite evidence for the two preparations, not a universal attractor
or asymptotic-convergence theorem.

Halving the step and radial spacing changes the Gaussian mean inner mass
by 7.3e-6. Moving the outer boundary from 64 to 96 changes it by .004781,
about .165%. The corresponding profile discrepancy changes by .00217.
Fine-run Hamilton-energy accounting errors remain below 3.4e-6, and the
independent inner-flux balance error is below 6e-6. Norm bookkeeping closes
below 2e-9. These finite comparisons are separate from relativistic and
ambient-field errors.

A stationary control remains close to its independently solved profile;
a no-gravity control agrees with the analytic dispersing Gaussian. The
six independent runs use six worker processes. Full radial fields are
hashed, not committed; the figure uses transient in-memory projections.

![Scalar cooling and the independently solved profile](figures/nsc-atlas-scalar-core.png)

## Admit the coherent core into the actual initial atlas

The Gaussian late inner mass selects scale lambda=.826621713254, hence
zeta_core=.0826621713254. Use that nodeless coherent representative at a
turning point, not the full radiating snapshot. The difference from the
time-averaged evolved core is the measured 3.77% above.

Add its real displacement f(R) to the inherited parent field phi0. In the
leading local profile f=sqrt(2/kappa) zeta_core^2 S(omega zeta_core R).
A quintic taper acts only between rho=15 and20, with value and its first
two derivatives matched. It reaches exactly zero before the material
wall; no additional surface source is inserted. Isotropic-to-areal radius
corrections belong to the envelope approximation, while the following
initial constraint equations are solved for the actual supplied f(R).

Retain homogeneous parent slice expansion v=H R and zero scalar momentum.
Write the geometric mass as m=M(R)+kappa U0 R^3/6, where U0=21/32 is the
parent's initial ambient density. Then

\[
\beta^2=1+H^2R^2-\frac{\kappa U_0R^2}{3}-\frac{2M(R)}R,
\]

\[
M'=\frac{\kappa R^2}{4}\beta^2f'^2
 +\frac{\kappa R^2}{2}U'_0f
 +\frac{\kappa R^2}{4}\omega^2f^2.
\]

The scalar cross term U'_0 f belongs to the same field, not an additional
fluid. M(0)=0 gives a regular pole. The momentum constraint vanishes
because Pi=0 and both slice extrinsic eigenvalues equal -H.
Outside the core, M is constant. Its value enters the unsquared outer-wall
Israel equation and determines H. Solving that coupled scalar root closes
the core and wall together. The existing common-middle interpolation uses
the new outer H and unchanged inner wall; both seams are checked.

The parent returns smoothly through its maximum areal radius and closes
at the core pole. The positive initial lapse collars realize both temporal
Israel jets. Removing f recovers the same-frequency shared reference, including its
extrinsic curvature, wall speed and proper length.

For the positive-phase coherent representative:

| Initial quantity | Result |
| --- | ---: |
| Isolated leading geometric mass | .0000724395411224 |
| Geometric mass including the inherited environment | .0000957140112889 |
| Inherited scalar cross contribution | .0000231104644327 |
| Proper slice energy excess | .00241200566080 |
| Outer-wall areal speed change | .000361840525751 |
| Parent maximum areal radius | 6.65064999805 |
| Total pole-to-pole proper length | 20.4717034369 |

The physical mass-unit conversion is 8pi M/kappa; it is not the proper
slice energy integral. Independent field quadrature and proper-distance
integration agree with the mass root within 5e-12. Sampled Hamilton
residuals are below 2e-8; middle and wall residuals are at floating precision.
This establishes initial-constraint compatibility, not the subsequent
nonlinear membrane evolution.
The finite cooling comparison also limits how precisely its coherent
representative is selected; the displayed solver digits are not a full
formation-history accuracy claim.
The changed wall speed compares two jointly prepared initial slices;
it is not a wall response propagated by the local cooling run.
The homogeneous reference's future shielding theorem is not asserted
for this new nonlinear core geometry.

## Relative phase has a geometric consequence

Changing f to -f keeps the isolated quadratic core energy fixed, but

\[
U(\phi_0+f)-U(\phi_0-f)=2U'_0f.
\]

The cross term reverses sign. Solving the same gravitational matching again
gives geometric mass .0000494933577268 and wall speed change .000187115382446.
Both configurations satisfy the initial constraints and full instantaneous
wall conditions. Their difference is coherent interference with the inherited
ambient field. A negative excess in part of the second profile is depletion
relative to that background; the canonical field's absolute energy remains
nonnegative. The gradient mass contribution also responds slightly to the
changed metric, so the total mass difference is not simply twice the cross
term evaluated on one fixed geometry.

This is why profile selection does not erase all state freedom. The local
cooling calculation can lose information to outgoing radiation while the
complete field state retains it, consistent with the earlier finite-flow
invertibility result. The environment and relative phase still matter.

## Reproduction and next connection

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_scalar_core.py --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_scalar_core.py
```

The [compact record](../results/summaries/nsc-atlas-scalar-core-v1.json)
contains the fixed cloud inputs, numerical comparisons, flux and energy
accounting, hierarchy, both phase embeddings and producer hashes. The
admitted producer budget is 240 CPU seconds including its workers.

An independent Grok audit checked the profile normalization, mass factors,
ambient cross term, both walls, zero-core recovery and short outflow
controls. The published representative uses the mass-selected amplitude
above; the unscaled amplitude .1 is also retained as a tested control.

The full formation history in the excited, shared evolving parent remains
unresolved. It must include ambient-field interaction, the outgoing
radiation, moving walls and eventual nonspherical disturbances. Cloud
origin, frequency selection, abundance and observational identification
also remain open. This entry supplies a physical profile-selection
candidate and compatible sourced initial geometry for that next problem;
it does not relabel P003 or generate a new cosmological constant.
In a single global formation history the outgoing scalar radiation and
the original cloud's gravitational contribution must also be included.
