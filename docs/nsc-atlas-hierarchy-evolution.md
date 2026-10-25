# The core remains localized while its inherited parent expands

The [scale hierarchy](nsc-atlas-scale-hierarchy.md) established compatible
initial data, including a reservoir that supplies the core's mass. The
next connection evolves that same high-frequency field preparation in a
leading classical envelope description, with the reservoir, self-gravity
and parent expansion retained. A finite localized mass persists; a control
without self-gravity disperses. The long run is not a carrier-resolved
Einstein--Klein--Gordon evolution.

## Rethinking and admitted calculation

The choices were a direct high-frequency relativistic grid, another static
supply estimate, and an envelope derived from the existing real scalar.
The first would resolve about a million carrier cycles per slow binding
time; the second would not answer persistence. Select the envelope and
check it against the full equations on a short, tractable weak-field family.
Keep a distinction between its finite numerical error and its model error.

The new state is the complex amplitude of the same real field and one parent
scale factor. The homogeneous amplitude is included in that total envelope;
it is not counted again as a separate field. Subtract it only for the local
Poisson source and the evolved radial perturbation. Source ownership stays
with the initial hierarchy; the original full relativistic solver remains
unchanged. Existing SciPy sine transforms advance radial dispersion; the
spherical gravitational kernel is implemented directly.

Three resolutions, separate time/outer controls, a no-gravity control,
two short envelope/full-equation comparisons and a full-equation refinement
use eleven independent cases in six processes. The initial budget is 1200
CPU seconds. Each completed run is checkpointed locally before validation.
The public record contains summaries and hashes, not full field histories.
Failure of conservation, weak-field comparison or outer-domain independence
would require a diagnosis before extending the evolution.

The real-scalar envelope reduction is attributed to the existing
[Fodor, Forgacs and Mezei source](https://arxiv.org/html/1007.0388).
[Schrodinger--Poisson descriptions in an expanding Friedmann background](https://arxiv.org/abs/1504.05968)
provide related prior art; its bibliographic record and abstract were
checked. The source normalization, excited closed parent, compensation and
fixed-physical-radius mass account below are derived for this atlas member.

## Same real field, a slow amplitude, and the actual excited parent

Fix the supplied member omega=1e9, zeta=.001 and half-depletion. Its initial
field and reservoir radius are those already solved by ScaleHierarchy.
The original parent has Lambda_v=9/32, excitation E_b=3/8 and spatial
curvature C=.0227765023676. Its normalized scale factor A starts at one;
the original curvature radius starts at 6.62607707610.

Use the total complex classical amplitude psi:

\[
\chi=\sqrt{2/\kappa}\,\zeta^2 A^{-3/2}
\operatorname{Re}[\psi(x,\tau)e^{-i\omega t}],\qquad
x=\frac{\omega\zeta R}{A},\quad\tau=\omega\zeta^2t.
\]

The background amplitude in this normalization is

\[
b=\frac{\sqrt{\kappa/2}\,d}{\omega\zeta^2}
=.000612372435696,\qquad d=\sqrt{3}/2.
\]

The field evolved numerically is u=x(psi-b). The initial psi equals the
existing tapered core plus b times the depleted ambient profile. This
matches the real scalar value on the exact constrained slice. The momentum
Pi=0 matches only at leading order: the A^(-3/2) prefactor has an initial
derivative -3H/2. The relative carrier-quadrature correction is O(H/omega).
The exact initial gradient contribution to geometric mass is also retained
as an explicit GR/envelope difference rather than tuned away.

Carrier averaging gives the parent equations

\[
H^2=\frac{\Lambda_v}{3}+\frac{\kappa E_b}{3A^3}-\frac{C}{A^2},\qquad
\dot H=\frac{C}{A^2}-\frac{\kappa E_b}{2A^3}.
\]

The excitation is initially 4/7 of total density, so replacing this parent
by pure de Sitter would lose a leading contribution. The local envelope is

\[
i\psi_\tau=-\frac{\Delta_x\psi}{2A^2}+V\psi,\qquad
\Delta_xV=\frac{|\psi|^2-b^2}{2A}.
\]

The factors of A follow from the physical volume and field normalization.
The full quadratic contrast is retained, including the coherent cross term:
|psi|^2-b^2=|psi-b|^2+2b Re(psi-b). Near this core the quadratic local term
dominates; a linear-density treatment would be a different approximation.
No probability postulate or additional quantum field has been introduced.

## An exact normalization check on the averaged background

For chi_b proportional to A^(-3/2) cos(omega t), differentiation gives

\[
(\partial_t^2+3H\partial_t+\omega^2)\chi_b
=\left(-\frac{3}{2}\dot H-\frac{9}{4}H^2\right)\chi_b
=\frac{3}{4}\left(\frac{C}{A^2}-\Lambda_v\right)\chi_b.
\]

The dust contribution cancels exactly. Thus constant b is the leading WKB
background after factoring out A^(-3/2), with a known curvature/vacuum
residual. It is not a second, frozen physical background added to chi.
This identity is on the averaged parent; it does not prove equality to the
unaveraged scalar-gravity solution. The coefficient divided by omega^2 is
about 2e-19 for this member. A short full-equation comparison below tests the
combined approximation in a less extreme regime.

## Mass transport and expansion work

The geometric excess mass in a physical ball is, at this order,

\[
\mathfrak{M}_{\rm env}(R)
=\frac{\zeta}{2\omega}
\int_0^{\omega\zeta R/A}x^2(|\psi|^2-b^2)\,dx.
\]

This factor of one-half follows directly from the spherical Einstein mass
normalization. Initially the envelope core mass differs from the exact
initial GR mass by -3.31923e-7 relatively. The residual total compensated
envelope mass has the same small fraction; it is preserved in the account.

Continuity is

\[
\partial_\tau|\psi|^2
+\frac{1}{x^2}\partial_x\left[
\frac{x^2}{A^2}\operatorname{Im}(\psi^*\psi_x)\right]=0.
\]

For a fixed physical radius, x_*=X/A moves inward. With h=A'/A,

\[
\frac{dN(x_*)}{d\tau}
=-\frac{x_*^2}{A^2}\operatorname{Im}(\psi^*\psi_x)
-h x_*^3(|\psi|^2-b^2).
\]

The second term is essential; a fixed comoving radius asks a different
question. The calculation independently integrates this current at physical
cuts X=20 and 40, then compares it with the density integral.

Define the radial energy functional

\[
T=\frac{1}{2A^2}\int x^2|\psi_x|^2dx,\qquad
W=\frac{1}{2}\int x^2V(|\psi|^2-b^2)dx.
\]

With the symmetric Newtonian Green kernel and fixed outer Dirichlet value
psi=b, the energy has the exact evolution law within these reduced equations:

\[
\frac{d(T+W)}{d\tau}=-\frac{A'}{A}(2T+W).
\]

Both kinetic and gravitational terms contribute to expansion work. The
boundary's fixed field has zero time derivative, so the associated energy
boundary work vanishes; its mass current need not vanish and is measured.
There is no absorber. A larger outer domain checks the finite boundary.

## Finite evolution and controls

The [compact record](../results/summaries/nsc-atlas-hierarchy-evolution-v1.json)
contains the complete protocol and comparisons. Evolve to tau 600, or atlas
time .6, and average the late interval tau 480--600. The parent expands to
A=1.27443638723. Mass within the original fixed physical core cutoff X=20
has late mean ratio

\[
\frac{\langle\mathfrak{M}_{\rm core}\rangle}
{\mathfrak{M}_{\rm core}(0)}=1.00374088\pm.00001314.
\]

The mass exceeds its initial value by a small inward transfer; it has not
been renormalized. With self-gravity removed, the same late ratio is
.000190975. The late mean physical half-mass radius is 2.23255821 in the
scaled length omega*zeta*R, with comparison envelope .00319767. The core
therefore stays localized while the reservoir and homogeneous parent expand. Its radius is not stationary: the early peak-to-peak range
is .0673906 in scaled length, increasing to .1399190 in the late interval.
The recorded growth of breathing prevents an asymptotic-stability claim.
These statements concern a finite spherical evolution within the envelope.

The mass-current discrepancy falls from 5.13e-6 to 1.26e-6 to 3.24e-7 of
initial core mass. The energy/expansion-work discrepancy falls from 1.75e-5
to 4.19e-6 to 1.04e-6 of initial energy. Total mass agrees with integrated
outer current to 3.15e-9 of initial core mass. Outer net transport is
8.66e-7 of initial core mass; enlarging the domain changes late core mass
by about 1.3e-11. Reported numerical ranges use twice the largest grid,
time-step or outer-domain change, not a rigorous interval enclosure.

## Comparison with the full relativistic equations

Hold omega*zeta^2=12.8 and evolve to atlas time .4 (slow time 5.12).
The initial scalar field and solved reservoir are supplied by the same
hierarchy constructor. The original full EKG functions remain unchanged.

| zeta | omega | Full EKG core-mass ratio | Envelope ratio | Difference |
| --- | ---: | ---: | ---: | ---: |
| .08 | 2000 | 1.007402025 | 1.007139302 | .000262723 |
| .04 | 8000 | 1.007190375 | 1.007121981 | .0000683932 |

Halving zeta reduces the discrepancy by about 3.84, consistent with the
expected weak-field ordering. Refining the smaller full-EKG case changes
its ratio by 2.15e-8, far less than the model discrepancy. The finite
boundary margins remain positive. Free expanding Gaussian propagation,
an independent Gaussian potential integral and moving-cut continuity
provide separate analytic checks.

This comparison does **not** establish a full-GR error bound at tau 600.
The long-run member has zeta^2=1e-6, initial H/omega=4.43e-10 and initial
parent spatial curvature C Rd^2=5.29e-9. These identify omitted scales;
their smallness alone is not a rigorous long-time error estimate.

## Conditional physical interpretation and next question

Under the previous, explicitly supplied Lambda_v=1e-52 m^-2 conversion,
the duration is about 3.36 billion years. The fixed mass cut is 34.4 kpc;
the late mean half-mass radius is about 3.84 kpc. No galaxy data were fitted.
The physical time and radii inherit the assumed ruler and approximation.

The result advances initial compatibility to finite radial localization.
It does not derive the initial cloud, select its abundance, prove asymptotic
or nonspherical stability, or transfer the older lens image predictions.
The next strong discriminator is a perturbation/persistence calculation
on this evolved state, or a joint observable calculation with the same
source law, rather than more refinement of this completed radial run.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_evolution.py --output results/summaries/nsc-atlas-hierarchy-evolution-v1.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_hierarchy_evolution.py
```

A focused Grok audit independently checked the total-envelope normalization,
the geometric-mass factor, both expansion-work terms, the constant-background
WKB residual and the initial-momentum limitation. The symbolic identities,
analytic reference tests and short full-equation comparison provide the
reproducible checks. The eleven underlying computations used about 629 CPU
seconds, within the 1200-second budget; cache-only postprocessing does not
repeat those computations.
